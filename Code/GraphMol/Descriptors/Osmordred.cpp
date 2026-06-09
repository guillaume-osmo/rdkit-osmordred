//  Copyright (c) 2025, Guillaume Godin Osmo Labs, PBC’s and others
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "Osmordred.h"
#include <boost/functional/hash.hpp> // For custom hashing of pairs
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Atom.h>
#include <RDGeneral/export.h>
#include <RDGeneral/types.h>


#include <boost/graph/adjacency_list.hpp>

#include <set>
#include <cmath> // For M_PI and pow
#include <tuple>
#include <map>
#include <string>
#include <utility> // for std::pair
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <climits>
#include <queue> // for fused rings
#include <stdexcept>
#include <iomanip>  // For std::fixed and std::setprecision
#include <sstream>  // For std::ostringstream
#include <iostream>
#include <cstring> // For memcpy
#include <functional>
#include <numeric>
#include <stack>

#ifdef RDK_BUILD_OSMORDRED_SUPPORT
#if defined(_MSC_VER) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#endif

#include <Eigen/Dense> // we should try to remove those...
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/PartialCharges/GasteigerParams.h>
#include <GraphMol/RingInfo.h>
#include <RDGeneral/RDThreads.h>
#include <lapacke.h>
#include <future>
#include <chrono>
#endif

// Define a custom hash function for std::pair<int, int>
namespace std {
    template <>
    struct hash<std::pair<int, int>> {
        size_t operator()(const std::pair<int, int>& p) const {
            // Combine the hash of the two elements of the pair
            return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
        }
    };
}

namespace RDKit {
namespace Descriptors {
namespace Osmordred {
#ifdef RDK_BUILD_OSMORDRED_SUPPORT
bool hasOsmordredSupport() { return true; }

// v2.0: Filter function to check if a molecule is too large (will cause hangs)
// Returns true if molecule has >10 rings OR >200 heavy atoms
// This prevents Osmordred.Calculate from hanging on very complex molecules
bool isMoleculeTooLarge(const RDKit::ROMol &mol) {
  const RDKit::RingInfo *ri = mol.getRingInfo();
  int numRings = 0;
  if (ri && ri->isInitialized()) {
    numRings = ri->numRings();
  }
  
  int numHeavyAtoms = mol.getNumHeavyAtoms();
  
  // Filter: >10 rings OR >200 heavy atoms
  return (numRings > 10 || numHeavyAtoms > 200);
}

namespace {    
void solveLinearSystem(const ROMol &mol, std::vector<double>& A, std::vector<double>& B,
		       int n, int nrhs, bool& success) {
    int lda = n; // Leading dimension of A
    int ldb = n; // Leading dimension of B
    int info;

    success = false; // Initialize success flag

    // CRITICAL FIX v2.0: Save original RHS before any LAPACK calls modify B
    // LAPACK routines modify B in-place, even when they fail!
    std::vector<double> B_original = B;

    // First, try dposv (Cholesky factorization for positive definite matrices)
    std::vector<double> A_copy = A; // Copy A because LAPACK modifies it
    info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', n, nrhs, A_copy.data(), lda, B.data(), ldb);

    if (info == 0) {
        success = true;
        return;
    } else {
        // dposv failed; fall back to dgesv (LU factorization)
        // CRITICAL FIX v2.0: Restore original RHS before calling dgesv
        // dposv modified B even though it failed!
        B = B_original;
        
        std::vector<int> ipiv(n); // Pivot array for dgesv
        A_copy = A; // Reset A because it was modified by dposv
        info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, A_copy.data(), lda, ipiv.data(), B.data(), ldb);

        if (info == 0) {
            success = true;
            return;
        } else {
            // dgesv failed (singular matrix); fall back to dgelss (pseudo-inverse via SVD)
            // CRITICAL FIX v2.0: Added dgelss fallback for singular matrices
            // This provides a minimum-norm least-squares solution when exact solution doesn't exist
            B = B_original; // Restore original RHS values
            
            std::vector<double> A_copy2 = A; // Fresh copy for dgelss
            std::vector<double> B_copy = B; // Copy B because dgelss modifies it
            
            // Allocate workspace for dgelss
            std::vector<double> s(n); // Singular values
            int rank; // Rank of matrix
            double rcond = 1e-15; // Condition number threshold
            
            // dgelss computes least-squares solution: min ||Ax - b||_2
            info = LAPACKE_dgelss(LAPACK_COL_MAJOR, n, n, nrhs,
                                   A_copy2.data(), lda, B_copy.data(), ldb,
                                   s.data(), rcond, &rank);
            
            if (info == 0) {
                // dgelss succeeded - copy solution back to B
                B = B_copy;
                success = true;
                return;
            } else {
                // All solvers failed - this is a true error
                std::string outputSmiles = RDKit::MolToSmiles(mol);
                std::cerr << "ERROR: All LAPACK solvers failed (dposv, dgesv, dgelss): info=" 
                          << info << ", Smiles:" << outputSmiles << "\n";
            }
        }
    }
}
}
  
// Function to count the number of endocyclic single bonds
int calcEndocyclicSingleBonds(const RDKit::ROMol &mol) {
  const RDKit::RingInfo *ri = mol.getRingInfo();
  if (!ri || !ri->isInitialized()) {
    return 0;  // No ring information available
  }

  std::unordered_set<int> bondIndices;

  // Collect all bond indices involved in rings
  for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
    if (ri->numBondRings(i) > 0) {  // Check if the bond is part of a ring
      bondIndices.insert(i);
    }
  }

  int nbonds = 0;
  for (const auto &bondIdx : bondIndices) {
    const RDKit::Bond *bond = mol.getBondWithIdx(bondIdx);
    if (bond->getBondType() == RDKit::Bond::SINGLE) {
      nbonds++;
    }
  }

  return nbonds;
}

static const std::vector<std::string> acidicSMARTS = {
    "[O;H1]-[C,S,P]=O", "[*;-;!$(*~[*;+])]", "[NH](S(=O)=O)C(F)(F)F",
    "n1nnnc1"};

static const std::vector<std::shared_ptr<RDKit::RWMol>> &GetAcidicSmarts() {
  static const std::vector<std::shared_ptr<RDKit::RWMol>> compiledAcidicSMARTS = [] {
    std::vector<std::shared_ptr<RDKit::RWMol>> res;
    for (const auto &smarts : acidicSMARTS) {
      auto mol = RDKit::SmartsToMol(smarts);
      if (mol) {
	res.emplace_back(std::shared_ptr<RDKit::RWMol>(mol));
      } else {
	std::cerr << "Invalid SMARTS: " << smarts << std::endl;
      }
    }
    return res;
  }();
  return compiledAcidicSMARTS;
}


static const std::vector<std::string> basicSMARTS = {
    "[NH2]-[CX4]",
    "[NH](-[CX4])-[CX4]",
    "N(-[CX4])(-[CX4])-[CX4]",
    "[*;+;!$(*~[*;-])]",
    "N=C-N",
    "N-C=N"
};

static const std::vector<std::shared_ptr<RDKit::RWMol>>& GetBasicSmarts() {
      static const std::vector<std::shared_ptr<RDKit::RWMol>>
          compiledBasicSMARTS = [] {
            std::vector<std::shared_ptr<RDKit::RWMol>> res;
            for (const auto &smarts : basicSMARTS) {
              auto mol = RDKit::SmartsToMol(smarts);
              if (mol) {
                res.emplace_back(std::shared_ptr<RDKit::RWMol>(mol));
              } else {
                std::cerr << "Invalid SMARTS: " << smarts << std::endl;
              }
            }
            return res;
          }();
  return compiledBasicSMARTS;
}
      

    // Function to count substructure matches
    int countMatches(const ROMol& mol, const std::vector<std::shared_ptr<RDKit::RWMol>>& patterns) {
        int count = 0;
        for (const auto& pattern : patterns) {
            if (pattern) {
                std::vector<RDKit::MatchVectType> matches;
                RDKit::SubstructMatch(mol, *pattern, matches);
                count += matches.size();
            }
        }
        return count;
    }



// additional features

    static const std::vector<std::string> alcoholsSMARTS = {
        "[#6;H2;!$(C=O)][OX2H]",
        "[#6;H1;!$(C=O)][OX2H]",
        "[#6;H0;!$(C=O)][OX2H]",
    };


 // Precompile SMARTS patterns for efficiency
    static const std::vector<std::shared_ptr<RDKit::RWMol>> &GetAlcoholSmarts() {
      static const std::vector<std::shared_ptr<RDKit::RWMol>>
          compiledalcoholsSMARTS = [] {
            std::vector<std::shared_ptr<RDKit::RWMol>> res;
            for (const auto &smarts : alcoholsSMARTS) {
              auto mol = RDKit::SmartsToMol(smarts);
              if (mol) {
                res.emplace_back(std::shared_ptr<RDKit::RWMol>(mol));
              } else {
                std::cerr << "Invalid SMARTS: " << smarts << std::endl;
              }
            }
            return res;
          }();
      return compiledalcoholsSMARTS;
    }

    // Function to count primary, secondary, and tertiary hydroxyl groups
std::vector<int> countHydroxylGroups(const RDKit::ROMol &mol) {
    std::vector<int> results(3, 0);

    for (size_t i = 0; i < GetAlcoholSmarts().size(); ++i) {
        std::vector<RDKit::MatchVectType> matches;
      RDKit::SubstructMatch(mol, *GetAlcoholSmarts()[i], matches);
        results[i] = matches.size();
    }

    return results;
}


// Define static SMARTS patterns for bridged bonds, polyacids, and polyalcohols
static const std::vector<std::string> smartsPatterns = {
    "[*x3,*x4,*x5,*x6]",            // Bridged bonds
    "[CX3](=O)[OX1H0-,OX2H1]",       // Polyacid
    "[#6;!$(C=O)][OX2H]"             // Polyalcohol
};

// Precompile SMARTS patterns for efficiency
static const std::vector<std::shared_ptr<RDKit::RWMol>> &GetSmarts() {
  static const std::vector<std::shared_ptr<RDKit::RWMol>> compiledSMARTS = [] {
    std::vector<std::shared_ptr<RDKit::RWMol>> res;
    for (const auto &smarts : smartsPatterns) {
      auto mol = RDKit::SmartsToMol(smarts);
      if (mol) {
        res.emplace_back(std::shared_ptr<RDKit::RWMol>(mol));
      } else {
        std::cerr << "Invalid SMARTS: " << smarts << std::endl;
      }
    }
    return res;
  }();
  return compiledSMARTS;
}

// Function to count the number of bridged bonds
int countBridgedBonds(const RDKit::ROMol &mol) {
    int nBridgeheads = RDKit::Descriptors::calcNumBridgeheadAtoms(mol);
    if (nBridgeheads > 0) {
        std::vector<RDKit::MatchVectType> matches;
        int nbonds = 0;

        if (RDKit::SubstructMatch(mol, *GetSmarts()[0], matches)) {  // Use precompiled pattern for bridged bonds
            for (const auto &match : matches) {
                for (const auto &atom : match) {
                    nbonds += std::distance(mol.getAtomNeighbors(mol.getAtomWithIdx(atom.second)).first,
                                            mol.getAtomNeighbors(mol.getAtomWithIdx(atom.second)).second);
                }
            }
        }
        return nbonds;
    }
    return 0;
}

// Function to check if a molecule is a polyacid
bool isPolyAcid(const RDKit::ROMol &mol) {
    return RDKit::SubstructMatch(mol, *GetSmarts()[1]).size() > 1;  // Use precompiled pattern for polyacid
}

// Function to check if a molecule is a polyalcohol
bool isPolyAlcohol(const RDKit::ROMol &mol) {
  return RDKit::SubstructMatch(mol, *GetSmarts()[2]).size() >
         1;  // Use precompiled pattern for polyalcohol
}



std::vector<double> calcAddFeatures(const RDKit::ROMol& mol) {
        std::vector<double> v(7,0.);
        auto hydroxylCounts = countHydroxylGroups(mol);
        v[0] = hydroxylCounts[0];
        v[1] = hydroxylCounts[1];
        v[2] = hydroxylCounts[2];
        v[3] = static_cast<double>(countBridgedBonds(mol));
        v[4] = static_cast<double>(isPolyAcid(mol));
        v[5] = static_cast<double>(isPolyAlcohol(mol));
        v[6] = static_cast<double>(calcEndocyclicSingleBonds(mol));

    return v;
}



    // Function to calculate the number of acidic groups in a molecule
    int calcAcidicGroupCount(const ROMol& mol) {
        // v3 fix: count acidic groups (carboxylic/sulfonic/phosphonic OH,
        // anions, triflylamide, tetrazole), not alcohols. The previous
        // GetAlcoholSmarts() counted phenols/alcohols (Cyanidin -> 5) and
        // missed COOH (Glutathione -> 0). GetAcidicSmarts() matches Mordred.
        return countMatches(mol, GetAcidicSmarts());
    }

    // Function to calculate the number of basic groups in a molecule
    int calcBasicGroupCount(const ROMol& mol) {
        return countMatches(mol, GetBasicSmarts());
    }

    // Function to calculate both acidic and basic group counts
    std::vector<int> calcAcidBase(const ROMol& mol) {
        return {calcAcidicGroupCount(mol), calcBasicGroupCount(mol)};
    }

    std::vector<double> calcABCIndex(const ROMol& mol) {
        std::vector<double> res(2,0.);
        double ggAbcIndex = 0.0;
        double abcIndex = 0.0;

        double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"
        unsigned int numAtoms = mol.getNumAtoms();

        for (const auto& bond : mol.bonds()) {
            int u = bond->getBeginAtomIdx();
            int v = bond->getEndAtomIdx();
            auto atom1 = bond->getBeginAtom();
            auto atom2 = bond->getEndAtom();
            int nu = 0, nv = 0;

            double du = static_cast<double>(atom1->getDegree());
            double dv = static_cast<double>(atom2->getDegree());

            abcIndex += std::sqrt((du + dv - 2.0) / (du * dv));

            for (size_t i = 0; i < numAtoms; ++i) {
                if (distances[u * numAtoms + i] < distances[v * numAtoms + i]) nu++;
                if (distances[v * numAtoms + i] < distances[u * numAtoms + i]) nv++;
            }

            ggAbcIndex += std::sqrt((nu + nv - 2.0) / (nu * nv));
        }
        res[0] = abcIndex;
        res[1] = ggAbcIndex;

        return res;
    }


    // Function to calculate the number of aromatic atoms in a molecule
    int countAromaticAtoms(const ROMol& mol) {
        int count = 0;
        for (const auto& atom : mol.atoms()) {
            if (atom->getIsAromatic()) {
                count++;
            }
        }
        return count;
    }

    // Function to calculate the number of aromatic bonds in a molecule
    int countAromaticBonds(const ROMol& mol) {
        int count = 0;
        for (const auto& bond : mol.bonds()) {
            if (bond->getIsAromatic()) {
                count++;
            }
        }
        return count;
    }


    std::vector<int> calcAromatic(const ROMol& mol) {
            return {countAromaticAtoms(mol), countAromaticBonds(mol)};
    }


    static const std::unordered_map<std::string, int> elementMapAtomCounts = {
        {"H", 1}, {"B", 2}, {"C", 3}, {"N", 4}, {"O", 5},
        {"S", 6}, {"P", 7}, {"F", 8}, {"Cl", 9}, {"Br", 10}, {"I", 11}
    };

    // Function to calculate the atom count descriptor
    std::vector<int> calcAtomCounts(const ROMol& mol) {
        std::unique_ptr<RDKit::ROMol> hmol(RDKit::MolOps::addHs(mol));

        // Initialize the counts for each atom type

        // nAtom,nHeavyAtom,nSpiro,nBridgehead,nHetero,nH,nB,nC,nN,nO,nS,nP,nF,nCl,nBr,nI,nX version higher then v1 else
        // nAtom,nHeavyAtom,nSpiro,nBridgehead,nH,nB,nC,nN,nO,nS,nP,nF,nCl,nBr,nI,nX for version 1

        //std::vector<int> counts(17, 0);

        int nAtoms  = hmol->getNumAtoms();
        int nHeavy = RDKit::Descriptors::calcNumHeavyAtoms(mol);
        int nSpiro = RDKit::Descriptors::calcNumSpiroAtoms(mol);
        int nBrigde = RDKit::Descriptors::calcNumBridgeheadAtoms(mol);


        // Halogen list (F, Cl, Br, I) not use in the logic ... also restrinction
        // std::vector<int> halogens = {9, 17, 35, 53};  // Atomic numbers of halogens (F, Cl, Br, I)  remark: also 85, 117 in Mordred, but honesttly we have rarely the case ...
        int nH = 0, nB = 0, nC = 0, nN = 0, nO = 0, nS = 0, nP = 0, nF = 0, nCl = 0, nBr = 0, nI = 0;

        // Iterate over all atoms in the molecule
        for (const auto& atom : hmol->atoms()) {
            std::string symbol = atom->getSymbol();

            auto it = elementMapAtomCounts.find(symbol);
            if (it != elementMapAtomCounts.end()) {
                switch (it->second) {
                    case 1:  nH++;  break;
                    case 2:  nB++;  break;
                    case 3:  nC++;  break;
                    case 4:  nN++;  break;
                    case 5:  nO++;  break;
                    case 6:  nS++;  break;
                    case 7:  nP++;  break;
                    case 8:  nF++;  break;
                    case 9:  nCl++; break;
                    case 10: nBr++; break;
                    case 11: nI++;  break;
                }
            }
        }

        int nX = nF+nCl+nBr+nI;
        int nHetero = RDKit::Descriptors::calcNumHeteroatoms(mol);
        return  {nAtoms,nHeavy,nSpiro,nBrigde,nHetero, nH,nB,nC,nN,nO,nS,nP,nF,nCl,nBr,nI,nX };
    }

    // return vector sum over rows
    std::vector<double> _VertexDegrees(const double* distances, const unsigned int numatoms) {
        std::vector<double> res(numatoms, 0.0);
        double sum;
        for (unsigned int i = 0; i < numatoms ; ++i) {
            sum = 0.0;
            for (unsigned int j = 0; j < numatoms ; ++j) {
                    sum += distances[ j * numatoms + i];
            }
            res[i] = sum;
        }
        return res;
    }

    double BalabanJ(const ROMol& mol) {


        double q = mol.getNumBonds();
        unsigned int n = mol.getNumAtoms();

        double* Topodistances = MolOps::getDistanceMat(mol, false, false, false, "Balaban"); // we need bond order diagonal!  based on the code mordred all false
        double* adjMat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO" );

        std::vector<double> s = _VertexDegrees(Topodistances, n);
        double mu = q - n + 1;

        double sum_ = 0.0;
        for (unsigned int i = 0; i < n; ++i) {
            double si = s[i];
            for (unsigned int j = i; j < n; ++j) {
                if (adjMat[i * n + j] == 1) {
                    sum_ += 1.0 / sqrt(si * s[j]);
                }
            }
        }

        double J = (mu + 1 != 0) ? (q / (mu + 1)) * sum_ : 0.0;

        return J;
    }


 std::vector<double> calcBalabanJ(const ROMol& mol) {

    std::vector<double> res(1,0.);
    res[0] = BalabanJ(mol);
    return res;

 }


    // bertyCT related functions (InfoGain can be found in rdkit ML/InfoGainFuncs part, but I have to change to input for matching python code)
    template <typename... Args>
    std::string makeKey(Args... args) {
        std::ostringstream oss;
        ((oss << args << "_"), ...);
        std::string key = oss.str();
        key.pop_back();  // Remove the trailing underscore
        return key;
    }

    template <class T>
    double InfoEntropy(const std::vector<T>& data) {
        T nInstances = 0;
        double accum = 0.0, d;

        for (const auto& val : data) {
            nInstances += val;
        }

        if (nInstances != 0) {
            for (const auto& val : data) {
                d = static_cast<double>(val) / nInstances;
                if (d != 0) {
                    accum += -d * std::log(d);
                }
            }
        }
        return accum / std::log(2.0);
    }



    template <class T>
    double WeightedInfoEntropy(const std::vector<T>& data,const std::vector<double>& w) {
        T nInstances = 0;
        double accum = 0.0, d;

        for (const auto& val : data) {
            nInstances += val;
        }

        if (nInstances != 0) {
            for (size_t i = 0; i < data.size(); ++i) {
                const auto& val = data[i];
                const auto& wi = w[i];
                d = static_cast<double>(val) / nInstances;
                if (d != 0) {
                    accum += -d * std::log(d) * wi;
                }
            }
        }
        return accum / std::log(2.0);
    }




    template <class T>
    double WeightedCrossInfoEntropy(const std::vector<T>& data,const std::vector<T>& w) {
        T nInstances = 0;
        double accum = 0.0, d;

        for (const auto& val : data) {
            nInstances += val;
        }

        if (nInstances != 0) {
            for (size_t i = 0; i < data.size(); ++i) {
                const auto& val = data[i];
                const auto& wi = w[i]*data[i];
                d = static_cast<double>(val) / nInstances;
                if (d != 0) {
                    accum += -d * std::log(d) * wi;
                }
            }
        }
        return accum / std::log(2.0);
    }





    // Function to assign symmetry classes to each atom based on the distance matrix
    std::vector<int> assignSymmetryClasses(const RDKit::ROMol& mol, const std::vector<std::vector<double>>&, int numAtoms, int cutoff) {
        std::vector<int> symList(numAtoms, 0);

        double* distances = MolOps::getDistanceMat(mol, true, false, true, "Balaban");
        std::vector<std::vector<double>> distMatrix(numAtoms, std::vector<double>(numAtoms, 0.0));

        // Fill the distance matrix
        for (int i = 0; i < numAtoms; ++i) {
            for (int j = i; j < numAtoms; ++j) {
                distMatrix[i][j] = distances[i * numAtoms + j];
                distMatrix[j][i] = distMatrix[i][j];
            }
        }

        // To store unique symmetry classes
        std::unordered_map<std::string, int> keysSeen;
        int currentClass = 1;

        // Assign symmetry classes based on distances
        for (int i = 0; i < numAtoms; ++i) {
            std::vector<double> tmpList(distMatrix[i].begin(), distMatrix[i].end());
            std::sort(tmpList.begin(), tmpList.end());

            // Take the first 'cutoff' distances to form the unique key
            std::string key = "";
            for (int j = 0; j < std::min(cutoff, static_cast<int>(tmpList.size())); ++j) {
                key += std::to_string(tmpList[j]) + ",";
            }

            if (keysSeen.find(key) == keysSeen.end()) {
                keysSeen[key] = currentClass++;
            }

            symList[i] = keysSeen[key];
        }

        return symList;
    }


    double lookUpBondOrder(int atom1Id, int atom2Id, const std::unordered_map<std::pair<int, int>, Bond::BondType>& bondDict) {
        std::pair<int, int> theKey = std::minmax(atom1Id, atom2Id);
        auto it = bondDict.find(theKey);

        if (it == bondDict.end()) {
            return 1.0;  // Default bond order if not found
        }

        Bond::BondType bondOrder = it->second;
        if (bondOrder == Bond::AROMATIC) {
            return 1.5;  // Aromatic bond order
        }
        return static_cast<double>(bondOrder);  // Convert bond type to numeric value
    }

    // Function to create bondDict, neighborList, and vdList efficiently
    std::tuple<std::unordered_map<std::pair<int, int>, Bond::BondType>, std::vector<std::vector<int>>, std::vector<int>> CreateBondDictEtc(const ROMol& mol, int numAtoms) {
        std::unordered_map<std::pair<int, int>, Bond::BondType> bondDict;
        std::vector<std::vector<int>> nList(numAtoms);  // List of neighbors for each atom
        std::vector<int> vdList(numAtoms, 0);  // Valency list, initialized to 0

        // Iterate over bonds in the molecule
        for (const auto& bond : mol.bonds()) {
            int atom1 = bond->getBeginAtomIdx();
            int atom2 = bond->getEndAtomIdx();

            // Ensure atom1 < atom2 for consistent key order in bondDict
            if (atom1 > atom2) std::swap(atom1, atom2);

            // Add bond to bondDict (aromatic bonds are marked as AROMATIC)
            bondDict[std::make_pair(atom1, atom2)] = bond->getBondType();

            // Update neighbors for atom1
            if (std::find(nList[atom1].begin(), nList[atom1].end(), atom2) == nList[atom1].end()) {
                nList[atom1].push_back(atom2);
            }

            // Update neighbors for atom2
            if (std::find(nList[atom2].begin(), nList[atom2].end(), atom1) == nList[atom2].end()) {
                nList[atom2].push_back(atom1);
            }
        }

        // Calculate the valency (number of neighbors) for each atom
        for (int i = 0; i < numAtoms; ++i) {
            vdList[i] = nList[i].size();
        }
        return std::make_tuple(bondDict, nList, vdList);
    }

    // Function to calculate the entropies of connections and atom types
    double CalculateEntropies(
        const std::unordered_map<std::string, double>& connectionDict,
        const std::unordered_map<int, double>& atomTypeDict,
        int numAtoms) {

        // Extract connection values into a list
        std::vector<double> connectionList;
        for (const auto& [key, value] : connectionDict) {
            connectionList.push_back(value);
        }

        if (connectionList.empty() || atomTypeDict.empty()) return 0.0;

        // Total connections
        double totConnections = std::accumulate(connectionList.begin(), connectionList.end(), 0.0);

        // Calculate connection entropy
        double connectionIE = totConnections * (InfoEntropy(connectionList) + std::log(totConnections) / std::log(2.0));

        // Extract atom type values into a list
        std::vector<double> atomTypeList;
        for (const auto& [key, value] : atomTypeDict) {
            atomTypeList.push_back(value);
        }


        // Calculate atom type entropy
        double atomTypeIE = numAtoms * InfoEntropy(atomTypeList);

        return atomTypeIE + connectionIE;
    }




    // Main BertzCT function (refactored)
    double BertzCT(const ROMol& mol) {
         int cutoff = 100;
        std::unordered_map<int, double> atomTypeDict;  // Maps atom type to count
        std::unordered_map<std::string, double> connectionDict;  // Maps bond classes to count

        int numAtoms = mol.getNumAtoms();

        double* distances = MolOps::getDistanceMat(mol, true, false, true, "Balaban");
        std::vector<std::vector<double>> dMat(numAtoms, std::vector<double>(numAtoms, 0.0));

        // Fill the distance matrix
        for (int i = 0; i < numAtoms; ++i) {
            for (int j = i; j < numAtoms; ++j) {
                dMat[i][j] = distances[i * numAtoms + j];
                dMat[j][i] = dMat[i][j];
            }
        }

        if (numAtoms < 2) return 0.0;

        // Create bondDict, neighborList, and vdList
        auto [bondDict, neighborList, vdList] = CreateBondDictEtc(mol, numAtoms);
        // Assign symmetry classes
        auto symmetryClasses = assignSymmetryClasses(mol, dMat, numAtoms, cutoff);

        // Iterate over atoms to compute atomTypeDict and connectionDict
        for (int atomIdx = 0; atomIdx < numAtoms; ++atomIdx) {
            int hingeAtomNumber = mol.getAtomWithIdx(atomIdx)->getAtomicNum();
            atomTypeDict[hingeAtomNumber]++;

            int hingeAtomClass = symmetryClasses[atomIdx];
            int numNeighbors = vdList[atomIdx];

            for (int i = 0; i < numNeighbors; ++i) {
                int neighbor_iIdx = neighborList[atomIdx][i];
                int NiClass = symmetryClasses[neighbor_iIdx];
                double bond_i_order = lookUpBondOrder(atomIdx, neighbor_iIdx, bondDict);

                if (bond_i_order > 1 && neighbor_iIdx > atomIdx) {
                    double numConnections = bond_i_order * (bond_i_order - 1) / 2;
                    std::string connectionKey = makeKey(std::min(hingeAtomClass, NiClass), std::max(hingeAtomClass, NiClass));
                    connectionDict[connectionKey] += numConnections;
                }

                for (int j = i + 1; j < numNeighbors; ++j) {
                    int neighbor_jIdx = neighborList[atomIdx][j];
                    int NjClass = symmetryClasses[neighbor_jIdx];
                    double bond_j_order = lookUpBondOrder(atomIdx, neighbor_jIdx, bondDict);
                    double numConnections = bond_i_order * bond_j_order;
                    std::string connectionKey = makeKey(std::min(NiClass, NjClass), hingeAtomClass, std::max(NiClass, NjClass));
                    connectionDict[connectionKey] += numConnections;

                }
            }
        }

        if (connectionDict.empty()) return 0.0;
        if (atomTypeDict.empty()) return 0.0;

        // Calculate and return the final entropy-based complexity value
        return CalculateEntropies(connectionDict, atomTypeDict, numAtoms);
    }

    std::vector<double> calcBertzCT(const RDKit::ROMol& mol) {
        std::vector<double> res(1,0.);
        res[0] = BertzCT(mol);
        return res;
    }



    // bondCount


    std::vector<int> calcBondCounts(const RDKit::ROMol& mol) {
        // Vector to hold bond counts: [Any, Single, Double, Triple, Aromatic, Multiple]
        std::vector<int> bondCounts(9, 0);



        std::unique_ptr<RDKit::ROMol> hmol(RDKit::MolOps::addHs(mol));




        bondCounts[0] = hmol->getNumBonds();


        // Loop through all bonds in the molecule
        for (const auto* bond : hmol->bonds()) {
            bool isAromatic = bond->getIsAromatic();
            auto bondType = bond->getBondType();

            // Count each bond type

            if (bondType == RDKit::Bond::BondType::SINGLE) {
                ++bondCounts[2]; // Single bond
            }

            if (bondType == RDKit::Bond::BondType::DOUBLE) {
                ++bondCounts[3]; // Double bond
            }

            if (bondType == RDKit::Bond::BondType::TRIPLE) {
                ++bondCounts[4]; // Triple bond
            }

            if (isAromatic || bondType == RDKit::Bond::BondType::AROMATIC) {
                ++bondCounts[5]; // Aromatic bond
            }

            if (isAromatic || bondType != RDKit::Bond::BondType::SINGLE) {
                ++bondCounts[6]; // Multiple bond
            }

            if (bond->getBeginAtom()->getAtomicNum() > 1 && bond->getEndAtom()->getAtomicNum() > 1) {
                ++bondCounts[1]; // Heavy bond
            }

        }


        RDKit::RWMol* kekulizedMol = new RDKit::RWMol(*hmol);
        RDKit::MolOps::Kekulize(*kekulizedMol, true);

        for (const auto* bond : kekulizedMol->bonds()) {
            auto bondType = bond->getBondType();

            if (bondType == RDKit::Bond::BondType::SINGLE) {
                ++bondCounts[7]; // Kekulize Single bond
            }

            if (bondType == RDKit::Bond::BondType::DOUBLE) {
                ++bondCounts[8]; // Kekulize Double bond
            }
        }

        delete kekulizedMol;

        return bondCounts;
    }





    // CarbonTypes there is an issue in the code not sure why this is not the same as in python code!
    // TODO: need debug
    std::vector<double> calcCarbonTypes(const RDKit::ROMol &mol) {
        int slots = 11;
        
        std::vector<double> carbonTypeCounts(slots, 0.);  // Store counts for C1SP1, C2SP1, C1SP2, ..., C4SP3
        unsigned int nSP3 = 0;  // Count of SP3 hybridized carbons
        unsigned int nSP2 = 0;  // Count of SP2 hybridized carbons

        // Make a copy of the molecule and Kekulize it
        // TO DO everywhere : must avoid segment fault so must replace all delete by this method !!!
        std::unique_ptr<RDKit::RWMol> kekulizedMol(new RDKit::RWMol(mol));
        RDKit::MolOps::Kekulize(*kekulizedMol, true);

        for (const RDKit::Atom* atom : kekulizedMol->atoms()) {
            // Check if the atom is a carbon (atomic number 6)
            if (atom->getAtomicNum() == 6) {


                // Count only neighboring carbons not the real degree (aka ignore any hetero atoms)
                int carbonNeighbors = 0;
                for (const auto& neighbor : kekulizedMol->atomNeighbors(atom)) {
                    if (neighbor->getAtomicNum() == 6) {
                        carbonNeighbors++;
                    }
                }

                RDKit::Atom::HybridizationType hybridization = atom->getHybridization();

                // Check SP1 hybridization
                if (hybridization == RDKit::Atom::SP) {
                    if (carbonNeighbors == 1) {
                        carbonTypeCounts[0]++; // C1SP1
                    } else if (carbonNeighbors == 2) {
                        carbonTypeCounts[1]++; // C2SP1
                    }
                }
                // Check SP2 hybridization
                else if (hybridization == RDKit::Atom::SP2) {
                    nSP2++;
                    if (carbonNeighbors == 1) {
                        carbonTypeCounts[2]++; // C1SP2
                    } else if (carbonNeighbors == 2) {
                        carbonTypeCounts[3]++; // C2SP2
                    } else if (carbonNeighbors == 3) {
                        carbonTypeCounts[4]++; // C3SP2
                    }
                }
                // Check SP3 hybridization
                else if (hybridization == RDKit::Atom::SP3) {
                    nSP3++;
                    if (carbonNeighbors == 1) {
                        carbonTypeCounts[5]++; // C1SP3
                    } else if (carbonNeighbors == 2) {
                        carbonTypeCounts[6]++; // C2SP3
                    } else if (carbonNeighbors == 3) {
                        carbonTypeCounts[7]++; // C3SP3
                    } else if (carbonNeighbors == 4) {
                        carbonTypeCounts[8]++; // C4SP3
                    }
                }
            }
        }

        // Calculate Hybridization Ratio (HybRatio)
        if (nSP2 + nSP3 > 0) {
             carbonTypeCounts[9] = static_cast<double>(nSP3) / (nSP2 + nSP3);
        }

        double fcps3 = RDKit::Descriptors::calcFractionCSP3(*kekulizedMol);
        carbonTypeCounts[10] = fcps3;
        return carbonTypeCounts;
    }






    // VertexAdjacencyInformation
    double VertexAdjacencyInformation(const RDKit::ROMol &mol)  {
        int m = 0;

        // Count the number of heavy-heavy bonds
        for (const auto &bond : mol.bonds()) {
            const auto *beginAtom = bond->getBeginAtom();
            const auto *endAtom = bond->getEndAtom();

            if (beginAtom->getAtomicNum() != 1 && endAtom->getAtomicNum() != 1) {
                ++m;
            }
        }


        // Calculate the descriptor value
        return 1.0 + std::log2(static_cast<double>(m));
    }




    std::vector<double> calcVertexAdjacencyInformation(const ROMol &mol) {
	
        std::vector<double> res(1,0.);
	res[0]=VertexAdjacencyInformation(mol);
        return res;
    }


    // WalkCount
    // TODO: not correct output not sure why to be investigation equation ...
    // Function to compute the adjacency matrix
    Eigen::MatrixXd computeAdjacencyMatrix(const RDKit::ROMol &mol, bool useBondOrder = false) {
        unsigned int nAtoms = mol.getNumAtoms();
        double* adjFlat = RDKit::MolOps::getAdjacencyMatrix(mol, useBondOrder, false, false, "noBO");

        Eigen::MatrixXd adjMat(nAtoms, nAtoms);

        for (unsigned int i = 0; i < nAtoms; ++i) {
            for (unsigned int j = 0; j < nAtoms; ++j) {
                adjMat(i, j) = adjFlat[i * nAtoms + j];
            }
        }
        return adjMat;
    }


    void upperTriangularSelfProduct(const Eigen::MatrixXd& A, Eigen::MatrixXd& result) {
        int n = A.rows();
        result = Eigen::MatrixXd::Zero(n, n);

        // Compute only upper triangle
        for (int i = 0; i < n; ++i) {
            for (int j = i; j < n; ++j) {
                for (int k = i; k <= j; ++k) {
                    result(i, j) += A(i, k) * A(j, k);
                }
                if (i != j) {
                    result(j, i) = result(i, j);  // Fill the symmetric part
                }
            }
        }
    }

    // Function to compute walk counts Need a Product Matrix order
    std::vector<double> calcWalkCounts(const RDKit::ROMol &mol) {
        const int maxOrder = 10;  // Maximum order of walks
        unsigned int nAtoms = mol.getNumAtoms();
        Eigen::MatrixXd adjMat = computeAdjacencyMatrix(mol);

        std::vector<double> results(21, 0.0);

        Eigen::MatrixXd powerMatrix = Eigen::MatrixXd::Identity(nAtoms, nAtoms);

        // we need to initialize at nAtoms both totals
        double totalMWC10 = nAtoms, totalSRW10 = nAtoms;

        //#pragma omp parallel for reduction(+:totalMWC10, totalSRW10)
        for (int order = 1; order <= maxOrder; ++order) {

            if (order == 1) {
                powerMatrix = adjMat;  // A^1
            } else {
                powerMatrix = powerMatrix * adjMat;  // A^order
            }


            // Compute MWC ie full
            double mwc = (order == 1) ? 0.5 * adjMat.sum() : std::log(powerMatrix.sum() + 1.0);
            results[order - 1] = mwc;  // Store MWC0n

            // Compute SRW ie trace
            double srw = std::log(powerMatrix.diagonal().sum() + 1.0);
            if (order > 1) {
                results[maxOrder + order -1] = srw;  // Store SRW0n
            }

            // Accumulate totals for MWC10 and SRW10
            totalMWC10 += mwc;
            totalSRW10 += srw;
        }

        results[maxOrder] = totalMWC10;  // TMWC10
        results[20] = totalSRW10;  // TSRW10

        return results;
    }



/*
    std::vector<double> calcWalkCountsBlas(const RDKit::ROMol &mol) {
        const int maxOrder = 10;  // Maximum order of walks
        unsigned int nAtoms = mol.getNumAtoms();

        // Get adjacency matrix from RDKit
        std::vector<double> adjMat(nAtoms * nAtoms, 0.0);
        double* adjFlat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "noBO");
        std::copy(adjFlat, adjFlat + nAtoms * nAtoms, adjMat.begin());

        std::vector<double> powerMatrix(adjMat); // Start with A^1
        std::vector<double> results(21, 0.0);

        // Initialize totals with the number of atoms
        double totalMWC10 = nAtoms, totalSRW10 = nAtoms;

        for (int order = 1; order <= maxOrder; ++order) {
            if (order > 1) {
                // Compute powerMatrix = powerMatrix * adjMat using BLAS (dgemm: C = alpha*A*B + beta*C)
                std::vector<double> tempMatrix(nAtoms * nAtoms, 0.0);
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                            nAtoms, nAtoms, nAtoms,
                            1.0, powerMatrix.data(), nAtoms,
                                adjMat.data(), nAtoms,
                            0.0, tempMatrix.data(), nAtoms);
                powerMatrix = tempMatrix;
            }

            // Compute MWC (full matrix sum)
            double mwc = (order == 1) ? 0.5 * std::accumulate(adjMat.begin(), adjMat.end(), 0.0)
                                    : std::log(std::accumulate(powerMatrix.begin(), powerMatrix.end(), 0.0) + 1.0);
            results[order - 1] = mwc;

            // Compute SRW (sum of diagonal elements)
            double srw = 0.0;
            for (unsigned int i = 0; i < nAtoms; ++i) {
                srw += powerMatrix[i * nAtoms + i];
            }
            srw = std::log(srw + 1.0);

            if (order > 1) {
                results[maxOrder + order - 1] = srw;
            }

            // Accumulate totals
            totalMWC10 += mwc;
            totalSRW10 += srw;
        }

        results[maxOrder] = totalMWC10;  // Store TMWC10
        results[20] = totalSRW10;  // Store TSRW10

        return results;
    }
*/
   // Weight
   // trick is to add the Hs for the average not only heavy atoms!
   // we need a function that can do this trick!!!
    std::vector<double>  calcWeight(const ROMol& mol) {
        std::vector<double> W(2, 0.);
        std::unique_ptr<RDKit::ROMol> hmol(RDKit::MolOps::addHs(mol));
        W[0] = RDKit::Descriptors::calcExactMW(mol);
        int fullatomsnumber = hmol->getNumAtoms();
        W[1] = W[0] / fullatomsnumber;
        return W;
    }

    // Wiener Index
    // working!
    std::vector<int>  calcWienerIndex(const ROMol& mol) {
        std::vector<int> WI(2, 0);

        double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"
        unsigned int numAtoms = mol.getNumAtoms();

        for (size_t i = 0; i < numAtoms; ++i) {
            for (size_t j = 0; j < numAtoms; ++j) {
                if ( distances[i * numAtoms + j] == 3) {
                    WI[1]+=1.;
                }
                WI[0]+=distances[i * numAtoms + j];
            }
        }
        // must return the half part (of the Distance matrix!) or only compute half matrix... to speed up!!!!
        WI[0] = static_cast<int>(WI[0]*0.5);
        WI[1] = static_cast<int>(WI[1]*0.5);

        return WI;
    }

    // Perform eigen decomposition on a symmetric matrix
    std::pair<std::vector<double>, std::vector<std::vector<double>>> eigenDecompositionSymmetric(
        const std::vector<std::vector<double>>& matrix) {

        int n = matrix.size();
        assert(matrix.size() == matrix[0].size() && "Matrix must be square");

        // Convert std::vector<std::vector<double>> to a 1D array in column-major order for LAPACK
        std::vector<double> matrixData(n * n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                matrixData[j * n + i] = matrix[i][j]; // Column-major order
            }
        }

        // Storage for eigenvalues
        std::vector<double> eigenValues(n);

        // Call LAPACK's dsyev to compute eigenvalues and eigenvectors
        int info;
        std::vector<double> work(1);
        int lwork = -1; // Request optimal workspace size
        info = LAPACKE_dsyev_work(LAPACK_COL_MAJOR, 'V', 'U', n, matrixData.data(), n, eigenValues.data(), work.data(), lwork);

        if (info != 0) {
            throw std::runtime_error("Error querying optimal workspace size for LAPACK dsyev");
        }

        lwork = static_cast<int>(work[0]);
        work.resize(lwork);

        // Perform eigen decomposition
        info = LAPACKE_dsyev_work(LAPACK_COL_MAJOR, 'V', 'U', n, matrixData.data(), n, eigenValues.data(), work.data(), lwork);
        if (info != 0) {
            throw std::runtime_error("LAPACK dsyev failed to compute eigen decomposition");
        }

        // Convert the eigenvectors back to std::vector<std::vector<double>>
        std::vector<std::vector<double>> eigenVectors(n, std::vector<double>(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                eigenVectors[i][j] = matrixData[j * n + i]; // Column-major to row-major
            }
        }

        return {eigenValues, eigenVectors};
    }



    // ZagrebIndex this is correct! we can use the Eigen too...
    std::vector<double> calcValence(const RDKit::ROMol &mol) {
        unsigned int nAtoms = mol.getNumAtoms();
        double* adjMat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "noBO");

        // Sum along each column to get the valence of each atom
        std::vector<double> valences(nAtoms, 0);
        for (unsigned int i = 0; i < nAtoms; ++i) {
            valences[i] = static_cast<double>(std::accumulate(&adjMat[i * nAtoms], &adjMat[(i + 1) * nAtoms], 0.0));
        }

        return valences;
    }






    // Zagreb
    // TODO: debug YES
    // this is correct but not matching the results myabe H's are added ???
    std::vector<double> calcZagrebIndex(const RDKit::ROMol &mol) {
        std::vector<double> zagrebIndex(4, 0.0);

        // Create a vector to store the degree of each atom
        std::vector<double> degrees = calcValence(mol);

        // Zagreb Index M1: sum of degrees raised to the power of 2 * lambda
        for (int degree : degrees) {
            zagrebIndex[0] += std::pow(degree, 2.);
            zagrebIndex[2] += std::pow(degree, -2.);
        }

        // Zagreb Index M2: sum of products of degrees of connected atoms raised to the power of lambda of (-1,1)

        for (const RDKit::Bond* bond : mol.bonds()) {
            zagrebIndex[1] += std::pow(degrees[bond->getBeginAtomIdx()] * degrees[bond->getEndAtomIdx()], 1.);
            zagrebIndex[3] += std::pow(degrees[bond->getBeginAtomIdx()] * degrees[bond->getEndAtomIdx()], -1.);
        }

        return zagrebIndex;
    }


    // Define Bondi radii for atomic contributions
    const std::unordered_map<int, double> bondiRadii = {
        {1, 1.20},  // H
        {5, 2.13},  // B
        {6, 1.70},  // C
        {7, 1.55},  // N
        {8, 1.52},  // O
        {9, 1.47},  // F
        {14, 2.10}, // Si
        {15, 1.80}, // P
        {16, 1.80}, // S
        {17, 1.75}, // Cl
        {35, 1.85}, // Br
        {33, 1.85}, // As
        {34, 1.90},  // Se
        {53, 1.98} // I missing in Mordred from paper/CDK implementation
    };

    // Precompute atomic contributions using 4/3 * pi * r^3
    const std::unordered_map<int, double> atomContributions = []() {
        std::unordered_map<int, double> contribs;
        for (const auto &[atomicNum, radius] : bondiRadii) {
            contribs[atomicNum] = (4.0 / 3.0) * M_PI * std::pow(radius, 3);
        }
        return contribs;
    }();

    // VdwVolumeABC
    // working "Need Hs explicit!"
    double VdwVolumeABC(const ROMol &mol) {


        std::unique_ptr<RDKit::ROMol> hmol(RDKit::MolOps::addHs(mol));

        //Nb is the number of bonds
        //NRa is the number of aromatic rings
        //NRA is the number of aliphatic rings
        int NRA = RDKit::Descriptors::calcNumAliphaticRings(*hmol);
        int NRa = RDKit::Descriptors::calcNumAromaticRings(*hmol);
        int Nb = hmol->getNumBonds();


        double ac = 0.0;

        // Calculate the sum of atomic contributions
        for (const auto &atom : hmol->atoms()) {
            int atomicNum = atom->getAtomicNum();
            auto it = atomContributions.find(atomicNum);
            if (it == atomContributions.end()) {
                throw std::runtime_error("Unknown atom type encountered in molecule.");
            }
            ac += it->second;
        }

        // Compute van der Waals volume
        return ac - 5.92 * Nb - 14.7 * NRa - 3.8 * NRA;
    }



    std::vector<double> calcVdwVolumeABC(const ROMol &mol) {
	
        std::vector<double> res(1,0.);
	res[0]=VdwVolumeABC(mol);
        return res;
    }

    // TopoPSA (adding S & P atoms effect to rdkit version!)
    // working for my molecules but need a S, P molecules to check ...
    // Helper function to calculate hydrogen count
    int hydrogenCount(const Atom* atom) {
        int nH = atom->getTotalNumHs();
        for (const auto& neighbor : atom->getOwningMol().atomNeighbors(atom)) {
            if (neighbor->getAtomicNum() == 1) {
                ++nH;  // Count the hydrogens attached to the atom
            }
        }
        return nH;
    }

    // Helper function to count bond types
    std::unordered_map<Bond::BondType, int> bondTypeCount(const Atom* atom) {
        std::unordered_map<Bond::BondType, int> bondCounts;
        for (const auto& bond : atom->getOwningMol().atomBonds(atom)) {
            if (bond->getIsAromatic()) {
                bondCounts[Bond::BondType::AROMATIC] += 1;
            } else {
                bondCounts[bond->getBondType()] += 1;
            }
        }
        return bondCounts;
    }



    // Function to calculate the phosphorus contribution to TPSA
    double getPhosphorusContribution(const Atom* atom) {
        int nH = hydrogenCount(atom);
        auto cnt = bondTypeCount(atom);

        if (atom->getFormalCharge() != 0 || atom->getIsAromatic()) {
            return 0.0;
        }

        if (nH == 1 && cnt == std::unordered_map<Bond::BondType, int>{{Bond::BondType::SINGLE, 3}, {Bond::BondType::DOUBLE, 1}}) {
            return 23.47;
        } else if (nH == 0) {
            if (cnt == std::unordered_map<Bond::BondType, int>{{Bond::BondType::SINGLE, 3}}) {
                return 13.59;
            } else if (cnt == std::unordered_map<Bond::BondType, int>{{Bond::BondType::SINGLE, 1}, {Bond::BondType::DOUBLE, 1}}) {
                return 34.14;
            } else if (cnt == std::unordered_map<Bond::BondType, int>{{Bond::BondType::SINGLE, 3}, {Bond::BondType::DOUBLE, 1}}) {
                return 9.81;
            }
        }
        return 0.0;
    }

    // Function to calculate the sulfur contribution to TPSA
    double getSulfurContribution(const Atom* atom) {
        int nH = hydrogenCount(atom);
        auto cnt = bondTypeCount(atom);

        if (atom->getFormalCharge() != 0) {
            return 0.0;
        }

        if (atom->getIsAromatic()) {
            if (nH == 0) {
                if (cnt == std::unordered_map<Bond::BondType, int>{{Bond::BondType::AROMATIC, 2}}) {
                    return 28.24;
                } else if (cnt == std::unordered_map<Bond::BondType, int>{{Bond::BondType::AROMATIC, 2}, {Bond::BondType::DOUBLE, 1}}) {
                    return 21.70;
                }
            }
        } else {
            if (nH == 1 && cnt == std::unordered_map<Bond::BondType, int>{{Bond::BondType::SINGLE, 2}}) {
                return 38.80;
            } else if (nH == 0) {
                if (cnt == std::unordered_map<Bond::BondType, int>{{Bond::BondType::SINGLE, 2}}) {
                    return 25.30;
                } else if (cnt == std::unordered_map<Bond::BondType, int>{{Bond::BondType::DOUBLE, 1}}) {
                    return 32.09;
                } else if (cnt == std::unordered_map<Bond::BondType, int>{{Bond::BondType::SINGLE, 2}, {Bond::BondType::DOUBLE, 1}}) {
                    return 19.21;
                } else if (cnt == std::unordered_map<Bond::BondType, int>{{Bond::BondType::SINGLE, 2}, {Bond::BondType::DOUBLE, 2}}) {
                    return 8.38;
                }
            }
        }
        return 0.0;
    }

    // Main function to calculate the Topological Polar Surface Area (TPSA)
    std::vector<double> calcTopoPSA(const ROMol& mol) {
        // v3 fix: res[0] = TopoPSA(NO) (N,O only); res[1] = TopoPSA (incl. S & P).
        // Use RDKit's calcTPSA includeSandP flag (Ertl), which correctly handles
        // thiols/sulfides/phosphorus on H-suppressed molecules. The previous
        // hand-rolled getSulfurContribution/getPhosphorusContribution assumed
        // explicit H and missed e.g. R-SH thiols (Glutathione lost the 38.80 term),
        // so res[1] silently equalled res[0] for thiols.
        // (getPhosphorusContribution/getSulfurContribution/hydrogenCount/
        //  bondTypeCount above are now unused and may be removed in cleanup.)
        std::vector<double> res(2, 0.0);
        res[0] = RDKit::Descriptors::calcTPSA(mol, false, false); // N,O only  -> TopoPSA(NO)
        res[1] = RDKit::Descriptors::calcTPSA(mol, false, true);  // incl. S,P -> TopoPSA
        return res;
    }


    // Main function to calculate the Topological Polar Surface Area (TPSA)
    std::vector<double> calcSLogP(const ROMol& mol) {
        std::vector<double> res(2, 0.0);
        double logP;
        double MR;
        RDKit::Descriptors::calcCrippenDescriptors(mol, logP, MR);

        res[0] = logP;
        res[1] = MR;

        return res;
    }



    // HygrogenBond from Rdkit code
    std::vector<double> calcHydrogenBond(const ROMol& mol) {
        std::vector<double> res(2, 0.);

        int nHBAcc = RDKit::Descriptors::calcNumHBA(mol);


        int nHBDon =  RDKit::Descriptors::calcNumHBD(mol);
        res[0] = static_cast<double>(nHBAcc);
        res[1] = static_cast<double>(nHBDon);

        return res;
    }

    // MOE need to implement EState here ;-)



    // Get principal quantum number for an atomic number
    int GetPrincipalQuantumNumber(int atomicNum) {
        if (atomicNum <= 2)   return 1;
        if (atomicNum <= 10)  return 2;
        if (atomicNum <= 18)  return 3;
        if (atomicNum <= 36)  return 4;
        if (atomicNum <= 54)  return 5;
        if (atomicNum <= 86)  return 6;
        return 7;
    }




    // rdkit version
    // caution use the Valence Method upper faster and cheaper!
    double getValenceElectrons(const Atom& atom) {
        unsigned int N = atom.getAtomicNum();
        if (N == 1) {
            return 0.0;  // Hydrogen has 0 valence electrons in this context
        }

        const PeriodicTable* tbl = PeriodicTable::getTable();
        double Zv = tbl->getNouterElecs(N) - atom.getFormalCharge();
        double Z = atom.getAtomicNum() - atom.getFormalCharge();
        unsigned int hi = atom.getTotalNumHs();


        unsigned int he = std::count_if(atom.getOwningMol().atomNeighbors(&atom).begin(), atom.getOwningMol().atomNeighbors(&atom).end(),
                                        [](const Atom* neighbor) { return neighbor->getAtomicNum() == 1; });
        unsigned int h = hi + he;

        return (Zv - h) / (Z - Zv - 1.);
    }



    double getSigmaElectrons(const Atom& atom) {
        // Retrieve the molecule owning the atom
        const ROMol& mol = atom.getOwningMol();

        // Get the neighbors of the atom
        auto neighbors = mol.atomNeighbors(&atom);

        // Count the number of neighbors that are not hydrogen
        unsigned int sigmaElectrons = std::count_if(neighbors.begin(), neighbors.end(),
            [](const Atom* neighbor) {
                return neighbor->getAtomicNum() != 1; // Exclude hydrogens
            });

        return static_cast<double>(sigmaElectrons);
    }

    // Function to get intrinsic state (equivalent to get_intrinsic_state)
    double getIntrinsicState(const Atom& atom) {
        unsigned int i = atom.getAtomicNum();
        double d = getSigmaElectrons(atom);
        double dv = getValenceElectrons(atom);

        if (d == 0) {
            return 0.;
        }

        // Use GetPrincipalQuantumNumber to get the principal quantum number
        int n = GetPrincipalQuantumNumber(i);

        return ((2.0 / n) * (2.0 / n) * dv + 1) / d;
    }




    double getRKHE(const Atom& atom) {
        const PeriodicTable* tbl = PeriodicTable::getTable();

        int hi = atom.getTotalNumHs();

        int Z = atom.getAtomicNum();

        int Sigma = atom.getDegree();

        int SigmaV = tbl->getNouterElecs(Z) - hi;

        int N = GetPrincipalQuantumNumber(Z);

        return (SigmaV - Sigma) / static_cast<double>(N * N);
    }



std::vector<double> calcIStateIndices(const RDKit::ROMol& mol){
        const RDKit::PeriodicTable* tbl = RDKit::PeriodicTable::getTable();
        int numAtoms = mol.getNumAtoms();

        std::vector<double> Is(numAtoms, 0.0);

        // Compute initial EState values
        for (int i = 0; i < numAtoms; ++i) {
            const auto* atom = mol.getAtomWithIdx(i);
            int degree = atom->getDegree();
            if (degree > 0) {
                int atomicNum = atom->getAtomicNum();
                int dv = tbl->getNouterElecs(atomicNum) - atom->getTotalNumHs();
                int N = GetPrincipalQuantumNumber(atomicNum);
                Is[i] = (4.0 / (N * N) * dv + 1.0) / degree;
            }
            mol.getAtomWithIdx(i)->setProp("IState", Is[i]);

        }
        return Is;
}


    // Calculate EState indices for a molecule
    std::vector<double> calcEStateIndices(const RDKit::ROMol& mol) {
       std::vector<double> Is = calcIStateIndices(mol);
        int numAtoms = mol.getNumAtoms();

        // Compute distance matrix
        double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"


        std::vector<double> accum(numAtoms, 0.0);

        // Compute accumulative EState contributions
        for (int i = 0; i < numAtoms; ++i) {
            for (int j = i + 1; j < numAtoms; ++j) {
                double p = distances[i*numAtoms +j] + 1.;  // Fix : use "p = distance + 1."
                if (p < 1e6) {  // Valid distance
                    double tmp = (Is[i] - Is[j]) / (p * p);
                    accum[i] += tmp;
                    accum[j] -= tmp;
                }
            }
        }

        // Combine initial EState values and accumulative contributions
        std::vector<double> res(numAtoms, 0.0);
        for (int i = 0; i < numAtoms; ++i) {
            res[i] = accum[i] + Is[i];
            mol.getAtomWithIdx(i)->setProp("EState", res[i]);

        }
        return res;
    }



    std::vector<double> CalcHEStateIndices(const RDKit::ROMol& mol) {
        size_t nAtoms = mol.getNumAtoms();


        int numAtoms = mol.getNumAtoms();

        std::vector<double> RKHE(numAtoms, 0.0);
        std::vector<double> hasHs(numAtoms, 0.0);

        for (int i = 0; i < numAtoms; ++i) {
            const auto* atom = mol.getAtomWithIdx(i);
            int hi = atom->getTotalNumHs();
            if (hi >0) {
                hasHs[i]= 1.;
                RKHE[i] = getRKHE(*atom);
            }

        }

        // Get distance matrix
       double* distMat = MolOps::getDistanceMat(mol, false, false, false); // no bond order, no weights, no hydrogens
;
        std::vector<std::vector<double>> invDists(nAtoms, std::vector<double>(nAtoms, 0.0));

        for (size_t i = 0; i < nAtoms; ++i) {
            for (size_t j = 0; j < nAtoms; ++j) {
                double dist = distMat[i*nAtoms+j]+1.;
                invDists[i][j] = dist > 0 ? 1.0 / (dist * dist) : 0.0;
            }
        }
        // Compute adjusted HEState
        std::vector<double> heStateIndices(nAtoms, 0.0);
        for (size_t i = 0; i < nAtoms; ++i) {
            double accum = 0.0;
            for (size_t j = 0; j < nAtoms; ++j) {
                if (i != j) {
                    accum += RKHE[j] * invDists[i][j];
                }
            }
            heStateIndices[i] = hasHs[i] * (accum + 2.0 * RKHE[i] + 0.2);
        }


        return heStateIndices;
    }




    // Helper function to calculate VSA_EState
    // Helper function for LabuteASA calculation
    std::vector<double> VSAcontrib(const RDKit::ROMol& mol, bool includeHs = true) {

        // Compute the Labute ASA contributions
        std::vector<double> atomContribs(mol.getNumAtoms(), 0.0);
        double totalASA = 0.0;
        RDKit::Descriptors::getLabuteAtomContribs(mol, atomContribs, totalASA, includeHs);


        // Include the total hydrogen contribution as the first element
        std::vector<double> Vi(atomContribs.size() + 1, 0.0);

        Vi[0] = totalASA;

        std::copy(atomContribs.begin(), atomContribs.end(), Vi.begin() + 1);

        // Cache the result as a property in the molecule
        // mol.setProp("_LabuteContribs", Vi);

        return Vi;
    }


    // Function to calculate EState-based VSA contributions
    std::vector<double> VSA_EState(const RDKit::ROMol& mol, const std::vector<double>& bins) {
        // Calculate EState indices
        std::vector<double> propContribs = calcEStateIndices(mol);

        // Calculate VSA contributions
        std::vector<double> volContribs = VSAcontrib(mol, true);


        // Initialize result vector
        std::vector<double> ans(bins.size() + 1, 0.0);

        // Assign contributions to bins
        for (size_t i = 0; i < propContribs.size(); ++i) {
            if (!std::isnan(propContribs[i])) {
                double volume = volContribs[i + 1];
                auto nbin = std::upper_bound(bins.begin(), bins.end(), volume) - bins.begin();
                ans[nbin] += propContribs[i];
            }
        }
        // Cache the result in the molecule
        // mol.setProp("_VSA_EState", ans);

        return ans;
    }


    // Function to calculate EState-based VSA contributions
    std::vector<double> EState_VSA(const RDKit::ROMol& mol, const std::vector<double>& bins) {
        // Calculate EState indices
        std::vector<double> propContribs = calcEStateIndices(mol);

        // Calculate VSA contributions
        std::vector<double> volContribs = VSAcontrib(mol, true);


        // Initialize result vector
        std::vector<double> ans(bins.size() + 1, 0.0);

        // Assign contributions to bins
        for (size_t i = 0; i < propContribs.size(); ++i) {
            if (!std::isnan(propContribs[i])) {
                auto nbin = std::upper_bound(bins.begin(), bins.end(), propContribs[i]) - bins.begin();
                ans[nbin] += volContribs[i + 1];
            }
        }
        // Cache the result in the molecule
        // mol.setProp("_EState_VSA", ans);

        return ans;
    }



    // this is fixed by p = distance + 1
    std::vector<double> calcVSA_EState(const RDKit::ROMol& mol) {
        const std::string customPropName = "EState";  // Assume EState values are stored in this property
        std::vector<double> vsaBins = {4.78, 5.00, 5.410, 5.740, 6.00, 6.07, 6.45, 7.00, 11.0};
        std::vector<double> vsaEstate  = VSA_EState( mol, vsaBins);
        return vsaEstate;
    }

    //
    std::vector<double> calcEState_VSA(const RDKit::ROMol& mol) {
        const std::string customPropName = "EState";  // Assume EState values are stored in this property
        std::vector<double> estateBins = {-0.390, 0.290, 0.717, 1.165, 1.540, 1.807, 2.05, 4.69, 9.17, 15.0};

        // Calculate the EState contributions for each atom
        std::vector<double> estateIndices = calcEStateIndices(mol);
        std::vector<double> eStateVSA  = EState_VSA( mol, estateBins);

        // Use calcCustomProp_VSA with EState bins
        //std::vector<double> eStateVSA = RDKit::Descriptors::calcCustomProp_VSA(
        //    mol, customPropName, estateBins);


        return eStateVSA;
    }

    std::vector<double> calcMoeType(const ROMol& mol) {
        // osmordredv3: 58 = LabuteASA + PEOE_VSA(14)+SMR_VSA(10)+SlogP_VSA(12)
        // +EState_VSA(11)+VSA_EState(10). Was res(54) with `p += size-1`, which
        // overlapped each family's first bin onto the previous family's last
        // bin (dropped 4 values). Now `p += size` (below) -> no overlap.
        std::vector<double> res(58, 0.0);
        double LabuteASA = RDKit::Descriptors::calcLabuteASA(mol);
        res[0] = LabuteASA;
        int p = 1;

        std::vector<double> PEOE_VSA = RDKit::Descriptors::calcPEOE_VSA(mol);
        //std::cout << PEOE_VSA.size() << " ";
        for (size_t i = 0; i < PEOE_VSA.size(); ++i) {
            res[p + i] = PEOE_VSA[i];  // Copy PEOE_VSA to res[1:13]
        }
        p += PEOE_VSA.size();
        //std::cout << "- new p: " << p << " | ";

        std::vector<double> SMR_VSA = RDKit::Descriptors::calcSMR_VSA(mol);
        //std::cout << SMR_VSA.size() << " ";

        for (size_t i = 0; i < SMR_VSA.size(); ++i) {
            res[p + i] = SMR_VSA[i];  // Copy SMR_VSA to res[14:23]
        }
        p += SMR_VSA.size();
        //std::cout << "- new p: " << p << " | ";

        std::vector<double> SlogP_VSA = RDKit::Descriptors::calcSlogP_VSA(mol);
        //std::cout << SlogP_VSA.size() << " ";
        for (size_t i = 0; i < SlogP_VSA.size(); ++i) {
            res[p + i] = SlogP_VSA[i];  // Copy SlogP_VSA to res[23:34]
        }
        p += SlogP_VSA.size();
        //std::cout << "- new p: " <<  p << " | ";


        std::vector<double> EState_VSA = calcEState_VSA(mol);
        //std::cout << EState_VSA.size() << " ";
        for (size_t i = 0; i < EState_VSA.size(); ++i) {
            res[p + i] = EState_VSA[i];  // Copy EState_VSA
        }
        p += EState_VSA.size();
        //std::cout << "- new p: " << p << " | ";

        // EState (mordred) VSA_EState 1-9 & EState_VSA 1-10
        std::vector<double> VSA_EState = calcVSA_EState(mol);
        //std::cout << VSA_EState.size() << "\n";
        for (size_t i = 0; i < VSA_EState.size(); ++i) {
            res[p + i] = VSA_EState[i];  // Copy SlogP_VSA to res[24:34]
	   // std::cout << " ( " << p+i;
	}
	//std::cout << ")\n";
        p += VSA_EState.size();
        //std::cout << "- new p: " << p << " end ";
        return res;
    }
/// logS


    // SMARTS patterns with their corresponding log contributions
    static const std::vector<std::pair<std::string, double>> smartsLogs = {
        {"[NH0;X3;v3]", 0.71535},
        {"[NH2;X3;v3]", 0.41056},
        {"[nH0;X3]", 0.82535},
        {"[OH0;X2;v2]", 0.31464},
        {"[OH0;X1;v2]", 0.14787},
        {"[OH1;X2;v2]", 0.62998},
        {"[CH2;!R]", -0.35634},
        {"[CH3;!R]", -0.33888},
        {"[CH0;R]", -0.21912},
        {"[CH2;R]", -0.23057},
        {"[ch0]", -0.37570},
        {"[ch1]", -0.22435},
        {"F", -0.21728},
        {"Cl", -0.49721},
        {"Br", -0.57982},
        {"I", -0.51547},
    };

    // Precompile SMARTS patterns to avoid repeated parsing
    static const std::vector<std::pair<std::shared_ptr<RDKit::RWMol>, double>> &GetSmartsLogs() {
      static const std::vector<std::pair<std::shared_ptr<RDKit::RWMol>, double>>
          compiledSmartsLogs = [] {
            std::vector<std::pair<std::shared_ptr<RDKit::RWMol>, double>> res;
            for (const auto &pair : smartsLogs) {
              auto mol = RDKit::SmartsToMol(pair.first);
              if (mol) {
                res.emplace_back(std::shared_ptr<RDKit::RWMol>(mol),
                                 pair.second);
              } else {
                std::cerr << "Invalid SMARTS: " << pair.first << std::endl;
              }
            }
            return res;
          }();
      return compiledSmartsLogs;
    }
    // Function to calculate LogS descriptor
    double LogS(const RDKit::ROMol& mol) {
        // Base formula contribution
        double molWeight = RDKit::Descriptors::calcAMW(mol, false); // Get molecular weight
        double logS = 0.89823 - 0.10369 * std::sqrt(molWeight);

        // Add contributions from precompiled SMARTS patterns
        for (const auto& pair : GetSmartsLogs()) {
            auto& smartsMol = pair.first;
            double logContribution = pair.second;

            if (!smartsMol) {
                continue;  // Skip invalid SMARTS
            }

            // Match SMARTS pattern
            std::vector<RDKit::MatchVectType> matches;
            RDKit::SubstructMatch(mol, *smartsMol, matches);

            // Add contributions for each match
            logS += matches.size() * logContribution;
        }

        return logS;
    }


    std::vector<double> calcLogS(const RDKit::ROMol& mol) {
        std::vector<double> res(1,0.);
        res[0] = LogS(mol);
        return res;
    }

    // Function to calculate Lipinski rule of 5
    int calculateLipinski(double LogP, double MW, double HBDon, double HBAcc) {
      bool L = (HBDon <= 5 && HBAcc <= 10 && MW <= 500 && LogP <= 5);

      if (L) {return 1;}
      else {return 0;}
    }

    // Function to calculate Ghose filter
    int calculateGhoseFilter(double MW, double LogP, double MR, int numAtoms) {
        bool G = (MW >= 160 && MW <= 480) &&
            (numAtoms >= 20 && numAtoms <= 70) &&
            (LogP >= -0.4 && LogP <= 5.6) &&
            (MR >= 40 && MR <= 130);

        if (G) {return 1;}
        else {return 0;}
    }

    // Main function to calculate Lipinski and Ghose filter
    std::vector<int> calcLipinskiGhose(const RDKit::ROMol& mol) {
        double MW =  RDKit::Descriptors::calcExactMW(mol);
        double HBDon = static_cast<double>(RDKit::Descriptors::calcNumHBD(mol));
        double HBAcc = static_cast<double>(RDKit::Descriptors::calcNumHBA(mol));
        double LogP;
        double MR;
        RDKit::Descriptors::calcCrippenDescriptors(mol, LogP, MR);

        int lipinski =  calculateLipinski(LogP, MW, HBDon, HBAcc);
        // must add Hs for Ghose

        std::unique_ptr<RDKit::ROMol> hmol(RDKit::MolOps::addHs(mol));


        int numAtoms = hmol->getNumAtoms();

        int ghoseFilter =  calculateGhoseFilter(MW, LogP, MR, numAtoms);

        return {lipinski, ghoseFilter};
    }

    static const std::map<int, double> VdWAtomicMap() {
        std::map<int, double> VdW_atomicRadii = {
            {1, 1.10}, {2, 1.40},
            {3, 1.82}, {4, 1.53}, {5, 1.92}, {6, 1.70}, {7, 1.55}, {8, 1.52}, {9, 1.47}, {10, 1.54},
            {11, 2.27}, {12, 1.73}, {13, 1.84}, {14, 2.10}, {15, 1.80}, {16, 1.80}, {17, 1.75}, {18, 1.88},
            {19, 2.75}, {20, 2.31}, {21, 2.15}, {22, 2.11}, {23, 2.07}, {24, 2.06}, {25, 2.05}, {26, 2.04}, {27, 2.00},
            {28, 1.97}, {29, 1.96}, {30, 2.01}, {31, 1.87}, {32, 2.11}, {33, 1.85}, {34, 1.90}, {35, 1.85}, {36, 2.02},
            {37, 3.03}, {38, 2.49}, {39, 2.32}, {40, 2.23}, {41, 2.18}, {42, 2.17}, {43, 2.16}, {44, 2.13}, {45, 2.10},
            {46, 2.10}, {47, 2.11}, {48, 2.18}, {49, 1.93}, {50, 2.17}, {51, 2.06}, {52, 2.06}, {53, 1.98}, {54, 2.16},
            {55, 3.43}, {56, 2.68},
            {57, 2.43}, {58, 2.42}, {59, 2.40}, {60, 2.39}, {61, 2.38}, {62, 2.36}, {63, 2.35}, {64, 2.34},
            {65, 2.33}, {66, 2.31}, {67, 2.30}, {68, 2.29}, {69, 2.27}, {70, 2.26}, {71, 2.24},
            {72, 2.23}, {73, 2.22}, {74, 2.18}, {75, 2.16}, {76, 2.16}, {77, 2.13}, {78, 2.13}, {79, 2.14},
            {80, 2.23}, {81, 1.96}, {82, 2.02}, {83, 2.07}, {84, 1.97}, {85, 2.02}, {86, 2.20},
            {87, 3.48}, {88, 2.83},
            {89, 2.47}, {90, 2.45}, {91, 2.43}, {92, 2.41}, {93, 2.39}, {94, 2.43}, {95, 2.44}, {96, 2.45},
            {97, 2.44}, {98, 2.45}, {99, 2.45}, {100, 2.45}, {101, 2.46}, {102, 2.46}, {103, 2.46}
        };

        return VdW_atomicRadii;
    }

    static const std::map<int, double> SandersonENAtomicMap() {
        std::map<int, double> SandersonElectronnegativityAtomicMap = {
            {1, 2.592}, {3, 0.670}, {4, 1.810}, {5, 2.275}, {6, 2.746}, {7, 3.194}, {8, 3.654}, {9, 4.000}, {10, 4.5},
            {11, 0.560}, {12, 1.318}, {13, 1.714}, {14, 2.138}, {15, 2.515}, {16, 2.957}, {17, 3.475}, {18, 3.31},
            {19, 0.445}, {20, 0.946}, {21, 1.02}, {22, 1.09}, {23, 1.39}, {24, 1.66}, {25, 2.2}, {26, 2.2}, {27, 2.56},
            {28, 1.94}, {29, 2.033}, {30, 2.223}, {31, 2.419}, {32, 2.618}, {33, 2.816}, {34, 3.014}, {35, 3.219},
            {36, 2.91}, {37, 0.312}, {38, 0.721}, {39, 0.65}, {40, 0.9}, {41, 1.42}, {42, 1.15}, {47, 1.826},
            {48, 1.978}, {49, 2.138}, {50, 2.298}, {51, 2.458}, {52, 2.618}, {53, 2.778}, {54, 2.34}, {55, 0.220},
            {56, 0.651}, {74, 0.98}, {80, 2.195}, {81, 2.246}, {82, 2.291}, {83, 2.342}
        };

      return SandersonElectronnegativityAtomicMap;
    }

    static const std::map<int, double> PaulingENAtomicMap() {
        std::map<int, double> PaulingelectronegativityAtomicMap = {
            {1, 2.2},  {3, 0.98}, {4, 1.57}, {5, 2.04},
            {6, 2.55}, {7, 3.04}, {8, 3.44}, {9, 3.98},
            {11, 0.93}, {12, 1.31}, {13, 1.61}, {14, 1.9}, {15, 2.19}, {16, 2.58}, {17, 3.16},
            {19, 0.82}, {20, 1.0}, {21, 1.36}, {22, 1.54},
            {23, 1.63}, {24, 1.66}, {25, 1.55}, {26, 1.83}, {27, 1.88}, {28, 1.91}, {29, 1.9}, {30, 1.65},
            {31, 1.81}, {32, 2.01}, {33, 2.18}, {34, 2.55}, {35, 2.96}, {36, 3.0}, {37, 0.82}, {38, 0.95},
            {39, 1.22}, {40, 1.33}, {41, 1.6}, {42, 2.16}, {43, 1.9}, {44, 2.2}, {45, 2.28}, {46, 2.2},
            {47, 1.93}, {48, 1.69}, {49, 1.78}, {50, 1.96}, {51, 2.05}, {52, 2.1}, {53, 2.66}, {54, 2.6},
            {55, 0.79}, {56, 0.89}, {57, 1.1}, {58, 1.12}, {59, 1.13}, {60, 1.14},
            {62, 1.17}, {64, 1.2},
            {66, 1.22}, {67, 1.23}, {68, 1.24}, {69, 1.25},
            {71, 1.27}, {72, 1.3}, {73, 1.5}, {74, 2.36}, {75, 1.9}, {76, 2.2}, {77, 2.2}, {78, 2.28},
            {79, 2.54}, {80, 2.0}, {81, 1.62}, {82, 2.33}, {83, 2.02}, {84, 2.0}, {85, 2.2},
            {87, 0.7}, {88, 0.9}, {89, 1.1}, {90, 1.3}, {91, 1.5}, {92, 1.38}, {93, 1.36}, {94, 1.28},
            {95, 1.3}, {96, 1.3}, {97, 1.3}, {98, 1.3}, {99, 1.3}, {100, 1.3}, {101, 1.3}, {102, 1.3}
        };

    return PaulingelectronegativityAtomicMap;
    }

    static const std::map<int, double> Allred_rocow_ENAtomicMap() {
        std::map<int, double> allred_rocow_electron_negativityAtomicMap = {
            {1, 2.20}, {3, 0.97}, {4, 1.47}, {5, 2.01}, {6, 2.50}, {7, 3.07}, {8, 3.50}, {9, 4.10},
            {11, 1.01}, {12, 1.23}, {13, 1.47}, {14, 1.74}, {15, 2.06}, {16, 2.44}, {17, 2.83},
            {19, 0.91}, {20, 1.04}, {21, 1.20}, {22, 1.32}, {23, 1.45}, {24, 1.56}, {25, 1.60},
            {26, 1.64}, {27, 1.70}, {28, 1.75}, {29, 1.75}, {30, 1.66}, {31, 1.82}, {32, 2.02},
            {33, 2.20}, {34, 2.48}, {35, 2.74}, {37, 0.89}, {38, 0.99}, {39, 1.11}, {40, 1.22},
            {41, 1.23}, {42, 1.30}, {43, 1.36}, {44, 1.42}, {45, 1.45}, {46, 1.35}, {47, 1.42},
            {48, 1.46}, {49, 1.49}, {50, 1.72}, {51, 1.82}, {52, 2.01}, {53, 2.21}, {55, 0.86},
            {56, 0.97}, {57, 1.08}, {72, 1.23}, {73, 1.33}, {74, 1.40}, {75, 1.46}, {76, 1.52},
            {77, 1.55}, {78, 1.44}, {79, 1.42}, {80, 1.44}, {81, 1.44}, {82, 1.55}, {83, 1.67},
            {84, 1.76}, {85, 1.90}
        };

    return allred_rocow_electron_negativityAtomicMap;
    }

    static const std::map<int, double> ionizationEnergyAtomicMap() {
        std::map<int, double> ionizationEnergyAtomicMap = {
            {1, 13.598443}, {2, 24.587387}, {3, 5.391719}, {4, 9.32270}, {5, 8.29802}, {6, 11.26030}, {7, 14.5341},
            {8, 13.61805}, {9, 17.4228}, {10, 21.56454}, {11, 5.139076}, {12, 7.646235}, {13, 5.985768}, {14, 8.15168},
            {15, 10.48669}, {16, 10.36001}, {17, 12.96763}, {18, 15.759610}, {19, 4.3406633}, {20, 6.11316},
            {21, 6.56149}, {22, 6.82812}, {23, 6.74619}, {24, 6.76651}, {25, 7.43402}, {26, 7.9024}, {27, 7.88101},
            {28, 7.6398}, {29, 7.72638}, {30, 9.394199}, {31, 5.999301}, {32, 7.89943}, {33, 9.7886}, {34, 9.75239},
            {35, 11.8138}, {36, 13.99961}, {37, 4.177128}, {38, 5.69485}, {39, 6.2173}, {40, 6.63390}, {41, 6.75885},
            {42, 7.09243}, {43, 7.28}, {44, 7.36050}, {45, 7.45890}, {46, 8.3369}, {47, 7.57623}, {48, 8.99382},
            {49, 5.78636}, {50, 7.34392}, {51, 8.60839}, {52, 9.0096}, {53, 10.45126}, {54, 12.12984}, {55, 3.893905},
            {56, 5.211664}, {57, 5.5769}, {58, 5.5387}, {59, 5.473}, {60, 5.5250}, {61, 5.582}, {62, 5.6437},
            {63, 5.67038}, {64, 6.14980}, {65, 5.8638}, {66, 5.9389}, {67, 6.0215}, {68, 6.1077}, {69, 6.18431},
            {70, 6.25416}, {71, 5.42586}, {72, 6.82507}, {73, 7.54957}, {74, 7.86403}, {75, 7.83352}, {76, 8.43823},
            {77, 8.96702}, {78, 8.9588}, {79, 9.22553}, {80, 10.4375}, {81, 6.108194}, {82, 7.41663}, {83, 7.2855},
            {84, 8.414}, {86, 10.7485}, {87, 4.072741}, {88, 5.278423}, {89, 5.17}, {90, 6.3067}, {91, 5.89},
            {92, 6.1941}, {93, 6.2657}, {94, 6.0260}, {95, 5.9738}, {96, 5.9914}, {97, 6.1979}, {98, 6.2817},
            {99, 6.42}, {100, 6.50}, {101, 6.58}, {102, 6.65}, {103, 4.9}, {104, 6.0}
            };

    return ionizationEnergyAtomicMap;
    }


    static const std::map<int, double> McGowanVolumAtomicMap() {
        std::map<int, double> atomicProperties = {
            {1, 8.71}, {2, 6.75}, {3, 22.23}, {4, 20.27}, {5, 18.31},
            {6, 16.35}, {7, 14.39}, {8, 12.43}, {9, 10.47}, {10, 8.51},
            {11, 32.71}, {12, 30.75}, {13, 28.79}, {14, 26.83}, {15, 24.87},
            {16, 22.91}, {17, 20.95}, {18, 18.99}, {19, 51.89}, {20, 50.28},
            {21, 48.68}, {22, 47.07}, {23, 45.47}, {24, 43.86}, {25, 42.26},
            {26, 40.65}, {27, 39.05}, {28, 37.44}, {29, 35.84}, {30, 34.23},
            {31, 32.63}, {32, 31.02}, {33, 29.42}, {34, 27.81}, {35, 26.21},
            {36, 24.60}, {37, 60.22}, {38, 58.61}, {39, 57.01}, {40, 55.40},
            {41, 53.80}, {42, 52.19}, {43, 50.59}, {44, 48.98}, {45, 47.38},
            {46, 45.77}, {47, 44.17}, {48, 42.56}, {49, 40.96}, {50, 39.35},
            {51, 37.75}, {52, 36.14}, {53, 34.54}, {54, 32.93}, {55, 77.25},
            {56, 76.00}, {57, 74.75}, {58, 73.49}, {59, 72.24}, {60, 70.99},
            {61, 69.74}, {62, 68.49}, {63, 67.23}, {64, 65.98}, {65, 64.73},
            {66, 63.48}, {67, 62.23}, {68, 60.97}, {69, 59.72}, {70, 58.47},
            {71, 57.22}, {72, 55.97}, {73, 54.71}, {74, 53.46}, {75, 52.21},
            {76, 50.96}, {77, 49.71}, {78, 48.45}, {79, 47.20}, {80, 45.95},
            {81, 44.70}, {82, 43.45}, {83, 42.19}, {84, 40.94}, {85, 39.69},
            {86, 38.44}, {87, 75.59}, {88, 74.34}, {89, 73.09}, {90, 71.83},
            {91, 70.58}, {92, 69.33}, {93, 68.08}, {94, 66.83}, {95, 65.57},
            {96, 64.32}, {97, 63.07}, {98, 61.82}, {99, 60.57}, {100, 59.31},
            {101, 58.06}, {102, 56.81}, {103, 55.56}
        };

        return atomicProperties;
    }



    double McGowanVolume(const RDKit::ROMol &mol) {
        // In Padel code this is /100 in order to match the Polarisability equation
        double res = 0.;
        std::unique_ptr<RDKit::ROMol> hmol(RDKit::MolOps::addHs(mol));

        std::map<int, double> mgvmap = McGowanVolumAtomicMap();

        for (const auto& atom : hmol->atoms()) {
            int atomicNum = atom->getAtomicNum();
               if (mgvmap.find(atomicNum) != mgvmap.end()) {
                res += mgvmap.at(atomicNum);
            }
        }
        double finalres = res - hmol->getNumBonds() * 6.56;

        return  finalres;
    }

  std::vector<double> calcMcGowanVolume(const RDKit::ROMol &mol) {
        std::vector<double> res(1,0.);
        res[0] = McGowanVolume(mol);
        return res;
 }


    // SMARTS patterns for fragments
    static const std::vector<std::string> PolFrags = {
        "[CX4H3]", "[CX4H2]", "[CX4H1]", "F", "Cl", "Br", "I",
        "[$([NX3](=O)=O),$([NX3+](=O)[O-])]", "[N,n]", "[O,o]",
        "[$([SX4](=[OX1])(=[OX1])([#6])[#6]),$([SX4+2]([OX1-])([OX1-])([#6])[#6])]", "[S,s]", "[P,p]"
    };

    // Corresponding coefficients for the fragments
    static const std::vector<double> coefPol = {
        10.152, 8.765, 5.702, 3.833, 16.557, 24.123, 38.506,
        10.488, 6.335, 4.307, 15.726, 22.366, 11.173
    };

    // Precompile SMARTS patterns for efficiency

    static const std::vector<std::shared_ptr<RDKit::RWMol>>& GetCompiledPolFrags() {
      static const std::vector<std::shared_ptr<RDKit::RWMol>> compiledPolFrags =
          [] {
            std::vector<std::shared_ptr<RDKit::RWMol>> res;
            for (const auto &smarts : PolFrags) {
              auto mol = RDKit::SmartsToMol(smarts);
              if (mol) {
                res.emplace_back(std::shared_ptr<RDKit::RWMol>(mol));
              } else {
                std::cerr << "Invalid SMARTS: " << smarts << std::endl;
                res.emplace_back(nullptr);  // Placeholder for invalid SMARTS
              }
            }
            return res;
          }();

      return compiledPolFrags;
    }

    // Function to count hydrogen atoms in the molecule
    int getNumHs(const RDKit::ROMol &mol) {
        int nHs = 0;
        for (const auto atom : mol.atoms()) {
            nHs += atom->getTotalNumHs();
        }
        return nHs;
    }

    // Function to calculate polarity descriptor
    double Polarity(const RDKit::ROMol &mol) {
        double res = -1.529;  // Intercept value

        // Add contributions from precompiled SMARTS patterns
        for (size_t i = 0; i < GetCompiledPolFrags().size(); ++i) {
            auto &pattern = GetCompiledPolFrags()[i];
            if (!pattern) continue;  // Skip invalid patterns

            std::vector<RDKit::MatchVectType> matches;
            RDKit::SubstructMatch(mol, *pattern, matches, true);  // uniquify = true
            res += matches.size() * coefPol[i];
        }

        // Add hydrogen contribution
        res += 3.391 * static_cast<double>(getNumHs(mol));

        return res;
    }


    std::vector<double> calcPol(const RDKit::ROMol &mol) {

        std::vector<double> res(1,0.);
        res[0] = Polarity(mol);
        return res;


    }

     double MRvalue(const RDKit::ROMol &mol) {

            return 4./3. * M_PI* Polarity(mol);
     }


    std::vector<double> calcMR(const RDKit::ROMol &mol) {
        std::vector<double> res(1,0.);
        res[0] =  MRvalue(mol);
        return res;

    }


    double ODT(const RDKit::ROMol &) {
        return 1;
    }


    std::vector<double> calcODT(const RDKit::ROMol &mol) {
        std::vector<double> res(1,0.);
        res[0] = ODT(mol);
        return res;


    }


// Function to calculate the Schultz descriptor
double Schultz(const ROMol &mol) {
    // Get the number of atoms in the molecule
    int nAtoms = mol.getNumAtoms();
    if (nAtoms == 0) return 0.0;

    // Get the distance matrix
    const double* distMat = MolOps::getDistanceMat(mol, false, false, false);

    // Get the adjacency matrix
    const double* adjMat = MolOps::getAdjacencyMatrix(mol, false, false, false, "noBO");
    // Calculate the vertex degree for each atom
    std::vector<double> vertexDegree(nAtoms, 0.0);
    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            vertexDegree[i] += adjMat[i * nAtoms + j];
        }
    }

    // Compute the Schultz descriptor
    double schultz = 0.0;
    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            double distance = distMat[i * nAtoms + j];
            double adjacency = adjMat[i * nAtoms + j];
            schultz += (distance + adjacency) * vertexDegree[j];
        }
    }

    return schultz;
}
std::vector<double> calcSchultz(const RDKit::ROMol &mol) {
    std::vector<double> res(1,0.);
    res[0] = Schultz(mol);
    return res;
}


    static const std::map<int, double> Polarizability78AtomicMap() {
        std::map<int, double> atomicProperties = {
            {1, 0.666793}, {2, 0.204956}, {3, 24.3}, {4, 5.6}, {5, 3.03},
            {6, 1.76}, {7, 1.1}, {8, 0.802}, {9, 0.557}, {10, 0.3956},
            {11, 23.6}, {12, 10.6}, {13, 6.8}, {14, 5.38}, {15, 3.63},
            {16, 2.9}, {17, 2.18}, {18, 1.6411}, {19, 43.4}, {20, 22.8},
            {21, 17.8}, {22, 14.6}, {23, 12.4}, {24, 11.6}, {25, 9.4},
            {26, 8.4}, {27, 7.5}, {28, 6.8}, {29, 6.1}, {30, 7.1},
            {31, 8.12}, {32, 6.07}, {33, 4.31}, {34, 3.77}, {35, 3.05},
            {36, 2.4844}, {37, 47.3}, {38, 27.6}, {39, 22.7}, {40, 17.9},
            {41, 15.7}, {42, 12.8}, {43, 11.4}, {44, 9.6}, {45, 8.6},
            {46, 4.8}, {47, 7.2}, {48, 7.2}, {49, 10.2}, {50, 7.7},
            {51, 6.6}, {52, 5.5}, {53, 5.35}, {54, 4.044}, {55, 59.6},
            {56, 39.7}, {57, 31.1}, {58, 29.6}, {59, 28.2}, {60, 31.4},
            {61, 30.1}, {62, 28.8}, {63, 27.7}, {64, 23.5}, {65, 25.5},
            {66, 24.5}, {67, 23.6}, {68, 22.7}, {69, 21.8}, {70, 21.0},
            {71, 21.9}, {72, 16.2}, {73, 13.1}, {74, 11.1}, {75, 9.7},
            {76, 8.5}, {77, 7.6}, {78, 6.5}, {79, 5.8}, {80, 5.7},
            {81, 7.6}, {82, 6.8}, {83, 7.4}, {84, 6.8}, {85, 6.0},
            {86, 5.3}, {87, 48.7}, {88, 38.3}, {89, 32.1}, {90, 32.1},
            {91, 25.4}, {92, 27.4}, {93, 24.8}, {94, 24.5}, {95, 23.3},
            {96, 23.0}, {97, 22.7}, {98, 20.5}, {99, 19.7}, {100, 23.8},
            {101, 18.2}, {102, 17.5}
        };

        return atomicProperties;
    }

    static const std::map<int, double> Polarizability94AtomicMap() {
        std::map<int, double> atomicProperties = {
            {1, 0.666793}, {2, 0.2050522}, {3, 24.33}, {4, 5.60}, {5, 3.03},
            {6, 1.67}, {7, 1.10}, {8, 0.802}, {9, 0.557}, {10, 0.39432},
            {11, 24.11}, {12, 10.6}, {13, 6.8}, {14, 5.53}, {15, 3.63},
            {16, 2.90}, {17, 2.18}, {18, 1.6411}, {19, 43.06}, {20, 22.8},
            {21, 17.8}, {22, 14.6}, {23, 12.4}, {24, 11.6}, {25, 9.4},
            {26, 8.4}, {27, 7.5}, {28, 6.8}, {29, 6.2}, {30, 5.75},
            {31, 8.12}, {32, 5.84}, {33, 4.31}, {34, 3.77}, {35, 3.05},
            {36, 2.4844}, {37, 47.24}, {38, 23.5}, {39, 22.7}, {40, 17.9},
            {41, 15.7}, {42, 12.8}, {43, 11.4}, {44, 9.6}, {45, 8.6},
            {46, 4.8}, {47, 6.78}, {48, 7.36}, {49, 10.2}, {50, 7.84},
            {51, 6.6}, {52, 5.5}, {53, 5.35}, {54, 4.044}, {55, 59.42},
            {56, 39.7}, {57, 31.1}, {58, 29.6}, {59, 28.2}, {60, 31.4},
            {61, 30.1}, {62, 28.8}, {63, 27.7}, {64, 23.5}, {65, 25.5},
            {66, 24.5}, {67, 23.6}, {68, 22.7}, {69, 21.8}, {70, 20.9},
            {71, 21.9}, {72, 16.2}, {73, 13.1}, {74, 11.1}, {75, 9.7},
            {76, 8.5}, {77, 7.6}, {78, 6.5}, {79, 5.8}, {80, 5.02},
            {81, 7.6}, {82, 7.01}, {83, 7.4}, {84, 6.8}, {85, 6.0},
            {86, 5.3}, {87, 48.6}, {88, 38.3}, {89, 32.1}, {90, 32.1},
            {91, 25.4}, {92, 24.9}, {93, 24.8}, {94, 24.5}, {95, 23.3},
            {96, 23.0}, {97, 22.7}, {98, 20.5}, {99, 19.7}, {100, 23.8},
            {101, 18.2}, {102, 16.4}, {112, 4.06}, {114, 4.59}, {119, 24.26}
        };

        return atomicProperties;
    }

    // Combined function for calculating both atomic and bond polarizability
    std::vector<double> calcPolarizability(const RDKit::ROMol &mol) {
        double atomicPol = 0.0;
        double bondPol = 0.0;
        const auto&  polmap = Polarizability94AtomicMap();
        std::unique_ptr<RDKit::ROMol> hmol(RDKit::MolOps::addHs(mol));



        for (const auto& atom : hmol->atoms()) {
            int atomicNum = atom->getAtomicNum();
            auto it = polmap.find(atomicNum);
            if (it != polmap.end()) {
                atomicPol += it->second;
            }
        }
        // Calculate bond polarizability
        for (const auto& bond : hmol->bonds()) {
            int atomicNum1 = bond->getBeginAtom()->getAtomicNum();
            int atomicNum2 = bond->getEndAtom()->getAtomicNum();

            auto it1 = polmap.find(atomicNum1);
            auto it2 = polmap.find(atomicNum2);

            if (it1 != polmap.end() && it2 != polmap.end()) {
                bondPol += std::abs(it1->second - it2->second);
            }
        }

        return {atomicPol,bondPol};
    }


    // Main Rotatabond
    std::vector<double> calcRotatableBond(const ROMol& mol) {
        double rot = RDKit::Descriptors::calcNumRotatableBonds(mol);

        std::vector<double> res(2, 0.0);
        res[0] = rot;

        int bondcountheavy = 0;
        for (const auto& bond : mol.bonds()) {
            if (bond->getBeginAtom()->getAtomicNum() != 1 && bond->getEndAtom()->getAtomicNum() != 1) {
                bondcountheavy +=1; // Heavy bond
            }
        }
        if (bondcountheavy>0) {
            res[1] = static_cast<double>(static_cast<double>(rot) / static_cast<double>(bondcountheavy));
        }
        return res;
    }

    double FragmentComplexity(const ROMol& mol) {
        // Number of atoms (A)
        int A = mol.getNumAtoms();

        // Number of bonds (B)
        int B = mol.getNumBonds();

        // Number of heteroatoms (H) - atoms that are not Carbon
        int H = 0;
        for (const auto& atom : mol.atoms()) {
            if (atom->getAtomicNum() != 6) {  // Exclude Carbon atoms
                H++;
            }
        }

        // Calculate fragment complexity: |B^2 - A^2 + A| + H / 100
        double fragCpx = std::abs(std::pow(B, 2) - std::pow(A, 2) + A) + H / 100.0;

        return fragCpx;
    }


    std::vector<double> calcFragmentComplexity(const ROMol& mol) {
	
        std::vector<double> res(1, 0.);
	    res[0]=FragmentComplexity(mol);
        return res;
    }

    // optimize construction of symetrics matrix only do j=i not j=0! Upper mat
     Eigen::MatrixXd calculateDistanceMatrix(const ROMol& mol) {
        unsigned int nAtoms = mol.getNumAtoms();

        // Get the distance matrix using RDKit's MolOps::getDistanceMat
        double* distanceMat = MolOps::getDistanceMat(mol, false, false, false); // no bond order, no weights, no hydrogens

        // Convert the raw pointer to an Eigen MatrixXd (distance matrix)
         Eigen::MatrixXd distMatrix(nAtoms, nAtoms);

        for (unsigned int i = 0; i < nAtoms; ++i) {
            for (unsigned int j = i; j < nAtoms; ++j) {
                distMatrix(i, j) = distanceMat[i * nAtoms + j];
                if (j>i) {
                    distMatrix(j,i) = distMatrix(i, j); // diagonal is already set only the symetrical part need to be clone
                }
            }
        }

        return distMatrix;
    }


    // Function to calculate eccentricity (max distance for each atom)
    Eigen::VectorXd calculateEccentricity(const ROMol& mol) {
        unsigned int nAtoms = mol.getNumAtoms();

        // Calculate the distance matrix using Eigen
        Eigen::MatrixXd distanceMatrix = calculateDistanceMatrix(mol);

        // For each atom, calculate the maximum distance to any other atom (eccentricity)
        Eigen::VectorXd eccentricity(nAtoms);

        // Eigen's rowwise max operation
        eccentricity = distanceMatrix.colwise().maxCoeff();

        return eccentricity;
    }

    // it is not correct!!!!
    std::vector<double> calcEccentricConnectivityIndex(const ROMol& mol) {
	std::vector<double> res(1,0.);
        Eigen::VectorXd E = calculateEccentricity(mol);
        std::vector<double> V = calcValence(mol);


        if ((size_t)E.size() != V.size()) {
            throw std::invalid_argument("Eccentricity and Valence vectors must have the same length.");
        }

        // Calculate element-wise product of E and V
        double productSum = 0.0;
        for (size_t i = 0; i < (size_t)E.size(); ++i) {
            productSum += static_cast<int>(E[i]) * V[i];  // Cast E[i] to int and multiply with V[i]
        }
	res[0] = productSum;
        // Return the result as an integer
        return res;//static_cast<int>(productSum);
    }


    // Function to find rings in a molecule using RDKit's ring detection
    std::vector<std::vector<int>> findRings(const ROMol& mol) {
        std::vector<std::vector<int>> rings;
        RingInfo* ri = mol.getRingInfo();
        if (ri) {
            for (size_t i = 0; i < ri->numRings(); ++i) {
                std::vector<int> ring = ri->bondRings()[i];
                rings.push_back(ring);
            }
        }
        return rings;
    }

    // Constitutional
    // code vs name
    //Z    a.GetAtomicNum()
    //m    mass[a.GetAtomicNum()]
    //v    vdw_volume[a.GetAtomicNum()]
    //se   sanderson[a.GetAtomicNum()]
    //pe   pauling[a.GetAtomicNum()]
    //are  allred_rocow[a.GetAtomicNum()]
    //p    polarizability94[a.GetAtomicNum()] (last as default!!!)
    //i    ionization_potentials[a.GetAtomicNum()]

    // convert radius to volume
    double vdw_volume (double r) {
     return (4.0 / 3.0) *  M_PI * std::pow(r, 3);
    }

    std::vector<double> calcConstitutional(const ROMol& mol) {
            double SZ = 0.;
            double Sm = 0.;
            double Sv = 0.;
            double Sse = 0.;
            double Spe = 0.;
            double Sare = 0.;
            double Sp = 0.;
            double Si = 0.;
            double MZ, Mm, Mv, Mse, Mpe, Mare, Mp, Mi;
            const PeriodicTable *tbl = PeriodicTable::getTable();
            std::unique_ptr<RDKit::ROMol> hmol(RDKit::MolOps::addHs(mol));

            double zcc = static_cast<double>(tbl->getAtomicNumber("C"));

            double mcc = static_cast<double>(tbl->getAtomicWeight("C"));

            std::map<int, double> vdwmap = VdWAtomicMap();
            double mvdwc = vdw_volume(vdwmap[6]);

            std::map<int, double> semap = SandersonENAtomicMap();
            double msec = semap[6];

            std::map<int, double> pemap = PaulingENAtomicMap();
            double mpec = pemap[6];

            std::map<int, double> aremap = Allred_rocow_ENAtomicMap();
            double marec = aremap[6];

            std::map<int, double> pmap = Polarizability94AtomicMap();
            double mpc = pmap[6];

            std::map<int, double> imap = ionizationEnergyAtomicMap();
            double mic = imap[6];

            for (const auto& atom : hmol->atoms()) {
                std::string symbol = atom->getSymbol();
                int atzi = atom->getAtomicNum();
                SZ += static_cast<double>(atzi) / zcc;
                Sm += static_cast<double>(tbl->getAtomicWeight(symbol)) / mcc;
                // need to convert radius to volume!!!
                if (vdwmap.find(atzi) != vdwmap.end()) {
                    Sv += vdw_volume(vdwmap.at(atzi)) / mvdwc;
                }
                if (semap.find(atzi) != semap.end()) {
                    Sse += semap.at(atzi) / msec;
                }
                if (pemap.find(atzi) != pemap.end()) {
                    Spe += pemap.at(atzi) / mpec;
                }
                if (aremap.find(atzi) != aremap.end()) {
                    Sare += aremap.at(atzi) / marec;
                }
                if (pmap.find(atzi) != pmap.end()) {
                    Sp += pmap.at(atzi) / mpc;
                }
                if (imap.find(atzi) != imap.end()) {
                    Si += imap.at(atzi) / mic;
                }

            }

            double natoms = static_cast<double>(hmol->getNumAtoms());
            MZ = SZ / natoms;  // correct
            Mm = Sm / natoms;  // correct
            Mv = Sv / natoms; // correct trick convert radius to volume!!!
            Mse = Sse / natoms; // correct
            Mpe = Spe / natoms; // correct
            Mare = Sare / natoms; // correct
            Mp = Sp / natoms; // correct
            Mi = Si / natoms; // correct

            return {SZ, Sm, Sv, Sse, Spe, Sare, Sp, Si, MZ, Mm, Mv, Mse, Mpe, Mare, Mp, Mi};
    }




////// Barysz Matrixes Eigen style


    // Function to compute eigenvalues and eigenvectors
    void compute_eigenvalues_and_eigenvectors(const Eigen::MatrixXd& matrix,
                                            Eigen::VectorXd& eigenvalues,
                                            Eigen::MatrixXd& eigenvectors) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix); // hermitian (because we have symetrical matrix! critial)
        eigenvalues = solver.eigenvalues().real();  // Take real parts if complex
        eigenvectors = solver.eigenvectors().real();

        // Canonicalize eigenvector signs: force first non-zero entry positive per column
        const double canonical_eps = 1e-12;
        for (int j = 0; j < eigenvectors.cols(); ++j) {
            for (int i = 0; i < eigenvectors.rows(); ++i) {
                double v = eigenvectors(i, j);
                if (std::abs(v) > canonical_eps) {
                    if (v < 0.0) {
                        eigenvectors.col(j) *= -1.0;
                    }
                    break;
                }
            }
        }
    }

    // Spectral Absolute Sum: sum of absolute eigenvalues
    double spAbs(const Eigen::VectorXd& eigenvalues) {
        return eigenvalues.array().abs().sum();
    }

    // Leading Eigenvalue: largest eigenvalue
    double spMax(const Eigen::VectorXd& eigenvalues) {
        return eigenvalues.maxCoeff();
    }

    // Spectral Diameter: difference between largest and smallest eigenvalues
    double spDiam(const Eigen::VectorXd& eigenvalues) {
        return eigenvalues.maxCoeff() - eigenvalues.minCoeff();
    }

    // Mean of Eigenvalues: average eigenvalue
    double spMean(const Eigen::VectorXd& eigenvalues) {
        return eigenvalues.mean();
    }

    // Spectral Absolute Deviation: sum of absolute deviations from the mean
    double spAD(const Eigen::VectorXd& eigenvalues, double mean) {
        return (eigenvalues.array() - mean).abs().sum();
    }

    // Estrada-like index: log sum of exponentials of eigenvalues
    double logEE(const Eigen::VectorXd& eigenvalues) {
        double a = std::max(eigenvalues.maxCoeff(), 0.0);
        double sx = (eigenvalues.array() - a).exp().sum() + std::exp(-a);
        return a + std::log(sx);
    }

    double SM1(const Eigen::MatrixXd& matrix) {
        return matrix.trace();
    }


    // Coefficient Sum of the Last Eigenvector (VE1)
    double VE1(const Eigen::MatrixXd& matrix, Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors) {
        compute_eigenvalues_and_eigenvectors(matrix, eigenvalues, eigenvectors);  // Reuse the eigenvalue computation
        Eigen::VectorXd eigenvector = eigenvectors.col(eigenvalues.size() - 1);  // Last eigenvector
        return eigenvector.array().abs().sum();  // Sum of the absolute values
    }

    // Average Coefficient of the Last Eigenvector (VE2)
    double VE2(const Eigen::MatrixXd& matrix, int numAtoms, Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors) {
        double VE1_value = VE1(matrix, eigenvalues, eigenvectors);
        return VE1_value / numAtoms;
    }

    // Logarithmic Coefficient Sum of the Last Eigenvector (VE3)
    double VE3(const Eigen::MatrixXd& matrix, int numAtoms, Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors) {
        double VE1_value = VE1(matrix, eigenvalues, eigenvectors);
        return std::log(0.1 * numAtoms * VE1_value);
    }


    // Randic-like Eigenvector-Based Index (VR1)
    double VR1(const Eigen::MatrixXd& matrix, const std::vector<std::pair<int, int>>& bonds,
            Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors) {
        compute_eigenvalues_and_eigenvectors(matrix, eigenvalues, eigenvectors);  // Reuse the eigenvalue computation
        Eigen::VectorXd eigenvector = eigenvectors.col(eigenvalues.size() - 1);  // Last eigenvector
        double result = 0.0;

        for (const auto& bond : bonds) {
            int i = bond.first;
            int j = bond.second;
            result += std::pow(std::abs(eigenvector(i)) * std::abs(eigenvector(j)), -0.5);
        }

        return result;
    }

    // Normalized Randic-like Eigenvector-Based Index (VR2)
    double VR2(const Eigen::MatrixXd& matrix, const std::vector<std::pair<int, int>>& bonds,
            int numAtoms, Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors) {
        double VR1_value = VR1(matrix, bonds, eigenvalues, eigenvectors);
        return VR1_value / numAtoms;
    }

    // Logarithmic Randic-like Eigenvector-Based Index (VR3)
    double VR3(const Eigen::MatrixXd& matrix, const std::vector<std::pair<int, int>>& bonds,
            int numAtoms, Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors) {
        double VR1_value = VR1(matrix, bonds, eigenvalues, eigenvectors);
        return std::log(0.1 * numAtoms * VR1_value);
    }


    // BaryszMatrix deps

    Eigen::MatrixXd floydWarshall(Eigen::MatrixXd& A) {
        int n = A.rows();

        // Set diagonal elements to 0
        A.diagonal().setZero();

        // Perform the relaxation to get shortest path BFS equivalent faster for complex case than BFS
        for (int i = 0; i < n; ++i) {
            // Broadcast the row i and column i sums and apply the min operation
            Eigen::RowVectorXd rowSum = A.row(i);  // i-th row
            Eigen::VectorXd colSum = A.col(i);    // i-th column

            // Ensure correct dimensions for broadcasting
            Eigen::MatrixXd rowpart = rowSum.replicate(n, 1);
            Eigen::MatrixXd colpart = colSum.replicate(1, n);

            Eigen::MatrixXd tempMatrix = rowpart + colpart;

            A = A.cwiseMin(tempMatrix);
        }

        return A;
    }



// option using inspired CDS paper https://www.sciencedirect.com/org/science/article/pii/S2635098X24001426
//  full step matrix implements

    std::vector<std::vector<float>>
    FastWeightedStepMatrix(int n_atom, std::vector<std::vector<float>> &Sa,
                        std::unordered_map<int, std::set<int>> &Ladj,
                        std::set<std::pair<int, int>> step_xy) {
        std::vector<std::vector<float>> SF = Sa;

        // Iterate step-wise from 1 to n_atom-2
        for (int m = 1; m < n_atom - 1; ++m) {
            if (step_xy.empty()) {
                break;
            }
            std::set<std::pair<int, int>> temp_step_xy{};
            for (const auto &item: step_xy) {
                int r = std::get<0>(item);
                for (int c: Ladj[std::get<1>(item)]) {
                    // Update path weight if a shorter weighted path is found
                    float newWeight = SF[r - 1][std::get<1>(item) - 1] + Sa[std::get<1>(item) - 1][c - 1];
                    if ((SF[r - 1][c - 1] == 0 || newWeight < SF[r - 1][c - 1]) && r != c) {
                        SF[r - 1][c - 1] = newWeight;
                        temp_step_xy.insert(std::make_pair(r, c));
                    }
                }
            }
            step_xy = temp_step_xy;
        }
        return SF;
    }



    std::vector<std::vector<float>> initializeWeightedAdjacency(
    const RDKit::ROMol &mol,
    const std::vector<double> &atomicProps) {
    int numAtoms = mol.getNumAtoms();
    std::vector<std::vector<float>> Sa(numAtoms, std::vector<float>(numAtoms, 0.0));

    for (const auto &bond : mol.bonds()) {
            unsigned int i = bond->getBeginAtomIdx();
            unsigned int j = bond->getEndAtomIdx();
            double pi = bond->getBondTypeAsDouble(); // Bond order
            Sa[i][j] = 1.0 / (atomicProps[i] * atomicProps[j] * pi); // Weighted edge
            Sa[j][i] = Sa[i][j];
        }

    return Sa;
    }


    std::tuple<std::vector<std::vector<float>>, std::unordered_map<int, std::set<int>>, std::set<std::pair<int, int>>>
    prepareAdjacencyData(const std::vector<std::vector<float>> &Sa) {
        int numAtoms = Sa.size();
        std::unordered_map<int, std::set<int>> Ladj;
        std::set<std::pair<int, int>> step_xy;

        for (int i = 0; i < numAtoms; ++i) {
            for (int j = 0; j < numAtoms; ++j) {
                if (Sa[i][j] > 0.0) { // Edge exists
                    Ladj[i + 1].insert(j + 1);
                    step_xy.insert({i + 1, j + 1});
                }
            }
        }

        return {Sa, Ladj, step_xy};
    }


    Eigen::MatrixXd computeBaryszMatrix0(
        const RDKit::ROMol& mol,
        const std::vector<double>& atomicProps) {

        const unsigned int numAtoms = mol.getNumAtoms();

        // Initialize Eigen matrix with a large value for non-edges
        double largeValue = 1e6;
        Eigen::MatrixXd baryszMatrix = Eigen::MatrixXd::Constant(numAtoms, numAtoms, largeValue);

        // Set diagonal elements to 0 initially already done in floydWarshallRelaxation
        for (unsigned int i = 0; i < numAtoms; ++i) {
            baryszMatrix(i, i) = 0.0;
        }


        std::vector<double> reciprocalAtomicProps(numAtoms);
        for (unsigned int i = 0; i < numAtoms; ++i) {
            reciprocalAtomicProps[i] = 1.0 / atomicProps[i];
        }

        // Fill the matrix with edge weights from bonds
        for (const auto& bond : mol.bonds()) {
            unsigned int i = bond->getBeginAtomIdx();
            unsigned int j = bond->getEndAtomIdx();

            double pi = bond->getBondTypeAsDouble();  // Bond order it is already normalized
            //double w = (cw * cw) / (atomicProps[i] * atomicProps[j] * pi); // already normalyzed so 1 instead of cw*cw
            double w = reciprocalAtomicProps[i] * reciprocalAtomicProps[j] / pi;
            baryszMatrix(i, j) = w;
            baryszMatrix(j, i) = w;
        }

        // Floyd-Warshall algorithm for shortest paths using Eigen
        Eigen::MatrixXd  baryszMatrix_ = floydWarshall(baryszMatrix);
        // can we use this instead for speed https://github.com/FangyouYan/Connectivity-Stepwise-Derivation-CSD-/blob/master/cpp_code/common_utils/MsUtils.cpp

        // Update diagonal with 1.0 - C / P[i] this is 1 / normP[i] in our case
        for (unsigned int i = 0; i < numAtoms; ++i) {
            baryszMatrix_(i, i) = 1.0 - 1 / atomicProps[i];
        }

        // Replace any remaining large values with 0 to ensure no invalid entries
        baryszMatrix_ = (baryszMatrix_.array() == largeValue).select(0.0, baryszMatrix_);


        return baryszMatrix_;
    }



Eigen::MatrixXd computeBaryszMatrix1(
    const RDKit::ROMol &mol,
    const std::vector<double> &atomicProps) {

    const unsigned int numAtoms = mol.getNumAtoms();

    auto Sa = initializeWeightedAdjacency(mol, atomicProps);

    // Prepare Ladj and step_xy structures for FastFullStepMatrix but this is the same for all descriptors or ???

    auto [SaMatrix, Ladj, step_xy] = prepareAdjacencyData(Sa);

    // Compute shortest path matrix using FastFullStepMatrix
    auto shortestPathMatrix = FastWeightedStepMatrix(mol.getNumAtoms(), SaMatrix, Ladj, step_xy);

    // Convert shortestPathMatrix to Eigen matrix
    Eigen::MatrixXd baryszMatrix = Eigen::MatrixXd::Zero(numAtoms, numAtoms);
        for (unsigned int i = 0; i < numAtoms; ++i) {
            for (unsigned int j = i+1; j < numAtoms; ++j) {
                baryszMatrix(i, j) = shortestPathMatrix[i][j];
                baryszMatrix(j, i) = baryszMatrix(i, j);
            }
        }
    // Update diagonal with 1.0 - (1 / P[i])
    for (unsigned int i = 0; i < numAtoms; ++i) {
        baryszMatrix(i, i) = 1.0 - (1.0 / atomicProps[i]);
    }

    return baryszMatrix;
}



// Function to calculate the Barysz matrix
Eigen::MatrixXd computeBaryszMatrix2(const RDKit::ROMol& mol, const std::vector<double>& w) {
    unsigned int natom = mol.getNumAtoms();
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(natom, natom);

    std::vector<double> reciprocalw(natom);
    for (unsigned int i = 0; i < natom; ++i) {
        reciprocalw[i] = 1.0 / w[i];
    }

    // Set diagonal entries
    for (unsigned int i = 0; i < natom; ++i) {
        matrix(i, i) = 1.0 - reciprocalw[i];
    }

    // Set off-diagonal entries
    for (unsigned int i = 0; i < natom; ++i) {
        for (unsigned int j = i + 1; j < natom; ++j) {

            std::list<int> path = RDKit::MolOps::getShortestPath(mol, i, j);

            for (auto it = path.begin(); std::next(it) != path.end(); ++it) {
                int atomIdx1 = *it;
                int atomIdx2 = *std::next(it);

                //const RDKit::Atom* atom1 = mol.getAtomWithIdx(atomIdx1);
                //const RDKit::Atom* atom2 = mol.getAtomWithIdx(atomIdx2);
                const RDKit::Bond* bond = mol.getBondBetweenAtoms(atomIdx1, atomIdx2);

                double weights = reciprocalw[atomIdx1] * reciprocalw[atomIdx2];
                if (bond) {
                    if (bond->getIsAromatic()) {
                        matrix(i, j) += weights / 1.5;
                    } else {
                        switch (bond->getBondType()) {
                            case RDKit::Bond::SINGLE:
                                matrix(i, j) += weights;
                                break;
                            case RDKit::Bond::DOUBLE:
                                matrix(i, j) += 0.5 * weights;
                                break;
                            case RDKit::Bond::TRIPLE:
                                matrix(i, j) += weights / 3.0;
                                break;
                            default:
                                break;
                        }
                    }
                }
            }
            matrix(j, i) = matrix(i, j);  // Ensure symmetry
        }
    }

    return matrix;
}


    std::vector<double> MatrixDescs(const RDKit::ROMol& mol, Eigen::MatrixXd Mat ) {

        int numAtoms = mol.getNumAtoms();

        Eigen::VectorXd eigenvalues;
        Eigen::MatrixXd eigenvectors;

        compute_eigenvalues_and_eigenvectors(Mat, eigenvalues, eigenvectors);


        std::vector<std::pair<int, int>> bonds;
            // Iterate over all bonds in the molecule
            for (const auto& bond : mol.bonds()) {
                int beginAtomIdx = bond->getBeginAtomIdx();
                int endAtomIdx = bond->getEndAtomIdx();
                bonds.push_back(std::make_pair(beginAtomIdx, endAtomIdx));
            }

        // Compute descriptors
        double Sp_Abs = spAbs(eigenvalues);
        double Sp_Max = spMax(eigenvalues);
        double Sp_Diam = spDiam(eigenvalues);
        double Sp_Mean = spMean(eigenvalues); // tmp values not needed to export as result
        double Sp_AD = spAD(eigenvalues, Sp_Mean);
        double Sp_MAD = Sp_AD / numAtoms;
        double Log_EE = logEE(eigenvalues);  // this one is not correct ??? do we need to add the bond
        // look correct per say! ve1... vr3

        double sm1 = SM1(Mat);
        double ve1 = VE1(Mat, eigenvalues, eigenvectors);
        double ve2 = VE2(Mat, numAtoms, eigenvalues, eigenvectors);
        double ve3 = VE3(Mat, numAtoms, eigenvalues, eigenvectors) ;
        double vr1 = VR1(Mat, bonds, eigenvalues, eigenvectors);
        double vr2 = VR2(Mat, bonds, numAtoms, eigenvalues, eigenvectors) ;
        double vr3 = VR3(Mat, bonds, numAtoms, eigenvalues, eigenvectors);

        return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE, sm1, ve1, ve2, ve3, vr1, vr2, vr3}; // 13 outputs

    }


    std::vector<double> calcBaryszMatrixDescs(const RDKit::ROMol& mol) {
        int method = 1;
        const PeriodicTable* tbl = PeriodicTable::getTable();

        std::map<int, double> vdwmap = VdWAtomicMap();
        std::map<int, double> semap = SandersonENAtomicMap();
        std::map<int, double> pemap = PaulingENAtomicMap();
        std::map<int, double> aremap = Allred_rocow_ENAtomicMap();
        std::map<int, double> pmap = Polarizability94AtomicMap();
        std::map<int, double> imap = ionizationEnergyAtomicMap();

        double zcc = static_cast<double>(tbl->getAtomicNumber("C"));
        double mcc = static_cast<double>(tbl->getAtomicWeight("C"));
        double mvdwc = vdw_volume(vdwmap[6]);
        double msec = semap[6];
        double mpec = pemap[6];
        double marec = aremap[6];
        double mpc = pmap[6];
        double mic = imap[6];

        // Prepare vectors to store the computed atomic properties
        std::vector<double> Zs_(mol.getNumAtoms(),0.);
        std::vector<double> ms_(mol.getNumAtoms(),0.);
        std::vector<double> vs_(mol.getNumAtoms(),0.);
        std::vector<double> pes_(mol.getNumAtoms(),0.);
        std::vector<double> ses_(mol.getNumAtoms(),0.);
        std::vector<double> ares_(mol.getNumAtoms(),0.);
        std::vector<double> ps_(mol.getNumAtoms(),0.);
        std::vector<double> is_(mol.getNumAtoms(),0.);

        for (const auto& atom : mol.atoms()) {
            int atzi = atom->getAtomicNum();
            // use autonormalized atomproperties so we simplify the code as it is already (P / cw not just P)
            double zValue = static_cast<double>(atzi) / zcc;
            double mValue = static_cast<double>(tbl->getAtomicWeight(atzi)) / mcc;
            double vValue = (vdwmap.find(atzi) != vdwmap.end()) ? vdw_volume(vdwmap.at(atzi)) / mvdwc : 0.0;
            double seValue = (semap.find(atzi) != semap.end()) ? semap.at(atzi) / msec : 0.0;
            double peValue = (pemap.find(atzi) != pemap.end()) ? pemap.at(atzi) / mpec : 0.0;
            double areValue = (aremap.find(atzi) != aremap.end()) ? aremap.at(atzi) / marec : 0.0;
            double pValue = (pmap.find(atzi) != pmap.end()) ? pmap.at(atzi) / mpc : 0.0;
            double iValue = (imap.find(atzi) != imap.end()) ? imap.at(atzi) / mic : 0.0;

            // You can use any of the computed values (zValue, mValue, etc.)
            // Here, choose the property you want to compute the Barysz matrix for:
            Zs_[atom->getIdx()] = zValue; //  zValue
            ms_[atom->getIdx()] = mValue; // mass
            vs_[atom->getIdx()] = vValue; // volume
            ses_[atom->getIdx()] = seValue; // sanderson
            pes_[atom->getIdx()] = peValue; // pauling
            ares_[atom->getIdx()] = areValue; // allred
            ps_[atom->getIdx()] = pValue; // polarisability
            is_[atom->getIdx()] = iValue; // ioni

        }

        Eigen::MatrixXd ZBaryszMat,mBaryszMat,vBaryszMat,seBaryszMat, peBaryszMat, areBaryszMat,pBaryszMat,  iBaryszMat;
        // no need to pass the carbon property as it is autonormalized
        if (method==0) {


            ZBaryszMat = computeBaryszMatrix0(mol, Zs_);
            mBaryszMat = computeBaryszMatrix0(mol, ms_);
            vBaryszMat = computeBaryszMatrix0(mol, vs_);
            seBaryszMat = computeBaryszMatrix0(mol, ses_);
            peBaryszMat = computeBaryszMatrix0(mol, pes_);
            areBaryszMat = computeBaryszMatrix0(mol, ares_);
            pBaryszMat = computeBaryszMatrix0(mol, ps_);
            iBaryszMat = computeBaryszMatrix0(mol, is_);

        }

        else if (method==1)
        {

            ZBaryszMat = computeBaryszMatrix1(mol, Zs_);
            mBaryszMat = computeBaryszMatrix1(mol, ms_);
            vBaryszMat = computeBaryszMatrix1(mol, vs_);
            seBaryszMat = computeBaryszMatrix1(mol, ses_);
            peBaryszMat = computeBaryszMatrix1(mol, pes_);
            areBaryszMat = computeBaryszMatrix1(mol, ares_);
            pBaryszMat = computeBaryszMatrix1(mol, ps_);
            iBaryszMat = computeBaryszMatrix1(mol, is_);

        }

        else if (method==2)
        {

            ZBaryszMat = computeBaryszMatrix2(mol, Zs_);
            mBaryszMat = computeBaryszMatrix2(mol, ms_);
            vBaryszMat = computeBaryszMatrix2(mol, vs_);
            seBaryszMat = computeBaryszMatrix2(mol, ses_);
            peBaryszMat = computeBaryszMatrix2(mol, pes_);
            areBaryszMat = computeBaryszMatrix2(mol, ares_);
            pBaryszMat = computeBaryszMatrix2(mol, ps_);
            iBaryszMat = computeBaryszMatrix2(mol, is_);

        }


        std::vector<double> zBaryszDesc = MatrixDescs(mol, ZBaryszMat);
        std::vector<double> mBaryszDesc = MatrixDescs(mol, mBaryszMat);
        std::vector<double> vBaryszDesc = MatrixDescs(mol, vBaryszMat);
        std::vector<double> seBaryszDesc = MatrixDescs(mol, seBaryszMat);
        std::vector<double> peBaryszDesc = MatrixDescs(mol, peBaryszMat);
        std::vector<double> areBaryszDesc = MatrixDescs(mol, areBaryszMat);
        std::vector<double> pBaryszDesc = MatrixDescs(mol, pBaryszMat);
        std::vector<double> iBaryszDesc = MatrixDescs(mol, iBaryszMat);


        // Concatenate all descriptors vector dimension is : 8*13
        std::vector<double> concatenatedDescriptors;
        concatenatedDescriptors.insert(concatenatedDescriptors.end(), zBaryszDesc.begin(), zBaryszDesc.end());
        concatenatedDescriptors.insert(concatenatedDescriptors.end(), mBaryszDesc.begin(), mBaryszDesc.end());
        concatenatedDescriptors.insert(concatenatedDescriptors.end(), vBaryszDesc.begin(), vBaryszDesc.end());
        concatenatedDescriptors.insert(concatenatedDescriptors.end(), seBaryszDesc.begin(), seBaryszDesc.end());
        concatenatedDescriptors.insert(concatenatedDescriptors.end(), peBaryszDesc.begin(), peBaryszDesc.end());
        concatenatedDescriptors.insert(concatenatedDescriptors.end(), areBaryszDesc.begin(), areBaryszDesc.end());
        concatenatedDescriptors.insert(concatenatedDescriptors.end(), pBaryszDesc.begin(), pBaryszDesc.end());
        concatenatedDescriptors.insert(concatenatedDescriptors.end(), iBaryszDesc.begin(), iBaryszDesc.end());


        return concatenatedDescriptors;



    }


    //// Baryz matrix lapack version



// Function to compute eigenvalues and eigenvectors using LAPACK
void compute_eigenvalues_and_eigenvectorsL(
    std::vector<std::vector<double>>& matrix,
    std::vector<double>& eigenvalues,
    std::vector<std::vector<double>>& eigenvectors) {
    int n = matrix.size();
    eigenvalues.resize(n);
    eigenvectors = matrix; // Copy matrix to preserve the original

    // Convert the 2D vector to a 1D array in column-major order for LAPACK
    std::vector<double> flatMatrix(n * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            flatMatrix[j * n + i] = eigenvectors[i][j];




    // Call LAPACKE_dsyev to compute eigenvalues and eigenvectors
    int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, flatMatrix.data(), n, eigenvalues.data());
    if (info != 0) {
        throw std::runtime_error("Error in LAPACKE_dsyev: " + std::to_string(info));
    }

    // Reshape the flatMatrix back into eigenvectors
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            eigenvectors[i][j] = flatMatrix[j * n + i];

    // Canonicalize eigenvector signs: force first non-zero entry positive per column
    const double canonical_eps = 1e-12;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            double v = eigenvectors[i][j];
            if (std::abs(v) > canonical_eps) {
                if (v < 0.0) {
                    for (int k = 0; k < n; ++k) {
                        eigenvectors[k][j] = -eigenvectors[k][j];
                    }
                }
                break;
            }
        }
    }
}


// Spectral Absolute Sum
double spAbsL(const std::vector<double>& eigenvalues) {
    return std::accumulate(eigenvalues.begin(), eigenvalues.end(), 0.0, [](double acc, double val) {
        return acc + std::abs(val);
    });
}

// Leading Eigenvalue
double spMaxL(const std::vector<double>& eigenvalues) {
    return *std::max_element(eigenvalues.begin(), eigenvalues.end());
}

// Spectral Diameter
double spDiamL(const std::vector<double>& eigenvalues) {
    auto [minIt, maxIt] = std::minmax_element(eigenvalues.begin(), eigenvalues.end());
    return *maxIt - *minIt;
}

// Mean of Eigenvalues
double spMeanL(const std::vector<double>& eigenvalues) {
    return std::accumulate(eigenvalues.begin(), eigenvalues.end(), 0.0) / eigenvalues.size();
}

// Spectral Absolute Deviation
double spADL(const std::vector<double>& eigenvalues, double mean) {
    return std::accumulate(eigenvalues.begin(), eigenvalues.end(), 0.0, [mean](double acc, double val) {
        return acc + std::abs(val - mean);
    });
}

double logEEL(const std::vector<double>& eigenvalues) {
    if (eigenvalues.empty()) {
        throw std::runtime_error("Eigenvalues vector is empty");
    }

    double maxVal = *std::max_element(eigenvalues.begin(), eigenvalues.end());
    double sumExp = 0.0;

    for (double val : eigenvalues) {
        sumExp += std::exp(val - maxVal);
    }

    // Avoid issues with log(0) if sumExp is very small
    if (sumExp < 1e-10) {
        return maxVal; // If sumExp is tiny, Log_EE reduces to maxVal
    }

    return maxVal + std::log(sumExp);
}

double logEE_stable(const std::vector<double>& eigenvalues, double threshold = 1e-10) {
    // Handle small eigenvalues close to zero
    std::vector<double> adjustedEigenvalues = eigenvalues;
    for (auto& val : adjustedEigenvalues) {
        if (std::abs(val) < threshold) {
            val = 0.0;
        }
    }

    double maxVal = *std::max_element(adjustedEigenvalues.begin(), adjustedEigenvalues.end());
    double sumExp = std::accumulate(adjustedEigenvalues.begin(), adjustedEigenvalues.end(), 0.0,
        [maxVal](double acc, double val) {
            return acc + std::exp(val - maxVal);
        });
    return maxVal + std::log(sumExp);
}

// Trace of Matrix
double SM1L(const std::vector<std::vector<double>>& matrix) {
    double trace = 0.0;
    for (size_t i = 0; i < matrix.size(); ++i) {
        trace += matrix[i][i];
    }
    return trace;
}

// Coefficient Sum of the Last Eigenvector
double VE1L(const std::vector<std::vector<double>>& eigenvectors) {
    size_t n = eigenvectors.size();
    //const auto& lastEigenvector = eigenvectors[n - 1];

    std::vector<double> lastEigenvector(n);
    for (unsigned int i = 0; i < n; ++i) {
        lastEigenvector[i] = eigenvectors[i][n - 1];
    }

    return std::accumulate(lastEigenvector.begin(), lastEigenvector.end(), 0.0, [](double acc, double val) {
        return acc + std::abs(val);
    });
}

// Average Coefficient of the Last Eigenvector
double VE2L(double ve1, int numAtoms) {
    return ve1 / numAtoms;
}

// Logarithmic Coefficient Sum of the Last Eigenvector
double VE3L(double ve1, int numAtoms) {
    return std::log(0.1 * numAtoms * ve1);
}

// Randic-like Eigenvector-Based Index
double VR1L(const std::vector<std::vector<double>>& eigenvectors,
           const std::vector<std::pair<int, int>>& bonds) {
    size_t n = eigenvectors.size();

    std::vector<double> lastEigenvector(n);
    for (size_t i = 0; i < n; ++i) {
        lastEigenvector[i] = eigenvectors[i][n - 1];
    }

    //const auto& lastEigenvector = eigenvectors[n - 1];
    double result = 0.0;

    for (const auto& bond : bonds) {
        int i = bond.first;
        int j = bond.second;
        double product = std::abs(lastEigenvector[i]) * std::abs(lastEigenvector[j]);
        result += 1.0 / std::sqrt(product);
    }

    return result;
}

// Normalized Randic-like Eigenvector-Based Index
double VR2L(double vr1, int numAtoms) {
    return vr1 / numAtoms;
}

// Logarithmic Randic-like Eigenvector-Based Index
double VR3L(double vr1, int numAtoms) {
    return std::log(0.1 * numAtoms * vr1);
}


std::vector<std::vector<double>> floydWarshallL(std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();

    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (matrix[i][k] < std::numeric_limits<double>::infinity() &&
                    matrix[k][j] < std::numeric_limits<double>::infinity()) {
                    matrix[i][j] = std::min(matrix[i][j], matrix[i][k] + matrix[k][j]);
                }
            }
        }
    }

    return matrix;
}



std::vector<std::vector<double>> computeBaryszMatrix0L(
    const RDKit::ROMol& mol,
    const std::vector<double>& atomicProps) {
    const unsigned int numAtoms = mol.getNumAtoms();
    double largeValue = 1e6;
    // int r = test();
    // Initialize adjacency matrix
    std::vector<std::vector<double>> baryszMatrix(numAtoms, std::vector<double>(numAtoms, largeValue));
    for (unsigned int i = 0; i < numAtoms; ++i) {
        baryszMatrix[i][i] = 0.0;
    }

    // Compute reciprocal atomic properties
    std::vector<double> reciprocalAtomicProps(numAtoms);
    for (unsigned int i = 0; i < numAtoms; ++i) {
        reciprocalAtomicProps[i] = 1.0 / atomicProps[i];
    }

    // Populate adjacency matrix
    for (const auto& bond : mol.bonds()) {
        unsigned int i = bond->getBeginAtomIdx();
        unsigned int j = bond->getEndAtomIdx();
        double pi = bond->getBondTypeAsDouble();
        double w = reciprocalAtomicProps[i] * reciprocalAtomicProps[j] / pi;
        baryszMatrix[i][j] = w;
        baryszMatrix[j][i] = w;
    }

    // Apply Floyd-Warshall
    baryszMatrix = floydWarshallL(baryszMatrix);

    // Update diagonal
    for (unsigned int i = 0; i < numAtoms; ++i) {
        baryszMatrix[i][i] = 1.0 - 1.0 / atomicProps[i];
    }

    // Replace large values with 0
    for (unsigned int i = 0; i < numAtoms; ++i) {
        for (unsigned int j = 0; j < numAtoms; ++j) {
            if (baryszMatrix[i][j] == largeValue) {
                baryszMatrix[i][j] = 0.0;
            }
        }
    }

    return baryszMatrix;
}

std::vector<double> MatrixDescsL(const RDKit::ROMol& mol, std::vector<std::vector<double>>& matrix) {
    int numAtoms = mol.getNumAtoms();

    std::vector<double> eigenvalues;
    std::vector<std::vector<double>> eigenvectors;

    compute_eigenvalues_and_eigenvectorsL(matrix, eigenvalues, eigenvectors);

    std::vector<std::pair<int, int>> bonds;
    for (const auto& bond : mol.bonds()) {
        bonds.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
    }

    double Sp_Abs = spAbsL(eigenvalues);
    double Sp_Max = spMaxL(eigenvalues);
    double Sp_Diam = spDiamL(eigenvalues);
    double Sp_Mean = spMeanL(eigenvalues);
    double Sp_AD = spADL(eigenvalues, Sp_Mean);
    double Sp_MAD = Sp_AD / numAtoms;
    double Log_EE = logEE_stable(eigenvalues);
    double sm1 = SM1L(matrix);
    double ve1 = VE1L(eigenvectors);
    double ve2 = VE2L(ve1, numAtoms);
    double ve3 = VE3L(ve1, numAtoms);
    double vr1 = VR1L(eigenvectors, bonds);
    double vr2 = VR2L(vr1, numAtoms);
    double vr3 = VR3L(vr1, numAtoms);

    return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE, sm1, ve1, ve2, ve3, vr1, vr2, vr3};
}



std::vector<double> calcBaryszMatrixDescsL(const RDKit::ROMol& mol) {
    const PeriodicTable* tbl = PeriodicTable::getTable();

    std::map<int, double> vdwmap = VdWAtomicMap();
    std::map<int, double> semap = SandersonENAtomicMap();
    std::map<int, double> pemap = PaulingENAtomicMap();
    std::map<int, double> aremap = Allred_rocow_ENAtomicMap();
    std::map<int, double> pmap = Polarizability94AtomicMap();
    std::map<int, double> imap = ionizationEnergyAtomicMap();

    double zcc = static_cast<double>(tbl->getAtomicNumber("C"));
    double mcc = static_cast<double>(tbl->getAtomicWeight("C"));
    double mvdwc = vdw_volume(vdwmap[6]);
    double msec = semap[6];
    double mpec = pemap[6];
    double marec = aremap[6];
    double mpc = pmap[6];
    double mic = imap[6];

    std::vector<double> Zs_(mol.getNumAtoms(), 0.0);
    std::vector<double> ms_(mol.getNumAtoms(), 0.0);
    std::vector<double> vs_(mol.getNumAtoms(), 0.0);
    std::vector<double> ses_(mol.getNumAtoms(), 0.0);
    std::vector<double> pes_(mol.getNumAtoms(), 0.0);
    std::vector<double> ares_(mol.getNumAtoms(), 0.0);
    std::vector<double> ps_(mol.getNumAtoms(), 0.0);
    std::vector<double> is_(mol.getNumAtoms(), 0.0);

    for (const auto& atom : mol.atoms()) {
        int atzi = atom->getAtomicNum();
        Zs_[atom->getIdx()] = static_cast<double>(atzi) / zcc;
        ms_[atom->getIdx()] = static_cast<double>(tbl->getAtomicWeight(atzi)) / mcc;
        vs_[atom->getIdx()] = (vdwmap.find(atzi) != vdwmap.end()) ? vdw_volume(vdwmap.at(atzi)) / mvdwc : 0.0;
        ses_[atom->getIdx()] = (semap.find(atzi) != semap.end()) ? semap.at(atzi) / msec : 0.0;
        pes_[atom->getIdx()] = (pemap.find(atzi) != pemap.end()) ? pemap.at(atzi) / mpec : 0.0;
        ares_[atom->getIdx()] = (aremap.find(atzi) != aremap.end()) ? aremap.at(atzi) / marec : 0.0;
        ps_[atom->getIdx()] = (pmap.find(atzi) != pmap.end()) ? pmap.at(atzi) / mpc : 0.0;
        is_[atom->getIdx()] = (imap.find(atzi) != imap.end()) ? imap.at(atzi) / mic : 0.0;
    }

    std::vector<std::vector<std::vector<double>>> matrices(8);
    matrices[0] = computeBaryszMatrix0L(mol, Zs_);
    matrices[1] = computeBaryszMatrix0L(mol, ms_);
    matrices[2] = computeBaryszMatrix0L(mol, vs_);
    matrices[3] = computeBaryszMatrix0L(mol, ses_);
    matrices[4] = computeBaryszMatrix0L(mol, pes_);
    matrices[5] = computeBaryszMatrix0L(mol, ares_);
    matrices[6] = computeBaryszMatrix0L(mol, ps_);
    matrices[7] = computeBaryszMatrix0L(mol, is_);

    std::vector<double> concatenatedDescriptors;
    for (auto& matrix : matrices) {
        auto desc = MatrixDescsL(mol, matrix);
        concatenatedDescriptors.insert(concatenatedDescriptors.end(), desc.begin(), desc.end());
    }

    return concatenatedDescriptors;
}

///TopologicalCharge


    // Function to get the Adjacency Matrix
    Eigen::MatrixXd calculateAdjacencyMatrix(const RDKit::ROMol& mol) {
        unsigned int nAtoms = mol.getNumAtoms();
        Eigen::MatrixXd adjMatrix(nAtoms, nAtoms);
        adjMatrix.setZero();

        // Populate the adjacency matrix using RDKit's getAdjacencyMatrix is this faster than calling
        for (unsigned int i = 0; i < nAtoms; ++i) {
            for (unsigned int j = 0; j < nAtoms; ++j) {
                const RDKit::Bond* bond = mol.getBondBetweenAtoms(i, j);
                adjMatrix(i, j) = (bond != nullptr) ? 1.0 : 0.0;
            }
        }

        return adjMatrix;
    }

    Eigen::MatrixXd calculateChargeTermMatrix(const Eigen::MatrixXd& A, const Eigen::MatrixXd& D) {
        // Step 1: Invert non-zero elements of D^2 and zero out the diagonal
        Eigen::MatrixXd D2 = (D.array() != 0).select(D.array().pow(-2), 0.0);
        D2.diagonal().setZero();

        // Step 2: Perform matrix multiplication and antisymmetric operation
        Eigen::MatrixXd M = A * D2;
        return M - M.transpose();
    }

    std::vector<double> calcTopologicalChargeDescs(const RDKit::ROMol& mol) {
        const int maxOrder = 10;
        std::vector<double> results(21, 0.0);

        // Calculate Distance and Adjacency Matrices
        Eigen::MatrixXd distanceMatrix = calculateDistanceMatrix(mol);
        Eigen::MatrixXd adjacencyMatrix = calculateAdjacencyMatrix(mol);

        // Step 1: Compute the Charge Term Matrix
        Eigen::MatrixXd CT = calculateChargeTermMatrix(adjacencyMatrix, distanceMatrix);

        // Step 3: Create a lower triangular distance matrix
        Eigen::MatrixXd lowerTriangularMask = Eigen::MatrixXd::Zero(distanceMatrix.rows(), distanceMatrix.cols());
        for (int i = 0; i < distanceMatrix.rows(); ++i) {
            for (int j = 0; j < i; ++j) {
                lowerTriangularMask(i, j) = 1.0;
            }
        }
        Eigen::MatrixXd D = distanceMatrix.cwiseProduct(lowerTriangularMask);
        D = D.unaryExpr([](double x) { return (x == 0) ? std::numeric_limits<double>::infinity() : x; });

        // Step 3: Calculate raw, mean, and global descriptors
        for (int order = 1; order <= maxOrder; ++order) {
            // Create mask for the current order
            Eigen::MatrixXd orderMask = (D.array() == order).cast<double>();

            // Apply mask to CT / D
            std::vector<double> filteredCT, filteredD;

            for (int i = 0; i < CT.rows(); ++i) {
                for (int j = 0; j < CT.cols(); ++j) {
                    if (orderMask(i, j) > 0) {
                        filteredCT.push_back(CT(i, j));
                        filteredD.push_back(D(i, j));
                    }
                }
            }

            // Raw descriptor: Absolute sum of filtered CT values
            double raw = 0.0;
            for (double val : filteredCT) {
                raw += std::abs(val);
            }
            results[order - 1] = raw;

            // Mean ie Normalize by frequencies
            if (!filteredCT.empty()) {
                // dict/Map to store frequencies of each distance value
                std::unordered_map<double, int> frequencies;
                for (double val : filteredD) {
                    frequencies[val]++;
                }

                // get mean descriptor
                double mean = 0.0;
                for (size_t i = 0; i < filteredCT.size(); ++i) {
                    mean += std::abs( filteredCT[i]) / frequencies[filteredD[i]];
                }
                results[10 + order - 1] = mean;
            }
        }

        // Global descriptor: full Sum of all mean values from 1 to maxOrder

        Eigen::MatrixXd orderMask = (D.array() <= maxOrder).cast<double>();

        // Apply mask to CT
        std::vector<double> gfilteredD, gfilteredCT;
        for (int i = 0; i < CT.rows(); ++i) {
            for (int j = 0; j < CT.cols(); ++j) {
                if (orderMask(i, j) > 0) {
                    gfilteredCT.push_back(CT(i, j));
                    gfilteredD.push_back(D(i, j));
                }
            }
        }


        // Mean descriptor: Normalize by frequencies
        if (!gfilteredCT.empty()) {
            // Map to store frequencies of each distance value
            std::unordered_map<double, int> gfrequencies;
            for (double val : gfilteredD) {
                gfrequencies[val]++;
            }
            double global = 0.0;
            // Calculate mean descriptor
            for (size_t i = 0; i < gfilteredCT.size(); ++i) {
                global += std::abs(gfilteredCT[i]) / gfrequencies[ gfilteredD[i]];
            }
            results[20] = global;
        }

        return results;
    }


    /// not tested TopologocalIndex (4 outputs => need eccentricities computation)

    // Function to calculate the radius (minimum of eccentricities)
    double calculateRadius(const Eigen::VectorXd& eccentricities) {
        if (eccentricities.size() == 0) {
            throw std::invalid_argument("Eccentricity data is empty");
        }
        return eccentricities.minCoeff();  // Eigen function to get the minimum value
    }

    // Function to calculate the diameter (maximum of eccentricities)
    double calculateDiameter(const Eigen::VectorXd& eccentricities) {
        if (eccentricities.size() == 0) {
            throw std::invalid_argument("Eccentricity data is empty");
        }
        return eccentricities.maxCoeff();  // Eigen function to get the maximum value
    }


    // Function to calculate the Topological Shape Index
    double calculateTopologicalShapeIndex(double& r, const double& d) {

        if (r == 0.0) {
            return 0.;
            //throw std::runtime_error("Division by zero: radius is 0");
        }

        return (d - r) / r;
    }

    // Function to calculate the Petitjean Index
    double calculatePetitjeanIndex(const double& r, const double& d) {

        if (d == 0.0) {
            return 0.;
            //throw std::runtime_error("Division by zero: diameter is 0");
        }

        return (d - r) / d;
    }

    std::vector<double> calcTopologicalIndex(const RDKit::ROMol& mol) {

        Eigen::VectorXd Eccentricity = calculateEccentricity(mol) ;
        double Radius = calculateRadius(Eccentricity);
        double Diameter = calculateDiameter(Eccentricity);
        double TSI = calculateTopologicalShapeIndex(Radius,Diameter);
        double PJI =  calculatePetitjeanIndex(Radius,Diameter);
        return {Diameter, Radius, TSI, PJI};
    }






///// CHI / Path deps

    // Enum for ChiType
    enum class ChiType {
        Path = 1,
        Cluster,
        PathCluster,
        Chain
    };


    // Function to convert ChiType to string
    std::string toString(ChiType type) {
        switch (type) {
            case ChiType::Path: return "Path";
            case ChiType::Cluster: return "Cluster";
            case ChiType::PathCluster: return "PathCluster";
            case ChiType::Chain: return "Chain";
            default: return "Unknown";
        }
    }

    void performDFS(const RDKit::ROMol& mol,
                    int startAtomIdx,
                    const std::vector<int>& path,
                    std::set<int>& visitedNodes,
                    std::set<std::pair<int, int>>& visitedEdges,
                    std::set<int>& degrees,
                    bool& isChain) {
        std::set<int> pathBonds(path.begin(), path.end()); // Bonds in the path for quick lookup
        std::unordered_map<int, std::set<int>> neighbors; // Neighbors in the subgraph

        // Populate neighbors for the subgraph
        for (int bondIdx : path) {
            const auto* bond = mol.getBondWithIdx(bondIdx);
            int begin = bond->getBeginAtomIdx();
            int end = bond->getEndAtomIdx();
            neighbors[begin].insert(end);
            neighbors[end].insert(begin);
        }

        // Perform DFS
        std::stack<int> stack;
        std::unordered_map<int, int> parent; // To track parent nodes in DFS
        stack.push(startAtomIdx);
        parent[startAtomIdx] = -1; // Root node has no parent

        while (!stack.empty()) {
            int node = stack.top();
            stack.pop();

            if (visitedNodes.count(node)) {
                continue;
            }

            visitedNodes.insert(node);

            // Calculate degree for this node based on subgraph neighbors
            int degree = neighbors[node].size();
            degrees.insert(degree); // Add degree to the set

            // Traverse neighbors in the subgraph
            for (int neighbor : neighbors[node]) {
                std::pair<int, int> edge = std::minmax(node, neighbor);

                if (!visitedNodes.count(neighbor)) {
                    stack.push(neighbor);
                    parent[neighbor] = node; // Set parent for the neighbor
                    visitedEdges.insert(edge);
                } else if (parent[node] != neighbor) { // Detect back edge
                    isChain = true; // Cycle detected
                }
            }
        }
    }

    bool allDegreesAreOneOrTwo(const std::set<int>& degrees) {
        return std::all_of(degrees.begin(), degrees.end(), [](int d) {
            return d == 1 || d == 2;
        });
    }

    // Main function to extract and classify subgraphs
    std::vector<std::tuple<std::vector<int>, std::set<int>, ChiType>> extractAndClassifyPaths(
        const RDKit::ROMol& mol, unsigned int targetLength, bool useHs) {

        std::vector<std::tuple<std::vector<int>, std::set<int>, ChiType>> results;

        // Get all subgraphs of the given length
        auto paths = findAllSubgraphsOfLengthN(mol, targetLength, useHs, -1); // return both path and complex subgraphs AllPath return only true Path ...! maybe we can leverage that except if it is too expensive...

        for (const auto& path : paths) {
            // Prepare sets for DFS traversal
            std::set<int> visitedNodes;
            std::set<std::pair<int, int>> visitedEdges;
            std::set<int> degrees;
            bool isChain = false;

            // Start DFS from the first bond in the path
            if (!path.empty()) {
                int startAtomIdx = mol.getBondWithIdx(path.front())->getBeginAtomIdx();
                performDFS(mol, startAtomIdx, path, visitedNodes, visitedEdges, degrees, isChain);
            }

            // If a cycle is detected, it's a Chain Path by definitin of the isChain bool flag from DFS code
            // this is a decision tree: Chain first than only 1 and 2 => Path than has 2 => Path Cluster else Cluster!
            ChiType type;
            if (isChain) {
                type = ChiType::Chain;
            } else if (allDegreesAreOneOrTwo(degrees)) {
                type = ChiType::Path;
            } else if (degrees.count(2)) {
                type = ChiType::PathCluster;
            } else {
                type = ChiType::Cluster;
            }
            results.emplace_back(path, visitedNodes, type);
        }
        return results;
    }


    // first code slow V1
    int calcPathsOfLengthN_(const RDKit::ROMol& mol, int order) {
        // Extract and classify subgraphs for the current radius order
        auto classifiedPaths = extractAndClassifyPaths(mol, order, false);
        int j = 0;
        for (const auto& [bonds, nodes, type] : classifiedPaths) {
            if (type == ChiType::Path) {
                j++;
            }
        }
        return j;
    }


    // second code  V2 faster than extractAndClassifyPaths
    int calcPathsOfLengthN(const RDKit::ROMol& mol, int order) {
        // Extract and classify subgraphs for the current radius order
        int j = 0;

        auto paths = findAllPathsOfLengthN(mol, order+1, false, false, -1, false); // Atoms indices we need +1 at it is linear path order bonds equal order+1 atoms paths!!!!

        for (const auto& atomPath : paths) {
            std::unordered_set<int> visitedAtoms; // Set to track visited atoms
            bool isDuplicate = false;

            for (size_t i = 0; i < atomPath.size(); ++i) {
                if (visitedAtoms.count(atomPath[i])) {
                    isDuplicate = true;
                    break;
                }
                visitedAtoms.insert(atomPath[i]);
            }

            if (!isDuplicate) {
                j++;         // Increment true path count and exclude Chain
            }
        }
        return j;
    }

    // thrid code no need for Iterator ... so very slightly faster then V2
    int calcPathsOfLengthN__(const RDKit::ROMol& mol, int order) {
        // Count for true paths
        int pathCount = 0;

        // Find linear atom paths (order + 1 atoms for a path of length `order`)
        auto atomPaths = findAllPathsOfLengthN(mol, order + 1, false, false, -1, false);

        for (const auto& atomPath : atomPaths) {
            // Check for duplicates in the atom path
            std::unordered_set<int> visitedAtoms(atomPath.begin(), atomPath.end());
            bool isChain = (visitedAtoms.size() < atomPath.size());

            if (!isChain) {
                pathCount++; // Increment path count for true linear paths
            }
        }

        return pathCount;
    }



// Function to convert bond paths into atom paths
std::vector<int> bondPathToAtomPath(const RDKit::ROMol& mol, const std::vector<int>& bondPath) {
    std::unordered_set<int> visitedAtoms;
    std::vector<int> atomPath;

    for (int bondIdx : bondPath) {
        const auto* bond = mol.getBondWithIdx(bondIdx);
        int begin = bond->getBeginAtomIdx();
        int end = bond->getEndAtomIdx();

        if (visitedAtoms.insert(begin).second) atomPath.push_back(begin);
        if (visitedAtoms.insert(end).second) atomPath.push_back(end);
    }

    return atomPath;
}


struct pathHash {
    std::size_t operator()(const std::vector<int>& path) const {
        std::size_t seed = 0;
        for (int i : path) {
            seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};


ChiType classifySubgraph(const std::set<int>& degrees, bool isChain) {
    // Classification logic based on degrees and cycle detection
    if (isChain) {
        return ChiType::Chain;
    } else if (allDegreesAreOneOrTwo(degrees)) {
        return ChiType::Path;
    } else if (degrees.count(2)) {
        return ChiType::PathCluster;
    } else {
        return ChiType::Cluster;
    }
}


bool detectCycle(const std::vector<int>& atomPath) {
    std::unordered_set<int> visited;
    for (int atomIdx : atomPath) {
        if (visited.count(atomIdx)) {
            return true; // Cycle detected
        }
        visited.insert(atomIdx);
    }
    return false; // No cycle
}

ChiType classifySubgraph(const RDKit::ROMol& mol, const std::vector<int>& bondPath) {
    // Map to store atom index mapping between subgraph and parent molecule
    INT_MAP_INT atomIdxMap;

    // Create a submolecule for the given bond path
    std::unique_ptr<RDKit::ROMol> subMol(Subgraphs::pathToSubmol(mol, bondPath, false, atomIdxMap));
    // Collect degrees in the subgraph
    std::set<int> degrees;
    for (const auto& atom : subMol->atoms()) {
        degrees.insert(atom->getDegree());
    }

    // Determine if the subgraph is a chain
    bool isChain = detectCycle(bondPath);

    // Apply classification logic
    if (isChain) {
        return ChiType::Chain;
    } else if (allDegreesAreOneOrTwo(degrees)) {
        return ChiType::Path;
    } else if (degrees.count(2)) {
        return ChiType::PathCluster;
    } else {
        return ChiType::Cluster;
    }
}


std::vector<std::pair<std::vector<int>, ChiType>> extractAndClassifyPathsAndSubgraphs(
    const RDKit::ROMol& mol,
    unsigned int targetLength,
    bool useHs
) {
    std::vector<std::pair<std::vector<int>, ChiType>> results;

    // Step 1: Classify linear atom paths
    auto atomPaths = findAllPathsOfLengthN(mol, targetLength + 1, useHs, false, -1, false);

    for (const auto& atomPath : atomPaths) {
        bool isChain = detectCycle(atomPath); // Detect cycles in atom path
        ChiType type = isChain ? ChiType::Chain : ChiType::Path;

        results.emplace_back(atomPath, type);
    }

    // Step 2: Classify bond subgraphs
    auto bondSubgraphs = findAllSubgraphsOfLengthN(mol, targetLength, useHs, -1);

    for (const auto& bondPath : bondSubgraphs) {
        // Convert bond subgraph to atom path
        std::vector<int> atomPath = bondPathToAtomPath(mol, bondPath);

        // Check if the atom path was already classified as Path or Chain
        auto it = std::find_if(results.begin(), results.end(), [&](const auto& pair) {
            return pair.first == atomPath && (pair.second == ChiType::Path || pair.second == ChiType::Chain);
        });

        if (it != results.end()) {
            continue; // Skip already classified paths
        }

        // Classify the bond subgraph
        ChiType type = classifySubgraph(mol, bondPath);

        // Only add the subgraph if it is not a Path or Chain
        if (type != ChiType::Path && type != ChiType::Chain) {
            results.emplace_back(atomPath, type);
        }
    }

    return results;
}

// Calculate path counts and weighted product
std::pair<int, double> calculatePathCount(const RDKit::ROMol& mol, int order) {
    int L = 0;          // Path count
    double piSum = 0.0; // Weighted bond product sum

    // Get all paths of the given length
    auto paths = findAllPathsOfLengthN(mol, order+1, false, false, -1, false); // Atoms indices we need +1 at it is linear path order bonds equal order+1 atoms paths!!!!

    for (const auto& atomPath : paths) {
        std::unordered_set<int> visitedAtoms; // Set to track visited atoms
        bool isDuplicate = false;

        double bondProduct = 1.0; // Initialize bond product

        for (size_t i = 0; i < atomPath.size(); ++i) {
            if (visitedAtoms.count(atomPath[i])) {
                isDuplicate = true;
                break;
            }
            visitedAtoms.insert(atomPath[i]);

            if (i > 0) {
                const auto* bond = mol.getBondBetweenAtoms(atomPath[i - 1], atomPath[i]);
                bondProduct *= bond->getBondTypeAsDouble();
            }
        }

        if (!isDuplicate) {
            L++;         // Increment path count
            piSum += bondProduct; // Add bond product to the sum
        }
    }

    return {L, piSum};
}

// Main function to calculate path descriptors
std::vector<double> calcPathCount(const RDKit::ROMol& mol) {
    std::vector<double> results(21, 0.0); // Output vector
    int totalMPC = mol.getNumAtoms();     // Initialize MPC1
    double cumulativePiSum = static_cast<double>(mol.getNumAtoms()); // Initialize piPC1

    for (int order = 1; order <= 10; ++order) {
        auto [L, piSum] = calculatePathCount(mol, order);

        if (order > 1) {
            results[order - 2] = L; // MPC2-MPC10 in indices 0-8
        }
        results[10 + order - 1] = std::log(piSum + 1.0); // piPC1-piPC10 in indices 10-19

        totalMPC += L;
        cumulativePiSum += piSum;
    }

    results[9] = totalMPC;                        // Total MPC in index 9
    results[20] = std::log(cumulativePiSum + 1); // Total piPC in index 20

    return results;
}

    // caution this is orignial kappa shape index so without the alpha term p 428 / 429!

    double kappa1Helperfix(double P1, double A) {
        double denom = P1 ;
        double kappa = 0.0;
        if (denom) {
            kappa = ( A ) * (A - 1) * (A - 1) / (denom * denom);
        }
        return kappa;
    }
    double calcKappa1fix(const ROMol &mol) {
        double P1 = mol.getNumBonds();
        double A = mol.getNumHeavyAtoms();
        double kappa = kappa1Helperfix(P1, A);
        return kappa;
    }

    double kappa2Helperfix(double P2, double A) {
        double denom = (P2 ) * (P2);
        double kappa = 0.0;
        if (denom) {
            kappa = (A - 1) * (A - 2) * (A - 2) / denom;
        }
        return kappa;
    }
    double calcKappa2fix(const ROMol &mol) {
        //double P2  = findAllPathsOfLengthN(mol, 2).size();
        double P2 = static_cast<double>(calcPathsOfLengthN(mol, 2));

        double A = mol.getNumHeavyAtoms();
        double kappa = kappa2Helperfix(P2, A);
        return kappa;
    }

    double kappa3Helperfix(double P3, int A) {
        double denom = P3  * P3;
        double kappa = 0.0;
        if (denom) {
            if (A % 2) {
                kappa = (A - 1) * (A - 3) * (A - 3) / denom;
            } else {
                kappa = (A - 2) * (A - 2) * (A - 3) / denom;
            }
        }
        return kappa;
        }


    double calcKappa3fix(const ROMol &mol) {
        //double P3 = findAllPathsOfLengthN(mol, 3).size();
        double P3 = static_cast<double>(calcPathsOfLengthN(mol, 3));

        double A = mol.getNumHeavyAtoms();
        double kappa = kappa3Helperfix(P3, A);
        return kappa;
    }

    std::vector<double> calcKappaShapeIndex(const RDKit::ROMol& mol) {
        std::vector<double> results(3, 0.0);
        results[0] = calcKappa1fix(mol);

        results[1] = calcKappa2fix(mol);

        results[2] = calcKappa3fix(mol);

    return results;
    }







    // it is not identical to "getRKHE"
    void hkDeltas(const ROMol& mol, std::vector<double>& deltas,  bool force) {
        PRECONDITION(deltas.size() >= mol.getNumAtoms(), "bad vector size");
        if (!force && mol.hasProp("_connectivityHKDeltas")) {
            mol.getProp("_connectivityHKDeltas", deltas);
            return;
        }

        // Compute the valence electrons for each atom
        ROMol::VERTEX_ITER atBegin, atEnd;
        boost::tie(atBegin, atEnd) = mol.getVertices();
        while (atBegin != atEnd) {
            const Atom* at = mol[*atBegin];
            double dv = getValenceElectrons(*at);  // Valence electrons

            // Replace previous delta calculation with valence electrons
            if (dv != 0.0) {
                deltas[at->getIdx()] = 1. / sqrt(dv);
            } else {
                deltas[at->getIdx()] = 0.0;
            }
            ++atBegin;
        }

        mol.setProp("_connectivityHKDeltas", deltas, true);
    }


    void nVals(const ROMol& mol, std::vector<double>& nVs, bool force) {
        PRECONDITION(nVs.size() >= mol.getNumAtoms(), "bad vector size");
        if (!force && mol.hasProp("_connectivityNVals")) {
            mol.getProp("_connectivityNVals", nVs);
            return;
        }

        // Compute the sigma electrons for each atom
        ROMol::VERTEX_ITER atBegin, atEnd;
        boost::tie(atBegin, atEnd) = mol.getVertices();
        while (atBegin != atEnd) {
            const Atom* at = mol[*atBegin];
            double sigma_electrons = getSigmaElectrons(*at);  // Sigma electrons

            if (sigma_electrons != 0.0) {
                nVs[at->getIdx()] = 1. / sqrt(sigma_electrons);
            } else {
                nVs[at->getIdx()] = 0.0;
            }
            ++atBegin;
        }

        mol.setProp("_connectivityNVals", nVs, true);
    }

    // optional alphaKappa version (for HeteroAtom differentiation)

    double getAlpha(const Atom &atom, bool &found) {
    double res = 0.0;
    found = false;
    switch (atom.getAtomicNum()) {
        case 1:
        res = 0.0;
        found = true;
        break;
        case 6:
        switch (atom.getHybridization()) {
            case Atom::SP:
            res = -0.22;
            found = true;
            break;
            case Atom::SP2:
            res = -0.13;
            found = true;
            break;
            default:
            res = 0.00;
            found = true;
        };
        break;
        case 7:
        switch (atom.getHybridization()) {
            case Atom::SP:
            res = -0.29;
            found = true;
            break;
            case Atom::SP2:
            res = -0.20;
            found = true;
            break;
            default:
            res = -0.04;
            found = true;
            break;
        };
        break;
        case 8:
        switch (atom.getHybridization()) {
            case Atom::SP2:
            res = -0.20;
            found = true;
            break;
            default:
            res = -0.04;
            found = true;
            break;
        };
        break;
        case 9:
        switch (atom.getHybridization()) {
            default:
            res = -0.07;
            found = true;
            break;
        };
        break;
        case 15:
        switch (atom.getHybridization()) {
            case Atom::SP2:
            res = 0.30;
            found = true;
            break;
            default:
            res = 0.43;
            found = true;
            break;
        };
        break;
        case 16:
        switch (atom.getHybridization()) {
            case Atom::SP2:
            res = 0.22;
            found = true;
            break;
            default:
            res = 0.35;
            found = true;
            break;
        };
        break;
        case 17:
        switch (atom.getHybridization()) {
            default:
            res = 0.29;
            found = true;
            break;
        };
        break;
        case 35:
        switch (atom.getHybridization()) {
            default:
            res = 0.48;
            found = true;
            break;
        };
        break;
        case 53:
        switch (atom.getHybridization()) {
            default:
            res = 0.73;
            found = true;
            break;
        };
        break;
        default:
        break;
    }
    return res;
    }



    // old style from connectivity rdkit deps file
    double calcHallKierAlpha(const ROMol &mol) {
    const PeriodicTable *tbl = PeriodicTable::getTable();
    double alphaSum = 0.0;
    double rC = tbl->getRb0(6);
    ROMol::VERTEX_ITER atBegin, atEnd;
    boost::tie(atBegin, atEnd) = mol.getVertices();
    while (atBegin != atEnd) {
        const Atom* at = mol[*atBegin];
        ++atBegin;
        unsigned int n = at->getAtomicNum();
        if (!n) continue;
        bool found;
        double alpha = getAlpha(*at, found);
        if (!found) {
        double rA = tbl->getRb0(n);
        alpha = rA / rC - 1.0;
        }
        alphaSum += alpha;
    }
    return alphaSum;
    };


    double alphakappa1Helper(double P1, double A, double alpha) {
        double denom = P1 + alpha;
        double kappa = 0.0;
        if (denom) {
            kappa = (A + alpha) * (A + alpha - 1) * (A + alpha - 1) / (denom * denom);
        }
        return kappa;
    }
    double alphakappa2Helper(double P2, double A, double alpha) {
        double denom = (P2 + alpha) * (P2 + alpha);
        double kappa = 0.0;
        if (denom) {
            kappa = (A + alpha - 1) * (A + alpha - 2) * (A + alpha - 2) / denom;
        }
        return kappa;
    }

    double alphakappa3Helperfix(double P3, int A, double alpha) {
        double denom = (P3 + alpha) * (P3 + alpha);
        double kappa = 0.0;
        if (denom) {
            if (A % 2) {
            kappa = (A + alpha - 1) * (A + alpha - 3) * (A + alpha - 3) / denom;
            } else {
            kappa = (A + alpha - 2) * (A + alpha - 2) * (A + alpha - 3) / denom;
            }
        }
        return kappa;
    }


    double calcalphaKappa1(const ROMol &mol) {
        double P1 = mol.getNumBonds();
        double A = mol.getNumHeavyAtoms();
        double alpha = calcHallKierAlpha(mol);
        double kappa = alphakappa1Helper(P1, A, alpha);
        return kappa;
    }
    double calcalphaKappa2(const ROMol &mol) {
        double P2 = static_cast<double>(calcPathsOfLengthN(mol, 2));
        double A = mol.getNumHeavyAtoms();
        double alpha = calcHallKierAlpha(mol);
        double kappa = alphakappa2Helper(P2, A, alpha);
        return kappa;
    }
    double calcalphaKappa3fix(const ROMol &mol) {
        double P3 = static_cast<double>(calcPathsOfLengthN(mol, 3));
        int A = mol.getNumHeavyAtoms();
        double alpha = calcHallKierAlpha(mol);
        double kappa = alphakappa3Helperfix(P3, A, alpha);
        return kappa;
    }


    double Flexibility(const RDKit::ROMol &mol) {

        double AK1 = calcalphaKappa1(mol);
        double AK2 = calcalphaKappa2(mol);
        int numHeavyAtom =  mol.getNumHeavyAtoms();
        return AK1*AK2/static_cast<double>(numHeavyAtom);

    }


    std::vector<double> calcFlexibility(const RDKit::ROMol &mol) {
        std::vector<double> res(1,0.);
        res[0] = Flexibility(mol);
        return res;


    }


    std::vector<double> calcAlphaKappaShapeIndex(const RDKit::ROMol& mol) {
        std::vector<double> results(3, 0.0);
        results[0] = calcalphaKappa1(mol);

        results[1] = calcalphaKappa2(mol);

        results[2] = calcalphaKappa3fix(mol);

    return results;
    }








    // Compute for Path (Xp)
    double calcXpDescriptor(const ROMol& mol, unsigned int order, bool useValenceElectrons, bool force) {
        std::vector<double> electronVals(mol.getNumAtoms());
        if (useValenceElectrons) {
            hkDeltas(mol, electronVals, force);  // Valence Electrons (dv)
        } else {
            nVals(mol, electronVals, force);  // Sigma Electrons (d)
        }

        PATH_LIST ps = findAllPathsOfLengthN(mol, order, false);
        double result = 0.0;
        for (const auto& p : ps) {
            double accum = 1.0;
            for (int aidx : p) {
                accum *= electronVals[aidx];  // Multiply the electron values along the path
            }
            result += accum;
        }
        return result;
    }

    // Compute for Chain (Xch)
    double calcXchDescriptor(const ROMol& mol, unsigned int order, bool useValenceElectrons, bool force) {
        std::vector<double> electronVals(mol.getNumAtoms());
        if (useValenceElectrons) {
            hkDeltas(mol, electronVals, force);  // Valence Electrons (dv)
        } else {
            nVals(mol, electronVals, force);  // Sigma Electrons (d)
        }

        PATH_LIST ps = findAllPathsOfLengthN(mol, order, true);  // Chain specific path finding (true means chain only ???)
        double result = 0.0;
        for (const auto& p : ps) {
            double accum = 1.0;
            for (int aidx : p) {
                accum *= electronVals[aidx];
            }
            result += accum;
        }
        return result;
    }


    // Function to compute all descriptors for different ranges
    double computeDescriptors(const ROMol& mol, bool force) {
        // Calculate for Xp (Path)
        for (unsigned int order = 0; order <= 7; ++order) {
            // Compute Xp-0d, Xp-1d, ..., Xp-7d
            double result_d = calcXpDescriptor(mol, order, false, force);
            std::cout << "Xp-" << order << "d: " << result_d << std::endl;

            // Compute Xp-0dv, Xp-1dv, ..., Xp-7dv
            double result_dv = calcXpDescriptor(mol, order, true, force);
            std::cout << "Xp-" << order << "dv: " << result_dv << std::endl;

            // Compute AXp-0d, AXp-1d, ..., AXp-7d
            // Assuming average means some form of averaging over the molecule
            double avg_result_d = result_d / mol.getNumAtoms();
            std::cout << "AXp-" << order << "d: " << avg_result_d << std::endl;

            // Compute AXp-0dv, AXp-1dv, ..., AXp-7dv
            double avg_result_dv = result_dv / mol.getNumAtoms();
            std::cout << "AXp-" << order << "dv: " << avg_result_dv << std::endl;
        }

        // Calculate for Xch (Chain)
        for (unsigned int order = 3; order <= 7; ++order) {
            double result_d = calcXchDescriptor(mol, order, false, force);
            std::cout << "Xch-" << order << "d: " << result_d << std::endl;

            double result_dv = calcXchDescriptor(mol, order, true, force);
            std::cout << "Xch-" << order << "dv: " << result_dv << std::endl;
        }

    return 1.0;
    }


// test function for Chi
    std::vector<int> calcChiPath(const RDKit::ROMol& mol) {
        std::vector<int> res(3, 0);

        for (int i = 1; i <= 3; i++) {
            // Extract and classify subgraphs for the current radius
            auto classifiedPaths = extractAndClassifyPaths(mol, i, false);

            //std::vector<std::vector<int>> atomPaths;

            // Filter only the paths classified as ChiType::Path
            int j = 0;
            for (const auto& [bonds, nodes, type] : classifiedPaths) {
                if (type == ChiType::Path) {
                    j++;
                    //atomPaths.emplace_back(nodes.begin(), nodes.end());
                }
            }

            // Count the number of unique atom paths for this radius
            res[i - 1] = j; //atomPaths.size();
        }

        return res;
    }

    // Function to compute Chi descriptor for a specific ChiType
    double computeChiDescriptor(const RDKit::ROMol&, const PATH_LIST& subgraphs, const std::vector<double>& P, bool averaged = false) {
        double chi_value = 0.0;

        for (const auto& nodes : subgraphs) {
            double c = 1.0;
            for (int node : nodes) {
                c *= P[node];  // Calculate the product of atomic properties
            }

            if (c <= 0) {
                std::cerr << "Error: Some properties are less than or equal to zero." << std::endl;
                return 0.0;  // Return 0 if any property is <= 0
            }

            chi_value += std::pow(c, -0.5);
        }

        // Optionally average the value by the number of subgraphs if `averaged` is true
        if (averaged && !subgraphs.empty()) {
            chi_value /= subgraphs.size();
        }

        return chi_value;
    }




//// DetourMatrix code here



    // Function to compute longest path using DFS
    void longestSimplePath(int u, const Eigen::MatrixXd& adjMatrix, std::vector<double>& result,
                        std::vector<bool>& visited, double distance, int N) {
        visited[u] = true;
        result[u] = std::max(result[u], distance);

        for (int v = 0; v < N; ++v) {
            if (adjMatrix(u, v) > 0 && !visited[v]) {
                longestSimplePath(v, adjMatrix, result, visited, distance + adjMatrix(u, v), N);
            }
        }

        visited[u] = false;
    }

    // Function to compute the detour matrix using an adjacency matrix
    Eigen::MatrixXd computeDetourMatrix(int N, const Eigen::MatrixXd& adjMatrix) {
        Eigen::MatrixXd detourMatrix(N, N);
        detourMatrix.setZero();  // Initialize the matrix with zeros

        // Compute longest path from each node
        for (int i = 0; i < N; ++i) {
            std::vector<double> result(N, 0.0);
            std::vector<bool> visited(N, false);
            longestSimplePath(i, adjMatrix, result, visited, 0.0, N);

            for (int j = 0; j < N; ++j) {
                detourMatrix(i, j) = result[j];
            }
        }

        return detourMatrix;
    }


    // Compute Detour Index
    int computeDetourIndex(const Eigen::MatrixXd& detourmatrix) {
        double sum = detourmatrix.sum();
        //int num_atoms = detourmatrix.rows();
        return static_cast<int>(0.5 * sum);
    }


    std::vector<double> calcDetourMatrixDescs(const RDKit::ROMol& mol) {
        Eigen::MatrixXd AdjMat = calculateAdjacencyMatrix( mol);

        int numAtoms = mol.getNumAtoms();
        Eigen::Matrix detourMatrix = computeDetourMatrix(numAtoms,  AdjMat);

        Eigen::VectorXd eigenvalues;

        Eigen::MatrixXd eigenvectors;

        compute_eigenvalues_and_eigenvectors(detourMatrix, eigenvalues, eigenvectors);




        // Compute descriptors
        double Sp_Abs = spAbs(eigenvalues);
        double Sp_Max = spMax(eigenvalues);
        double Sp_Diam = spDiam(eigenvalues);
        double Sp_Mean = spMean(eigenvalues);
        double Sp_AD = spAD(eigenvalues, Sp_Mean);
        double Sp_MAD = Sp_AD / numAtoms;
        double Log_EE = logEE(eigenvalues);

        std::vector<std::pair<int, int>> bonds;

            // Iterate over all bonds in the molecule
            for (const auto& bond : mol.bonds()) {
                int beginAtomIdx = bond->getBeginAtomIdx();
                int endAtomIdx = bond->getEndAtomIdx();
                bonds.push_back(std::make_pair(beginAtomIdx, endAtomIdx));
            }


        // Calculate the descriptors
        double ve1 = VE1(detourMatrix, eigenvalues, eigenvectors);
        double ve2 = VE2(detourMatrix, numAtoms, eigenvalues, eigenvectors);
        double ve3 = VE3(detourMatrix, numAtoms, eigenvalues, eigenvectors) ;
        double vr1 = VR1(detourMatrix, bonds, eigenvalues, eigenvectors);
        double vr2 = VR2(detourMatrix, bonds, numAtoms, eigenvalues, eigenvectors) ;
        double vr3 = VR3(detourMatrix, bonds, numAtoms, eigenvalues, eigenvectors);
        double detourindex = static_cast<double>(computeDetourIndex(detourMatrix));

        double sm1 = SM1(detourMatrix);
        return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE, sm1, ve1,ve2,ve3,vr1,vr2,vr3,detourindex};


    }

/// much faster lapack version

// Function to compute longest path using DFS
void longestSimplePathL(int u, const std::vector<std::vector<double>>& adjMatrix,
                        std::vector<double>& result, std::vector<bool>& visited,
                        double distance, int N) {
    visited[u] = true;
    result[u] = std::max(result[u], distance);

    for (int v = 0; v < N; ++v) {
        if (adjMatrix[u][v] > 0 && !visited[v]) {
            longestSimplePathL(v, adjMatrix, result, visited, distance + adjMatrix[u][v], N);
        }
    }

    visited[u] = false;
}

// Function to compute the detour matrix using an adjacency matrix
std::vector<std::vector<double>> computeDetourMatrixL(int N, const std::vector<std::vector<double>>& adjMatrix) {
    std::vector<std::vector<double>> detourMatrix(N, std::vector<double>(N, 0.0));

    // Compute longest path from each node
    for (int i = 0; i < N; ++i) {
        std::vector<double> result(N, 0.0);
        std::vector<bool> visited(N, false);
        longestSimplePathL(i, adjMatrix, result, visited, 0.0, N);

        for (int j = 0; j < N; ++j) {
            detourMatrix[i][j] = result[j];
        }
    }

    return detourMatrix;
}

// Compute Detour Index
int computeDetourIndexL(const std::vector<std::vector<double>>& detourMatrix) {
    double sum = 0.0;
    for (size_t i = 0; i < detourMatrix.size(); ++i) {
        for (size_t j = 0; j < detourMatrix[i].size(); ++j) {
            sum += detourMatrix[i][j];
        }
    }
    return static_cast<int>(0.5 * sum);
}



// Function to get the Adjacency Matrix using standard C++ structures
std::vector<std::vector<double>> calculateAdjacencyMatrixL(const RDKit::ROMol& mol) {
    unsigned int nAtoms = mol.getNumAtoms();

    // Initialize a 2D vector with zeros
    std::vector<std::vector<double>> adjMatrix(nAtoms, std::vector<double>(nAtoms, 0.0));

    // Populate the adjacency matrix using RDKit's bond information
    for (unsigned int i = 0; i < nAtoms; ++i) {
        for (unsigned int j = 0; j < nAtoms; ++j) {
            const RDKit::Bond* bond = mol.getBondBetweenAtoms(i, j);
            adjMatrix[i][j] = (bond != nullptr) ? 1.0 : 0.0;
        }
    }

    return adjMatrix;
}

std::vector<double> calcDetourMatrixDescsL(const RDKit::ROMol& mol) {
    auto adjMatrix = calculateAdjacencyMatrixL(mol); // Convert adjacency matrix computation to LAPACK
    int numAtoms = mol.getNumAtoms();

    auto detourMatrix = computeDetourMatrixL(numAtoms, adjMatrix);

    // Convert detourMatrix to 1D array in column-major order for LAPACK
    std::vector<double> flatMatrix(numAtoms * numAtoms);
    for (int i = 0; i < numAtoms; ++i) {
        for (int j = 0; j < numAtoms; ++j) {
            flatMatrix[j * numAtoms + i] = detourMatrix[i][j];
        }
    }

    // Compute eigenvalues and eigenvectors using LAPACK
    std::vector<double> eigenvalues(numAtoms);
    std::vector<std::vector<double>> eigenvectors(numAtoms, std::vector<double>(numAtoms));

    compute_eigenvalues_and_eigenvectorsL(detourMatrix, eigenvalues, eigenvectors);

    // Compute descriptors
    double Sp_Abs = spAbsL(eigenvalues);
    double Sp_Max = spMaxL(eigenvalues);
    double Sp_Diam = spDiamL(eigenvalues);
    double Sp_Mean = spMeanL(eigenvalues);
    double Sp_AD = spADL(eigenvalues, Sp_Mean);
    double Sp_MAD = Sp_AD / numAtoms;
    double Log_EE = logEE_stable(eigenvalues);

    std::vector<std::pair<int, int>> bonds;

    // Iterate over all bonds in the molecule
    for (const auto& bond : mol.bonds()) {
        int beginAtomIdx = bond->getBeginAtomIdx();
        int endAtomIdx = bond->getEndAtomIdx();
        bonds.push_back(std::make_pair(beginAtomIdx, endAtomIdx));
    }

    // Calculate descriptors
    double ve1 = VE1L(eigenvectors);
    double ve2 = VE2L(ve1, numAtoms);
    double ve3 = VE3L(ve1, numAtoms);
    double vr1 = VR1L(eigenvectors, bonds);
    double vr2 = VR2L(vr1, numAtoms);
    double vr3 = VR3L(vr1, numAtoms);
    double detourindex = static_cast<double>(computeDetourIndexL(detourMatrix));


    double sm1 = SM1L(detourMatrix);
    return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE, sm1, ve1,ve2,ve3,vr1,vr2,vr3,detourindex};

}


 //// AdjacencyMatrix


    std::vector<double> calcAdjMatrixDescs(const RDKit::ROMol& mol) {
        Eigen::MatrixXd AdjMat = calculateAdjacencyMatrix( mol);

        int numAtoms = mol.getNumAtoms();

        Eigen::VectorXd eigenvalues;

        Eigen::MatrixXd eigenvectors;

        compute_eigenvalues_and_eigenvectors(AdjMat, eigenvalues, eigenvectors);

        std::vector<std::pair<int, int>> bonds;

        // Iterate over all bonds in the molecule
        for (const auto& bond : mol.bonds()) {
            int beginAtomIdx = bond->getBeginAtomIdx();
            int endAtomIdx = bond->getEndAtomIdx();
            bonds.push_back(std::make_pair(beginAtomIdx, endAtomIdx));
        }
    /*
        // Compute descriptors
        std::cout << "eigenvalues from Eigen: ";
        for (auto ei : eigenvalues) {
            std::cout << ei << " ";
        }
        std::cout << std::endl;
    */
        double Sp_Abs = spAbs(eigenvalues);
        double Sp_Max = spMax(eigenvalues);
        double Sp_Diam = spDiam(eigenvalues);
        double Sp_Mean = spMean(eigenvalues); // tmp values not needed to export as result
        double Sp_AD = spAD(eigenvalues, Sp_Mean);
        double Sp_MAD = Sp_AD / numAtoms;
        double Log_EE = logEE(eigenvalues);  // this one is not correct ??? do we need to add the bond
        double ve1 = VE1(AdjMat, eigenvalues, eigenvectors);
        double ve2 = VE2(AdjMat, numAtoms, eigenvalues, eigenvectors);
        double ve3 = VE3(AdjMat, numAtoms, eigenvalues, eigenvectors) ;
        double vr1 = VR1(AdjMat, bonds, eigenvalues, eigenvectors);
        double vr2 = VR2(AdjMat, bonds, numAtoms, eigenvalues, eigenvectors) ;
        double vr3 = VR3(AdjMat, bonds, numAtoms, eigenvalues, eigenvectors);
        return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE,  ve1, ve2, ve3, vr1, vr2, vr3};

    }

    std::vector<double> calcAdjMatrixDescsL(const RDKit::ROMol& mol) {
        auto AdjMat = calculateAdjacencyMatrixL(mol);  // Use LAPACK-compatible adjacency matrix

        int numAtoms = mol.getNumAtoms();

        std::vector<double> eigenvalues;
        std::vector<std::vector<double>> eigenvectors;

        compute_eigenvalues_and_eigenvectorsL(AdjMat, eigenvalues, eigenvectors);  // LAPACK eigen computation

        std::vector<std::pair<int, int>> bonds;
        for (const auto& bond : mol.bonds()) {
            bonds.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
        }
/*
        // Compute descriptors
        std::cout << "eigenvalues from Lapack: ";
        for (auto ei : eigenvalues) {
            std::cout << ei << " ";
        }
        std::cout << std::endl;
*/
        // Compute descriptors
        double Sp_Abs = spAbsL(eigenvalues);
        double Sp_Max = spMaxL(eigenvalues);
        double Sp_Diam = spDiamL(eigenvalues);
        double Sp_Mean = spMeanL(eigenvalues);
        double Sp_AD = spADL(eigenvalues, Sp_Mean);
        double Sp_MAD = Sp_AD / numAtoms;
        double Log_EE = logEE_stable(eigenvalues);

        double ve1 = VE1L(eigenvectors);
        double ve2 = VE2L(ve1, numAtoms);
        double ve3 = VE3L(ve1, numAtoms);
        double vr1 = VR1L(eigenvectors, bonds);
        double vr2 = VR2L(vr1, numAtoms);
        double vr3 = VR3L(vr1, numAtoms);
        return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE,  ve1, ve2, ve3, vr1, vr2, vr3};
    }




    std::vector<double> calcDistMatrixDescs(const RDKit::ROMol& mol) {
        Eigen::MatrixXd DMat = calculateDistanceMatrix( mol);

        int numAtoms = mol.getNumAtoms();

        Eigen::VectorXd eigenvalues;

        Eigen::MatrixXd eigenvectors;

        compute_eigenvalues_and_eigenvectors(DMat, eigenvalues, eigenvectors);

        std::vector<std::pair<int, int>> bonds;

            // Iterate over all bonds in the molecule
            for (const auto& bond : mol.bonds()) {
                int beginAtomIdx = bond->getBeginAtomIdx();
                int endAtomIdx = bond->getEndAtomIdx();
                bonds.push_back(std::make_pair(beginAtomIdx, endAtomIdx));
            }

        // Compute descriptors
        double Sp_Abs = spAbs(eigenvalues);
        double Sp_Max = spMax(eigenvalues);
        double Sp_Diam = spDiam(eigenvalues);
        double Sp_Mean = spMean(eigenvalues); // tmp values not needed to export as result
        double Sp_AD = spAD(eigenvalues, Sp_Mean);
        double Sp_MAD = Sp_AD / numAtoms;
        double Log_EE = logEE(eigenvalues);  // this one is not correct ??? do we need to add the bond
        double ve1 = VE1(DMat, eigenvalues, eigenvectors);
        double ve2 = VE2(DMat, numAtoms, eigenvalues, eigenvectors);
        double ve3 = VE3(DMat, numAtoms, eigenvalues, eigenvectors) ;
        double vr1 = VR1(DMat, bonds, eigenvalues, eigenvectors);
        double vr2 = VR2(DMat, bonds, numAtoms, eigenvalues, eigenvectors) ;
        double vr3 = VR3(DMat, bonds, numAtoms, eigenvalues, eigenvectors);

        // SM1(DMat) intentionally omitted for parity with L variant
        return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE, ve1, ve2, ve3, vr1, vr2, vr3};

    }

std::vector<std::vector<double>> calculateDistanceMatrixL(const RDKit::ROMol& mol) {
    unsigned int nAtoms = mol.getNumAtoms();

    // Get the distance matrix using RDKit's MolOps::getDistanceMat
    double* distanceMat = MolOps::getDistanceMat(mol, false, false, false);  // No bond order, no weights, no hydrogens

    // Convert the raw pointer to a 2D vector
    std::vector<std::vector<double>> distMatrix(nAtoms, std::vector<double>(nAtoms, 0.0));
    for (unsigned int i = 0; i < nAtoms; ++i) {
        for (unsigned int j = 0; j < nAtoms; ++j) {
            distMatrix[i][j] = distanceMat[i * nAtoms + j];
        }
    }

    return distMatrix;
}

std::vector<double> calcDistMatrixDescsL(const RDKit::ROMol& mol) {
    auto DMat = calculateDistanceMatrixL(mol);  // Use LAPACK-compatible distance matrix

    int numAtoms = mol.getNumAtoms();

    std::vector<double> eigenvalues;
    std::vector<std::vector<double>> eigenvectors;

    compute_eigenvalues_and_eigenvectorsL(DMat, eigenvalues, eigenvectors);  // LAPACK eigen computation

    std::vector<std::pair<int, int>> bonds;
    for (const auto& bond : mol.bonds()) {
        bonds.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
    }

    // Compute descriptors
    double Sp_Abs = spAbsL(eigenvalues);
    double Sp_Max = spMaxL(eigenvalues);
    double Sp_Diam = spDiamL(eigenvalues);
    double Sp_Mean = spMeanL(eigenvalues);
    double Sp_AD = spADL(eigenvalues, Sp_Mean);
    double Sp_MAD = Sp_AD / numAtoms;
    double Log_EE = logEE_stable(eigenvalues);
    double ve1 = VE1L(eigenvectors);
    double ve2 = VE2L(ve1, numAtoms);
    double ve3 = VE3L(ve1, numAtoms);
    double vr1 = VR1L(eigenvectors, bonds);
    double vr2 = VR2L(vr1, numAtoms);
    double vr3 = VR3L(vr1, numAtoms);
    return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE, ve1, ve2, ve3, vr1, vr2, vr3};
}



    // Molecular Distance Edge Calculation
    std::vector<double> calcMolecularDistanceEdgeDescs(const RDKit::ROMol& mol) {
        using namespace RDKit;

        // Initialize result container
        std::vector<double> results;

        // Supported atomic numbers and valence pairs
        std::vector<int> atomic_nums = {6, 8, 7}; // Carbon, Oxygen, Nitrogen

        std::map<int, std::vector<std::pair<int, int>>> valence_pairs = {
            {6, {{1, 1}, {1, 2}, {1, 3}, {1, 4},
                {2, 2}, {2, 3}, {2, 4},
                {3, 3}, {3, 4}, {4, 4}}},
            {8, {{1, 1}, {1, 2}, {2, 2}}},
            {7, {{1, 1}, {1, 2}, {1, 3},
                {2, 2}, {2, 3}, {3, 3}}}};

        // Compute distance matrix
        const int num_atoms = mol.getNumAtoms();
        Eigen::MatrixXd distance_matrix = calculateDistanceMatrix(mol);
        std::vector<double> valences = calcValence(mol);

        // Iterate over atomic numbers
        for (const auto& atomic_num : atomic_nums) {
            const auto& pairs = valence_pairs[atomic_num];

            // Iterate over valence pairs
            for (const auto& [valence1, valence2] : pairs) {
                std::vector<double> Dv;

                // Filter distances based on valence and atomic number
                for (int i = 0; i < num_atoms; ++i) {
                    const auto atom_i = mol.getAtomWithIdx(i);
                    for (int j = i + 1; j < num_atoms; ++j) {
                        const auto atom_j = mol.getAtomWithIdx(j);
                        if ((valences[i] == valence1 && valences[j] == valence2) || (valences[j] == valence1 && valences[i] == valence2)) {
                            if (atom_j->getAtomicNum() == atomic_num && atom_i->getAtomicNum() == atomic_num) {
                                Dv.push_back(distance_matrix(i, j));
                            }
                        }
                    }
                }

                // Compute descriptor if valid distances are found
                int n = Dv.size();
                if (n > 0) {
                    double log_sum = 0.0;
                    for (const auto& d : Dv) {
                        log_sum += std::log(d);
                    }
                    double log_dx = log_sum / (2.0 * n);
                    double dx = std::exp(log_dx);
                    double descriptor_value = n / (dx * dx);

                    results.push_back(descriptor_value);
                } else {
                    results.push_back(0.0); // Default value when no pairs are found
                }
            }
        }

        return results;
    }

////////////////////////// ExtendedTopochemicalAtom  //////////////////////////


    // ExtendedTopochemicalAtom  all codes
    RDKit::RWMol* modifyMolecule(const RDKit::ROMol& mol, bool explicitHydrogens, bool saturated) {
        RDKit::RWMol* newMol = new RDKit::RWMol();
        std::map<unsigned int, unsigned int> atomMap;

        // Add atoms to the new molecule
        for (const auto& atom : mol.atoms()) {
            if (atom->getAtomicNum() == 1) continue; // Skip hydrogens

            RDKit::Atom* newAtom = nullptr;

            if (saturated) {
                // Preserve atomic number and formal charge
                newAtom = new RDKit::Atom(atom->getAtomicNum());
                newAtom->setFormalCharge(atom->getFormalCharge());
                newAtom->setIsAromatic(atom->getIsAromatic()); // Retain aromaticity
            } else {
                // Replace all non-hydrogen atoms with carbon
                newAtom = new RDKit::Atom(6); // Carbon
                newAtom->setIsAromatic(false); // Reset aromaticity
            }

            unsigned int newIdx = newMol->addAtom(newAtom, false, true);
            atomMap[atom->getIdx()] = newIdx;
        }

        // Add bonds to the new molecule
        for (const auto& bond : mol.bonds()) {
            unsigned int startIdx = bond->getBeginAtomIdx();
            unsigned int endIdx = bond->getEndAtomIdx();

            if (atomMap.find(startIdx) != atomMap.end() && atomMap.find(endIdx) != atomMap.end()) {
                RDKit::Bond::BondType bondType;

                if (saturated) {
                    // Retain original bond type if at least one atom is not carbon
                    const auto* startAtom = mol.getAtomWithIdx(startIdx);
                    const auto* endAtom = mol.getAtomWithIdx(endIdx);
                    if (startAtom->getAtomicNum() != 6 || endAtom->getAtomicNum() != 6) {
                        bondType = bond->getBondType();
                    } else {
                        bondType = RDKit::Bond::SINGLE; // Default to single bond
                    }
                } else {
                    bondType = RDKit::Bond::SINGLE; // Force single bonds for reference molecules
                }

                newMol->addBond(atomMap[startIdx], atomMap[endIdx], bondType);
            }
        }

        // Attempt sanitization
        unsigned int failedOps = 0;
        try {
            RDKit::MolOps::sanitizeMol(*newMol, failedOps, RDKit::MolOps::SANITIZE_ALL ^ RDKit::MolOps::SANITIZE_PROPERTIES);
        } catch (const RDKit::MolSanitizeException& e) {
            std::cerr << "Sanitization failed: " << e.what() << "\n";
        }

        if (failedOps > 0) {
            std::cerr << "Sanitization failed on some operations, continuing. Failed ops: " << failedOps << "\n";
        }

        // Add explicit hydrogens if requested
        if (explicitHydrogens) {
            RDKit::MolOps::addHs(*newMol);
        }

        // Ensure the molecule is Kekulized
        try {
            RDKit::MolOps::Kekulize(*newMol, false);
        } catch (const RDKit::KekulizeException& e) {
            std::cerr << "Kekulization failed: " << e.what() << "\n";
        }

        return newMol;
    }


    // Function to check if any heteroatom has more than 4 neighbors
    bool hasExcessiveNeighbors(const RDKit::ROMol& mol) {
        for (auto atom : mol.atoms()) {
            int total_neighbors = atom->getDegree() + atom->getTotalNumHs();
            if (atom->getAtomicNum() != 6 && total_neighbors > 4) {
                std::cout << "Skipping molecule due to heteroatom with more than 4 neighbors (Atom Index: "
                        << atom->getIdx() << " Atomic Num: " << atom->getAtomicNum() << " Neighbors: "
                        << total_neighbors << ")\n";
                return true; // A heteroatom with more than 4 neighbors found
            }
        }
        return false;  // Safe to proceed
    }


    // CAUTION : wrong for type 4 ie saturated True use modifyMolecule instead
    RDKit::RWMol* cloneAndModifyMolecule(const RDKit::ROMol& originalMol, bool explicitHydrogens, bool saturated) {
        try {

            // Check before proceeding (only meaningful when atoms keep their real
            // identity; short-circuit on `saturated` so the carbon-skeleton path
            // does not emit spurious warnings for these hypothetical objects).
            if (saturated && hasExcessiveNeighbors(originalMol)) {
                return nullptr;
            }


            // Clone the molecule
            RDKit::RWMol* clonedMol = new RDKit::RWMol(originalMol);

            // Iterate over atoms and modify them
            for (auto atom : clonedMol->atoms()) {
                if (atom->getAtomicNum() == 1) {
                    continue;
                }

                if (saturated) {
                    // Retain the atomic number but adjust formal charges if necessary ...
                    atom->setFormalCharge(atom->getFormalCharge());
                } else {
                    // Change heavy atoms (non-hydrogen) to Carbon (atomic number 6)
                    if (atom->getAtomicNum() != 6) {
                        atom->setAtomicNum(6);
                    }
                }
            }

            // Iterate over bonds and modify them
            for (auto bond : clonedMol->bonds()) {
                // Adjust bond types if not saturated
                if (!saturated ) {
                    bond->setBondType(RDKit::Bond::SINGLE);
                }
            }

            // v3 fix: this is an artificial graph (heavy atoms abstracted to
            // carbon, bonds forced to single), not a real molecule, so a
            // carried-over hypervalent atom (e.g. pentavalent P -> C5) must NOT
            // raise a valence error. Compute valences non-strictly and skip the
            // strict valence check (SANITIZE_PROPERTIES).
            unsigned int failedOps = 0;
            clonedMol->updatePropertyCache(false);
            RDKit::MolOps::sanitizeMol(*clonedMol, failedOps,
                                       RDKit::MolOps::SANITIZE_ALL ^ RDKit::MolOps::SANITIZE_PROPERTIES);

            if (failedOps > 0) {
                std::cerr << "Sanitization failed on properties, but continuing. Failed ops: " << failedOps << "\n";
            }

            if (explicitHydrogens) {
                RDKit::MolOps::addHs(*clonedMol);
            }

            // Ensure the molecule is Kekulized
            RDKit::MolOps::Kekulize(*clonedMol, false); // same error as before !!! keep aromatic flags please!!!

            return clonedMol; // Return the modified clone
        } catch (const RDKit::MolSanitizeException& e) {
            std::cerr << "Sanitization failed: " << e.what() << "\n";
            return nullptr; // Return nullptr if the modification fails
        } catch (const std::exception& e) {
            std::cerr << "Error cloning and modifying molecule: " << e.what() << "\n";
            return nullptr;
        }
    }


    // Helper function to get the core count
    double getCoreCount(const Atom& atom) {
        int Z = atom.getAtomicNum();
        if (Z == 1) {
            return 0.0;
        }

        const RDKit::PeriodicTable* tbl = RDKit::PeriodicTable::getTable();
        double Zv = tbl->getNouterElecs(Z);
        int PN = GetPrincipalQuantumNumber(Z);

        return (Z - Zv) / (Zv * (PN - 1.));
    }

    // Helper function to calculate eta_epsilon
    double getEtaEpsilon(const Atom& atom) {
        const RDKit::PeriodicTable* tbl = RDKit::PeriodicTable::getTable();
        double Zv = tbl->getNouterElecs(atom.getAtomicNum());
        return 0.3 * Zv - getCoreCount(atom);
    }

    // Helper function to calculate eta_beta_sigma
    double getEtaBetaSigma(const Atom& atom) {
        double e = getEtaEpsilon(atom);
        double betaSigma = 0.0;

        for (const auto& neighbor : atom.getOwningMol().atomNeighbors(&atom)) {
            if (neighbor->getAtomicNum() != 1) {
                double neighborEpsilon = getEtaEpsilon(*neighbor);
                betaSigma += (std::abs(neighborEpsilon - e) <= 0.3) ? 0.5 : 0.75;
            }
        }
        return betaSigma;
    }

    // Helper function to calculate non-sigma contribution
    double getEtaNonSigmaContribute(const Bond& bond) {
        if (bond.getBondType() == Bond::BondType::SINGLE) {
            return 0.0;
        }
        double f = 1.0;
        if (bond.getBondType() == Bond::BondType::TRIPLE) {
            f = 2.0;
        }
        const Atom* a = bond.getBeginAtom();
        const Atom* b = bond.getEndAtom();
        double dEps = std::abs(getEtaEpsilon(*a) - getEtaEpsilon(*b));
        double y = 1.0;
        if (bond.getIsAromatic()) {
            y = 2.0;
        } else if (dEps > 0.3) {
            y = 1.5;
        }
        return y * f;
    }

    bool isAtomInRing(const RDKit::Atom& atom) {
        const RDKit::RingInfo* ringInfo = atom.getOwningMol().getRingInfo();
        if (ringInfo && ringInfo->isInitialized()) {
            return ringInfo->numAtomRings(atom.getIdx()) > 0;
        }
        return false; // Not in a ring if no ring info is available
    }

    // Helper function to calculate eta_beta_delta
    double getEtaBetaDelta(const RDKit::Atom& atom) {
        const RDKit::PeriodicTable* tbl = RDKit::PeriodicTable::getTable();
        if (atom.getIsAromatic() || isAtomInRing(atom) || (tbl->getNouterElecs(atom.getAtomicNum()) - atom.getTotalValence() <= 0)) {
            return 0.0;
        }

        for (const auto& neighbor : atom.getOwningMol().atomNeighbors(&atom)) {
            if (neighbor->getIsAromatic()) {
                return 0.5;
            }
        }
        return 0.0;
    }

    // Helper function to get the other atom in a bond
    const RDKit::Atom* getOtherAtom(const RDKit::Bond& bond, const RDKit::Atom& atom) {
        if (bond.getBeginAtom() == &atom) {
            return bond.getEndAtom();
        } else {
            return bond.getBeginAtom();
        }
    }

    // Helper function to calculate eta_beta_non_sigma
    double getEtaBetaNonSigma(const RDKit::Atom& atom) {
        double betaNonSigma = 0.0;

        for (const auto& bond : atom.getOwningMol().atomBonds(&atom)) {
            const RDKit::Atom* otherAtom = getOtherAtom(*bond, atom);
            if (otherAtom->getAtomicNum() != 1) {
                betaNonSigma += getEtaNonSigmaContribute(*bond);
            }
        }
        return betaNonSigma;
    }

    // Helper function to calculate eta_gamma
    double getEtaGamma(const RDKit::Atom& atom) {
        double beta = getEtaBetaSigma(atom) + getEtaBetaNonSigma(atom) + getEtaBetaDelta(atom);
        if (beta == 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return getCoreCount(atom) / beta;
    }

    //  Etacorecount of reference can be merge with next function with an additional parameter if needed...
    double calculateEtaCoreCountRef(const RDKit::ROMol& mol, bool averaged) {
        const RDKit::ROMol* targetMol = &mol; // Default to the input molecule
        //RDKit::ROMol* molWithHs = nullptr;
        targetMol = cloneAndModifyMolecule(mol, false, false); // this is false, false based on the python source code

        // Handle potential failure in cloning
        if (!targetMol) {
            // Suppress noisy warnings: return NaN silently
            delete targetMol;
            return std::numeric_limits<double>::quiet_NaN();  // Return NaN instead of crashing
        }

        double coreCount = 0.0;
        for (const auto& atom : targetMol->atoms()) {
            coreCount += getCoreCount(*atom);
        }
        if (averaged) {
            coreCount /= mol.getNumHeavyAtoms();
        }
        delete targetMol;
        return coreCount;
    }

    //  Etacorecount
    std::vector<double> calculateEtaCoreCount(const RDKit::ROMol& mol) {
        double coreCount = 0.0;
        for (const auto& atom : mol.atoms()) {
            coreCount += getCoreCount(*atom);
        }

        return {coreCount, coreCount / mol.getNumHeavyAtoms()};
    }

    std::vector<double> calculateEtaShapeIndex(const RDKit::ROMol& mol, double alpha) {
        double shapeAlphaP = 0.0;
        double shapeAlphaY = 0.0;
        double shapeAlphaX = 0.0;

        for (const auto& atom : mol.atoms()) {
            switch (atom->getDegree()) {
                case 1:
                    shapeAlphaP += getCoreCount(*atom);
                    break;
                case 3:
                    shapeAlphaY += getCoreCount(*atom);
                    break;
                case 4:
                    shapeAlphaX += getCoreCount(*atom);
                    break;
                default:
                    break;  // Ignore atoms with other degrees
            }
        }

        return {shapeAlphaP / alpha, shapeAlphaY / alpha, shapeAlphaX / alpha};
    }

    std::vector<double> computeEtaBetaDescriptors(const RDKit::Atom& atom) {
        const RDKit::PeriodicTable* tbl = RDKit::PeriodicTable::getTable();

        double epsilon = getEtaEpsilon(atom);
        double betaSigma = 0.0;
        double betaNonSigma = 0.0;
        double betaDelta = 0.0;

        bool isAromatic = atom.getIsAromatic();
        bool inRing = isAtomInRing(atom);

        // Check if the atom satisfies delta condition
        if (!isAromatic && !inRing && (tbl->getNouterElecs(atom.getAtomicNum()) - atom.getTotalValence() > 0)) {
            for (const auto& neighbor : atom.getOwningMol().atomNeighbors(&atom)) {
                if (neighbor->getIsAromatic()) {
                    betaDelta = 0.5;
                    break;
                }
            }
        }

        // Iterate through bonds
        for (const auto& bond : atom.getOwningMol().atomBonds(&atom)) {
            const RDKit::Atom* neighbor = getOtherAtom(*bond, atom);
            if (neighbor->getAtomicNum() != 1) {
                double neighborEpsilon = getEtaEpsilon(*neighbor);

                // Compute betaSigma (sigma contribution)
                betaSigma += (std::abs(neighborEpsilon - epsilon) <= 0.3) ? 0.5 : 0.75;

                // Compute betaNonSigma (non-sigma contribution)
                if (bond->getBondType() != RDKit::Bond::SINGLE) {
                    double bondFactor = (bond->getBondType() == RDKit::Bond::TRIPLE) ? 2.0 : 1.0;
                    double dEps = std::abs(epsilon - neighborEpsilon);
                    double bondContribution = bond->getIsAromatic() ? 2.0 : (dEps > 0.3 ? 1.5 : 1.0);
                    betaNonSigma += bondContribution * bondFactor;
                }
            }
        }

        return {betaSigma, betaNonSigma, betaDelta};
    }

    std::vector<double> calculateEtaVEMCount(const RDKit::ROMol& mol) {
        double beta = 0.0;
        double beta_s = 0.0;
        double beta_ns = 0.0;
        double beta_ns_d = 0.0;

        int numAtoms = mol.getNumAtoms();

        for (const auto& atom : mol.atoms()) {

            std::vector<double> EBD = computeEtaBetaDescriptors(*atom);
            double betaDelta = EBD[2];
            double betaSigma = EBD[0] / 2.0;
            double betaNonSigma = EBD[1] / 2.0 + betaDelta;

            beta_s += betaSigma;
            beta_ns += betaNonSigma;
            beta_ns_d += betaDelta;
            beta += betaSigma + betaNonSigma;
        }

        // Calculate averaged values
        double avg_beta = beta / numAtoms;
        double avg_beta_s = beta_s / numAtoms;
        double avg_beta_ns = beta_ns / numAtoms;
        double avg_beta_ns_d = beta_ns_d / numAtoms;

        return {
            beta,          // ETA_beta
            avg_beta,      // AETA_beta
            beta_s,        // ETA_beta_s
            avg_beta_s,    // AETA_beta_s
            beta_ns,       // ETA_beta_ns
            avg_beta_ns,   // AETA_beta_ns
            beta_ns_d,     // ETA_beta_ns_d
            avg_beta_ns_d  // AETA_beta_ns_d
        };
    }



    // Function to Calculate ETA Composite Index
    double calculateEtaCompositeIndex(const RDKit::ROMol& mol, bool useReference, bool local, bool averaged) {
        // Fetch molecule (reference or input)

        const RDKit::ROMol* targetMol = &mol; // Default to the input molecule
        //RDKit::ROMol* molWithHs = nullptr;

        if (useReference) {
            targetMol = cloneAndModifyMolecule(mol, false, false); // this is false, false based on the python source code

            // Handle potential failure in cloning
            if (!targetMol) {
                // Suppress noisy warnings: return NaN silently
                delete targetMol;
                return std::numeric_limits<double>::quiet_NaN();  // Return NaN instead of crashing
            }
        }



        // Calculate distance matrix

        Eigen::MatrixXd distanceMatrix = calculateDistanceMatrix(*targetMol);
        int numAtoms = targetMol->getNumAtoms();


        // Define gamma values for each atom
        std::vector<double> gamma(numAtoms, 0.0);
        for (const auto& atom : targetMol->atoms()) {
            gamma[atom->getIdx()] = getEtaGamma(*atom);
        }

        // ETA calculation "triangle computation" ie  j = i + 1 trick to go faster
        double eta = 0.0;
        for (int i = 0; i < numAtoms; ++i) {
            for (int j = i + 1; j < numAtoms; ++j) {
                if (local && distanceMatrix(i,j) != 1.0) continue;
                if (!local && distanceMatrix(i,j) == 0.0) continue;
                eta += std::sqrt(gamma[i] * gamma[j] / (distanceMatrix(i,j) * distanceMatrix(i,j)));
            }
        }

        // Averaged value if needed
        if (averaged) {
            eta /= numAtoms;
        }
        if (useReference) {
            delete  targetMol;
        }
        return eta;
    }


    std::vector<double> calculateEtaCompositeIndices(const RDKit::ROMol& mol) {
        std::vector<double> etaValues(8, 0.0);

        // Define options for different cases
        std::vector<std::tuple<bool, bool, bool>> options = {
            {false, false, false}, // ETA_eta
            {false, false, true},  // AETA_eta
            {false, true, false},  // ETA_eta_L
            {false, true, true},   // AETA_eta_L
            {true, false, false},  // ETA_eta_R
            {true, false, true},   // AETA_eta_R
            {true, true, false},   // ETA_eta_RL
            {true, true, true}     // AETA_eta_RL
        };

        for (size_t idx = 0; idx < options.size(); ++idx) {
            bool useReference = std::get<0>(options[idx]);
            bool local = std::get<1>(options[idx]);
            bool averaged = std::get<2>(options[idx]);

            // Fetch molecule (reference or input)
            const RDKit::ROMol* targetMol = &mol;  // Default to the input molecule
            if (useReference) {
                targetMol = cloneAndModifyMolecule(mol, false, false);  // Clone with modifications (why not "True")

                // Handle potential failure in cloning
                if (!targetMol) {
                    // Suppress noisy warnings: return NaN silently
                    delete targetMol;
                    return std::vector<double>(8, std::numeric_limits<double>::quiet_NaN());  // Return vector with NaN
                }

            }

            // Calculate distance matrix
            Eigen::MatrixXd distanceMatrix = calculateDistanceMatrix(*targetMol);
            int numAtoms = targetMol->getNumAtoms();

            // Define gamma values for each atom
            std::vector<double> gamma(numAtoms, 0.0);
            for (const auto& atom : targetMol->atoms()) {
                gamma[atom->getIdx()] = getEtaGamma(*atom);
            }

            // ETA calculation using the optimized "triangle computation"
            double eta = 0.0;
            for (int i = 0; i < numAtoms; ++i) {
                for (int j = i + 1; j < numAtoms; ++j) {
                    if (local && distanceMatrix(i, j) != 1.0) continue;
                    if (!local && distanceMatrix(i, j) == 0.0) continue;

                    eta += std::sqrt(gamma[i] * gamma[j] / (distanceMatrix(i, j) * distanceMatrix(i, j)));
                }
            }

            // Averaged value if required
            if (averaged) {
                eta /= numAtoms;
            }

            etaValues[idx] = eta;

            if (useReference) {
                delete targetMol;  // Clean up dynamically allocated molecule
            }
        }

        return etaValues;
    }


    std::vector<double> calculateEtaFunctionalityIndices(const RDKit::ROMol& mol) {
        std::vector<double> etaFunctionalityValues(4, 0.0);

        // Define options for different cases
        std::vector<std::tuple<bool, bool>> options = {
            {false, false}, // ETA_eta_F
            {false, true},  // AETA_eta_F
            {true, false},  // ETA_eta_FL
            {true, true}    // AETA_eta_FL
        };

        for (size_t idx = 0; idx < options.size(); ++idx) {
            bool local = std::get<0>(options[idx]);
            bool averaged = std::get<1>(options[idx]);

            // Calculate eta without reference and with reference
            double eta = calculateEtaCompositeIndex(mol, false, local, false);
            double etaRef = calculateEtaCompositeIndex(mol, true, local, false);

            // Compute functionality index
            double etaF = etaRef - eta;

            // Apply averaging if needed
            if (averaged) {
                etaF /= mol.getNumAtoms();
            }

            etaFunctionalityValues[idx] = etaF;
        }

        return etaFunctionalityValues;
    }


    std::vector<double> calculateEtaBranchingIndices(const RDKit::ROMol& mol) {
        std::vector<double> etaBranchingValues(4, 0.0);
        int atomCount = mol.getNumAtoms();

        if (atomCount <= 1) {
            return etaBranchingValues;  // Return zeros if the molecule has only one atom
        }

        // Calculate non-local branching term
        double eta_NL = (atomCount == 2) ? 1.0 : (std::sqrt(2.0) + 0.5 * (atomCount - 3));

        double eta_RL = calculateEtaCompositeIndex(mol, true, true, false);


        // Calculate ring count once
        double ringCount = RDKit::Descriptors::calcNumRings(mol);

        // Define options for different cases
        std::vector<std::tuple<bool, bool>> options = {
            {false, false}, // ETA_eta_B
            {false, true},  // AETA_eta_B
            {true, false},  // ETA_eta_BR
            {true, true}    // AETA_eta_BR
        };

        for (size_t idx = 0; idx < options.size(); ++idx) {
            bool useRingCount = std::get<0>(options[idx]);
            bool averaged = std::get<1>(options[idx]);

            double etaBranch = eta_NL - eta_RL;
            if (useRingCount) {
                etaBranch += 0.086 * ringCount;
            }

            if (averaged) {
                etaBranch /= atomCount;
            }

            etaBranchingValues[idx] = etaBranch;
        }

        return etaBranchingValues;
    }

    std::vector<double> calculateEtaDeltaBetaAll(const RDKit::ROMol& mol, double beta_ns, double beta_s) {
        std::vector<double> results(2, 0.0);

        // Calculate ETA_beta (delta_beta without averaging)
        double delta_beta = beta_ns - beta_s;
        results[0] = delta_beta;  // ETA_beta

        // Calculate AETA_beta (averaged delta_beta)
        if (mol.getNumAtoms() > 0) {
            results[1] = delta_beta / mol.getNumAtoms();  // AETA_beta
        } else {
            results[1] = 0.0;  // Handle division by zero case
        }

        return results;
    }


    // Function to calculate ETA Psi
    double calculateEtaPsi(const RDKit::ROMol& mol, double alpha, double epsilon) {
        if (epsilon == 0) {
            throw std::runtime_error("Epsilon cannot be zero");
        }
        return alpha / (mol.getNumAtoms() * epsilon);
    }


    std::vector<double> calculateEtaDeltaPsiAll(const RDKit::ROMol&, double psi) {
        std::vector<double> results(2, 0.0);

        // Constants for reference values
        const double L = 0.714;
        const double R = psi;

        // Calculate ETA_dPsi_A (L - R)
        results[0] = std::max(L - R, 0.0);  // ETA_dPsi_A

        // Calculate ETA_dPsi_B (R - L)
        results[1] = std::max(R - L, 0.0);  // ETA_dPsi_B

        return results;
    }


    std::vector<double> calculateEtaDeltaAlpha(const RDKit::ROMol& mol, double alpha, double alpha_R) {
        std::vector<double> etaDeltaAlphaValues(2, 0.0);
        int numAtoms = mol.getNumAtoms();

        if (numAtoms <= 0) {
            return etaDeltaAlphaValues;  // Return zeros if the molecule has no atoms
        }

        // Calculate dAlpha_A and dAlpha_B
        etaDeltaAlphaValues[0] = std::max((alpha - alpha_R) / numAtoms, 0.0);  // ETA_dAlpha_A
        etaDeltaAlphaValues[1] = std::max((alpha_R - alpha) / numAtoms, 0.0);  // ETA_dAlpha_B

        return etaDeltaAlphaValues;
    }

std::vector<double> calculateEtaEpsilonAll(const RDKit::ROMol& mol) {
    std::vector<double> etaEps(5, 0.0);

    // Step 1: Calculate ETA epsilon for type 2 (non-hydrogen atoms)
    double epsilon_sum_type2 = 0.0;
    for (const auto& atom : mol.atoms()) {
        epsilon_sum_type2 += getEtaEpsilon(*atom);
    }
    etaEps[1] = epsilon_sum_type2 / mol.getNumAtoms();  // ETA_epsilon_2

    // Step 2: Add hydrogens for types 1 and 5 calculations
    std::unique_ptr<RDKit::ROMol> hmol(RDKit::MolOps::addHs(mol));


    // Step 3: Calculate ETA epsilon for type 1 and type 5 (hydrogen inclusion)
    double epsilon_sum_type1 = 0.0, epsilon_sum_type5 = 0.0;
    int count_type1 = 0, count_type5 = 0;
    // Debug counters
    int debug_num_h_in_hmol = 0;
    int debug_excluded_ch = 0;

    for (const auto& atom : hmol->atoms()) {
        // type 1
        double EtaEps =  getEtaEpsilon(*atom);
        epsilon_sum_type1 += EtaEps;
        count_type1++;
        // type 5
        bool hasCarbonNeighbor = false;
        for (const auto& neighbor : hmol->atomNeighbors(atom)) {
            if (neighbor->getAtomicNum() == 6) {
                hasCarbonNeighbor = true;
                break;
            }
        }
        if (atom->getAtomicNum() == 1) {
            // count hydrogens in hmol
            debug_num_h_in_hmol++;
        }
        if (atom->getAtomicNum() != 1 || !hasCarbonNeighbor) {
            epsilon_sum_type5 += EtaEps;
            count_type5++;
        } else {
            // Hydrogen attached to carbon is excluded from epsilon_5
            if (atom->getAtomicNum() == 1 && hasCarbonNeighbor) {
                debug_excluded_ch++;
            }
        }
    }
    etaEps[0] = epsilon_sum_type1 / count_type1;  // ETA_epsilon_1
    etaEps[4] = count_type5 > 0 ? epsilon_sum_type5 / count_type5 : 0.0;  // ETA_epsilon_5

    // Optional debug logging (enabled when OSMO_ETA_DEBUG is set)
    const char* eta_dbg = std::getenv("OSMO_ETA_DEBUG");
    if (eta_dbg) {
        std::cerr << "[ETA_DEBUG] hmol_atoms=" << hmol->getNumAtoms()
                  << " hydrogens_in_hmol=" << debug_num_h_in_hmol
                  << " excluded_CH_for_eps5=" << debug_excluded_ch
                  << " epsilon_sum_type1=" << epsilon_sum_type1
                  << " epsilon_sum_type2=" << epsilon_sum_type2
                  << " epsilon_sum_type5=" << epsilon_sum_type5
                  << " count_type1=" << count_type1
                  << " count_type5=" << count_type5
                  << std::endl;
    }

    // Step 4: Calculate ETA epsilon for types 3 and 4 (modified molecules)
    const RDKit::ROMol* targetMol3 = modifyMolecule(mol, true, false);
    const RDKit::ROMol* targetMol4 = modifyMolecule(mol, true, true);

    double epsilon_sum_type3 = 0.0;
    for (const auto& atom : targetMol3->atoms()) {
        epsilon_sum_type3 += getEtaEpsilon(*atom);
    }
    etaEps[2] = epsilon_sum_type3 / targetMol3->getNumAtoms();  // ETA_epsilon_3

    double epsilon_sum_type4 = 0.0;
    for (const auto& atom : targetMol4->atoms()) {
        epsilon_sum_type4 += getEtaEpsilon(*atom);
    }
    etaEps[3] = epsilon_sum_type4 / targetMol4->getNumAtoms();  // ETA_epsilon_4

    // Clean up dynamically allocated molecules
    delete targetMol3;
    delete targetMol4;

    return etaEps;
}


    std::vector<double> calcExtendedTopochemicalAtom(const RDKit::ROMol& mol) {
        std::vector<double> results;


        std::unique_ptr<RDKit::RWMol> kekulizedMol(new RDKit::RWMol(mol));
        RDKit::MolOps::Kekulize(*kekulizedMol, false);

        //  Eta alpha  Descriptors correct  "EtaCoreCount" : alpha & shape can be group in one function
        std::vector<double> myalphas = calculateEtaCoreCount(*kekulizedMol);
        double alpha = myalphas[0];

        results.push_back(myalphas[0]);  // Eta alpha
        results.push_back(myalphas[1]);  // average Eta alpha
        // Shape Descriptors correct     "EtaShapeIndex",

        // working for all atoms types
        std::vector<double> ESP = calculateEtaShapeIndex(*kekulizedMol, myalphas[0]);

        results.push_back(ESP[0]);  // ETA_shape_p
        results.push_back(ESP[1]);  // ETA_shape_y
        results.push_back(ESP[2]);  // ETA_shape_x

        // Beta Descriptors ie EtaVEMCount correct : can be group in one function for optimizatoin
        std::vector<double> EtaVEM = calculateEtaVEMCount(*kekulizedMol);

        results.push_back(EtaVEM[0]);  // ETA_beta
        results.push_back(EtaVEM[1]);  // AETA_beta
        double etabetas = EtaVEM[2];
        results.push_back(EtaVEM[2]);  // ETA_beta_s
        results.push_back(EtaVEM[3]);  // AETA_beta_s
        double etabetans = EtaVEM[4];
        results.push_back(EtaVEM[4]);  // ETA_beta_ns
        results.push_back(EtaVEM[5]);  // AETA_beta_ns
        results.push_back(EtaVEM[6]);  // ETA_beta_ns_d
        results.push_back(EtaVEM[7]);  // AETA_beta_ns_d

        // ETA Descriptors   ==  "EtaCompositeIndex"
        std::vector<double> ECI = calculateEtaCompositeIndices(*kekulizedMol);
        results.push_back(ECI[0]); // ETA_eta
        results.push_back(ECI[1]); // AETA_eta
        results.push_back(ECI[2]); // ETA_eta_L
        results.push_back(ECI[3]); // AETA_eta_L
        results.push_back(ECI[4]); // ETA_eta_R
        results.push_back(ECI[5]); // AETA_eta_R
        results.push_back(ECI[6]); // ETA_eta_RL
        results.push_back(ECI[7]); // AETA_eta_RL

       // Functionality and Branching EtaFunctionalityIndex not working for heteroatom molecule...
        std::vector<double> EFI = calculateEtaFunctionalityIndices(*kekulizedMol);
        results.push_back(EFI[0]); // ETA_eta_F
        results.push_back(EFI[1]); // AETA_eta_F
        results.push_back(EFI[2]); // ETA_eta_FL
        results.push_back(EFI[3]); // AETA_eta_FL

        // EtaBranchingIndex  working
        std::vector<double> EBI = calculateEtaBranchingIndices(*kekulizedMol);
        results.push_back(EBI[0]); // ETA_eta_B
        results.push_back(EBI[1]); // AETA_eta_B
        results.push_back(EBI[2]); // ETA_eta_BR
        results.push_back(EBI[3]); // AETA_eta_BR

        //"EtaDeltaAlpha" :  dAlpha_A, dAlpha_B  correct
        double alpha_R = calculateEtaCoreCountRef(*kekulizedMol, false);
        std::vector<double> EDA = calculateEtaDeltaAlpha(*kekulizedMol,alpha, alpha_R);
        results.push_back(EDA[0]);
        results.push_back(EDA[1]);

        // "EtaEpsilon":
        std::vector<double> EtaEps = calculateEtaEpsilonAll(*kekulizedMol);
        results.push_back(EtaEps[0]);
        results.push_back(EtaEps[1]);
        results.push_back(EtaEps[2]);
        results.push_back(EtaEps[3]);
        results.push_back(EtaEps[4]);

        // "EtaDeltaEpsilon", A,B,C,D use EtaEps diff cases
        results.push_back(EtaEps[1-1]-EtaEps[3-1]);
        results.push_back(EtaEps[1-1]-EtaEps[4-1]);
        results.push_back(EtaEps[3-1]-EtaEps[4-1]);
        results.push_back(EtaEps[2-1]-EtaEps[5-1]);

        //"EtaDeltaBeta":
        std::vector<double> EDBA = calculateEtaDeltaBetaAll(*kekulizedMol, etabetans, etabetas );
        results.push_back(EDBA[0]);
        results.push_back(EDBA[1]);

        // EtaPsi:
        double EtaPsi = calculateEtaPsi(*kekulizedMol, alpha, EtaEps[1]);
        results.push_back(EtaPsi);                        // ETA_psi_1

       //EtaDeltaPsi:
        std::vector<double> EDPA = calculateEtaDeltaPsiAll(*kekulizedMol, EtaPsi);
        results.push_back(EDPA[0]);
        results.push_back(EDPA[1]);


        return results;
    }

    // Fast aggregate: compute all descriptors in C++ in one pass
    std::vector<double> calcOsmordred(const RDKit::ROMol& mol) {
        // Silence RDKit warnings locally
        RDLog::LogStateSetter guard;

        std::vector<double> out;
        out.reserve(4200); // approximate capacity to avoid many reallocations

        // Always compute the full v2 set; keep signature for ABI stability
        const bool doExEstate = true;

        // Precompute/cached intermediates where safe
        // Note: We do not change algorithms; just reuse intermediates across calls
        std::unique_ptr<RDKit::RWMol> kekulizedMol(new RDKit::RWMol(mol));
        try {
            RDKit::MolOps::Kekulize(*kekulizedMol, false);
        } catch (...) {
            // leave kekulizedMol as-is when kekulization fails
        }

        // No shared matrix caching in baseline fast path

        // Some families already build needed matrices internally; where public
        // helpers exist (e.g., Adj/Dist matrices), we call the descriptor that
        // accepts version flags so we do not duplicate logic.

        // Collect results with minimal inserts
        auto append = [&out](const std::vector<double>& v){ out.insert(out.end(), v.begin(), v.end()); };
        auto appendInt = [&out](const std::vector<int>& v){ out.reserve(out.size()+v.size()); for(int x: v) out.push_back(static_cast<double>(x)); };

        append(calcABCIndex(mol));
        appendInt(calcAcidBase(mol));
        append(calcAdjMatrixDescsL(mol));
        appendInt(calcAromatic(mol));
        appendInt(calcAtomCounts(mol));
        append(calcAutoCorrelation(mol));
        append(calcBCUTs(mol));
        append(calcBalabanJ(mol));
        // Use L variant (faster), to match Python CalcBaryszMatrix
        append(calcBaryszMatrixDescsL(mol));
        append(calcBertzCT(mol));
        appendInt(calcBondCounts(mol));
        append(calcRNCG_RPCG(mol));
        append(calcCarbonTypes(mol));
        append(calcAllChiDescriptors(mol));
        append(calcConstitutional(mol));
        append(calcDetourMatrixDescsL(mol));
        append(calcDistMatrixDescsL(mol));
        append(calcEStateDescs(mol, doExEstate));
        append(calcEccentricConnectivityIndex(mol));
        append(calcExtendedTopochemicalAtom(mol));
        append(calcFragmentComplexity(mol));
        append(calcFramework(mol));
        append(calcHydrogenBond(mol));
        append(calcLogS(mol));
        append(calcInformationContent(mol, 5));
        append(calcKappaShapeIndex(mol));
        appendInt(calcLipinskiGhose(mol));
        append(calcMcGowanVolume(mol));
        append(calcMoeType(mol));
        append(calcMolecularDistanceEdgeDescs(mol));
        append(calcMolecularId(mol));
        append(calcPathCount(mol));
        append(calcPolarizability(mol));
        appendInt(calcRingCount(mol));
        append(calcRotatableBond(mol));
        append(calcSLogP(mol));
        append(calcTopoPSA(mol));
        append(calcTopologicalChargeDescs(mol));
        append(calcTopologicalIndex(mol));
        append(calcVdwVolumeABC(mol));
        append(calcVertexAdjacencyInformation(mol));
        append(calcWalkCounts(mol));
        append(calcWeight(mol));
        appendInt(calcWienerIndex(mol));
        append(calcZagrebIndex(mol));
        append(calcPol(mol));
        append(calcMR(mol));
        append(calcFlexibility(mol));
        append(calcSchultz(mol));
        append(calcAlphaKappaShapeIndex(mol));
        append(calcHEStateDescs(mol));
        append(calcBEStateDescs(mol));
        append(calcAbrahams(mol));
        append(calcANMat(mol));
        append(calcASMat(mol));
        append(calcAZMat(mol));
        append(calcDSMat(mol));
        append(calcDN2Mat(mol));
        append(calcFrags(mol));
        append(calcAddFeatures(mol));

        return out;
    }

    // v2.0: Timeout constant for Osmordred computation (60 seconds = 1 minute)
    constexpr int OSMORDRED_TIMEOUT_SECONDS = 60;
    
    // v2.0: Single molecule with timeout protection (all-or-nothing)
    // Returns NaN vector if computation exceeds timeout_seconds
    RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcOsmordredWithTimeout(
        const RDKit::ROMol& mol, int timeout_seconds) {
        
        auto future = std::async(std::launch::async, [&mol]() {
            return calcOsmordred(mol);
        });
        
        int actual_timeout = timeout_seconds > 0 ? timeout_seconds : OSMORDRED_TIMEOUT_SECONDS;
        auto status = future.wait_for(std::chrono::seconds(actual_timeout));
        
        if (status == std::future_status::ready) {
            return future.get();
        } else {
            // Timeout - return NaN vector (3585 NaN values)
            return std::vector<double>(3585, std::numeric_limits<double>::quiet_NaN());
        }
    }
    
    // v2.0: Batch version from SMILES: parses each SMILES -> NEW mol (tautomer canonical LOST).
    // For tautomer-canonical mols use calcOsmordredBatchFromMols(mols) with mols from ToBinary.
    RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>> calcOsmordredBatch(
        const std::vector<std::string>& smiles_list, int n_jobs) {

        std::vector<std::vector<double>> results;
        results.reserve(smiles_list.size());

        unsigned int nThreads = getNumThreadsToUse(n_jobs);

        if (nThreads <= 1 || smiles_list.size() < 10) {
            for (const auto& smi : smiles_list) {
                auto future = std::async(std::launch::async, [&smi]() {
                    ROMol* mol = SmilesToMol(smi);
                    if (mol) {
                        auto desc = calcOsmordred(*mol);
                        delete mol;
                        return desc;
                    }
                    return std::vector<double>();
                });
                
                auto status = future.wait_for(std::chrono::seconds(OSMORDRED_TIMEOUT_SECONDS));
                if (status == std::future_status::ready) {
                    results.push_back(future.get());
                } else {
                    results.push_back(std::vector<double>(3585, std::numeric_limits<double>::quiet_NaN()));
                }
            }
            return results;
        }
        
        // Parallel processing using std::async with timeout
        std::vector<std::future<std::vector<double>>> futures;
        futures.reserve(smiles_list.size());
        
        for (size_t idx = 0; idx < smiles_list.size(); ++idx) {
            const auto& smi = smiles_list[idx];
            
            futures.emplace_back(std::async(std::launch::async, [smi]() {
                try {
                    ROMol* mol = SmilesToMol(smi);
                    if (mol) {
                        try {
                            std::vector<double> descriptors = calcOsmordred(*mol);
                            delete mol;
                            return descriptors;
                        } catch (...) {
                            delete mol;
                            return std::vector<double>();
                        }
                    } else {
                        return std::vector<double>();
                    }
                } catch (...) {
                    return std::vector<double>();
                }
            }));
        }
        
        for (auto& f : futures) {
            auto status = f.wait_for(std::chrono::seconds(OSMORDRED_TIMEOUT_SECONDS));
            if (status == std::future_status::ready) {
                results.push_back(f.get());
            } else {
                results.push_back(std::vector<double>(3585, std::numeric_limits<double>::quiet_NaN()));
            }
        }
        
        return results;
    }

    // v2.0: Batch version from mol objects: PRESERVES tautomer canonical.
    // Python binding uses mol.ToBinary() -> MolPickler::molFromPickle -> these mols.
    RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>> calcOsmordredBatchFromMols(
        const std::vector<const ROMol*>& mols, int n_jobs) {

        std::vector<std::vector<double>> results;
        results.reserve(mols.size());

        unsigned int nThreads = getNumThreadsToUse(n_jobs);
        const size_t nFeatures = 3585;
        const std::vector<double> nanRow(nFeatures, std::numeric_limits<double>::quiet_NaN());

        if (nThreads <= 1 || mols.size() < 10) {
            for (const ROMol* mol : mols) {
                if (mol) {
                    try {
                        results.push_back(calcOsmordred(*mol));
                    } catch (...) {
                        results.push_back(nanRow);
                    }
                } else {
                    results.push_back(nanRow);
                }
            }
            return results;
        }

        std::vector<std::future<std::vector<double>>> futures;
        futures.reserve(mols.size());
        for (size_t idx = 0; idx < mols.size(); ++idx) {
            const ROMol* mol = mols[idx];
            futures.emplace_back(std::async(std::launch::async, [mol, nFeatures]() {
                if (mol) {
                    try {
                        return calcOsmordred(*mol);
                    } catch (...) {
                        return std::vector<double>(nFeatures, std::numeric_limits<double>::quiet_NaN());
                    }
                }
                return std::vector<double>(nFeatures, std::numeric_limits<double>::quiet_NaN());
            }));
        }
        for (auto& f : futures) {
            auto status = f.wait_for(std::chrono::seconds(OSMORDRED_TIMEOUT_SECONDS));
            if (status == std::future_status::ready) {
                results.push_back(f.get());
            } else {
                results.push_back(nanRow);
            }
        }
        return results;
    }

    // v2.0: Get descriptor names in the same order as calcOsmordred returns values
    std::vector<std::string> getOsmordredDescriptorNames() {
        std::vector<std::string> names;
        names.reserve(3589);
        
// osmordredv3: real Mordred-style names (was addNames placeholders)
        // ABCIndex (2)
        names.insert(names.end(), {"ABC", "ABCGG"});
        // AcidBase (2)
        names.insert(names.end(), {"nAcid", "nBase"});
        // AdjacencyMatrix (12)
        names.insert(names.end(), {
            "SpAbs_A", "SpMax_A", "SpDiam_A", "SpAD_A", "SpMAD_A", "LogEE_A", "VE1_A", "VE2_A", "VE3_A", "VR1_A",
            "VR2_A", "VR3_A"
        });
        // Aromatic (2)
        names.insert(names.end(), {"nAromAtom", "nAromBond"});
        // AtomCount (17)
        names.insert(names.end(), {
            "nAtom", "nHeavyAtom", "nSpiro", "nBridgehead", "nHetero", "nH", "nB", "nC", "nN", "nO",
            "nS", "nP", "nF", "nCl", "nBr", "nI", "nX"
        });
        // Autocorrelation (606)
        names.insert(names.end(), {
            "ATS0dv", "ATS1dv", "ATS2dv", "ATS3dv", "ATS4dv", "ATS5dv", "ATS6dv", "ATS7dv", "ATS8dv", "ATS0d",
            "ATS1d", "ATS2d", "ATS3d", "ATS4d", "ATS5d", "ATS6d", "ATS7d", "ATS8d", "ATS0s", "ATS1s",
            "ATS2s", "ATS3s", "ATS4s", "ATS5s", "ATS6s", "ATS7s", "ATS8s", "ATS0Z", "ATS1Z", "ATS2Z",
            "ATS3Z", "ATS4Z", "ATS5Z", "ATS6Z", "ATS7Z", "ATS8Z", "ATS0m", "ATS1m", "ATS2m", "ATS3m",
            "ATS4m", "ATS5m", "ATS6m", "ATS7m", "ATS8m", "ATS0v", "ATS1v", "ATS2v", "ATS3v", "ATS4v",
            "ATS5v", "ATS6v", "ATS7v", "ATS8v", "ATS0se", "ATS1se", "ATS2se", "ATS3se", "ATS4se", "ATS5se",
            "ATS6se", "ATS7se", "ATS8se", "ATS0pe", "ATS1pe", "ATS2pe", "ATS3pe", "ATS4pe", "ATS5pe", "ATS6pe",
            "ATS7pe", "ATS8pe", "ATS0are", "ATS1are", "ATS2are", "ATS3are", "ATS4are", "ATS5are", "ATS6are", "ATS7are",
            "ATS8are", "ATS0p", "ATS1p", "ATS2p", "ATS3p", "ATS4p", "ATS5p", "ATS6p", "ATS7p", "ATS8p",
            "ATS0i", "ATS1i", "ATS2i", "ATS3i", "ATS4i", "ATS5i", "ATS6i", "ATS7i", "ATS8i", "AATS0dv",
            "AATS1dv", "AATS2dv", "AATS3dv", "AATS4dv", "AATS5dv", "AATS6dv", "AATS7dv", "AATS8dv", "AATS0d", "AATS1d",
            "AATS2d", "AATS3d", "AATS4d", "AATS5d", "AATS6d", "AATS7d", "AATS8d", "AATS0s", "AATS1s", "AATS2s",
            "AATS3s", "AATS4s", "AATS5s", "AATS6s", "AATS7s", "AATS8s", "AATS0Z", "AATS1Z", "AATS2Z", "AATS3Z",
            "AATS4Z", "AATS5Z", "AATS6Z", "AATS7Z", "AATS8Z", "AATS0m", "AATS1m", "AATS2m", "AATS3m", "AATS4m",
            "AATS5m", "AATS6m", "AATS7m", "AATS8m", "AATS0v", "AATS1v", "AATS2v", "AATS3v", "AATS4v", "AATS5v",
            "AATS6v", "AATS7v", "AATS8v", "AATS0se", "AATS1se", "AATS2se", "AATS3se", "AATS4se", "AATS5se", "AATS6se",
            "AATS7se", "AATS8se", "AATS0pe", "AATS1pe", "AATS2pe", "AATS3pe", "AATS4pe", "AATS5pe", "AATS6pe", "AATS7pe",
            "AATS8pe", "AATS0are", "AATS1are", "AATS2are", "AATS3are", "AATS4are", "AATS5are", "AATS6are", "AATS7are", "AATS8are",
            "AATS0p", "AATS1p", "AATS2p", "AATS3p", "AATS4p", "AATS5p", "AATS6p", "AATS7p", "AATS8p", "AATS0i",
            "AATS1i", "AATS2i", "AATS3i", "AATS4i", "AATS5i", "AATS6i", "AATS7i", "AATS8i", "ATSC0c", "ATSC1c",
            "ATSC2c", "ATSC3c", "ATSC4c", "ATSC5c", "ATSC6c", "ATSC7c", "ATSC8c", "ATSC0dv", "ATSC1dv", "ATSC2dv",
            "ATSC3dv", "ATSC4dv", "ATSC5dv", "ATSC6dv", "ATSC7dv", "ATSC8dv", "ATSC0d", "ATSC1d", "ATSC2d", "ATSC3d",
            "ATSC4d", "ATSC5d", "ATSC6d", "ATSC7d", "ATSC8d", "ATSC0s", "ATSC1s", "ATSC2s", "ATSC3s", "ATSC4s",
            "ATSC5s", "ATSC6s", "ATSC7s", "ATSC8s", "ATSC0Z", "ATSC1Z", "ATSC2Z", "ATSC3Z", "ATSC4Z", "ATSC5Z",
            "ATSC6Z", "ATSC7Z", "ATSC8Z", "ATSC0m", "ATSC1m", "ATSC2m", "ATSC3m", "ATSC4m", "ATSC5m", "ATSC6m",
            "ATSC7m", "ATSC8m", "ATSC0v", "ATSC1v", "ATSC2v", "ATSC3v", "ATSC4v", "ATSC5v", "ATSC6v", "ATSC7v",
            "ATSC8v", "ATSC0se", "ATSC1se", "ATSC2se", "ATSC3se", "ATSC4se", "ATSC5se", "ATSC6se", "ATSC7se", "ATSC8se",
            "ATSC0pe", "ATSC1pe", "ATSC2pe", "ATSC3pe", "ATSC4pe", "ATSC5pe", "ATSC6pe", "ATSC7pe", "ATSC8pe", "ATSC0are",
            "ATSC1are", "ATSC2are", "ATSC3are", "ATSC4are", "ATSC5are", "ATSC6are", "ATSC7are", "ATSC8are", "ATSC0p", "ATSC1p",
            "ATSC2p", "ATSC3p", "ATSC4p", "ATSC5p", "ATSC6p", "ATSC7p", "ATSC8p", "ATSC0i", "ATSC1i", "ATSC2i",
            "ATSC3i", "ATSC4i", "ATSC5i", "ATSC6i", "ATSC7i", "ATSC8i", "AATSC0c", "AATSC1c", "AATSC2c", "AATSC3c",
            "AATSC4c", "AATSC5c", "AATSC6c", "AATSC7c", "AATSC8c", "AATSC0dv", "AATSC1dv", "AATSC2dv", "AATSC3dv", "AATSC4dv",
            "AATSC5dv", "AATSC6dv", "AATSC7dv", "AATSC8dv", "AATSC0d", "AATSC1d", "AATSC2d", "AATSC3d", "AATSC4d", "AATSC5d",
            "AATSC6d", "AATSC7d", "AATSC8d", "AATSC0s", "AATSC1s", "AATSC2s", "AATSC3s", "AATSC4s", "AATSC5s", "AATSC6s",
            "AATSC7s", "AATSC8s", "AATSC0Z", "AATSC1Z", "AATSC2Z", "AATSC3Z", "AATSC4Z", "AATSC5Z", "AATSC6Z", "AATSC7Z",
            "AATSC8Z", "AATSC0m", "AATSC1m", "AATSC2m", "AATSC3m", "AATSC4m", "AATSC5m", "AATSC6m", "AATSC7m", "AATSC8m",
            "AATSC0v", "AATSC1v", "AATSC2v", "AATSC3v", "AATSC4v", "AATSC5v", "AATSC6v", "AATSC7v", "AATSC8v", "AATSC0se",
            "AATSC1se", "AATSC2se", "AATSC3se", "AATSC4se", "AATSC5se", "AATSC6se", "AATSC7se", "AATSC8se", "AATSC0pe", "AATSC1pe",
            "AATSC2pe", "AATSC3pe", "AATSC4pe", "AATSC5pe", "AATSC6pe", "AATSC7pe", "AATSC8pe", "AATSC0are", "AATSC1are", "AATSC2are",
            "AATSC3are", "AATSC4are", "AATSC5are", "AATSC6are", "AATSC7are", "AATSC8are", "AATSC0p", "AATSC1p", "AATSC2p", "AATSC3p",
            "AATSC4p", "AATSC5p", "AATSC6p", "AATSC7p", "AATSC8p", "AATSC0i", "AATSC1i", "AATSC2i", "AATSC3i", "AATSC4i",
            "AATSC5i", "AATSC6i", "AATSC7i", "AATSC8i", "MATS1c", "MATS2c", "MATS3c", "MATS4c", "MATS5c", "MATS6c",
            "MATS7c", "MATS8c", "MATS1dv", "MATS2dv", "MATS3dv", "MATS4dv", "MATS5dv", "MATS6dv", "MATS7dv", "MATS8dv",
            "MATS1d", "MATS2d", "MATS3d", "MATS4d", "MATS5d", "MATS6d", "MATS7d", "MATS8d", "MATS1s", "MATS2s",
            "MATS3s", "MATS4s", "MATS5s", "MATS6s", "MATS7s", "MATS8s", "MATS1Z", "MATS2Z", "MATS3Z", "MATS4Z",
            "MATS5Z", "MATS6Z", "MATS7Z", "MATS8Z", "MATS1m", "MATS2m", "MATS3m", "MATS4m", "MATS5m", "MATS6m",
            "MATS7m", "MATS8m", "MATS1v", "MATS2v", "MATS3v", "MATS4v", "MATS5v", "MATS6v", "MATS7v", "MATS8v",
            "MATS1se", "MATS2se", "MATS3se", "MATS4se", "MATS5se", "MATS6se", "MATS7se", "MATS8se", "MATS1pe", "MATS2pe",
            "MATS3pe", "MATS4pe", "MATS5pe", "MATS6pe", "MATS7pe", "MATS8pe", "MATS1are", "MATS2are", "MATS3are", "MATS4are",
            "MATS5are", "MATS6are", "MATS7are", "MATS8are", "MATS1p", "MATS2p", "MATS3p", "MATS4p", "MATS5p", "MATS6p",
            "MATS7p", "MATS8p", "MATS1i", "MATS2i", "MATS3i", "MATS4i", "MATS5i", "MATS6i", "MATS7i", "MATS8i",
            "GATS1c", "GATS2c", "GATS3c", "GATS4c", "GATS5c", "GATS6c", "GATS7c", "GATS8c", "GATS1dv", "GATS2dv",
            "GATS3dv", "GATS4dv", "GATS5dv", "GATS6dv", "GATS7dv", "GATS8dv", "GATS1d", "GATS2d", "GATS3d", "GATS4d",
            "GATS5d", "GATS6d", "GATS7d", "GATS8d", "GATS1s", "GATS2s", "GATS3s", "GATS4s", "GATS5s", "GATS6s",
            "GATS7s", "GATS8s", "GATS1Z", "GATS2Z", "GATS3Z", "GATS4Z", "GATS5Z", "GATS6Z", "GATS7Z", "GATS8Z",
            "GATS1m", "GATS2m", "GATS3m", "GATS4m", "GATS5m", "GATS6m", "GATS7m", "GATS8m", "GATS1v", "GATS2v",
            "GATS3v", "GATS4v", "GATS5v", "GATS6v", "GATS7v", "GATS8v", "GATS1se", "GATS2se", "GATS3se", "GATS4se",
            "GATS5se", "GATS6se", "GATS7se", "GATS8se", "GATS1pe", "GATS2pe", "GATS3pe", "GATS4pe", "GATS5pe", "GATS6pe",
            "GATS7pe", "GATS8pe", "GATS1are", "GATS2are", "GATS3are", "GATS4are", "GATS5are", "GATS6are", "GATS7are", "GATS8are",
            "GATS1p", "GATS2p", "GATS3p", "GATS4p", "GATS5p", "GATS6p", "GATS7p", "GATS8p", "GATS1i", "GATS2i",
            "GATS3i", "GATS4i", "GATS5i", "GATS6i", "GATS7i", "GATS8i"
        });
        // BCUT (24)
        names.insert(names.end(), {
            "BCUTc-1l", "BCUTc-1h", "BCUTdv-1l", "BCUTdv-1h", "BCUTd-1l", "BCUTd-1h", "BCUTs-1l", "BCUTs-1h", "BCUTZ-1l", "BCUTZ-1h",
            "BCUTm-1l", "BCUTm-1h", "BCUTv-1l", "BCUTv-1h", "BCUTse-1l", "BCUTse-1h", "BCUTpe-1l", "BCUTpe-1h", "BCUTare-1l", "BCUTare-1h",
            "BCUTp-1l", "BCUTp-1h", "BCUTi-1l", "BCUTi-1h"
        });
        // BalabanJ (1)
        names.insert(names.end(), {"J"});
        // BaryszMatrix (104)
        names.insert(names.end(), {
            "SpAbs_DzZ", "SpMax_DzZ", "SpDiam_DzZ", "SpAD_DzZ", "SpMAD_DzZ", "LogEE_DzZ", "SM1_DzZ", "VE1_DzZ", "VE2_DzZ", "VE3_DzZ",
            "VR1_DzZ", "VR2_DzZ", "VR3_DzZ", "SpAbs_Dzm", "SpMax_Dzm", "SpDiam_Dzm", "SpAD_Dzm", "SpMAD_Dzm", "LogEE_Dzm", "SM1_Dzm",
            "VE1_Dzm", "VE2_Dzm", "VE3_Dzm", "VR1_Dzm", "VR2_Dzm", "VR3_Dzm", "SpAbs_Dzv", "SpMax_Dzv", "SpDiam_Dzv", "SpAD_Dzv",
            "SpMAD_Dzv", "LogEE_Dzv", "SM1_Dzv", "VE1_Dzv", "VE2_Dzv", "VE3_Dzv", "VR1_Dzv", "VR2_Dzv", "VR3_Dzv", "SpAbs_Dzse",
            "SpMax_Dzse", "SpDiam_Dzse", "SpAD_Dzse", "SpMAD_Dzse", "LogEE_Dzse", "SM1_Dzse", "VE1_Dzse", "VE2_Dzse", "VE3_Dzse", "VR1_Dzse",
            "VR2_Dzse", "VR3_Dzse", "SpAbs_Dzpe", "SpMax_Dzpe", "SpDiam_Dzpe", "SpAD_Dzpe", "SpMAD_Dzpe", "LogEE_Dzpe", "SM1_Dzpe", "VE1_Dzpe",
            "VE2_Dzpe", "VE3_Dzpe", "VR1_Dzpe", "VR2_Dzpe", "VR3_Dzpe", "SpAbs_Dzare", "SpMax_Dzare", "SpDiam_Dzare", "SpAD_Dzare", "SpMAD_Dzare",
            "LogEE_Dzare", "SM1_Dzare", "VE1_Dzare", "VE2_Dzare", "VE3_Dzare", "VR1_Dzare", "VR2_Dzare", "VR3_Dzare", "SpAbs_Dzp", "SpMax_Dzp",
            "SpDiam_Dzp", "SpAD_Dzp", "SpMAD_Dzp", "LogEE_Dzp", "SM1_Dzp", "VE1_Dzp", "VE2_Dzp", "VE3_Dzp", "VR1_Dzp", "VR2_Dzp",
            "VR3_Dzp", "SpAbs_Dzi", "SpMax_Dzi", "SpDiam_Dzi", "SpAD_Dzi", "SpMAD_Dzi", "LogEE_Dzi", "SM1_Dzi", "VE1_Dzi", "VE2_Dzi",
            "VE3_Dzi", "VR1_Dzi", "VR2_Dzi", "VR3_Dzi"
        });
        // BertzCT (1)
        names.insert(names.end(), {"BertzCT"});
        // BondCount (9)
        names.insert(names.end(), {"nBonds", "nBondsO", "nBondsS", "nBondsD", "nBondsT", "nBondsA", "nBondsM", "nBondsKS", "nBondsKD"});
        // RNCGRPCG (2)
        names.insert(names.end(), {"RNCG", "RPCG"});
        // CarbonTypes (11)
        names.insert(names.end(), {
            "C1SP1", "C2SP1", "C1SP2", "C2SP2", "C3SP2", "C1SP3", "C2SP3", "C3SP3", "C4SP3", "HybRatio",
            "FCSP3"
        });
        // Chi (56)
        names.insert(names.end(), {
            "Xch-3d", "Xch-4d", "Xch-5d", "Xch-6d", "Xch-7d", "Xch-3dv", "Xch-4dv", "Xch-5dv", "Xch-6dv", "Xch-7dv",
            "Xc-3d", "Xc-4d", "Xc-5d", "Xc-6d", "Xc-3dv", "Xc-4dv", "Xc-5dv", "Xc-6dv", "Xpc-4d", "Xpc-5d",
            "Xpc-6d", "Xpc-4dv", "Xpc-5dv", "Xpc-6dv", "Xp-0d", "Xp-1d", "Xp-2d", "Xp-3d", "Xp-4d", "Xp-5d",
            "Xp-6d", "Xp-7d", "AXp-0d", "AXp-1d", "AXp-2d", "AXp-3d", "AXp-4d", "AXp-5d", "AXp-6d", "AXp-7d",
            "Xp-0dv", "Xp-1dv", "Xp-2dv", "Xp-3dv", "Xp-4dv", "Xp-5dv", "Xp-6dv", "Xp-7dv", "AXp-0dv", "AXp-1dv",
            "AXp-2dv", "AXp-3dv", "AXp-4dv", "AXp-5dv", "AXp-6dv", "AXp-7dv"
        });
        // Constitutional (16)
        names.insert(names.end(), {
            "SZ", "Sm", "Sv", "Sse", "Spe", "Sare", "Sp", "Si", "MZ", "Mm",
            "Mv", "Mse", "Mpe", "Mare", "Mp", "Mi"
        });
        // DetourMatrix (14)
        names.insert(names.end(), {
            "SpAbs_Dt", "SpMax_Dt", "SpDiam_Dt", "SpAD_Dt", "SpMAD_Dt", "LogEE_Dt", "SM1_Dt", "VE1_Dt", "VE2_Dt", "VE3_Dt",
            "VR1_Dt", "VR2_Dt", "VR3_Dt", "DetourIndex"
        });
        // DistanceMatrix (12)
        names.insert(names.end(), {
            "SpAbs_D", "SpMax_D", "SpDiam_D", "SpAD_D", "SpMAD_D", "LogEE_D", "VE1_D", "VE2_D", "VE3_D", "VR1_D",
            "VR2_D", "VR3_D"
        });
        // EState (404)
        names.insert(names.end(), {
            "EState_1", "EState_2", "EState_3", "EState_4", "EState_5", "EState_6", "NsCH3", "EState_8", "EState_9", "EState_10",
            "EState_11", "NaaCH", "NsssCH", "EState_14", "EState_15", "NdssC", "NaasC", "EState_18", "EState_19", "EState_20",
            "NsNH2", "EState_22", "EState_23", "EState_24", "NaaNH", "EState_26", "EState_27", "EState_28", "NaaN", "EState_30",
            "EState_31", "EState_32", "EState_33", "NsOH", "NdO", "EState_36", "EState_37", "EState_38", "EState_39", "EState_40",
            "EState_41", "EState_42", "EState_43", "EState_44", "EState_45", "EState_46", "EState_47", "EState_48", "EState_49", "EState_50",
            "EState_51", "EState_52", "EState_53", "EState_54", "EState_55", "EState_56", "EState_57", "EState_58", "EState_59", "EState_60",
            "EState_61", "EState_62", "EState_63", "EState_64", "EState_65", "EState_66", "EState_67", "EState_68", "EState_69", "EState_70",
            "EState_71", "EState_72", "EState_73", "EState_74", "EState_75", "EState_76", "EState_77", "EState_78", "EState_79", "EState_80",
            "EState_81", "EState_82", "EState_83", "EState_84", "EState_85", "EState_86", "EState_87", "EState_88", "EState_89", "EState_90",
            "EState_91", "EState_92", "EState_93", "EState_94", "EState_95", "EState_96", "EState_97", "EState_98", "EState_99", "EState_100",
            "EState_101", "EState_102", "EState_103", "EState_104", "EState_105", "EState_106", "EState_107", "SsCH3", "EState_109", "EState_110",
            "EState_111", "EState_112", "EState_113", "SsssCH", "EState_115", "EState_116", "SdssC", "SaasC", "EState_119", "EState_120",
            "EState_121", "SsNH2", "EState_123", "EState_124", "EState_125", "SaaNH", "EState_127", "EState_128", "EState_129", "SaaN",
            "EState_131", "EState_132", "EState_133", "EState_134", "SsOH", "SdO", "EState_137", "EState_138", "EState_139", "EState_140",
            "EState_141", "EState_142", "EState_143", "EState_144", "EState_145", "EState_146", "EState_147", "EState_148", "EState_149", "EState_150",
            "EState_151", "EState_152", "EState_153", "EState_154", "EState_155", "EState_156", "EState_157", "EState_158", "EState_159", "EState_160",
            "EState_161", "EState_162", "EState_163", "EState_164", "EState_165", "EState_166", "EState_167", "EState_168", "EState_169", "EState_170",
            "EState_171", "EState_172", "EState_173", "EState_174", "EState_175", "EState_176", "EState_177", "EState_178", "EState_179", "EState_180",
            "EState_181", "EState_182", "EState_183", "EState_184", "EState_185", "EState_186", "EState_187", "EState_188", "EState_189", "EState_190",
            "EState_191", "EState_192", "EState_193", "EState_194", "EState_195", "EState_196", "EState_197", "EState_198", "EState_199", "EState_200",
            "EState_201", "EState_202", "EState_203", "EState_204", "EState_205", "EState_206", "EState_207", "EState_208", "EState_209", "EState_210",
            "EState_211", "EState_212", "EState_213", "EState_214", "EState_215", "EState_216", "EState_217", "EState_218", "EState_219", "EState_220",
            "EState_221", "EState_222", "EState_223", "EState_224", "EState_225", "EState_226", "EState_227", "EState_228", "EState_229", "EState_230",
            "EState_231", "EState_232", "EState_233", "EState_234", "EState_235", "EState_236", "EState_237", "EState_238", "EState_239", "EState_240",
            "EState_241", "EState_242", "EState_243", "EState_244", "EState_245", "EState_246", "EState_247", "EState_248", "EState_249", "EState_250",
            "EState_251", "EState_252", "EState_253", "EState_254", "EState_255", "EState_256", "EState_257", "EState_258", "EState_259", "EState_260",
            "EState_261", "EState_262", "EState_263", "EState_264", "EState_265", "EState_266", "EState_267", "EState_268", "EState_269", "EState_270",
            "EState_271", "EState_272", "EState_273", "EState_274", "EState_275", "EState_276", "EState_277", "EState_278", "EState_279", "EState_280",
            "EState_281", "EState_282", "EState_283", "EState_284", "EState_285", "EState_286", "EState_287", "EState_288", "EState_289", "EState_290",
            "EState_291", "EState_292", "EState_293", "EState_294", "EState_295", "EState_296", "EState_297", "EState_298", "EState_299", "EState_300",
            "EState_301", "EState_302", "EState_303", "EState_304", "EState_305", "EState_306", "EState_307", "EState_308", "EState_309", "EState_310",
            "EState_311", "EState_312", "EState_313", "EState_314", "EState_315", "EState_316", "EState_317", "EState_318", "EState_319", "EState_320",
            "EState_321", "EState_322", "EState_323", "EState_324", "EState_325", "EState_326", "EState_327", "EState_328", "EState_329", "EState_330",
            "EState_331", "EState_332", "EState_333", "EState_334", "EState_335", "EState_336", "EState_337", "EState_338", "EState_339", "EState_340",
            "EState_341", "EState_342", "EState_343", "EState_344", "EState_345", "EState_346", "EState_347", "EState_348", "EState_349", "EState_350",
            "EState_351", "EState_352", "EState_353", "EState_354", "EState_355", "EState_356", "EState_357", "EState_358", "EState_359", "EState_360",
            "EState_361", "EState_362", "EState_363", "EState_364", "EState_365", "EState_366", "EState_367", "EState_368", "EState_369", "EState_370",
            "EState_371", "EState_372", "EState_373", "EState_374", "EState_375", "EState_376", "EState_377", "EState_378", "EState_379", "EState_380",
            "EState_381", "EState_382", "EState_383", "EState_384", "EState_385", "EState_386", "EState_387", "EState_388", "EState_389", "EState_390",
            "EState_391", "EState_392", "EState_393", "EState_394", "EState_395", "EState_396", "EState_397", "EState_398", "EState_399", "EState_400",
            "EState_401", "EState_402", "EState_403", "EState_404"
        });
        // EccentricConnectivityIndex (1)
        names.insert(names.end(), {"ECIndex"});
        // ExtendedTopochemicalAtom (45)
        names.insert(names.end(), {
            "ETA_alpha", "AETA_alpha", "ETA_shape_p", "ETA_shape_y", "ETA_shape_x", "ETA_beta", "AETA_beta", "ETA_beta_s", "AETA_beta_s", "ETA_beta_ns",
            "AETA_beta_ns", "ETA_beta_ns_d", "AETA_beta_ns_d", "ETA_eta", "AETA_eta", "ETA_eta_L", "AETA_eta_L", "ETA_eta_R", "AETA_eta_R", "ETA_eta_RL",
            "AETA_eta_RL", "ETA_eta_F", "AETA_eta_F", "ETA_eta_FL", "AETA_eta_FL", "ETA_eta_B", "AETA_eta_B", "ETA_eta_BR", "AETA_eta_BR", "ETA_dAlpha_A",
            "ETA_dAlpha_B", "ETA_epsilon_1", "ETA_epsilon_2", "ETA_epsilon_3", "ETA_epsilon_4", "ETA_epsilon_5", "ETA_dEpsilon_A", "ETA_dEpsilon_B", "ETA_dEpsilon_C", "ETA_dEpsilon_D",
            "ETA_dBeta", "AETA_dBeta", "ETA_psi_1", "ETA_dPsi_A", "ETA_dPsi_B"
        });
        // FragmentComplexity (1)
        names.insert(names.end(), {"fragCpx"});
        // Framework (1)
        names.insert(names.end(), {"fMF"});
        // HydrogenBond (2)
        names.insert(names.end(), {"nHBAcc", "nHBDon"});
        // LogS (1)
        names.insert(names.end(), {"FilterItLogS"});
        // InformationContent (42)
        names.insert(names.end(), {
            "IC0", "IC1", "IC2", "IC3", "IC4", "IC5", "TIC0", "TIC1", "TIC2", "TIC3",
            "TIC4", "TIC5", "SIC0", "SIC1", "SIC2", "SIC3", "SIC4", "SIC5", "BIC0", "BIC1",
            "BIC2", "BIC3", "BIC4", "BIC5", "CIC0", "CIC1", "CIC2", "CIC3", "CIC4", "CIC5",
            "MIC0", "MIC1", "MIC2", "MIC3", "MIC4", "MIC5", "ZMIC0", "ZMIC1", "ZMIC2", "ZMIC3",
            "ZMIC4", "ZMIC5"
        });
        // KappaShapeIndex (3)
        names.insert(names.end(), {"Kier1", "Kier2", "Kier3"});
        // Lipinski (2)
        names.insert(names.end(), {"Lipinski", "GhoseFilter"});
        // McGowanVolume (1)
        names.insert(names.end(), {"VMcGowan"});
        // MoeType (58): LabuteASA + VSA descriptors (RDKit binning); was 54 MoeType_ placeholders
        names.insert(names.end(), {
            "LabuteASA", "PEOE_VSA1", "PEOE_VSA2", "PEOE_VSA3", "PEOE_VSA4", "PEOE_VSA5", "PEOE_VSA6", "PEOE_VSA7", "PEOE_VSA8", "PEOE_VSA9",
            "PEOE_VSA10", "PEOE_VSA11", "PEOE_VSA12", "PEOE_VSA13", "PEOE_VSA14", "SMR_VSA1", "SMR_VSA2", "SMR_VSA3", "SMR_VSA4", "SMR_VSA5",
            "SMR_VSA6", "SMR_VSA7", "SMR_VSA8", "SMR_VSA9", "SMR_VSA10", "SlogP_VSA1", "SlogP_VSA2", "SlogP_VSA3", "SlogP_VSA4", "SlogP_VSA5",
            "SlogP_VSA6", "SlogP_VSA7", "SlogP_VSA8", "SlogP_VSA9", "SlogP_VSA10", "SlogP_VSA11", "SlogP_VSA12", "EState_VSA1", "EState_VSA2", "EState_VSA3",
            "EState_VSA4", "EState_VSA5", "EState_VSA6", "EState_VSA7", "EState_VSA8", "EState_VSA9", "EState_VSA10", "EState_VSA11", "VSA_EState1", "VSA_EState2",
            "VSA_EState3", "VSA_EState4", "VSA_EState5", "VSA_EState6", "VSA_EState7", "VSA_EState8", "VSA_EState9", "VSA_EState10"
        });
        // MolecularDistanceEdge (19)
        names.insert(names.end(), {
            "MDEC-11", "MDEC-12", "MDEC-13", "MDEC-14", "MDEC-22", "MDEC-23", "MDEC-24", "MDEC-33", "MDEC-34", "MDEC-44",
            "MDEO-11", "MDEO-12", "MDEO-22", "MDEN-11", "MDEN-12", "MDEN-13", "MDEN-22", "MDEN-23", "MDEN-33"
        });
        // MolecularId (12)
        names.insert(names.end(), {
            "MID", "AMID", "MID_h", "AMID_h", "MID_C", "AMID_C", "MID_N", "AMID_N", "MID_O", "AMID_O",
            "MID_X", "AMID_X"
        });
        // PathCount (21)
        names.insert(names.end(), {
            "MPC2", "MPC3", "MPC4", "MPC5", "MPC6", "MPC7", "MPC8", "MPC9", "MPC10", "TMPC10",
            "piPC1", "piPC2", "piPC3", "piPC4", "piPC5", "piPC6", "piPC7", "piPC8", "piPC9", "piPC10",
            "TpiPC10"
        });
        // Polarizability (2)
        names.insert(names.end(), {"apol78", "bpol78"});
        // RingCount (138)
        names.insert(names.end(), {
            "nRing", "n3Ring", "n4Ring", "n5Ring", "n6Ring", "n7Ring", "n8Ring", "n9Ring", "n10Ring", "n11Ring",
            "n12Ring", "nG12Ring", "nHRing", "n3HRing", "n4HRing", "n5HRing", "n6HRing", "n7HRing", "n8HRing", "n9HRing",
            "n10HRing", "n11HRing", "n12HRing", "nG12HRing", "naRing", "n3aRing", "n4aRing", "n5aRing", "n6aRing", "n7aRing",
            "n8aRing", "n9aRing", "n10aRing", "n11aRing", "n12aRing", "nG12aRing", "naHRing", "n3aHRing", "n4aHRing", "n5aHRing",
            "n6aHRing", "n7aHRing", "n8aHRing", "n9aHRing", "n10aHRing", "n11aHRing", "n12aHRing", "nG12aHRing", "nARing", "n3ARing",
            "n4ARing", "n5ARing", "n6ARing", "n7ARing", "n8ARing", "n9ARing", "n10ARing", "n11ARing", "n12ARing", "nG12ARing",
            "nAHRing", "n3AHRing", "n4AHRing", "n5AHRing", "n6AHRing", "n7AHRing", "n8AHRing", "n9AHRing", "n10AHRing", "n11AHRing",
            "n12AHRing", "nG12AHRing", "nFRing", "n4FRing", "n5FRing", "n6FRing", "n7FRing", "n8FRing", "n9FRing", "n10FRing",
            "n11FRing", "n12FRing", "nG12FRing", "nFHRing", "n4FHRing", "n5FHRing", "n6FHRing", "n7FHRing", "n8FHRing", "n9FHRing",
            "n10FHRing", "n11FHRing", "n12FHRing", "nG12FHRing", "nFaRing", "n4FaRing", "n5FaRing", "n6FaRing", "n7FaRing", "n8FaRing",
            "n9FaRing", "n10FaRing", "n11FaRing", "n12FaRing", "nG12FaRing", "nFaHRing", "n4FaHRing", "n5FaHRing", "n6FaHRing", "n7FaHRing",
            "n8FaHRing", "n9FaHRing", "n10FaHRing", "n11FaHRing", "n12FaHRing", "nG12FaHRing", "nFARing", "n4FARing", "n5FARing", "n6FARing",
            "n7FARing", "n8FARing", "n9FARing", "n10FARing", "n11FARing", "n12FARing", "nG12FARing", "nFAHRing", "n4FAHRing", "n5FAHRing",
            "n6FAHRing", "n7FAHRing", "n8FAHRing", "n9FAHRing", "n10FAHRing", "n11FAHRing", "n12FAHRing", "nG12FAHRing"
        });
        // RotatableBond (2)
        names.insert(names.end(), {"nRot", "RotRatio"});
        // SLogP (2)
        names.insert(names.end(), {"SLogP", "SMR"});
        // TopoPSA (2)
        names.insert(names.end(), {"TopoPSA(NO)", "TopoPSA"});
        // TopologicalCharge (21)
        names.insert(names.end(), {
            "GGI1", "GGI2", "GGI3", "GGI4", "GGI5", "GGI6", "GGI7", "GGI8", "GGI9", "GGI10",
            "JGI1", "JGI2", "JGI3", "JGI4", "JGI5", "JGI6", "JGI7", "JGI8", "JGI9", "JGI10",
            "JGT10"
        });
        // TopologicalIndex (4)
        names.insert(names.end(), {"Diameter", "Radius", "TopoShapeIndex", "PetitjeanIndex"});
        // VdwVolumeABC (1)
        names.insert(names.end(), {"Vabc"});
        // VertexAdjacencyInformation (1)
        names.insert(names.end(), {"VAdjMat"});
        // WalkCount (21)
        names.insert(names.end(), {
            "MWC01", "MWC02", "MWC03", "MWC04", "MWC05", "MWC06", "MWC07", "MWC08", "MWC09", "MWC10",
            "TMWC10", "SRW02", "SRW03", "SRW04", "SRW05", "SRW06", "SRW07", "SRW08", "SRW09", "SRW10",
            "TSRW10"
        });
        // Weight (2)
        names.insert(names.end(), {"MW", "AMW"});
        // WienerIndex (2)
        names.insert(names.end(), {"WPath", "WPol"});
        // ZagrebIndex (4)
        names.insert(names.end(), {"Zagreb1", "Zagreb2", "mZagreb1", "mZagreb2"});
        // Pol (1)
        names.insert(names.end(), {"Pol"});
        // MR (1)
        names.insert(names.end(), {"MR"});
        // Flexibility (1)
        names.insert(names.end(), {"Flexibility"});
        // Schultz (1)
        names.insert(names.end(), {"Schultz"});
        // AlphaKappaShapeIndex (3)
        names.insert(names.end(), {"KierAlpha1", "KierAlpha2", "KierAlpha3"});
        // HEState (88)
        names.insert(names.end(), {
            "HEState_1", "HEState_2", "HEState_3", "HEState_4", "HEState_5", "HEState_6", "HEState_7", "HEState_8", "HEState_9", "HEState_10",
            "HEState_11", "HEState_12", "HEState_13", "HEState_14", "HEState_15", "HEState_16", "HEState_17", "HEState_18", "HEState_19", "HEState_20",
            "HEState_21", "HEState_22", "HEState_23", "HEState_24", "HEState_25", "HEState_26", "HEState_27", "HEState_28", "HEState_29", "HEState_30",
            "HEState_31", "HEState_32", "HEState_33", "HEState_34", "HEState_35", "HEState_36", "HEState_37", "HEState_38", "HEState_39", "HEState_40",
            "HEState_41", "HEState_42", "HEState_43", "HEState_44", "HEState_45", "HEState_46", "HEState_47", "HEState_48", "HEState_49", "HEState_50",
            "HEState_51", "HEState_52", "HEState_53", "HEState_54", "HEState_55", "HEState_56", "HEState_57", "HEState_58", "HEState_59", "HEState_60",
            "HEState_61", "HEState_62", "HEState_63", "HEState_64", "HEState_65", "HEState_66", "HEState_67", "HEState_68", "HEState_69", "HEState_70",
            "HEState_71", "HEState_72", "HEState_73", "HEState_74", "HEState_75", "HEState_76", "HEState_77", "HEState_78", "HEState_79", "HEState_80",
            "HEState_81", "HEState_82", "HEState_83", "HEState_84", "HEState_85", "HEState_86", "HEState_87", "HEState_88"
        });
        // BEState (1460)
        names.insert(names.end(), {
            "BEState_1", "BEState_2", "BEState_3", "BEState_4", "BEState_5", "BEState_6", "BEState_7", "BEState_8", "BEState_9", "BEState_10",
            "BEState_11", "BEState_12", "BEState_13", "BEState_14", "BEState_15", "BEState_16", "BEState_17", "BEState_18", "BEState_19", "BEState_20",
            "BEState_21", "BEState_22", "BEState_23", "BEState_24", "BEState_25", "BEState_26", "BEState_27", "BEState_28", "BEState_29", "BEState_30",
            "BEState_31", "BEState_32", "BEState_33", "BEState_34", "BEState_35", "BEState_36", "BEState_37", "BEState_38", "BEState_39", "BEState_40",
            "BEState_41", "BEState_42", "BEState_43", "BEState_44", "BEState_45", "BEState_46", "BEState_47", "BEState_48", "BEState_49", "BEState_50",
            "BEState_51", "BEState_52", "BEState_53", "BEState_54", "BEState_55", "BEState_56", "BEState_57", "BEState_58", "BEState_59", "BEState_60",
            "BEState_61", "BEState_62", "BEState_63", "BEState_64", "BEState_65", "BEState_66", "BEState_67", "BEState_68", "BEState_69", "BEState_70",
            "BEState_71", "BEState_72", "BEState_73", "BEState_74", "BEState_75", "BEState_76", "BEState_77", "BEState_78", "BEState_79", "BEState_80",
            "BEState_81", "BEState_82", "BEState_83", "BEState_84", "BEState_85", "BEState_86", "BEState_87", "BEState_88", "BEState_89", "BEState_90",
            "BEState_91", "BEState_92", "BEState_93", "BEState_94", "BEState_95", "BEState_96", "BEState_97", "BEState_98", "BEState_99", "BEState_100",
            "BEState_101", "BEState_102", "BEState_103", "BEState_104", "BEState_105", "BEState_106", "BEState_107", "BEState_108", "BEState_109", "BEState_110",
            "BEState_111", "BEState_112", "BEState_113", "BEState_114", "BEState_115", "BEState_116", "BEState_117", "BEState_118", "BEState_119", "BEState_120",
            "BEState_121", "BEState_122", "BEState_123", "BEState_124", "BEState_125", "BEState_126", "BEState_127", "BEState_128", "BEState_129", "BEState_130",
            "BEState_131", "BEState_132", "BEState_133", "BEState_134", "BEState_135", "BEState_136", "BEState_137", "BEState_138", "BEState_139", "BEState_140",
            "BEState_141", "BEState_142", "BEState_143", "BEState_144", "BEState_145", "BEState_146", "BEState_147", "BEState_148", "BEState_149", "BEState_150",
            "BEState_151", "BEState_152", "BEState_153", "BEState_154", "BEState_155", "BEState_156", "BEState_157", "BEState_158", "BEState_159", "BEState_160",
            "BEState_161", "BEState_162", "BEState_163", "BEState_164", "BEState_165", "BEState_166", "BEState_167", "BEState_168", "BEState_169", "BEState_170",
            "BEState_171", "BEState_172", "BEState_173", "BEState_174", "BEState_175", "BEState_176", "BEState_177", "BEState_178", "BEState_179", "BEState_180",
            "BEState_181", "BEState_182", "BEState_183", "BEState_184", "BEState_185", "BEState_186", "BEState_187", "BEState_188", "BEState_189", "BEState_190",
            "BEState_191", "BEState_192", "BEState_193", "BEState_194", "BEState_195", "BEState_196", "BEState_197", "BEState_198", "BEState_199", "BEState_200",
            "BEState_201", "BEState_202", "BEState_203", "BEState_204", "BEState_205", "BEState_206", "BEState_207", "BEState_208", "BEState_209", "BEState_210",
            "BEState_211", "BEState_212", "BEState_213", "BEState_214", "BEState_215", "BEState_216", "BEState_217", "BEState_218", "BEState_219", "BEState_220",
            "BEState_221", "BEState_222", "BEState_223", "BEState_224", "BEState_225", "BEState_226", "BEState_227", "BEState_228", "BEState_229", "BEState_230",
            "BEState_231", "BEState_232", "BEState_233", "BEState_234", "BEState_235", "BEState_236", "BEState_237", "BEState_238", "BEState_239", "BEState_240",
            "BEState_241", "BEState_242", "BEState_243", "BEState_244", "BEState_245", "BEState_246", "BEState_247", "BEState_248", "BEState_249", "BEState_250",
            "BEState_251", "BEState_252", "BEState_253", "BEState_254", "BEState_255", "BEState_256", "BEState_257", "BEState_258", "BEState_259", "BEState_260",
            "BEState_261", "BEState_262", "BEState_263", "BEState_264", "BEState_265", "BEState_266", "BEState_267", "BEState_268", "BEState_269", "BEState_270",
            "BEState_271", "BEState_272", "BEState_273", "BEState_274", "BEState_275", "BEState_276", "BEState_277", "BEState_278", "BEState_279", "BEState_280",
            "BEState_281", "BEState_282", "BEState_283", "BEState_284", "BEState_285", "BEState_286", "BEState_287", "BEState_288", "BEState_289", "BEState_290",
            "BEState_291", "BEState_292", "BEState_293", "BEState_294", "BEState_295", "BEState_296", "BEState_297", "BEState_298", "BEState_299", "BEState_300",
            "BEState_301", "BEState_302", "BEState_303", "BEState_304", "BEState_305", "BEState_306", "BEState_307", "BEState_308", "BEState_309", "BEState_310",
            "BEState_311", "BEState_312", "BEState_313", "BEState_314", "BEState_315", "BEState_316", "BEState_317", "BEState_318", "BEState_319", "BEState_320",
            "BEState_321", "BEState_322", "BEState_323", "BEState_324", "BEState_325", "BEState_326", "BEState_327", "BEState_328", "BEState_329", "BEState_330",
            "BEState_331", "BEState_332", "BEState_333", "BEState_334", "BEState_335", "BEState_336", "BEState_337", "BEState_338", "BEState_339", "BEState_340",
            "BEState_341", "BEState_342", "BEState_343", "BEState_344", "BEState_345", "BEState_346", "BEState_347", "BEState_348", "BEState_349", "BEState_350",
            "BEState_351", "BEState_352", "BEState_353", "BEState_354", "BEState_355", "BEState_356", "BEState_357", "BEState_358", "BEState_359", "BEState_360",
            "BEState_361", "BEState_362", "BEState_363", "BEState_364", "BEState_365", "BEState_366", "BEState_367", "BEState_368", "BEState_369", "BEState_370",
            "BEState_371", "BEState_372", "BEState_373", "BEState_374", "BEState_375", "BEState_376", "BEState_377", "BEState_378", "BEState_379", "BEState_380",
            "BEState_381", "BEState_382", "BEState_383", "BEState_384", "BEState_385", "BEState_386", "BEState_387", "BEState_388", "BEState_389", "BEState_390",
            "BEState_391", "BEState_392", "BEState_393", "BEState_394", "BEState_395", "BEState_396", "BEState_397", "BEState_398", "BEState_399", "BEState_400",
            "BEState_401", "BEState_402", "BEState_403", "BEState_404", "BEState_405", "BEState_406", "BEState_407", "BEState_408", "BEState_409", "BEState_410",
            "BEState_411", "BEState_412", "BEState_413", "BEState_414", "BEState_415", "BEState_416", "BEState_417", "BEState_418", "BEState_419", "BEState_420",
            "BEState_421", "BEState_422", "BEState_423", "BEState_424", "BEState_425", "BEState_426", "BEState_427", "BEState_428", "BEState_429", "BEState_430",
            "BEState_431", "BEState_432", "BEState_433", "BEState_434", "BEState_435", "BEState_436", "BEState_437", "BEState_438", "BEState_439", "BEState_440",
            "BEState_441", "BEState_442", "BEState_443", "BEState_444", "BEState_445", "BEState_446", "BEState_447", "BEState_448", "BEState_449", "BEState_450",
            "BEState_451", "BEState_452", "BEState_453", "BEState_454", "BEState_455", "BEState_456", "BEState_457", "BEState_458", "BEState_459", "BEState_460",
            "BEState_461", "BEState_462", "BEState_463", "BEState_464", "BEState_465", "BEState_466", "BEState_467", "BEState_468", "BEState_469", "BEState_470",
            "BEState_471", "BEState_472", "BEState_473", "BEState_474", "BEState_475", "BEState_476", "BEState_477", "BEState_478", "BEState_479", "BEState_480",
            "BEState_481", "BEState_482", "BEState_483", "BEState_484", "BEState_485", "BEState_486", "BEState_487", "BEState_488", "BEState_489", "BEState_490",
            "BEState_491", "BEState_492", "BEState_493", "BEState_494", "BEState_495", "BEState_496", "BEState_497", "BEState_498", "BEState_499", "BEState_500",
            "BEState_501", "BEState_502", "BEState_503", "BEState_504", "BEState_505", "BEState_506", "BEState_507", "BEState_508", "BEState_509", "BEState_510",
            "BEState_511", "BEState_512", "BEState_513", "BEState_514", "BEState_515", "BEState_516", "BEState_517", "BEState_518", "BEState_519", "BEState_520",
            "BEState_521", "BEState_522", "BEState_523", "BEState_524", "BEState_525", "BEState_526", "BEState_527", "BEState_528", "BEState_529", "BEState_530",
            "BEState_531", "BEState_532", "BEState_533", "BEState_534", "BEState_535", "BEState_536", "BEState_537", "BEState_538", "BEState_539", "BEState_540",
            "BEState_541", "BEState_542", "BEState_543", "BEState_544", "BEState_545", "BEState_546", "BEState_547", "BEState_548", "BEState_549", "BEState_550",
            "BEState_551", "BEState_552", "BEState_553", "BEState_554", "BEState_555", "BEState_556", "BEState_557", "BEState_558", "BEState_559", "BEState_560",
            "BEState_561", "BEState_562", "BEState_563", "BEState_564", "BEState_565", "BEState_566", "BEState_567", "BEState_568", "BEState_569", "BEState_570",
            "BEState_571", "BEState_572", "BEState_573", "BEState_574", "BEState_575", "BEState_576", "BEState_577", "BEState_578", "BEState_579", "BEState_580",
            "BEState_581", "BEState_582", "BEState_583", "BEState_584", "BEState_585", "BEState_586", "BEState_587", "BEState_588", "BEState_589", "BEState_590",
            "BEState_591", "BEState_592", "BEState_593", "BEState_594", "BEState_595", "BEState_596", "BEState_597", "BEState_598", "BEState_599", "BEState_600",
            "BEState_601", "BEState_602", "BEState_603", "BEState_604", "BEState_605", "BEState_606", "BEState_607", "BEState_608", "BEState_609", "BEState_610",
            "BEState_611", "BEState_612", "BEState_613", "BEState_614", "BEState_615", "BEState_616", "BEState_617", "BEState_618", "BEState_619", "BEState_620",
            "BEState_621", "BEState_622", "BEState_623", "BEState_624", "BEState_625", "BEState_626", "BEState_627", "BEState_628", "BEState_629", "BEState_630",
            "BEState_631", "BEState_632", "BEState_633", "BEState_634", "BEState_635", "BEState_636", "BEState_637", "BEState_638", "BEState_639", "BEState_640",
            "BEState_641", "BEState_642", "BEState_643", "BEState_644", "BEState_645", "BEState_646", "BEState_647", "BEState_648", "BEState_649", "BEState_650",
            "BEState_651", "BEState_652", "BEState_653", "BEState_654", "BEState_655", "BEState_656", "BEState_657", "BEState_658", "BEState_659", "BEState_660",
            "BEState_661", "BEState_662", "BEState_663", "BEState_664", "BEState_665", "BEState_666", "BEState_667", "BEState_668", "BEState_669", "BEState_670",
            "BEState_671", "BEState_672", "BEState_673", "BEState_674", "BEState_675", "BEState_676", "BEState_677", "BEState_678", "BEState_679", "BEState_680",
            "BEState_681", "BEState_682", "BEState_683", "BEState_684", "BEState_685", "BEState_686", "BEState_687", "BEState_688", "BEState_689", "BEState_690",
            "BEState_691", "BEState_692", "BEState_693", "BEState_694", "BEState_695", "BEState_696", "BEState_697", "BEState_698", "BEState_699", "BEState_700",
            "BEState_701", "BEState_702", "BEState_703", "BEState_704", "BEState_705", "BEState_706", "BEState_707", "BEState_708", "BEState_709", "BEState_710",
            "BEState_711", "BEState_712", "BEState_713", "BEState_714", "BEState_715", "BEState_716", "BEState_717", "BEState_718", "BEState_719", "BEState_720",
            "BEState_721", "BEState_722", "BEState_723", "BEState_724", "BEState_725", "BEState_726", "BEState_727", "BEState_728", "BEState_729", "BEState_730",
            "BEState_731", "BEState_732", "BEState_733", "BEState_734", "BEState_735", "BEState_736", "BEState_737", "BEState_738", "BEState_739", "BEState_740",
            "BEState_741", "BEState_742", "BEState_743", "BEState_744", "BEState_745", "BEState_746", "BEState_747", "BEState_748", "BEState_749", "BEState_750",
            "BEState_751", "BEState_752", "BEState_753", "BEState_754", "BEState_755", "BEState_756", "BEState_757", "BEState_758", "BEState_759", "BEState_760",
            "BEState_761", "BEState_762", "BEState_763", "BEState_764", "BEState_765", "BEState_766", "BEState_767", "BEState_768", "BEState_769", "BEState_770",
            "BEState_771", "BEState_772", "BEState_773", "BEState_774", "BEState_775", "BEState_776", "BEState_777", "BEState_778", "BEState_779", "BEState_780",
            "BEState_781", "BEState_782", "BEState_783", "BEState_784", "BEState_785", "BEState_786", "BEState_787", "BEState_788", "BEState_789", "BEState_790",
            "BEState_791", "BEState_792", "BEState_793", "BEState_794", "BEState_795", "BEState_796", "BEState_797", "BEState_798", "BEState_799", "BEState_800",
            "BEState_801", "BEState_802", "BEState_803", "BEState_804", "BEState_805", "BEState_806", "BEState_807", "BEState_808", "BEState_809", "BEState_810",
            "BEState_811", "BEState_812", "BEState_813", "BEState_814", "BEState_815", "BEState_816", "BEState_817", "BEState_818", "BEState_819", "BEState_820",
            "BEState_821", "BEState_822", "BEState_823", "BEState_824", "BEState_825", "BEState_826", "BEState_827", "BEState_828", "BEState_829", "BEState_830",
            "BEState_831", "BEState_832", "BEState_833", "BEState_834", "BEState_835", "BEState_836", "BEState_837", "BEState_838", "BEState_839", "BEState_840",
            "BEState_841", "BEState_842", "BEState_843", "BEState_844", "BEState_845", "BEState_846", "BEState_847", "BEState_848", "BEState_849", "BEState_850",
            "BEState_851", "BEState_852", "BEState_853", "BEState_854", "BEState_855", "BEState_856", "BEState_857", "BEState_858", "BEState_859", "BEState_860",
            "BEState_861", "BEState_862", "BEState_863", "BEState_864", "BEState_865", "BEState_866", "BEState_867", "BEState_868", "BEState_869", "BEState_870",
            "BEState_871", "BEState_872", "BEState_873", "BEState_874", "BEState_875", "BEState_876", "BEState_877", "BEState_878", "BEState_879", "BEState_880",
            "BEState_881", "BEState_882", "BEState_883", "BEState_884", "BEState_885", "BEState_886", "BEState_887", "BEState_888", "BEState_889", "BEState_890",
            "BEState_891", "BEState_892", "BEState_893", "BEState_894", "BEState_895", "BEState_896", "BEState_897", "BEState_898", "BEState_899", "BEState_900",
            "BEState_901", "BEState_902", "BEState_903", "BEState_904", "BEState_905", "BEState_906", "BEState_907", "BEState_908", "BEState_909", "BEState_910",
            "BEState_911", "BEState_912", "BEState_913", "BEState_914", "BEState_915", "BEState_916", "BEState_917", "BEState_918", "BEState_919", "BEState_920",
            "BEState_921", "BEState_922", "BEState_923", "BEState_924", "BEState_925", "BEState_926", "BEState_927", "BEState_928", "BEState_929", "BEState_930",
            "BEState_931", "BEState_932", "BEState_933", "BEState_934", "BEState_935", "BEState_936", "BEState_937", "BEState_938", "BEState_939", "BEState_940",
            "BEState_941", "BEState_942", "BEState_943", "BEState_944", "BEState_945", "BEState_946", "BEState_947", "BEState_948", "BEState_949", "BEState_950",
            "BEState_951", "BEState_952", "BEState_953", "BEState_954", "BEState_955", "BEState_956", "BEState_957", "BEState_958", "BEState_959", "BEState_960",
            "BEState_961", "BEState_962", "BEState_963", "BEState_964", "BEState_965", "BEState_966", "BEState_967", "BEState_968", "BEState_969", "BEState_970",
            "BEState_971", "BEState_972", "BEState_973", "BEState_974", "BEState_975", "BEState_976", "BEState_977", "BEState_978", "BEState_979", "BEState_980",
            "BEState_981", "BEState_982", "BEState_983", "BEState_984", "BEState_985", "BEState_986", "BEState_987", "BEState_988", "BEState_989", "BEState_990",
            "BEState_991", "BEState_992", "BEState_993", "BEState_994", "BEState_995", "BEState_996", "BEState_997", "BEState_998", "BEState_999", "BEState_1000",
            "BEState_1001", "BEState_1002", "BEState_1003", "BEState_1004", "BEState_1005", "BEState_1006", "BEState_1007", "BEState_1008", "BEState_1009", "BEState_1010",
            "BEState_1011", "BEState_1012", "BEState_1013", "BEState_1014", "BEState_1015", "BEState_1016", "BEState_1017", "BEState_1018", "BEState_1019", "BEState_1020",
            "BEState_1021", "BEState_1022", "BEState_1023", "BEState_1024", "BEState_1025", "BEState_1026", "BEState_1027", "BEState_1028", "BEState_1029", "BEState_1030",
            "BEState_1031", "BEState_1032", "BEState_1033", "BEState_1034", "BEState_1035", "BEState_1036", "BEState_1037", "BEState_1038", "BEState_1039", "BEState_1040",
            "BEState_1041", "BEState_1042", "BEState_1043", "BEState_1044", "BEState_1045", "BEState_1046", "BEState_1047", "BEState_1048", "BEState_1049", "BEState_1050",
            "BEState_1051", "BEState_1052", "BEState_1053", "BEState_1054", "BEState_1055", "BEState_1056", "BEState_1057", "BEState_1058", "BEState_1059", "BEState_1060",
            "BEState_1061", "BEState_1062", "BEState_1063", "BEState_1064", "BEState_1065", "BEState_1066", "BEState_1067", "BEState_1068", "BEState_1069", "BEState_1070",
            "BEState_1071", "BEState_1072", "BEState_1073", "BEState_1074", "BEState_1075", "BEState_1076", "BEState_1077", "BEState_1078", "BEState_1079", "BEState_1080",
            "BEState_1081", "BEState_1082", "BEState_1083", "BEState_1084", "BEState_1085", "BEState_1086", "BEState_1087", "BEState_1088", "BEState_1089", "BEState_1090",
            "BEState_1091", "BEState_1092", "BEState_1093", "BEState_1094", "BEState_1095", "BEState_1096", "BEState_1097", "BEState_1098", "BEState_1099", "BEState_1100",
            "BEState_1101", "BEState_1102", "BEState_1103", "BEState_1104", "BEState_1105", "BEState_1106", "BEState_1107", "BEState_1108", "BEState_1109", "BEState_1110",
            "BEState_1111", "BEState_1112", "BEState_1113", "BEState_1114", "BEState_1115", "BEState_1116", "BEState_1117", "BEState_1118", "BEState_1119", "BEState_1120",
            "BEState_1121", "BEState_1122", "BEState_1123", "BEState_1124", "BEState_1125", "BEState_1126", "BEState_1127", "BEState_1128", "BEState_1129", "BEState_1130",
            "BEState_1131", "BEState_1132", "BEState_1133", "BEState_1134", "BEState_1135", "BEState_1136", "BEState_1137", "BEState_1138", "BEState_1139", "BEState_1140",
            "BEState_1141", "BEState_1142", "BEState_1143", "BEState_1144", "BEState_1145", "BEState_1146", "BEState_1147", "BEState_1148", "BEState_1149", "BEState_1150",
            "BEState_1151", "BEState_1152", "BEState_1153", "BEState_1154", "BEState_1155", "BEState_1156", "BEState_1157", "BEState_1158", "BEState_1159", "BEState_1160",
            "BEState_1161", "BEState_1162", "BEState_1163", "BEState_1164", "BEState_1165", "BEState_1166", "BEState_1167", "BEState_1168", "BEState_1169", "BEState_1170",
            "BEState_1171", "BEState_1172", "BEState_1173", "BEState_1174", "BEState_1175", "BEState_1176", "BEState_1177", "BEState_1178", "BEState_1179", "BEState_1180",
            "BEState_1181", "BEState_1182", "BEState_1183", "BEState_1184", "BEState_1185", "BEState_1186", "BEState_1187", "BEState_1188", "BEState_1189", "BEState_1190",
            "BEState_1191", "BEState_1192", "BEState_1193", "BEState_1194", "BEState_1195", "BEState_1196", "BEState_1197", "BEState_1198", "BEState_1199", "BEState_1200",
            "BEState_1201", "BEState_1202", "BEState_1203", "BEState_1204", "BEState_1205", "BEState_1206", "BEState_1207", "BEState_1208", "BEState_1209", "BEState_1210",
            "BEState_1211", "BEState_1212", "BEState_1213", "BEState_1214", "BEState_1215", "BEState_1216", "BEState_1217", "BEState_1218", "BEState_1219", "BEState_1220",
            "BEState_1221", "BEState_1222", "BEState_1223", "BEState_1224", "BEState_1225", "BEState_1226", "BEState_1227", "BEState_1228", "BEState_1229", "BEState_1230",
            "BEState_1231", "BEState_1232", "BEState_1233", "BEState_1234", "BEState_1235", "BEState_1236", "BEState_1237", "BEState_1238", "BEState_1239", "BEState_1240",
            "BEState_1241", "BEState_1242", "BEState_1243", "BEState_1244", "BEState_1245", "BEState_1246", "BEState_1247", "BEState_1248", "BEState_1249", "BEState_1250",
            "BEState_1251", "BEState_1252", "BEState_1253", "BEState_1254", "BEState_1255", "BEState_1256", "BEState_1257", "BEState_1258", "BEState_1259", "BEState_1260",
            "BEState_1261", "BEState_1262", "BEState_1263", "BEState_1264", "BEState_1265", "BEState_1266", "BEState_1267", "BEState_1268", "BEState_1269", "BEState_1270",
            "BEState_1271", "BEState_1272", "BEState_1273", "BEState_1274", "BEState_1275", "BEState_1276", "BEState_1277", "BEState_1278", "BEState_1279", "BEState_1280",
            "BEState_1281", "BEState_1282", "BEState_1283", "BEState_1284", "BEState_1285", "BEState_1286", "BEState_1287", "BEState_1288", "BEState_1289", "BEState_1290",
            "BEState_1291", "BEState_1292", "BEState_1293", "BEState_1294", "BEState_1295", "BEState_1296", "BEState_1297", "BEState_1298", "BEState_1299", "BEState_1300",
            "BEState_1301", "BEState_1302", "BEState_1303", "BEState_1304", "BEState_1305", "BEState_1306", "BEState_1307", "BEState_1308", "BEState_1309", "BEState_1310",
            "BEState_1311", "BEState_1312", "BEState_1313", "BEState_1314", "BEState_1315", "BEState_1316", "BEState_1317", "BEState_1318", "BEState_1319", "BEState_1320",
            "BEState_1321", "BEState_1322", "BEState_1323", "BEState_1324", "BEState_1325", "BEState_1326", "BEState_1327", "BEState_1328", "BEState_1329", "BEState_1330",
            "BEState_1331", "BEState_1332", "BEState_1333", "BEState_1334", "BEState_1335", "BEState_1336", "BEState_1337", "BEState_1338", "BEState_1339", "BEState_1340",
            "BEState_1341", "BEState_1342", "BEState_1343", "BEState_1344", "BEState_1345", "BEState_1346", "BEState_1347", "BEState_1348", "BEState_1349", "BEState_1350",
            "BEState_1351", "BEState_1352", "BEState_1353", "BEState_1354", "BEState_1355", "BEState_1356", "BEState_1357", "BEState_1358", "BEState_1359", "BEState_1360",
            "BEState_1361", "BEState_1362", "BEState_1363", "BEState_1364", "BEState_1365", "BEState_1366", "BEState_1367", "BEState_1368", "BEState_1369", "BEState_1370",
            "BEState_1371", "BEState_1372", "BEState_1373", "BEState_1374", "BEState_1375", "BEState_1376", "BEState_1377", "BEState_1378", "BEState_1379", "BEState_1380",
            "BEState_1381", "BEState_1382", "BEState_1383", "BEState_1384", "BEState_1385", "BEState_1386", "BEState_1387", "BEState_1388", "BEState_1389", "BEState_1390",
            "BEState_1391", "BEState_1392", "BEState_1393", "BEState_1394", "BEState_1395", "BEState_1396", "BEState_1397", "BEState_1398", "BEState_1399", "BEState_1400",
            "BEState_1401", "BEState_1402", "BEState_1403", "BEState_1404", "BEState_1405", "BEState_1406", "BEState_1407", "BEState_1408", "BEState_1409", "BEState_1410",
            "BEState_1411", "BEState_1412", "BEState_1413", "BEState_1414", "BEState_1415", "BEState_1416", "BEState_1417", "BEState_1418", "BEState_1419", "BEState_1420",
            "BEState_1421", "BEState_1422", "BEState_1423", "BEState_1424", "BEState_1425", "BEState_1426", "BEState_1427", "BEState_1428", "BEState_1429", "BEState_1430",
            "BEState_1431", "BEState_1432", "BEState_1433", "BEState_1434", "BEState_1435", "BEState_1436", "BEState_1437", "BEState_1438", "BEState_1439", "BEState_1440",
            "BEState_1441", "BEState_1442", "BEState_1443", "BEState_1444", "BEState_1445", "BEState_1446", "BEState_1447", "BEState_1448", "BEState_1449", "BEState_1450",
            "BEState_1451", "BEState_1452", "BEState_1453", "BEState_1454", "BEState_1455", "BEState_1456", "BEState_1457", "BEState_1458", "BEState_1459", "BEState_1460"
        });
        // Abrahams (6)
        names.insert(names.end(), {"A", "B", "S", "E", "L", "V"});
        // ANMat (25)
        names.insert(names.end(), {
            "ANMat_1", "ANMat_2", "ANMat_3", "ANMat_4", "ANMat_5", "ANMat_6", "ANMat_7", "ANMat_8", "ANMat_9", "ANMat_10",
            "ANMat_11", "ANMat_12", "ANMat_13", "ANMat_14", "ANMat_15", "ANMat_16", "ANMat_17", "ANMat_18", "ANMat_19", "ANMat_20",
            "ANMat_21", "ANMat_22", "ANMat_23", "ANMat_24", "ANMat_25"
        });
        // ASMat (20)
        names.insert(names.end(), {
            "ASMat_1", "ASMat_2", "ASMat_3", "ASMat_4", "ASMat_5", "ASMat_6", "ASMat_7", "ASMat_8", "ASMat_9", "ASMat_10",
            "ASMat_11", "ASMat_12", "ASMat_13", "ASMat_14", "ASMat_15", "ASMat_16", "ASMat_17", "ASMat_18", "ASMat_19", "ASMat_20"
        });
        // AZMat (15)
        names.insert(names.end(), {
            "AZMat_1", "AZMat_2", "AZMat_3", "AZMat_4", "AZMat_5", "AZMat_6", "AZMat_7", "AZMat_8", "AZMat_9", "AZMat_10",
            "AZMat_11", "AZMat_12", "AZMat_13", "AZMat_14", "AZMat_15"
        });
        // DSMat (20)
        names.insert(names.end(), {
            "DSMat_1", "DSMat_2", "DSMat_3", "DSMat_4", "DSMat_5", "DSMat_6", "DSMat_7", "DSMat_8", "DSMat_9", "DSMat_10",
            "DSMat_11", "DSMat_12", "DSMat_13", "DSMat_14", "DSMat_15", "DSMat_16", "DSMat_17", "DSMat_18", "DSMat_19", "DSMat_20"
        });
        // DN2Mat (20)
        names.insert(names.end(), {
            "DN2Mat_1", "DN2Mat_2", "DN2Mat_3", "DN2Mat_4", "DN2Mat_5", "DN2Mat_6", "DN2Mat_7", "DN2Mat_8", "DN2Mat_9", "DN2Mat_10",
            "DN2Mat_11", "DN2Mat_12", "DN2Mat_13", "DN2Mat_14", "DN2Mat_15", "DN2Mat_16", "DN2Mat_17", "DN2Mat_18", "DN2Mat_19", "DN2Mat_20"
        });
        // Frags (215)
        names.insert(names.end(), {
            "Frags_1", "Frags_2", "Frags_3", "Frags_4", "Frags_5", "Frags_6", "Frags_7", "Frags_8", "Frags_9", "Frags_10",
            "Frags_11", "Frags_12", "Frags_13", "Frags_14", "Frags_15", "Frags_16", "Frags_17", "Frags_18", "Frags_19", "Frags_20",
            "Frags_21", "Frags_22", "Frags_23", "Frags_24", "Frags_25", "Frags_26", "Frags_27", "Frags_28", "Frags_29", "Frags_30",
            "Frags_31", "Frags_32", "Frags_33", "Frags_34", "Frags_35", "Frags_36", "Frags_37", "Frags_38", "Frags_39", "Frags_40",
            "Frags_41", "Frags_42", "Frags_43", "Frags_44", "Frags_45", "Frags_46", "Frags_47", "Frags_48", "Frags_49", "Frags_50",
            "Frags_51", "Frags_52", "Frags_53", "Frags_54", "Frags_55", "Frags_56", "Frags_57", "Frags_58", "Frags_59", "Frags_60",
            "Frags_61", "Frags_62", "Frags_63", "Frags_64", "Frags_65", "Frags_66", "Frags_67", "Frags_68", "Frags_69", "Frags_70",
            "Frags_71", "Frags_72", "Frags_73", "Frags_74", "Frags_75", "Frags_76", "Frags_77", "Frags_78", "Frags_79", "Frags_80",
            "Frags_81", "Frags_82", "Frags_83", "Frags_84", "Frags_85", "Frags_86", "Frags_87", "Frags_88", "Frags_89", "Frags_90",
            "Frags_91", "Frags_92", "Frags_93", "Frags_94", "Frags_95", "Frags_96", "Frags_97", "Frags_98", "Frags_99", "Frags_100",
            "Frags_101", "Frags_102", "Frags_103", "Frags_104", "Frags_105", "Frags_106", "Frags_107", "Frags_108", "Frags_109", "Frags_110",
            "Frags_111", "Frags_112", "Frags_113", "Frags_114", "Frags_115", "Frags_116", "Frags_117", "Frags_118", "Frags_119", "Frags_120",
            "Frags_121", "Frags_122", "Frags_123", "Frags_124", "Frags_125", "Frags_126", "Frags_127", "Frags_128", "Frags_129", "Frags_130",
            "Frags_131", "Frags_132", "Frags_133", "Frags_134", "Frags_135", "Frags_136", "Frags_137", "Frags_138", "Frags_139", "Frags_140",
            "Frags_141", "Frags_142", "Frags_143", "Frags_144", "Frags_145", "Frags_146", "Frags_147", "Frags_148", "Frags_149", "Frags_150",
            "Frags_151", "Frags_152", "Frags_153", "Frags_154", "Frags_155", "Frags_156", "Frags_157", "Frags_158", "Frags_159", "Frags_160",
            "Frags_161", "Frags_162", "Frags_163", "Frags_164", "Frags_165", "Frags_166", "Frags_167", "Frags_168", "Frags_169", "Frags_170",
            "Frags_171", "Frags_172", "Frags_173", "Frags_174", "Frags_175", "Frags_176", "Frags_177", "Frags_178", "Frags_179", "Frags_180",
            "Frags_181", "Frags_182", "Frags_183", "Frags_184", "Frags_185", "Frags_186", "Frags_187", "Frags_188", "Frags_189", "Frags_190",
            "Frags_191", "Frags_192", "Frags_193", "Frags_194", "Frags_195", "Frags_196", "Frags_197", "Frags_198", "Frags_199", "Frags_200",
            "Frags_201", "Frags_202", "Frags_203", "Frags_204", "Frags_205", "Frags_206", "Frags_207", "Frags_208", "Frags_209", "Frags_210",
            "Frags_211", "Frags_212", "Frags_213", "Frags_214", "Frags_215"
        });
        // AddFeatures (7)
        names.insert(names.end(), {"nBridgedBonds", "nHydroxylPrimary", "nHydroxylSecondary", "nHydroxylTertiary", "isPolyAcid", "isPolyAlcohol", "nEndocyclicSingleBonds"});
        
        return names;
    }

//// new version faster but not correct!!!



std::vector<std::vector<int>> calculateTopologicalMatrix(const RDKit::ROMol& mol) {
    int numAtoms = mol.getNumAtoms();
    std::vector<std::vector<int>> distanceMatrix(numAtoms, std::vector<int>(numAtoms, -1));

    for (int i = 0; i < numAtoms; ++i) {
        distanceMatrix[i][i] = 0;  // Distance to self is always zero
        for (int j = i + 1; j < numAtoms; ++j) {
            auto path = RDKit::MolOps::getShortestPath(mol, i, j);
            distanceMatrix[i][j] = path.size();
            distanceMatrix[j][i] = distanceMatrix[i][j];  // Symmetric assignment
        }
    }
    return distanceMatrix;
}




// Function to calculate atomic properties
void calculateAtomicDescriptors(const RDKit::ROMol &mol, std::vector<double> &alpha, std::vector<double> &alphaR, std::vector<double> &epsilon,  double &Alpha_P, double &Alpha_Y, double &Alpha_X) {
    int numAtoms = mol.getNumAtoms();
    alpha.resize(numAtoms, 0.0);
    alphaR.resize(numAtoms, 0.5);
    epsilon.resize(numAtoms, 0.3);
    const RDKit::PeriodicTable* tbl = RDKit::PeriodicTable::getTable();
    Alpha_P = 0.0;
    Alpha_Y = 0.0;
    Alpha_X = 0.0;

    for (int i = 0; i < numAtoms; ++i) {
        const auto atom = mol.getAtomWithIdx(i);
        int atomicNum = atom->getAtomicNum();

        if (atomicNum != 1) {
            int Z = atomicNum;
            int Zv = tbl->getNouterElecs(Z);
            int period = GetPrincipalQuantumNumber(Z);

            alpha[i] = ((Z - Zv) / static_cast<double>(Zv)) * (1.0 / (period - 1));
            epsilon[i] = -alpha[i] + 0.3 * Zv;

            int nonHNeighbors = 0;
            for (const auto &bond : mol.atomBonds(atom)) {
                const auto neighbor = bond->getOtherAtom(atom);
                if (neighbor->getAtomicNum() != 1) {
                    nonHNeighbors++;
                }
            }

            if (nonHNeighbors == 1) {
                Alpha_P += alpha[i];
            } else if (nonHNeighbors == 3) {
                Alpha_Y += alpha[i];
            } else if (nonHNeighbors == 4) {
                Alpha_X += alpha[i];
            }

        }
    }
}

void computeVEMContributions(const RDKit::ROMol &mol, const std::vector<double>& epsilon,
                             std::vector<double>& betaS, std::vector<double>& betaNS,
                             std::vector<double>& betaNSd, std::vector<double>& beta,
                             std::vector<double>& gamma, std::vector<double>& gammaRef,
                             const std::vector<double>& alpha, const std::vector<double>& alphaR) {

    int numAtoms = mol.getNumAtoms();
    betaS.resize(numAtoms, 0.0);
    betaNS.resize(numAtoms, 0.0);
    betaNSd.resize(numAtoms, 0.0);
    beta.resize(numAtoms, 0.0);
    gamma.resize(numAtoms, 0.0);
    gammaRef.resize(numAtoms, 0.0);

    const RDKit::PeriodicTable* tbl = RDKit::PeriodicTable::getTable();

    for (int i = 0; i < numAtoms; ++i) {
        const auto atom = mol.getAtomWithIdx(i);
        betaS[i] = 0.0;
        betaNS[i] = 0.0;
        betaNSd[i] = 0.0;
        double betaR = 0.0;
        bool aromatic = false;
        double nonconjugatedSum = 0.0;
        bool conjugated = false;
        bool hasConjugatedSingleBond = false;
        bool hasDoubleTripleBond = false;
        int conjugatedMultiplier = 1;
        bool delta = false;
        for (const auto &bond : mol.atomBonds(atom)) {
            const auto neighbor = bond->getOtherAtom(atom);
            if (neighbor->getAtomicNum() == 1) continue;
            betaR += 0.5;
            int j = neighbor->getIdx();
            double epsilonDiff = std::abs(epsilon[i] - epsilon[j]);
            if (epsilonDiff <= 0.3) {
                betaS[i] += 0.5/2;  // error in code have to divide by 2??? look like count twice ???
            } else {
                betaS[i] += 0.75/2; // error in code have to divide by 2??? look like count twice ???
            }

            if (bond->getIsAromatic()) {
                aromatic = true;
            } else if (bond->getBondType() == RDKit::Bond::SINGLE) {
                if (!conjugated) {
                    if (neighbor->getAtomicNum() == 6) {  // Check if neighbor is carbon
                        for (const auto& nbrBond : mol.atomBonds(neighbor)) {
                            if (nbrBond->getBondType() == RDKit::Bond::DOUBLE || nbrBond->getIsAromatic()) {
                                hasConjugatedSingleBond = true;
                                break;
                            } else if (nbrBond->getBondType() == RDKit::Bond::TRIPLE) {
                                conjugatedMultiplier = 2;
                                hasConjugatedSingleBond = true;
                                break;
                            }
                        }
                    }
                    else if (tbl->getNouterElecs(neighbor->getAtomicNum()) - neighbor->getFormalCharge() - mol.getAtomDegree(neighbor) >= 2) {
                        hasConjugatedSingleBond = true;
                }
                    if (hasDoubleTripleBond && hasConjugatedSingleBond) conjugated = true;
            }
            } else if (bond->getBondType() == RDKit::Bond::DOUBLE || bond->getBondType() == RDKit::Bond::TRIPLE) {
                int multiplier = (bond->getBondType() == RDKit::Bond::DOUBLE) ? 1 : 2;
                conjugatedMultiplier = multiplier;

                hasDoubleTripleBond = true;
                if (!conjugated) {
                    if (hasConjugatedSingleBond) {
                        conjugated = true;
                    } else if (neighbor->getAtomicNum() == 6) {
                        for (const auto& bond2 : mol.atomBonds(neighbor)) {
                            if (bond2->getBondType() == RDKit::Bond::SINGLE) {
                                const auto atom3 = bond2->getOtherAtom(neighbor);
                                if (atom3->getAtomicNum() == 6) {
                                    for (const auto& bond3 : mol.atomBonds(atom3)) {
                                        if (bond3->getBondType() == RDKit::Bond::DOUBLE || bond3->getIsAromatic()) {
                                            hasConjugatedSingleBond = true;
                                            break;
                                        } else if (bond3->getBondType() == RDKit::Bond::TRIPLE) {
                                            conjugatedMultiplier = 2;
                                            hasConjugatedSingleBond = true;
                                            break;
                                        }
                                    }
                                } else if (tbl->getNouterElecs(atom3->getAtomicNum()) - atom3->getFormalCharge() - mol.getAtomDegree(atom3) >= 2) {
                                    hasConjugatedSingleBond = true;
                                }
                                if (hasConjugatedSingleBond) {
                                    conjugated = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                double g = (epsilonDiff <= 0.3) ? 1.0 : 1.5;
                nonconjugatedSum += g * multiplier;
            }
            if (tbl->getNouterElecs(atom->getAtomicNum()) - atom->getFormalCharge() - mol.getAtomDegree(atom) >= 2 &&
                !atom->getIsAromatic() &&
                !isAtomInRing(*atom) &&
                neighbor->getIsAromatic() &&
                isAtomInRing(*neighbor)) {
                delta = true;
            }
        }

        betaNS[i] += (aromatic ? 2.0 : 0.0) + (conjugated ? 1.5 * conjugatedMultiplier : nonconjugatedSum) + (delta ? 0.5 : 0.0);
        betaNSd[i] = (delta ? 0.5 : 0.0);
        beta[i] = betaS[i] + betaNS[i];
        // Compute gamma values
        gamma[i] = alpha[i] / beta[i];
        gammaRef[i] = alphaR[i] / betaR;
    }
}

// Main function to compute ETA descriptors
std::vector<double> calculateETADescriptors(const RDKit::ROMol& mol) {
    int numAtoms = mol.getNumAtoms();

    std::unique_ptr<RDKit::RWMol> kekulizedMol(new RDKit::RWMol(mol));
    RDKit::MolOps::Kekulize(*kekulizedMol, false);


    // Step 1: Calculate atomic descriptors
    std::vector<double> alpha, alphaR, epsilon;
    double Alpha_P, Alpha_Y, Alpha_X;
    calculateAtomicDescriptors(*kekulizedMol, alpha, alphaR, epsilon, Alpha_P, Alpha_Y, Alpha_X);

    // Step 2: Compute VEM contributions fixed by /2 the values why mistary ????
    std::vector<double> betaS, betaNS, betaNSd, beta, gamma, gammaRef;
    computeVEMContributions(*kekulizedMol, epsilon, betaS, betaNS, betaNSd, beta, gamma, gammaRef, alpha, alphaR);

    // Step 3: Calculate topological matrix (distance matrix)
    auto pathLengths = calculateTopologicalMatrix(mol);

    // Step 4: Count rings
    int maxRings = RDKit::MolOps::findSSSR(mol);

    // Step 5: Compute descriptor sums
    double alphaSum = std::accumulate(alpha.begin(), alpha.end(), 0.0);
    double alphaRSum = std::accumulate(alphaR.begin(), alphaR.end(), 0.0);
    double epsilonSum = std::accumulate(epsilon.begin(), epsilon.end(), 0.0);
    double betaSSum = std::accumulate(betaS.begin(), betaS.end(), 0.0);
    double betaNSSum = std::accumulate(betaNS.begin(), betaNS.end(), 0.0);
    double betaNSdSum = std::accumulate(betaNSd.begin(), betaNSd.end(), 0.0);

    // Step 6: Compute various ETA values based on atomic and bond properties
    double etaSum = 0.0, etaRSum = 0.0, etaLSum = 0.0, etaRLSum = 0.0;

    for (int i = 0; i < numAtoms; ++i) {
        for (int j = i + 1; j < numAtoms; ++j) {
            if (pathLengths[i][j] > 0) {
                etaSum += std::sqrt(gamma[i] * gamma[j] / (pathLengths[i][j] * pathLengths[i][j]));
                etaRSum += std::sqrt(gammaRef[i] * gammaRef[j] / (pathLengths[i][j] * pathLengths[i][j]));

                if (pathLengths[i][j] == 1) {
                    etaLSum += std::sqrt(gamma[i] * gamma[j]);
                    etaRLSum += std::sqrt(gammaRef[i] * gammaRef[j]);
                }
            }
        }
    }



    // Step 7: Compute ETA indices

    int N = numAtoms; // Total atoms count
    int Nv = numAtoms - std::count_if(mol.atoms().begin(), mol.atoms().end(), [](const RDKit::Atom* a){ return a->getAtomicNum() == 1; });
    int Nr = 5 * N; // Placeholder for number of reference hydrogens, adjust as needed
    int Nss = N;    // Total number of saturated atoms
    int Nxh = 0;    // Count heteroatom-bound hydrogens (modify logic if needed)
    double epsilonEHSum = epsilonSum;
    double epsilonRSum = epsilonSum;
    double epsilonSSSum = epsilonSum;
    double epsilonXHSum = 0.0;  // Placeholder for heteroatom-bonded H sum

    double PsiTemp = 0.71429;  // Fixed constant


    // EtaCoreCount
    double ETA_Alpha = alphaSum;
    double ETA_AlphaP = alphaSum / Nv;

    // EtaShapeIndex
    double ETA_Shape_P = Alpha_P / ETA_Alpha;
    double ETA_Shape_Y = Alpha_Y / ETA_Alpha;
    double ETA_Shape_X = Alpha_X / ETA_Alpha;

    // VEMCount
    double ETA_Beta = betaSSum + betaNSSum;
    double ETA_BetaP = ETA_Beta / Nv;
    double ETA_Beta_s = betaSSum;
    double ETA_BetaP_s = betaSSum / Nv;
    double ETA_Beta_ns = betaNSSum;
    double ETA_BetaP_ns = betaNSSum / Nv;
    double ETA_Beta_ns_d = betaNSdSum;
    double ETA_BetaP_ns_d = ETA_Beta_ns_d / Nv;

    //  "EtaCompositeIndex"
    double ETA_Eta = etaSum;
    double ETA_EtaP = etaSum / Nv;
    double ETA_Eta_L = etaLSum;
    double ETA_EtaP_L = etaLSum / Nv;
    double ETA_Eta_R = etaRSum;
    double ETA_EtaP_R = etaRSum / Nv;
    double ETA_Eta_R_L = etaRLSum;
    double ETA_EtaP_R_L = etaRLSum / Nv;

    // EtaFunctionalityIndex
    double ETA_Eta_F = etaRSum - etaSum;
    double ETA_EtaP_F = ETA_Eta_F / Nv;
    double ETA_Eta_F_L = etaRLSum - etaLSum;
    double ETA_EtaP_F_L = ETA_Eta_F_L / Nv;

    // EtaBranchingIndex
    double ETA_Eta_B = (Nv > 3 ? (1.414 + (Nv - 3) * 0.5) - etaRLSum : 0);
    double ETA_EtaP_B = ETA_Eta_B / Nv;
    double ETA_Eta_B_RC = ETA_Eta_B + 0.086 * maxRings;
    double ETA_EtaP_B_RC = ETA_Eta_B_RC / Nv;

    //"EtaDeltaAlpha" :  dAlpha_A, dAlpha_B
    double ETA_dAlpha_A = std::max((alphaSum - alphaRSum) / Nv, 0.0);
    double ETA_dAlpha_B = std::max((alphaRSum - alphaSum) / Nv, 0.0);

    // "EtaEpsilon":
    double ETA_Epsilon_1 = epsilonSum / N;
    double ETA_Epsilon_2 = epsilonEHSum / Nv;
    double ETA_Epsilon_3 = epsilonRSum / Nr;
    double ETA_Epsilon_4 = epsilonSSSum / Nss;
    double ETA_Epsilon_5 = (epsilonEHSum + epsilonXHSum) / (Nv + Nxh);

    // "EtaDeltaEpsilon", A,B,C,D use EtaEps diff cases
    double ETA_dEpsilon_A = ETA_Epsilon_1 - ETA_Epsilon_3;
    double ETA_dEpsilon_B = ETA_Epsilon_1 - ETA_Epsilon_4;
    double ETA_dEpsilon_C = ETA_Epsilon_3 - ETA_Epsilon_4;
    double ETA_dEpsilon_D = ETA_Epsilon_2 - ETA_Epsilon_5;

    //"EtaDeltaBeta" 2 Values
    double ETA_dBeta = betaNSSum - betaSSum;
    double ETA_dBetaP = ETA_dBeta / Nv;

    // EtaPsi
    double ETA_Psi_1 = alphaSum / epsilonEHSum;

    //EtaDeltaPsi
    double ETA_dPsi_A = std::max(PsiTemp - ETA_Psi_1, 0.0);
    double ETA_dPsi_B = std::max(ETA_Psi_1 - PsiTemp, 0.0);
    // Step 8: Return all calculated values
    return {
        ETA_Alpha, ETA_AlphaP,
        ETA_Shape_P, ETA_Shape_Y, ETA_Shape_X,
        ETA_Beta, ETA_BetaP, ETA_Beta_s, ETA_BetaP_s, ETA_Beta_ns, ETA_BetaP_ns, ETA_Beta_ns_d, ETA_BetaP_ns_d,
        ETA_Eta, ETA_EtaP,ETA_Eta_L, ETA_EtaP_L, ETA_Eta_R,   ETA_EtaP_R,  ETA_Eta_R_L, ETA_EtaP_R_L,
        ETA_Eta_F, ETA_EtaP_F, ETA_Eta_F_L, ETA_EtaP_F_L,
        ETA_Eta_B, ETA_EtaP_B, ETA_Eta_B_RC, ETA_EtaP_B_RC,
        ETA_dAlpha_A, ETA_dAlpha_B,
        ETA_Epsilon_1, ETA_Epsilon_2, ETA_Epsilon_3, ETA_Epsilon_4, ETA_Epsilon_5,
        ETA_dEpsilon_A, ETA_dEpsilon_B, ETA_dEpsilon_C, ETA_dEpsilon_D,
        ETA_dBeta ,ETA_dBetaP,  ETA_Psi_1, ETA_dPsi_A, ETA_dPsi_B
    };

}


    // start ringcount
    std::vector<std::vector<int>>  GetSSSR(const RDKit::ROMol& mol) {
        RDKit::ROMol mol_copy(mol); // Create a non-const copy of the molecule
        std::vector<std::vector<int>> sssr;
        RDKit::MolOps::symmetrizeSSSR(mol_copy, sssr);
        return sssr;
    }

    // Get all rings in the molecule
    std::vector<std::set<int>> GetRings(const RDKit::ROMol&,std::vector<std::vector<int>> sssr ) {

        std::vector<std::set<int>> rings;
        for (const auto& ring : sssr) {
            rings.emplace_back(ring.begin(), ring.end());
        }
        return rings;
    }

    // Type aliases for readability
    using Ring = std::set<int>;
    using Rings = std::vector<Ring>;

    // Helper function to find connected components in an undirected graph
    std::vector<std::vector<int>> findConnectedComponents(const std::unordered_map<int, std::unordered_set<int>>& graph) {
        std::vector<std::vector<int>> components;
        std::unordered_set<int> visited;

        for (const auto& [node, _] : graph) {
            if (visited.count(node) == 0) {
                // Perform BFS to find all nodes in the current component
                std::vector<int> component;
                std::queue<int> toVisit;
                toVisit.push(node);
                visited.insert(node);

                while (!toVisit.empty()) {
                    int current = toVisit.front();
                    toVisit.pop();
                    component.push_back(current);

                    for (int neighbor : graph.at(current)) {
                        if (visited.count(neighbor) == 0) {
                            visited.insert(neighbor);
                            toVisit.push(neighbor);
                        }
                    }
                }

                components.push_back(component);
            }
        }

        return components;
    }


    // Function to calculate fused rings
    std::vector<std::vector<int>> calcFusedRings(const Rings& rings) {
        if (rings.size() < 2) {
            return {};
        }

        // Graph adjacency list representation
        std::unordered_map<int, std::unordered_set<int>> graph;

        size_t numRings = rings.size();
        for (size_t i = 0; i < numRings; ++i) {
            for (size_t j = i + 1; j < numRings; ++j) {
                std::vector<int> intersection;
                std::set_intersection(
                    rings[i].begin(), rings[i].end(),
                    rings[j].begin(), rings[j].end(),
                    std::back_inserter(intersection)
                );

                if (intersection.size() >= 2) {
                    graph[i].insert(j);
                    graph[j].insert(i);
                }
            }
        }

        // Find connected components in the graph
        std::vector<std::vector<int>> components = findConnectedComponents(graph);

        // Map connected ring indices to their atom sets
        std::vector<std::vector<int>> fusedRings;
        for (const auto& component : components) {
            std::unordered_set<int> fusedRingAtoms;
            for (int ringIdx : component) {
                fusedRingAtoms.insert(rings[ringIdx].begin(), rings[ringIdx].end());
            }

            fusedRings.emplace_back(fusedRingAtoms.begin(), fusedRingAtoms.end());
        }

        return fusedRings;
    }

    // Calculate ring count descriptors
    std::vector<int> calcRingCount(const ROMol& mol) {
        std::vector<int> descriptors(138, 0); // Placeholder for descriptor vector

        std::vector<std::vector<int>>  sssr = GetSSSR(mol);


        // part one make rings inventory

        descriptors[0] = sssr.size();
        // Iterate over all rings and count descriptors based on size, aromaticity, fused
        for (const auto& ring : sssr) {

            size_t ring_size = ring.size();

            bool is_aromatic = true;

            bool has_hetero = false;


            for (int atom_idx : ring) {
                const Atom* atom = mol.getAtomWithIdx(atom_idx);
                if (!atom->getIsAromatic()) {
                    is_aromatic = false;
                }
                if (atom->getAtomicNum() != 6) {
                    has_hetero = true;
                }
            }

            // Update ring counts for the respective features
            if (ring_size >= 3 && ring_size <= 12) {
                descriptors[ring_size - 3+1]++; // n3Ring to n12Ring
                    if (has_hetero ) {
                        descriptors[12]++;  // nHRing sum
                        descriptors[12 + (ring_size - 3)+1]++;
                    }

                    if (is_aromatic ) {
                        descriptors[24]++; // naRing

                        descriptors[24 + (ring_size - 3)+1]++; // n3aRing to n12aRing
                    }

                    if (is_aromatic && has_hetero) {
                        descriptors[36]++; // naHRing

                        descriptors[36 + (ring_size - 3)+1]++; // n3aHRing to n12aHRing
                    }

                    if (!is_aromatic) {
                        descriptors[48]++; // nARing
                        descriptors[48 + (ring_size - 3)+1]++; // n3ARing to n12ARing
                    }

                    if (!is_aromatic && has_hetero) {
                        descriptors[60]++; // nAHRing
                        descriptors[60 + (ring_size - 3)+1]++; // n3AHRing to n12AHRing
                    }
                }


            if (ring_size > 12) {

                // greater
                if (has_hetero) {
                    descriptors[12]++; // nHRing sum
                    descriptors[23]++; // nG12 HRing

                }
                if (is_aromatic) {
                    descriptors[24]++; // naRing  sum
                    descriptors[35]++; // nG12 aRing
                }

                if (is_aromatic && has_hetero) {
                    descriptors[36]++; // naHRing  sum
                    descriptors[47]++; // nG12 aHRing
                }

                if (!is_aromatic) {
                    descriptors[48]++; // nARing  sum
                    descriptors[59]++; // nG12 ARing
                }

                if (!is_aromatic &&  has_hetero) {
                    descriptors[60]++; // nAHRing  sum
                    descriptors[71]++; // nG12AHRing
                }
            }

        }

        // part two make fused ring inventory
        // maybe we don't need the rings step...
        auto rings = GetRings(mol, sssr);
        auto fusedRings = calcFusedRings(rings);


        descriptors[72] = fusedRings.size();

        for (const auto& fusedring : fusedRings) {

            size_t fused_ring_size = fusedring.size();

            bool is_aromatic = true;

            bool has_hetero = false;
            // we could also break it...
            for (int atom_idx : fusedring) {
                const Atom* atom = mol.getAtomWithIdx(atom_idx);
                if (!atom->getIsAromatic()) {
                    is_aromatic = false;
                }
                if (atom->getAtomicNum() != 6) {
                    has_hetero = true;
                }
            }

            // Update ring counts for the respective features
            if (fused_ring_size >= 4 && fused_ring_size <= 12) {
                descriptors[72 + fused_ring_size - 4+1]++; // n4FRing to n12FRing
                if (has_hetero  ) {
                    descriptors[83]++;  // nFHRing sum
                    descriptors[83 + (fused_ring_size - 4)+1]++;  //n4FHRing to n12FHRing
                }
                if (is_aromatic ) {
                    descriptors[94]++;  // nFHRing sum
                    descriptors[94 + (fused_ring_size - 4)+1]++;  //n4FaRing to n12FaRing
                }

                if (is_aromatic && has_hetero) {
                    descriptors[105]++; // naHRing

                    descriptors[105 + (fused_ring_size - 4)+1]++; // n3aHRing to n12aHRing
                }

                if (!is_aromatic ) {
                    descriptors[116]++; // nARing
                    descriptors[116 + (fused_ring_size - 4)+1]++; // n3ARing to n12ARing
                }

                if (!is_aromatic && has_hetero) {
                    descriptors[127]++; // nAHRing
                    descriptors[127 + (fused_ring_size - 4)+1]++; // n3AHRing to n12AHRing
                }
            }
            if (fused_ring_size > 12) {
                    // greater
                descriptors[82]++; // v3 fix: nG12FRing (plain) was never incremented -> always 0
                if (has_hetero) {
                    descriptors[83]++; // nFHRing sum
                    descriptors[93]++; // nG12FHRing
                }
                if (is_aromatic) {
                    descriptors[94]++; // nFaRing sum
                    descriptors[104]++; // nG12FaRing
                }

                if (is_aromatic && has_hetero) {
                    descriptors[105]++; // nFaHRing sum
                    descriptors[115]++; // nG12FaHRing
                }

                if (!is_aromatic) {
                    descriptors[116]++; // nFARing sum
                    descriptors[126]++; // nG12FARing
                }

                if (!is_aromatic &&  has_hetero) {
                    descriptors[127]++; // nFAHRing sum
                    descriptors[137]++; // nG12FAHRing
                }
            }
        }

        return descriptors;
    }






    static const std::vector<std::pair<std::string, std::string>>  hsPatterns = {
            {"HsOH", "[OD1H]-*"}, {"HdNH", "[ND1H]=*"},  {"HsSH", "[SD1H]-*"}, {"HsNH2", "[ND1H2]-*"},
            {"HssNH", "[ND2H](-*)-*"}, {"HaaNH", "[nD2H](:*):*"}, {"HsNH3p", "[ND1H3]-*"},{"HssNH2p", "[ND2H2](-*)-*"},
            {"HsssNHp", "[ND3H](-*)(-*)-*"},  {"HtCH", "[CD1H]#*"}, {"HdCH2", "[CD1H2]=*"},
            {"HdsCH", "[CD2H](=*)-*"}, {"HaaCH", "[#6D2H](:*):*"}, {"HCHnX", "[CX4;!H0]-[F,Cl,Br,I]"}, {"HCsats", "[CX4;!H0]~[!a]"},
            {"HCsatu", "[CX4;!H0]-[*]:,=,#[*]"}, {"HAvin", "[CX3H](=C)-[c]"}, {"Hall", "[*;!H0]"}, {"Hother", "[*;!H0]:,=,#[*]"},
            {"Hmisc", "[*;!H0]~[B,Si,P,Ge,As,Se,Sn,Pb]"}
    };



    // Define the EState atom types and their SMARTS patterns
    // NOTE (osmordred vs Mordred -- MINOR, intentionally left as-is): the
    // aromatic-N E-state types below (aaN / aaNH, and aasC) can disagree with
    // Mordred by ~0.05-0.11 on a single borderline heteroaromatic (Histidine's
    // imidazole), where one ring N straddles the aaN/aaNH boundary under RDKit
    // aromaticity. This is a tolerance-level atom-typing/tautomer effect, not a
    // clear bug -- not changed in v3.
    static const std::vector<std::pair<std::string, std::string>> esPatterns = {
    {"sLi", "[LiD1]-*"},{"ssBe", "[BeD2](-*)-*"},{"ssssBe", "[BeD4](-*)(-*)(-*)-*"},
    {"ssBH", "[BD2H](-*)-*"},{"sssB", "[BD3](-*)(-*)-*"},{"ssssB", "[BD4](-*)(-*)(-*)-*"},
    {"sCH3", "[CD1H3]-*"},{"dCH2", "[CD1H2]=*"},{"ssCH2", "[CD2H2](-*)-*"},{"tCH", "[CD1H]#*"},
    {"dsCH", "[CD2H](=*)-*"},{"aaCH", "[C,c;D2H](:*):*"},{"sssCH", "[CD3H](-*)(-*)-*"},{"ddC", "[CD2H0](=*)=*"},
    {"tsC", "[CD2H0](#*)-*"},{"dssC", "[CD3H0](=*)(-*)-*"},{"aasC", "[C,c;D3H0](:*)(:*)-*"},{"aaaC", "[C,c;D3H0](:*)(:*):*"},
    {"ssssC", "[CD4H0](-*)(-*)(-*)-*"},{"sNH3", "[ND1H3]-*"},{"sNH2", "[ND1H2]-*"},{"ssNH2", "[ND2H2](-*)-*"},{"dNH", "[ND1H]=*"},
    {"ssNH", "[ND2H](-*)-*"},{"aaNH", "[N,nD2H](:*):*"},{"tN", "[ND1H0]#*"},{"sssNH", "[ND3H](-*)(-*)-*"},{"dsN", "[ND2H0](=*)-*"},
    {"aaN", "[N,nD2H0](:*):*"},{"sssN", "[ND3H0](-*)(-*)-*"},{"ddsN", "[ND3H0](~[OD1H0])(~[OD1H0])-,:*"},
    {"aasN", "[N,nD3H0](:*)(:*)-,:*"}, {"ssssN", "[ND4H0](-*)(-*)(-*)-*"},{"sOH", "[OD1H]-*"},{"dO", "[OD1H0]=*"},{"ssO", "[OD2H0](-*)-*"},
    {"aaO", "[O,oD2H0](:*):*"},{"sF", "[FD1]-*"},{"sSiH3", "[SiD1H3]-*"},{"ssSiH2", "[SiD2H2](-*)-*"},{"sssSiH", "[SiD3H1](-*)(-*)-*"},
    {"ssssSi", "[SiD4H0](-*)(-*)(-*)-*"},{"sPH2", "[PD1H2]-*"},{"ssPH", "[PD2H1](-*)-*"},{"sssP", "[PD3H0](-*)(-*)-*"},
    {"dsssP", "[PD4H0](=*)(-*)(-*)-*"},{"sssssP", "[PD5H0](-*)(-*)(-*)(-*)-*"},{"sSH", "[SD1H1]-*"},{"dS", "[SD1H0]=*"},
    {"ssS", "[SD2H0](-*)-*"},{"aaS", "[S,sD2H0](:*):*"},{"dssS", "[SD3H0](=*)(-*)-*"},{"ddssS", "[SD4H0](~[OD1H0])(~[OD1H0])(-*)-*"},
    {"sCl", "[ClD1]-*"},{"sGeH3", "[GeD1H3](-*)"},{"ssGeH2", "[GeD2H2](-*)-*"},{"sssGeH", "[GeD3H1](-*)(-*)-*"},
    {"ssssGe", "[GeD4H0](-*)(-*)(-*)-*"},{"sAsH2", "[AsD1H2]-*"},{"ssAsH", "[AsD2H1](-*)-*"},{"sssAs", "[AsD3H0](-*)(-*)-*"},
    {"sssdAs", "[AsD4H0](=*)(-*)(-*)-*"},{"sssssAs", "[AsD5H0](-*)(-*)(-*)(-*)-*"},{"sSeH", "[SeD1H1]-*"},{"dSe", "[SeD1H0]=*"},
    {"ssSe", "[SeD2H0](-*)-*"},{"aaSe", "[SeD2H0](:*):*"},{"dssSe", "[SeD3H0](=*)(-*)-*"},{"ddssSe", "[SeD4H0](=*)(=*)(-*)-*"},
    {"sBr", "[BrD1]-*"},{"sSnH3", "[SnD1H3]-*"},{"ssSnH2", "[SnD2H2](-*)-*"},{"sssSnH", "[SnD3H1](-*)(-*)-*"},{"ssssSn", "[SnD4H0](-*)(-*)(-*)-*"},
    {"sI", "[ID1]-*"},{"sPbH3", "[PbD1H3]-*"},{"ssPbH2", "[PbD2H2](-*)-*"},{"sssPbH", "[PbD3H1](-*)(-*)-*"},{"ssssPb", "[PbD4H0](-*)(-*)(-*)-*"}
    };

    // Define the EState atom types and their SMARTS patterns
    static const std::vector<std::pair<std::string, std::string>> esPatternsFromOEState = {
    {"sLi", "[LiD1]-*"},{"ssBe", "[BeD2](-*)-*"},{"ssssBe", "[BeD4](-*)(-*)(-*)-*"},
    {"ssBH", "[BD2H](-*)-*"},{"sssB", "[BD3](-*)(-*)-*"},{"ssssB", "[BD4](-*)(-*)(-*)-*"},
    {"sCH3", "[CD1H3]-*"},{"dCH2", "[CD1H2]=*"},{"ssCH2", "[CD2H2](-*)-*"},{"tCH", "[CD1H]#*"},
    {"dsCH", "[CD2H](=*)-*"},{"aaCH", "[C,c;D2H](:*):*"},{"sssCH", "[CD3H](-*)(-*)-*"},{"ddC", "[CD2H0](=*)=*"},
    {"tsC", "[CD2H0](#*)-*"},{"dssC", "[CD3H0](=*)(-*)-*"},{"aasC", "[C,c;D3H0](:*)(:*)-*"},{"aaaC", "[C,c;D3H0](:*)(:*):*"},
    {"ssssC", "[CD4H0](-*)(-*)(-*)-*"},{"sNH3", "[ND1H3]-*"},{"sNH2", "[ND1H2]-*"},{"ssNH2", "[ND2H2](-*)-*"},{"dNH", "[ND1H]=*"},
    {"ssNH", "[ND2H](-*)-*"},{"aaNH", "[N,nD2H](:*):*"},{"tN", "[ND1H0]#*"},{"sssNH", "[ND3H](-*)(-*)-*"},{"dsN", "[ND2H0](=*)-*"},
    {"aaN", "[N,nD2H0](:*):*"},{"sssN", "[ND3H0](-*)(-*)-*"},{"ddsN", "[ND3H0](~[OD1H0])(~[OD1H0])-,:*"},
    {"aasN", "[N,nD3H0](:*)(:*)-,:*"}, {"ssssN", "[ND4H0](-*)(-*)(-*)-*"},{"sOH", "[OD1H]-*"},{"dO", "[OD1H0]=*"},{"ssO", "[OD2H0](-*)-*"},
    {"aaO", "[O,oD2H0](:*):*"},{"sF", "[FD1]-*"},{"sSiH3", "[SiD1H3]-*"},{"ssSiH2", "[SiD2H2](-*)-*"},{"sssSiH", "[SiD3H1](-*)(-*)-*"},
    {"ssssSi", "[SiD4H0](-*)(-*)(-*)-*"},{"sPH2", "[PD1H2]-*"},{"ssPH", "[PD2H1](-*)-*"},{"sssP", "[PD3H0](-*)(-*)-*"},
    {"dsssP", "[PD4H0](=*)(-*)(-*)-*"},{"sssssP", "[PD5H0](-*)(-*)(-*)(-*)-*"},{"sSH", "[SD1H1]-*"},{"dS", "[SD1H0]=*"},
    {"ssS", "[SD2H0](-*)-*"},{"aaS", "[S,sD2H0](:*):*"},{"dssS", "[SD3H0](=*)(-*)-*"},{"ddssS", "[SD4H0](~[OD1H0])(~[OD1H0])(-*)-*"},
    {"sCl", "[ClD1]-*"},{"sGeH3", "[GeD1H3](-*)"},{"ssGeH2", "[GeD2H2](-*)-*"},{"sssGeH", "[GeD3H1](-*)(-*)-*"},
    {"ssssGe", "[GeD4H0](-*)(-*)(-*)-*"},{"sAsH2", "[AsD1H2]-*"},{"ssAsH", "[AsD2H1](-*)-*"},{"sssAs", "[AsD3H0](-*)(-*)-*"},
    {"sssdAs", "[AsD4H0](=*)(-*)(-*)-*"},{"sssssAs", "[AsD5H0](-*)(-*)(-*)(-*)-*"},{"sSeH", "[SeD1H1]-*"},{"dSe", "[SeD1H0]=*"},
    {"ssSe", "[SeD2H0](-*)-*"},{"aaSe", "[SeD2H0](:*):*"},{"dssSe", "[SeD3H0](=*)(-*)-*"},{"ddssSe", "[SeD4H0](=*)(=*)(-*)-*"},
    {"sBr", "[BrD1]-*"},{"sSnH3", "[SnD1H3]-*"},{"ssSnH2", "[SnD2H2](-*)-*"},{"sssSnH", "[SnD3H1](-*)(-*)-*"},{"ssssSn", "[SnD4H0](-*)(-*)(-*)-*"},
    {"sI", "[ID1]-*"},{"sPbH3", "[PbD1H3]-*"},{"ssPbH2", "[PbD2H2](-*)-*"},{"sssPbH", "[PbD3H1](-*)(-*)-*"},{"ssssPb", "[PbD4H0](-*)(-*)(-*)-*"},
    {"sNH2(A)","[#7;D1;X3;H2][CX4;A]"},{"sNH2(a)","[#7;D1;X3;H2][c]"},    {"sNH2(oth)","[$([#7;D1;X3;H2][CX3])]"},
    {"ssNH(A)","[$([#7;X3;D2;H]([CX4;A][CX4;A]))]"},   {"ssNH(a)","[$([#7;X3;D2;H]([c;a])-*)]"},   {"ssNH(oth)","[$([#7;X3;D2;H]([CX3])-*)]"},
    {"sssN(A)","[$([#7;D3;X3;H0]([CX4;A])([CX4;A])[CX4;A])]"},   {"sssN(a)","[$([#7;D3;X3;H0]([c;a])(-*)-*)]"},   {"sssN(oth)","[$([#7D3;X3;H0]([CX3])(-*)-*)]"},
    {"ddsN(nitro)","[$([#7X3](=O)=O),$([#7X3+](=O)[O-])][!#8]"},   {"sOH(A)","[$([#8;X2;D1;H][CX4;A])]"},   {"sOH(a)","[$([#8;X2;D1;H][c;a])]"},
    {"sOH(acid)","[$([#8;X2;D1;H][CX3](=[OX1]))]"},   {"sOH(zwit)","[$([#8;X2;D1;H,OX1-][CX3](=[OX1])[NX3,NX4+])]"},
    {"ssO(ester)","[$([#8;X2;D2;H0]([CX3]=[OX1H0])[#6])]"},   {"dOfix","[#8;D1;X1;H0]~*"},   {"dO(keto)","[$([#8;X1;D1;H0]=[#6X3]([#6])[#6])]"},
    {"dO(acid)","[$([#8;X1;D1;H0]=[#6X3]([OX2H1]))]"},   {"dO(ester)","[$([#8;X1;D1;H0]=[#6X3]([OX2H0])[#6])]"},
    {"dO(amid)","[$([#8;X1;D1;H0]=[#6X3]([#7])[#6])]"},   {"dO(nitro)","[$([#8;X1;D1;H0]~[#7X3]~[#8;X1;D1;H0])]"},
    {"dO(sulfo)","[$([#8;X1]=[#16;X4]=[#8;X1]),$([#8;X1-][#16;X4+2][#8;X1-])]"}
    };


    // Function to add SMARTS queries safely to the vector of pairs
    void addToQueries(std::vector<std::pair<std::string, std::shared_ptr<RDKit::RWMol>>>& queries,
                      const std::string& key,
                      std::shared_ptr<RDKit::RWMol> mol) {
        for (auto& pair : queries) {
            if (pair.first == key) {  // Key already exists, update the value
                pair.second = mol;
                return;
            }
        }
        queries.emplace_back(key, mol);  // Key not found, add a new pair
    }


    // Precompile SMARTS patterns for efficiency
    static const std::vector<std::pair<std::string, std::shared_ptr<RDKit::RWMol>>>& GetHsQueries() {
      static const std::vector<
          std::pair<std::string, std::shared_ptr<RDKit::RWMol>>>
          hsQueries = [] {
            std::vector<std::pair<std::string, std::shared_ptr<RDKit::RWMol>>>
                queries;
            for (const auto &entry : hsPatterns) {
              auto mol = RDKit::SmartsToMol(entry.second);
              if (mol) {
                addToQueries(queries, entry.first,
                             std::shared_ptr<RDKit::RWMol>(mol));
              } else {
                std::cerr << "Invalid SMARTS: " << entry.second << std::endl;
              }
            }
            return queries;
          }();
      return hsQueries;
    }

    // Precompile SMARTS patterns for efficiency
    static const std::vector<std::pair<std::string, std::shared_ptr<RDKit::RWMol>>>& GetesQueries() {
      static const std::vector<
          std::pair<std::string, std::shared_ptr<RDKit::RWMol>>>
          esQueries = [] {
            std::vector<std::pair<std::string, std::shared_ptr<RDKit::RWMol>>>
                queries;
            for (const auto &entry : esPatterns) {
              auto mol = RDKit::SmartsToMol(entry.second);
              if (mol) {
                addToQueries(queries, entry.first,
                             std::shared_ptr<RDKit::RWMol>(mol));
              } else {
                std::cerr << "Invalid SMARTS: " << entry.second << std::endl;
              }
            }
            return queries;
          }();
      return esQueries;
    }

    // Precompile SMARTS patterns for efficiency
    static const std::vector<std::pair<std::string, std::shared_ptr<RDKit::RWMol>>>&
    GetesExtQueries() {
      static const std::vector<
          std::pair<std::string, std::shared_ptr<RDKit::RWMol>>>
          esExtQueries = [] {
            std::vector<std::pair<std::string, std::shared_ptr<RDKit::RWMol>>>
                queries;
            for (const auto &entry : esPatternsFromOEState) {
              auto mol = RDKit::SmartsToMol(entry.second);
              if (mol) {
                addToQueries(queries, entry.first,
                             std::shared_ptr<RDKit::RWMol>(mol));
              } else {
                std::cerr << "Invalid SMARTS: " << entry.second << std::endl;
              }
            }
            return queries;
          }();
      return esExtQueries;
    }


    // Function to switch between standard and extended SMARTS queries
    const std::vector<std::pair<std::string, std::shared_ptr<RDKit::RWMol>>> & getEStateQueries(bool extended) {
      return extended ? GetesExtQueries() : GetesQueries();
    }


    // Function to calculate EState fingerprints
    std::vector<double> calcEStateDescs(const RDKit::ROMol& mol, bool extended) {
        const auto& queries = getEStateQueries(extended);
        //const std::vector<std::pair<std::string, std::string>> esPat = extended ? esPatternsFromOEState : esPatterns;
        size_t nPatts = queries.size();
        std::vector<int> counts(nPatts, 0);
        std::vector<double> sums(nPatts, 0.0);
        std::vector<double> maxValues(nPatts, std::numeric_limits<double>::lowest());
        std::vector<double> minValues(nPatts, std::numeric_limits<double>::max());
        // Parse SMARTS strings into RDKit molecule objects

        // Calculate EState indices for the molecule
        std::vector<double> esIndices = calcEStateIndices(mol);

        size_t i = 0;

        for (const auto& [name, pattern] : queries) {

            // Find all substructure matches
            std::vector<RDKit::MatchVectType> matches;
            RDKit::SubstructMatch(mol, *pattern, matches, true);

            // Update counts, sums, max, and min
            counts[i] = static_cast<int>(matches.size());
            for (const auto& match : matches) {
                int atomIdx = match[0].second; // Atom index from the match
                double value = esIndices[atomIdx];
                sums[i] += value;
                maxValues[i] = std::max(maxValues[i], value);
                minValues[i] = std::min(minValues[i], value);
            }

            // Handle cases where there are no matches
            if (counts[i] == 0) {
                maxValues[i] = 0.0;
                minValues[i] = 0.0;
            }

            ++i;  // Increment the index

        }

        // Concatenate counts, sums, maxValues, and minValues into a single vector
        std::vector<double> results;
        results.reserve(4 * nPatts);
        results.insert(results.end(), counts.begin(), counts.end());    // Counts
        results.insert(results.end(), sums.begin(), sums.end());        // Sums
        results.insert(results.end(), maxValues.begin(), maxValues.end()); // Max values
        results.insert(results.end(), minValues.begin(), minValues.end()); // Min values

        return results;
    }

    std::vector<double> calcHBDHBAtDescs(const RDKit::ROMol& mol,
                                        const std::vector<double>& esIndices) {
        // Indices for HBD, wHBD, HBA, and wHBA patterns
        // const std::vector<int> HBD = {10, 21, 23, 24, 25, 34, 48};
        // const std::vector<int> wHBD = {38, 54};
        // const std::vector<int> HBA = {21, 23, 24, 25, 29, 30, 34, 35, 36, 37, 38, 50, 51, 54, 70};
        // const std::vector<int> wHBA = {18, 10, 11, 12, 14, 15, 16, 17, 18};

        // Indices for HBD, wHBD, HBA, and wHBA patterns
        const std::vector<std::string> HBD = {"dNH", "sNH2", "ssNH2", "aaNH", "sssNH", "ddsN", "aasN"};
        const std::vector<std::string> wHBD = {"sssN", "ssssN"};
        const std::vector<std::string> HBA = {"sOH", "ssO", "dO", "aaO", "ssO(ester)", "dO(acid)", "dO(ester)", "dO(amid)", "dO(nitro)", "dO(sulfo)"};
        const std::vector<std::string> wHBA = {"dO", "sOH", "aaO"};

        // Initialize accumulators for the groups
        int nHBd = 0, nwHBd = 0, nHBa = 0, nwHBa = 0;
        double SHBd = 0.0, SwHBd = 0.0, SHBa = 0.0, SwHBa = 0.0;


        // Function to find index of a SMARTS pattern in `esQueries`
        auto findIndex = [&](const std::string& smarts) -> int {
            auto esQueries = GetesQueries();
            for (size_t i = 0; i < esQueries.size(); ++i) {
                if (esQueries[i].first == smarts) return i;  // Found index
            }
            return -1;  // Not found
        };




        // Helper function to process matches
        auto processMatches = [&](const std::vector<std::string>& patterns, int& count, double& sum) {
            for (const auto& smarts : patterns) {
                int idx = findIndex(smarts);
                if (idx == -1 || !GetesQueries()[idx].second) continue;  // Skip if not found or null pointer

                std::vector<RDKit::MatchVectType> matches;
                RDKit::SubstructMatch(mol, *GetesQueries()[idx].second, matches, true);

                for (const auto& match : matches) {
                    int atomIdx = match[0].second;
                    double value = esIndices[atomIdx];
                    count++;
                    sum += value;
                }
            }
        };

        /*
        for (const auto& smarts : HBD) {
            auto it = esQueries.find(smarts);
            if (it == esQueries.end() || !it->second) continue;

            std::vector<RDKit::MatchVectType> matches;
            RDKit::SubstructMatch(mol, *it->second, matches, true);

            for (const auto& match : matches) {
                int atomIdx = match[0].second;
                double value = esIndices[atomIdx];
                nHBd++;
                SHBd += value;
            }
        }

        for (const auto& smarts : wHBD) {
            auto it = esQueries.find(smarts);
            if (it == esQueries.end() || !it->second) continue;

            std::vector<RDKit::MatchVectType> matches;
            RDKit::SubstructMatch(mol, *it->second, matches, true);

            for (const auto& match : matches) {
                int atomIdx = match[0].second;
                double value = esIndices[atomIdx];
                nwHBd++;
                SwHBd += value;
            }
        }

        for (const auto& smarts : HBA) {
            auto it = esQueries.find(smarts);
            if (it == esQueries.end() || !it->second) continue;

            std::vector<RDKit::MatchVectType> matches;
            RDKit::SubstructMatch(mol, *it->second, matches, true);

            for (const auto& match : matches) {
                int atomIdx = match[0].second;
                double value = esIndices[atomIdx];
                nHBa++;
                SHBa += value;
            }
        }

        for (const auto& smarts : wHBA) {
            auto it = esQueries.find(smarts);
            if (it == esQueries.end() || !it->second) continue;

            std::vector<RDKit::MatchVectType> matches;
            RDKit::SubstructMatch(mol, *it->second, matches, true);

            for (const auto& match : matches) {
                int atomIdx = match[0].second;
                double value = esIndices[atomIdx];
                nwHBa++;
                SwHBa += value;
            }
        }


        // Combine results into a single vector
        std::vector<double> results = {
            static_cast<double>(nHBd), SHBd,
            static_cast<double>(nwHBd), SwHBd,
            static_cast<double>(nHBa), SHBa,
            static_cast<double>(nwHBa), SwHBa
        };



        return results;
        */

        // Process each descriptor group
        processMatches(HBD, nHBd, SHBd);
        processMatches(wHBD, nwHBd, SwHBd);
        processMatches(HBA, nHBa, SHBa);
        processMatches(wHBA, nwHBa, SwHBa);

        // Combine results into a single vector
        return {static_cast<double>(nHBd), SHBd,
                static_cast<double>(nwHBd), SwHBd,
                static_cast<double>(nHBa), SHBa,
                static_cast<double>(nwHBa), SwHBa};

    }


    // Function to calculate HEState fingerprints + need to add the HBD, wHDBm HBA and wHBA patterns
    std::vector<double> calcHEStateDescs(const RDKit::ROMol& mol) {
        auto hsQueries = GetHsQueries();
        size_t nPatts = hsQueries.size();
        std::vector<int> counts(nPatts, 0);
        std::vector<double> sums(nPatts, 0.0);
        std::vector<double> maxValues(nPatts, 0.);
        std::vector<double> minValues(nPatts, 0.);
        // Parse SMARTS strings into RDKit molecule objects

        // Calculate EState indices for the molecule
        std::vector<double> hesIndices = CalcHEStateIndices(mol);

        int i = 0;
        for (const auto& [name, pattern] : hsQueries) {
            if (!pattern) continue; // Skip invalid SMARTS patterns

            // Find all substructure matches
            std::vector<RDKit::MatchVectType> matches;
            RDKit::SubstructMatch(mol, *pattern, matches, true);

            // Update counts, sums, max, and min
            counts[i] = static_cast<int>(matches.size());
            for (const auto& match : matches) {
                int atomIdx = match[0].second; // Atom index from the match
                double value = hesIndices[atomIdx];
                sums[i] += value;
                maxValues[i] = std::max(maxValues[i], value);
                minValues[i] = std::min(minValues[i], value);
            }

            // Handle cases where there are no matches
            if (counts[i] == 0) {
                maxValues[i] = 0.0;
                minValues[i] = 0.0;
            }
            i++;
        }

        // Concatenate counts, sums, maxValues, and minValues into a single vector
        std::vector<double> esIndices = CalcHEStateIndices(mol);
        std::vector<double> HDADescs = calcHBDHBAtDescs(mol,esIndices);

        std::vector<double> results;
        results.reserve(4 * nPatts+HDADescs.size());


        results.insert(results.end(), counts.begin(), counts.end());    // Counts
        results.insert(results.end(), sums.begin(), sums.end());        // Sums
        results.insert(results.end(), maxValues.begin(), maxValues.end()); // Max values
        results.insert(results.end(), minValues.begin(), minValues.end()); // Min values
        results.insert(results.end(), HDADescs.begin(), HDADescs.end()); // HDADescs values

        return results;
    }



// Chi computation Chain dv and d from 3 to 7 orders
    std::vector<double> calcChichain(const RDKit::ROMol& mol) {
        // Output vector: First 10 elements are MPC 2-10 and Total, next 11 are piPC 1-10 and Total
        std::vector<double> results(10, 0.0);

        for (int order = 3; order <= 7; ++order) {
            auto classifiedPaths = extractAndClassifyPaths(mol, order, false);
            double xd = 0.0, xdv = 0.0;
            for (const auto& [path, nodes, type] : classifiedPaths) {
                double cd = 1.0, cdv = 1.0;
                if (type==ChiType::Chain) {
                    for (const auto& node : nodes) {
                        const Atom* at = mol.getAtomWithIdx(node);
                        double d = getSigmaElectrons(*at);  // d
                        cd *= d;
                        double dv = getValenceElectrons(*at);  // dv
                        cdv *= dv;
                    }
                    xd += 1. / std::sqrt(cd);
                    xdv += 1./ std::sqrt(cdv);
                }
            }
            results[order-3]=xd;
            results[5+order-3]=xdv;


        }

        return results;
    }


// Chi computation CLuster ie c for  dv and d from 3 to 6 orders
    std::vector<double> calcChicluster(const RDKit::ROMol& mol) {
        // Output vector: First 10 elements are MPC 2-10 and Total, next 11 are piPC 1-10 and Total
        std::vector<double> results(8, 0.0);

        for (int order = 3; order <= 6; ++order) {
            auto classifiedPaths = extractAndClassifyPaths(mol, order, false);
            double xd = 0.0, xdv = 0.0;
            for (const auto& [path, nodes, type] : classifiedPaths) {
                double cd = 1.0, cdv = 1.0;
                if (type==ChiType::Cluster) {
                    for (const auto& node : nodes) {
                        const Atom* at = mol.getAtomWithIdx(node);
                        double d = getSigmaElectrons(*at);  // d
                        cd *= d;
                        double dv = getValenceElectrons(*at);  // dv
                        cdv *= dv;
                    }
                    xd += 1. / std::sqrt(cd);
                    xdv += 1./ std::sqrt(cdv);
                }
            }
            results[order-3]=xd;
            results[4+order-3]=xdv;


        }

        return results;
    }


// Chi computation Chain dv and d from 3 to 7 orders
    std::vector<double> calcChipathcluster(const RDKit::ROMol& mol) {
        // Output vector: First 10 elements are MPC 2-10 and Total, next 11 are piPC 1-10 and Total
        std::vector<double> results(6, 0.0);

        for (int order = 4; order <= 6; ++order) {
            auto classifiedPaths = extractAndClassifyPaths(mol, order, false);
            double xd = 0.0, xdv = 0.0;
            for (const auto& [path, nodes, type] : classifiedPaths) {
                double cd = 1.0, cdv = 1.0;
                if (type==ChiType::PathCluster) {
                    for (const auto& node : nodes) {
                        const Atom* at = mol.getAtomWithIdx(node);
                        double d = getSigmaElectrons(*at);  // d
                        cd *= d;
                        double dv = getValenceElectrons(*at);  // dv
                        cdv *= dv;
                    }
                    xd += 1. / std::sqrt(cd);
                    xdv += 1./ std::sqrt(cdv);
                }
            }
            results[order-4]=xd;
            results[3+order-4]=xdv;


        }

        return results;
    }


// Chi computation Chain dv and d from 3 to 7 orders
    std::vector<double> calcChipath(const RDKit::ROMol& mol) {
        // Output vector: First 10 elements are MPC 2-10 and Total, next 11 are piPC 1-10 and Total
        std::vector<double> results(32, 0.0);

        double xd = 0.0, xdv = 0.0;

        for (const auto& at : mol.atoms()) {
            double cd = 1.0, cdv = 1.0;

            double d = getSigmaElectrons(*at);  // d
            cd *= d;
            double dv = getValenceElectrons(*at);  // dv
            cdv *= dv;

            xd += 1. / std::sqrt(cd);
            xdv += 1./ std::sqrt(cdv);
        }

        results[0] = xd;
        results[8] = xd/mol.getNumAtoms();
        results[16] = xdv;
        results[24] = xdv/mol.getNumAtoms();

        for (int order = 1; order <= 7; ++order) {
            auto classifiedPaths = extractAndClassifyPaths(mol, order, false);

            double xd = 0.0, xdv = 0.0;

            int nbnodes = 0;
            for (const auto& [bonds, nodes, type] : classifiedPaths) {
                double cd = 1.0, cdv = 1.0;
                if (type==ChiType::Path) {
                    nbnodes+=1;
                    for (const auto& node : nodes) {
                        const Atom* at = mol.getAtomWithIdx(node);
                        double d = getSigmaElectrons(*at);  // d
                        cd *= d;
                        double dv = getValenceElectrons(*at);  // dv
                        cdv *= dv;
                    }
                    xd += 1. / std::sqrt(cd);
                    xdv += 1./ std::sqrt(cdv);
                }
            }
            results[order]=xd;
            results[8+order]=xd/nbnodes;
            results[16+order]=xdv;
            results[24+order]=xdv/nbnodes;


        }

        return results;
    }

std::vector<double> calcAllChiDescriptors(const RDKit::ROMol& mol) {
    std::vector<double> results(56, 0.0); // Full results vector for all descriptors

    // Precompute sigma and valence electrons for all atoms
    std::vector<double> sigmaElectrons(mol.getNumAtoms(), 0.0);
    std::vector<double> valenceElectrons(mol.getNumAtoms(), 0.0);

    double path_0_xd = 0.0, path_0_xdv = 0.0;

    for (const auto& atom : mol.atoms()) {
        auto idx = atom->getIdx();
        sigmaElectrons[idx] = getSigmaElectrons(*atom);
        valenceElectrons[idx] = getValenceElectrons(*atom);
        path_0_xd += 1.0 / std::sqrt(sigmaElectrons[idx]);
        path_0_xdv += 1.0 / std::sqrt(valenceElectrons[idx]);
    }

    // Path-specific accumulators for totals and averages
    double path_xd_total = 0.0, path_xdv_total = 0.0;
    std::vector<double> path_xd(8, 0.0), path_xdv(8, 0.0);
    std::vector<int> path_node_counts(8, 0);

    // Store Path results for order 0
    results[24] = path_0_xd;                     // Total Xp-0d
    results[32] = path_0_xd / mol.getNumAtoms(); // Average Xp-0d
    results[40] = path_0_xdv;                    // Total Xp-0dv
    results[48] = path_0_xdv / mol.getNumAtoms(); // Average Xp-0dv

    for (int order = 1; order <= 7; ++order) {
        auto classifiedPaths = extractAndClassifyPaths(mol, order, false);

        double chain_xd = 0.0, chain_xdv = 0.0;
        double cluster_xd = 0.0, cluster_xdv = 0.0;
        double pathcluster_xd = 0.0, pathcluster_xdv = 0.0;

        for (const auto& [bonds, nodes, type] : classifiedPaths) {
            if (type == ChiType::Chain && order >= 3 && order <= 7) {
                double cd = 1.0, cdv = 1.0;
                for (const auto& node : nodes) {
                    const Atom* at = mol.getAtomWithIdx(node);
                    double d = getSigmaElectrons(*at); // d
                    cd *= d;
                    double dv = getValenceElectrons(*at); // dv
                    cdv *= dv;
                }
                chain_xd += 1.0 / std::sqrt(cd);
                chain_xdv += 1.0 / std::sqrt(cdv);
            } else if (type == ChiType::Cluster && order >= 3 && order <= 6) {
                double cd = 1.0, cdv = 1.0;
                for (const auto& node : nodes) {
                    const Atom* at = mol.getAtomWithIdx(node);
                    double d = getSigmaElectrons(*at); // d
                    cd *= d;
                    double dv = getValenceElectrons(*at); // dv
                    cdv *= dv;
                }
                cluster_xd += 1.0 / std::sqrt(cd);
                cluster_xdv += 1.0 / std::sqrt(cdv);
            } else if (type == ChiType::PathCluster && order >= 4 && order <= 6) {
                double cd = 1.0, cdv = 1.0;
                for (const auto& node : nodes) {
                    const Atom* at = mol.getAtomWithIdx(node);
                    double d = getSigmaElectrons(*at); // d
                    cd *= d;
                    double dv = getValenceElectrons(*at); // dv
                    cdv *= dv;
                }
                pathcluster_xd += 1.0 / std::sqrt(cd);
                pathcluster_xdv += 1.0 / std::sqrt(cdv);
            } else if (type == ChiType::Path) {
                double cd = 1.0, cdv = 1.0;
                for (const auto& node : nodes) {
                    const Atom* at = mol.getAtomWithIdx(node);
                    double d = getSigmaElectrons(*at); // d
                    cd *= d;
                    double dv = getValenceElectrons(*at); // dv
                    cdv *= dv;
                }
                path_node_counts[order] += 1;
                path_xd[order] += 1.0 / std::sqrt(cd);
                path_xdv[order] += 1.0 / std::sqrt(cdv);
            }
        }

        // Update total Path values
        if (order <= 7) {
            path_xd_total += path_xd[order];
            path_xdv_total += path_xdv[order];
        }

        // Store Chain results (order 3-7)
        if (order >= 3 && order <= 7) {
            results[order - 3] = chain_xd;
            results[5 + (order - 3)] = chain_xdv;
        }

        // Store Cluster results (order 3-6)
        if (order >= 3 && order <= 6) {
            results[10 + (order - 3)] = cluster_xd;
            results[14 + (order - 3)] = cluster_xdv;
        }

        // Store PathCluster results (order 4-6)
        if (order >= 4 && order <= 6) {
            results[18 + (order - 4)] = pathcluster_xd;
            results[21 + (order - 4)] = pathcluster_xdv;
        }
    }

    // Finalize Path results
    for (int order = 1; order <= 7; ++order) {
        results[24 + order] = path_xd[order];                        // Total path xd
        results[32 + order] = path_xd[order] / path_node_counts[order]; // Average path xd
        results[40 + order] = path_xdv[order];                       // Total path xdv
        results[48 + order] = path_xdv[order] / path_node_counts[order]; // Average path xdv
    }



    return results;
}


///////
    // Define the graph as an adjacency list
    using Graph = std::unordered_map<int, std::vector<std::pair<int, double>>>;


    // Build the molecular graph
    Graph buildGraph(const RDKit::ROMol& mol) {
        Graph graph;
        for (const auto& bond : mol.bonds()) {
            int start = bond->getBeginAtomIdx();
            int end = bond->getEndAtomIdx();
            double weight = static_cast<double>(mol.getAtomWithIdx(start)->getDegree() * mol.getAtomWithIdx(end)->getDegree());

            graph[start].emplace_back(end, weight);
            graph[end].emplace_back(start, weight);
        }
        return graph;
    }

    // Recursive DFS for atomic ID computation
    double computeAtomicId(const Graph& graph, int atomIdx, double epsilon, double currentWeight,
                        std::unordered_set<int>& visited, double limit) {
        double id = 0.0;

        visited.insert(atomIdx);

        for (const auto& [nextAtom, edgeWeight] : graph.at(atomIdx)) {
            if (visited.count(nextAtom)) continue;

            double combinedWeight = currentWeight * edgeWeight;

            id += 1.0 / std::sqrt(combinedWeight);  // Contribution to ID

            if (combinedWeight < limit) {
                id += computeAtomicId(graph, nextAtom, epsilon, combinedWeight, visited, limit);
            }
        }

        visited.erase(atomIdx);  // Backtrack
        return id;
    }


    // Compute atomic IDs for all atoms in the molecule
    std::vector<double> computeAtomicIds(const RDKit::ROMol& mol, double epsilon) {
        int natoms = mol.getNumAtoms();
        std::vector<double> atomicIds(natoms, 0.0);
        if (natoms>1) {
            Graph graph = buildGraph(mol);
            double limit = 1.0 / (epsilon * epsilon);

            for (int atomIdx = 0; atomIdx < natoms; ++atomIdx) {
                std::unordered_set<int> visited;
                double id = computeAtomicId(graph, atomIdx, epsilon, 1.0, visited, limit);
                atomicIds[atomIdx] = 1.0 + id / 2.0;  // Normalize
            }
        }
        return atomicIds;
    }

    std::vector<double> calcMolecularId(const RDKit::ROMol& mol) {
        double epsilon = 1e-10;
        std::vector<double> atomicIds = computeAtomicIds(mol, epsilon);
        size_t numAtoms = mol.getNumAtoms();

        // 12 results:
        // MID (all), AMID (all averaged), MID_h, AMID_h, MID_C, AMID_C, MID_N, AMID_N, MID_O, AMID_O, MID_X, AMID_X
        std::vector<double> results(12, 0.0);

        for (size_t i = 0; i < atomicIds.size(); ++i) {
            const auto* atom = mol.getAtomWithIdx(i);
            int atomicNum = atom->getAtomicNum();

            // All atoms
            results[0] += atomicIds[i]; // MID
            if (atomicNum > 1 && atomicNum != 6) { // Heavy atoms
                results[2] += atomicIds[i]; // MID_h
            }

            // Specific element sums
            if (atomicNum == 6) results[4] += atomicIds[i]; // MID_C
            if (atomicNum == 7) results[6] += atomicIds[i]; // MID_N
            if (atomicNum == 8) results[8] += atomicIds[i]; // MID_O
            if (atomicNum == 9 || atomicNum == 17 || atomicNum == 35 || atomicNum == 53 || atomicNum == 85) {
                // Halogens: F, Cl, Br, I, At
                results[10] += atomicIds[i]; // MID_X
            }
        }

        // Compute averages of each type
        results[1] = results[0] / numAtoms; // AMID
        results[3] = results[2] / numAtoms; // AMID_h
        results[5] = results[4] / numAtoms; // AMID_C
        results[7] = results[6] / numAtoms; // AMID_N
        results[9] = results[8] / numAtoms; // AMID_O
        results[11] = results[10] / numAtoms; // AMID_X

        return results;
    }

    // Function to build the Burden matrix
    Eigen::MatrixXd buildBurdenMatrix(const RDKit::ROMol& mol) {
        size_t numAtoms = mol.getNumAtoms();
        Eigen::MatrixXd burdenMatrix = Eigen::MatrixXd::Constant(numAtoms, numAtoms, 0.001);

        // Iterate through bonds in the molecule
        for (const auto& bond : mol.bonds()) {
            const auto* a = bond->getBeginAtom();
            const auto* b = bond->getEndAtom();
            int i = a->getIdx();
            int j = b->getIdx();

            double w = bond->getBondTypeAsDouble() / 10.0;
            if (a->getDegree() == 1 || b->getDegree() == 1) {
                w += 0.01;
            }
            burdenMatrix(i, j) = w;
            burdenMatrix(j, i) = w;
        }

        return burdenMatrix;
    }

    // Function to compute the highest and lowest eigenvalues of the Burden matrix
    Eigen::VectorXd calcEigenValues(const Eigen::MatrixXd& burdenMatrix) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(burdenMatrix);

        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Eigenvalue computation failed");
        }

        Eigen::VectorXd eigenvalues = solver.eigenvalues().real();

        return eigenvalues;
    }



// v2.0: Control function to check if Gasteiger parameters exist for all atoms
// BEFORE calling computeGasteigerCharges. This prevents crashes on unusual elements
// like Hg (mercury) or certain phosphorus environments.
//
// Returns:
//   - true: All atoms have parameters -> safe to call computeGasteigerCharges
//   - false: At least one atom lacks parameters -> return NaN instead
bool checkGasteigerParameters(const RDKit::ROMol& mol) {
    PeriodicTable *table = PeriodicTable::getTable();
    const GasteigerParams *params = GasteigerParams::getParams();
    
    // Check each atom: if ONE atom doesn't have parameters, the ENTIRE calculation is invalid
    for (auto &atom : mol.atoms()) {
        std::string elem = table->getElementSymbol(atom->getAtomicNum());
        std::string mode;
        
        // Reproduce EXACTLY the logic from computeGasteigerCharges to determine mode
        switch (atom->getHybridization()) {
            case Atom::SP3:
                mode = "sp3";
                break;
            case Atom::SP2:
                mode = "sp2";
                break;
            case Atom::SP:
                mode = "sp";
                break;
            default:
                // For atoms with unspecified hybridization, determine mode based on element
                if (atom->getAtomicNum() == 1) {
                    mode = "*"; // Hydrogen always uses "*" mode
                } else if (atom->getAtomicNum() == 16) {
                    // Sulfur: check oxygen neighbors
                    ROMol::ADJ_ITER nbrIdx, endIdx;
                    boost::tie(nbrIdx, endIdx) = mol.getAtomNeighbors(atom);
                    int no = 0;
                    while (nbrIdx != endIdx) {
                        if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() == 8) {
                            no++;
                        }
                        nbrIdx++;
                    }
                    if (no == 2) {
                        mode = "so2";
                    } else if (no == 1) {
                        mode = "so";
                    } else {
                        mode = "sp3";
                    }
                } else {
                    mode = "sp3"; // Default
                }
        }
        
        // Check if parameters exist for this element + mode
        try {
            params->getParams(elem, mode, true);
        } catch (...) {
            // Parameters don't exist (e.g., Hg, certain P environments)
            return false;
        }
    }
    
    return true;
}

std::vector<double> calcRNCG_RPCG(const RDKit::ROMol& mol){
        // subclass of CPSA using only 2D descriptors available in Mordred v1
        // v2.0: Check Gasteiger parameters first (on original mol, before adding H)
        if (!checkGasteigerParameters(mol)) {
            return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
        }
        
         std::unique_ptr<RDKit::ROMol> hmol(RDKit::MolOps::addHs(mol));

        double maxpos = 0;
        double maxneg = 0;
        double totalneg =0;
        double totalpos =0;

        // v2.0: Try-catch as safety net
        try {
            computeGasteigerCharges(*hmol, 12, true);
        } catch (...) {
            return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
        }

        for (auto &atom : hmol->atoms()) {
                double ch = atom->getProp<double>(common_properties::_GasteigerCharge);

                if (atom->hasProp(common_properties::_GasteigerHCharge)) {
                    ch += atom->getProp<double>(common_properties::_GasteigerHCharge);
                }

                if (ch < 0) {
                    totalneg += -ch;
                    if (-ch > maxneg) {
                        maxneg = -ch;
                    }


                }
                else if ( ch > 0) {
                    totalpos += ch;
                    if (ch > maxpos) {
                        maxpos = ch;
                    }

                }
        }
        if (totalneg == 0 || totalpos == 0){
            return {0.,0.};
        }
        return {maxneg/totalneg, maxpos/totalpos};



}

    // Function to compute BCUT descriptors for multiple properties
    std::vector<double> calcBCUTsEigen(const RDKit::ROMol& mol) {
        // v2.0: Check Gasteiger parameters first
        bool gasteiger_ok = checkGasteigerParameters(mol);
        
        // Atomic properties to compute
        auto* tbl = RDKit::PeriodicTable::getTable();
        std::map<int, double> vdwMap = VdWAtomicMap();
        std::map<int, double> sandersonENMap = SandersonENAtomicMap();
        std::map<int, double> paulingENMap = PaulingENAtomicMap();
        std::map<int, double> allredENMap = Allred_rocow_ENAtomicMap();
        std::map<int, double> polarizabilityMap = Polarizability94AtomicMap();
        std::map<int, double> ionizationMap = ionizationEnergyAtomicMap();

        size_t numAtoms = mol.getNumAtoms();
        std::vector<double> gasteigerCharges(numAtoms, 0.0);
        
        // v2.0: Only compute Gasteiger if parameters exist
        if (gasteiger_ok) {
            try {
                computeGasteigerCharges(mol, 12, true);
                for (size_t i = 0; i < numAtoms; ++i) {
                    const auto* atom = mol.getAtomWithIdx(i);
                    double ch = atom->getProp<double>(common_properties::_GasteigerCharge);
                    if (atom->hasProp(common_properties::_GasteigerHCharge)) {
                        ch += atom->getProp<double>(common_properties::_GasteigerHCharge);
                    }
                    gasteigerCharges[i] = ch;
                }
            } catch (...) {
                // Gasteiger failed - charges stay at 0.0
                gasteiger_ok = false;
            }
        }
        // If gasteiger_ok is false, gasteigerCharges remains all zeros

        // Vector to store BCUT results
        std::vector<double> results;
        results.reserve(24);

        // Build the Burden matrix
        Eigen::MatrixXd burdenMatrix = buildBurdenMatrix(mol);

        // List of atomic property vectors
        std::vector<std::vector<double>> atomicProperties(12, std::vector<double>(numAtoms, 0.0));


        // Populate atomic property vectors
        for (size_t i = 0; i < numAtoms; ++i) {
            const auto* atom = mol.getAtomWithIdx(i);
            int atomNumber = atom->getAtomicNum();
            atomicProperties[0][i] = gasteigerCharges[i];                                    // Gasteiger charge (c)
            atomicProperties[1][i] = getValenceElectrons(*atom);                             // Valence electrons (dv)
            atomicProperties[2][i] = getSigmaElectrons(*atom);                              // Sigma electrons (d)
            atomicProperties[3][i] = getIntrinsicState(*atom);                              // Intrinsic state (s)
            atomicProperties[4][i] = static_cast<double>(atomNumber);                         // Atomic number (Z)
            atomicProperties[5][i] = tbl->getAtomicWeight(atomNumber);                       // Atomic weight (m)
            atomicProperties[6][i] = vdw_volume(vdwMap[atomNumber]);                          // Van der Waals volume need vdw_volume to go from (r) to (v)
            atomicProperties[7][i] = sandersonENMap[atomNumber];                             // Sanderson electronegativity (se)
            atomicProperties[8][i] = paulingENMap[atomNumber];                               // Pauling electronegativity (pe)
            atomicProperties[9][i] = allredENMap[atomNumber];                                // Allred-Rocow electronegativity (are)
            atomicProperties[10][i] = polarizabilityMap[atomNumber];                          // Polarizability (p)
            atomicProperties[11][i] = ionizationMap[atomNumber];                              // Ionization energy (i)
        }
        // Compute eigenvalues for each atomic property
        for (const auto& prop : atomicProperties) {
            Eigen::MatrixXd adjustedMatrix = burdenMatrix;
            for (size_t i = 0; i < numAtoms; ++i) {
                adjustedMatrix(i, i) = prop[i];
            }

            auto eigenValues = calcEigenValues(adjustedMatrix);

            results.push_back(eigenValues.maxCoeff());
            results.push_back(eigenValues.minCoeff());
        }

        return results;
    }

// Function to build the Burden matrix
std::vector<std::vector<double>> buildBurdenMatrixL(const RDKit::ROMol& mol) {
    size_t numAtoms = mol.getNumAtoms();
    std::vector<std::vector<double>> burdenMatrix(numAtoms, std::vector<double>(numAtoms, 0.001));

    // Iterate through bonds in the molecule
    for (const auto& bond : mol.bonds()) {
        const auto* a = bond->getBeginAtom();
        const auto* b = bond->getEndAtom();
        int i = a->getIdx();
        int j = b->getIdx();

        double w = bond->getBondTypeAsDouble() / 10.0;
        if (a->getDegree() == 1 || b->getDegree() == 1) {
            w += 0.01;
        }
        burdenMatrix[i][j] = w;
        burdenMatrix[j][i] = w;
    }

    return burdenMatrix;
}

// Function to compute eigenvalues using LAPACK ...
std::vector<double> calcEigenValuesLAPACK(const std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    std::vector<double> flatMatrix(n * n);
    std::vector<double> eigenvalues(n);

    // Convert 2D matrix to a flat array (column-major order)
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            flatMatrix[j * n + i] = matrix[i][j];

    // LAPACKE_dsyev computes eigenvalues and eigenvectors
    int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'N', 'U', n, flatMatrix.data(), n, eigenvalues.data());
    if (info != 0) {
        throw std::runtime_error("Error in LAPACKE_dsyev: " + std::to_string(info));
    }

    return eigenvalues;
}

// Function to compute BCUT descriptors for multiple properties
std::vector<double> calcBCUTs(const RDKit::ROMol& mol) {
    // v2.0: Check Gasteiger parameters first
    bool gasteiger_ok = checkGasteigerParameters(mol);
    
    // Atomic properties to compute
    auto* tbl = RDKit::PeriodicTable::getTable();
    std::map<int, double> vdwMap = VdWAtomicMap();
    std::map<int, double> sandersonENMap = SandersonENAtomicMap();
    std::map<int, double> paulingENMap = PaulingENAtomicMap();
    std::map<int, double> allredENMap = Allred_rocow_ENAtomicMap();
    std::map<int, double> polarizabilityMap = Polarizability94AtomicMap();
    std::map<int, double> ionizationMap = ionizationEnergyAtomicMap();

    size_t numAtoms = mol.getNumAtoms();
    std::vector<double> gasteigerCharges(numAtoms, 0.0);
    
    // v2.0: Only compute Gasteiger if parameters exist
    if (gasteiger_ok) {
        try {
            computeGasteigerCharges(mol, 12, true);
            for (size_t i = 0; i < numAtoms; ++i) {
                const auto* atom = mol.getAtomWithIdx(i);
                double ch = atom->getProp<double>(common_properties::_GasteigerCharge);
                if (atom->hasProp(common_properties::_GasteigerHCharge)) {
                    ch += atom->getProp<double>(common_properties::_GasteigerHCharge);
                }
                gasteigerCharges[i] = ch;
            }
        } catch (...) {
            // Gasteiger failed - charges stay at 0.0
        }
    }
    // If gasteiger_ok is false, gasteigerCharges remains all zeros

    // Vector to store BCUT results
    std::vector<double> results;
    results.reserve(24);

    // Build the Burden matrix
    auto burdenMatrix = buildBurdenMatrixL(mol);

    // List of atomic property vectors
    std::vector<std::vector<double>> atomicProperties(12, std::vector<double>(numAtoms, 0.0));

    // Populate atomic property vectors
    for (size_t i = 0; i < numAtoms; ++i) {
        const auto* atom = mol.getAtomWithIdx(i);
        int atomNumber = atom->getAtomicNum();
        atomicProperties[0][i] = gasteigerCharges[i];                                    // Gasteiger charge (c)
        atomicProperties[1][i] = getValenceElectrons(*atom);                             // Valence electrons (dv)
        atomicProperties[2][i] = getSigmaElectrons(*atom);                              // Sigma electrons (d)
        atomicProperties[3][i] = getIntrinsicState(*atom);                              // Intrinsic state (s)
        atomicProperties[4][i] = static_cast<double>(atomNumber);                         // Atomic number (Z)
        atomicProperties[5][i] = tbl->getAtomicWeight(atomNumber);                       // Atomic weight (m)
        atomicProperties[6][i] = vdw_volume(vdwMap[atomNumber]);                          // Van der Waals volume need vdw_volume to go from (r) to (v)
        atomicProperties[7][i] = sandersonENMap[atomNumber];                             // Sanderson electronegativity (se)
        atomicProperties[8][i] = paulingENMap[atomNumber];                               // Pauling electronegativity (pe)
        atomicProperties[9][i] = allredENMap[atomNumber];                                // Allred-Rocow electronegativity (are)
        atomicProperties[10][i] = polarizabilityMap[atomNumber];                          // Polarizability (p)
        atomicProperties[11][i] = ionizationMap[atomNumber];                              // Ionization energy (i)
    }

    // Compute eigenvalues for each atomic property
    for (const auto& prop : atomicProperties) {
        // Adjust diagonal with atomic property
        auto adjustedMatrix = burdenMatrix;
        for (size_t i = 0; i < numAtoms; ++i) {
            adjustedMatrix[i][i] = prop[i];
        }

        auto eigenValues = calcEigenValuesLAPACK(adjustedMatrix);

        // v3 fix: descriptor names are <prop>-1l (low) then <prop>-1h (high),
        // so emit min then max. Previously emitted max then min, swapping the
        // -1l/-1h label of every BCUT descriptor.
        results.push_back(*std::min_element(eigenValues.begin(), eigenValues.end())); // <prop>-1l : min eigenvalue
        results.push_back(*std::max_element(eigenValues.begin(), eigenValues.end())); // <prop>-1h : max eigenvalue
    }

    return results;
}

    // Constitutional
    // code vs name
    //Z    a.GetAtomicNum()
    //m    mass[a.GetAtomicNum()]
    //v    vdw_volume[a.GetAtomicNum()]
    //se   sanderson[a.GetAtomicNum()]
    //pe   pauling[a.GetAtomicNum()]
    //are  allred_rocow[a.GetAtomicNum()]
    //p    polarizability94[a.GetAtomicNum()] (last as default!!!)
    //i    ionization_potentials[a.GetAtomicNum()]
    //c    gasteiger charge


// autocorrelation


    // Function to compute the ATS descriptors
    std::vector<double> calcAutoCorrelationEigen(const RDKit::ROMol &mol) {
        std::unique_ptr<RDKit::ROMol> hmol(RDKit::MolOps::addHs(mol));

        double *dist = RDKit::MolOps::getDistanceMat(*hmol, false); // Topological matrix
        unsigned int numAtoms = hmol->getNumAtoms();

        // Lookup tables
        std::map<int, double> vdwmap = VdWAtomicMap();
        std::map<int, double> semap = SandersonENAtomicMap();
        std::map<int, double> pemap = PaulingENAtomicMap();
        std::map<int, double> aremap = Allred_rocow_ENAtomicMap();
        std::map<int, double> pmap = Polarizability94AtomicMap();
        std::map<int, double> imap = ionizationEnergyAtomicMap();
        const auto* tbl = RDKit::PeriodicTable::getTable();

        // Eigen vector for atomic properties
        Eigen::MatrixXd propertyMatrix(12, numAtoms);
        propertyMatrix.setZero();

        // Compute Gasteiger charges
        computeGasteigerCharges(*hmol, 12, true);
        std::vector<double> gasteigerCharges(numAtoms, 0.0);

        for (unsigned int i = 0; i < numAtoms; ++i) {
            const auto* atom = hmol->getAtomWithIdx(i);
            double charge = atom->getProp<double>(RDKit::common_properties::_GasteigerCharge);
            if (atom->hasProp(RDKit::common_properties::_GasteigerHCharge)) {
                charge += atom->getProp<double>(RDKit::common_properties::_GasteigerHCharge);
            }
            gasteigerCharges[i] = charge;

            int atomNumber = atom->getAtomicNum();
            propertyMatrix(0, i) = getValenceElectrons(*atom);
            propertyMatrix(1, i) = getSigmaElectrons(*atom);
            propertyMatrix(2, i) = getIntrinsicState(*atom);
            propertyMatrix(3, i) = static_cast<double>(atomNumber);
            propertyMatrix(4, i) = tbl->getAtomicWeight(atomNumber);
            propertyMatrix(5, i) = vdw_volume(vdwmap[atomNumber]);
            propertyMatrix(6, i) = semap[atomNumber];
            propertyMatrix(7, i) = pemap[atomNumber];
            propertyMatrix(8, i) = aremap[atomNumber];
            propertyMatrix(9, i) = pmap[atomNumber];
            propertyMatrix(10, i) = imap[atomNumber];
            propertyMatrix(11, i) = gasteigerCharges[i]; // need to change order ...

        }

        // Initialize the topological distance matrix as an Eigen matrix
        Eigen::MatrixXd distanceMatrix(numAtoms, numAtoms);
        for (unsigned int i = 0; i < numAtoms; ++i) {
            for (unsigned int j = 0; j < numAtoms; ++j) {
                distanceMatrix(i, j) = dist[i * numAtoms + j];
            }
        }

        // Compute the ATS descriptors
        const int maxDistance = 9;
        std::vector<double> ATSDescriptors(11 * maxDistance, 0.0); // not charges descriptors
        std::vector<double> AATSDescriptors(11 * maxDistance, 0.0); // not charges descriptors
        std::vector<double> ATSCDescriptors(12 * maxDistance, 0.0); // we have now Charges too
        std::vector<double> AATSCDescriptors(12 * maxDistance, 0.0); // we have now Charges too
        std::vector<double> MATSDescriptors(12 * (maxDistance-1), 0.0); // we have now Charges too but not zeros order
        std::vector<double> GATSDescriptors(12 * (maxDistance-1), 0.0); // we have now Charges too but not zeros order

        // ATSC for centered ie avec - avec.mean()

        for (int k = 0; k < maxDistance; ++k) {


            Eigen::MatrixXd binaryMatrix = (distanceMatrix.array() == (k)).cast<double>();
            double gsum = binaryMatrix.array().sum();
            for (int t = 0; t < 12; ++t) {
                Eigen::VectorXd weights = propertyMatrix.row(t).transpose();
                Eigen::VectorXd weights_centered = weights.array() - weights.mean();
                Eigen::MatrixXd weights_col = weights.replicate(1, weights.size());
                Eigen::MatrixXd weights_row = weights.transpose().replicate(weights.size(), 1);
                Eigen::MatrixXd diffSquared = (weights_row.array() - weights_col.array()).square();

                double gatsdenominator = weights_centered.array().square().sum() / (weights_centered.size() - 1);
                if (k > 0) {
                    if (t<11) {
                        ATSDescriptors[t * maxDistance + k] = 0.5 * (weights.transpose() * binaryMatrix * weights).value();
                        AATSDescriptors[t * maxDistance + k] = ATSDescriptors[t * maxDistance + k] / (0.5 * gsum);
                        ATSCDescriptors[(t+1) * maxDistance + k] =  0.5 * (weights_centered.transpose() * binaryMatrix * weights_centered).value();
                        AATSCDescriptors[(t+1) * maxDistance + k] =  ATSCDescriptors[(t+1) * maxDistance + k] / (0.5 * gsum) ;
                        MATSDescriptors [(t+1) * (maxDistance-1) + k-1] = weights.size() * AATSCDescriptors[(t+1) * maxDistance + k] / weights_centered.array().square().sum();

                        // Compute n
                        GATSDescriptors [(t+1) * (maxDistance-1) + k-1] =  ((binaryMatrix.array() * diffSquared.array()).sum() / (2.0 * gsum)) / gatsdenominator;




                    } else {
                        ATSCDescriptors[k] =  0.5 * (weights_centered.transpose() * binaryMatrix * weights_centered).value();
                        AATSCDescriptors[k] = ATSCDescriptors[k] / (0.5*gsum) ;
                        MATSDescriptors [k-1] =  weights.size() * AATSCDescriptors[k]  / weights_centered.array().square().sum();
                        GATSDescriptors [k-1] =  ((binaryMatrix.array() * diffSquared.array()).sum() / (2.0 * gsum)) / gatsdenominator;

                    }
                } else {
                    if (t<11) {
                        ATSDescriptors[t * maxDistance + k] = weights.array().square().sum();
                        AATSDescriptors[t * maxDistance + k] = ATSDescriptors[t * maxDistance + k] / gsum;
                        ATSCDescriptors[(t+1) * maxDistance + k] =  weights_centered.array().square().sum();
                        AATSCDescriptors[(t+1) * maxDistance + k] = ATSCDescriptors[(t+1) * maxDistance + k]  / gsum;
                    }
                    else {
                        ATSCDescriptors[k] =   weights_centered.array().square().sum();
                        AATSCDescriptors[k] =   ATSCDescriptors[k] / gsum;

                    }


                }
            }
        }

        std::vector<double> autocorrDescriptors;
        autocorrDescriptors.reserve(606); // Pre-allocate memory for 606 elements

        // Append ATSDescriptors
        for (const auto &val : ATSDescriptors) {
            autocorrDescriptors.push_back(val);
        }

        // Append AATSDescriptors
        for (const auto &val : AATSDescriptors) {
            autocorrDescriptors.push_back(val);
        }

        // Append ATSCDescriptors
        for (const auto &val : ATSCDescriptors) {
            autocorrDescriptors.push_back(val);
        }

        // Append AATSCDescriptors
        for (const auto &val : AATSCDescriptors) {
            autocorrDescriptors.push_back(val);
        }

        // Append MATSDescriptors
        for (const auto &val : MATSDescriptors) {
            autocorrDescriptors.push_back(val);
        }

        // Append GATSDescriptors
        for (const auto &val : GATSDescriptors) {
            autocorrDescriptors.push_back(val);
        }

        return autocorrDescriptors;
    }



    // Function to compute the ATS descriptors without Eigen
    std::vector<double> calcAutoCorrelation(const RDKit::ROMol &mol) {
        // v2.0: Check Gasteiger parameters first (on original mol, before adding H)
        bool gasteiger_ok = checkGasteigerParameters(mol);
        
        std::unique_ptr<RDKit::ROMol> hmol(RDKit::MolOps::addHs(mol));
        double *dist = RDKit::MolOps::getDistanceMat(*hmol, false);  // Topological matrix
        const unsigned int numAtoms = hmol->getNumAtoms();
        const unsigned int numProperties = 12;
        // Lookup tables
        std::map<int, double> vdwmap = VdWAtomicMap();
        std::map<int, double> semap = SandersonENAtomicMap();
        std::map<int, double> pemap = PaulingENAtomicMap();
        std::map<int, double> aremap = Allred_rocow_ENAtomicMap();
        std::map<int, double> pmap = Polarizability94AtomicMap();
        std::map<int, double> imap = ionizationEnergyAtomicMap();
        const auto *tbl = RDKit::PeriodicTable::getTable();
        // Property matrix as a vector of vectors
        std::vector<std::vector<double>> propertyMatrix(numProperties, std::vector<double>(numAtoms, 0.0)); // correct

        // v2.0: Compute Gasteiger charges only if parameters exist
        std::vector<double> gasteigerCharges(numAtoms, 0.0);
        if (gasteiger_ok) {
            try {
                computeGasteigerCharges(*hmol, 12, true);
                for (unsigned int i = 0; i < numAtoms; ++i) {
                    const auto *atom = hmol->getAtomWithIdx(i);
                    double charge = atom->getProp<double>(RDKit::common_properties::_GasteigerCharge);
                    if (atom->hasProp(RDKit::common_properties::_GasteigerHCharge)) {
                        charge += atom->getProp<double>(RDKit::common_properties::_GasteigerHCharge);
                    }
                    gasteigerCharges[i] = charge;
                }
            } catch (...) {
                // Gasteiger failed - charges stay at 0.0
            }
        }

        for (unsigned int i = 0; i < numAtoms; ++i) {
            const auto *atom = hmol->getAtomWithIdx(i);
            int atomNumber = atom->getAtomicNum();
            propertyMatrix[0][i] = gasteigerCharges[i]; // not for ATS & AATS output
            propertyMatrix[1][i] = getValenceElectrons(*atom);
            propertyMatrix[2][i] = getSigmaElectrons(*atom);
            propertyMatrix[3][i] = getIntrinsicState(*atom);
            propertyMatrix[4][i] = static_cast<double>(atomNumber);
            propertyMatrix[5][i] = tbl->getAtomicWeight(atomNumber);
            propertyMatrix[6][i] = vdw_volume(vdwmap[atomNumber]);
            propertyMatrix[7][i] = semap[atomNumber];
            propertyMatrix[8][i] = pemap[atomNumber];
            propertyMatrix[9][i] = aremap[atomNumber];
            propertyMatrix[10][i] = pmap[atomNumber];
            propertyMatrix[11][i] = imap[atomNumber];
        }

        // Initialize the topological symetric distance matrix without diagonal
        std::vector<std::vector<double>> distanceMatrix(numAtoms, std::vector<double>(numAtoms, 0.0));
        for (unsigned int i = 0; i < numAtoms; ++i) {
            for (unsigned int j = i+1; j < numAtoms; ++j) {
                distanceMatrix[i][j] = dist[i * numAtoms + j];
            }
        }

        // Compute the ATS descriptors
        const int maxDistance = 8;

        std::vector<std::vector<double>> ATS_(maxDistance + 1, std::vector<double>(numProperties-1, 0.0));
        std::vector<std::vector<double>> AATS_(maxDistance + 1, std::vector<double>(numProperties-1, 0.0));
        std::vector<std::vector<double>> ATSC(maxDistance + 1, std::vector<double>(numProperties, 0.0));
        std::vector<std::vector<double>> AATSC(maxDistance + 1, std::vector<double>(numProperties, 0.0));
        std::vector<std::vector<double>> MATS(maxDistance + 1, std::vector<double>(numProperties, 0.0)); // we skip the first one at the end!
        std::vector<std::vector<double>> GATS(maxDistance + 1, std::vector<double>(numProperties, 0.0)); // we skip the first one at the end!

        // Centered property values
        std::vector<std::vector<double>> centeredProperties(numProperties, std::vector<double>(numAtoms, 0.0));
        for (unsigned int t = 0; t < numProperties; ++t) {
            double sum = 0.0;
            for (unsigned int i = 0; i < numAtoms; ++i) {
                sum += propertyMatrix[t][i];
            }
            double mean = sum / numAtoms;
            for (unsigned int i = 0; i < numAtoms; ++i) {
                centeredProperties[t][i] = propertyMatrix[t][i] - mean;
            }
        }

        // Lag 0: self-correlations
        for (unsigned int t = 0; t < numProperties; ++t) {
            for (unsigned int i = 0; i < numAtoms; ++i) {
                if (t >0) {
                    // we skip Charges
                    ATS_[0][t-1] += propertyMatrix[t][i] * propertyMatrix[t][i];
                }
                ATSC[0][t] += centeredProperties[t][i] * centeredProperties[t][i];
            }
            if (t>0) {
                // we skip Charges
                AATS_[0][t-1] = ATS_[0][t-1] / numAtoms;
            }
            AATSC[0][t] = ATSC[0][t] / numAtoms;
            MATS[0][t] = AATSC[0][t];
            GATS[0][t] = ATSC[0][t] / (numAtoms - 1);
        }

        // Lags 1 to maxLag: pairwise correlations
        for (int k = 1; k <= maxDistance; ++k) {
            int maxkVertexPairs = 0;
            for (unsigned int i = 0; i < numAtoms; ++i) {
                for (unsigned int j = i + 1; j < numAtoms; ++j) {
                    if (distanceMatrix[i][j] == k) {
                        ++maxkVertexPairs;
                        for (unsigned int t = 0; t < numProperties; ++t) {
                            double diff = propertyMatrix[t][i] - propertyMatrix[t][j];
                            if (t>0) {
                                ATS_[k][t-1] += propertyMatrix[t][i] * propertyMatrix[t][j];
                            }
                            ATSC[k][t] += centeredProperties[t][i] * centeredProperties[t][j];
                            GATS[k][t] += diff * diff;
                        }
                    }
                }
            }

            if (maxkVertexPairs > 0) {
                for (unsigned int t = 0; t < numProperties; ++t) {
                    if (t>0) {
                        AATS_[k][t-1] = ATS_[k][t-1] / maxkVertexPairs;
                    }
                    AATSC[k][t] = ATSC[k][t] / maxkVertexPairs;
                    if (MATS[0][t] > 0.0) {
                        MATS[k][t] = AATSC[k][t] / MATS[0][t];
                    }
                    if (GATS[0][t] > 0.0) {
                        GATS[k][t] /= (2 * maxkVertexPairs * GATS[0][t]);
                    }
                }
            }
        }
        std::vector<double> descriptors;


        // Flatten the descriptors into a single vector
        for (const auto& vec : {ATS_, AATS_}) {
            for (unsigned int t = 0; t < numProperties-1; ++t) {
                for (unsigned int k = 0; k <= maxDistance; ++k) {
                    descriptors.push_back(vec[k][t]);
                }
            }
        }


        // Flatten the descriptors into a single vector
        for (const auto& vec : { ATSC, AATSC}) {
            for (unsigned int t = 0; t < numProperties; ++t) {
                for (unsigned int k = 0; k <= maxDistance; ++k) {
                    descriptors.push_back(vec[k][t]);
                }
            }
        }

        // Flatten the descriptors into a single vector
        for (const auto& vec : { MATS, GATS}) {
            for (unsigned int t = 0; t < numProperties; ++t) {
                for (unsigned int k = 1; k <= maxDistance; ++k) {
                    descriptors.push_back(vec[k][t]);
                }
            }
        }


        return descriptors;
    }



    std::unordered_set<int> findLinkersWithBFS(const ROMol& mol, const std::unordered_set<int>& ringAtoms) {
        const RingInfo* ringInfo = mol.getRingInfo();
        std::unordered_set<int> nonRingAtoms;
        std::unordered_set<int> linkers;

        // Collect all non-ring atoms
        int numAtoms = mol.getNumAtoms();
        for (int i = 0; i < numAtoms; ++i) {
            if (ringAtoms.find(i) == ringAtoms.end()) {
                nonRingAtoms.insert(i);
            }
        }

        // Perform BFS for each pair of fused rings
        std::vector<std::vector<int>> atomRings = ringInfo->atomRings();
        for (unsigned int i = 0; i < atomRings.size(); ++i) {
            for (unsigned int j = i + 1; j < atomRings.size(); ++j) {
                // Identify start atoms (non-ring neighbors of the first ring)
                std::queue<std::pair<int, std::vector<int>>> q; // Queue of (current atom, path)
                std::unordered_set<int> visited;

                for (int ringAtom : atomRings[i]) {
                    const Atom* atom = mol.getAtomWithIdx(ringAtom);
                    for (const Atom* neighbor : mol.atomNeighbors(atom)) {
                        int neighborIdx = neighbor->getIdx();
                        if (nonRingAtoms.find(neighborIdx) != nonRingAtoms.end()) {
                            q.push({neighborIdx, {neighborIdx}});
                            visited.insert(neighborIdx);
                        }
                    }
                }

                // BFS traversal
                while (!q.empty()) {
                    auto [current, path] = q.front();
                    q.pop();

                    const Atom* currentAtom = mol.getAtomWithIdx(current);
                    for (const Atom* neighbor : mol.atomNeighbors(currentAtom)) {
                        int neighborIdx = neighbor->getIdx();

                        if (visited.find(neighborIdx) != visited.end()) continue;

                        // Check if we've reached the second ring
                        if (ringAtoms.find(neighborIdx) != ringAtoms.end() &&
                            std::find(atomRings[j].begin(), atomRings[j].end(), neighborIdx) != atomRings[j].end()) {
                            // Valid path found
                            linkers.insert(path.begin(), path.end());
                            break;
                        }

                        // Continue BFS through non-ring atoms
                        if (nonRingAtoms.find(neighborIdx) != nonRingAtoms.end()) {
                            std::vector<int> newPath = path;
                            newPath.push_back(neighborIdx);
                            q.push({neighborIdx, newPath});
                            visited.insert(neighborIdx);
                        }
                    }
                }
            }
        }

        return linkers;
    }

    // Function to calculate the FMF ratio
    double Framework(const ROMol& mol) {
        const RingInfo* ringInfo = mol.getRingInfo();
        std::unordered_set<int> ringAtoms;

        // Collect all atoms that are part of rings
        std::vector<std::vector<int>> atomRings = ringInfo->atomRings();
        for (const auto& ring : atomRings) {
            for (int atom : ring) {
                ringAtoms.insert(atom);
            }
        }

        // Find linkers using BFS on non-ring atoms
        std::unordered_set<int> linkers = findLinkersWithBFS(mol, ringAtoms);

        // Total number of atoms (including hydrogens)
        std::unique_ptr<RDKit::ROMol> hmol(RDKit::MolOps::addHs(mol));
        int totalAtoms = hmol->getNumAtoms();

        // Number of framework atoms: linkers + ring atoms
        int frameworkAtoms = linkers.size() + ringAtoms.size();

        // Calculate FMF
        double FMF = static_cast<double>(frameworkAtoms) / totalAtoms;

        return FMF;
    }



     std::vector<double> calcFramework(const ROMol& mol) {
        std::vector<double> res(1,0.);
        res[0] = Framework(mol);
        return res;
     }

    // BRStates: Tetko version only organis !

    const std::vector<std::string>& getOrganicBondKeys() {
    static const std::vector<std::string> organicbondkeys = {"7-S-7","9-S-7","11-S-7","13-S-7","15-S-7","16-S-7","17-S-7","19-S-7","20-S-7","21-S-7","22-S-7","24-S-7",
    "27-S-7","28-S-7","30-S-7","31-S-7","33-S-7","34-S-7","36-S-7","38-S-7","48-S-7","50-S-7","52-S-7","53-S-7","54-S-7","70-S-7",
    "75-S-7","39-S-7","8-D-8","11-D-8","14-D-8","16-D-8","23-D-8","28-D-8","35-D-8","49-D-8","52-D-8","9-S-9","11-S-9","13-S-9",
    "15-S-9","16-S-9","17-S-9","19-S-9","20-S-9","21-S-9","22-S-9","24-S-9","27-S-9","28-S-9","30-S-9","31-S-9","33-S-9","34-S-9",
    "36-S-9","38-S-9","48-S-9","50-S-9","52-S-9","53-S-9","54-S-9","70-S-9","75-S-9","39-S-9","10-T-10","15-T-10","26-T-10","11-D-11",
    "11-S-11","13-S-11","14-D-11","15-S-11","16-D-11","16-S-11","17-S-11","19-S-11","20-S-11","21-S-11","22-S-11","23-D-11","24-S-11",
    "27-S-11","28-D-11","28-S-11","30-S-11","31-S-11","33-S-11","34-S-11","35-D-11","36-S-11","38-S-11","48-S-11","49-D-11","50-S-11",
    "52-D-11","52-S-11","54-S-11","70-S-11","75-S-11","39-S-11","12-A-12","17-A-12","18-A-12","25-A-12","29-A-12","32-A-12","37-A-12",
    "13-S-13","15-S-13","16-S-13","17-S-13","19-S-13","20-S-13","21-S-13","22-S-13","24-S-13","27-S-13","28-S-13","30-S-13","31-S-13",
    "33-S-13","34-S-13","36-S-13","38-S-13","48-S-13","50-S-13","52-S-13","53-S-13","54-S-13","70-S-13","75-S-13","39-S-13","14-D-14",
    "16-D-14","23-D-14","28-D-14","35-D-14","49-D-14","52-D-14","53-D-14","15-T-15","15-S-15","17-S-15","19-S-15","20-S-15","21-S-15",
    "22-S-15","24-S-15","26-T-15","27-S-15","28-S-15","30-S-15","31-S-15","33-S-15","34-S-15","36-S-15","38-S-15","48-S-15","50-S-15",
    "52-S-15","53-S-15","54-S-15","70-S-15","75-S-15","39-S-15","16-D-16","16-S-16","17-S-16","19-S-16","20-S-16","21-S-16","22-S-16",
    "23-D-16","24-S-16","27-S-16","28-D-16","28-S-16","30-S-16","31-D-16","34-S-16","35-D-16","36-S-16","38-D-16","48-S-16","49-D-16",
    "50-S-16","52-D-16","52-S-16","53-D-16","53-S-16","54-S-16","70-S-16","75-S-16","39-S-16","17-A-17","17-S-17","18-A-17","19-S-17",
    "20-S-17","21-S-17","22-S-17","24-S-17","28-S-17","30-S-17","31-S-17","33-S-17","34-S-17","36-S-17","37-A-17","38-S-17","48-S-17",
    "50-S-17","51-A-17","52-S-17","53-S-17","54-S-17","70-S-17","75-S-17","39-S-17","18-A-18","25-A-18","29-A-18","32-A-18","37-A-18",
    "51-A-18","19-S-19","20-S-19","21-S-19","22-S-19","24-S-19","27-S-19","28-S-19","30-S-19","31-S-19","33-S-19","34-S-19","36-S-19",
    "38-S-19","48-S-19","50-S-19","52-S-19","53-S-19","54-S-19","70-S-19","75-S-19","39-S-19","24-S-20","28-S-20","30-S-20","50-S-20",
    "52-S-20","53-S-20","21-S-21","24-S-21","28-S-21","52-S-21","53-S-21","24-S-22","28-S-22","30-S-22","50-S-22","52-S-22","53-S-22",
    "28-D-23","28-S-23","52-D-23","52-S-23","53-D-23","53-D-23","24-S-24","28-S-24","52-S-24","53-S-24","25-A-25","29-S-25","37-A-25",
    "51-A-25","28-S-27","30-S-27","50-S-27","52-S-27","53-S-27","28-D-28","28-S-28","30-S-28","31-S-28","35-D-28","36-S-28","38-S-28",
    "49-D-28","52-D-28","52-S-28","53-D-28","53-S-28","29-A-29","37-A-29","51-A-29","30-S-30","31-S-30","34-S-30","36-S-30","37-S-30",
    "52-S-30","53-S-30","31-S-31","35-D-31","35-S-32","50-S-33","52-S-33","53-S-33","34-S-34","36-S-34","52-S-34","53-S-34","52-D-35",
    "53-D-35","39-S-36","50-S-36","52-S-36","53-S-36","54-S-36","70-S-36","75-S-36","37-A-37","51-A-37","38-S-38","39-S-38","50-S-38",
    "52-S-38","53-S-38","54-S-38","70-S-38","75-S-38","39-S-39","48-S-39","50-S-39","52-S-39","53-S-39","54-S-39","70-S-39","75-S-39",
    "48-S-48","50-S-48","52-S-48","53-S-48","52-D-49","53-D-49","50-S-50","52-S-50","53-S-50","54-S-50","70-S-50","75-S-50","52-S-52",
    "53-S-52","54-S-52","70-S-52","75-S-52","54-S-54","70-S-54","75-S-54","70-S-70","75-S-70","75-S-75"};
            return organicbondkeys;
    }

    struct BondEStateResult {
        std::vector<std::string> SumKeys;
        std::vector<double> BEStotal;
        std::vector<double> SumBES;
        std::vector<double> nBES;
        std::vector<double> minBES;
        std::vector<double> maxBES;
    };


    // retreive the positional of the
    std::vector<int> NamePosES(const RDKit::ROMol &mol, const std::vector<std::pair<std::string, std::shared_ptr<RDKit::RWMol>>> &queries) {
        size_t nAtoms = mol.getNumAtoms();
        std::vector<int> pos(nAtoms, 0);  // Initialize positions with 0
        for (unsigned int idx=0; idx<queries.size(); idx++) {
            auto entry = queries[idx];
            if (!entry.second) continue;  // Skip invalid SMARTS patterns

            std::vector<RDKit::MatchVectType> qmatches;
            if ( RDKit::SubstructMatch(mol, *entry.second, qmatches, true) )
            {
                for(unsigned int i = 0 ; i < qmatches.size() ; ++i )
                {       int atomIdx = qmatches[i][0].second;
                        //std::cout << entry.first << ":" << atomIdx << ": " << idx+1 << "\n";
                        pos[atomIdx] = idx+1;
                }
            }
        }

        /*std::cout << "C++ nPatts: " << queries.size() << std::endl;
        std::cout << "pos : ";
        for (int i=0; i< pos.size();i++) {
            std::cout << pos[i] << " ";
        }
        std::cout << "\n";
        */
        return pos;
    }

    // Tetko version
    BondEStateResult getBEStateFeatures(const RDKit::ROMol &mol, bool extended) {
        size_t nBonds = mol.getNumBonds();
        size_t nAtoms = mol.getNumAtoms();

        // Precompute atomic EState indices
        std::vector<double> Is = calcIStateIndices(mol);

        /*std::cout << "IStatesIndices : ";
        for (int i=0; i< Is.size();i++) {
            std::cout << Is[i] << " ";
        }
        std::cout << "\n";
        */
        // Get the distance matrix using RDKit's MolOps::getDistanceMat
        double* dists = MolOps::getDistanceMat(mol, false, false, false); // no bond order, no weights, no hydrogens

        /*
        for (int i=0; i<nAtoms; i++){
            for (int j=0; j<nAtoms; j++){
                std::cout << dists[j * nAtoms + i]+1. << ", ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
        */


        // Bond-specific indices
        std::vector<double> Iij(nBonds, 0.0);
        std::vector<std::string> BEScode(nBonds);

        const auto &queries = extended ? GetesExtQueries() : GetesQueries();
        const std::vector<int> pos = NamePosES(mol, queries);

        // Compute bond contributions
        for (size_t t = 0; t < nBonds; ++t) {
            const RDKit::Bond *bt = mol.getBondWithIdx(t);
            int i = bt->getBeginAtomIdx();
            int j = bt->getEndAtomIdx();

            // Bond type and code
            std::string btype;
            switch (bt->getBondType()) {
                case RDKit::Bond::SINGLE: btype = "S"; break;
                case RDKit::Bond::DOUBLE: btype = "D"; break;
                case RDKit::Bond::TRIPLE: btype = "T"; break;
                default: btype = "A"; break;
            }

            // Calculate Iij
            Iij[t] = 0.5*(Is[i] + Is[j]);

            // Generate bond code
            int posi = pos[i];
            int posj = pos[j];
            if (posi >= posj) {
                BEScode[t] = std::to_string(posi) + "-" + btype + "-" + std::to_string(posj);
            } else {
                BEScode[t] = std::to_string(posj) + "-" + btype + "-" + std::to_string(posi);
            }
        }

        // Initialize BES
        std::vector<double> BES(nBonds, 0.0);

        // Compute edge EState contributions
        for (size_t i = 0; i < nBonds; ++i) {
            const RDKit::Bond *bt1 = mol.getBondWithIdx(i);
            int i1 = bt1->getBeginAtomIdx();
            int i2 = bt1->getEndAtomIdx();

            for (size_t j = 0; j < i; ++j) {
                const RDKit::Bond *bt2 = mol.getBondWithIdx(j);
                int j1 = bt2->getBeginAtomIdx();
                int j2 = bt2->getEndAtomIdx();
                double d11 = dists[i1 * nAtoms + j1];
                double d12 = dists[i1 * nAtoms + j2];
                double d21 = dists[i2 * nAtoms + j1];
                double d22 = dists[i2 * nAtoms + j2];

                // Topological distance adding one or ???
                double p = 1+(d11+d12+d21+d22)/4;
                double dI = (Iij[i] - Iij[j]) / (p * p);
                BES[i] += dI;
                BES[j] -= dI;
            }
        }



        // Total EState contributions
        std::vector<double> BEStotal(nBonds);
        for (size_t i = 0; i < nBonds; ++i) {
            BEStotal[i] = BES[i] + Iij[i];
        }
        /*
        std::cout << "Iij : ";
        for (int i=0; i< Iij.size();i++) {
            std::cout << Iij[i] << " ";
        }
        std::cout << "\n";


        std::cout << "BES : ";
        for (int i=0; i< BES.size();i++) {
            std::cout << BES[i] << " ";
        }
        std::cout << "\n";


        std::cout << "BEStotal : ";
        for (int i=0; i< BEStotal.size();i++) {
            std::cout << BEStotal[i] << " ";
        }
        std::cout << "\n";
        */

        // Summ contributions by bond code
        std::unordered_map<std::string, std::vector<double>> codeMap;
        for (size_t i = 0; i < nBonds; ++i) {
            codeMap[BEScode[i]].push_back(BEStotal[i]);
        }

        // Aggregated results
        std::vector<std::string> SumKeys;
        std::vector<double> SumBES, minBES, maxBES, nBES;

        for (const auto &[key, values] : codeMap) {
            SumKeys.push_back("S" + key);
            SumBES.push_back(std::accumulate(values.begin(), values.end(), 0.0));
            nBES.push_back(static_cast<double>(values.size()));
            minBES.push_back(*std::min_element(values.begin(), values.end()));
            maxBES.push_back(*std::max_element(values.begin(), values.end()));


        }
        /*
        std::sort(SumKeys.begin(), SumKeys.end());
        for (const auto &k : SumKeys) {
            std::cout << k << " ";
        }
        std::cout << "\n";
        */
        return {SumKeys, BEStotal, SumBES, nBES, minBES, maxBES};
    }




    // Function to compute Bond E-State fingerprints
    std::vector<double> calcBEStateDescs(const RDKit::ROMol &mol) {
        // Call the function to calculate Bond E-State descriptors using the extended patterns (aka true)
        auto [SumKeys_i, BEStotal_i, SumBES_i, nBES_i, minBES_i, maxBES_i] = getBEStateFeatures(mol, true);


        const auto&  orgbondkeys = getOrganicBondKeys();


        // Initialize results with size equal to organicbondkeys plus one (for unmatched patterns)
        size_t nKeys = orgbondkeys.size();

        std::vector<double> sumsBES(nKeys + 1, 0.);
        std::vector<double> nBES(nKeys + 1, 0.);
        std::vector<double> minBES(nKeys + 1, 99999);
        std::vector<double> maxBES(nKeys + 1, 0.);

        // Process each descriptor
        for (size_t i = 0; i < SumBES_i.size(); ++i) {

            //std::cout << "Key: " << SumKeys_i[i] << " BEStotal: " << BEStotal_i[i] << " SumBES:" << SumBES_i[i] << " nBES:" << nBES_i[i] << " minBES:" << minBES_i[i] << " maxBES:" << maxBES_i[i] << "\n";

            const auto &pattern = SumKeys_i[i];
            const auto &descriptor = pattern.substr(1); // Extract the key after "S" only one to remove!


            // Check if the descriptor exists in organicbondkeys
            auto it = std::find(orgbondkeys.begin(), orgbondkeys.end(), descriptor);
            if (it != orgbondkeys.end()) {

                // Get the position index
                size_t posidx = std::distance(orgbondkeys.begin(), it);

                //std::cout << "found at position: " << posidx <<"\n";

                if (posidx < sumsBES.size()) {

                    // Update the corresponding values
                    sumsBES[posidx] += SumBES_i[i];
                    nBES[posidx] += nBES_i[i];
                    minBES[posidx] = std::min(minBES[posidx], minBES_i[i]);
                    maxBES[posidx] = std::max(maxBES[posidx], maxBES_i[i]);
                }
                else {
                    std::cerr << "Out-of-bounds access detected for posIdx: " << posidx << " : " << sumsBES.size() << "\n";

                }

            } else {
                // Update the last index (unmatched patterns)
                size_t unmatchedIdx = nKeys;
                sumsBES[unmatchedIdx] += SumBES_i[i];
                nBES[unmatchedIdx] += nBES_i[i];
                minBES[unmatchedIdx] = std::min(minBES[unmatchedIdx], minBES_i[i]) ;
                maxBES[unmatchedIdx] = std::max(maxBES[unmatchedIdx], maxBES_i[i]);
            }
        }


        for (unsigned int i=0; i<nKeys+1; i++) {
            if (minBES[i]==99999) {
                minBES[i] = 0.;
            }
        }


        // Concatenate all vectors into a single vector<double>
        std::vector<double> concatenatedResult;
        concatenatedResult.reserve(sumsBES.size() + nBES.size() + minBES.size() + maxBES.size());

        concatenatedResult.insert(concatenatedResult.end(), sumsBES.begin(), sumsBES.end());
        concatenatedResult.insert(concatenatedResult.end(), nBES.begin(), nBES.end());
        concatenatedResult.insert(concatenatedResult.end(), minBES.begin(), minBES.end());
        concatenatedResult.insert(concatenatedResult.end(), maxBES.begin(), maxBES.end());

        return concatenatedResult;
    }


    static const std::vector<std::string> AFragments = {
        "[C][OX2H]","[c][OX2H]","[C][NX3;H2]","[c][NX3;H2;!$(NC=O)]","[C][NX3;H1;!R][C]","[C][NX3;H1;R][C]","[c][NX3;H1;!$(NC=O)][C]","[c][nX3;H1][c]","[CX3](=O)[OX1H0-,OX2H1]",
        "[CX3](=[OX1])[NX3;H2]","[CX3](=[OX1])[NX3;H1][C]","[CX3](=[OX1])[NX3;H1][c]","[$([SX4](=[OX1])(=[OX1])([!O])[NH,NH2,NH3+]),$([SX4+2]([OX1-])([OX1-])([!O])[NH,NH2,NH3+])]",
        "[NX3;H1]C(=[OX1])[NX3;H1]","[NX3;H0]C(=[OX1])[NX3;H1]","[NX3;H1]C(=[OX1])O","[NX3;H1]C(=N)[NX3;H0]","[C]#[CH]","P[OH,O-]","[CH][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]",
        "[CH]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]",
        "[CX4]([CX3](=O)[OX1H0-,OX2H1])[CX4][CX3](=O)[OX1H0-,OX2H1]","[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX3](=O)[OX1H0-,OX2H1]",
        "[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[OH]","[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][OH]",
        "[nX3;H1]:n","[nX3;H1]:c:n","[OX2;H1]CC[O,N]","[OX2;H1]C[C,N]=[O,S]","[OX2;H1]c1ccccc1[O,NX3]","[OX2;H1]c1ccccc1C=[O,S]","[OX2;H1]c1ccccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
        "[NH,NH2,NH3+]CC[O,N]","[NH,NH2,NH3+]c1ccccc1[O,N]","[NH,NH2,NH3+]c1ccccc1[C,N]=[O,S]","[OX2H]c1ccccc1[Cl,Br,I]","[OX1]=[C,c]~[C,c]C[OH]","[OH]c1cccc2cccnc12",
        "[OH]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1","[OH]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1",
        "[NH,NH2,NH3+]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1","[NH,NH2,NH3+]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1",
        "[CX3](=O)([OX1H0-,OX2H1])c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1",
        "[CX3](=O)([OX1H0-,OX2H1])c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1","[OH]c1c([CX4])cccc1[CX4]","[NH,NH2,NH3+]c1c([CX4])cccc1[CX4]",
        "[OH]c1c(C[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cccc1","[OH]c1cc([CX3](=O)[OX1H0-,OX2H1])ccc1","[OH]c1ccc([CX3](=O)[OX1H0-,OX2H1])cc1",
        "[OH]c1cc([$([CH](=O)),$(C(=O)C)])ccc1","[OH]c1ccc([$([CH](=O)),$(C(=O)C)])cc1"};

    static const std::vector<std::string> BSELFragments = {
        "[CX4H3]","[CX4H2]","[CX4H1]","[CX4H0]","*=[CX3H2]","[$(*=[CX3H1]),$([cX3H1](a)a)]","[$(*=[CX3H0]),$([cX3H0](a)(a)A)]","c(a)(a)a","*#C","[C][NX3;H2]","[c][NX3;H2]","[C][NX3;H1][C]",
        "[c][NX3;H1]","[c][nX3;H1][c]","[C][NX3;H0](C)[C]","[c][NX3;H0](C)[C]","[c][nX3;H0][c]","*=[Nv3;!R]","*=[Nv3;R]","[nX2H0,nX3H1+](a)a","N#C[A;!#1]","N#C[a;!#1]",
        "[$([A;!#1][NX3](=O)=O),$([A;!#1][NX3+](=O)[O-])]","[$([a;!#1][NX3](=O)=O),$([a;!#1][NX3+](=O)[O-])]","[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]","[OH]","[OX2;H0;!R]",
        "[OX2;H0;R]","[oX2](a)a","*=O","[SX2](*)*","[sX2](a)a","*=[SX1]","*=[SX3]","[$([#16X4](=[OX1])(=[OX1])([!#8])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([!#8])[OX2H0])]","[S,s]","[P,p]",
        "FA","Fa","Cl","Br","I","[CX3;!R](=[OX1])[OX2H0]","[CX3;R](=[OX1])[OX2H0;R]","P(=[OX1])(O)(O)O","[CX3](=[OX1])([OX2H0])[OX2H0]","[CX3](=O)[OX1H0-,OX2H1]","nC=[OX1]","[N;!R]C=[OX1]",
        "[N;R][C;R]=[OX1]","[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]","NC(=[OX1])N","[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]","[CX3](=[OX1])[NX3][CX3](=[OX1])",
        "C1(=[OX1])C=CC(=[OX1])C=C1","[$([CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])]",
        "[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]",
        "*1~*2~*(~*3~*(~*~*~*~*3)~*1)~*~*~*1~*2~*~*~*1","[OX2H]CC[O,N]","[OX2H]C[C,N]=[O,S]","[OX2H]c1ccccc1[O,Nv3]","[OX2H]c1ccccc1C=[O,S]","[OX2H]c1ccccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
        "[NH,NH2,NH3+]CC[O,N]","[NH,NH2,NH3+]c1ccccc1[O,N]","[NH,NH2,NH3+]c1ccccc1[C,N]=[O,S]","[OX2H]c1ccccc1[Cl,Br,I]","[CX4]([OH])[CX4][OH]","n:n","o:n","n:c:n","o:c:n","n:c:c:n",
        "[F,Cl,Br,I,N,O,S]-c:c-[F,Cl,Br,I,N,O,S]","[F,Cl,Br,I,N,O,S]-c:c:c-[F,Cl,Br,I,N,O,S]","[F,Cl,Br,I,N,O,S]-c:c:c:c-[F,Cl,Br,I,N,O,S]","P(=[OX1])N","Nc:n","[$(cC[OH]);!$(c[CX3](=O)[OX1H0-,OX2H1])]",
        "[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]","[OX2]-c:c-[OX2]"};


    // Precompile SMARTS patterns for efficiency
    static const std::vector<std::shared_ptr<RDKit::RWMol>> &GetQueriesA() {
      static const std::vector<std::shared_ptr<RDKit::RWMol>> queriesA = [] {
        std::vector<std::shared_ptr<RDKit::RWMol>> res;
        for (const auto &smi : AFragments) {
          auto mol = RDKit::SmartsToMol(smi);
          if (mol) {
            res.emplace_back(std::move(mol));
          } else {
            std::cerr << "Invalid SMARTS: " << smi << std::endl;
          }
        }
        return res;
      }();
      return queriesA;
    }

    static const std::vector<std::shared_ptr<RDKit::RWMol>> &GetQueriesB() {
      static const std::vector<std::shared_ptr<RDKit::RWMol>> queriesB = [] {
        std::vector<std::shared_ptr<RDKit::RWMol>> res;
        for (const auto &smi : BSELFragments) {
          auto mol = RDKit::SmartsToMol(smi);
          if (mol) {
            res.emplace_back(std::move(mol));
          } else {
            std::cerr << "Invalid SMARTS: " << smi << std::endl;
          }
        }
        return res;
      }();
      return queriesB;
    }

    static const std::vector<double> coefAFragments = {
        0.345,0.543,0.177,0.247,0.087,0.321,0.194,0.371,0.243,0.275,0.281,-0.091,0.356,-0.165,-0.119,-0.105,0.170,0.082,0.493,0.019,0.050,-0.362,0.118,0.1,0.051,0.194,0.042,-0.089,-0.161,
        -0.251,-0.418,-0.45,-0.155,0.,-0.093,-0.11,-0.601,-0.475,0.119,0.176,0.08,0.084,0.085,0.055,-0.162,-0.181,0.195,-0.203,0.096,0.185,0.203,0.003 };

    static const std::vector<double> coefBHFragments = {
        0.007,0.,0.011,0.037,0.019,0.011,0.,0.019,0.028,0.481,0.275,0.541,0.415,0.316,0.653,0.321,0.392,0.200,0.596,0.321,0.242,0.103,-0.476,-0.525,-0.204,0.307,0.211,0.331,0.047,0.334,
        0.168,0.043,0.071,0.448,-0.188,0.,1.183,-0.036,0.,0.,-0.011,0.,-0.206,-0.214,-0.394,-0.267,-0.308,-0.095,-0.287,-0.231,-0.446,-0.076,-0.252,-0.148,-0.051,-0.014,0.013,0.267,0.,
        -0.068,-0.079,-0.387,-0.126,0.,-0.059,-0.045,-0.130,0.,-0.132,-0.157,-0.098,-0.170,-0.089,0.031,-0.035,-0.023,-0.668,-0.042,0.131,-0.408,-0.216,0.071};

    static const std::vector<double> coefBOFragments = {
        0.000,0.000,0.020,0.047,0.024,0.012,0.000,0.018,0.032,0.486,0.326,0.543,0.426,0.267,0.655,0.338,0.338,0.202,0.589,0.300,0.245,0.093,-0.595,-0.533,-0.202,0.311,0.226,0.330,0.060,
        0.339,0.175,0.083,0.069,0.319,-0.190,0.000,1.189,-0.033,0.000,0.000,0.000,0.000,-0.223,-0.169,-0.408,-0.298,-0.312,-0.038,-0.292,-0.242,-0.443,-0.054,-0.251,-0.149,-0.050,-0.016,
        0.010,0.218,0.000,-0.090,-0.122,-0.403,-0.120,0.000,-0.027,-0.069,-0.130,-0.018,-0.094,-0.141,-0.113,-0.184,-0.073,0.025,-0.033,-0.025,-0.668,-0.057,0.129,-0.405,-0.218,0.064 };

    static const std::vector<double> coefSFragments = {
        -0.075,0.,0.036,0.071,-0.085,0.050,0.101,0.121,0.034,0.175,0.383,0.265,0.311,0.221,0.323,0.295,0.265,0.125,0.254,0.223,0.694,0.390,0.,-0.231,-0.476,0.247,0.185,0.185,0.,
        0.370,0.189,0.,0.618,1.065,-0.505,0.643,0.703,-0.042,0.,0.082,0.161,0.198,-0.225,0.360,-0.240,-0.190,-0.412,-0.076,0.175,-0.1,-0.569,-0.553,-0.588,-0.510,-0.411,-0.050,
        0.000,1.029,-0.067,-0.095,-0.237,-0.344,-0.276,-0.102,0.,-0.140,-0.120,0.052,0.024,0.047,-0.040,0.087,-0.051,-0.043,-0.038,0.,-0.452,0.098,0.,0.434,0.380,0.277 };

    static const std::vector<double> coefEFragments = {
        -0.104, 0.,0.089,0.187,-0.045,0.068,0.18,0.3,0.04,0.085,0.163,0.138,0.192,-0.03,0.22,0.346,0.083,0.117,0.121,0.046,0.,0.,0.2,0.21,0.,0.061,0.014,0.013,-0.125,
        -0.041,0.33,0.116,0.364,0.413,0.,0.465,0.295,-0.18,-0.23,0.023,0.196,0.533,-0.113,0.,-0.1,0.,-0.192,0.221,0.,0.061,-0.111,-0.11,0.,0.,0.,-0.017,0.012,
        0.285,0.029,0.,-0.069,0.,0.,0.,0.,0.,-0.1,-0.043,0.092,-0.113,0.,0.052,0.,0.,0.,0.,-0.08,0.185,0.,0.,0.,0.248 };

    static const std::vector<double> coefLFragments = {
        0.321,0.499,0.449,0.443,0.244,0.469,0.624,0.744,0.332,0.781,0.949,0.568,0.912,1.25,0.4,0.869,0.794,-0.235,-0.24,0.574,0.757,0.732,0.278,0.347,0.,0.672,0.36,0.359,
        0.057,0.495,1.258,0.848,0.954,2.196,0.,0.554,2.051,-0.143,-0.147,0.669,1.097,1.590,-0.39,0.406,-0.483,0.,-0.369,0.,0.603,0.583,0.,0.,0.,0.,0.,-0.111,
        0.054,0.488,-0.072,-0.337,0.,-0.303,-0.364,0.062,0.,0.169,-0.4,0.1,-0.179,0.,0.042,0.209,-0.058,-0.081,-0.026,0.,0.,0.149,-0.145,0.,0.,0.13 };


    std::vector<double> calcAbrahams(const RDKit::ROMol& mol) {
        std::vector<double> retval(6, 0.0);

        try {
            // Calculate A descriptor
            auto queriesA = GetQueriesA();
            for (size_t i = 0; i < queriesA.size(); ++i) {

                std::vector<RDKit::MatchVectType> matches;
                RDKit::SubstructMatch(mol, *queriesA[i], matches, true);  // uniquify = true
                retval[0] += matches.size() * coefAFragments[i];
            }

            // Calculate BSEL descriptors
            int sulphurCount = 0;
            auto queriesB = GetQueriesB();
            for (size_t i = 0; i < queriesB.size(); ++i) {

                std::vector<RDKit::MatchVectType> matches;
                RDKit::SubstructMatch(mol, *queriesB[i], matches, true);  // uniquify = true

                int uniqueMatches = matches.size();
                if (30 <= i && i <= 34) {
                    sulphurCount += uniqueMatches;
                } else if (i == 35) {
                    uniqueMatches -= sulphurCount;
                }

                retval[1] += uniqueMatches * coefBHFragments[i];  // BH
                retval[2] += uniqueMatches * coefBOFragments[i];  // BO
                retval[3] += uniqueMatches * coefSFragments[i];   // S
                retval[4] += uniqueMatches * coefEFragments[i];   // E
                retval[5] += uniqueMatches * coefLFragments[i];   // L
            }

            // Add intercepts
            if (coefAFragments.size() > queriesA.size()) {
                retval[0] += coefAFragments.back();  // A intercept
            }
            if (coefBHFragments.size() > queriesB.size()) {
                retval[1] += coefBHFragments.back(); // BH intercept
                retval[2] += coefBOFragments.back(); // BO intercept
                retval[3] += coefSFragments.back();  // S intercept
                retval[4] += coefEFragments.back();  // E intercept
                retval[5] += coefLFragments.back();  // L intercept
            }



        } catch (const std::exception& e) {
            std::cerr << "Error in SMARTS matching: " << e.what() << std::endl;
            throw std::runtime_error("Error in SMARTSQueryTool");
        }

        return retval;
    }



    //// IC based on the Roy, Basak, Harriss, Magnuson paper Neighorhoo complexities and symmetry of chemical graphs and their biological applications
    // in this algorithm we use the cluster list to iterate per radius not the whole atoms accross clusters.
    // so we don't need to compare all the keys over all atoms but only in a cluster the logic is to split the cluster until we get a singleton cluster or we rich radius 5
    // We don't need to generate a special key for a radius because we compute all the radius in one iterative step.
    // this make the code simpler and sorting and comparison more intuitive.

    // Function to generate a key for an atom's environment based on the paper we don't need the Neighbor degree for a delta radius key cluster separation.
    // but we may need the "truncated leaf" how to encode them ? A truncated leaf dictionnary ?
    // for me the logic would be: truncated leaf key*100 + radius of the truncated leaf to be included in the remaining extension radius...
    // TODO : check if "-2" case is properly used in cluster because by definition an empty key is a key too to discrimitate!

    int generateKey(int rootNum, int rootDeg, int bondOrder, int neighNum, int neighDeg) {
        // osmordredv3: fold in the NEIGHBOR degree too (Mordred parity). The radius-0
        // partition is atomic-number only (built before the r>=1 loop), so adding
        // neighbor degree here only affects r>=1 -> no radius-0 regressions.
        return ((rootNum * 10 + rootDeg) * 1000 + bondOrder * 100 + neighNum) * 10 + neighDeg;
    }


    int getbondtypeint(const Bond::BondType &bd) {
        if (bd == Bond::BondType::AROMATIC ) {
            return 4;
        } else {
            return static_cast<int>(bd);
        }
    }

    // Function to get bond type as double
    double getbondtypeindouble(const RDKit::Bond::BondType &bd) {
        switch (bd) {
            case RDKit::Bond::BondType::SINGLE: return 1.0;
            case RDKit::Bond::BondType::DOUBLE: return 2.0;
            case RDKit::Bond::BondType::TRIPLE: return 3.0;
            case RDKit::Bond::BondType::QUADRUPLE: return 4.0;
            case RDKit::Bond::BondType::AROMATIC: return 1.5;
            default: return static_cast<double>(bd);
        }
    }

    // Initialize adjacency and shortest-path matrices
    std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>>
    initializeMatrixAndSP(int nAtoms, int maxradius) {
        std::vector<std::vector<int>> M(nAtoms, std::vector<int>(nAtoms, -1)); // not sure of N+1 here ???
        std::vector<std::vector<int>> SP(nAtoms, std::vector<int>(maxradius + 1, -1));

        for (int i = 0; i < nAtoms; ++i) {
            M[i][0] = i;
            SP[i][0] = 0;
        }

        return {M, SP};
    }

    // Find the last occupied position in a matrix row
    int findLastOccupied(const std::vector<std::vector<int>>& M, int atomIdx) {
        auto it = std::find(M[atomIdx].begin(), M[atomIdx].end(), -1);
        return it == M[atomIdx].begin() ? 0 : std::distance(M[atomIdx].begin(), it) - 1;
    }

    // Update clusters by grouping atoms with identical keys
    std::vector<std::vector<int>> updateClustersWithKeys(const std::vector<std::pair<int, std::vector<int>>>& keys) {
        std::vector<std::pair<int, std::vector<int>>> sortedKeys = keys;
        std::sort(sortedKeys.begin(), sortedKeys.end(), [](auto& a, auto& b) {
            return a.second < b.second;
        });

        std::vector<std::vector<int>> newClusters;
        std::vector<int> currentCluster;
        std::vector<int> currentKey;

        for (auto& [atomIdx, key] : sortedKeys) {
            if (currentCluster.empty() || key != currentKey) {
                if (!currentCluster.empty()) {
                    newClusters.push_back(currentCluster);
                }
                currentCluster = {atomIdx};
                currentKey = key;
            } else {
                currentCluster.push_back(atomIdx);
            }
        }
        if (!currentCluster.empty()) newClusters.push_back(currentCluster);
        return newClusters;
    }

    // Main pipeline
    std::map<int, std::vector<std::vector<int>>>  computePipeline(RDKit::RWMol& mol, int maxRadius, bool addDeadKeys = false) {
        int nAtoms = rdcast<int>(mol.getNumAtoms());

        if (nAtoms == 0) {
            std::cerr << "Error: Molecule has no atoms." << std::endl;
            return {};
        }

        std::string smi = RDKit::MolToSmiles(mol);

        bool debug= false; // (smi=="FP(F)F" || smi=="BrBr");


        if (debug) {
            std::cerr << "Debugging enabled for molecule: " << smi << "n & NumAtoms: " << nAtoms << std::endl;
        }

        auto [M, SP] = initializeMatrixAndSP(nAtoms, maxRadius);



        if (debug) {
            std::cerr << "Initial M matrix:\n";
            for (const auto& row : M) {
                for (int val : row) {
                    std::cerr << val << " ";
                }
                std::cerr << "\n";
            }

            std::cerr << "Initial SP matrix:\n";
            for (const auto& row : SP) {
                for (int val : row) {
                    std::cerr << val << " ";
                }
                std::cerr << "\n";
            }
        }

        std::map<int, std::vector<std::vector<int>>> CN; // Combined CN and AN
        std::map<int, std::vector<int>> clustersByAN;

        // Radius 0: Group atoms by atomic number
        for (auto& atom : mol.atoms()) {
            clustersByAN[atom->getAtomicNum()].push_back(atom->getIdx());
        }

        std::vector<std::vector<int>> clusters;
        for (auto& [_, atoms] : clustersByAN) {
            clusters.push_back(atoms);
        }

        // Initialize CN[0]
        CN[0].resize(2); // Two vectors: sizes and last atomic values
        for (auto& cluster : clusters) {
            if (cluster.empty()) continue;  // Skip empty clusters
            CN[0][0].push_back(static_cast<int>(cluster.size()));  // Cluster size
            int atomIdx = cluster.back();
            if (atomIdx < 0 || atomIdx >= nAtoms) {
                std::cerr << "Error: Invalid atom index: " << atomIdx << std::endl;
                continue;
            }
            CN[0][1].push_back(mol.getAtomWithIdx(atomIdx)->getAtomicNum());   // Last atomic number


        }

        if (debug) {
            std::cerr << "\nM matrix before radius 1:\n";
            for (const auto& row : M) {
                for (int val : row) {
                    std::cerr << val << " ";
                }
                std::cerr << "\n";
            }

            std::cerr << "\nSP matrix before radius 1:\n";
            for (const auto& row : SP) {
                for (int val : row) {
                    std::cerr << val << " ";
                }
                std::cerr << "\n";
            }
        }

        for (int r = 1; r <= maxRadius; ++r) {
            bool stopExpansion = true;

            std::vector<std::vector<int>> newClusters;

            for (auto& cluster : clusters) {
                if (cluster.size() == 1) {
                    newClusters.push_back(cluster);
                    continue;
                }

                std::vector<std::pair<int, std::vector<int>>> clusterKeys;

                for (int atomIdx : cluster) {
                    int start = SP[atomIdx][r - 1]; // get the start neighbor to visite at this radius
                    int stop = findLastOccupied(M, atomIdx); // get the end


                    // osmordredv3 FIX (keep-all): when an atom's BFS frontier is
                    // exhausted (SP[atomIdx][r-1] == -2) it has no new neighbors, but it
                    // must STAY in the partition. Previously this `continue` dropped it
                    // BEFORE it was appended to clusterKeys, so atoms vanished as the
                    // radius grew, classes merged and IC/TIC/CIC/SIC/BIC came out too low
                    // (e.g. Benzene IC5 collapsed to 0). Keep it with a stable empty key so
                    // CN[r] sizes always sum to N. Proven in Python on a faithful
                    // computePipeline port (matches the binary 840/840): fixes the
                    // collapse, 0 regressions vs the Mordred references.
                    if (start == -2) {
                        clusterKeys.emplace_back(atomIdx, std::vector<int>{});  // keep atom; no new frontier
                        continue;
                    }

                    std::vector<int> eqKeys;
                    std::vector<int> neighbors;

                    if (debug) {
                        std::cerr << "Radius " << r << ", Atom " << atomIdx << " .Symbol: " << mol.getAtomWithIdx(atomIdx)->getSymbol()
                               << ", Start " << start << ", Stop " << stop << std::endl;
                    }

                    for (int pos = start; pos <= stop; ++pos) {
                        int rootIdx = M[atomIdx][pos];
                        if (rootIdx<0 || rootIdx >= nAtoms) {
                            std::cout << "Atom index out of boundaries:" << rootIdx << "smile:" << smi << "\n" ;
                            continue;
                        }
                        const Atom* rootAtom = mol.getAtomWithIdx(rootIdx);
                        int rootNum = rootAtom->getAtomicNum();
                        int rootDeg = rootAtom->getDegree();

                        for (const auto& nb : mol.atomNeighbors(rootAtom)) {
                            int nbIdx = nb->getIdx();
                            if (std::find(M[atomIdx].begin(), M[atomIdx].begin() + stop + 1, nbIdx) != M[atomIdx].begin() + stop + 1) {
                                continue; // Already visited no effect of adding truncated branch this was already seen in previous radius!
                            }

                            neighbors.push_back(nbIdx);
                            const Bond* bond = mol.getBondBetweenAtoms(rootIdx, nbIdx);
                            int bondOrder = getbondtypeint(bond->getBondType()); // don't need kekulize like in Mordred
                            int neighNum = mol.getAtomWithIdx(nbIdx)->getAtomicNum();
                            int neighDeg = mol.getAtomWithIdx(nbIdx)->getDegree();  // osmordredv3: neighbor degree (Mordred parity)
                            eqKeys.push_back(generateKey(rootNum, rootDeg, bondOrder, neighNum, neighDeg));
                        }
                    }

                    // option one for the deadkeys
                    if (neighbors.empty()) {
                        if (SP[atomIdx][r] == -1) {
                            SP[atomIdx][r] = -2;  // Mark as exhausted
                        }
                        // Add dead key if enabled
                        if (addDeadKeys) {
                            eqKeys.push_back(-2);
                        }
                    }

                    std::sort(eqKeys.begin(), eqKeys.end());
                    clusterKeys.emplace_back(atomIdx, eqKeys);

                    // Appends only validated visited new neighbors, not all neighbors (i.e., exclude backward, leaf, or cycle paths)
                    if (!neighbors.empty()) {
                        stopExpansion = false;  // Continue expansion if new neighbors are found
                        for (int i = 0; i < static_cast<int>(neighbors.size()); ++i) {
                            int freeSlot = findLastOccupied(M, atomIdx) + 1;
                            if (i == 0) SP[atomIdx][r] = (freeSlot < nAtoms) ? freeSlot : -2;
                            if (freeSlot < nAtoms) {
                                M[atomIdx][freeSlot] = neighbors[i];
                            }
                        }
                    }
                    // option two for the deadkeys
                    //else {
                    //    if (SP[atomIdx][r] == -1) {
                    //        SP[atomIdx][r] = -2;  // Mark as exhausted, no further expansion
                    //    }
                    //}
                }

                auto subClusters = updateClustersWithKeys(clusterKeys);
                newClusters.insert(newClusters.end(), subClusters.begin(), subClusters.end());
            }

            clusters = newClusters;

            if (debug) {
                std::cerr << "M Matrix after radius " << r << ":\n";
                for (const auto& row : M) {
                    for (int val : row) {
                        std::cerr << val << " ";
                    }
                    std::cerr << std::endl;
                }

                std::cerr << "SP Matrix after radius " << r << ":\n";
                for (const auto& row : SP) {
                    for (int val : row) {
                        std::cerr << val << " ";
                    }
                    std::cerr << std::endl;
                }
            }


            CN[r].resize(2); // Two vectors: sizes and last atomic values
            for (auto& cluster : clusters) {
                if (cluster.empty()) continue;  // Skip empty clusters
                int atomIdx = cluster.back();  // Atom index

                if (atomIdx<0 || atomIdx >= nAtoms) {
                    std::cout << "Atom index out of boundaries:" << atomIdx << "\n" ;
                    continue;
                }

                int atomicNum = mol.getAtomWithIdx(atomIdx)->getAtomicNum();   // Last atomic mass

                CN[r][0].push_back(static_cast<int>(cluster.size()));  // store Cluster size
                CN[r][1].push_back(atomicNum);   // store atomic mass

            }

            if (stopExpansion ||
                std::all_of(clusters.begin(), clusters.end(), [](auto& c) { return c.size() == 1; })) {

                if (debug) {
                        std::cerr << "Stopping expansion at radius " << r
                                << " - Reason: "
                                << (stopExpansion ? "No new neighbors" : "All clusters are singletons")
                                << std::endl;
                    }

                break;
            }
        }


        // Finalize CN by padding remaining radii

        if (!CN.empty()) {
            auto lastRadius = CN.rbegin()->first;
            for (int rr = lastRadius + 1; rr <= maxRadius; ++rr) {
                CN[rr] = CN[lastRadius];
            }
        }

        if (debug) {
            std::cerr << "Final CN values for " << smi << ":\n";
            for (const auto& [r, values] : CN) {
                std::cerr << "Radius " << r << ": ";
                for (const auto& val : values[0]) {
                    std::cerr << val << " ";
                }
                std::cerr << std::endl;
            }
        }

        // clean up the memory!
        //M.clear();
        //SP.clear();
        //clustersByAN.clear();
        //clusters.clear();

        return CN;
    }


 std::vector<double> ShannonEntropies(std::map<int, std::vector<std::vector<int>>> CN, int maxradius, double log2nA, double log2nB, int nAtoms) {

        const auto* tbl = RDKit::PeriodicTable::getTable();

        int singleoutputsize = maxradius+1;

        std::vector<double> icvalues(7*singleoutputsize, 0.);

        for (auto& [r, data]: CN) {

            icvalues[r] =   InfoEntropy(data[0]);// IC
            icvalues[r+singleoutputsize ] = icvalues[r] * nAtoms; // TIC
            if (log2nA > 0) {
                icvalues[r+2*singleoutputsize] =icvalues[r] / log2nA; // SIC
            }

            if (log2nB>0) {

                icvalues[r+3*singleoutputsize] = icvalues[r] / log2nB; // BIC

            }

            icvalues[r+4*singleoutputsize] = log2nA - icvalues[r]; // CIC

            std::vector<double> w(data[1].size(),0.);


            for (unsigned int j=0; j < data[1].size(); j++) {
                w[j] = tbl->getAtomicWeight(data[1][j]); // catch Atomic mass
            }

            icvalues[r+5*singleoutputsize] = WeightedInfoEntropy(data[0],w); // MIC  use the atomic mass weighted Shannon Entropy
            icvalues[r+6*singleoutputsize] = WeightedCrossInfoEntropy(data[0],data[1]); // ZMIC  use the atomic mass weighted cross Shannon Entropy
            w.clear();

        }

        return icvalues;

 }

    std::vector<double> calcInformationContent(const RDKit::ROMol &mol, int maxradius)
    {

        std::unique_ptr<RDKit::RWMol> hmol(new RDKit::RWMol(mol));
        RDKit::MolOps::addHs(*hmol);

        int nAtoms = hmol->getNumAtoms();

        if (nAtoms == 0) {
            std::cerr << "Error: Molecule has no atoms after adding hydrogens." << std::endl;
            return {};
        }

        double nBonds = 0.;
        for (auto& bond : hmol->bonds()) {
            nBonds += getbondtypeindouble(bond->getBondType());
        }

        double log2nA = std::log(static_cast<double>(nAtoms)) / std::log(2);
        double log2nB = std::log(static_cast<double>(nBonds)) / std::log(2);


        auto CN = computePipeline(*hmol, maxradius);

        if (CN.empty()) {
            std::cerr << "Error: ComputePipeline returned empty CN." << std::endl;
            return {};
        }

        return ShannonEntropies(CN, maxradius, log2nA, log2nB,  nAtoms);
    }





std::vector<double> calcInformationContent_(const RDKit::ROMol& mol) {
    int maxradius = 5;
    // Dynamically allocate RWMol using new
    RDKit::RWMol* hmol = new RDKit::RWMol(mol);

    try {
        // Add hydrogens
        RDKit::MolOps::addHs(*hmol);

        int nAtoms = hmol->getNumAtoms();
        if (nAtoms == 0) {
            std::cerr << "Error: Molecule has no atoms after adding hydrogens." << std::endl;
            delete hmol; // Clean up memory
            return {};
        }

        double nBonds = 0.0;
        for (auto& bond : hmol->bonds()) {
            nBonds += getbondtypeindouble(bond->getBondType());
        }

        double log2nA = std::log(static_cast<double>(nAtoms)) / std::log(2);
        double log2nB = (nBonds > 0) ? std::log(static_cast<double>(nBonds)) / std::log(2) : 0;

        auto CN = computePipeline(*hmol, maxradius);
        if (CN.empty()) {
            std::cerr << "Error: ComputePipeline returned empty CN." << std::endl;
            delete hmol; // Clean up memory
            return {};
        }

        // Calculate Shannon Entropies
        std::vector<double> icvalues = ShannonEntropies(CN, maxradius, log2nA, log2nB, nAtoms);

        delete hmol; // Clean up memory
        return icvalues;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        delete hmol; // Clean up memory
        return {};
    }
}

// triplet example AZ

std::vector<double> TIn(const RDKit::ROMol &mol, const std::vector<double> b) {
    double ti1 = std::accumulate(b.begin(), b.end(), 0.0);
    double ti2 = std::accumulate(b.begin(), b.end(), 0.0, [](double acc, double yi) { return acc + yi * yi; });
    double ti3 = std::accumulate(b.begin(), b.end(), 0.0, [](double acc, double yi) { return acc + std::sqrt(yi); });

    double ti4 = 0;
    for (unsigned int i = 0; i < b.size(); ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            const auto* bond = mol.getBondBetweenAtoms(i, j);
            if (bond != nullptr) {
                ti4 += std::pow(b[i] * b[j], -0.5);
            }
        }
    }

    double ti5 = b.size() * std::pow(std::accumulate(b.begin(), b.end(), 1.0, std::multiplies<double>()), 1.0 / b.size());

    return   {ti1, ti2, ti3, ti4, ti5};
}



// triplet example AZ

std::vector<double> TIn(const RDKit::ROMol &mol, const std::vector<double> b, int nhrs=1) {
    int n = mol.getNumAtoms();
    std::vector res(5 * nhrs, 0.0);
    for (int i = 0; i < nhrs; i++) {
        std::vector<double> bi(b.begin() + i * n, b.begin() + (i + 1) * n);

        double ti1 = std::accumulate(bi.begin(), bi.end(), 0.0);
        double ti2 = std::accumulate(bi.begin(), bi.end(), 0.0, [](double acc, double yi) { return acc + yi * yi; });
        double ti3 = std::accumulate(bi.begin(), bi.end(), 0.0, [](double acc, double yi) { return acc + std::sqrt(yi); });

        double ti4 = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < i; ++j) {
                const auto* bond = mol.getBondBetweenAtoms(i, j);
                if (bond != nullptr) {
                    ti4 += std::pow(bi[i] * bi[j], -0.5);
                }
            }
        }

        double ti5 = n * std::pow(std::accumulate(bi.begin(), bi.end(), 1.0, std::multiplies<double>()), 1.0 / n);
        res[i*5 + 0] = ti1;
        res[i*5 + 1] = ti2;
        res[i*5 + 2] = ti3;
        res[i*5 + 3] = ti4;
        res[i*5 + 4] = ti5;

    }
    return res;
}





// triplet AN*x = V :  S,V,Z,I,N
std::vector<double> calcANMat(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A
    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"

    std::vector<double> DistSum (nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<double>(distances[i * nAtoms + j]); // "S"
        }
    }

    // Modify diagonal to include atomic numbers

    for (int i = 0; i < nAtoms; ++i) {
        adjMatVec[i * nAtoms + i] = static_cast<double>(nAtoms); // Atomic number on diagonal  == Z
    }
    int n = nAtoms, nrhs = 5;

    std::vector<double> V(nAtoms * nrhs, 0.0); // B matrix for 5 right-hand sides

    for (int j = 0; j < nrhs; ++j) { // Iterate over columns
        for (int i = 0; i < nAtoms; ++i) { // Iterate over rows
            switch (j) {
                case 0: V[i + j * nAtoms] = DistSum[i]; break;  // S
                case 1: V[i + j * nAtoms] = static_cast<double>(mol.getAtomWithIdx(i)->getDegree());  break; // V
                case 2: V[i + j * nAtoms] = static_cast<double>(mol.getAtomWithIdx(i)->getAtomicNum());  break; // Z
                case 3: V[i + j * nAtoms] = static_cast<double>(1); break; // I
                case 4: V[i + j * nAtoms] = static_cast<double>(nAtoms);  break; // N
            }
        }
    }


    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);
    bool success = false;

    std::vector<double> B(V);

    solveLinearSystem(mol, A_flat, B, n, nrhs, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,B,nrhs );
}



// triplet AZ*x = V : V,S,N
std::vector<double> calcAZMat(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A
    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"

    std::vector<double> DistSum (nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<double>(distances[i * nAtoms + j]); // "S"
        }
    }

    // Modify diagonal to include atomic numbers

    for (int i = 0; i < nAtoms; ++i) {
        const auto* atom = mol.getAtomWithIdx(i);
        adjMatVec[i * nAtoms + i] = static_cast<double>(atom->getAtomicNum()); // Atomic number on diagonal  == Z
    }
    int n = nAtoms, nrhs = 3;

    std::vector<double> V(nAtoms * nrhs, 0.0); // B matrix for 5 right-hand sides

    for (int j = 0; j < nrhs; ++j) { // Iterate over columns
        for (int i = 0; i < nAtoms; ++i) { // Iterate over rows
            switch (j) {
                case 0: V[i + j * nAtoms] = static_cast<double>(mol.getAtomWithIdx(i)->getDegree());  break; // V
                case 1: V[i + j * nAtoms] = DistSum[i]; break;  // S
                case 2: V[i + j * nAtoms] = static_cast<double>(nAtoms);  break; // N
            }
        }
    }


    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);
    bool success = false;

    std::vector<double> B(V);

    solveLinearSystem(mol, A_flat, B, n, nrhs, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,B,nrhs );
}



// triplet AS*x = V: N,V,Z,I
std::vector<double> calcASMat(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A
    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"

    std::vector<double> DistSum (nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<double>(distances[i * nAtoms + j]); // "S"
        }
    }

    // Modify diagonal to include atomic numbers

    for (int i = 0; i < nAtoms; ++i) {
        adjMatVec[i * nAtoms + i] = static_cast<double>(DistSum[i]); // Atomic number on diagonal  == Z
    }
    int n = nAtoms, nrhs = 4;

    std::vector<double> V(nAtoms * nrhs, 0.0); // B matrix for 5 right-hand sides

    for (int j = 0; j < nrhs; ++j) { // Iterate over columns
        for (int i = 0; i < nAtoms; ++i) { // Iterate over rows
            switch (j) {
                case 0: V[i + j * nAtoms] = static_cast<double>(mol.getAtomWithIdx(i)->getDegree());  break; // V
                case 1: V[i + j * nAtoms] = static_cast<double>(mol.getAtomWithIdx(i)->getAtomicNum()); break;  // Z
                case 2: V[i + j * nAtoms] = static_cast<double>(1); break;  // I
                case 3: V[i + j * nAtoms] = static_cast<double>(nAtoms);  break; // N
            }
        }
    }


    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);
    bool success = false;

    std::vector<double> B(V);

    solveLinearSystem(mol, A_flat, B, n, nrhs, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,B,nrhs );
}



// triplet DS*x = V: V,I,N,Z
std::vector<double> calcDSMat(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0); // "D"
    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"
    std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

    std::vector<double> DistSum (nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<double>(distances[i * nAtoms + j]); // "S"
        }
    }

    // Modify diagonal to include atomic numbers

    for (int i = 0; i < nAtoms; ++i) {
        DistMatVec[i * nAtoms + i] = static_cast<double>(DistSum[i]); // Atomic number on diagonal  == Z
    }
    int n = nAtoms, nrhs = 4;

    std::vector<double> V(nAtoms * nrhs, 0.0); // B matrix for 5 right-hand sides

    for (int j = 0; j < nrhs; ++j) { // Iterate over columns
        for (int i = 0; i < nAtoms; ++i) { // Iterate over rows
            switch (j) {
                case 0: V[i + j * nAtoms] = static_cast<double>(mol.getAtomWithIdx(i)->getDegree());  break; // V
                case 1: V[i + j * nAtoms] = static_cast<double>(1); break;  // I
                case 2: V[i + j * nAtoms] = static_cast<double>(nAtoms);  break; // N
                case 3: V[i + j * nAtoms] = static_cast<double>(mol.getAtomWithIdx(i)->getAtomicNum()); break;  // Z
            }
        }
    }


    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(DistMatVec);
    bool success = false;

    std::vector<double> B(V);

    solveLinearSystem(mol, A_flat, B, n, nrhs, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,B,nrhs );
}



// triplet DN2*x = V: S,I,N,Z
std::vector<double> calcDN2Mat(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0); // "D"
    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"
    std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

    std::vector<double> DistSum (nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<double>(distances[i * nAtoms + j]); // "S"
        }
    }

    // Modify diagonal to include atomic numbers

    for (int i = 0; i < nAtoms; ++i) {
        DistMatVec[i * nAtoms + i] = nAtoms*nAtoms; // Atomic number on diagonal  == Z
    }
    int n = nAtoms, nrhs = 4;

    std::vector<double> V(nAtoms * nrhs, 0.0); // B matrix for 5 right-hand sides

    for (int j = 0; j < nrhs; ++j) { // Iterate over columns
        for (int i = 0; i < nAtoms; ++i) { // Iterate over rows
            switch (j) {
                case 0: V[i + j * nAtoms] = DistSum[i]; break;  // S
                case 1: V[i + j * nAtoms] = static_cast<double>(1); break;  // I
                case 2: V[i + j * nAtoms] = static_cast<double>(nAtoms);  break; // N
                case 3: V[i + j * nAtoms] = static_cast<double>(mol.getAtomWithIdx(i)->getAtomicNum()); break;  // Z

            }
        }
    }


    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(DistMatVec);
    bool success = false;

    std::vector<double> B(V);

    solveLinearSystem(mol, A_flat, B, n, nrhs, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,B,nrhs );
}




// triplet AZ*x = V

std::vector<double> calcAZV(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A
    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());


    // Modify diagonal to include atomic numbers

    for (int i = 0; i < nAtoms; ++i) {
        const auto* atom = mol.getAtomWithIdx(i);
        adjMatVec[i * nAtoms + i] = static_cast<double>(atom->getAtomicNum()); // Atomic number on diagonal  == Z
    }


    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        const auto* atom = mol.getAtomWithIdx(i);
        V[i] = static_cast<double>(atom->getDegree());                        // Degree of the atom  == V
    }

    int n = nAtoms, nrhs = 1;

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);

    std::vector<double> b(V);
    bool success = false;

    solveLinearSystem(mol, A_flat, b, n, nrhs, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b, 1);
}





std::vector<double> calcASV(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A

    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());


    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"

    std::vector<int> DistSum (nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<int>(distances[i * nAtoms + j]); // "S"
        }
    }

    //std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);


    for (int i = 0; i < nAtoms; ++i) {
        const auto* atom = mol.getAtomWithIdx(i);
        adjMatVec[i * nAtoms + i] = static_cast<double>(DistSum[i]); // Distance Sum == S
        V[i] = static_cast<double>(atom->getDegree());   // Degree of the atom  == V
    }

    int n = nAtoms;

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);

    std::vector<double> b(V);
    bool success = false;


    solveLinearSystem(mol, A_flat, b, n, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}


std::vector<double> calcDSV(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0); // "D"
    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"
    std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

    std::vector<int> DistSum (nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<int>(distances[i * nAtoms + j]); // "S"
        }
    }
    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        const auto* atom = mol.getAtomWithIdx(i);
        DistMatVec[i * nAtoms + i] = static_cast<double>(DistSum[i]); // Distance Sum == S
        V[i] = static_cast<double>(atom->getDegree());   // Degree of the atom  == V
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(DistMatVec);

    std::vector<double> b(V);
    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}


std::vector<double> calcAZS(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A

    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());


    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"

    std::vector<double> DistSum(nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<double>(distances[i * nAtoms + j]); // "S"
        }
    }

    //std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        const auto* atom = mol.getAtomWithIdx(i);
        adjMatVec[i * nAtoms + i] = static_cast<double>(atom->getAtomicNum()); // Z
        //V[i] = static_cast<double>(atom->getDegree());   // Degree of the atom  == V
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);

    std::vector<double> b(DistSum); // S

    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}


std::vector<double> calcASZ(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A

    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());


    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"

    std::vector<double> DistSum(nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<double>(distances[i * nAtoms + j]); // "S"
        }
    }

    //std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        const auto* atom = mol.getAtomWithIdx(i);
        adjMatVec[i * nAtoms + i] = static_cast<double>(DistSum[i] ); // S
        V[i] = static_cast<double>(atom->getAtomicNum());   // Z
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);

    std::vector<double> b(V); // Z
    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}



std::vector<double> calcDN2S(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0); // D
    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"
    std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

    std::vector<double> DistSum (nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<double>(distances[i * nAtoms + j]); // "S"
        }
    }
    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        DistMatVec[i * nAtoms + i] = nAtoms*nAtoms; // N2
        //V[i] = static_cast<double>(atom->getDegree());   // Degree of the atom  == V
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(DistMatVec);

    std::vector<double> b(DistSum); // S
    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}


std::vector<double> calcDN2I(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0); // D
    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"
    std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());


    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        DistMatVec[i * nAtoms + i] = nAtoms*nAtoms; // N2
        V[i] = static_cast<double>(1);   // I
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(DistMatVec);

    std::vector<double> b(V); // I
    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}




std::vector<double> calcASI(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A

    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());


    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"

    std::vector<double> DistSum(nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<double>(distances[i * nAtoms + j]); // "S"
        }
    }

    //std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        adjMatVec[i * nAtoms + i] = static_cast<double>(DistSum[i] ); // S
        V[i] = static_cast<double>(1);   // I
    }

     // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);

    std::vector<double> b(V); // Z
    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}


std::vector<double> calcDSI(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0); // "D"
    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"
    std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

    std::vector<int> DistSum (nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<int>(distances[i * nAtoms + j]); // "S"
        }
    }
    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        DistMatVec[i * nAtoms + i] = static_cast<double>(DistSum[i]); // Distance Sum == S
        V[i] = static_cast<double>(1);   // I
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(DistMatVec);

    std::vector<double> b(V);
    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}




std::vector<double> calcASN(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A

    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());


    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"

    std::vector<double> DistSum(nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<double>(distances[i * nAtoms + j]); // "S"
        }
    }

    //std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        adjMatVec[i * nAtoms + i] = static_cast<double>(DistSum[i] ); // S
        V[i] = static_cast<double>(nAtoms);   // N
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);

    std::vector<double> b(V); // Z
    bool success = false;


    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}


std::vector<double> calcDSN(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0); // "D"
    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"
    std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

    std::vector<int> DistSum (nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<int>(distances[i * nAtoms + j]); // "S"
        }
    }
    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        DistMatVec[i * nAtoms + i] = static_cast<double>(DistSum[i]); // Distance Sum == S
        V[i] = static_cast<double>(nAtoms);   // N
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(DistMatVec);

    std::vector<double> b(V);
    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}



std::vector<double> calcDN2N(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0); // D
    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"
    std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());


    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        DistMatVec[i * nAtoms + i] = nAtoms*nAtoms; // N2
        V[i] = static_cast<double>(nAtoms);   // N
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(DistMatVec);

    std::vector<double> b(V); // I
    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}




std::vector<double> calcANS(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A

    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());


    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"

    std::vector<double> DistSum(nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<double>(distances[i * nAtoms + j]); // "S"
        }
    }

    //std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        adjMatVec[i * nAtoms + i] = static_cast<double>( nAtoms ); // N
        V[i] = static_cast<double>(DistSum[i] );   // S
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);

    std::vector<double> b(V); // Z
    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}




std::vector<double> calcANV(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A

    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());


    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);


    for (int i = 0; i < nAtoms; ++i) {
        const auto* atom = mol.getAtomWithIdx(i);
        adjMatVec[i * nAtoms + i] = static_cast<double>( nAtoms ); // N
        V[i] = static_cast<double>( atom->getDegree() );   //  V
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);

    std::vector<double> b(V); // V

    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}


std::vector<double> calcAZN(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A
    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());


    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
      const auto* atom = mol.getAtomWithIdx(i);
        adjMatVec[i * nAtoms + i] = static_cast<double>(atom->getAtomicNum()); // Z
        V[i] = static_cast<double>(nAtoms);   //  N
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);

    std::vector<double> b(V); // N
    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}



std::vector<double> calcANZ(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A

    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());


    //std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        const auto* atom = mol.getAtomWithIdx(i);
        adjMatVec[i * nAtoms + i] = static_cast<double>(nAtoms); // N
        V[i] = static_cast<double>(atom->getAtomicNum());   // Z
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);

    std::vector<double> b(V); // Z
    bool success = false;


    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}




std::vector<double> calcANI(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A

    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());


    //std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        adjMatVec[i * nAtoms + i] = static_cast<double>(nAtoms); // N
        V[i] = static_cast<double>(1);   // I
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);

    std::vector<double> b(V); // I
    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}




std::vector<double> calcDSZ(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0); // "D"
    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"
    std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

    std::vector<int> DistSum (nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<int>(distances[i * nAtoms + j]); // "S"
        }
    }
    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        const auto* atom = mol.getAtomWithIdx(i);
        DistMatVec[i * nAtoms + i] = static_cast<double>(DistSum[i]); // Distance Sum == S
        V[i] = static_cast<double>(atom->getAtomicNum());   // Z
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(DistMatVec);

    std::vector<double> b(V);
    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}



std::vector<double> calcANN(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0); // A

    double* Mat = RDKit::MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
    std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());


    //std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
        adjMatVec[i * nAtoms + i] = static_cast<double>(nAtoms); // N
        V[i] = static_cast<double>(nAtoms);   // N
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(adjMatVec);

    std::vector<double> b(V); // I
    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}



std::vector<double> calcDN2Z(const RDKit::ROMol &mol) {

    int nAtoms = rdcast<int>(mol.getNumAtoms());

    // Use RDKit's built-in function to get the adjacency matrix
    std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0); // "D"
    double* distances = MolOps::getDistanceMat(mol, false, false, false); // no need for "Bond order"
    std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

    std::vector<int> DistSum (nAtoms, 0.);

    for (int i = 0; i < nAtoms; ++i) {
        for (int j = 0; j < nAtoms; ++j) {
            DistSum[i] += static_cast<int>(distances[i * nAtoms + j]); // "S"
        }
    }
    // Modify diagonal to include atomic numbers
    std::vector<double> V(nAtoms, 0.0);

    for (int i = 0; i < nAtoms; ++i) {
      const auto* atom = mol.getAtomWithIdx(i);
        DistMatVec[i * nAtoms + i] = static_cast<double>(nAtoms*nAtoms); // Distance Sum == S
        V[i] = static_cast<double>(atom->getAtomicNum());   // Z
    }

    // Copy the adjacency matrix into a working buffer
    std::vector<double> A_flat(DistMatVec);

    std::vector<double> b(V);
    bool success = false;

    solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

    if (!success) {return {0.,0.,0.,0.,0.};}

    return TIn(mol,b,1);
}


    static const std::vector<std::string> frags = {
        "[CX4H3]","[CX4H2]","[CX4H1]","[CX4H0]","*=[CX3H2]","[$(*=[CX3H1]),$([cX3H1](a)a)]","[$(*=[CX3H0]),$([cX3H0](a)(a)A)]","c(a)(a)a","*#C","[C][NX3;H2]","[c][NX3;H2]","[C][NX3;H1][C]",
        "[c][NX3;H1]","[c][nX3;H1][c]","[C][NX3;H0](C)[C]","[c][NX3;H0](C)[C]","[c][nX3;H0][c]","*=[Nv3;!R]","*=[Nv3;R]",
        "[nX2H0,nX3H1+](a)a","N#C[A;!#1]","N#C[a;!#1]","[$([A;!#1][NX3](=O)=O),$([A;!#1][NX3+](=O)[O-])]",
        "[$([a;!#1][NX3](=O)=O),$([a;!#1][NX3+](=O)[O-])]","[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]","[OH]","[OX2;H0;!R]","[OX2;H0;R]","[oX2](a)a",
        "*=O","[SX2](*)*","[sX2](a)a","*=[SX1]","*=[SX3]","[$([#16X4](=[OX1])(=[OX1])([!#8])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([!#8])[OX2H0])]","[S,s]",
        "[P,p]","FA","Fa","Cl","Br","I","[CX3;!R](=[OX1])[OX2H0]","[CX3;R](=[OX1])[OX2H0;R]","P(=[OX1])(O)(O)O","[CX3](=[OX1])([OX2H0])[OX2H0]","[CX3](=O)[OX1H0-,OX2H1]",
        "nC=[OX1]","[N;!R]C=[OX1]","[N;R][C;R]=[OX1]","[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]","NC(=[OX1])N","[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]","[CX3](=[OX1])[NX3][CX3](=[OX1])",
        "C1(=[OX1])C=CC(=[OX1])C=C1","[$([CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])]",
        "[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]","*1~*2~*(~*3~*(~*~*~*~*3)~*1)~*~*~*1~*2~*~*~*1",
        "[OX2;H1]CC[O,N]","[OX2;H1]C[C,N]=[O,S]","[OX2;H1]c1ccccc1[O,NX3]","[OX2;H1]c1ccccc1C=[O,S]","[OX2;H1]c1ccccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]","[NH,NH2,NH3+]CC[O,N]","[NH,NH2,NH3+]c1ccccc1[O,N]",
        "[NH,NH2,NH3+]c1ccccc1[C,N]=[O,S]","[OX2H]c1ccccc1[Cl,Br,I]","[CX4]([OH])[CX4][OH]","n:n","o:n","n:c:n","o:c:n","n:c:c:n","[F,Cl,Br,I,N,O,S]-c:c-[F,Cl,Br,I,N,O,S]","[F,Cl,Br,I,N,O,S]-c:c:c-[F,Cl,Br,I,N,O,S]",
        "[F,Cl,Br,I,N,O,S]-c:c:c:c-[F,Cl,Br,I,N,O,S]","P(=[OX1])N","Nc:n","[$(cC[OH]);!$(c[CX3](=O)[OX1H0-,OX2H1])]","[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]",
        "[OX2]-c:c-[OX2]","[C][OX2H]","[c][OX2H]","[C][NX3;H1;!R][C]","[C][NX3;H1;R][C]","[c][NX3;H1;!$(NC=O)][C]","[CX3](=[OX1])[NX3H2]","[CX3](=[OX1])[NX3;H1][C]","[CX3](=[OX1])[NX3;H1][c]",
        "[$([SX4](=[OX1])(=[OX1])([!O])[NH,NH2,NH3+]),$([SX4+2]([OX1-])([OX1-])([!O])[NH,NH2,NH3+])]","[NX3;H1]C(=[OX1])[NX3;H1]","[NX3;H0]C(=[OX1])[NX3;H1]","[NX3;H1]C(=[OX1])O","[NX3;H1]C(=N)[NX3;H0]",
        "[C]#[CH]","P[OH,O-]","[CH][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]",
        "[CH]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]","[CX4]([CX3](=O)[OX1H0-,OX2H1])[CX4][CX3](=O)[OX1H0-,OX2H1]"
        "[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX3](=O)[OX1H0-,OX2H1]","[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[OH]",
        "[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][OH]","[nX3;H1]:n","[nX3;H1]:c:n","[OX1]=[C,c]~[C,c]C[OH]","[OH]c1cccc2cccnc12",
        "[OH]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1","[OH]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1",
        "[NH,NH2,NH3+]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1","[NH,NH2,NH3+]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1",
        "[CX3](=O)([OX1H0-,OX2H1])c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1","[CX3](=O)([OX1H0-,OX2H1])c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1",
        "[OH]c1c([CX4])cccc1[CX4]","[NH,NH2,NH3+]c1c([CX4])cccc1[CX4]","[OH]c1c(C[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cccc1","[OH]c1cc([CX3](=O)[OX1H0-,OX2H1])ccc1",
        "[OH]c1ccc([CX3](=O)[OX1H0-,OX2H1])cc1","[OH]c1cc([$([CH](=O)),$(C(=O)C)])ccc1","[OH]c1ccc([$([CH](=O)),$(C(=O)C)])cc1","[OX2H0+0]-[cX3H0;$(*-A)]:[cX3H0;$(*-a)]-[cX3H0;$(*-a)]:[cX3H1]:[cX3H1]",
        "[CX4H3]-[nX3H0+0;$(*-A)]","[FX1H0]-[CX4H1](-[nX3H0+0;$(*-A)]1:[cX3H0;$(*-A)](-[CX4H3]):[nX2H0+0;!$(*-a);!$(*~A)]:[nX3H0+0;$(*-a)](:[cX3H0;$(*=A)]:1=[OX1H0+0])-[cX3H0;$(*-a)]1:[cX3H1]:[cX3H0;$(*-A)](:[cX3H0;$(*-A)]:[cX3H1]:[cX3H0;$(*-A)]:1-[ClX1H0])-[NX3H1+0]-[SX4H0](=[OX1H0+0])(=[OX1H0+0])-[CX4H3])-[FX1H0]",
        "[cX3H1]:[cX3H1]:[cX3H1]:[cX3H0;$(*-A)]-[CX3H1]=[CX3H1]-[CX2H0]#[NX1H0+0]","[OX2H1+0]-[cX3H0;$(*-A)](:[cX3H1]):[cX3H1]","[CX4H3]-[cX3H0;$(*-A)]:[cX3H0;$(*-A)]","[CX4H3]-[cX3H0;$(*-A)]:[cX3H0;!$(*-a);!$(*~A)]:[cX3H0;!$(*-a);!$(*~A)]",
        "[CX4H3]-[OX2H0+0]-[cX3H0;$(*-A)]","[CX4H2]-[SX2H0]","[CX4H2]-[CX4H3]","[cX3H1]:[cX3H1]:[cX3H0;!$(*-a);!$(*~A)](:[cX3H0;!$(*-a);!$(*~A)]):[cX3H0;!$(*-a);!$(*~A)]","[cX3H0;!$(*-a);!$(*~A)]:[cX3H0;!$(*-a);!$(*~A)]:[cX3H1]:[cX3H0;!$(*-a);!$(*~A)]:[cX3H0;!$(*-a);!$(*~A)]:[cX3H1]:[cX3H1]",
        "[nX2H0+0;!$(*-a);!$(*~A)]:[nX3H0+0;$(*-A)]","[cX3H0;$(*-A)]:[cX3H0;$(*-A)]:[nX2H0+0;!$(*-a);!$(*~A)]","[CX4H3]-[SX2H0]","[CX4H3]-[NX3H0+0]-[CX4H3]","[CX4H2]-[CX4H2]-[cX3H0;$(*-A)]",
        "[cX3H1]:[cX3H0;$(*-A)]-[CX2H0]#[NX1H0+0]","[OX2H0+0]-[cX3H0;$(*-A)]:[cX3H0;$(*-A)]:[cX3H1]:[cX3H1]","[CX4H1]-[CX4H2]","[CX4H1]-[cX3H0;$(*-A)]","[cX3H1]:[cX3H0;$(*-A)]:[nX2H0+0;!$(*-a);!$(*~A)]","[cX3H1]:[cX3H1]:[cX3H0;$(*-A)]-[CX3H0]=[OX1H0+0]",
        "[NX3H1+0]-[CX3H0]=[OX1H0+0]","[FX1H0]-[CX4H0](-[FX1H0])-[FX1H0]","[cX3H1]:[cX3H0;$(*-a)](:[cX3H1])-[cX3H0;$(*-a)](:[cX3H1]):[cX3H1]","[NX3H1+0]-[cX3H0;$(*-A)](:[cX3H1]):[cX3H1]","[cX3H1]:[cX3H1]:[cX3H0;$(*-A)]-[NX3H0+0](=[OX1H0+0])=[OX1H0+0]",
        "[CX4H3]-[cX3H0;$(*-A)]:[cX3H1]:[cX3H1]:[cX3H1]","[CX4H3]-[CX4H0]-[CX4H3]","[CX4H1]-[CX4H1]","[CX4H2]-[CX3H0](=[OX1H0+0])-[NX3H0+0]","[NX3H0+0]-[cX3H0;$(*-A)](:[cX3H0;$(*-A)]):[cX3H0;$(*-A)]","[OX2H1+0]-[cX3H0;$(*-A)](:[cX3H1]):[cX3H0;$(*-A)]-[ClX1H0]","[CX4H2]-[OX2H0+0]-[CX3H0]=[OX1H0+0]",
        "[CX4H3]-[OX2H0+0]-[PX4H0](=[SX1H0])-[OX2H0+0]-[CX4H3]","[OX2H0+0]-[CX3H0](=[OX1H0+0])-[cX3H0;$(*-A)]:[cX3H0;$(*-A)]","[NX3H0+0]-[cX3H0;$(*-A)]:[cX3H0;$(*-A)]:[cX3H1]:[cX3H0;$(*-A)]","[cX3H0;$(*-A)]:[cX3H1]:[cX3H0;$(*-A)]:[cX3H1]:[cX3H0;$(*-A)]",
        "[CX4H2]-[cX3H0;$(*-A)]:[cX3H1]:[cX3H0;$(*-A)]","[CX4H0]-[CX4H0]-[ClX1H0]","[cX3H0;$(*-A)]1:[cX3H0;$(*-A)]:[cX3H0;$(*-A)]:[cX3H0;$(*-A)]:[cX3H0;$(*-A)]:[cX3H0;$(*-A)]:1","[A!#1x0+0]#[A!#1x0+0]","[SX2H0]","[ax3+0;$(*-[A!#1])]","[NX3H0+0]","[CX3H1]=[CX3H2]","[CX3H1]=[OX1H0+0]","[cX3H0;$(*-A)]",
        "[ax3+0;$(*-a)]","[nX3H1+0]","[OX2H0+0]","[A!#1x0+0]","[#8]","[cX3H0;!$(*-a);!$(*~A)]","[#7]","[#6]","[SX2H1]","[CX3](=O)[OX2H1]","[$([CX3H][#6]),$([CX3H2])]=[OX1]","[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]",
        "[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]","[CX4](F)(F)F","[NX3]=[CX3]","[NX3][CX3]=[NX3]","[NX1]#[CX2]","[CX3]=[OX1]","[#6][CX3](=O)[#6]","[CX3H1](=O)[#6]","[NX3][CX3](=[OX1])[#6]",
        "[NX3][CX3](=[OX1])[#5]","[NX3][CX3](=[OX1])[OX2H0]","[NX3,NX4+][CX3](=[OX1])[OX2H,OX1-]","[#6][CX3](=[OX1])[OX2H0][#6]","[CX3](=[OX1])[OX1-]","[CX3](=O)[OX2H1]","[OX1]=[CX3]([OX2])[OX2]","[CX3]=[SX1]","[NX3][NX3]","[NX2]=N",
        "[NX2]=[OX1]","[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]","[OX1]=[NX2][OX2]","[OX2,OX1-][OX2,OX1-]","[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]","[$([#16X4](=[OX1])(=[OX1])([#6])[#6]),$([#16X4+2]([OX1-])([OX1-])([#6])[#6])]",
        "[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]","[$([#16X3](=[OX1])[OX2H0]),$([#16X3+]([OX1-])[OX2H0])]","[$([#16X3](=[OX1])[OX2H,OX1H0-]),$([#16X3+]([OX1-])[OX2H,OX1H0-])]",
        "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]","[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]",
        "[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]","[$([SX4](=O)(=O)(O)O),$([SX4+2]([O-])([O-])(O)O)]","[#16X2H0][#16X2H0]","[PX5](=[OX1])([OX1-])[OX1-]",
        "[PX5](=[OX1])([OX2H])[OX2H]","[PX6](=[OX1])([OX1-])([OX1-])[OX1-]"};


     static const std::vector<std::shared_ptr<RDKit::RWMol>> GetQueriesFrags() {
      static const std::vector<std::shared_ptr<RDKit::RWMol>> queriesFrags =
          [] {
            std::vector<std::shared_ptr<RDKit::RWMol>> res;
            for (const auto &smi : frags) {
              auto mol = RDKit::SmartsToMol(smi);
              if (mol) {
                res.emplace_back(std::move(mol));
              } else {
                std::cerr << "Invalid SMARTS: " << smi << std::endl;
              }
            }
            return res;
          }();
       return queriesFrags;
    }


    std::vector<double> calcFrags(const RDKit::ROMol& mol) {
        auto queriesFrags = GetQueriesFrags();
        std::vector<double> retval(queriesFrags.size(), 0.0);

        try {
            // Calculate A descriptor
            for (size_t i = 0; i < queriesFrags.size(); ++i) {

                std::vector<RDKit::MatchVectType> matches;
                RDKit::SubstructMatch(mol, *queriesFrags[i], matches, true);  // uniquify = true
                retval[i] = matches.size();
            }


        } catch (const std::exception& e) {
            std::cerr << "Error in SMARTS matching: " << e.what() << std::endl;
            throw std::runtime_error("Error in SMARTSQueryTool");
        }

        return retval;
    }
#else
bool hasOsmordredSupport() {
    return false;
}

std::vector<double> calcABCIndex(const ROMol& mol) { return std::vector<double>(); }
std::vector<int> calcAcidBase(const ROMol& mol) { return std::vector<int>(); }
std::vector<int> calcAromatic(const ROMol& mol) { return std::vector<int>(); }
std::vector<int> calcAtomCounts(const ROMol& mol) { return std::vector<int>(); }
std::vector<double> calcBalabanJ(const ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcBertzCT(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<int> calcBondCounts(const RDKit::ROMol& mol) { return std::vector<int>(); }
std::vector<double> calcVertexAdjacencyInformation(const ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcWeight(const ROMol& mol) { return std::vector<double>(); }
std::vector<int> calcWienerIndex(const ROMol& mol) { return std::vector<int>(); }
std::vector<double> calcVdwVolumeABC(const ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcTopoPSA(const ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcSLogP(const ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcHydrogenBond(const ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcLogS(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<int> calcLipinskiGhose(const RDKit::ROMol& mol) { return std::vector<int>(); }
std::vector<double> calcMcGowanVolume(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcPolarizability(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcRotatableBond(const ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcFragmentComplexity(const ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcConstitutional(const ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcTopologicalIndex(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcDetourMatrixDescs(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcDetourMatrixDescsL(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcDistMatrixDescs(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcDistMatrixDescsL(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcAdjMatrixDescs(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcAdjMatrixDescsL(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcCarbonTypes(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcEccentricConnectivityIndex(const ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcBaryszMatrixDescsL(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcBaryszMatrixDescs(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcZagrebIndex(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcMoeType(const ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcMolecularDistanceEdgeDescs(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcEStateDescs(const RDKit::ROMol& mol, bool extended) { return std::vector<double>(); }
std::vector<double> calcWalkCounts(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcTopologicalChargeDescs(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcAllChiDescriptors(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcPathCount(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcKappaShapeIndex(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<int> calcRingCount(const ROMol& mol) { return std::vector<int>(); }
std::vector<double> calcMolecularId(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcBCUTs(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcAutoCorrelation(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcFramework(const ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcExtendedTopochemicalAtom(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calculateETADescriptors(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcChipath(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcChichain(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcChicluster(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcChipathcluster(const RDKit::ROMol& mol) { return std::vector<double>(); }
int calcAcidicGroupCount(const ROMol& mol) { return 0; }
int calcBasicGroupCount(const ROMol& mol) { return 0; }
int countAromaticAtoms(const ROMol& mol) { return 0; }
int countAromaticBonds(const ROMol& mol) { return 0; }
std::vector<double> calcBEStateDescs(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcHEStateDescs(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcAlphaKappaShapeIndex(const RDKit::ROMol& mol) { return std::vector<double>(); } // closer "missing" k3 path count not correct on few cases
std::vector<double> calcAbrahams(const RDKit::ROMol& mol) { return std::vector<double>(); } // Platts, Butina, Abraham, Hersey  paper J Chem Inf Comput Sci. 1999 30/8/01;39(5):835-45
std::vector<double> calcPol(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcMR(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcFlexibility(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcODT(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcSchultz(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcRNCG_RPCG(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcAZV(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcASV(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcDSV(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcAZS(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcASZ(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcDN2S(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcDN2I(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcASI(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcDSI(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcASN(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcDSN(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcDN2N(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcANS(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcANV(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcAZN(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcANZ(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcANI(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcDSZ(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcANN(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcDN2Z(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcANMat(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcAZMat(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcASMat(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcDSMat(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcDN2Mat(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcFrags(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcAddFeatures(const RDKit::ROMol& mol) { return std::vector<double>(); }
std::vector<double> calcInformationContent(const RDKit::ROMol& mol, int maxradius) { return std::vector<double>(); } // Inspired by 1984 Basak paper

// Aggregated fast path that calls all Osmordred descriptors in C++
std::vector<double> calcOsmordred(const RDKit::ROMol& mol) { return std::vector<double>(); }
  
#endif
} // end namespae descriptors
} // end namespace osmordred
} // end rdkit namespace
