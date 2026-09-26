// RDKit 217 Descriptors Implementation
// Extracts all 217 RDKit descriptors in exact order matching Python's Descriptors._descList
// This is a comprehensive implementation for cascade models

#include "RDKit217Descriptors.h"
#include "../MolDescriptors.h"
#include "../Crippen.h"
#include "../MolSurf.h"
#include "../Lipinski.h"
#include "../ConnectivityDescriptors.h"
#include "../BCUT.h"
#include "../Property.h"
#include "../Osmordred.h"

#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <DataStructs/SparseIntVect.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDThreads.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <string>
#include <set>
#include <map>
#include <utility>
#include <memory>  // For std::shared_ptr
#include <future>  // For std::async
#include <array>
#include <cstdio>
#include <stdexcept>
#include <ML/InfoTheory/InfoGainFuncs.h>  // RDInfoTheory::InfoEntropy (rdInfoTheory)

// Python parity: this file's own arithmetic ports Python/numpy expressions,
// which are evaluated one IEEE operation at a time. Disable FMA contraction
// (clang's default -ffp-contract=on) for everything below; code in the headers
// above, like the C++ that Python calls, keeps the default.
#if defined(__clang__)
#pragma clang fp contract(off)
#endif

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// Defined (and exported) in OsmordredBasicPhyschemCountsRules.cpp but not declared in
// Osmordred.h; declared here so this file does not modify the v3 sources.
std::vector<double> calcEState_VSA(const ROMol &mol);
std::vector<double> calcVSA_EState(const ROMol &mol);

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

// ---------------------------------------------------------------------------
// Python-arithmetic helpers.
//
// Everything below the file-level "fp contract(off)" pragma (see top of file)
// is evaluated one IEEE operation at a time, the way CPython and numpy do:
// clang's default -ffp-contract=on would otherwise fuse a*b+c into an FMA and
// change last bits. C++ code that Python itself calls (rdMolDescriptors,
// rdInfoTheory.InfoEntropy, ...) is called here too, compiled with the same
// default flags as the Python extension modules.
// ---------------------------------------------------------------------------

// math.log: ValueError for x <= 0 (nan -> nan, inf -> inf).
double pyLog(double x) {
    if (x <= 0.0) {
        throw std::domain_error("math domain error");
    }
    return std::log(x);
}

// math.exp: OverflowError when a finite argument overflows.
double pyExp(double x) {
    const double r = std::exp(x);
    if (std::isinf(r) && std::isfinite(x)) {
        throw std::range_error("math range error");
    }
    return r;
}

// numpy's pairwise summation (pairwise_sum_DOUBLE, numpy 1.26), used by
// add.reduce and therefore by numpy.trace. Verified bit-for-bit against
// numpy.trace on random diagonals.
double numpyPairwiseSum(const double* a, std::size_t n, std::size_t stride) {
    if (n < 8) {
        double res = 0.;
        for (std::size_t i = 0; i < n; ++i) {
            res += a[i * stride];
        }
        return res;
    }
    if (n <= 128) {
        double r[8];
        for (std::size_t j = 0; j < 8; ++j) {
            r[j] = a[j * stride];
        }
        std::size_t i = 8;
        for (; i < n - (n % 8); i += 8) {
            for (std::size_t j = 0; j < 8; ++j) {
                r[j] += a[(i + j) * stride];
            }
        }
        double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; ++i) {
            res += a[i * stride];
        }
        return res;
    }
    std::size_t n2 = n / 2;
    n2 -= n2 % 8;
    return numpyPairwiseSum(a, n2, stride) + numpyPairwiseSum(a + n2 * stride, n - n2, stride);
}

// rdkit.ML.InfoTheory.entropy.InfoEntropy is the C++ RDInfoTheory::InfoEntropy
// (rdInfoTheory module) on a float64 array. Call the same template.
double infoEntropy(std::vector<double> v) {
    return RDInfoTheory::InfoEntropy(v.data(), static_cast<long int>(v.size()));
}

// numpy.dot(A, B) for an n x n 0/1 matrix A and float64 B, as computed by the
// OpenBLAS dgemm numpy links (openblas64 0.3.23): each C[i][j] accumulates
// fma(A[i][k], B[k][j], acc) sequentially over k inside a K block, and the
// blocks are added to C in order. K blocks follow OpenBLAS level3.c with
// GEMM_Q = 128: take 128 while >= 256 remain, split the last (128, 256) in two
// halves (first = ceil(r/2)). Verified bit-for-bit against numpy.dot for
// n = 2..700. With A in {0, 1}, fma(a, b, acc) == acc + a*b exactly; a zero
// entry is skipped unless the B row holds inf/nan (0*inf = nan in the kernel).
void openblasDotAdjacency(const std::vector<double>& A, const std::vector<double>& B,
                          std::vector<double>& C, unsigned int n) {
    std::vector<char> rowFinite(n, 1);
    for (unsigned int k = 0; k < n; ++k) {
        for (unsigned int j = 0; j < n; ++j) {
            if (!std::isfinite(B[k * n + j])) {
                rowFinite[k] = 0;
                break;
            }
        }
    }
    std::vector<double> acc(n);
    unsigned int k0 = 0;
    bool first = true;
    while (k0 < n) {
        unsigned int kl = n - k0;
        if (kl >= 256) {
            kl = 128;
        } else if (kl > 128) {
            kl = (kl + 1) / 2;
        }
        for (unsigned int i = 0; i < n; ++i) {
            std::fill(acc.begin(), acc.end(), 0.0);
            for (unsigned int k = k0; k < k0 + kl; ++k) {
                const double a = A[i * n + k];
                if (a == 0.0 && rowFinite[k]) {
                    continue;
                }
                const double* b = &B[k * n];
                for (unsigned int j = 0; j < n; ++j) {
                    acc[j] = acc[j] + a * b[j];
                }
            }
            double* c = &C[i * n];
            if (first) {
                for (unsigned int j = 0; j < n; ++j) {
                    c[j] = acc[j];
                }
            } else {
                for (unsigned int j = 0; j < n; ++j) {
                    c[j] = c[j] + acc[j];
                }
            }
        }
        first = false;
        k0 += kl;
    }
}

// abs(Graphs.CharacteristicPolynomial(mol, adjMat)) with
// adjMat = numpy.equal(Chem.GetDistanceMatrix(mol, False), 1)  (GraphDescriptors.Ipc).
// Le Verrier-Faddeev-Frame, literally:
//   I = identity; An = A; res[0] = 1
//   for n in 1..N: res[n] = 1./n * trace(An); Bn = An - res[n]*I; An = dot(A, Bn)
//   res[1:] *= -1
std::vector<double> ipcAbsCharPoly(const ROMol& mol) {
    const unsigned int n = mol.getNumAtoms();
    const double* dMat = MolOps::getDistanceMat(mol, false, false, false);
    std::vector<double> A(static_cast<std::size_t>(n) * n);
    for (std::size_t i = 0; i < A.size(); ++i) {
        A[i] = (dMat[i] == 1.0) ? 1.0 : 0.0;
    }
    std::vector<double> res(n + 1, 0.0);
    res[0] = 1.0;
    std::vector<double> An = A;
    std::vector<double> Bn(A.size());
    for (unsigned int m = 1; m <= n; ++m) {
        res[m] = 1. / m * numpyPairwiseSum(An.data(), n, n + 1);
        for (unsigned int i = 0; i < n; ++i) {
            for (unsigned int j = 0; j < n; ++j) {
                Bn[i * n + j] = An[i * n + j] - res[m] * (i == j ? 1.0 : 0.0);
            }
        }
        openblasDotAdjacency(A, Bn, An, n);
    }
    for (unsigned int m = 1; m <= n; ++m) {
        res[m] *= -1;
    }
    for (auto& v : res) {
        v = std::fabs(v);
    }
    return res;
}

// GraphDescriptors.BalabanJ, literally. The distance matrix is
// GetDistanceMatrix(useBO=True, prefix="Balaban"); the adjacency test is done
// on mol._adjMat, which AvgIpc (called first by CalcMolDescriptors) set to the
// topological distance matrix: "adjMat[i, j] == 1" is then "bonded".
// s = sum(dMat) is Python's builtin sum over rows (unreachable pairs count as
// 1e8, as in Python).
double pyBalabanJ(const ROMol& mol) {
    const unsigned int n = mol.getNumAtoms();
    const double* dMat = MolOps::getDistanceMat(mol, true, false, false, "Balaban");
    const double* adjMat = MolOps::getDistanceMat(mol, false, false, false);
    std::vector<double> s(n, 0.0);
    for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = 0; j < n; ++j) {
            s[j] = s[j] + dMat[i * n + j];
        }
    }
    const int q = static_cast<int>(mol.getNumBonds());
    const int mu = q - static_cast<int>(n) + 1;
    double sum_ = 0.;
    for (unsigned int i = 0; i < n; ++i) {
        const double si = s[i];
        for (unsigned int j = i; j < n; ++j) {
            if (adjMat[i * n + j] == 1) {
                sum_ += 1. / std::sqrt(si * s[j]);
            }
        }
    }
    if (mu + 1 != 0) {
        return static_cast<double>(q) / static_cast<double>(mu + 1) * sum_;
    }
    return 0.0;
}

// GraphDescriptors.BertzCT(mol, cutoff=100, forceDMat=1), literally, including
// dict insertion order (which fixes the summation order of the entropies).
double pyBertzCT(const ROMol& mol) {
    const int cutoff = 100;
    const unsigned int numAtoms = mol.getNumAtoms();
    if (numAtoms < 2) {
        return 0.0;
    }
    // _CreateBondDictEtc
    std::map<std::pair<unsigned int, unsigned int>, double> bondDict;
    std::vector<std::vector<unsigned int>> nList(numAtoms);
    for (const auto bond : mol.bonds()) {
        unsigned int atom1 = bond->getBeginAtomIdx();
        unsigned int atom2 = bond->getEndAtomIdx();
        if (atom1 > atom2) {
            std::swap(atom1, atom2);
        }
        bondDict[{atom1, atom2}] = bond->getIsAromatic() ? 1.5 : bond->getBondTypeAsDouble();
        if (std::find(nList[atom1].begin(), nList[atom1].end(), atom2) == nList[atom1].end()) {
            nList[atom1].push_back(atom2);
        }
        if (std::find(nList[atom2].begin(), nList[atom2].end(), atom1) == nList[atom2].end()) {
            nList[atom2].push_back(atom1);
        }
    }
    for (auto& element : nList) {
        std::sort(element.begin(), element.end());
    }
    auto lookUpBondOrder = [&bondDict](unsigned int a1, unsigned int a2) {
        return a1 < a2 ? bondDict.at({a1, a2}) : bondDict.at({a2, a1});
    };

    // _AssignSymmetryClasses: rows of the (forced) Balaban BO distance matrix,
    // sorted, cut at `cutoff` columns, keyed by '%.4f,' * ncols.
    const double* bdMat = MolOps::getDistanceMat(mol, true, false, true, "Balaban");
    const unsigned int nCols = std::min<unsigned int>(numAtoms, cutoff);
    std::map<std::string, int> keysSeen;
    std::vector<int> symmetryClasses(numAtoms);
    std::vector<double> row(numAtoms);
    char buf[64];
    for (unsigned int i = 0; i < numAtoms; ++i) {
        std::copy(bdMat + i * numAtoms, bdMat + (i + 1) * numAtoms, row.begin());
        std::sort(row.begin(), row.end());
        std::string key;
        for (unsigned int j = 0; j < nCols; ++j) {
            std::snprintf(buf, sizeof(buf), "%.4f,", row[j]);
            key += buf;
        }
        const int next = static_cast<int>(keysSeen.size());
        symmetryClasses[i] = keysSeen.emplace(key, next).first->second + 1;
    }

    // atomTypeDict / connectionDict with Python dict insertion order.
    std::vector<std::pair<int, int>> atomTypes;
    std::vector<std::pair<std::array<int, 3>, double>> connections;
    std::map<std::array<int, 3>, std::size_t> connectionIdx;
    auto addConnection = [&](const std::array<int, 3>& key, double numConnections) {
        auto it = connectionIdx.find(key);
        if (it == connectionIdx.end()) {
            connectionIdx.emplace(key, connections.size());
            connections.emplace_back(key, 0 + numConnections);
        } else {
            connections[it->second].second = connections[it->second].second + numConnections;
        }
    };
    for (unsigned int atomIdx = 0; atomIdx < numAtoms; ++atomIdx) {
        const int hingeAtomNumber = mol.getAtomWithIdx(atomIdx)->getAtomicNum();
        auto at = std::find_if(atomTypes.begin(), atomTypes.end(),
                               [hingeAtomNumber](const auto& p) { return p.first == hingeAtomNumber; });
        if (at == atomTypes.end()) {
            atomTypes.emplace_back(hingeAtomNumber, 1);
        } else {
            ++at->second;
        }
        const int hingeAtomClass = symmetryClasses[atomIdx];
        const auto& neighbors = nList[atomIdx];
        const std::size_t numNeighbors = neighbors.size();
        for (std::size_t i = 0; i < numNeighbors; ++i) {
            const unsigned int neighbor_iIdx = neighbors[i];
            const int NiClass = symmetryClasses[neighbor_iIdx];
            const double bond_i_order = lookUpBondOrder(atomIdx, neighbor_iIdx);
            if (bond_i_order > 1 && neighbor_iIdx > atomIdx) {
                const double numConnections = bond_i_order * (bond_i_order - 1) / 2;
                // 2-tuple key (min, max); third slot 0 never occurs in a 3-tuple
                addConnection({std::min(hingeAtomClass, NiClass), std::max(hingeAtomClass, NiClass), 0},
                              numConnections);
            }
            for (std::size_t j = i + 1; j < numNeighbors; ++j) {
                const unsigned int neighbor_jIdx = neighbors[j];
                const int NjClass = symmetryClasses[neighbor_jIdx];
                const double bond_j_order = lookUpBondOrder(atomIdx, neighbor_jIdx);
                const double numConnections = bond_i_order * bond_j_order;
                addConnection({std::min(NiClass, NjClass), hingeAtomClass, std::max(NiClass, NjClass)},
                              numConnections);
            }
        }
    }
    std::vector<double> connectionList;
    if (connections.empty()) {
        connectionList.push_back(1.0);  // connectionDict = {'a': 1}
    } else {
        for (const auto& c : connections) {
            connectionList.push_back(c.second);
        }
    }

    // _CalculateEntropies
    double totConnections = 0;
    for (double v : connectionList) {
        totConnections = totConnections + v;
    }
    const double log2val = std::log(2.0);
    const double connectionIE =
        totConnections * (infoEntropy(connectionList) + pyLog(totConnections) / log2val);
    std::vector<double> atomTypeList;
    for (const auto& p : atomTypes) {
        atomTypeList.push_back(static_cast<double>(p.second));
    }
    const double atomTypeIE = static_cast<double>(numAtoms) * infoEntropy(atomTypeList);
    return atomTypeIE + connectionIE;
}

}  // namespace

// Helper functions for cached SMARTS queries (EXACT Osmordred pattern with IIFE)
const std::vector<std::shared_ptr<RWMol>>& GetQEDAcceptorQueries() {
    static const std::vector<std::shared_ptr<RWMol>> queries = [] {
        std::vector<std::shared_ptr<RWMol>> res;
        std::vector<std::string> patterns = {
            "[oH0;X2]", "[OH1;X2;v2]", "[OH0;X2;v2]", "[OH0;X1;v2]", "[O-;X1]",
            "[SH0;X2;v2]", "[SH0;X1;v2]", "[S-;X1]", "[nH0;X2]", "[NH0;X1;v3]",
            "[$([N;+0;X3;v3]);!$(N[C,S]=O)]"
        };
        for (const auto& pattern : patterns) {
            auto mol = RDKit::SmartsToMol(pattern);
            if (mol) {
                res.emplace_back(std::shared_ptr<RWMol>(mol));
            } else {
                res.emplace_back(nullptr);
            }
        }
        return res;
    }();  // Immediately invoked lambda expression (IIFE)
    return queries;
}

const std::vector<std::shared_ptr<RWMol>>& GetQEDAlertQueries() {
    static const std::vector<std::shared_ptr<RWMol>> queries = [] {
        std::vector<std::shared_ptr<RWMol>> res;
        std::vector<std::string> patterns = {
            "*1[O,S,N]*1", "[S,C](=[O,S])[F,Br,Cl,I]", "[CX4][Cl,Br,I]", "[#6]S(=O)(=O)O[#6]",
            "[$([CH]),$(CC)]#CC(=O)[#6]", "[$([CH]),$(CC)]#CC(=O)O[#6]", "n[OH]",
            "[$([CH]),$(CC)]#CS(=O)(=O)[#6]", "C=C(C=O)C=O", "n1c([F,Cl,Br,I])cccc1", "[CH1](=O)", "[#8][#8]",
            "[C;!R]=[N;!R]", "[N!R]=[N!R]", "[#6](=O)[#6](=O)", "[#16][#16]", "[#7][NH2]", "C(=O)N[NH2]",
            "[#6]=S", "[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]",
            "C1(=[O,N])C=CC(=[O,N])C=C1", "C1(=[O,N])C(=[O,N])C=CC=C1", "a21aa3a(aa1aaaa2)aaaa3",
            "a31a(a2a(aa1)aaaa2)aaaa3", "a1aa2a3a(a1)A=AA=A3=AA=A2", "c1cc([NH2])ccc1",
            "[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd,Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]",
            "I", "OS(=O)(=O)[O-]", "[N+](=O)[O-]", "C(=O)N[OH]", "C1NC(=O)NC(=O)1", "[SH]", "[S-]",
            "c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]", "c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]",
            "[CR1]1[CR1][CR1][CR1][CR1][CR1][CR1]1", "[CR1]1[CR1][CR1]cc[CR1][CR1]1",
            "[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2][CR2]1", "[CR2]1[CR2][CR2]cc[CR2][CR2][CR2]1",
            "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1", "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
            "C#C", "[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]", "[$([N+R]),$([n+R]),$([N+]=C)][O-]",
            "[#6]=N[OH]", "[#6]=NOC=O", "[#6](=O)[CX4,CR0X3,O][#6](=O)", "c1ccc2c(c1)ccc(=O)o2",
            "[O+,o+,S+,s+]", "N=C=O", "[NX3,NX4][F,Cl,Br,I]", "c1ccccc1OC(=O)[#6]", "[CR0]=[CR0][CR0]=[CR0]",
            "[C+,c+,C-,c-]", "N=[N+]=[N-]", "C12C(NC(N1)=O)CSC2", "c1c([OH])c([OH,NH2,NH])ccc1", "P",
            "[N,O,S]C#N", "C=C=O", "[Si][F,Cl,Br,I]", "[SX2]O", "[SiR0,CR0](c1ccccc1)(c2ccccc2)(c3ccccc3)",
            "O1CCCCC1OC2CCC3CCCCC3C2", "N=[CR0][N,n,O,S]",
            "[cR2]1[cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2][cR2]1[cR2]2[cR2][cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2]2",
            "C=[C!r]C#N", "[cR2]1[cR2]c([N+0X3R0,nX3R0])c([N+0X3R0,nX3R0])[cR2][cR2]1",
            "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2]c([N+0X3R0,nX3R0])[cR2]1",
            "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2][cR2]c1([N+0X3R0,nX3R0])", "[OH]c1ccc([OH,NH2,NH])cc1",
            "c1ccccc1OC(=O)O", "[SX2H0][N]", "c12ccccc1(SC(S)=N2)", "c12ccccc1(SC(=S)N2)", "c1nnnn1C=O",
            "s1c(S)nnc1NC=O", "S1C=CSC1=S", "C(=O)Onnn", "OS(=O)(=O)C(F)(F)F", "N#CC[OH]", "N#CC(=O)",
            "S(=O)(=O)C#N", "N[CH2]C#N", "C1(=O)NCC1", "S(=O)(=O)[O-,OH]", "NC[F,Cl,Br,I]", "C=[C!r]O",
            "[NX2+0]=[O+0]", "[OR0,NR0][OR0,NR0]", "C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]", "[CX2R0][NX3R0]",
            "c1ccccc1[C;!R]=[C;!R]c2ccccc2", "[NX3R0,NX4R0,OR0,SX2R0][CX4][NX3R0,NX4R0,OR0,SX2R0]",
            "[s,S,c,C,n,N,o,O]~[n+,N+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,C,n,N,o,O]",
            "[s,S,c,C,n,N,o,O]~[nX3+,NX3+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]", "[*]=[N+]=[*]", "[SX3](=O)[O-,OH]",
            "N#N", "F.F.F.F", "[R0;D2][R0;D2][R0;D2][R0;D2]", "[cR,CR]~C(=O)NC(=O)~[cR,CR]", "C=!@CC=[O,S]",
            "[#6,#8,#16][#6](=O)O[#6]", "c[C;R0](=[O,S])[#6]", "c[SX2][C;!R]", "C=C=C",
            "c1nc([F,Cl,Br,I,S])ncc1", "c1ncnc([F,Cl,Br,I,S])c1", "c1nc(c2c(n1)nc(n2)[F,Cl,Br,I])",
            "[#6]S(=O)(=O)c1ccc(cc1)F", "[15N]", "[13C]", "[18O]", "[34S]"
        };
        for (const auto& pattern : patterns) {
            auto mol = RDKit::SmartsToMol(pattern);
            if (mol) {
                res.emplace_back(std::shared_ptr<RWMol>(mol));
            } else {
                res.emplace_back(nullptr);
            }
        }
        return res;
    }();  // Immediately invoked lambda expression (IIFE)
    return queries;
}

// The 85 fr_* descriptors, in Descriptors._descList order, with the SMARTS
// copied verbatim from Data/FragmentDescriptors.csv (column 3; Fragments.py
// keeps the trailing newline, which the SMARTS parser ignores).
// Fragments.py counts len(mol.GetSubstructMatches(patt, uniquify=True)),
// i.e. uniquified matches with the default maxMatches=1000.
struct FragmentDef {
    const char* name;
    const char* smarts;
};
constexpr FragmentDef kFragmentDefs[] = {
    {"fr_Al_COO", "C-C(=O)[O;H1,-]"},
    {"fr_Al_OH", "[C!$(C=O)]-[OH]"},
    {"fr_Al_OH_noTert", "[$(C-[OX2H]);!$([CX3](-[OX2H])=[OX1]);!$([CD4]-[OX2H])]"},
    {"fr_ArN", "[$(a-[NX3H2]),$(a-[NH1][NH2]),$(a-C(=[OX1])[NH1][NH2]),$(a-C(=[NH])[NH2])]"},
    {"fr_Ar_COO", "c-C(=O)[O;H1,-]"},
    {"fr_Ar_N", "n"},
    {"fr_Ar_NH", "[nH]"},
    {"fr_Ar_OH", "c[OH1]"},
    {"fr_COO", "[#6]C(=O)[O;H,-1]"},
    {"fr_COO2", "[CX3](=O)[OX1H0-,OX2H1]"},
    {"fr_C_O", "[CX3]=[OX1]"},
    {"fr_C_O_noCOO", "[C!$(C-[OH])]=O"},
    {"fr_C_S", "C=[SX1]"},
    {"fr_HOCCN", "[$([OX2H1][CX4][CX4H2][NX3&R1]),$([OH1][CX4][CX4H2][NX3][CX4](C)(C)C)]"},
    {"fr_Imine", "[Nv3](=C)-[#6]"},
    {"fr_NH0", "[NH0,nH0]"},
    {"fr_NH1", "[NH1,nH1]"},
    {"fr_NH2", "[NH2,nH2]"},
    {"fr_N_O", "[N!$(N=O)](-O)-C"},
    {"fr_Ndealkylation1", "[$(N(-[CH3])-C-[$(C~O),$(C-a),$(C-N),$(C=C)]),$(N(-[CH2][CH3])-C-[$(C~O),$(C-a),$(C-N),$(C=C)])]"},
    {"fr_Ndealkylation2", "[$([N&R1]1(-C)CCC1),$([N&R1]1(-C)CCCC1),$([N&R1]1(-C)CCCCC1),$([N&R1]1(-C)CCCCCC1),$([N&R1]1(-C)CCCCCCC1)]"},
    {"fr_Nhpyrrole", "[nH]"},
    {"fr_SH", "[SH]"},
    {"fr_aldehyde", "[CX3H1](=O)[#6]"},
    {"fr_alkyl_carbamate", "C[NH1]C(=O)OC"},
    {"fr_alkyl_halide", "[CX4]-[Cl,Br,I,F]"},
    {"fr_allylic_oxid", "[$(C=C-C);!$(C=C-C-[N,O,S]);!$(C=C-C-C-[N,O]);!$(C12=CC(=O)CCC1C3C(C4C(CCC4)CC3)CC2)]"},
    {"fr_amide", "C(=O)-N"},
    {"fr_amidine", "C(=N)(-N)-[!#7]"},
    {"fr_aniline", "c-[NX3;!$(N=*)]"},
    {"fr_aryl_methyl", "[$(a-[CH3]),$(a-[CH2]-[CH3]),$(a-[CH2]-[CH2]~[!N;!O]);!$(a(:a!:*):a!:*)]"},
    {"fr_azide", "[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]"},
    {"fr_azo", "[#6]-N=N-[#6]"},
    {"fr_barbitur", "C1C(=O)NC(=O)NC1=O"},
    {"fr_benzene", "c1ccccc1"},
    {"fr_benzodiazepine", "[c&R2]12[c&R1][c&R1][c&R1][c&R1][c&R2]1[N&R1][C&R1][C&R1][N&R1]=[C&R1]2"},
    {"fr_bicyclic", "[R2][R2]"},
    {"fr_diazo", "[N+]#N"},
    {"fr_dihydropyridine", "[$([NX3H1]1-C=C-C-C=C1),$([Nv3]1=C-C-C=C-C1),$([Nv3]1=C-C=C-C-C1),$([NX3H1]1-C-C=C-C=C1)]"},
    {"fr_epoxide", "O1CC1"},
    {"fr_ester", "[#6][CX3](=O)[OX2H0][#6]"},
    {"fr_ether", "[OD2]([#6])[#6]"},
    {"fr_furan", "o1cccc1"},
    {"fr_guanido", "C(=N)(N)N"},
    {"fr_halogen", "[#9,#17,#35,#53]"},
    {"fr_hdrzine", "[NX3]-[NX3]"},
    {"fr_hdrzone", "C=N-[NX3]"},
    {"fr_imidazole", "n1cncc1"},
    {"fr_imide", "N(-C(=O))-C=O"},
    {"fr_isocyan", "N=C=O"},
    {"fr_isothiocyan", "N=C=S"},
    {"fr_ketone", "[#6][CX3](=O)[#6]"},
    {"fr_ketone_Topliss", "[$([CX3](=[OX1])(C)([c,C]));!$([CX3](=[OX1])([CH1]=C)[c,C])]"},
    {"fr_lactam", "N1C(=O)CC1"},
    {"fr_lactone", "[C&R1](=O)[O&R1][C&R1]"},
    {"fr_methoxy", "[OX2](-[#6])-[CH3]"},
    {"fr_morpholine", "O1CCNCC1"},
    {"fr_nitrile", "[NX1]#[CX2]"},
    {"fr_nitro", "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]"},
    {"fr_nitro_arom", "[$(c1(-[$([NX3](=O)=O),$([NX3+](=O)[O-])])ccccc1)]"},
    {"fr_nitro_arom_nonortho", "[$(c1(-[$([NX3](=O)=O),$([NX3+](=O)[O-])])ccccc1);!$(cc-!:*)]"},
    {"fr_nitroso", "[N!$(N-O)]=O"},
    {"fr_oxazole", "c1ocnc1"},
    {"fr_oxime", "[CX3]=[NX2]-[OX2]"},
    {"fr_para_hydroxylation", "[$([cH]1[cH]cc(c[cH]1)~[$([#8,$([#8]~[H,c,C])])]),$([cH]1[cH]cc(c[cH]1)~[$([#7X3,$([#7](~[H,c,C])~[H,c,C])])]),$([cH]1[cH]cc(c[cH]1)-!:[$([NX3H,$(NC(=O)[H,c,C])])])]"},
    {"fr_phenol", "[OX2H]-c1ccccc1"},
    {"fr_phenol_noOrthoHbond", "[$(c1(-[OX2H])ccccc1);!$(cc-!:[CH2]-[OX2H]);!$(cc-!:C(=O)[O;H1,-]);!$(cc-!:C(=O)-[NH2])]"},
    {"fr_phos_acid", "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]"},
    {"fr_phos_ester", "[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]"},
    {"fr_piperdine", "N1CCCCC1"},
    {"fr_piperzine", "N1CCNCC1"},
    {"fr_priamide", "C(=O)-[NH2]"},
    {"fr_prisulfonamd", "[NH2]-S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])-[#6]"},
    {"fr_pyridine", "n1ccccc1"},
    {"fr_quatN", "[$([NX4+]),$([NX4]=*)]"},
    {"fr_sulfide", "[SX2](-[#6])-C"},
    {"fr_sulfonamd", "N-S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])-[#6]"},
    {"fr_sulfone", "S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]"},
    {"fr_term_acetylene", "C#[CH]"},
    {"fr_tetrazole", "c1nnnn1"},
    {"fr_thiazole", "c1scnc1"},
    {"fr_thiocyan", "S-C#N"},
    {"fr_thiophene", "s1cccc1"},
    {"fr_unbrch_alkane", "[CR0;D2,D1][CR0;D2][CR0;D2][CR0;D2,D1]"},
    {"fr_urea", "C(=O)(-N)-N"},
};

const std::vector<std::shared_ptr<RWMol>>& GetFragmentQueries() {
    static const std::vector<std::shared_ptr<RWMol>> queries = [] {
        std::vector<std::shared_ptr<RWMol>> res;
        res.reserve(std::size(kFragmentDefs));
        for (const auto& def : kFragmentDefs) {
            res.emplace_back(RDKit::SmartsToMol(def.smarts));
        }
        return res;
    }();  // Immediately invoked lambda expression (IIFE)
    return queries;
}

// Constants for molecule size limits to prevent hanging
constexpr unsigned int MAX_HEAVY_ATOMS_RDKIT = 200;
constexpr unsigned int MAX_RINGS_RDKIT = 30;

// Extract all 217 RDKit descriptors in exact order matching Python's Descriptors._descList
std::vector<double> extractRDKitDescriptors(const ROMol& mol) {
    std::vector<double> descriptors;
    descriptors.reserve(217);
    
    // Early exit for molecules that are too large (prevents hanging/timeout)
    unsigned int nHeavyCheck = mol.getNumHeavyAtoms();
    unsigned int nRingsCheck = RDKit::Descriptors::calcNumRings(mol);
    if (nHeavyCheck > MAX_HEAVY_ATOMS_RDKIT || nRingsCheck > MAX_RINGS_RDKIT) {
        // Return vector of NaN values (217 features)
        return std::vector<double>(217, std::numeric_limits<double>::quiet_NaN());
    }
    
    // Get some precomputed values
    double MW = calcAMW(mol);
    double exactMW = calcExactMW(mol);
    double heavyMW = calcAMW(mol, true);
    unsigned int nHeavy = calcNumHeavyAtoms(mol);
    double logP, MR;
    calcCrippenDescriptors(mol, logP, MR);
    
    // Get EState indices (needed for several descriptors)
    // Implement EState indices calculation matching Python's EStateIndices function
    // Reference: Hall, Mohney and Kier. JCICS _31_ 76-81 (1991)
    std::vector<double> estateIndices;
    try {
        const PeriodicTable* tbl = PeriodicTable::getTable();
        unsigned int nAtoms = mol.getNumAtoms();
        estateIndices.resize(nAtoms, 0.0);
        
        // Step 1: Calculate initial I-state values (Is)
        std::vector<double> Is(nAtoms, 0.0);
        for (unsigned int i = 0; i < nAtoms; ++i) {
            const Atom* atom = mol.getAtomWithIdx(i);
            unsigned int d = atom->getDegree();
            if (d > 0) {
                unsigned int atNum = atom->getAtomicNum();
                int dv = tbl->getNouterElecs(atNum) - atom->getTotalNumHs();
                // Get principal quantum number N (period number)
                int N = 1;
                if (atNum <= 2) N = 1;
                else if (atNum <= 10) N = 2;
                else if (atNum <= 18) N = 3;
                else if (atNum <= 36) N = 4;
                else if (atNum <= 54) N = 5;
                else if (atNum <= 86) N = 6;
                else N = 7;
                
                Is[i] = (4.0 / (N * N) * dv + 1.0) / d;
            }
        }
        
        // Step 2: Get distance matrix (useBO=0, useAtomWts=0 as in Python)
        double* distances = MolOps::getDistanceMat(mol, false, false, false);
        
        // Step 3: Calculate accumulative contributions
        std::vector<double> accum(nAtoms, 0.0);
        for (unsigned int i = 0; i < nAtoms; ++i) {
            for (unsigned int j = i + 1; j < nAtoms; ++j) {
                double p = distances[i * nAtoms + j] + 1.0;  // p = distance + 1
                if (p < 1e6) {  // Valid distance
                    double tmp = (Is[i] - Is[j]) / (p * p);
                    accum[i] += tmp;
                    accum[j] -= tmp;
                }
            }
        }
        
        // Step 4: Combine Is and accum
        for (unsigned int i = 0; i < nAtoms; ++i) {
            estateIndices[i] = accum[i] + Is[i];
        }
        
        // NOTE: Do NOT delete[] distances - RDKit caches it in molecule properties
        // The documentation says "The caller should NOT delete this pointer"
    } catch (...) {
        estateIndices = std::vector<double>(mol.getNumAtoms(), 0.0);
    }

    // 0: MaxAbsEStateIndex
    if (!estateIndices.empty()) {
        double maxAbs = 0.0;
        for (double v : estateIndices) {
            double absV = std::abs(v);
            if (absV > maxAbs) maxAbs = absV;
        }
        descriptors.push_back(maxAbs);
    } else {
        descriptors.push_back(0.0);
    }
    
    // 1: MaxEStateIndex
    if (!estateIndices.empty()) {
        descriptors.push_back(*std::max_element(estateIndices.begin(), estateIndices.end()));
    } else {
        descriptors.push_back(0.0);
    }
    
    // 2: MinAbsEStateIndex
    if (!estateIndices.empty()) {
        double minAbs = std::numeric_limits<double>::max();
        for (double v : estateIndices) {
            double absV = std::abs(v);
            if (absV < minAbs) minAbs = absV;
        }
        descriptors.push_back(minAbs == std::numeric_limits<double>::max() ? 0.0 : minAbs);
    } else {
        descriptors.push_back(0.0);
    }
    
    // 3: MinEStateIndex
    if (!estateIndices.empty()) {
        descriptors.push_back(*std::min_element(estateIndices.begin(), estateIndices.end()));
    } else {
        descriptors.push_back(0.0);
    }
    
    // 4: qed - QED descriptor (Quantitative Estimation of Drug-likeness)
    // Python: qed(mol, w=WEIGHT_MEAN) - uses ADS transformation on properties
    double qed_value = 0.0;
    try {
        // Remove hydrogens like Python does
        std::unique_ptr<RDKit::RWMol> molNoH_ptr(new RDKit::RWMol(mol));
        RDKit::RWMol& molNoH = *molNoH_ptr;
        RDKit::MolOps::removeHs(molNoH);
        
        // Calculate properties (matching Python's properties() function)
        double MW_qed = calcExactMW(molNoH);
        double logP_qed, MR_qed;
        calcCrippenDescriptors(molNoH, logP_qed, MR_qed);
        double ALOGP = logP_qed;
        
        // HBA: count acceptors using QED acceptor queries
        unsigned int HBA = 0;
        auto& acceptorQueries = GetQEDAcceptorQueries();
        for (const auto& query : acceptorQueries) {
            if (query) {
                std::vector<RDKit::MatchVectType> matches;
                RDKit::SubstructMatch(molNoH, *query, matches);
                HBA += matches.size();
            }
        }
        
        // HBD: hydrogen bond donors
        unsigned int HBD = calcNumHBD(molNoH);
        
        // PSA: topological polar surface area
        double PSA = calcTPSA(molNoH);
        
        // ROTB: rotatable bonds (strict mode like Python)
        unsigned int ROTB = calcNumRotatableBonds(molNoH);
        
        // AROM: aromatic rings (matching Python: len(Chem.GetSSSR(Chem.DeleteSubstructs(Chem.Mol(mol), AliphaticRings))))
        unsigned int AROM = 0;
        try {
            auto aliphaticRings = RDKit::SmartsToMol("[$([A;R][!a])]");
            if (aliphaticRings) {
                RDKit::ROMol* molDeleted = RDKit::deleteSubstructs(molNoH, *aliphaticRings);
                if (molDeleted) {
                    std::vector<std::vector<int>> sssr;
                    RDKit::MolOps::findSSSR(*molDeleted, sssr);
                    AROM = sssr.size();
                    delete molDeleted;
                }
                delete aliphaticRings;
            } else {
                std::vector<std::vector<int>> sssr;
                RDKit::MolOps::findSSSR(molNoH, sssr);
                AROM = sssr.size();
            }
        } catch (...) {
            AROM = 0;
        }
        
        // ALERTS: structural alerts count
        unsigned int ALERTS = 0;
        auto& alertQueries = GetQEDAlertQueries();
        for (const auto& query : alertQueries) {
            if (query) {
                std::vector<RDKit::MatchVectType> matches;
                RDKit::SubstructMatch(molNoH, *query, matches);
                if (!matches.empty()) {
                    ALERTS++;
                }
            }
        }
        
        // ADS parameters (from Python QED.py)
        struct ADSparam {
            double A, B, C, D, E, F, DMAX;
        };
        ADSparam adsParams[8] = {
            {2.817065973, 392.5754953, 290.7489764, 2.419764353, 49.22325677, 65.37051707, 104.9805561},  // MW
            {3.172690585, 137.8624751, 2.534937431, 4.581497897, 0.822739154, 0.576295591, 131.3186604},  // ALOGP
            {2.948620388, 160.4605972, 3.615294657, 4.435986202, 0.290141953, 1.300669958, 148.7763046},  // HBA
            {1.618662227, 1010.051101, 0.985094388, 0.000000001, 0.713820843, 0.920922555, 258.1632616},  // HBD
            {1.876861559, 125.2232657, 62.90773554, 87.83366614, 12.01999824, 28.51324732, 104.5686167},  // PSA
            {0.010000000, 272.4121427, 2.558379970, 1.565547684, 1.271567166, 2.758063707, 105.4420403},  // ROTB
            {3.217788970, 957.7374108, 2.274627939, 0.000000001, 1.317690384, 0.375760881, 312.3372610},  // AROM
            {0.010000000, 1199.094025, -0.09002883, 0.000000001, 0.185904477, 0.875193782, 417.7253140},  // ALERTS
        };
        
        // WEIGHT_MEAN (from Python)
        double weights[8] = {0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95};
        double properties[8] = {MW, ALOGP, static_cast<double>(HBA), static_cast<double>(HBD), PSA, static_cast<double>(ROTB), static_cast<double>(AROM), static_cast<double>(ALERTS)};
        
        // Apply ADS transformation: ads(x, p) = (p.A + p.B / exp1 * (1 - 1 / exp2)) / p.DMAX
        // where exp1 = 1 + exp(-(x - p.C + p.D/2) / p.E)
        //       exp2 = 1 + exp(-(x - p.C - p.D/2) / p.F)
        double sum_weighted_log = 0.0;
        double sum_weights = 0.0;
        for (int i = 0; i < 8; ++i) {
            double x = properties[i];
            const ADSparam& p = adsParams[i];
            
            double exp1 = 1.0 + std::exp(-(x - p.C + p.D / 2.0) / p.E);
            double exp2 = 1.0 + std::exp(-(x - p.C - p.D / 2.0) / p.F);
            double dx = (p.A + p.B / exp1 * (1.0 - 1.0 / exp2)) / p.DMAX;
            
            if (dx > 0.0) {
                sum_weighted_log += weights[i] * std::log(dx);
                sum_weights += weights[i];
            }
        }
        
        if (sum_weights > 0.0) {
            qed_value = std::exp(sum_weighted_log / sum_weights);
        }
    } catch (...) {
        qed_value = 0.0;
    }
    descriptors.push_back(qed_value);
    
    // 5: SPS - SPS descriptor (SpacialScore)
    // Pattern from Python SpacialScore.py: molCp = Chem.Mol(mol); rdmolops.FindPotentialStereoBonds(molCp)
    // Python does NOT sanitize - just creates a copy and calls FindPotentialStereoBonds
    double sps_value = 0.0;
    try {
        // Create a deep copy exactly like Python's Chem.Mol(mol)
        std::unique_ptr<RDKit::RWMol> molCopy_ptr(new RDKit::RWMol(mol));
        RDKit::RWMol& molCopy = *molCopy_ptr;
        
        // First: Find potential stereo bonds (like Python: rdmolops.FindPotentialStereoBonds(molCp))
        RDKit::MolOps::findPotentialStereoBonds(molCopy);
        
        // Find stereo centers - Python uses:
        //   Chem.FindMolChiralCenters(molCp, includeUnassigned=True, includeCIP=False, useLegacyImplementation=False)
        // That corresponds to tetrahedral stereo centers (not square-planar, octahedral, etc.).
        std::set<unsigned int> chiral_idxs;
        std::vector<RDKit::Chirality::StereoInfo> stereoInfo = RDKit::Chirality::findPotentialStereo(molCopy, false, true);  // cleanIt=false, flagPossible=true
        for (const auto& info : stereoInfo) {
            if (info.type == RDKit::Chirality::StereoType::Atom_Tetrahedral) {
                // Include all potential tetrahedral centers (matching includeUnassigned=True)
                chiral_idxs.insert(info.centeredOn);
            }
        }
        
        // Find stereo double bonds (E/Z) - read from bonds after findPotentialStereoBonds
        std::map<std::pair<unsigned int, unsigned int>, RDKit::Bond::BondStereo> doublebonds_stereo;
        for (auto bond : molCopy.bonds()) {
            if (bond->getBondType() == RDKit::Bond::DOUBLE) {
                auto pair = std::make_pair(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
                doublebonds_stereo[pair] = bond->getStereo();
            }
        }
        
        unsigned int nHeavy = molCopy.getNumHeavyAtoms();
        if (nHeavy == 0) {
            sps_value = 0.0;
        } else {
            double total_score = 0.0;
            
            for (auto atom : molCopy.atoms()) {
                unsigned int atom_idx = atom->getIdx();
                
                // Hybridization score (from SpacialScore.py _hybridisations dict)
                int hyb_score = 4;  // default
                RDKit::Atom::HybridizationType hyb = atom->getHybridization();
                if (hyb == RDKit::Atom::SP) hyb_score = 1;
                else if (hyb == RDKit::Atom::SP2) hyb_score = 2;
                else if (hyb == RDKit::Atom::SP3) hyb_score = 3;
                
                // Stereo score (from _accountForStereo)
                int stereo_score = 1;
                if (chiral_idxs.find(atom_idx) != chiral_idxs.end()) {
                    stereo_score = 2;
                } else {
                    // Check if atom is part of a stereo double bond
                    for (const auto& db_pair : doublebonds_stereo) {
                        if (db_pair.second != RDKit::Bond::STEREONONE) {
                            if (atom_idx == db_pair.first.first || atom_idx == db_pair.first.second) {
                                stereo_score = 2;
                                break;
                            }
                        }
                    }
                }
                
                // Ring score (from _accountForRing)
                // Python: if atom.GetIsAromatic(): return 1; elif atom.IsInRing(): return 2; else: return 1
                int ring_score = 1;
                if (atom->getIsAromatic()) {
                    ring_score = 1;  // Aromatic rings not promoted
                } else {
                    // Check if atom is in a ring (like Osmordred does)
                    const RDKit::RingInfo *ri = molCopy.getRingInfo();
                    if (ri && ri->isInitialized() && ri->numAtomRings(atom_idx) > 0) {
                        ring_score = 2;  // Non-aromatic rings
                    }
                }
                
                // Bond score (neighbor score) - squared degree
                int bond_score = atom->getDegree();
                bond_score = bond_score * bond_score;
                
                // Total score for this atom (from _calculateScoreForAtom)
                double atom_score = static_cast<double>(hyb_score * stereo_score * ring_score * bond_score);
                total_score += atom_score;
            }
            
            // Normalize by number of heavy atoms (nSPS) - default normalize=True
            sps_value = total_score / static_cast<double>(nHeavy);
        }
    } catch (...) {
        sps_value = 0.0;
    }
    descriptors.push_back(sps_value);
    
    // 6: MolWt
    descriptors.push_back(MW);
    
    // 7: HeavyAtomMolWt
    descriptors.push_back(heavyMW);
    
    // 8: ExactMolWt
    descriptors.push_back(exactMW);
    
    // 9: NumValenceElectrons - Sum of valence electrons for all atoms
    // Python: sum(tbl.GetNOuterElecs(atom.GetAtomicNum()) - atom.GetFormalCharge() + atom.GetTotalNumHs() for atom in mol.GetAtoms())
    unsigned int totalValence = 0;
    const PeriodicTable* tbl = PeriodicTable::getTable();
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        const Atom* atom = mol.getAtomWithIdx(i);
        unsigned int nOuterElecs = tbl->getNouterElecs(atom->getAtomicNum());
        int formalCharge = atom->getFormalCharge();
        unsigned int totalNumHs = atom->getTotalNumHs();
        totalValence += nOuterElecs - formalCharge + totalNumHs;
    }
    descriptors.push_back(static_cast<double>(totalValence));
    
    // 10: NumRadicalElectrons - Sum of radical electrons for all atoms
    unsigned int totalRadicals = 0;
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        const Atom* atom = mol.getAtomWithIdx(i);
        totalRadicals += atom->getNumRadicalElectrons();
    }
    descriptors.push_back(static_cast<double>(totalRadicals));
    
    // 11-14: Partial charge descriptors
    // Literal port of Descriptors.py:_ChargeDescriptors:
    //   ComputeGasteigerCharges(mol)            (nIter=12, throwOnParamFailure=False)
    //   minChg, maxChg = 500., -500.
    //   for at in mol.GetAtoms():
    //     chg = float(at.GetProp('_GasteigerCharge'))
    //     minChg = min(chg, minChg); maxChg = max(chg, maxChg)
    // Charges are always recomputed (Python never reuses existing props). The
    // GetProp string round trip is exact (17 significant digits), so the double
    // is read directly. Python's builtin min(a, b) returns b if b < a else a
    // (likewise max with >); with a NaN charge the comparison is false, so the
    // running value is reset to the NaN atom and recovers at the next atom.
    // This comparison order is reproduced exactly, including NaN propagation.
    double maxCharge = std::numeric_limits<double>::quiet_NaN();
    double minCharge = std::numeric_limits<double>::quiet_NaN();
    double maxAbsCharge = std::numeric_limits<double>::quiet_NaN();
    double minAbsCharge = std::numeric_limits<double>::quiet_NaN();
    try {
        RDKit::computeGasteigerCharges(mol, 12, false);
        double minChg = 500.0;
        double maxChg = -500.0;
        for (const auto atom : mol.atoms()) {
            const double chg = atom->getProp<double>(common_properties::_GasteigerCharge);
            minChg = (minChg < chg) ? minChg : chg;  // min(chg, minChg)
            maxChg = (maxChg > chg) ? maxChg : chg;  // max(chg, maxChg)
        }
        minCharge = minChg;
        maxCharge = maxChg;
        const double a1 = std::fabs(minChg);
        const double a2 = std::fabs(maxChg);
        maxAbsCharge = (a2 > a1) ? a2 : a1;  // max(abs(v1), abs(v2))
        minAbsCharge = (a2 < a1) ? a2 : a1;  // min(abs(v1), abs(v2))
    } catch (...) {
        // Python raises -> CalcMolDescriptors missing value (NaN)
    }

    descriptors.push_back(maxCharge);      // MaxPartialCharge
    descriptors.push_back(minCharge);      // MinPartialCharge
    descriptors.push_back(maxAbsCharge);   // MaxAbsPartialCharge
    descriptors.push_back(minAbsCharge);   // MinAbsPartialCharge
    
    // 15-17: FpDensityMorgan - Fingerprint density = (num nonzero elements) / (num heavy atoms)
    // Python: _FingerprintDensity(mol, _getMorganCountFingerprint, radius)
    // Returns: len(fp.GetNonzeroElements()) / mol.GetNumHeavyAtoms()
    unsigned int numHeavy = mol.getNumHeavyAtoms();
    if (numHeavy == 0) {
        descriptors.push_back(0.0);  // FpDensityMorgan1
        descriptors.push_back(0.0);  // FpDensityMorgan2
        descriptors.push_back(0.0);  // FpDensityMorgan3
    } else {
        // Use RDKit's Morgan fingerprint generator for radius 1, 2, 3
        // unique_ptr ensures the generator is freed (fixes leak from getMorganGenerator)
        for (unsigned int radius = 1; radius <= 3; ++radius) {
            try {
                std::unique_ptr<FingerprintGenerator<std::uint32_t>> mgen(
                    MorganFingerprint::getMorganGenerator<std::uint32_t>(radius));
                if (!mgen) {
                    descriptors.push_back(0.0);
                    continue;
                }
                FingerprintFuncArguments args;
                auto fp = mgen->getSparseCountFingerprint(mol, args);
                if (!fp) {
                    descriptors.push_back(0.0);
                    continue;
                }
                // Get number of nonzero elements (like Python's len(fp.GetNonzeroElements()))
                unsigned int numNonzero = fp->getNonzeroElements().size();
                double density = static_cast<double>(numNonzero) / static_cast<double>(numHeavy);
                descriptors.push_back(density);
            } catch (...) {
                // If fingerprint generation fails for this radius, use zero
                descriptors.push_back(0.0);
            }
        }
    }
    
    // 18-25: BCUT2D descriptors
    // Use existing BCUT2D implementation from BCUT.cpp; on exception set NaN (match Python CalcMolDescriptors)
    // Reference: cpp/snn_features_bcut_fix.cpp — wrap BCUT in try/catch; on exception fill all BCUT slots with NaN.
    // Returns: [MWHI, MWLOW, CHGHI, CHGLO, LOGPHI, LOGPLOW, MRHI, MRLOW]
    {
        const double bcutNaN = std::numeric_limits<double>::quiet_NaN();
        try {
            std::vector<double> bcut_values = BCUT2D(mol);
            if (bcut_values.size() == 8) {
                descriptors.push_back(bcut_values[0]);  // BCUT2D_MWHI
                descriptors.push_back(bcut_values[1]);  // BCUT2D_MWLOW
                descriptors.push_back(bcut_values[2]);  // BCUT2D_CHGHI
                descriptors.push_back(bcut_values[3]);  // BCUT2D_CHGLO
                descriptors.push_back(bcut_values[4]);  // BCUT2D_LOGPHI
                descriptors.push_back(bcut_values[5]);  // BCUT2D_LOGPLOW
                descriptors.push_back(bcut_values[6]);  // BCUT2D_MRHI
                descriptors.push_back(bcut_values[7]);  // BCUT2D_MRLOW
            } else {
                for (int i = 0; i < 8; ++i) descriptors.push_back(bcutNaN);
            }
        } catch (const std::exception&) {
            for (int i = 0; i < 8; ++i) descriptors.push_back(bcutNaN);
        } catch (...) {
            for (int i = 0; i < 8; ++i) descriptors.push_back(bcutNaN);
        }
    }
    // 26-28: AvgIpc, BalabanJ, BertzCT -- literal ports of GraphDescriptors.py
    // (pure Python there), evaluated in CalcMolDescriptors order.
    double avgIpcValue = kNaN;
    double ipcValue = kNaN;
    try {
        const std::vector<double> cPoly = ipcAbsCharPoly(mol);
        const double entropy = infoEntropy(cPoly);
        double cPolySum = 0;
        for (double v : cPoly) {
            cPolySum = cPolySum + v;
        }
        avgIpcValue = entropy;
        ipcValue = cPolySum * entropy;
    } catch (...) {
    }
    descriptors.push_back(avgIpcValue);  // 26: AvgIpc

    double balabanJ = kNaN;
    try {
        balabanJ = pyBalabanJ(mol);
    } catch (...) {
    }
    descriptors.push_back(balabanJ);  // 27: BalabanJ

    double bertzCT = kNaN;
    try {
        bertzCT = pyBertzCT(mol);
    } catch (...) {
    }
    descriptors.push_back(bertzCT);  // 28: BertzCT
    
    // 29: Chi0 - Python uses sum(sqrt(1/degree)) for all atoms with degree > 0
    // From equations (1),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    double chi0 = 0.0;
    try {
        for (auto atom : mol.atoms()) {
            unsigned int degree = atom->getDegree();
            if (degree > 0) {
                chi0 += std::sqrt(1.0 / degree);
            }
        }
    } catch (...) {
        chi0 = 0.0;
    }
    descriptors.push_back(chi0);
    
    // 30: Chi0n
    descriptors.push_back(calcChi0n(mol));
    
    // 31: Chi0v
    descriptors.push_back(calcChi0v(mol));
    
    // 32: Chi1 - Python uses sum(sqrt(1/(deg1*deg2))) for all bonds with deg1*deg2 > 0
    // From equations (1),(11) and (12) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    double chi1 = 0.0;
    try {
        for (auto bond : mol.bonds()) {
            unsigned int deg1 = bond->getBeginAtom()->getDegree();
            unsigned int deg2 = bond->getEndAtom()->getDegree();
            unsigned int product = deg1 * deg2;
            if (product > 0) {
                chi1 += std::sqrt(1.0 / product);
            }
        }
    } catch (...) {
        chi1 = 0.0;
    }
    descriptors.push_back(chi1);
    
    // 33: Chi1n
    descriptors.push_back(calcChi1n(mol));
    
    // 34: Chi1v
    descriptors.push_back(calcChi1v(mol));
    
    // 35: Chi2n
    descriptors.push_back(calcChi2n(mol));
    
    // 36: Chi2v
    descriptors.push_back(calcChi2v(mol));
    
    // 37: Chi3n
    descriptors.push_back(calcChi3n(mol));
    
    // 38: Chi3v
    descriptors.push_back(calcChi3v(mol));
    
    // 39: Chi4n
    descriptors.push_back(calcChi4n(mol));
    
    // 40: Chi4v
    descriptors.push_back(calcChi4v(mol));
    
    // 41: HallKierAlpha
    descriptors.push_back(calcHallKierAlpha(mol));
    
    // 42: Ipc = sum(cPoly) * InfoEntropy(cPoly)  (cPoly computed with AvgIpc above)
    descriptors.push_back(ipcValue);
    
    // 43: Kappa1
    descriptors.push_back(calcKappa1(mol));
    
    // 44: Kappa2
    descriptors.push_back(calcKappa2(mol));
    
    // 45: Kappa3
    descriptors.push_back(calcKappa3(mol));
    
    // 46: LabuteASA
    descriptors.push_back(calcLabuteASA(mol));
    
    // 47-60: PEOE_VSA1-14
    // Python order: 1, 10, 11, 12, 13, 14, 2, 3, 4, 5, 6, 7, 8, 9
    // C++ indices:  0,  9, 10, 11, 12, 13, 1, 2, 3, 4, 5, 6, 7, 8
    try {
        std::vector<double> peoeVSA = RDKit::Descriptors::calcPEOE_VSA(mol);
        // Python _descList order: PEOE_VSA1, PEOE_VSA10, PEOE_VSA11, PEOE_VSA12, PEOE_VSA13, PEOE_VSA14, PEOE_VSA2-9
        int pythonOrder[14] = {0, 9, 10, 11, 12, 13, 1, 2, 3, 4, 5, 6, 7, 8};
        for (int i = 0; i < 14; ++i) {
            int cppIdx = pythonOrder[i];
            descriptors.push_back(cppIdx < static_cast<int>(peoeVSA.size()) ? peoeVSA[cppIdx] : 0.0);
        }
    } catch (...) {
        for (int i = 0; i < 14; ++i) {
            descriptors.push_back(0.0);
        }
    }
    
    // 61-70: SMR_VSA1-10
    // Python order: 1, 10, 2, 3, 4, 5, 6, 7, 8, 9
    // C++ indices:  0,  9, 1, 2, 3, 4, 5, 6, 7, 8
    try {
        std::vector<double> smrVSA = RDKit::Descriptors::calcSMR_VSA(mol);
        int pythonOrder[10] = {0, 9, 1, 2, 3, 4, 5, 6, 7, 8};
        for (int i = 0; i < 10; ++i) {
            int cppIdx = pythonOrder[i];
            descriptors.push_back(cppIdx < static_cast<int>(smrVSA.size()) ? smrVSA[cppIdx] : 0.0);
        }
    } catch (...) {
        for (int i = 0; i < 10; ++i) {
            descriptors.push_back(0.0);
        }
    }
    
    // 71-82: SlogP_VSA1-12
    // Python order: 1, 10, 11, 12, 2, 3, 4, 5, 6, 7, 8, 9
    // C++ indices:  0,  9, 10, 11, 1, 2, 3, 4, 5, 6, 7, 8
    try {
        std::vector<double> slogpVSA = RDKit::Descriptors::calcSlogP_VSA(mol);
        int pythonOrder[12] = {0, 9, 10, 11, 1, 2, 3, 4, 5, 6, 7, 8};
        for (int i = 0; i < 12; ++i) {
            int cppIdx = pythonOrder[i];
            descriptors.push_back(cppIdx < static_cast<int>(slogpVSA.size()) ? slogpVSA[cppIdx] : 0.0);
        }
    } catch (...) {
        for (int i = 0; i < 12; ++i) {
            descriptors.push_back(0.0);
        }
    }
    
    // 83: TPSA
    descriptors.push_back(calcTPSA(mol));
    
    // 84-94: EState_VSA1-11 - Use exported calcEState_VSA from Osmordred
    // Python order: 1, 10, 11, 2, 3, 4, 5, 6, 7, 8, 9
    // C++ indices:  0,  9, 10, 1, 2, 3, 4, 5, 6, 7, 8
    try {
        std::vector<double> estateVSA = calcEState_VSA(mol);
        int pythonOrder[11] = {0, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8};
        for (int i = 0; i < 11; ++i) {
            int cppIdx = pythonOrder[i];
            descriptors.push_back(cppIdx < static_cast<int>(estateVSA.size()) ? estateVSA[cppIdx] : 0.0);
        }
    } catch (...) {
        for (int i = 0; i < 11; ++i) {
            descriptors.push_back(0.0);
        }
    }
    
    // 95-104: VSA_EState1-10 - Use exported calcVSA_EState from Osmordred
    // Python order: 1, 10, 2, 3, 4, 5, 6, 7, 8, 9
    // C++ indices:  0,  9, 1, 2, 3, 4, 5, 6, 7, 8
    try {
        std::vector<double> vsaEState = calcVSA_EState(mol);
        int pythonOrder[10] = {0, 9, 1, 2, 3, 4, 5, 6, 7, 8};
        for (int i = 0; i < 10; ++i) {
            int cppIdx = pythonOrder[i];
            descriptors.push_back(cppIdx < static_cast<int>(vsaEState.size()) ? vsaEState[cppIdx] : 0.0);
        }
    } catch (...) {
        for (int i = 0; i < 10; ++i) {
            descriptors.push_back(0.0);
        }
    }
    
    // 105: FractionCSP3
    descriptors.push_back(calcFractionCSP3(mol));
    
    // 106: HeavyAtomCount
    descriptors.push_back(static_cast<double>(nHeavy));
    
    // 107-131: Count descriptors
    // NHOHCount - Python uses rdMolDescriptors.CalcNumLipinskiHBD
    unsigned int nhohCount = calcLipinskiHBD(mol);  // Use calcLipinskiHBD, not calcNumHBD
    descriptors.push_back(static_cast<double>(nhohCount));  // NHOHCount
    
    // NOCount - Count N or O atoms
    unsigned int noCount = 0;
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        const Atom* atom = mol.getAtomWithIdx(i);
        unsigned int atomicNum = atom->getAtomicNum();
        if (atomicNum == 7 || atomicNum == 8) {  // N or O
            noCount++;
        }
    }
    descriptors.push_back(static_cast<double>(noCount));     // NOCount
    descriptors.push_back(calcNumAliphaticCarbocycles(mol));
    descriptors.push_back(calcNumAliphaticHeterocycles(mol));
    descriptors.push_back(calcNumAliphaticRings(mol));
    descriptors.push_back(calcNumAmideBonds(mol));
    descriptors.push_back(calcNumAromaticCarbocycles(mol));
    descriptors.push_back(calcNumAromaticHeterocycles(mol));
    descriptors.push_back(calcNumAromaticRings(mol));
    descriptors.push_back(static_cast<double>(numAtomStereoCenters(mol)));
    descriptors.push_back(calcNumBridgeheadAtoms(mol));
    descriptors.push_back(static_cast<double>(calcNumHBA(mol)));  // NumHAcceptors
    descriptors.push_back(static_cast<double>(calcNumHBD(mol)));  // NumHDonors
    descriptors.push_back(static_cast<double>(calcNumHeteroatoms(mol)));  // NumHeteroatoms
    descriptors.push_back(static_cast<double>(calcNumHeterocycles(mol)));  // NumHeterocycles
    descriptors.push_back(static_cast<double>(calcNumRotatableBonds(mol)));  // NumRotatableBonds
    descriptors.push_back(calcNumSaturatedCarbocycles(mol));
    descriptors.push_back(calcNumSaturatedHeterocycles(mol));
    descriptors.push_back(calcNumSaturatedRings(mol));
    descriptors.push_back(calcNumSpiroAtoms(mol));
    descriptors.push_back(static_cast<double>(numUnspecifiedAtomStereoCenters(mol)));
    
    // 128: Phi
    descriptors.push_back(calcPhi(mol));
    
    // 129: RingCount
    descriptors.push_back(static_cast<double>(calcNumRings(mol)));
    
    // 130: MolLogP
    descriptors.push_back(logP);
    
    // 131: MolMR
    descriptors.push_back(MR);
    
    // 132-216: Fragment counts (fr_*) - 85 fragment descriptors
    // These are functional group counts using SMARTS pattern matching
    // Order matches Python's Descriptors._descList: fr_Al_COO, fr_Al_OH, ..., fr_urea
    // SMARTS patterns from FragmentDescriptors.csv
    // Use shared_ptr like Osmordred does for safe memory management
    auto& fragmentQueries = GetFragmentQueries();
    SubstructMatchParameters fragParams;  // uniquify=true, maxMatches=1000 (Python defaults)
    for (const auto& query : fragmentQueries) {
        double count = std::numeric_limits<double>::quiet_NaN();
        if (query) {
            try {
                count = static_cast<double>(RDKit::SubstructMatch(mol, *query, fragParams).size());
            } catch (...) {
            }
        }
        descriptors.push_back(count);
    }
    
    // Ensure exactly 217 descriptors
    if (descriptors.size() < 217) {
        descriptors.resize(217, 0.0);
    } else if (descriptors.size() > 217) {
        descriptors.resize(217);
    }
    
    return descriptors;
}

// Batch version: Extract RDKit descriptors for multiple molecules using parallel processing
// UNIFIED PARALLEL IMPLEMENTATION: Uses same approach as calcPhysChemPropBatch
RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>> extractRDKitDescriptorsBatch(
    const std::vector<std::string>& smiles_list, int n_jobs) {
    
    std::vector<std::vector<double>> results;
    results.reserve(smiles_list.size());
    
    // Determine number of threads to use
    unsigned int nThreads = getNumThreadsToUse(n_jobs);
    
    // For small batches, use sequential processing (avoid async overhead)
    if (nThreads <= 1 || smiles_list.size() < 10) {
        for (const auto& smi : smiles_list) {
            ROMol* mol = SmilesToMol(smi);
            if (mol) {
                results.push_back(extractRDKitDescriptors(*mol));
                delete mol;
            } else {
                // Return vector of zeros for invalid SMILES
                results.push_back(std::vector<double>(217, 0.0));
            }
        }
        return results;
    }
    
    // Parallel processing using std::async (same approach as calcPhysChemPropBatch)
    std::vector<std::future<std::vector<double>>> futures;
    futures.reserve(smiles_list.size());
    
    // Process molecules in parallel
    for (size_t idx = 0; idx < smiles_list.size(); ++idx) {
        const auto& smi = smiles_list[idx];
        
        futures.emplace_back(std::async(std::launch::async, [smi]() {
            try {
                ROMol* mol = SmilesToMol(smi);
                if (mol) {
                    try {
                        std::vector<double> descriptors = extractRDKitDescriptors(*mol);
                        delete mol;
                        return descriptors;
                    } catch (...) {
                        // If extractRDKitDescriptors fails, return zeros and continue processing
                        delete mol;
                        return std::vector<double>(217, 0.0);
                    }
                } else {
                    return std::vector<double>(217, 0.0);
                }
            } catch (...) {
                // Catch any other exceptions (e.g., memory errors) and return zeros
                return std::vector<double>(217, 0.0);
            }
        }));
    }
    
    // Collect results
    for (auto& f : futures) {
        results.push_back(f.get());
    }
    
    return results;
}

// Batch version: Extract RDKit descriptors from mol objects directly (SNN support)
// Handles nullptr gracefully (returns NaN row for that molecule — same as Python missing value)
RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>> extractRDKitDescriptorsFromMolsBatch(
    const std::vector<const ROMol*>& mols, int n_jobs) {
    
    const double kNaN = std::numeric_limits<double>::quiet_NaN();
    const std::vector<double> nanRow(217, kNaN);

    std::vector<std::vector<double>> results;
    results.reserve(mols.size());

    // Determine number of threads to use
    unsigned int nThreads = getNumThreadsToUse(n_jobs);

    // For small batches, use sequential processing (avoid async overhead)
    if (nThreads <= 1 || mols.size() < 10) {
        for (const ROMol* mol : mols) {
            if (mol) {
                try {
                    results.push_back(extractRDKitDescriptors(*mol));
                    if (results.back().size() != 217u) results.back() = nanRow;
                } catch (...) {
                    results.push_back(nanRow);
                }
            } else {
                results.push_back(nanRow);
            }
        }
        return results;
    }

    // Parallel processing using std::async (same approach as extractRDKitDescriptorsBatch)
    std::vector<std::future<std::vector<double>>> futures;
    futures.reserve(mols.size());

    for (size_t idx = 0; idx < mols.size(); ++idx) {
        const ROMol* mol = mols[idx];

        futures.emplace_back(std::async(std::launch::async, [mol, kNaN]() {
            if (mol) {
                try {
                    std::vector<double> row = extractRDKitDescriptors(*mol);
                    if (row.size() != 217u) return std::vector<double>(217, kNaN);
                    return row;
                } catch (...) {
                    return std::vector<double>(217, kNaN);
                }
            }
            return std::vector<double>(217, kNaN);
        }));
    }

    for (auto& f : futures) {
        results.push_back(f.get());
    }

    return results;
}

// Get descriptor names in the same order as extractRDKitDescriptors returns values
// This matches Python's Descriptors._descList order exactly
RDKIT_DESCRIPTORS_EXPORT std::vector<std::string> getRDKit217DescriptorNames() {
    // Return the 217 descriptor names in exact order matching extractRDKitDescriptors
    // This order matches Python's Descriptors._descList
    return {
        "MaxAbsEStateIndex", "MaxEStateIndex", "MinAbsEStateIndex", "MinEStateIndex",
        "qed", "SPS",
        "MolWt", "HeavyAtomMolWt", "ExactMolWt",
        "NumValenceElectrons", "NumRadicalElectrons",
        "MaxPartialCharge", "MinPartialCharge", "MaxAbsPartialCharge", "MinAbsPartialCharge",
        "FpDensityMorgan1", "FpDensityMorgan2", "FpDensityMorgan3",
        "BCUT2D_MWHI", "BCUT2D_MWLOW", "BCUT2D_CHGHI", "BCUT2D_CHGLO",
        "BCUT2D_LOGPHI", "BCUT2D_LOGPLOW", "BCUT2D_MRHI", "BCUT2D_MRLOW",
        "AvgIpc", "BalabanJ", "BertzCT",
        "Chi0", "Chi0n", "Chi0v", "Chi1", "Chi1n", "Chi1v",
        "Chi2n", "Chi2v", "Chi3n", "Chi3v", "Chi4n", "Chi4v",
        "HallKierAlpha", "Ipc", "Kappa1", "Kappa2", "Kappa3", "LabuteASA",
        "PEOE_VSA1", "PEOE_VSA10", "PEOE_VSA11", "PEOE_VSA12", "PEOE_VSA13", "PEOE_VSA14",
        "PEOE_VSA2", "PEOE_VSA3", "PEOE_VSA4", "PEOE_VSA5", "PEOE_VSA6", "PEOE_VSA7", "PEOE_VSA8", "PEOE_VSA9",
        "SMR_VSA1", "SMR_VSA10", "SMR_VSA2", "SMR_VSA3", "SMR_VSA4", "SMR_VSA5",
        "SMR_VSA6", "SMR_VSA7", "SMR_VSA8", "SMR_VSA9",
        "SlogP_VSA1", "SlogP_VSA10", "SlogP_VSA11", "SlogP_VSA12",
        "SlogP_VSA2", "SlogP_VSA3", "SlogP_VSA4", "SlogP_VSA5", "SlogP_VSA6", "SlogP_VSA7", "SlogP_VSA8", "SlogP_VSA9",
        "TPSA",
        "EState_VSA1", "EState_VSA10", "EState_VSA11", "EState_VSA2", "EState_VSA3", "EState_VSA4",
        "EState_VSA5", "EState_VSA6", "EState_VSA7", "EState_VSA8", "EState_VSA9",
        "VSA_EState1", "VSA_EState10", "VSA_EState2", "VSA_EState3", "VSA_EState4",
        "VSA_EState5", "VSA_EState6", "VSA_EState7", "VSA_EState8", "VSA_EState9",
        "FractionCSP3", "HeavyAtomCount",
        "NHOHCount", "NOCount",
        "NumAliphaticCarbocycles", "NumAliphaticHeterocycles", "NumAliphaticRings", "NumAmideBonds",
        "NumAromaticCarbocycles", "NumAromaticHeterocycles", "NumAromaticRings",
        "NumAtomStereoCenters", "NumBridgeheadAtoms",
        "NumHAcceptors", "NumHDonors", "NumHeteroatoms", "NumHeterocycles", "NumRotatableBonds",
        "NumSaturatedCarbocycles", "NumSaturatedHeterocycles", "NumSaturatedRings",
        "NumSpiroAtoms", "NumUnspecifiedAtomStereoCenters",
        "Phi", "RingCount",
        "MolLogP", "MolMR",
        "fr_Al_COO", "fr_Al_OH", "fr_Al_OH_noTert", "fr_ArN", "fr_Ar_COO", "fr_Ar_N", "fr_Ar_NH", "fr_Ar_OH",
        "fr_COO", "fr_COO2", "fr_C_O", "fr_C_O_noCOO", "fr_C_S", "fr_HOCCN", "fr_Imine",
        "fr_NH0", "fr_NH1", "fr_NH2", "fr_N_O", "fr_Ndealkylation1", "fr_Ndealkylation2", "fr_Nhpyrrole", "fr_SH",
        "fr_aldehyde", "fr_alkyl_carbamate", "fr_alkyl_halide", "fr_allylic_oxid", "fr_amide", "fr_amidine", "fr_aniline",
        "fr_aryl_methyl", "fr_azide", "fr_azo", "fr_barbitur", "fr_benzene", "fr_benzodiazepine", "fr_bicyclic",
        "fr_diazo", "fr_dihydropyridine", "fr_epoxide", "fr_ester", "fr_ether", "fr_furan", "fr_guanido", "fr_halogen",
        "fr_hdrzine", "fr_hdrzone", "fr_imidazole", "fr_imide", "fr_isocyan", "fr_isothiocyan", "fr_ketone", "fr_ketone_Topliss",
        "fr_lactam", "fr_lactone", "fr_methoxy", "fr_morpholine", "fr_nitrile", "fr_nitro", "fr_nitro_arom", "fr_nitro_arom_nonortho",
        "fr_nitroso", "fr_oxazole", "fr_oxime", "fr_para_hydroxylation", "fr_phenol", "fr_phenol_noOrthoHbond",
        "fr_phos_acid", "fr_phos_ester", "fr_piperdine", "fr_piperzine", "fr_priamide", "fr_prisulfonamd", "fr_pyridine", "fr_quatN",
        "fr_sulfide", "fr_sulfonamd", "fr_sulfone", "fr_term_acetylene", "fr_tetrazole", "fr_thiazole", "fr_thiocyan", "fr_thiophene",
        "fr_unbrch_alkane", "fr_urea"
    };
}

}  // namespace Osmordred
}  // namespace Descriptors
}  // namespace RDKit

