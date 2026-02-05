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

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// Helper to safely compute descriptor, return 0.0 on error
template<typename Func>
double safeCompute(Func func, const ROMol& mol) {
    try {
        return func(mol);
    } catch (...) {
        return 0.0;
    }
}

// Helper to get first element of vector descriptor, return 0.0 if empty
double getFirst(const std::vector<double>& vec) {
    return vec.empty() ? 0.0 : vec[0];
}

// Helper: Information Entropy (Shannon entropy in base 2)
// From rdkit.ML.InfoTheory.entropy.InfoEntropy
template <class T>
double calcInfoEntropy(const std::vector<T>& data) {
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

// Helper: Characteristic Polynomial using Le Verrier-Faddeev-Frame method
// From rdkit.Chem.Graphs.CharacteristicPolynomial
// Returns coefficients of characteristic polynomial: [c0, c1, c2, ..., cn]
// where det(A - λI) = c0 + c1*λ + c2*λ² + ... + cn*λⁿ
std::vector<double> characteristicPolynomial(const ROMol& mol, const std::vector<std::vector<double>>& adjMat) {
    unsigned int nAtoms = mol.getNumAtoms();
    std::vector<double> res(nAtoms + 1, 0.0);
    res[0] = 1.0;
    
    if (nAtoms == 0) return res;
    
    // Identity matrix
    std::vector<std::vector<double>> I(nAtoms, std::vector<double>(nAtoms, 0.0));
    for (unsigned int i = 0; i < nAtoms; ++i) {
        I[i][i] = 1.0;
    }
    
    // An = A (adjacency matrix)
    std::vector<std::vector<double>> An = adjMat;
    
    // Le Verrier-Faddeev-Frame method
    for (unsigned int n = 1; n <= nAtoms; ++n) {
        // Calculate trace of An
        double trace = 0.0;
        for (unsigned int i = 0; i < nAtoms; ++i) {
            trace += An[i][i];
        }
        
        res[n] = trace / static_cast<double>(n);
        
        // Bn = An - res[n] * I
        std::vector<std::vector<double>> Bn(nAtoms, std::vector<double>(nAtoms, 0.0));
        for (unsigned int i = 0; i < nAtoms; ++i) {
            for (unsigned int j = 0; j < nAtoms; ++j) {
                Bn[i][j] = An[i][j] - res[n] * I[i][j];
            }
        }
        
        // An = A * Bn (matrix multiplication)
        std::vector<std::vector<double>> AnNew(nAtoms, std::vector<double>(nAtoms, 0.0));
        for (unsigned int i = 0; i < nAtoms; ++i) {
            for (unsigned int j = 0; j < nAtoms; ++j) {
                for (unsigned int k = 0; k < nAtoms; ++k) {
                    AnNew[i][j] += adjMat[i][k] * Bn[k][j];
                }
            }
        }
        An = AnNew;
    }
    
    // Negate coefficients (except c0)
    for (unsigned int i = 1; i <= nAtoms; ++i) {
        res[i] *= -1.0;
    }
    
    return res;
}

// Helper: Get adjacency matrix from molecule (1 if bonded, 0 otherwise)
std::vector<std::vector<double>> getAdjacencyMatrix(const ROMol& mol) {
    unsigned int nAtoms = mol.getNumAtoms();
    std::vector<std::vector<double>> adjMat(nAtoms, std::vector<double>(nAtoms, 0.0));
    
    for (auto bond : mol.bonds()) {
        unsigned int i = bond->getBeginAtomIdx();
        unsigned int j = bond->getEndAtomIdx();
        adjMat[i][j] = 1.0;
        adjMat[j][i] = 1.0;
    }
    
    return adjMat;
}

// Helper: Calculate Ipc (Information Content of characteristic polynomial coefficients)
double calcIpc(const ROMol& mol, bool avg = false) {
    try {
        // Get adjacency matrix (1 if bonded, 0 otherwise)
        std::vector<std::vector<double>> adjMat = getAdjacencyMatrix(mol);
        
        // Calculate characteristic polynomial
        std::vector<double> cPoly = characteristicPolynomial(mol, adjMat);
        
        // Take absolute values
        std::vector<double> absCPoly;
        for (double val : cPoly) {
            absCPoly.push_back(std::abs(val));
        }
        
        // Calculate information entropy
        double entropy = calcInfoEntropy(absCPoly);
        
        if (avg) {
            return entropy;
        } else {
            // Sum of coefficients * entropy
            double sum = 0.0;
            for (double val : absCPoly) {
                sum += val;
            }
            return sum * entropy;
        }
    } catch (...) {
        return 0.0;
    }
}

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

const std::vector<std::shared_ptr<RWMol>>& GetFragmentQueries() {
    static const std::vector<std::shared_ptr<RWMol>> queries = [] {
        std::vector<std::shared_ptr<RWMol>> res;
        std::vector<std::string> fragment_smarts = {
            "C-C(=O)[O;H1,-]", "[C!$(C=O)]-[OH]", "[$(C-[OX2H]);!$([CX3](-[OX2H])=[OX1]);!$([CD4]-[OX2H])]",
            "[$(a-[NX3H2]),$(a-[NH1][NH2]),$(a-C(=[OX1])[NH1][NH2]),$(a-C(=[NH])[NH2])]", "c-C(=O)[O;H1,-]",
            "n", "[nH]", "c[OH1]", "[#6]C(=O)[O;H,-1]", "[CX3](=O)[OX1H0-,OX2H1]", "[CX3]=[OX1]",
            "[C!$(C-[OH])]=O", "C=[SX1]", "[$([OX2H1][CX4][CX4H2][NX3&R1]),$([OH1][CX4][CX4H2][NX3][CX4](C)(C)C)]",
            "[Nv3](=C)-[#6]", "[NH0,nH0]", "[NH1,nH1]", "[NH2,nH2]", "[N!$(N=O)](-O)-C",
            "[$(N(-[CH3])-C-[$(C~O),$(C-a),$(C-N),$(C=C)]),$(N(-[CH2][CH3])-C-[$(C~O),$(C-a),$(C-N),$(C=C)])]",
            "[$([N&R1]1(-C)CCC1),$([N&R1]1(-C)CCCC1),$([N&R1]1(-C)CCCCC1),$([N&R1]1(-C)CCCCCC1),$([N&R1]1(-C)CCCCCCC1)]",
            "[nH]", "[SH]", "[CX3H1](=O)[#6]", "C[NH1]C(=O)OC", "[CX4]-[Cl,Br,I,F]",
            "[$(C=C-C);!$(C=C-C-[N,O,S]);!$(C=C-C-C-[N,O]);!$(C12=CC(=O)CCC1C3C(C4C(CCC4)CC3)CC2)]",
            "C(=O)-N", "C(=N)(-N)-[!#7]", "c-[NX3;!$(N=*)]",
            "[$(a-[CH3]),$(a-[CH2]-[CH3]),$(a-[CH2]-[CH2]~[!N;!O]);!$(a(:a!:*):a!:*)]",
            "[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]", "[#6]-N=N-[#6]", "C1C(=O)NC(=O)NC1=O",
            "c1ccccc1", "[c&R2]12[c&R1][c&R1][c&R1][c&R1][c&R2]1[N&R1][C&R1][C&R1][N&R1]=[C&R1]2", "[R2][R2]",
            "[N+]#N", "[$([NX3H1]1-C=C-C-C=C1),$([Nv3]1=C-C-C=C-C1),$([Nv3]1=C-C=C-C-C1),$([NX3H1]1-C-C=C-C=C1)]",
            "O1CC1", "[#6][CX3](=O)[OX2H0][#6]", "[OD2]([#6])[#6]", "o1cccc1", "C(=N)(N)N",
            "[#9,#17,#35,#53]", "[NX3]-[NX3]", "C=N-[NX3]", "n1cncc1", "N(-C(=O))-C=O", "N=C=O", "N=C=S",
            "[#6][CX3](=O)[#6]", "[$([CX3](=[OX1])(C)([c,C]));!$([CX3](=[OX1])([CH1]=C)[c,C])]", "N1C(=O)CC1",
            "[C&R1](=O)[O&R1][C&R1]", "[OX2](-[#6])-[CH3]", "O1CCNCC1", "[NX1]#[CX2]",
            "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]", "[$(c1(-[$([NX3](=O)=O),$([NX3+](=O)[O-])])ccccc1)]",
            "[$(c1(-[$([NX3](=O)=O),$([NX3+](=O)[O-])])ccccc1);!$(cc-!:*)]", "[N!$(N-O)]=O", "c1ocnc1",
            "[CX3]=[NX2]-[OX2]", "[$([cH]1[cH]cc(c[cH]1)~[$([#8,$([#8]~[H,c,C])])]),$([cH]1[cH]cc(c[cH]1)~[$([#7X3,$([#7](~[H,c,C])~[H,c,C])])]),$([cH]1[cH]cc(c[cH]1)-!:[$([NX3H,$(NC(=O)[H,c,C])])])]",
            "[OX2H]-c1ccccc1", "[$(c1(-[OX2H])ccccc1);!$(cc-!:[CH2]-[OX2H]);!$(cc-!:C(=O)[O;H1,-]);!$(cc-!:C(=O)-[NH2])]",
            "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]",
            "[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]",
            "N1CCCCC1", "N1CCNCC1", "C(=O)-[NH2]", "[NH2]-S(=,-[OX1;+0;-1])(=,-[OX1;+0;-1])-[#6]", "n1ccccc1",
            "[$([NX4+]),$([NX4]=*)]", "[SX2](-[#6])-C", "N-S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])-[#6]",
            "S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]", "C#[CH]", "c1nnnn1", "c1scnc1", "S-C#N",
            "s1cccc1", "[CR0;D2,D1][CR0;D2][CR0;D2][CR0;D2,D1]", "C(=O)(-N)-N"
        };
        res.reserve(fragment_smarts.size());
        for (const auto& smarts : fragment_smarts) {
            auto mol = RDKit::SmartsToMol(smarts);
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
    
    // Get BalabanJ - Python uses useBO=1, so we implement it directly here
    // Python: GetDistanceMatrix(mol, useBO=1, useAtomWts=0)
    double balabanJ = 0.0;
    try {
        double q = static_cast<double>(mol.getNumBonds());
        unsigned int n = mol.getNumAtoms();
        
        // Get distance matrix with useBO=1 (bond order) to match Python
        double* dMat = MolOps::getDistanceMat(mol, true, false, false);  // useBO=1, useAtomWts=0
        // Get adjacency matrix (useBO=0 for adjacency)
        double* adjMat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
        
        // Calculate vertex degrees s from distance matrix
        // Python's _VertexDegrees: row sum of distance matrix (sum of distances for each atom)
        std::vector<double> s(n, 0.0);
        for (unsigned int i = 0; i < n; ++i) {
            for (unsigned int j = 0; j < n; ++j) {
                double dist = dMat[i * n + j];
                if (dist < 1e6) {  // Valid distance
                    s[i] += dist;
                }
            }
        }
        
        double mu = q - n + 1;
        double sum_ = 0.0;
        for (unsigned int i = 0; i < n; ++i) {
            double si = s[i];
            for (unsigned int j = i; j < n; ++j) {
                if (adjMat[i * n + j] == 1) {  // Adjacent atoms
                    sum_ += 1.0 / std::sqrt(si * s[j]);
                }
            }
        }
        
        if (mu + 1 != 0) {
            balabanJ = (q / (mu + 1)) * sum_;
        }
        
        // NOTE: Do NOT delete[] dMat or adjMat - RDKit caches them in molecule properties
        // The documentation says "The caller should NOT delete this pointer"
    } catch (...) {
        balabanJ = 0.0;
    }
    
    // Get BertzCT (returns vector, need first element)
    std::vector<double> bertzCT_vec = calcBertzCT(mol);
    double bertzCT = getFirst(bertzCT_vec);
    
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
    // Python: MaxPartialCharge, MinPartialCharge, MaxAbsPartialCharge, MinAbsPartialCharge
    // Python's _ChargeDescriptors returns (minCharge, maxCharge), then:
    //   MaxAbsPartialCharge = max(abs(minCharge), abs(maxCharge))
    //   MinAbsPartialCharge = min(abs(minCharge), abs(maxCharge))
    //
    // CRITICAL: When mol comes from Python (e.g. ToBinary after create_rdkit_descriptors),
    // Python has already set _GasteigerCharge on every atom. Reusing those values (instead
    // of calling computeGasteigerCharges again) gives identical min/max and avoids ~3e-12
    // float differences from C++ vs Python Gasteiger iteration order.
    double maxCharge = std::numeric_limits<double>::quiet_NaN();
    double minCharge = std::numeric_limits<double>::quiet_NaN();
    double maxAbsCharge = std::numeric_limits<double>::quiet_NaN();
    double minAbsCharge = std::numeric_limits<double>::quiet_NaN();
    try {
        // Check if every atom already has _GasteigerCharge (e.g. mol from Python with charges set)
        unsigned int nAtoms = mol.getNumAtoms();
        bool allHaveCharge = (nAtoms > 0);
        for (unsigned int i = 0; i < nAtoms && allHaveCharge; ++i) {
            const Atom* atom = mol.getAtomWithIdx(i);
            double d;
            std::string s;
            if (!atom->getPropIfPresent("_GasteigerCharge", d) && !atom->getPropIfPresent("_GasteigerCharge", s))
                allHaveCharge = false;
        }
        if (!allHaveCharge)
            RDKit::computeGasteigerCharges(mol, 12, false);

        double minChg = 500.0;
        double maxChg = -500.0;
        for (unsigned int i = 0; i < nAtoms; ++i) {
            const Atom* atom = mol.getAtomWithIdx(i);
            double chg;
            bool haveChg = false;
            if (atom->getPropIfPresent("_GasteigerCharge", chg)) {
                haveChg = true;
            } else {
                std::string s;
                if (atom->getPropIfPresent("_GasteigerCharge", s)) {
                    chg = std::stod(s);
                    haveChg = true;
                }
            }
            if (haveChg && !std::isnan(chg)) {
                minChg = std::min(chg, minChg);
                maxChg = std::max(chg, maxChg);
            }
        }
        // Fallback: same as Gasteiger splitChargeConjugated — when no valid Gasteiger charge (e.g. metals),
        // use formal charge so [Hg+2], [Fe+2], etc. get min=max=formal (matches Python).
        if (minChg >= maxChg && nAtoms > 0) {
            minChg = 500.0;
            maxChg = -500.0;
            for (unsigned int i = 0; i < nAtoms; ++i) {
                double fc = static_cast<double>(mol.getAtomWithIdx(i)->getFormalCharge());
                if (fc < minChg) minChg = fc;
                if (fc > maxChg) maxChg = fc;
            }
        }
        // Only NaN when no atom had a valid charge (sentinels unchanged: minChg > maxChg). When minChg == maxChg (e.g. one atom, II) we have a valid result.
        if (minChg > maxChg) {
            minCharge = maxCharge = maxAbsCharge = minAbsCharge = std::numeric_limits<double>::quiet_NaN();
        } else {
            minCharge = minChg;
            maxCharge = maxChg;
            double absMin = std::abs(minCharge);
            double absMax = std::abs(maxCharge);
            maxAbsCharge = std::max(absMin, absMax);
            minAbsCharge = std::min(absMin, absMax);
        }
    } catch (...) {
        // If something fails, keep NaNs (so Python/C++ can agree on missingness).
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
    
    // 26: AvgIpc - Average Information Content
    // Uses characteristic polynomial of adjacency matrix (proper implementation)
    descriptors.push_back(calcIpc(mol, true));  // avg = true
    
    // 27: BalabanJ
    descriptors.push_back(balabanJ);
    
    // 28: BertzCT
    descriptors.push_back(bertzCT);
    
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
    
    // 42: Ipc - Information Content
    // Uses characteristic polynomial of adjacency matrix (proper implementation)
    descriptors.push_back(calcIpc(mol, false));  // avg = false
    
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
    // Use const ROMol& directly like Osmordred does (uniquify=true works with const ROMol&)
    for (const auto& query : fragmentQueries) {
        if (query) {
            try {
                std::vector<RDKit::MatchVectType> matches;
                RDKit::SubstructMatch(mol, *query, matches, true);  // uniquify = true, like Osmordred
                descriptors.push_back(static_cast<double>(matches.size()));
            } catch (...) {
                descriptors.push_back(0.0);
            }
        } else {
            descriptors.push_back(0.0);
        }
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

