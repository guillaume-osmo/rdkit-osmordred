//
//  C++ port of dmpnn-cpp-fp16/tools/featurize.py.
//
//  The script is the authority on the trained input format and this file follows
//  it line by line. Three things are load-bearing:
//    * The graph is built on AddHs(mol); the 217 descriptors on mol itself.
//    * MMFFMolProperties(work) MUTATES `work` (MMFF aromaticity model). The
//      script builds it before reading any atom, bond or whole-molecule feature,
//      so every one of those sees the MMFF-perceived molecule. Same order here.
//    * The descriptor block goes float64 -> float32, non-finite -> NaN,
//      float32 inf -> NaN, NaN -> 0, then arcsinh IN FLOAT32.
//
#include "DmpnnFeaturizer.h"

#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/Lipinski.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/MolSurf.h>
#include <GraphMol/Descriptors/rdkit217/RDKit217Descriptors.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <utility>

namespace RDKit {
namespace DmpnnFeaturizer {
namespace {

constexpr int MMFF_TYPE_BINS = 128;
const char *const MMFF_VARIANT = "MMFF94s";

// _FALLBACK_BY_ELEMENT, _HALOGEN_TYPES and the boron codes of featurize.py.
constexpr int H_TYPE = 5;
constexpr int B_SP2_C = 120, B_SP2_O = 121, B_SP2_N = 122;
constexpr int B_SP2_HALOGEN = 123, B_SP2_S = 124, B_TETRA = 125;
constexpr int B_RING = 126, UNKNOWN = 127;

int fallbackByElement(int z) {
  switch (z) {
    case 1: return 5;
    case 6: return 1;
    case 7: return 8;
    case 8: return 6;
    case 9: return 11;
    case 16: return 15;
    case 17: return 12;
    case 35: return 13;
    case 53: return 14;
    default: return UNKNOWN;
  }
}

int halogenType(int z) {  // 0 when z is not one of the four halogens
  switch (z) {
    case 9: return 11;
    case 17: return 12;
    case 35: return 13;
    case 53: return 14;
    default: return 0;
  }
}

bool atomInRing(const Atom *atom) {
  return atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx()) > 0;
}

int fallbackType(const Atom *atom) {
  const int z = atom->getAtomicNum();
  if (z != 5) {
    return fallbackByElement(z);
  }
  std::vector<int> neighbours;
  for (const auto nbr : atom->getOwningMol().atomNeighbors(atom)) {
    neighbours.push_back(nbr->getAtomicNum());
  }
  if (atom->getFormalCharge() < 0 || atom->getDegree() >= 4) {
    return B_TETRA;
  }
  for (int n : neighbours) {
    if (halogenType(n)) {
      return B_SP2_HALOGEN;
    }
  }
  if (atomInRing(atom)) {
    return B_RING;
  }
  for (const auto &[element, code] :
       {std::pair<int, int>{16, B_SP2_S}, {7, B_SP2_N}, {8, B_SP2_O}}) {
    if (std::find(neighbours.begin(), neighbours.end(), element) != neighbours.end()) {
      return code;
    }
  }
  return B_SP2_C;
}

// _atom_types(work). Constructing MMFFMolProperties mutates `work`, whether or
// not typing succeeds, exactly as AllChem.MMFFGetMoleculeProperties does.
std::vector<int> atomTypes(RWMol &work) {
  MMFF::MMFFMolProperties props(work, MMFF_VARIANT);
  std::vector<int> types;
  types.reserve(work.getNumAtoms());
  if (props.isValid()) {
    for (const auto atom : work.atoms()) {
      types.push_back(static_cast<int>(props.getMMFFAtomType(atom->getIdx())));
    }
    return types;
  }
  // A lone hydrogen halide: two atoms, one bond.
  if (work.getNumAtoms() == 2 && work.getNumBonds() == 1) {
    const int z0 = work.getAtomWithIdx(0)->getAtomicNum();
    const int z1 = work.getAtomWithIdx(1)->getAtomicNum();
    if (z0 == 1 || z1 == 1) {
      const int hIdx = (z0 == 1) ? 0 : 1;  // numbers.index(1): the first hydrogen
      const int other = 1 - hIdx;
      const int otherZ = other == 0 ? z0 : z1;
      if (halogenType(otherZ)) {
        types.assign(2, 0);
        types[hIdx] = H_TYPE;
        types[other] = halogenType(otherZ);
        return types;
      }
    }
  }
  for (const auto atom : work.atoms()) {
    types.push_back(fallbackType(atom));
  }
  return types;
}

void nodes(RWMol &work, Graph &g) {
  const auto types = atomTypes(work);
  const auto table = PeriodicTable::getTable();
  g.nNodes = static_cast<int>(work.getNumAtoms());
  g.x.assign(static_cast<size_t>(g.nNodes) * NODE_DIM, 0.0f);
  for (const auto atom : work.atoms()) {
    float *row = &g.x[static_cast<size_t>(atom->getIdx()) * NODE_DIM];
    const int mmffType = types[atom->getIdx()];
    row[std::min(std::max(mmffType, 0), MMFF_TYPE_BINS - 1)] = 1.0f;
    const int z = atom->getAtomicNum();
    const double numeric[13] = {
        static_cast<double>(z),
        atom->getMass(),
        table->getRcovalent(z),
        table->getRvdw(z),
        static_cast<double>(atom->getDegree()),
        static_cast<double>(atom->getTotalDegree()),
        static_cast<double>(atom->getTotalValence()),
        static_cast<double>(atom->getFormalCharge()),
        static_cast<double>(atom->getTotalNumHs()),
        static_cast<double>(atom->getIsAromatic()),
        static_cast<double>(atomInRing(atom)),
        static_cast<double>(static_cast<int>(atom->getHybridization())),
        static_cast<double>(mmffType),
    };
    for (int k = 0; k < 13; ++k) {
      row[MMFF_TYPE_BINS + k] = static_cast<float>(numeric[k]);
    }
  }
}

void edges(const RWMol &work, Graph &g) {
  const auto ri = work.getRingInfo();
  for (const auto bond : work.bonds()) {
    const int i = static_cast<int>(bond->getBeginAtomIdx());
    const int j = static_cast<int>(bond->getEndAtomIdx());
    const double row[EDGE_DIM] = {
        bond->getBondTypeAsDouble(),
        static_cast<double>(bond->getIsAromatic()),
        static_cast<double>(bond->getIsConjugated()),
        static_cast<double>(ri->numBondRings(bond->getIdx()) > 0),
        0.0,
        static_cast<double>(work.getAtomWithIdx(i)->getAtomicNum() == 1 ||
                            work.getAtomWithIdx(j)->getAtomicNum() == 1),
    };
    g.src.push_back(i);
    g.dst.push_back(j);
    g.src.push_back(j);
    g.dst.push_back(i);
    for (int twice = 0; twice < 2; ++twice) {
      for (unsigned int k = 0; k < EDGE_DIM; ++k) {
        g.edgeAttr.push_back(static_cast<float>(row[k]));
      }
    }
  }
  g.nEdges = static_cast<int>(g.src.size());
  // position = {(src, dst): e}; the last duplicate wins, as with a dict.
  std::map<std::pair<int, int>, int> position;
  for (int e = 0; e < g.nEdges; ++e) {
    position[{g.src[e], g.dst[e]}] = e;
  }
  g.rev.resize(g.nEdges);
  for (int e = 0; e < g.nEdges; ++e) {
    const auto it = position.find({g.dst[e], g.src[e]});
    g.rev[e] = it == position.end() ? -1 : it->second;
  }
}

// Lipinski.NumRotatableBonds as RDKit 2025.09.x computes it, the release the weights
// were trained with. RDKit 2026 (#9096, NumRotatableBondsVersion 3.1.0 -> 3.2.0) added
// &!$([CH3]) to both ends of the strict pattern. On a hydrogen-explicit molecule that
// stops counting every methyl rotor, which moves this column on ~77% of molecules
// (the hydrogen-suppressed 217-descriptor block is unaffected). Counting follows
// Lipinski.cpp's ss_matcher: uniquified SubstructMatch on a copy of the query,
// because recursive queries are not thread safe.
unsigned int numRotatableBonds2025(const ROMol &mol) {
  static const std::unique_ptr<const ROMol> pattern(SmartsToMol(
      "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])("
      "[CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]="
      "[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-,:;!@[!$"
      "(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])(["
      "CH3])[CH3])]"));
  const ROMol query(*pattern, true);
  std::vector<MatchVectType> matches;
  SubstructMatch(mol, query, matches);
  return static_cast<unsigned int>(matches.size());
}

// Lipinski.NumHAcceptors as RDKit 2025.09.x computes it: RDKit 2026 (#9060,
// version 2.0.1 -> 2.0.2) narrowed $([nH0,o,s;+0]) to $([nH0X2,o,s;+0]), so a
// degree-three aromatic N (pyrrole-type, pyridone) no longer counts. Same counting
// as Lipinski.cpp's SMARTSCOUNTFUNC.
unsigned int numHBA2025(const ROMol &mol) {
  static const std::unique_ptr<const ROMol> pattern(SmartsToMol(
      "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$("
      "[N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0,o,s;+0])]"));
  const ROMol query(*pattern, true);
  std::vector<MatchVectType> matches;
  SubstructMatch(mol, query, matches);
  return static_cast<unsigned int>(matches.size());
}

// Index of NumHAcceptors in the 217 (Descriptors._descList order).
int numHAcceptorsIndex() {
  static const int idx = [] {
    const auto names = Descriptors::Osmordred::getRDKit217DescriptorNames();
    const auto it = std::find(names.begin(), names.end(), "NumHAcceptors");
    return it == names.end() ? -1 : static_cast<int>(it - names.begin());
  }();
  return idx;
}

// _whole_molecule(work): the 11 columns on the hydrogen-explicit molecule.
void wholeMolecule(const RWMol &work, Definitions defs, std::vector<float> &out) {
  const bool legacy = defs == Definitions::RDKit2025;
  int heavy = 0, charge = 0;
  for (const auto atom : work.atoms()) {
    heavy += atom->getAtomicNum() > 1;
    charge += atom->getFormalCharge();
  }
  double logp = 0.0, mr = 0.0;
  Descriptors::calcCrippenDescriptors(work, logp, mr, true, false);
  const double v[WHOLE_MOL_DIM] = {
      static_cast<double>(work.getNumAtoms()),
      static_cast<double>(heavy),
      static_cast<double>(work.getNumBonds()),
      Descriptors::calcAMW(work, false),
      logp,
      Descriptors::calcTPSA(work, false, false),
      static_cast<double>(Descriptors::calcNumRings(work)),
      static_cast<double>(legacy ? numRotatableBonds2025(work)
                                 : Descriptors::calcNumRotatableBonds(work, Descriptors::Strict)),
      static_cast<double>(legacy ? numHBA2025(work) : Descriptors::calcNumHBA(work)),
      static_cast<double>(Descriptors::calcNumHBD(work)),
      static_cast<double>(charge),
  };
  for (double d : v) {
    out.push_back(static_cast<float>(d));
  }
}

// _descriptors(mol): float64 -> float32, non-finite -> NaN, float32 inf -> NaN,
// NaN -> 0, then np.arcsinh on the float32 array (float32 arithmetic).
void descriptors(const ROMol &mol, Definitions defs, std::vector<float> &out) {
  const auto raw = rdkit217(mol, defs);
  for (unsigned int i = 0; i < DESCRIPTOR_DIM; ++i) {
    const double v = i < raw.size() ? raw[i] : std::numeric_limits<double>::quiet_NaN();
    float f = std::isfinite(v) ? static_cast<float>(v)
                               : std::numeric_limits<float>::quiet_NaN();
    if (std::isinf(f) || std::isnan(f)) {
      f = 0.0f;
    }
    out.push_back(std::asinh(f));
  }
}

template <typename T>
void appendRaw(std::string &buf, const T *data, size_t n) {
  buf.append(reinterpret_cast<const char *>(data), n * sizeof(T));
}

}  // namespace

std::vector<double> rdkit217(const ROMol &mol, Definitions defs) {
  auto values = Descriptors::Osmordred::extractRDKitDescriptors(mol);
  const int k = numHAcceptorsIndex();
  if (defs == Definitions::RDKit2025 && k >= 0 && k < static_cast<int>(values.size())) {
    values[k] = static_cast<double>(numHBA2025(mol));
  }
  return values;
}

std::vector<float> molBlock(const ROMol &mol, Definitions defs) {
  RWMol work(mol);
  MolOps::addHs(work);
  std::vector<int> unused = atomTypes(work);  // the MMFF mutation happens first
  (void)unused;
  std::vector<float> block;
  block.reserve(MOL_DIM);
  wholeMolecule(work, defs, block);
  descriptors(mol, defs, block);
  return block;
}

Graph featurize(const ROMol &mol, Definitions defs) {
  Graph g;
  RWMol work(mol);
  MolOps::addHs(work);
  nodes(work, g);
  edges(work, g);
  g.mol.reserve(MOL_DIM);
  wholeMolecule(work, defs, g.mol);
  descriptors(mol, defs, g.mol);
  return g;
}

std::string graphBytes(const Graph &g) {
  static_assert(sizeof(float) == 4 && sizeof(std::int32_t) == 4, "layout");
  std::string buf;
  const std::int32_t header[5] = {g.nNodes, g.nEdges, static_cast<std::int32_t>(NODE_DIM),
                                  static_cast<std::int32_t>(EDGE_DIM),
                                  static_cast<std::int32_t>(g.mol.size())};
  appendRaw(buf, header, 5);
  appendRaw(buf, g.x.data(), g.x.size());
  appendRaw(buf, g.edgeAttr.data(), g.edgeAttr.size());
  appendRaw(buf, g.src.data(), g.src.size());
  appendRaw(buf, g.dst.data(), g.dst.size());
  appendRaw(buf, g.rev.data(), g.rev.size());
  appendRaw(buf, g.mol.data(), g.mol.size());
  return buf;
}

}  // namespace DmpnnFeaturizer
}  // namespace RDKit
