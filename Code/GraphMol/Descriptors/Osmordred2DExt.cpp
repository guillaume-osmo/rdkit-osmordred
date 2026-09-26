//
//  Copyright (C) 2026 Osmo Labs
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//! \file Osmordred2DExt.cpp
//! Osmordred v4 — 2D Dragon/alvaDesc families absent from v3. See Osmordred2DExt.h.

#include "Osmordred2DExt.h"

#include <GraphMol/MolOps.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <Eigen/Dense>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/biconnected_components.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace RDKit {
namespace Descriptors {
namespace Osmordred2DExt {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

//! Dragon's "conventional bond order": aromatic counts as 1.5.
double conventionalBondOrder(const Bond *b) {
  switch (b->getBondType()) {
    case Bond::SINGLE: return 1.0;
    case Bond::DOUBLE: return 2.0;
    case Bond::TRIPLE: return 3.0;
    case Bond::AROMATIC: return 1.5;
    default: return 1.0;
  }
}

// --------------------------------------------------------------------------
// Atomic weighting tables, all scaled on carbon.
// --------------------------------------------------------------------------
double bondiRadius(int z) {
  static const std::unordered_map<int, double> t = {
      {1, 1.10},  {6, 1.70},  {7, 1.55},  {8, 1.52},  {9, 1.47},  {14, 2.10},
      {15, 1.80}, {16, 1.80}, {17, 1.75}, {33, 1.85}, {34, 1.90}, {35, 1.85},
      {52, 2.06}, {53, 1.98}};
  auto it = t.find(z);
  return it == t.end() ? 1.70 : it->second;
}

//! \note NOT the classic tabulated vdW volumes — those give only 137/582 bit-exact.
//! The Bondi radius cube ratio gives 576/582.
double wVolume(const Atom *a) {
  const double r = bondiRadius(a->getAtomicNum()) / bondiRadius(6);
  return r * r * r;
}
double wMass(const Atom *a) { return a->getMass() / 12.011; }
double wSanderson(const Atom *a) {
  static const std::unordered_map<int, double> t = {
      {1, 2.59},  {6, 2.75},  {7, 3.19},  {8, 3.65},  {9, 4.00},  {14, 2.14},
      {15, 2.52}, {16, 2.96}, {17, 3.48}, {35, 3.22}, {53, 2.78}};
  auto it = t.find(a->getAtomicNum());
  return (it == t.end() ? 2.75 : it->second) / 2.75;
}
double wPolarizability(const Atom *a) {
  static const std::unordered_map<int, double> t = {
      {1, 0.667}, {6, 1.76},  {7, 1.10},  {8, 0.802}, {9, 0.557}, {14, 5.38},
      {15, 3.63}, {16, 2.90}, {17, 2.18}, {35, 3.05}, {53, 5.35}};
  auto it = t.find(a->getAtomicNum());
  return (it == t.end() ? 1.76 : it->second) / 1.76;
}
double wIonization(const Atom *a) {
  static const std::unordered_map<int, double> t = {
      {1, 13.598}, {6, 11.260},  {7, 14.534}, {8, 13.618}, {9, 17.423},
      {14, 8.152}, {15, 10.487}, {16, 10.360}, {17, 12.968}, {35, 11.814},
      {53, 10.451}};
  auto it = t.find(a->getAtomicNum());
  return (it == t.end() ? 11.260 : it->second) / 11.260;
}
//! Kier-Hall intrinsic state. Partially validated (298/582 bit-exact) — the form is
//! right but some element/hybridisation classes still disagree with alvaDesc.
double wIState(const Atom *a) {
  const int z = a->getAtomicNum();
  const int n = z <= 2 ? 1 : z <= 10 ? 2 : z <= 18 ? 3 : 4;
  const int dv = PeriodicTable::getTable()->getNouterElecs(z) -
                 static_cast<int>(a->getTotalNumHs());
  const int d = std::max<int>(a->getDegree(), 1);
  const double f = 2.0 / static_cast<double>(n);
  return (f * f * static_cast<double>(dv) + 1.0) / static_cast<double>(d);
}

// --------------------------------------------------------------------------
// Osmordred-parity atomic properties.
//
// 🔴 These are the tables from Code/Osmordred/Osmordred.cpp, NOT Dragon's. They are
// deliberately separate from burdenWeights() above, which is Dragon-parity: the two
// blocks answer different questions and must not be collapsed into one. Values are RAW
// (unnormalised) — normalisation is what decides whether the Burden spectrum is
// sign-split, so changing it here silently changes the extraction convention.
// --------------------------------------------------------------------------
double propLookup(const std::unordered_map<int, double> &t, int z, double fb) {
  const auto it = t.find(z);
  return it == t.end() ? fb : it->second;
}
const std::unordered_map<int, double> &osmoVdwRadius() {
  static const std::unordered_map<int, double> t{
      {1, 1.10},  {2, 1.40},  {5, 1.92},  {6, 1.70},  {7, 1.55},  {8, 1.52},
      {9, 1.47},  {14, 2.10}, {15, 1.80}, {16, 1.80}, {17, 1.75}, {35, 1.85},
      {53, 1.98}};
  return t;
}
const std::unordered_map<int, double> &osmoSandersonEN() {
  static const std::unordered_map<int, double> t{
      {1, 2.592},  {5, 2.275},  {6, 2.746},  {7, 3.194},  {8, 3.654},  {9, 4.000},
      {14, 2.138}, {15, 2.515}, {16, 2.957}, {17, 3.475}, {35, 3.219}, {53, 2.778}};
  return t;
}
const std::unordered_map<int, double> &osmoPaulingEN() {
  static const std::unordered_map<int, double> t{
      {1, 2.2},  {5, 2.04},  {6, 2.55},  {7, 3.04},  {8, 3.44},  {9, 3.98},
      {14, 1.9}, {15, 2.19}, {16, 2.58}, {17, 3.16}, {35, 2.96}, {53, 2.66}};
  return t;
}
const std::unordered_map<int, double> &osmoAllredEN() {
  static const std::unordered_map<int, double> t{
      {1, 2.20},  {5, 2.01},  {6, 2.50},  {7, 3.07},  {8, 3.50},  {9, 4.10},
      {14, 1.74}, {15, 2.06}, {16, 2.44}, {17, 2.83}, {35, 2.74}, {53, 2.21}};
  return t;
}
const std::unordered_map<int, double> &osmoPolarizability() {
  static const std::unordered_map<int, double> t{
      {1, 0.666793}, {5, 3.03},  {6, 1.67},  {7, 1.10},  {8, 0.802}, {9, 0.557},
      {14, 5.53},    {15, 3.63}, {16, 2.90}, {17, 2.18}, {35, 3.05}, {53, 5.35}};
  return t;
}
const std::unordered_map<int, double> &osmoIonisation() {
  static const std::unordered_map<int, double> t{
      {1, 13.598443}, {5, 8.29802},   {6, 11.26030},  {7, 14.5341},
      {8, 13.61805},  {9, 17.4228},   {14, 8.15168},  {15, 10.48669},
      {16, 10.36001}, {17, 12.96763}, {35, 11.8138},  {53, 10.45126}};
  return t;
}
int osmoPQN(int z) {
  if (z <= 2) return 1;
  if (z <= 10) return 2;
  if (z <= 18) return 3;
  if (z <= 36) return 4;
  if (z <= 54) return 5;
  if (z <= 86) return 6;
  return 7;
}
//! osmordred's getValenceElectrons: hydrogen is defined as 0, not 1.
double osmoValenceElectrons(const Atom *a) {
  const int z = a->getAtomicNum();
  if (z == 1) return 0.0;
  const double zv = PeriodicTable::getTable()->getNouterElecs(z) - a->getFormalCharge();
  const double zt = z - a->getFormalCharge();
  const double denom = zt - zv - 1.0;
  return denom == 0.0 ? 0.0 : (zv - a->getTotalNumHs()) / denom;
}
//! osmordred's getSigmaElectrons: count of NON-hydrogen neighbours.
double osmoSigmaElectrons(const Atom *a) {
  double n = 0;
  for (const auto nbr : a->getOwningMol().atomNeighbors(a)) {
    if (nbr->getAtomicNum() != 1) n += 1.0;
  }
  return n;
}
double osmoIntrinsicState(const Atom *a) {
  const double d = osmoSigmaElectrons(a);
  if (d == 0.0) return 0.0;
  const double n = osmoPQN(a->getAtomicNum());
  return ((2.0 / n) * (2.0 / n) * osmoValenceElectrons(a) + 1.0) / d;
}

//! The 12 property vectors, in osmordred's own order: c dv d s Z m v se pe are p i.
std::vector<std::vector<double>> osmoAtomProperties(const ROMol &mol) {
  const auto *tbl = PeriodicTable::getTable();
  const unsigned int n = mol.getNumAtoms();
  ROMol copy(mol);
  computeGasteigerCharges(copy, 12, false);
  std::vector<std::vector<double>> p(12, std::vector<double>(n, 0.0));
  for (unsigned int i = 0; i < n; ++i) {
    const Atom *a = mol.getAtomWithIdx(i);
    const Atom *ac = copy.getAtomWithIdx(i);
    const int z = a->getAtomicNum();
    double ch = 0.0, hch = 0.0;
    ac->getPropIfPresent(common_properties::_GasteigerCharge, ch);
    if (ac->getPropIfPresent(common_properties::_GasteigerHCharge, hch)) ch += hch;
    p[0][i] = std::isfinite(ch) ? ch : 0.0;
    p[1][i] = osmoValenceElectrons(a);
    p[2][i] = osmoSigmaElectrons(a);
    p[3][i] = osmoIntrinsicState(a);
    p[4][i] = static_cast<double>(z);
    p[5][i] = tbl->getAtomicWeight(z);
    const double r = propLookup(osmoVdwRadius(), z, 2.0);
    p[6][i] = (4.0 / 3.0) * M_PI * r * r * r;
    p[7][i] = propLookup(osmoSandersonEN(), z, 0.0);
    p[8][i] = propLookup(osmoPaulingEN(), z, 0.0);
    p[9][i] = propLookup(osmoAllredEN(), z, 0.0);
    p[10][i] = propLookup(osmoPolarizability(), z, 0.0);
    p[11][i] = propLookup(osmoIonisation(), z, 0.0);
  }
  return p;
}

constexpr int kMatFunctionals = 12;
//! SpMax SpMin SpDiam SpAD SpMAD EE VE1 VE2 VE3 VR1 VR2 VR3 of a symmetric matrix.
std::vector<double> matrixFunctionals(const Eigen::MatrixXd &m) {
  std::vector<double> out(kMatFunctionals, 0.0);
  if (m.rows() < 2 || !m.allFinite()) return out;
  const Eigen::MatrixXd sym = 0.5 * (m + m.transpose());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(sym);
  if (es.info() != Eigen::Success) return out;
  const Eigen::VectorXd ev = es.eigenvalues();
  const Eigen::VectorXd lead = es.eigenvectors().col(ev.size() - 1).cwiseAbs();
  const double n = static_cast<double>(ev.size());
  const double spAD = (ev.array() - ev.mean()).abs().sum();
  double ee = 0.0;
  for (int i = 0; i < ev.size(); ++i) {
    ee += std::exp(std::min(std::max(ev(i), -50.0), 50.0));   // clamp: exp overflows
  }
  const double ve1 = lead.sum();
  double vr1 = 0.0;
  for (int i = 0; i < lead.size(); ++i) {
    for (int j = 0; j < lead.size(); ++j) {
      if (lead(i) > 1e-12 && lead(j) > 1e-12) vr1 += 1.0 / std::sqrt(lead(i) * lead(j));
    }
  }
  out = {ev(ev.size() - 1), ev(0), ev(ev.size() - 1) - ev(0), spAD, spAD / n,
         std::log1p(ee), ve1, ve1 / n, ve1 > 0 ? std::log1p(n * ve1) : 0.0,
         vr1, vr1 / n, vr1 > 0 ? std::log1p(n * vr1) : 0.0};
  return out;
}

using WeightFn = double (*)(const Atom *);
const std::vector<std::pair<const char *, WeightFn>> &burdenWeights() {
  static const std::vector<std::pair<const char *, WeightFn>> w = {
      {"m", wMass}, {"v", wVolume},     {"e", wSanderson},
      {"p", wPolarizability}, {"i", wIonization}, {"s", wIState}};
  return w;
}

// --------------------------------------------------------------------------
// Augmented edge adjacency matrix. Rows/cols are bonds; off-diagonal 1 when two
// bonds share an atom; diagonal is the bond weighting.
// --------------------------------------------------------------------------
enum EdgeWeight { EW_BO = 0, EW_ED };

Eigen::MatrixXd augmentedEdgeAdjacency(const ROMol &mol, EdgeWeight w) {
  const auto nb = static_cast<Eigen::Index>(mol.getNumBonds());
  if (nb == 0) {
    return Eigen::MatrixXd::Zero(1, 1);
  }
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(nb, nb);
  for (Eigen::Index i = 0; i < nb; ++i) {
    const Bond *bi = mol.getBondWithIdx(static_cast<unsigned int>(i));
    if (w == EW_BO) {
      A(i, i) = conventionalBondOrder(bi);
    } else {
      // Dragon "edge degree": delta_i + delta_j - 2 over heavy-atom degrees.
      A(i, i) = static_cast<double>(bi->getBeginAtom()->getDegree()) +
                static_cast<double>(bi->getEndAtom()->getDegree()) - 2.0;
    }
    for (Eigen::Index j = i + 1; j < nb; ++j) {
      const Bond *bj = mol.getBondWithIdx(static_cast<unsigned int>(j));
      const bool share = bi->getBeginAtomIdx() == bj->getBeginAtomIdx() ||
                         bi->getBeginAtomIdx() == bj->getEndAtomIdx() ||
                         bi->getEndAtomIdx() == bj->getBeginAtomIdx() ||
                         bi->getEndAtomIdx() == bj->getEndAtomIdx();
      if (share) {
        A(i, j) = A(j, i) = 1.0;
      }
    }
  }
  return A;
}

//! Burden matrix on the hydrogen-INCLUDED graph.
Eigen::MatrixXd burdenMatrix(const ROMol &molH, WeightFn w) {
  const auto n = static_cast<Eigen::Index>(molH.getNumAtoms());
  Eigen::MatrixXd B = Eigen::MatrixXd::Constant(n, n, 0.001);
  for (const auto b : molH.bonds()) {
    const auto i = static_cast<Eigen::Index>(b->getBeginAtomIdx());
    const auto j = static_cast<Eigen::Index>(b->getEndAtomIdx());
    double v = std::sqrt(conventionalBondOrder(b));
    if (b->getBeginAtom()->getDegree() == 1 || b->getEndAtom()->getDegree() == 1) {
      v += 0.1;  // terminal bond augmentation: +0.1, not +0.01
    }
    B(i, j) = B(j, i) = v;
  }
  for (Eigen::Index i = 0; i < n; ++i) {
    B(i, i) = w(molH.getAtomWithIdx(static_cast<unsigned int>(i)));
  }
  return B;
}

// --------------------------------------------------------------------------
// CATS2D pharmacophore typing (SMARTS-based, self-contained).
// --------------------------------------------------------------------------
enum PPP { P_DON = 0, P_ACC, P_POS, P_NEG, P_LIP, P_N };

//! Parsed once; the lambda-initialised static makes the first call thread-safe
//! (the previous `if (pats.empty()) push_back(...)` raced when two threads
//! computed CATS2D for the first time concurrently).
const std::vector<ROMol *> &pppPatterns() {
  static const std::vector<ROMol *> pats = [] {
    std::vector<ROMol *> res;
    static const char *sm[P_N] = {
        "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),$([n;H1;+0])]",
        "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),"
        "$([N;v3;!$(N-*=[O,N,P,S])]),$([n;H0;+0])]",
        "[$([*;+]),$([N;H2,H1,H0;X3;!$(N-[!#6]);!$(N-*=[O,N,P,S])])]",
        "[$([*;-]),$([C,S,P](=O)[O;H1,-])]",
        // Lipophilic. The `c` branch and the halogen/thioether branch were added after
        // measurement: uppercase [C] is ALIPHATIC-ONLY in SMARTS, so the previous pattern
        // typed no aromatic carbon at all. Against a fresh 575-molecule alvaDesc run the
        // old pattern fired on 9.8 nonzero CATS3D bins per molecule where the reference
        // has 15.0; with aromatic carbon it reaches 12.9, and reconstruction of the
        // reference block rises from R^2 0.255 to 0.291.
        //
        // Aromatic carbon belongs in L rather than in a separate type because the
        // five-type scheme is the ORIGINAL one (Schneider et al., Angew Chem Int Ed 38:
        // 2894, 1999). The `R` (aromatic) type in every open-source CATS — Guha 2007,
        // RDKit port by Arthur 2015, inherited by iwatobipen/CATS2D, PyBioMed and
        // molfeat — is a later EXTENSION to six types. With no R there is nowhere else
        // for it to go.
        "[$([C;!$(C=[O,N,S]);!$(C#N);!$(C[O,N])]),$([c]),$([Cl,Br,I]),"
        "$([S;D2;$(S(C)(C))])]"};
    for (int i = 0; i < P_N; ++i) {
      res.push_back(SmartsToMol(sm[i]));
    }
    return res;
  }();
  return pats;
}

const char *pppShort(int t) {
  static const char *s[P_N] = {"D", "A", "P", "N", "L"};
  return s[t];
}

}  // namespace

// ===========================================================================
// public
// ===========================================================================

const unsigned int OSMORDRED_2D_EXT_NUM_DESCRIPTORS =
    2533;  // 310 + matrix2d 608 + MDE 19 + atom pairs 1596

unsigned int calcNumCircuits(const ROMol &mol) {
  const RingInfo *ri = mol.getRingInfo();
  if (!ri->isInitialized()) {
    MolOps::findSSSR(const_cast<ROMol &>(mol));
  }
  const auto &bondRings = mol.getRingInfo()->bondRings();
  const size_t nr = bondRings.size();
  if (nr == 0) {
    return 0;
  }
  // GF(2) sums over the SSSR basis; keep each result that is a single connected cycle.
  // Bond sets are bit masks (one bit per bond), so the symmetric difference is an XOR.
  const unsigned int nBonds = mol.getNumBonds();
  const size_t nWords = (nBonds + 63) / 64;
  std::vector<std::vector<std::uint64_t>> ringBits(
      nr, std::vector<std::uint64_t>(nWords, 0));
  for (size_t r = 0; r < nr; ++r) {
    for (int bidx : bondRings[r]) {
      ringBits[r][bidx / 64] ^= std::uint64_t(1) << (bidx % 64);
    }
  }
  const unsigned int nAtoms = mol.getNumAtoms();
  std::vector<int> deg(nAtoms, 0);
  std::vector<std::array<int, 2>> nbrs(nAtoms);
  std::vector<int> touched;
  std::vector<std::uint64_t> cyc(nWords);
  std::set<std::vector<std::uint64_t>> seen;
  unsigned int count = 0;
  const size_t combos = (nr < 20) ? (1u << nr) : (1u << 20);
  for (size_t mask = 1; mask < combos; ++mask) {
    std::fill(cyc.begin(), cyc.end(), 0);
    for (size_t r = 0; r < nr && r < 20; ++r) {
      if (!(mask & (1u << r))) {
        continue;
      }
      for (size_t w = 0; w < nWords; ++w) {
        cyc[w] ^= ringBits[r][w];
      }
    }
    // every vertex must have degree 2 ...
    touched.clear();
    bool allTwo = true;
    for (size_t w = 0; w < nWords && allTwo; ++w) {
      for (std::uint64_t bits = cyc[w]; bits && allTwo; bits &= bits - 1) {
        const int bidx = static_cast<int>(w * 64) + std::countr_zero(bits);
        const Bond *b = mol.getBondWithIdx(static_cast<unsigned int>(bidx));
        const int ends[2] = {static_cast<int>(b->getBeginAtomIdx()),
                             static_cast<int>(b->getEndAtomIdx())};
        for (int e = 0; e < 2; ++e) {
          const int u = ends[e];
          if (deg[u] == 0) {
            touched.push_back(u);
          }
          if (deg[u] == 2) {
            allTwo = false;  // a third bond at u
            break;
          }
          nbrs[u][deg[u]++] = ends[1 - e];
        }
      }
    }
    bool isCycle = allTwo && !touched.empty();
    for (int u : touched) {
      if (deg[u] != 2) {
        isCycle = false;
      }
    }
    // ... AND the edge set must be connected, or it is two disjoint cycles, not one.
    if (isCycle) {
      size_t visited = 1;
      int prev = touched.front(), cur = nbrs[prev][0];
      while (cur != touched.front()) {
        const int next = nbrs[cur][0] == prev ? nbrs[cur][1] : nbrs[cur][0];
        prev = cur;
        cur = next;
        ++visited;
      }
      isCycle = visited == touched.size();
    }
    for (int u : touched) {
      deg[u] = 0;
    }
    if (isCycle && seen.insert(cyc).second) {
      ++count;
    }
  }
  return count;
}

std::vector<double> calcNarumi(const ROMol &mol) {
  const unsigned int n = mol.getNumAtoms();
  if (n == 0) {
    return {kNaN, kNaN, kNaN};
  }
  double logSum = 0.0, invSum = 0.0;
  bool zeroDegree = false;
  for (const auto a : mol.atoms()) {
    const unsigned int d = a->getDegree();
    if (d == 0) {
      zeroDegree = true;
      continue;
    }
    logSum += std::log(static_cast<double>(d));
    invSum += 1.0 / static_cast<double>(d);
  }
  const double gnar = zeroDegree ? kNaN : std::exp(logSum / static_cast<double>(n));
  const double hnar = (invSum > 0.0) ? static_cast<double>(n) / invSum : kNaN;
  return {gnar, zeroDegree ? kNaN : logSum, hnar};
}

std::vector<double> calcEdgeAdjacency(const ROMol &mol) {
  std::vector<double> out;
  out.reserve(60);
  for (EdgeWeight w : {EW_BO, EW_ED}) {
    const Eigen::MatrixXd A = augmentedEdgeAdjacency(mol, w);
    // SM0k = ln(trace(A^k) + 1). Computed from the spectrum, which is numerically
    // steadier than repeated matrix products at k = 15.
    //
    // The +1 is NOT cosmetic and was missing here. Measured against fresh Dragon 6 and
    // alvaDesc runs over 575 molecules, both agreeing exactly: at SM02_AEA(bo) plain
    // ln(trace) is bit-exact on 0.0% of molecules and ln(trace+1) on 99.1%; at
    // SM02_AEA(ed) it is 4.2% against 100.0%. The offset is 1/trace, so it is 1/15 at
    // k=2 and 1/23759 by k=8 — the two forms converge to ~96% agreement at SM08, which
    // is why a high-order spot check passes and the low orders stay wrong. Reading the
    // residual back through exp() gives exactly 1 at every order.
    // only eigenvalues are used; skipping the eigenvectors leaves them unchanged
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A, Eigen::EigenvaluesOnly);
    const Eigen::VectorXd ev =
        (es.info() == Eigen::Success) ? es.eigenvalues() : Eigen::VectorXd::Zero(A.rows());
    for (int k = 1; k <= 15; ++k) {
      double tr = 0.0;
      for (Eigen::Index i = 0; i < ev.size(); ++i) {
        tr += std::pow(ev(i), static_cast<double>(k));
      }
      out.push_back(tr > -1.0 ? std::log1p(tr) : 0.0);
    }
    for (int k = 1; k <= 15; ++k) {
      const Eigen::Index idx = ev.size() - k;
      out.push_back(idx >= 0 ? ev(idx) : 0.0);
    }
  }
  return out;
}

std::vector<double> calcBurdenEigenvalues(const ROMol &mol) {
  std::vector<double> out;
  out.reserve(96);
  const std::unique_ptr<ROMol> molH(MolOps::addHs(mol));
  for (const auto &wp : burdenWeights()) {
    const Eigen::MatrixXd B = burdenMatrix(*molH, wp.second);
    // only eigenvalues are used; skipping the eigenvectors leaves them unchanged
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(B, Eigen::EigenvaluesOnly);
    if (es.info() != Eigen::Success) {
      out.insert(out.end(), 16, kNaN);
      continue;
    }
    const Eigen::VectorXd ev = es.eigenvalues();  // ascending
    for (int k = 1; k <= 8; ++k) {                // SpMax_k = k-th largest
      const Eigen::Index idx = ev.size() - k;
      out.push_back(idx >= 0 ? ev(idx) : kNaN);
    }
    // SpMin_k = |k-th most negative|, over NEGATIVE eigenvalues only. Using the
    // k-th smallest overall instead decays agreement to 0.693 by k=8.
    std::vector<double> neg;
    for (Eigen::Index i = 0; i < ev.size(); ++i) {
      if (ev(i) < 0.0) {
        neg.push_back(ev(i));
      }
    }
    std::sort(neg.begin(), neg.end());  // most negative first
    for (int k = 1; k <= 8; ++k) {
      out.push_back(neg.size() >= static_cast<size_t>(k)
                        ? std::fabs(neg[static_cast<size_t>(k) - 1])
                        : kNaN);
    }
  }
  return out;
}

std::vector<double> calcCATS2D(const ROMol &mol) {
  const unsigned int n = mol.getNumAtoms();
  std::vector<std::vector<char>> lab(P_N, std::vector<char>(n, 0));
  const auto &pats = pppPatterns();
  for (int t = 0; t < P_N; ++t) {
    if (!pats[t]) {
      continue;
    }
    std::vector<MatchVectType> hits;
    SubstructMatch(mol, *pats[t], hits, true);
    for (const auto &h : hits) {
      for (const auto &p : h) {
        if (p.second >= 0 && static_cast<unsigned int>(p.second) < n) {
          lab[t][p.second] = 1;
        }
      }
    }
  }
  const double *D = MolOps::getDistanceMat(mol, false, false, false);
  std::vector<double> out;
  out.reserve(150);
  for (int a = 0; a < P_N; ++a) {
    for (int b = a; b < P_N; ++b) {
      std::vector<double> cnt(10, 0.0);
      double na = 0.0, nb = 0.0;
      for (unsigned int i = 0; i < n; ++i) {
        na += lab[a][i];
        nb += lab[b][i];
      }
      for (unsigned int i = 0; i < n; ++i) {
        if (!lab[a][i]) continue;
        for (unsigned int j = (a == b ? i + 1 : 0); j < n; ++j) {
          if (!lab[b][j] || i == j) continue;
          const int k = static_cast<int>(std::lround(D[i * n + j]));
          if (k >= 0 && k < 10) cnt[k] += 1.0;
        }
      }
      const double scale = std::max(na + nb, 1.0);
      for (int k = 0; k < 10; ++k) out.push_back(cnt[k] / scale);
    }
  }
  return out;
}

// --------------------------------------------------------------------------
// 2D matrix-based descriptors (604).
//
// Reimplemented from the primary literature. The family is the set of graph-invariant
// functionals of a weighted molecular matrix catalogued by Todeschini & Consonni
// (Handbook of Molecular Descriptors; Consonni & Todeschini, "Multivariate Analysis of
// Molecular Descriptors", Statistical Modelling of Molecular Descriptors in QSAR/QSPR,
// vol. 2, 2012). Individual descriptors trace to their own sources: the Wiener index
// (Wiener 1947), the Randic connectivity and hyper-Wiener indices (Randic), the Balaban
// J index, the Estrada index, the Mohar indices and algebraic connectivity of the
// Laplacian, the Barysz weighted-distance matrix, and the Burden matrix (Burden, F. R.,
// "Molecular identification number for substructure searches", J. Chem. Inf. Model.
// 29(3), 225-227, 1989).
//
// Supersedes the 144-value approximation this file carried before, which was never a
// faithful implementation: it called the reciprocal-distance matrix "Dt" when Dt denotes
// the DETOUR matrix (reciprocal-squared distance is H2), emitted an SpMin_* family and
// Dz(se|pe|are|s) weightings that the published catalogue does not define, and used
// ad-hoc VE3/VR1/VR3 formulas.
//
// Several conventions in the published definitions are ambiguous in prose and were
// settled empirically, by checking which reading reproduces published descriptor values
// on a 4763-molecule set. They are recorded here so they are not "simplified" away:
//   * "log" and "ln" are distinct: EE/SM*/Ho/HyWi are natural logarithms;
//     VE3/VE3sign/VR3/SpPosLog are base-10.
//   * l, the "last eigenvector", belongs to the largest-MAGNITUDE negative eigenvalue,
//     i.e. the last under a descending sort. It must carry NO sign test: the Laplacian
//     is positive semi-definite, so an ev<0 guard fires on floating-point noise and made
//     VE*/VR*_L a coin flip (right on 50.7% of molecules).
//   * HyWi is the classic hyper-Wiener, 0.5*sum over the UPPER TRIANGLE — half of what a
//     literal reading of the doubly-indexed sum gives.
//   * Barysz takes the MINIMUM-SUM path among tied SHORTEST paths.
//   * The polarizability weighting is the 1978 tabulation, not the 1994 one.
//
// Python reference implementation and its validation harness:
// src/sandbox/guillaume/osmordred_ext/matrix2d_v2.py.
// --------------------------------------------------------------------------

//! 1978 atomic polarizabilities. Distinct from osmoPolarizability(), which carries the
//! 1994 tabulation for the osmordred-parity BCUT block; the matrix-based block needs the
//! 1978 values to reproduce published descriptor values (33/67 columns agreeing exactly
//! against 6/67). Both tables are kept because they serve different descriptor families.
const std::unordered_map<int, double> &osmoPolarizability78() {
  static const std::unordered_map<int, double> t{
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
      {101, 18.2}, {102, 17.5}};
  return t;
}

constexpr int kMatFnCount = 34;
static const char *kMatFnNames[kMatFnCount] = {
      "Wi", "WiA", "AVS", "H", "Chi", "ChiA",
      "J", "HyWi", "SpAbs", "SpPos", "SpPosA", "SpPosLog",
      "SpMax", "SpMaxA", "SpDiam", "SpAD", "SpMAD", "Ho",
      "EE", "SM1", "SM2", "SM3", "SM4", "SM5",
      "SM6", "VE1", "VE2", "VE3", "VE1sign", "VE2sign",
      "VE3sign", "VR1", "VR2", "VR3"};

//! I-state for the matrix-based family: ((2/L)^2 * delta_v + 1) / delta with
//! **delta_v = Zv - h**.
//! 🔴 NOT osmoIntrinsicState(), which uses Kier-Hall's higher-row form
//! delta_v = (Zv - h)/(Z - Zv - 1). For second-row atoms the denominator is 1 and the
//! two coincide; they diverge only for S, P and the halogens. Solved from published
//! values rather than assumed -- every B(w) shares the same weight-independent
//! off-diagonal, so AVS_B(s) - AVS_B(m) isolates the diagonals and recovers the
//! reference I-state sum: the ratio form scores r = 0.998388 against it, Zv - h scores
//! r = 1.000000. osmoIntrinsicState() is left alone; it feeds the osmordred-parity
//! Burden/BCUT block, which is a different convention on purpose.
double osmoIntrinsicStateSimple(const Atom *a) {
  const double d = osmoSigmaElectrons(a);
  if (d == 0.0) return 0.0;
  const auto *tbl = PeriodicTable::getTable();
  const int z = a->getAtomicNum();
  const double L = osmoPQN(z);
  const double zv = tbl->getNouterElecs(z) - a->getFormalCharge();
  const double dv = zv - a->getTotalNumHs();
  return ((2.0 / L) * (2.0 / L) * dv + 1.0) / d;
}

//! Per-atom values of one weighting scheme (Z e i m p v s).
std::vector<double> osmoMatrixWeight(const ROMol &mol, char key) {
  const auto *tbl = PeriodicTable::getTable();
  const unsigned int n = mol.getNumAtoms();
  std::vector<double> w(n, 0.0);
  for (unsigned int i = 0; i < n; ++i) {
    const Atom *a = mol.getAtomWithIdx(i);
    const int z = a->getAtomicNum();
    switch (key) {
      case 'Z': w[i] = static_cast<double>(z); break;
      case 'm': w[i] = tbl->getAtomicWeight(z); break;
      case 's': w[i] = osmoIntrinsicStateSimple(a); break;
      case 'v': {
        const double r = propLookup(osmoVdwRadius(), z, 2.0);
        w[i] = (4.0 / 3.0) * M_PI * r * r * r;
        break;
      }
      case 'e': w[i] = propLookup(osmoSandersonEN(), z, 0.0); break;
      case 'p': w[i] = propLookup(osmoPolarizability78(), z, 0.0); break;
      default:  w[i] = propLookup(osmoIonisation(), z, 0.0); break;  // 'i'
    }
  }
  return w;
}

//! w_C -- the weight on a plain sp3 carbon, taken from ethane so that the graph-derived
//! intrinsic state resolves exactly as it does for any other carbon.
double osmoCarbonReference(char key) {
  static std::unique_ptr<ROMol> ethane(SmilesToMol("CC"));
  const auto w = osmoMatrixWeight(*ethane, key);
  return w.empty() ? 1.0 : w[0];
}

// All-pairs longest simple path WITHIN one biconnected component (a single ring system,
// hence small): naive backtracking DFS from every node, confined to the BCC so the
// exponential cost stays bounded. adj, visited and the rows of lsp are indexed by global
// atom index; lsp(s, g) receives the longest s-g path for every pair of block nodes.
void allPairsLongestPathBCC(const std::vector<int> &nodes,
                            const std::vector<std::vector<int>> &adj,
                            std::vector<char> &visited,
                            std::vector<std::vector<double>> &lsp) {
  std::vector<double> result(adj.size(), 0.0);
  std::function<void(int, double)> dfs = [&](int u, double dist) {
    visited[u] = 1;
    if (dist > result[u]) result[u] = dist;
    for (int v : adj[u]) {
      if (!visited[v]) dfs(v, dist + 1.0);
    }
    visited[u] = 0;
  };
  for (int s : nodes) {
    for (int n : nodes) result[n] = 0.0;
    dfs(s, 0.0);
    for (int g : nodes) {
      lsp[s][g] = std::max(lsp[s][g], result[g]);
      lsp[g][s] = lsp[s][g];
    }
  }
}

//! Detour matrix (all-pairs LONGEST simple path) via biconnected decomposition, ported
//! from osmordred v3's computeDetourMatrixBCC (rdkit-brian 99e253646, "fix DetourMatrix
//! hang on polycyclic molecules"), which follows Mordred's CalcDetour. Longest simple
//! path is NP-hard, so the exponential
//! work is confined to each biconnected block (one ring system, hence small) and the
//! blocks are combined along the block-cut tree, where any cross-block path must pass
//! through the single shared cut atom. Returns an empty matrix when one block is
//! genuinely intractable (circuit rank > 12 — a fullerene-like cage); callers emit NaN.
//! 🔴 The rank test is PER BLOCK: many separate rings stay cheap.
//! The per-block longest paths and the merge use dense N x N tables indexed by atom;
//! all entries are sums of unit bond lengths, hence exact.
std::vector<std::vector<double>> osmoDetourMatrix(const ROMol &mol) {
  using namespace boost;
  const int N = static_cast<int>(mol.getNumAtoms());
  if (N == 0) return {};
  if (N == 1) return {{0.0}};

  typedef adjacency_list<vecS, vecS, undirectedS, no_property,
                         property<edge_index_t, std::size_t>>
      BGraph;
  BGraph g(N);
  std::size_t eidx = 0;
  for (const auto &bond : mol.bonds()) {
    add_edge(bond->getBeginAtomIdx(), bond->getEndAtomIdx(), eidx++, g);
  }
  if (num_edges(g) == 0) return std::vector<std::vector<double>>(N, std::vector<double>(N, 0.0));

  std::vector<std::size_t> comp(num_edges(g));
  auto compMap = make_iterator_property_map(comp.begin(), get(edge_index, g));
  std::size_t nBcc = biconnected_components(g, compMap);

  std::vector<std::vector<std::pair<int, int>>> bccEdges(nBcc);
  graph_traits<BGraph>::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
    bccEdges[compMap[*ei]].emplace_back(static_cast<int>(source(*ei, g)),
                                        static_cast<int>(target(*ei, g)));
  }

  // Longest paths inside each block. Every pair of nodes of a block gets an entry in
  // blockLsp (row-major, global indices); pairs in two different blocks are never
  // read from it.
  std::vector<std::vector<int>> blockNodes;
  std::vector<std::vector<double>> blockLsp(N, std::vector<double>(N, 0.0));
  std::vector<std::vector<int>> adj(N);
  std::vector<char> visited(N, 0);
  std::vector<char> inBlock(N, 0);
  for (const auto &edgesInBcc : bccEdges) {
    if (edgesInBcc.empty()) continue;
    std::vector<int> bnodes;
    for (const auto &[a, b] : edgesInBcc) {
      for (int n : {a, b}) {
        if (!inBlock[n]) {
          inBlock[n] = 1;
          bnodes.push_back(n);
        }
      }
      adj[a].push_back(b);
      adj[b].push_back(a);
    }
    std::sort(bnodes.begin(), bnodes.end());
    const int rank =
        static_cast<int>(edgesInBcc.size()) - static_cast<int>(bnodes.size()) + 1;
    if (rank > 12) return {};  // intractable single ring system -> caller emits NaN
    allPairsLongestPathBCC(bnodes, adj, visited, blockLsp);
    for (int n : bnodes) {
      adj[n].clear();
      inBlock[n] = 0;
    }
    blockNodes.push_back(std::move(bnodes));
  }
  if (blockNodes.empty()) return std::vector<std::vector<double>>(N, std::vector<double>(N, 0.0));

  // Merge blocks along the block-cut tree (Mordred CalcDetour.merge / calc_weight): a new
  // block shares exactly one (cut) atom with the atoms merged so far; paths between an old
  // atom i and a new atom j go through that atom.
  std::vector<std::vector<double>> D(N, std::vector<double>(N, 0.0));
  std::vector<char> merged(N, 0);
  std::vector<int> nodes = blockNodes.back();
  for (int i : nodes) {
    merged[i] = 1;
    for (int j : nodes) D[i][j] = blockLsp[i][j];
  }
  blockNodes.pop_back();
  while (!blockNodes.empty()) {
    int found = -1, common = -1;
    for (int i = static_cast<int>(blockNodes.size()) - 1; i >= 0; --i) {
      int inter = -1, nInter = 0;
      for (int n : blockNodes[i])
        if (merged[n]) { inter = n; if (++nInter > 1) break; }
      if (nInter == 0) continue;
      if (nInter > 1) return {};  // block-cut property violated (shouldn't happen)
      found = i; common = inter; break;
    }
    if (found < 0) return {};  // disconnected (shouldn't happen for a valid molecule)
    std::vector<int> block = std::move(blockNodes[found]);
    blockNodes.erase(blockNodes.begin() + found);
    for (int j : block) {
      if (j == common) continue;
      for (int i : block) D[i][j] = D[j][i] = blockLsp[i][j];
      for (int i : nodes) {
        if (i == common) continue;
        D[i][j] = D[j][i] = D[i][common] + blockLsp[j][common];
      }
    }
    for (int j : block) {
      if (!merged[j]) {
        merged[j] = 1;
        nodes.push_back(j);
      }
    }
  }
  return D;
}

//! Barysz matrix. Off-diagonal is the sum over bonds on the shortest i->j path of
//! wC^2/(pi_b*w_a*w_b); diagonal is 1 - wC/w_i.
//! 🔴 Where several shortest paths tie, the one with the MINIMUM SUM is taken. That
//! tie-break is load-bearing — rings routinely give several equally short paths, and
//! taking an arbitrary one costs 116 exactly-agreeing columns. Implemented as a lexicographic
//! shortest path (minimise hop count, then bond-term sum) with a BFS-layered DP: relax
//! u->v only when v is one BFS layer beyond u, visiting nodes in hop order so every
//! predecessor is final before use.
//! For every source atom s, all atoms stably sorted by hop distance from s (the
//! visiting order of osmoBaryszMatrix's BFS-layered DP). It depends only on the
//! distance matrix, so it is computed once for all six weightings.
std::vector<std::vector<int>> osmoHopOrders(int n, const double *dm) {
  std::vector<std::vector<int>> orders(n, std::vector<int>(n));
  for (int s = 0; s < n; ++s) {
    const double *hops = dm + static_cast<size_t>(s) * n;
    std::vector<int> &order = orders[s];
    for (int i = 0; i < n; ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [hops](int a, int b) {
      const double ha = std::isfinite(hops[a]) ? hops[a] : 1e18;
      const double hb = std::isfinite(hops[b]) ? hops[b] : 1e18;
      return ha < hb;
    });
  }
  return orders;
}

Eigen::MatrixXd osmoBaryszMatrix(const ROMol &mol, char key, const double *dm,
                                 const std::vector<std::vector<int>> &hopOrders) {
  const int n = static_cast<int>(mol.getNumAtoms());
  const auto w = osmoMatrixWeight(mol, key);
  const double wc = osmoCarbonReference(key);
  std::vector<std::vector<std::pair<int, double>>> adj(n);
  for (const auto b : mol.bonds()) {
    const int a = b->getBeginAtomIdx(), c = b->getEndAtomIdx();
    const double denom = b->getBondTypeAsDouble() * w[a] * w[c];
    const double cost = std::fabs(denom) > 1e-12 ? wc * wc / denom : 0.0;
    adj[a].emplace_back(c, cost);
    adj[c].emplace_back(a, cost);
  }
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, n);
  const double kInf = std::numeric_limits<double>::infinity();
  for (int s = 0; s < n; ++s) {
    const double *hops = dm + static_cast<size_t>(s) * n;
    const std::vector<int> &order = hopOrders[s];
    std::vector<double> best(n, kInf);
    best[s] = 0.0;
    for (int u : order) {
      if (best[u] == kInf) continue;
      for (const auto &vc : adj[u]) {
        if (hops[vc.first] == hops[u] + 1.0 && best[u] + vc.second < best[vc.first]) {
          best[vc.first] = best[u] + vc.second;
        }
      }
    }
    for (int j = 0; j < n; ++j) {
      if (j != s && best[j] < kInf) M(s, j) = best[j];
    }
  }
  for (int i = 0; i < n; ++i) {
    M(i, i) = 1.0 - (std::fabs(w[i]) > 1e-12 ? wc / w[i] : 0.0);
  }
  return M;
}

//! Burden matrix in the form used by the matrix-based family (Burden 1989): diagonal
//! w_i/wC, bonded sqrt(pi_b) with a +0.1 bonus when either atom is terminal, 0.001 for
//! non-bonded pairs.
//! 🔴 NOT osmordred's Burden convention (pi_b/10 with a +0.01 terminal bonus) — the two
//! blocks deliberately differ. The intrinsic state is the one weighting that is NOT
//! carbon-scaled, since it is a graph invariant rather than an element constant.
Eigen::MatrixXd osmoBurdenMatrix2D(const ROMol &mol, char key) {
  const int n = static_cast<int>(mol.getNumAtoms());
  const auto w = osmoMatrixWeight(mol, key);
  const double wc = osmoCarbonReference(key);
  Eigen::MatrixXd M = Eigen::MatrixXd::Constant(n, n, 0.001);
  for (const auto b : mol.bonds()) {
    const int i = b->getBeginAtomIdx(), j = b->getEndAtomIdx();
    double v = std::sqrt(b->getBondTypeAsDouble());
    if (mol.getAtomWithIdx(i)->getDegree() == 1 || mol.getAtomWithIdx(j)->getDegree() == 1) {
      v += 0.1;
    }
    M(i, j) = M(j, i) = v;
  }
  for (int i = 0; i < n; ++i) {
    M(i, i) = (key == 's') ? w[i] : (std::fabs(wc) > 1e-12 ? w[i] / wc : w[i]);
  }
  return M;
}

//! The 34 matrix functionals, in the order of kMatFnNames. n is nSK, nBO the bond count
//! (H excluded), nCIC the circuit count, VS_i the i-th row sum, a_ij the adjacency,
//! lambda the eigenvalues and l the last eigenvector.
//! When \c eigenvaluesOut is given and the decomposition succeeds, it receives
//! the (ascending) eigenvalues of the symmetrised matrix.
std::array<double, kMatFnCount> osmoMatrixFunctionals(
    const Eigen::MatrixXd &M, const Eigen::MatrixXd &A, int nBO, int nCIC,
    bool isLaplace, Eigen::VectorXd *eigenvaluesOut = nullptr) {
  std::array<double, kMatFnCount> f;
  f.fill(0.0);
  const int n = static_cast<int>(M.rows());
  if (n < 2 || !M.allFinite()) return f;

  const Eigen::MatrixXd sym = 0.5 * (M + M.transpose());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(sym);
  if (es.info() != Eigen::Success) return f;
  const Eigen::VectorXd ev = es.eigenvalues();          // ascending
  if (eigenvaluesOut) {
    *eigenvaluesOut = ev;
  }
  const Eigen::MatrixXd evec = es.eigenvectors();

  double diagSum = 0.0, offSum = 0.0, offRecip = 0.0, hyper = 0.0;
  bool zeroDiag = true;
  for (int i = 0; i < n; ++i) {
    diagSum += M(i, i);
    if (std::fabs(M(i, i)) > 1e-12) zeroDiag = false;
    for (int j = i + 1; j < n; ++j) {
      offSum += M(i, j);
      if (std::fabs(M(i, j)) > 1e-12) offRecip += 1.0 / M(i, j);
      hyper += M(i, j) * M(i, j) + M(i, j);
    }
  }
  const Eigen::VectorXd vs = M.rowwise().sum();
  double chi = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (A(i, j) > 0) {
        const double prod = vs(i) * vs(j);
        if (prod > 0) chi += 1.0 / std::sqrt(prod);
      }
    }
  }
  const double wi = diagSum + offSum;
  f[0] = wi;                                                          // Wi
  f[1] = zeroDiag ? 2.0 * wi / (n * (n - 1.0)) : 2.0 * wi / (n * (n + 1.0));  // WiA
  f[2] = M.sum() / n;                                                 // AVS
  f[3] = diagSum + offRecip;                                          // H
  f[4] = chi;                                                         // Chi
  f[5] = nBO ? chi / nBO : 0.0;                                       // ChiA
  f[6] = (nBO / (nCIC + 1.0)) * chi;                                  // J
  f[7] = std::log1p(0.5 * hyper);                                     // HyWi

  double spAbs = 0.0, pos = 0.0;
  for (int i = 0; i < n; ++i) {
    spAbs += std::fabs(ev(i));
    if (ev(i) > 0) pos += ev(i);
  }
  f[8] = spAbs;                                                       // SpAbs
  f[9] = pos;                                                         // SpPos
  f[10] = pos / n;                                                    // SpPosA
  f[11] = pos > 0 ? (n / 10.0) * std::log10(pos) : 0.0;               // SpPosLog
  f[12] = ev(n - 1);                                                  // SpMax
  f[13] = ev(n - 1) / n;                                              // SpMaxA
  f[14] = ev(n - 1) - ((isLaplace && n > 1) ? ev(1) : ev(0));         // SpDiam
  const double mean = ev.mean();
  double spAD = 0.0;
  for (int i = 0; i < n; ++i) spAD += std::fabs(ev(i) - mean);
  f[15] = spAD;                                                       // SpAD
  f[16] = spAD / n;                                                   // SpMAD

  // Ho: characteristic-polynomial coefficients, expanded from the eigenvalues.
  std::vector<double> c(n + 1, 0.0);
  c[0] = 1.0;
  for (int i = 0; i < n; ++i) {
    for (int k = i + 1; k >= 1; --k) c[k] -= ev(i) * c[k - 1];
  }
  double coefSum = 0.0;
  for (int k = 0; k <= n; ++k) coefSum += std::fabs(c[k]);
  f[17] = std::log1p(coefSum);                                        // Ho

  double ee = 0.0;
  for (int i = 0; i < n; ++i) ee += std::exp(std::max(-700.0, std::min(700.0, ev(i))));
  f[18] = std::log1p(ee);                                             // EE

  for (int k = 1; k <= 6; ++k) {                                      // SM1..SM6
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += std::pow(ev(i), k);
    const double sign = (s > 0) - (s < 0);
    f[18 + k] = sign * std::log1p(std::fabs(s));
  }

  // l is the eigenvector of the largest-magnitude negative eigenvalue -- the last one
  // under a descending sort, i.e. column 0 here. NO sign test (see the block comment).
  // Undefined, and therefore not emitted, when that eigenvalue is not unique.
  if (n > 1 && std::fabs(ev(0) - ev(1)) < 1e-9) return f;
  double ve1 = 0.0, ve1s = 0.0;
  for (int i = 0; i < n; ++i) {
    ve1 += std::fabs(evec(i, 0));
    ve1s += evec(i, 0);
  }
  // 🔴 On a symmetric molecule the last eigenvector is antisymmetric, so sum(l) cancels
  // to floating-point noise (~1e-16) rather than to a small real value, and log10 of
  // that noise is meaningless (~-15). Clamp it: published values carry VE1sign exactly 0
  // for such molecules (466 of a 4763-molecule set), their smallest genuine nonzero is
  // 8.4e-4, and where VE1sign is 0 they report VE2sign 0 and VE3sign NaN.
  ve1s = std::fabs(ve1s);
  if (ve1s <= 1e-9) ve1s = 0.0;
  f[25] = ve1;                                                        // VE1
  f[26] = ve1 / n;                                                    // VE2
  f[27] = ve1 > 0 ? (n / 10.0) * std::log10(ve1) : 0.0;               // VE3
  f[28] = ve1s;                                                       // VE1sign
  f[29] = ve1s / n;                                                   // VE2sign
  f[30] = ve1s > 0 ? (n / 10.0) * std::log10(ve1s)
                   : std::numeric_limits<double>::quiet_NaN();        // VE3sign
  double vr1 = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (A(i, j) > 0) {
        vr1 += 1.0 / std::sqrt(1.0 + std::fabs(evec(i, 0)) * std::fabs(evec(j, 0)));
      }
    }
  }
  f[31] = vr1;                                                        // VR1
  f[32] = vr1 / n;                                                    // VR2
  f[33] = vr1 > 0 ? (n / 10.0) * std::log10(vr1) : 0.0;               // VR3
  for (int i = 0; i < kMatFnCount; ++i) {
    if (i != 30 && !std::isfinite(f[i])) f[i] = 0.0;  // 30 = VE3sign, legitimately NaN
  }
  return f;
}

//! The (matrix, functional) grid actually emitted. Not every functional is defined for
//! every matrix -- e.g. SpAbs_A is omitted because it would duplicate SpAD_A, and Wi_A
//! would duplicate nBO. Order matches the Python reference's family_names() exactly.
struct MatrixSpec {
  const char *name;
  std::vector<int> fns;
};
const std::vector<MatrixSpec> &osmoMatrixGrid() {
  static const std::vector<MatrixSpec> g{
    {"A", {6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33}},
    {"D", {0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 20, 21, 22, 23, 24}},
    {"Dt", {0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 20, 21, 22, 23, 24}},
    {"D/Dt", {0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 20, 21, 22, 23, 24}},
    {"H2", {0, 1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 20, 21, 22, 23, 24}},
    {"L", {9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 20, 21, 22, 23, 24}},
    {"X", {2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 20, 21, 22, 23, 24}},
    {"Dz(Z)", {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 19, 20, 21, 22, 23, 24}},
    {"Dz(e)", {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 19, 20, 21, 22, 23, 24}},
    {"Dz(i)", {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 19, 20, 21, 22, 23, 24}},
    {"Dz(m)", {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 19, 20, 21, 22, 23, 24}},
    {"Dz(p)", {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 19, 20, 21, 22, 23, 24}},
    {"Dz(v)", {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 19, 20, 21, 22, 23, 24}},
    {"B(e)", {0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 19, 20, 21, 22, 23, 24}},
    {"B(i)", {0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 19, 20, 21, 22, 23, 24}},
    {"B(m)", {0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 19, 20, 21, 22, 23, 24}},
    {"B(p)", {0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 19, 20, 21, 22, 23, 24}},
    {"B(s)", {0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 19, 20, 21, 22, 23, 24}},
    {"B(v)", {0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 31, 32, 33, 19, 20, 21, 22, 23, 24}},
  };
  return g;
}

namespace {
//! The heavy-atom molecule (removeHs) and its topological distance matrix
//! (MolOps::getDistanceMat defaults; null below two atoms) that Matrix2D, MDE
//! and AtomPairs2D all start from, built once by calcOsmordred2DExt.
struct HeavyAtomGraph {
  explicit HeavyAtomGraph(const ROMol &mol) : work(mol) {
    MolOps::removeHs(work);
    if (work.getNumAtoms() >= 2) {
      dm = MolOps::getDistanceMat(work);
    }
  }
  RWMol work;
  const double *dm = nullptr;
};

std::vector<double> calcMatrix2D(const HeavyAtomGraph &hg) {
  const auto &grid = osmoMatrixGrid();
  int total = 5;  // + the five Laplace-only descriptors appended below
  for (const auto &spec : grid) total += static_cast<int>(spec.fns.size());

  const ROMol &work = hg.work;
  const int n = static_cast<int>(work.getNumAtoms());
  std::vector<double> out;
  out.reserve(total);
  if (n < 2) {
    out.assign(total, 0.0);
    return out;
  }

  const double *dm = hg.dm;
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n), D(n, n), H2(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      D(i, j) = dm[static_cast<size_t>(i) * n + j];
      H2(i, j) = D(i, j) > 0 ? 1.0 / (D(i, j) * D(i, j)) : 0.0;
    }
  }
  for (const auto b : work.bonds()) {
    A(b->getBeginAtomIdx(), b->getEndAtomIdx()) = 1.0;
    A(b->getEndAtomIdx(), b->getBeginAtomIdx()) = 1.0;
  }
  const Eigen::VectorXd deg = A.rowwise().sum();
  Eigen::MatrixXd L = Eigen::MatrixXd(deg.asDiagonal()) - A;
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (A(i, j) > 0 && deg(i) > 0 && deg(j) > 0) X(i, j) = 1.0 / std::sqrt(deg(i) * deg(j));
    }
  }

  // Dt and D/Dt; empty when the detour guard trips on an intractable ring system.
  const auto detour = osmoDetourMatrix(work);
  const bool haveDetour = !detour.empty();
  Eigen::MatrixXd Dt(n, n), DDt(n, n);
  if (haveDetour) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        Dt(i, j) = detour[i][j];
        DDt(i, j) = (i != j && detour[i][j] > 0) ? D(i, j) / detour[i][j] : 0.0;
      }
    }
  }

  const int nBO = static_cast<int>(work.getNumBonds());
  // 🔴 nCIC here is the CYCLOMATIC number (nBO - nSK + components), NOT osmordred's
  // calcNumCircuits(), which counts every circuit in the GF(2) span of the SSSR basis.
  // They differ as soon as rings are fused -- naphthalene is 2 against 3 -- and that
  // propagates into J for every one of the 19 matrices as a flat 4/3 factor.
  std::vector<int> frag;
  const int nComp = static_cast<int>(MolOps::getMolFrags(work, frag));
  const int nCIC = nBO - n + nComp;

  std::unordered_map<std::string, Eigen::MatrixXd> built;
  built.emplace("A", A);
  built.emplace("D", D);
  built.emplace("H2", H2);
  built.emplace("L", L);
  built.emplace("X", X);
  if (haveDetour) {
    built.emplace("Dt", Dt);
    built.emplace("D/Dt", DDt);
  }
  static const char kBarysz[6] = {'Z', 'e', 'i', 'm', 'p', 'v'};
  const auto hopOrders = osmoHopOrders(n, dm);
  for (char key : kBarysz) {
    built.emplace(std::string("Dz(") + key + ")",
                  osmoBaryszMatrix(work, key, dm, hopOrders));
  }
  static const char kBurden[6] = {'e', 'i', 'm', 'p', 's', 'v'};
  for (char key : kBurden) {
    built.emplace(std::string("B(") + key + ")", osmoBurdenMatrix2D(work, key));
  }

  Eigen::VectorXd laplaceEigenvalues;  // filled by the "L" grid entry
  for (const auto &spec : grid) {
    auto it = built.find(spec.name);
    if (it == built.end()) {                      // detour guard tripped
      out.insert(out.end(), spec.fns.size(), std::numeric_limits<double>::quiet_NaN());
      continue;
    }
    const bool isLaplace = std::string(spec.name) == "L";
    const auto f = osmoMatrixFunctionals(it->second, A, nBO, nCIC, isLaplace,
                                         isLaplace ? &laplaceEigenvalues : nullptr);
    for (int idx : spec.fns) out.push_back(f[idx]);
  }

  // Five descriptors defined only on the Laplacian, outside the generic grid:
  //   Acon   algebraic connectivity (Fiedler value), the second-smallest eigenvalue
  //   QW_L   quasi-Wiener / Kirchhoff number, nSK * sum 1/lambda over the NONZERO
  //          eigenvalues -- the Laplacian always carries one structural zero
  //   TI1_L  first Mohar index, 2*log10(nBO/nSK) * QW_L
  //   TI2_L  second Mohar index, 4/(nSK*Acon)
  //   STN_L  spanning-tree number, ln(prod(nonzero lambda) / nSK)
  // L is exactly symmetric, so 0.5 * (L + L^T) == L bit for bit and the grid's
  // decomposition of the "L" entry already produced these eigenvalues.
  if (laplaceEigenvalues.size() != n) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> les(L);
    laplaceEigenvalues = les.eigenvalues();
  }
  const Eigen::VectorXd &lev = laplaceEigenvalues;  // ascending; lev(0) is the zero
  const double acon = lev(1);
  double qw = 0.0, logProd = 0.0;
  for (int i = 1; i < n; ++i) {
    if (lev(i) > 1e-9) qw += 1.0 / lev(i);
    if (lev(i) > 1e-12) logProd += std::log(lev(i));
  }
  qw *= n;
  out.push_back(acon);                                                   // Acon
  out.push_back(qw);                                                     // QW_L
  out.push_back(2.0 * std::log10(static_cast<double>(nBO) / n) * qw);    // TI1_L
  out.push_back(acon > 1e-12 ? 4.0 / (n * acon) : 0.0);                  // TI2_L
  out.push_back(logProd - std::log(static_cast<double>(n)));             // STN_L
  return out;
}

// --------------------------------------------------------------------------
// Molecular Distance Edge (19) — Liu, Cao & Li, J Chem Inf Comput Sci 38:387 (1998).
// MDE(a,b) = n / g, g the GEOMETRIC mean topological distance over pairs of atoms of
// the same element with degrees a and b. The published denominator mde_km is the
// 1/(2*n_km)-th root of the distance product, so mde_km^2 IS g -- not g^2. Squaring g a
// second time held this block to |r| = 0.816; dividing by g reproduces the published
// values exactly (|r| = 1.0000, max |delta| ~1e-14 over a 4763-molecule reference set).
// --------------------------------------------------------------------------
std::vector<double> calcMDE(const HeavyAtomGraph &hg) {
  static const std::vector<std::pair<int, std::vector<int>>> kClasses{
      {6, {1, 2, 3, 4}}, {8, {1, 2}}, {7, {1, 2, 3}}};
  const ROMol &work = hg.work;
  const unsigned int n = work.getNumAtoms();
  const double *dm = hg.dm;
  std::vector<double> out;
  out.reserve(19);
  for (const auto &cl : kClasses) {
    for (size_t a = 0; a < cl.second.size(); ++a) {
      for (size_t b = a; b < cl.second.size(); ++b) {
        std::vector<unsigned int> ia, ib;
        for (const auto at : work.atoms()) {
          if (at->getAtomicNum() != cl.first) continue;
          if (static_cast<int>(at->getDegree()) == cl.second[a]) ia.push_back(at->getIdx());
          if (static_cast<int>(at->getDegree()) == cl.second[b]) ib.push_back(at->getIdx());
        }
        double logsum = 0.0;
        int cnt = 0;
        for (unsigned int x : ia) {
          for (unsigned int y : ib) {
            if (a == b ? x >= y : x == y) continue;
            const double d = dm ? dm[x * n + y] : 0.0;
            if (d > 0 && std::isfinite(d)) { logsum += std::log(d); ++cnt; }
          }
        }
        if (cnt == 0) { out.push_back(0.0); continue; }
        const double gm = std::exp(logsum / cnt);   // geometric, not arithmetic
        out.push_back(gm > 0 ? cnt / gm : 0.0);
      }
    }
  }
  return out;
}

// --------------------------------------------------------------------------
// 2D Atom Pairs (1596) — T(X..Y) heteroatom topological distance sums (36), then
// B01..B10 presence flags and F01..F10 counts over 78 element pairs.
//
// Validated exact against alvaDesc: median |r| = 1.0000 across the block.
// --------------------------------------------------------------------------
std::vector<double> calcAtomPairs2D(const HeavyAtomGraph &hg) {
  static const char *kBF[12] = {"C", "N", "O", "S", "P", "F",
                                "Cl", "Br", "I", "B", "Si", "X"};
  static const int kTZ[8] = {7, 8, 16, 15, 9, 17, 35, 53};
  constexpr int kDist = 10;
  const size_t nT = 36, nP = 78;
  std::vector<double> tOut(nT, 0.0), fOut(kDist * nP, 0.0);

  const ROMol &work = hg.work;
  const unsigned int n = work.getNumAtoms();
  if (n >= 2) {
    const double *dm = hg.dm;
    std::vector<int> cls(n, 11), tcls(n, -1);   // 11 = "X", the catch-all bucket
    for (unsigned int i = 0; i < n; ++i) {
      const std::string sym = work.getAtomWithIdx(i)->getSymbol();
      for (int k = 0; k < 12; ++k) {
        if (sym == kBF[k]) { cls[i] = k; break; }
      }
      for (int k = 0; k < 8; ++k) {
        if (work.getAtomWithIdx(i)->getAtomicNum() == kTZ[k]) tcls[i] = k;
      }
    }
    auto triIndex = [](int a, int b, int m) {
      if (a > b) std::swap(a, b);
      return a * m - (a * (a - 1)) / 2 + (b - a);
    };
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = i + 1; j < n; ++j) {
        const double d = dm[i * n + j];
        if (!std::isfinite(d)) continue;
        const int k = triIndex(cls[i], cls[j], 12);
        const int di = static_cast<int>(d);
        if (di >= 1 && di <= kDist) fOut[(di - 1) * nP + k] += 1.0;
        if (tcls[i] >= 0 && tcls[j] >= 0) tOut[triIndex(tcls[i], tcls[j], 8)] += d;
      }
    }
  }
  std::vector<double> out;
  out.reserve(nT + 2 * kDist * nP);
  out.insert(out.end(), tOut.begin(), tOut.end());
  for (double v : fOut) out.push_back(v > 0 ? 1.0 : 0.0);   // B: presence
  out.insert(out.end(), fOut.begin(), fOut.end());          // F: frequency
  return out;
}
}  // namespace

std::vector<double> calcMatrix2D(const ROMol &mol) {
  return calcMatrix2D(HeavyAtomGraph(mol));
}

std::vector<double> calcMDE(const ROMol &mol) {
  return calcMDE(HeavyAtomGraph(mol));
}

std::vector<double> calcAtomPairs2D(const ROMol &mol) {
  return calcAtomPairs2D(HeavyAtomGraph(mol));
}

std::vector<double> calcOsmordred2DExt(const ROMol &mol) {
  std::vector<double> out;
  out.reserve(OSMORDRED_2D_EXT_NUM_DESCRIPTORS);
  out.push_back(static_cast<double>(calcNumCircuits(mol)));
  auto append = [&out](const std::vector<double> &v) {
    out.insert(out.end(), v.begin(), v.end());
  };
  append(calcNarumi(mol));
  append(calcEdgeAdjacency(mol));
  append(calcBurdenEigenvalues(mol));
  append(calcCATS2D(mol));
  // Matrix2D, MDE and AtomPairs2D share one heavy-atom graph
  const HeavyAtomGraph heavy(mol);
  append(calcMatrix2D(heavy));
  append(calcMDE(heavy));
  append(calcAtomPairs2D(heavy));
  out.resize(OSMORDRED_2D_EXT_NUM_DESCRIPTORS, kNaN);
  return out;
}

std::vector<std::string> getOsmordred2DExtDescriptorNames() {
  std::vector<std::string> n;
  n.reserve(OSMORDRED_2D_EXT_NUM_DESCRIPTORS);
  n.emplace_back("nCIR");
  n.insert(n.end(), {"GNar", "SNar", "HNar"});
  for (const char *w : {"bo", "ed"}) {
    for (int k = 1; k <= 15; ++k) {
      char buf[32];
      snprintf(buf, sizeof(buf), "SM%02d_AEA(%s)", k, w);
      n.emplace_back(buf);
    }
    for (int k = 1; k <= 15; ++k) {
      char buf[32];
      snprintf(buf, sizeof(buf), "Eig%02d_AEA(%s)", k, w);
      n.emplace_back(buf);
    }
  }
  for (const auto &wp : burdenWeights()) {
    for (int k = 1; k <= 8; ++k) {
      n.push_back("SpMax" + std::to_string(k) + "_Bh(" + wp.first + ")");
    }
    for (int k = 1; k <= 8; ++k) {
      n.push_back("SpMin" + std::to_string(k) + "_Bh(" + wp.first + ")");
    }
  }
  for (int a = 0; a < P_N; ++a) {
    for (int b = a; b < P_N; ++b) {
      for (int k = 0; k < 10; ++k) {
        char buf[32];
        snprintf(buf, sizeof(buf), "CATS2D_%02d_%s%s", k, pppShort(a), pppShort(b));
        n.emplace_back(buf);
      }
    }
  }
    {
    // 🔴 v3 already ships four matrix blocks (AdjacencyMatrix 12, BaryszMatrix 104,
    // DetourMatrix 14, DistanceMatrix 12 = 142), so 30 of the names below would appear
    // TWICE in a v3 + extension descriptor vector -- with different values, because v3
    // uses its own EE/VE conventions where this block uses log10 and the most-negative
    // eigenvector. Both are kept deliberately; the extension's copies take an "ext_"
    // prefix so the combined vector has no ambiguous column.
    // (v3's Barysz names are DzZ/Dzm/... without parentheses, so Dz(Z) etc. do not clash.)
    static const std::set<std::string> kClashesWithV3{
      "SpAD_A", "SpAD_D", "SpAD_Dt", "SpDiam_A", "SpDiam_D",
      "SpDiam_Dt", "SpMAD_A", "SpMAD_D", "SpMAD_Dt", "SpMax_A",
      "SpMax_D", "SpMax_Dt", "VE1_A", "VE1_D", "VE1_Dt",
      "VE2_A", "VE2_D", "VE2_Dt", "VE3_A", "VE3_D",
      "VE3_Dt", "VR1_A", "VR1_D", "VR1_Dt", "VR2_A",
      "VR2_D", "VR2_Dt", "VR3_A", "VR3_D", "VR3_Dt"};
    for (const auto &spec : osmoMatrixGrid()) {
      for (int idx : spec.fns) {
        std::string nm = std::string(kMatFnNames[idx]) + "_" + spec.name;
        n.push_back(kClashesWithV3.count(nm) ? "ext_" + nm : nm);
      }
    }
    for (const char *e : {"Acon", "QW_L", "TI1_L", "TI2_L", "STN_L"}) n.push_back(e);
    static const char *kMdeTag[3] = {"C", "O", "N"};
    static const int kMdeDeg[3][4] = {{1, 2, 3, 4}, {1, 2, 0, 0}, {1, 2, 3, 0}};
    static const int kMdeN[3] = {4, 2, 3};
    for (int c = 0; c < 3; ++c) {
      for (int a = 0; a < kMdeN[c]; ++a) {
        for (int b = a; b < kMdeN[c]; ++b) {
          n.push_back("MDE" + std::string(kMdeTag[c]) + "-" +
                      std::to_string(kMdeDeg[c][a]) + std::to_string(kMdeDeg[c][b]));
        }
      }
    }
    static const char *kHet[8] = {"N", "O", "S", "P", "F", "Cl", "Br", "I"};
    for (int i = 0; i < 8; ++i) {
      for (int j = i; j < 8; ++j) {
        n.push_back(std::string("T(") + kHet[i] + ".." + kHet[j] + ")");
      }
    }
    static const char *kBF[12] = {"C", "N", "O", "S", "P", "F",
                                  "Cl", "Br", "I", "B", "Si", "X"};
    for (const char *pre : {"B", "F"}) {
      for (int k = 1; k <= 10; ++k) {
        for (int i = 0; i < 12; ++i) {
          for (int j = i; j < 12; ++j) {
            char buf[48];
            snprintf(buf, sizeof(buf), "%s%02d[%s-%s]", pre, k, kBF[i], kBF[j]);
            n.emplace_back(buf);
          }
        }
      }
    }
  }
  return n;
}

}  // namespace Osmordred2DExt
}  // namespace Descriptors
}  // namespace RDKit
