//
//  Copyright (C) 2026 Osmo Labs
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//! \file Osmordred2DExt.h
//! \brief Osmordred v4 — 2D Dragon/alvaDesc families that v3 does not implement.
//!
//! Companion to Osmordred3D.h. Everything here is **topological**: no conformer is
//! read, so it belongs with `calcOsmordred()` rather than `calcOsmordred3D()`.
//!
//! | family | n   | validation                                                        |
//! |--------|-----|-------------------------------------------------------------------|
//! | nCIR   |  1  | **474/474 bit-exact** vs Dragon                                    |
//! | Narumi |  3  | GNar / SNar / HNar, **473-474/474 bit-exact**                       |
//! | AEA    | 60  | SM01-15 + Eig01-15 for (bo) and (ed); **r=1.0000 at every order**  |
//! | Burden | 96  | SpMax1-8 / SpMin1-8 x 6 weightings; **r=1.000, 576/582 bit-exact** |
//! | MAT2D  | 144 | spectral functionals of A, D, 1/D + 9 Barysz-weighted matrices  |
//! | MDE    |  19 | molecular distance edge, C/O/N by degree pair                   |
//! | AP2D   |1596 | T(X..Y) sums, B01-B10 presence, F01-F10 counts                  |
//! | CATS2D | 150 | 15 pharmacophore pairs x 10 topological bins                       |
//!
//! Validated against 474 molecules of Dragon output and 582 of alvaDesc output
//! gathered from public repositories.
//!
//! ## Definitions that are easy to get wrong
//!
//! **nCIR** is *not* the SSSR ring count. It counts every **connected** simple cycle in
//! the cycle space (GF(2) sums over the SSSR basis). The connectivity test is
//! load-bearing: a degree-2 check alone accepts a *disconnected* union of two cycles, so
//! diphenyl ether scores 3 instead of 2. SSSR alone matches Dragon on only 451/474.
//!
//! **SM0k is `ln(trace(A^k) + 1)`.** Two errors are possible and both look benign. The
//! raw trace ranks perfectly (rho = 1.0) but is off by ~94 absolute units. Plain
//! `ln(trace)` — the obvious correction, and what this file did until the +1 was
//! measured — is bit-exact on **0.0%** of molecules at SM02_AEA(bo) where `ln(trace+1)`
//! is exact on **99.1%**; at SM02_AEA(ed) it is 4.2% against 100.0%. Because the offset
//! is 1/trace it decays from 1/15 at k=2 to 1/23759 at k=8, so the two forms agree to
//! ~96% by SM08 and a high-order spot check passes while SM02/SM03 stay wrong.
//! Verified against fresh Dragon 6 and alvaDesc runs over 575 molecules, agreeing
//! exactly; pushing the residual back through exp() gives 1 at every order.
//!
//! **The AEA `(ed)` diagonal is `delta_i + delta_j - 2`** over heavy-atom degrees. Every
//! alternative loses: `delta_i+delta_j` 0.995, `delta_i*delta_j` 0.978,
//! `sqrt(delta_i*delta_j)` 0.996, `(delta_i-1)(delta_j-1)` 0.973, H-included 0.73.
//!
//! **The Burden matrix** is built on the **hydrogen-INCLUDED** graph; bonded off-diagonals
//! are `sqrt(conventional bond order)` (not 0.1x bond order); terminal bonds get **+0.1**
//! (not +0.01); every other pair is 0.001; the diagonal is the atomic property scaled on
//! carbon. Getting any one of those wrong yields a plausible-looking r ~ -0.6.
//!
//! **SpMin_k is the absolute value of the k-th most NEGATIVE eigenvalue**, considering
//! negative eigenvalues only — not the k-th smallest overall. With the wrong convention
//! agreement decays from 0.999 at k=1 to 0.693 at k=8; with the right one it is 1.000 at
//! every k. NaN when the matrix has fewer than k negative eigenvalues.
//!
//! **The `(v)` weighting is `(r_Bondi / r_carbon)^3`**, not the classic tabulated van der
//! Waals volumes. The tabulated values give r = 0.998 but only 137/582 bit-exact; the
//! Bondi cube ratio gives 576/582.
//!
//! ## Deliberately absent
//!
//! `SM*_AEA(dm)` and `SM*_AEA(ri)` are **not** implemented. Both Dragon and alvaDesc emit
//! corrupt values for them — a monotonicity screen (a spectral moment can never decrease
//! in k) shows `(dm)` and `(ri)` breaking at SM02->SM03 in 75% and 83% of molecules, and
//! `(bo)` breaking at SM08->SM09 in 100%. Only `(ed)` is genuine across all orders
//! (0 decreasing steps in 8632). There is therefore no reference anywhere in either
//! product against which a dipole-moment or resonance-integral weighting could be fitted.
//!
#include <RDGeneral/export.h>
#ifndef OSMORDRED_2D_EXT_H
#define OSMORDRED_2D_EXT_H

#include <GraphMol/RDKitBase.h>

#include <string>
#include <vector>

namespace RDKit {
namespace Descriptors {
namespace Osmordred2DExt {

//! Number of descriptors returned by calcOsmordred2DExt().
RDKIT_DESCRIPTORS_EXPORT extern const unsigned int OSMORDRED_2D_EXT_NUM_DESCRIPTORS;

//! \brief Dragon `nCIR` — the number of connected simple cycles in the cycle space.
//! \note Not the SSSR count. See the header notes.
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumCircuits(const ROMol &mol);

//! Narumi topological indices: GNar (geometric), SNar (harmonic-sum), HNar (harmonic).
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcNarumi(const ROMol &mol);

//! \brief Augmented edge-adjacency descriptors: `ln(trace(A^k))` for k=1..15 then the
//! 15 largest eigenvalues, for the `(bo)` and `(ed)` weightings (60 values).
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcEdgeAdjacency(const ROMol &mol);

//! \brief Burden eigenvalues: SpMax1-8 then SpMin1-8, for weightings m, v, e, p, i, s
//! (96 values). SpMin_k is `|k-th most negative eigenvalue|`; NaN when there are fewer
//! than k negative eigenvalues.
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBurdenEigenvalues(const ROMol &mol);

//! 150 CATS2D pharmacophore-pair counts: 15 point-type pairs x 10 topological bins.
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcCATS2D(const ROMol &mol);

//! All 2533 values, in the order of getOsmordred2DExtDescriptorNames().
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcOsmordred2DExt(const ROMol &mol);

//! 144 spectral descriptors of A, D, 1/D and nine Barysz property-weighted matrices.
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMatrix2D(const ROMol &mol);
//! 19 Molecular Distance Edge descriptors (Liu, Cao & Li 1998).
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMDE(const ROMol &mol);
//! 1596 2D atom pairs: T(X..Y) distance sums, then B/F presence and counts.
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAtomPairs2D(const ROMol &mol);
//! Descriptor names, in the order calcOsmordred2DExt() returns values.
RDKIT_DESCRIPTORS_EXPORT std::vector<std::string>
getOsmordred2DExtDescriptorNames();

}  // namespace Osmordred2DExt
}  // namespace Descriptors
}  // namespace RDKit

#endif
