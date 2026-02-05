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
#include <RDGeneral/export.h>

#ifndef _OSMORDRED_H
#define _OSMORDRED_H

#include <GraphMol/ROMol.h>
#include <vector>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

RDKIT_DESCRIPTORS_EXPORT  bool hasOsmordredSupport();

// Control function to check if Gasteiger parameters exist for all atoms BEFORE calling computeGasteigerCharges
// Returns true if all atoms have parameters for their specific environment, false otherwise
RDKIT_DESCRIPTORS_EXPORT bool checkGasteigerParameters(const RDKit::ROMol& mol);

// Filter function to check if a molecule is too large (will cause hangs during descriptor calculation)
// Returns true if molecule has >10 rings OR >200 heavy atoms
// This prevents Osmordred.Calculate from hanging on very complex molecules
// Should be called BEFORE processing molecules in batch operations
RDKIT_DESCRIPTORS_EXPORT bool isMoleculeTooLarge(const RDKit::ROMol& mol);
  
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcABCIndex(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcAcidBase(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcAromatic(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcAtomCounts(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBalabanJ(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBertzCT(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcBondCounts(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcVertexAdjacencyInformation(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcWeight(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcWienerIndex(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcVdwVolumeABC(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcTopoPSA(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcSLogP(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcHydrogenBond(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcLogS(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcLipinskiGhose(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMcGowanVolume(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcPolarizability(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcRotatableBond(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcFragmentComplexity(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcConstitutional(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcTopologicalIndex(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDetourMatrixDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDetourMatrixDescsL(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDistMatrixDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDistMatrixDescsL(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAdjMatrixDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAdjMatrixDescsL(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcCarbonTypes(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcEccentricConnectivityIndex(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBaryszMatrixDescsL(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBaryszMatrixDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcZagrebIndex(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMoeType(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMolecularDistanceEdgeDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcEStateDescs(const RDKit::ROMol& mol, bool extended = false);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcEState_VSA(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcVSA_EState(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcWalkCounts(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcTopologicalChargeDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAllChiDescriptors(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcPathCount(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcKappaShapeIndex(const RDKit::ROMol& mol); // closer "missing" k3 path count not correct on few cases
RDKIT_DESCRIPTORS_EXPORT std::vector<int> calcRingCount(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMolecularId(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBCUTs(const RDKit::ROMol& mol); // 10x faster the
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAutoCorrelation(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcFramework(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcExtendedTopochemicalAtom(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calculateETADescriptors(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChipath(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChichain(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChicluster(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcChipathcluster(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT int calcAcidicGroupCount(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT int calcBasicGroupCount(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT int countAromaticAtoms(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT int countAromaticBonds(const ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcBEStateDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcHEStateDescs(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAlphaKappaShapeIndex(const RDKit::ROMol& mol); // closer "missing" k3 path count not correct on few cases
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAbrahams(const RDKit::ROMol& mol); // Platts, Butina, Abraham, Hersey  paper J Chem Inf Comput Sci. 1999 30/8/01;39(5):835-45
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAbrahamsV2(const RDKit::ROMol& mol); // Abraham V2: GBT for A,S; Ridge for B,E,L,V
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAbrahamsV2Features(const RDKit::ROMol& mol); // Return the 291 features used by V2
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAbrahamsV2Cascade(const RDKit::ROMol& mol); // Abraham V2 Cascade: XGBoost for A,S; Ridge for V,L,E,B with LSFER golden features

// Meta37 property calculation from features (deprecated - use calcMeta37Simple instead)
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMeta37FromFeatures(const std::vector<double>& osmo, 
                                                                     const std::vector<double>& av2,
                                                                     double MW);

RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcPol(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMR(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcFlexibility(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcODT(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcSchultz(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcRNCG_RPCG(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAZV(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASV(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSV(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAZS(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASZ(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2S(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2I(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASI(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSI(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASN(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSN(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2N(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANS(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANV(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAZN(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANZ(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANI(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSZ(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANN(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2Z(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcANMat(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAZMat(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcASMat(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDSMat(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcDN2Mat(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcFrags(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAddFeatures(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcInformationContent(const RDKit::ROMol& mol, int maxradius=5); // Inspired by 1984 Basak paper

// Aggregated fast path that calls all Osmordred descriptors in C++
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcOsmordred(const RDKit::ROMol& mol);

// Batch version: Calculate Osmordred descriptors for multiple molecules using parallel processing
// Single molecule with timeout protection (default 60 seconds)
// Returns NaN vector (3585 NaN values) if computation exceeds timeout
// This is the RECOMMENDED function for production use to prevent hanging
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcOsmordredWithTimeout(
    const RDKit::ROMol& mol, int timeout_seconds = 60);

// NEW: Osmordred with PER-FAMILY timeout (default 5 seconds per family)
// If a descriptor family times out, ONLY that family = NaN, rest continues
// This gives partial results instead of all NaN on slow molecules
// family_timeout_ms: timeout in MILLISECONDS per descriptor family (default 5000 = 5s)
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcOsmordredWithPartialTimeout(
    const RDKit::ROMol& mol, int family_timeout_ms = 5000);

// Batch from SMILES: parses each SMILES with SmilesToMol() -> NEW mol. Tautomer canonical LOST.
// Use CalcOsmordredBatchFromMols when you have mol objects (e.g. Python ToBinary after tautomer canonical).
RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>> calcOsmordredBatch(
    const std::vector<std::string>& smiles_list, int n_jobs = 0);

// Batch from mol objects (Python Mol via ToBinary/MolPickler). PRESERVES tautomer canonical.
// Same object for all 3585 feature steps. Use when mols come from tautomer-canonical pipeline.
RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>> calcOsmordredBatchFromMols(
    const std::vector<const RDKit::ROMol*>& mols, int n_jobs = 0);

// Get descriptor names in the same order as calcOsmordred returns values
RDKIT_DESCRIPTORS_EXPORT std::vector<std::string> getOsmordredDescriptorNames();

// ========================================================================
// UNIFIED PHYSICOCHEMICAL PROPERTY CALCULATOR
// Computes all 22 properties in one optimized call
// ========================================================================

// Property indices for CalcPhysChemProp result
enum PhysChemPropIndex {
    PROP_V = 0, PROP_POLARIZABILITY = 1, PROP_L = 2, PROP_E = 3, PROP_B = 4,
    PROP_DENSITY = 5, PROP_RI = 6, PROP_S = 7, PROP_A = 8,
    PROP_DD = 9, PROP_DH = 10, PROP_DP = 11,
    PROP_DELTAHVAP = 12,  // Enthalpy of Vaporization (BEFORE BP/logVP in cascade)
    PROP_BP = 13, PROP_LOGVP = 14, PROP_LOGPOW = 15, PROP_LOGWS = 16,
    PROP_DELTAHF = 17, PROP_DELTAHC = 18, PROP_MP = 19, PROP_FLASHPOINT = 20,
    PROP_LOGHENRYCC = 21, PROP_DIPOLEMOMENT = 22, PROP_LOGODT = 23,
    PROP_LOGVISCOSITY = 24,  // log10(viscosity in cP) at 25°C
    NUM_PHYSCHEM_PROPS = 25
};

// Calculate ALL physicochemical properties in one optimized call
// Returns: [V, Polarizability, L, E, B, Density, RI, S, A, dD, dH, dP,
//           deltaHvap, BP, logVP, logPow, logWS, deltaHf, deltaHc, MP, Flashpoint, logHenrycc, DipoleMoment, logODT, logViscosity]
// Note: deltaHvap is calculated BEFORE BP and logVP in the cascade order
// logViscosity is log10(viscosity in cP) at 25°C
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcPhysChemProp(const RDKit::ROMol& mol);

// calcPhysChemProp2 v3.3_full: Calculate ALL 27 cascade properties using v3.3_full models
// Returns: [V, Polarizability, L, E, B, Density, RI, S, A, Modularity, HansenTotal, 
//           dD, dH, dP, BP, deltaHvap, logVP, logPow, logWS, deltaHf, deltaHc, MP, 
//           Flashpoint, logHenrycc, DipoleMoment, logViscosity, logODT]
// Uses validated Python FeatureOrderManager.build_features_for_inference() logic
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcPhysChemProp2(const RDKit::ROMol& mol);

// Batch version: Calculate properties for multiple molecules using OpenMP parallelization
// Input: vector of SMILES strings
// Output: vector of property vectors (one per molecule)
// Uses C++ OpenMP for parallelization (n_jobs parameter sets OMP_NUM_THREADS)
RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>> calcPhysChemPropBatch(
    const std::vector<std::string>& smiles_list, int n_jobs = 0);

// Get property names
RDKIT_DESCRIPTORS_EXPORT const char* getPhysChemPropName(int index);
RDKIT_DESCRIPTORS_EXPORT std::vector<std::string> getPhysChemPropNames();

// Expose feature extraction for debugging/comparison
RDKIT_DESCRIPTORS_EXPORT std::vector<double> buildPolarizabilityV2Features(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> buildPolarizabilityFeatures(const RDKit::ROMol& mol, double V, double Density, double RI);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> buildLogODTFeatures(const RDKit::ROMol& mol, const std::vector<double>& cascade_props);

// ========================================================================
// OPTIMIZED VERSIONS - Compute base features ONCE and reuse
// ========================================================================

// Struct to hold all cached base features (computed once, reused everywhere)
struct CachedMolFeatures {
    std::vector<double> osmordred;      // 3585 features
    std::vector<double> rdkit;          // 217 features
    std::vector<double> smarts;         // 291 features
    std::vector<double> base_508;       // rdkit + smarts concatenated
    std::vector<double> abraham;        // 6 Abraham parameters
    double MW;                          // Molecular weight
    double logP;                        // logP from Crippen
    double MR;                          // Molar refractivity from Crippen
};

// Extract ALL base features once (for optimal reuse)
// This is the key optimization - call this ONCE and reuse the result
RDKIT_DESCRIPTORS_EXPORT CachedMolFeatures extractAllBaseFeatures(const RDKit::ROMol& mol);

// OPTIMIZED version of calcPhysChemProp - computes base features ONCE
// Returns same 25 properties as calcPhysChemProp but 5-7x faster
// Key optimization: osmordred (3585), rdkit (217), smarts (291) computed only ONCE
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcPhysChemPropOptimal(const RDKit::ROMol& mol);

// Batch version of optimal calculation (with 10-second per-molecule timeout)
RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>> calcPhysChemPropOptimalBatch(
    const std::vector<std::string>& smiles_list, int n_jobs = 0);

// Single molecule with TRUE timeout - returns NaN if computation exceeds timeout_seconds
// Default timeout is 60 seconds (1 minute)
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcPhysChemPropWithTimeout(
    const RDKit::ROMol& mol, int timeout_seconds = 60);

// PARTIAL TIMEOUT: Returns NaN only for properties that fail/timeout
// Other properties still get valid values even if some stages fail.
// This is the recommended function for production use.
// Stages have proportional timeouts: osmordred(50%), rdkit(30%), smarts(10%), rest(10%)
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcPhysChemPropPartialTimeout(
    const RDKit::ROMol& mol, int total_timeout_seconds = 60);

// Get ALL features at once: osmordred (3585) + rdkit (217) + physchem (25) = 3827 total
// This is the most efficient way to get everything in a single call
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAllFeaturesOnce(const RDKit::ROMol& mol);

// Single molecule with timeout protection (default 60 seconds)
// Returns NaN vector (3827 NaN values) if computation exceeds timeout
// This is the RECOMMENDED function for production use
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcAllFeaturesOnceWithTimeout(
    const RDKit::ROMol& mol, int timeout_seconds = 60);

// Batch version with 60-second per-molecule timeout
RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>> calcAllFeaturesOnceBatch(
    const std::vector<std::string>& smiles_list, int n_jobs = 0);

} // namespace Osmordred
} // namespace Descriptors
} // namespace RDKit

#endif //_OSMORDRED_H
