//  Copyright (c) 2025, Guillaume Godin Osmo Labs, PBC's and others
//  All rights reserved.
//
// SMARTS291 - Abraham V2 SMARTS-based Features
// 
// This module provides 291 SMARTS-based features for molecular property prediction.
// The features consist of:
//   - 241 base features: SMARTS pattern counts (sorted alphabetically)
//   - 50 golden features: Ratio features derived from base features
//
// These features are designed for Abraham parameter prediction (A, B, E, L, S, V)
// and are used as input to machine learning models for physicochemical property prediction.

#ifndef SMARTS291_H
#define SMARTS291_H

#include <RDGeneral/export.h>
#include <GraphMol/ROMol.h>
#include <vector>
#include <string>

namespace RDKit {
namespace Descriptors {
namespace SMARTS291 {

// Check if SMARTS291 support is available
RDKIT_DESCRIPTORS_EXPORT bool hasSMARTS291Support();

// Extract 241 base SMARTS features
// These are SMARTS pattern match counts, sorted alphabetically by feature name
// Returns: vector of 241 double values (count of matches for each SMARTS pattern)
RDKIT_DESCRIPTORS_EXPORT std::vector<double> extractBaseFeatures(const RDKit::ROMol& mol);

// Generate 50 golden features from base features
// These are ratio features: baseFeatures[i] / baseFeatures[j]
// Different Abraham parameters (A, B, E, L, S, V) use different golden feature definitions
// Returns: vector of 50 double values (0.0 if denominator is 0)
RDKIT_DESCRIPTORS_EXPORT std::vector<double> generateGoldenFeaturesA(const std::vector<double>& baseFeatures);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> generateGoldenFeaturesS(const std::vector<double>& baseFeatures);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> generateGoldenFeaturesB(const std::vector<double>& baseFeatures);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> generateGoldenFeaturesE(const std::vector<double>& baseFeatures);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> generateGoldenFeaturesL(const std::vector<double>& baseFeatures);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> generateGoldenFeaturesV(const std::vector<double>& baseFeatures);

// Extract all 291 SMARTS features for a given Abraham parameter
// Returns: vector of 291 double values (241 base + 50 golden)
RDKIT_DESCRIPTORS_EXPORT std::vector<double> extractSMARTS291_A(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> extractSMARTS291_S(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> extractSMARTS291_B(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> extractSMARTS291_E(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> extractSMARTS291_L(const RDKit::ROMol& mol);
RDKIT_DESCRIPTORS_EXPORT std::vector<double> extractSMARTS291_V(const RDKit::ROMol& mol);

// Batch extraction from SMILES list
RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>> extractSMARTS291Batch(
    const std::vector<std::string>& smiles_list, char param = 'A', int n_jobs = 0);

// Get feature names
RDKIT_DESCRIPTORS_EXPORT std::vector<std::string> getBaseFeatureNames();
RDKIT_DESCRIPTORS_EXPORT std::vector<std::string> getGoldenFeatureNames(char param = 'A');
RDKIT_DESCRIPTORS_EXPORT std::vector<std::string> getSMARTS291FeatureNames(char param = 'A');

}  // namespace SMARTS291
}  // namespace Descriptors
}  // namespace RDKit

#endif  // SMARTS291_H
