#pragma once
// Auto-generated Cascade Integration
// MW is at index 6 (MolWt), MR is at index 131 (MolMR) in base_features
// DO NOT add them separately!

#include <vector>
#include <cmath>

#include "../OsmordredCascadeV.h"
#include "../OsmordredCascadeE.h"
#include "../OsmordredCascadeL.h"
#include "../OsmordredCascadeB.h"
#include "../OsmordredCascadeS.h"
#include "../OsmordredCascadeA.h"
#include "../OsmordredCascadeDensity.h"
#include "../OsmordredCascadeRI.h"

namespace Osmordred {

// Golden features: raw predictions, squares, pairwise products, MW/MR interactions
inline void computeGoldenFeatures(
    const std::vector<double>& preds,
    double MW, double MR,
    std::vector<double>& golden) {
    
    golden.clear();
    int n = preds.size();
    if (n == 0) return;
    
    double eps = 1e-6;
    
    // 1. Raw predictions
    for (int i = 0; i < n; ++i) {
        golden.push_back(preds[i]);
    }
    
    // 2. Squared predictions
    for (int i = 0; i < n; ++i) {
        golden.push_back(preds[i] * preds[i]);
    }
    
    // 3. Pairwise products
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            golden.push_back(preds[i] * preds[j]);
        }
    }
    
    // 4. MW/MR interactions
    for (int i = 0; i < n; ++i) {
        golden.push_back(preds[i] * MW);
        golden.push_back(preds[i] / (MW + eps));
        golden.push_back(preds[i] * MR);
        golden.push_back(preds[i] / (MR + eps));
    }
}

// Input: base_features (508 features: 217 RDKit + 291 SMARTS)
// MW is at index 6, MR is at index 131 - DO NOT add them again!
// Output: [V, E, L, B, S, A, Density, RI]
inline std::vector<double> predictCascade(const std::vector<double>& base_features) {
    // Extract MW and MR from base_features
    double MW = base_features[6];   // MolWt
    double MR = base_features[131]; // MolMR
    
    std::vector<double> results(8, 0.0);
    std::vector<double> preds;
    std::vector<double> golden;
    std::vector<double> features;
    
    // V - index 0
    features = base_features;
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    results[0] = CascadeV::predict(features);
    preds.push_back(results[0]);
    
    // E - index 1
    features = base_features;
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    results[1] = CascadeE::predict(features);
    preds.push_back(results[1]);
    
    // L - index 2
    features = base_features;
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    results[2] = CascadeL::predict(features);
    preds.push_back(results[2]);
    
    // B - index 3
    features = base_features;
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    results[3] = CascadeB::predict(features);
    preds.push_back(results[3]);
    
    // S - index 4
    features = base_features;
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    results[4] = CascadeS::predict(features);
    preds.push_back(results[4]);
    
    // A - index 5
    features = base_features;
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    results[5] = CascadeA::predict(features);
    preds.push_back(results[5]);
    
    // Density - index 6
    features = base_features;
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    results[6] = CascadeDensity::predict(features);
    preds.push_back(results[6]);
    
    // RI - index 7
    features = base_features;
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    results[7] = CascadeRI::predict(features);
    preds.push_back(results[7]);
    
    return results;
}

// Get Abraham parameters only (V, E, L, B, S, A)
inline std::vector<double> predictAbraham(const std::vector<double>& base_features) {
    auto full = predictCascade(base_features);
    return std::vector<double>(full.begin(), full.begin() + 6);
}

} // namespace Osmordred
