#pragma once
// Osmordred Cascade Integration
// Cascade: V → E → L → B → S → A → Density → RI
// All models use XGBoost with golden features from previous predictions

#include <vector>
#include <cmath>
#include <string>

#include "OsmordredCascadeV.h"
#include "OsmordredCascadeE.h"
#include "OsmordredCascadeL.h"
#include "OsmordredCascadeB.h"
#include "OsmordredCascadeS.h"
#include "OsmordredCascadeA.h"
#include "OsmordredCascadeDensity.h"
#include "OsmordredCascadeRI.h"


// Golden feature computation for cascade
inline void computeGoldenFeatures(
    const std::vector<double>& preds,  // Previous predictions [V, E, L, B, S, A, Density]
    double MW, double MR,
    std::vector<double>& golden) {
    
    golden.clear();
    const double eps = 0.01;
    int n = preds.size();
    
    if (n == 0) return;
    
    // 1) Raw values + squares + sqrt
    for (int i = 0; i < n; ++i) {
        golden.push_back(preds[i]);
        golden.push_back(preds[i] * preds[i]);
        golden.push_back(std::sqrt(std::max(0.0, preds[i])));
    }
    
    // 2) Pairwise products
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            golden.push_back(preds[i] * preds[j]);
        }
    }
    
    // 3) Ratios (first 3 vs all)
    for (int i = 0; i < std::min(3, n); ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                golden.push_back(preds[i] / (preds[j] + eps));
            }
        }
    }
    
    // 4) Physics: MR_A = E + 2.83195*V - 0.52553 (if V and E available)
    if (n >= 2) {
        double V = preds[0];
        double E = preds[1];
        double MR_A = E + 2.83195 * V - 0.52553;
        golden.push_back(MR_A);
        golden.push_back(MR / (MR_A + 1.0));
    }
    
    // 5) MW/MR interactions
    for (int i = 0; i < n; ++i) {
        golden.push_back(preds[i] * MW);
        golden.push_back(preds[i] / MW);
        golden.push_back(preds[i] * MR);
    }
}

// Physics features for Density and RI
inline void computePhysicsFeatures(
    const std::vector<double>& abraham,  // [V, E, L, B, S, A]
    double MW, double MR, double Density,
    std::vector<double>& physics) {
    
    physics.clear();
    const double eps = 0.01;
    
    double V = abraham[0];
    double E = abraham[1];
    
    // MR_A (Abraham's molar refractivity)
    double MR_A = E + 2.83195 * V - 0.52553;
    physics.push_back(MR_A);
    physics.push_back(MR / (MR_A + 1.0));
    
    // Vm proxy from McGowan volume
    double Vm_mcgowan = V * 100.0;  // cm³/mol
    physics.push_back(MW / Vm_mcgowan);
    
    // If Density is predicted
    if (Density > 0) {
        double Vm = MW / (Density + eps);
        physics.push_back(Vm);
        
        // Lorentz-Lorenz factor proxy
        double f_n_approx = MR * Density / MW;
        physics.push_back(f_n_approx);
        physics.push_back(f_n_approx * f_n_approx);
        physics.push_back(Density);
    }
}

namespace Osmordred {

// Full cascade prediction
// Input: base_features (217 RDKit descriptors), MW, MR
// Output: [V, E, L, B, S, A, Density, RI]
inline std::vector<double> predictCascade(
    const std::vector<double>& base_features,
    double MW, double MR) {
    
    std::vector<double> results(8, 0.0);
    std::vector<double> preds;  // Accumulates predictions
    std::vector<double> golden;
    std::vector<double> features;
    
    // V (index 0)
    features = base_features;
    features.push_back(MW);
    features.push_back(MR);
    results[0] = CascadeV::predict(features);
    preds.push_back(results[0]);
    
    // E (index 1)
    features = base_features;
    features.push_back(MW);
    features.push_back(MR);
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    results[1] = CascadeE::predict(features);
    preds.push_back(results[1]);
    
    // L (index 2)
    features = base_features;
    features.push_back(MW);
    features.push_back(MR);
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    results[2] = CascadeL::predict(features);
    preds.push_back(results[2]);
    
    // B (index 3)
    features = base_features;
    features.push_back(MW);
    features.push_back(MR);
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    results[3] = CascadeB::predict(features);
    preds.push_back(results[3]);
    
    // S (index 4)
    features = base_features;
    features.push_back(MW);
    features.push_back(MR);
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    results[4] = CascadeS::predict(features);
    preds.push_back(results[4]);
    
    // A (index 5)
    features = base_features;
    features.push_back(MW);
    features.push_back(MR);
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    results[5] = CascadeA::predict(features);
    preds.push_back(results[5]);
    
    // Density (index 6) - uses all Abraham + physics
    features = base_features;
    features.push_back(MW);
    features.push_back(MR);
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    // Add MR_A physics
    double V = results[0], E = results[1];
    double MR_A = E + 2.83195 * V - 0.52553;
    features.push_back(MR_A);
    features.push_back(MR / (MR_A + 1.0));
    features.push_back(MW / (V * 100.0));  // MW/Vm proxy
    results[6] = CascadeDensity::predict(features);
    preds.push_back(results[6]);
    
    // RI (index 7) - uses all Abraham + Density + Lorentz-Lorenz
    features = base_features;
    features.push_back(MW);
    features.push_back(MR);
    computeGoldenFeatures(preds, MW, MR, golden);
    features.insert(features.end(), golden.begin(), golden.end());
    // Physics features
    double Density = results[6];
    features.push_back(MR_A);
    features.push_back(MR / (MR_A + 1.0));
    features.push_back(Density);
    double f_n_approx = MR * Density / MW;
    features.push_back(f_n_approx);
    features.push_back(f_n_approx * f_n_approx);
    features.push_back(MW / (Density + 0.01));  // Vm
    results[7] = CascadeRI::predict(features);
    
    return results;
}

} // namespace Osmordred
