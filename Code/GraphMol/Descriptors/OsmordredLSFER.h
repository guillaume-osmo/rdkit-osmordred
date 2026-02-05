// Auto-generated: LSFER Correction Equations for Abraham V2
// Uses V and L residuals to correct A and S predictions
// Trained on 80,284 molecules

#pragma once
#include <vector>
#include <cmath>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// ============================================================================
// LSFER CONSISTENCY EQUATIONS
// ============================================================================

// V = f(A, B, S, E, L) - R² = 0.986092
const std::vector<double> V_MEAN = {0.3261210834, 0.9526256069, 1.5070779651, 1.3076007338, 8.3217785674};
const std::vector<double> V_SCALE = {0.3388214503, 0.5887213789, 0.6863085728, 0.7760665971, 3.6862989500};
const std::vector<double> V_COEF = {-0.0359761231, -0.0373312741, -0.1154352534, -0.3319457315, 1.0502913877};
const double V_INTERCEPT = 1.7258821756;

// L = f(A, B, S, E, V) - R² = 0.993525
const std::vector<double> L_MEAN = {0.3261210834, 0.9526256069, 1.5070779651, 1.3076007338, 1.7258821756};
const std::vector<double> L_SCALE = {0.3388214503, 0.5887213789, 0.6863085728, 0.7760665971, 0.7258940098};
const std::vector<double> L_COEF = {0.1116715959, 0.1695472651, 0.4057246967, 1.1786579678, 2.4833076605};
const double L_INTERCEPT = 8.3217785674;

// A_corrected = f(B, S, E, L, V, V_residual) - R² = 0.999996
const std::vector<double> A_CORR_MEAN = {0.9526256069, 1.5070779651, 1.3076007338, 8.3217785674, 1.7258821756, -0.0000000000};
const std::vector<double> A_CORR_SCALE = {0.5887213789, 0.6863085728, 0.7760665971, 3.6862989500, 0.7258940098, 0.0856064117};
const std::vector<double> A_CORR_COEF = {-0.3503435378, -1.0843961609, -3.1188262235, 9.8680388496, -6.8203834049, 0.8043430839};
const double A_CORR_INTERCEPT = 0.3261210834;

// S_corrected = f(A, B, E, L, V, V_residual, L_residual) - R² = 1.000000
const std::vector<double> S_CORR_MEAN = {0.3261210834, 0.9526256069, 1.3076007338, 8.3217785674, 1.7258821756, -0.0000000000, 0.0000000000};
const std::vector<double> S_CORR_SCALE = {0.3388214503, 0.5887213789, 0.7760665971, 3.6862989500, 0.7258940098, 0.0856064117, 0.2966364504};
const std::vector<double> S_CORR_COEF = {-0.1777819902, -0.3148589131, -2.0010905552, 6.2276476177, -4.1474596032, -0.2231104200, -0.7214088641};
const double S_CORR_INTERCEPT = 1.5070779651;

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

// Predict V from A, B, S, E, L
inline double predictV(double A, double B, double S, double E, double L) {
    std::vector<double> input = {A, B, S, E, L};
    double result = V_INTERCEPT;
    for (size_t i = 0; i < 5; ++i) {
        double scaled = (input[i] - V_MEAN[i]) / V_SCALE[i];
        result += V_COEF[i] * scaled;
    }
    return result;
}

// Predict L from A, B, S, E, V
inline double predictL(double A, double B, double S, double E, double V) {
    std::vector<double> input = {A, B, S, E, V};
    double result = L_INTERCEPT;
    for (size_t i = 0; i < 5; ++i) {
        double scaled = (input[i] - L_MEAN[i]) / L_SCALE[i];
        result += L_COEF[i] * scaled;
    }
    return result;
}

// Compute V residual
inline double computeVResidual(double A, double B, double S, double E, double L, double V) {
    return V - predictV(A, B, S, E, L);
}

// Compute L residual
inline double computeLResidual(double A, double B, double S, double E, double L, double V) {
    return L - predictL(A, B, S, E, V);
}

// Correct A using LSFER self-consistency
// Input: B, S, E, L, V (from original Abraham), V_residual
// Returns: Corrected A value
inline double correctA(double B, double S, double E, double L, double V, double V_residual) {
    std::vector<double> input = {B, S, E, L, V, V_residual};
    double result = A_CORR_INTERCEPT;
    for (size_t i = 0; i < 6; ++i) {
        double scaled = (input[i] - A_CORR_MEAN[i]) / A_CORR_SCALE[i];
        result += A_CORR_COEF[i] * scaled;
    }
    return result;
}

// Correct S using LSFER self-consistency
// Input: A, B, E, L, V (from original Abraham), V_residual, L_residual
// Returns: Corrected S value
inline double correctS(double A, double B, double E, double L, double V, 
                       double V_residual, double L_residual) {
    std::vector<double> input = {A, B, E, L, V, V_residual, L_residual};
    double result = S_CORR_INTERCEPT;
    for (size_t i = 0; i < 7; ++i) {
        double scaled = (input[i] - S_CORR_MEAN[i]) / S_CORR_SCALE[i];
        result += S_CORR_COEF[i] * scaled;
    }
    return result;
}

// Compute corrected Abraham V2 parameters
// Returns: [A_corrected, B, S_corrected, E, L, V]
inline std::vector<double> calcAbrahamsV2Corrected(double A, double B, double S, 
                                                    double E, double L, double V) {
    // Compute residuals
    double V_residual = computeVResidual(A, B, S, E, L, V);
    double L_residual = computeLResidual(A, B, S, E, L, V);
    
    // Correct A and S
    double A_corrected = correctA(B, S, E, L, V, V_residual);
    double S_corrected = correctS(A, B, E, L, V, V_residual, L_residual);
    
    return {A_corrected, B, S_corrected, E, L, V};
}

// Compute outlier score (for quality assessment)
inline double computeOutlierScore(double A, double B, double S, double E, double L, double V) {
    double V_residual = computeVResidual(A, B, S, E, L, V);
    double L_residual = computeLResidual(A, B, S, E, L, V);
    return std::sqrt(V_residual * V_residual + L_residual * L_residual);
}

// Check if molecule is an outlier (threshold = 0.585, 95th percentile)
inline bool isOutlier(double A, double B, double S, double E, double L, double V,
                      double threshold = 0.585) {
    return computeOutlierScore(A, B, S, E, L, V) > threshold;
}

}  // namespace Osmordred
}  // namespace Descriptors
}  // namespace RDKit
