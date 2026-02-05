// Auto-generated: Top 1000 feature indices for logODT model
// Generated from Python feature importance analysis
// Model: XGBoost with R² = 0.6703, MAE = 0.8390

#ifndef LOGODT_TOP1000_FEATURES_H
#define LOGODT_TOP1000_FEATURES_H

#include <vector>
#include <cstddef>

namespace Osmordred {
namespace LogODTFeatures {

// Total features available: 4115
// Selected features: 1000
// NOTE: Functional group indicators are already included in RDKit/SMARTS/Osmordred
constexpr size_t TOTAL_FEATURES = 4115;
constexpr size_t SELECTED_FEATURES = 1000;

// Feature type ranges in full feature vector
constexpr size_t SMARTS_START = 0;
constexpr size_t SMARTS_END = 291;
constexpr size_t RDKIT_START = 291;
constexpr size_t RDKIT_END = 508;
constexpr size_t OSMORDRED_START = 508;
constexpr size_t OSMORDRED_END = 4093;
constexpr size_t CASCADE_START = 4093;
constexpr size_t CASCADE_END = 4115;

// Top 1000 feature indices (sorted)
extern const std::vector<size_t> TOP_1000_FEATURE_INDICES;

// Feature mask (boolean array, true for selected features)
extern const std::vector<bool> FEATURE_MASK;

// Helper function to check if a feature index is selected
inline bool isFeatureSelected(size_t idx) {
    if (idx >= TOTAL_FEATURES) return false;
    return FEATURE_MASK[idx];
}

} // namespace LogODTFeatures
} // namespace Osmordred

#endif // LOGODT_TOP1000_FEATURES_H
