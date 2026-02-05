// Auto-generated header for A_abraham Gradient Boosting model
// Generated from Python scikit-learn GradientBoostingRegressor

#ifndef ABRAHAM_GB_AABRAHAM_H
#define ABRAHAM_GB_AABRAHAM_H

#include <vector>
#include <string>
#include "AbrahamGBCommon.h"

namespace RDKit {
namespace Descriptors {

struct GBModel_AABRAHAM {
    static constexpr int n_estimators = 200;
    static constexpr double learning_rate = 0.0500000000;
    static constexpr double init_mean = 0.2386943874;
    static constexpr int n_features = 291;
    
    static std::vector<GBTree> get_trees();
    static std::vector<std::string> get_feature_names();
    
    // Predict function
    static double predict(const std::vector<double>& features);
};

} // namespace Descriptors
} // namespace RDKit

#endif // ABRAHAM_GB_AABRAHAM_H
