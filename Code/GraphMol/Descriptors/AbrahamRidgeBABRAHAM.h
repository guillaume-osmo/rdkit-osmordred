// Auto-generated header for B_abraham Ridge regression model
// Generated from Python scikit-learn Ridge

#ifndef ABRAHAM_RIDGE_BABRAHAM_H
#define ABRAHAM_RIDGE_BABRAHAM_H

#include <vector>
#include <string>

namespace RDKit {
namespace Descriptors {

struct RidgeModel_BABRAHAM {
    static constexpr double intercept = 0.6930157997;
    static constexpr double alpha = 1.0000000000;
    static constexpr int n_features = 291;
    
    static std::vector<double> get_coefficients();
    static std::vector<double> get_mean();  // Feature means for scaling
    static std::vector<double> get_scale();  // Feature scales for scaling
    static std::vector<std::string> get_feature_names();
    
    // Predict function (features should be scaled)
    static double predict(const std::vector<double>& features_scaled);
    
    // Predict function (features NOT scaled - will scale internally)
    static double predict_unscaled(const std::vector<double>& features);
};

} // namespace Descriptors
} // namespace RDKit

#endif // ABRAHAM_RIDGE_BABRAHAM_H
