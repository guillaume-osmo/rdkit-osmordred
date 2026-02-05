// Auto-generated: logWS GBT model (trained on full data)
#ifndef OSMORDRED_GBT_LOGWS_H
#define OSMORDRED_GBT_LOGWS_H

#include <vector>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

struct GBTModel_LOGWS {
    static constexpr double learning_rate = 0.0300000000;
    static constexpr double init_value = -2.9023535372;
    static constexpr int n_features = 200;
    static constexpr int n_trees = 300;
    static double predict(const std::vector<double>& features);
};

}}}

#endif
