// Auto-generated: BP GBT model (trained on full data)
#ifndef OSMORDRED_GBT_BP_H
#define OSMORDRED_GBT_BP_H

#include <vector>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

struct GBTModel_BP {
    static constexpr double learning_rate = 0.0300000000;
    static constexpr double init_value = 234.0837020240;
    static constexpr int n_features = 200;
    static constexpr int n_trees = 300;
    static double predict(const std::vector<double>& features);
};

}}}

#endif
