// Auto-generated: dH GBT model (trained on full data)
#ifndef OSMORDRED_GBT_DH_H
#define OSMORDRED_GBT_DH_H

#include <vector>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

struct GBTModel_DH {
    static constexpr double learning_rate = 0.0300000000;
    static constexpr double init_value = 5.9132186712;
    static constexpr int n_features = 200;
    static constexpr int n_trees = 300;
    static double predict(const std::vector<double>& features);
};

}}}

#endif
