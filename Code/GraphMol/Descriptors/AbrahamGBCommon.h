// Common structures for Gradient Boosting models
#ifndef ABRAHAM_GB_COMMON_H
#define ABRAHAM_GB_COMMON_H

#include <vector>

namespace RDKit {
namespace Descriptors {

struct GBTreeNode {
    int feature_idx;
    double threshold;
    double value;
    int left_child;
    int right_child;
};

struct GBTree {
    std::vector<GBTreeNode> nodes;
};

} // namespace Descriptors
} // namespace RDKit

#endif // ABRAHAM_GB_COMMON_H


















