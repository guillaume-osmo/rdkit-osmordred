// Osmordred Parallel Inference Support
// Uses RDKit threading for parallel tree prediction
#pragma once

#include <vector>
#include <future>
#include <RDGeneral/RDThreads.h>
#include <GraphMol/ROMol.h>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// Parallel tree inference utilities

// Predict sum of trees in range [start, end) in parallel
template<typename TreePredictFunc>
inline double predictTreesParallel(
    TreePredictFunc predictTree,
    int nTrees,
    const std::vector<double>& features,
    double baseScore,
    double learningRate,
    int numThreads = 0
) {
    unsigned int nThreads = getNumThreadsToUse(numThreads);
    
    // For small tree counts, sequential is faster
    if (nThreads <= 1 || nTrees < 50) {
        double pred = baseScore;
        for (int t = 0; t < nTrees; ++t) {
            pred += learningRate * predictTree(t, features);
        }
        return pred;
    }
    
    // Parallel sum of tree predictions
    std::vector<std::future<double>> futures;
    futures.reserve(nThreads);
    
    int treesPerThread = nTrees / nThreads;
    int remainder = nTrees % nThreads;
    
    int start = 0;
    for (unsigned int t = 0; t < nThreads; ++t) {
        int count = treesPerThread + (t < static_cast<unsigned int>(remainder) ? 1 : 0);
        int end = start + count;
        
        futures.emplace_back(std::async(std::launch::async,
            [&predictTree, &features, start, end]() {
                double sum = 0.0;
                for (int i = start; i < end; ++i) {
                    sum += predictTree(i, features);
                }
                return sum;
            }));
        start = end;
    }
    
    // Collect results
    double pred = baseScore;
    for (auto& f : futures) {
        pred += learningRate * f.get();
    }
    return pred;
}


} // namespace Osmordred
} // namespace Descriptors
} // namespace RDKit
