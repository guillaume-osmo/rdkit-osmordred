#pragma once

#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <vector>
#include <memory>
#include <string>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// Helper functions for cached SMARTS queries (like Osmordred pattern)
const std::vector<std::shared_ptr<RWMol>>& GetQEDAcceptorQueries();
const std::vector<std::shared_ptr<RWMol>>& GetQEDAlertQueries();
const std::vector<std::shared_ptr<RWMol>>& GetFragmentQueries();

// Extract all 217 RDKit descriptors in exact order matching Python's Descriptors._descList
// This matches the order used in Python training: desc_names = [d[0] for d in Descriptors._descList]
RDKIT_DESCRIPTORS_EXPORT std::vector<double> extractRDKitDescriptors(const ROMol& mol);

// Batch version: Extract RDKit descriptors for multiple molecules using parallel processing
// Input: vector of SMILES strings
// Output: vector of descriptor vectors (one per molecule, each with 217 descriptors)
// Uses std::async for parallelization (n_jobs parameter controls thread count)
// For small batches (< 10), uses sequential processing to avoid overhead
RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>> extractRDKitDescriptorsBatch(
    const std::vector<std::string>& smiles_list, int n_jobs = 0);

// Batch version: Extract RDKit descriptors from mol objects directly (SNN support)
// Input: vector of ROMol pointers (can contain nullptr for invalid molecules)
// Output: vector of descriptor vectors (one per molecule, each with 217 descriptors)
// Handles nullptr gracefully (returns zeros for that molecule)
// Uses std::async for parallelization (n_jobs parameter controls thread count)
RDKIT_DESCRIPTORS_EXPORT std::vector<std::vector<double>> extractRDKitDescriptorsFromMolsBatch(
    const std::vector<const ROMol*>& mols, int n_jobs = 0);

// Get descriptor names in the same order as extractRDKitDescriptors returns values
// Returns exactly 217 names matching Python's Descriptors._descList order
RDKIT_DESCRIPTORS_EXPORT std::vector<std::string> getRDKit217DescriptorNames();

}  // namespace Osmordred
}  // namespace Descriptors
}  // namespace RDKit
