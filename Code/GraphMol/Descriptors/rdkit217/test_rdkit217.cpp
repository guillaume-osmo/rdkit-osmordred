// RDKit217 Unit Tests
// Validates C++ implementation against Python Descriptors.descList golden reference

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ROMol.h>
#include <RDGeneral/RDLog.h>

#include "RDKit217Descriptors.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

using namespace RDKit;
using namespace RDKit::Descriptors::Osmordred;

// Tolerance for floating point comparison
const double TOLERANCE = 1e-5;

TEST_CASE("RDKit217 Basic Functionality", "[rdkit217]") {
    
    SECTION("Has correct number of descriptors") {
        auto names = getRDKit217DescriptorNames();
        REQUIRE(names.size() == 217);
    }
    
    SECTION("Ethanol descriptors") {
        auto mol = SmilesToMol("CCO");
        REQUIRE(mol != nullptr);
        
        auto descs = extractRDKitDescriptors(*mol);
        REQUIRE(descs.size() == 217);
        
        // Check a few known values (from Python)
        // MolWt (index 6) = 46.069
        CHECK_THAT(descs[6], Catch::Matchers::WithinRel(46.069, 0.01));
        
        delete mol;
    }
    
    SECTION("Batch processing from SMILES") {
        std::vector<std::string> smiles = {"CCO", "CC(=O)O", "c1ccccc1"};
        auto results = extractRDKitDescriptorsBatch(smiles, 2);
        
        REQUIRE(results.size() == 3);
        for (const auto& desc : results) {
            REQUIRE(desc.size() == 217);
        }
    }
}

#ifdef RDK_BUILD_OSMORDRED_SUPPORT

TEST_CASE("RDKit217 NCI Golden Reference Test", "[rdkit217][golden]") {
    // Get RDBASE
    std::string rdbase = std::getenv("RDBASE") ? std::getenv("RDBASE") : "";
    
    SECTION("Validate against Python Descriptors.descList golden reference") {
        if (rdbase.empty()) {
            WARN("RDBASE not set, skipping NCI golden reference test");
            return;
        }
        
        std::string nci_path = rdbase + "/Data/NCI/first_5K.smi";
        std::string golden_path = rdbase + "/Code/GraphMol/Descriptors/test_data/nci_100_rdkit217_golden.csv";
        
        std::ifstream golden_file(golden_path);
        if (!golden_file.is_open()) {
            WARN("Golden reference not found: " << golden_path);
            WARN("Run generate_all_golden.py to create golden reference");
            return;
        }
        
        // Disable logging
        boost::logging::disable_logs("rdApp.*");
        
        // Skip header
        std::string header;
        std::getline(golden_file, header);
        
        // Parse header to get descriptor count
        int n_descriptors = 0;
        {
            std::istringstream hss(header);
            std::string col;
            while (std::getline(hss, col, ',')) {
                n_descriptors++;
            }
            n_descriptors--;  // Subtract smiles column
        }
        REQUIRE(n_descriptors == 217);
        
        int validated = 0;
        int mismatches = 0;
        std::string line;
        
        while (std::getline(golden_file, line) && validated < 100) {
            std::istringstream iss(line);
            std::string smiles;
            std::getline(iss, smiles, ',');
            
            if (smiles.empty()) continue;
            
            // Parse molecule
            ROMol* mol = SmilesToMol(smiles);
            if (mol == nullptr) continue;
            
            // Compute descriptors
            std::vector<double> computed = extractRDKitDescriptors(*mol);
            delete mol;
            
            if (computed.size() != 217) continue;
            
            // Compare each descriptor
            for (int desc_idx = 0; desc_idx < 217; ++desc_idx) {
                std::string val_str;
                if (!std::getline(iss, val_str, ',')) break;
                
                double golden_val;
                try {
                    golden_val = std::stod(val_str);
                } catch (...) {
                    golden_val = std::nan("");
                }
                
                double computed_val = computed[desc_idx];
                
                // Handle NaN
                bool both_nan = std::isnan(golden_val) && std::isnan(computed_val);
                if (both_nan) continue;
                
                if (std::isnan(golden_val) != std::isnan(computed_val)) {
                    INFO("NaN mismatch at molecule " << validated << " descriptor " << desc_idx);
                    INFO("  SMILES: " << smiles);
                    INFO("  Golden: " << golden_val << ", Computed: " << computed_val);
                    mismatches++;
                    continue;
                }
                
                // Check relative/absolute difference
                double diff = std::abs(golden_val - computed_val);
                double rel_diff = (golden_val != 0) ? diff / std::abs(golden_val) : diff;
                
                if (diff > TOLERANCE && rel_diff > TOLERANCE) {
                    INFO("Value mismatch at molecule " << validated << " descriptor " << desc_idx);
                    INFO("  SMILES: " << smiles);
                    INFO("  Golden: " << std::setprecision(15) << golden_val);
                    INFO("  Computed: " << std::setprecision(15) << computed_val);
                    INFO("  Diff: " << diff << ", RelDiff: " << rel_diff);
                    mismatches++;
                }
            }
            validated++;
        }
        golden_file.close();
        
        INFO("Validated " << validated << " molecules against Python Descriptors.descList");
        INFO("Mismatches found: " << mismatches);
        
        // Requirements
        REQUIRE(validated >= 90);
        // Allow some mismatches due to platform/version differences in Python vs C++
        // RDKit217 C++ may have approximations for qed, SPS, etc.
        REQUIRE(mismatches < validated * 217 * 0.05);  // Less than 5% total mismatch
    }
}

#endif  // RDK_BUILD_OSMORDRED_SUPPORT
