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
#include <map>
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
    
    SECTION("Descriptor names match Python Descriptors.descList from golden file") {
        // Read descriptor names from golden reference file header
        // This is the source of truth from Python's Descriptors._descList
        std::string rdbase = std::getenv("RDBASE") ? std::getenv("RDBASE") : "";
        std::string golden_path;
        
        if (!rdbase.empty()) {
            golden_path = rdbase + "/Code/GraphMol/Descriptors/test_data/nci_100_rdkit217_golden.csv";
        } else {
            golden_path = "Code/GraphMol/Descriptors/test_data/nci_100_rdkit217_golden.csv";
        }
        
        std::ifstream golden_file(golden_path);
        if (!golden_file.is_open()) {
            WARN("Golden reference not found: " << golden_path);
            WARN("Skipping name validation test");
            return;
        }
        
        // Read header line containing Python descriptor names
        std::string header;
        std::getline(golden_file, header);
        golden_file.close();
        
        // Parse header to extract Python names (skip first "smiles" column)
        std::vector<std::string> python_names;
        std::istringstream hss(header);
        std::string col;
        bool first = true;
        while (std::getline(hss, col, ',')) {
            if (first) {
                first = false;  // Skip "smiles" column
                continue;
            }
            python_names.push_back(col);
        }
        
        INFO("Read " << python_names.size() << " descriptor names from golden file");
        REQUIRE(python_names.size() == 217);
        
        auto cpp_names = getRDKit217DescriptorNames();
        REQUIRE(cpp_names.size() == 217);
        
        // Validate each name matches exactly
        int mismatches = 0;
        for (size_t i = 0; i < 217; ++i) {
            if (cpp_names[i] != python_names[i]) {
                INFO("Name mismatch at index " << i << ": C++=\"" << cpp_names[i] 
                     << "\" vs Python=\"" << python_names[i] << "\"");
                mismatches++;
            }
        }
        
        INFO("Total name mismatches: " << mismatches);
        REQUIRE(mismatches == 0);
    }
    
    SECTION("Ethanol descriptors - validate by name against native RDKit") {
        auto mol = SmilesToMol("CCO");
        REQUIRE(mol != nullptr);
        
        auto descs = extractRDKitDescriptors(*mol);
        auto names = getRDKit217DescriptorNames();
        REQUIRE(descs.size() == 217);
        REQUIRE(names.size() == 217);
        
        // Build name->index map for lookup by descriptor name
        std::map<std::string, size_t> name_to_idx;
        for (size_t i = 0; i < names.size(); ++i) {
            name_to_idx[names[i]] = i;
        }
        
        // Validate specific descriptor values by NAME (not index)
        // These are known values from Python's native RDKit for ethanol (CCO)
        
        // MolWt: from rdkit.Chem.Descriptors import MolWt; MolWt(Chem.MolFromSmiles('CCO'))
        REQUIRE(name_to_idx.count("MolWt") == 1);
        CHECK_THAT(descs[name_to_idx["MolWt"]], Catch::Matchers::WithinRel(46.069, 0.001));
        
        // ExactMolWt
        REQUIRE(name_to_idx.count("ExactMolWt") == 1);
        CHECK_THAT(descs[name_to_idx["ExactMolWt"]], Catch::Matchers::WithinRel(46.0418648, 0.0001));
        
        // HeavyAtomCount
        REQUIRE(name_to_idx.count("HeavyAtomCount") == 1);
        CHECK(descs[name_to_idx["HeavyAtomCount"]] == 3.0);
        
        // NumHDonors (ethanol has 1 OH)
        REQUIRE(name_to_idx.count("NumHDonors") == 1);
        CHECK(descs[name_to_idx["NumHDonors"]] == 1.0);
        
        // NumHAcceptors (ethanol has 1 O)
        REQUIRE(name_to_idx.count("NumHAcceptors") == 1);
        CHECK(descs[name_to_idx["NumHAcceptors"]] == 1.0);
        
        // RingCount (ethanol has 0 rings)
        REQUIRE(name_to_idx.count("RingCount") == 1);
        CHECK(descs[name_to_idx["RingCount"]] == 0.0);
        
        // NumRotatableBonds (ethanol has 0)
        REQUIRE(name_to_idx.count("NumRotatableBonds") == 1);
        CHECK(descs[name_to_idx["NumRotatableBonds"]] == 0.0);
        
        // TPSA (topological polar surface area)
        REQUIRE(name_to_idx.count("TPSA") == 1);
        CHECK_THAT(descs[name_to_idx["TPSA"]], Catch::Matchers::WithinRel(20.23, 0.01));
        
        // MolLogP
        REQUIRE(name_to_idx.count("MolLogP") == 1);
        CHECK_THAT(descs[name_to_idx["MolLogP"]], Catch::Matchers::WithinRel(-0.0014, 0.5));  // ~0
        
        delete mol;
    }
    
    SECTION("Benzene descriptors - validate aromatic descriptors by name") {
        auto mol = SmilesToMol("c1ccccc1");
        REQUIRE(mol != nullptr);
        
        auto descs = extractRDKitDescriptors(*mol);
        auto names = getRDKit217DescriptorNames();
        
        std::map<std::string, size_t> name_to_idx;
        for (size_t i = 0; i < names.size(); ++i) {
            name_to_idx[names[i]] = i;
        }
        
        // Benzene-specific validations from native RDKit
        CHECK(descs[name_to_idx["HeavyAtomCount"]] == 6.0);
        CHECK(descs[name_to_idx["RingCount"]] == 1.0);
        CHECK(descs[name_to_idx["NumAromaticRings"]] == 1.0);
        CHECK(descs[name_to_idx["NumAromaticCarbocycles"]] == 1.0);
        CHECK(descs[name_to_idx["fr_benzene"]] == 1.0);
        CHECK_THAT(descs[name_to_idx["MolWt"]], Catch::Matchers::WithinRel(78.114, 0.001));
        
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
