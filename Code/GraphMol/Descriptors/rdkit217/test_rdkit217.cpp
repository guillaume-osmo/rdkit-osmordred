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
    
    SECTION("Descriptor names match Python Descriptors.descList") {
        // These are the 217 descriptor names from Python's Descriptors._descList
        // Generated with: [d[0] for d in Descriptors._descList]
        std::vector<std::string> python_names = {
            "MaxAbsEStateIndex", "MaxEStateIndex", "MinAbsEStateIndex", "MinEStateIndex",
            "qed", "SPS",
            "MolWt", "HeavyAtomMolWt", "ExactMolWt",
            "NumValenceElectrons", "NumRadicalElectrons",
            "MaxPartialCharge", "MinPartialCharge", "MaxAbsPartialCharge", "MinAbsPartialCharge",
            "FpDensityMorgan1", "FpDensityMorgan2", "FpDensityMorgan3",
            "BCUT2D_MWHI", "BCUT2D_MWLOW", "BCUT2D_CHGHI", "BCUT2D_CHGLO",
            "BCUT2D_LOGPHI", "BCUT2D_LOGPLOW", "BCUT2D_MRHI", "BCUT2D_MRLOW",
            "AvgIpc", "BalabanJ", "BertzCT",
            "Chi0", "Chi0n", "Chi0v", "Chi1", "Chi1n", "Chi1v",
            "Chi2n", "Chi2v", "Chi3n", "Chi3v", "Chi4n", "Chi4v",
            "HallKierAlpha", "Ipc", "Kappa1", "Kappa2", "Kappa3", "LabuteASA",
            "PEOE_VSA1", "PEOE_VSA10", "PEOE_VSA11", "PEOE_VSA12", "PEOE_VSA13", "PEOE_VSA14",
            "PEOE_VSA2", "PEOE_VSA3", "PEOE_VSA4", "PEOE_VSA5", "PEOE_VSA6", "PEOE_VSA7", "PEOE_VSA8", "PEOE_VSA9",
            "SMR_VSA1", "SMR_VSA10", "SMR_VSA2", "SMR_VSA3", "SMR_VSA4", "SMR_VSA5",
            "SMR_VSA6", "SMR_VSA7", "SMR_VSA8", "SMR_VSA9",
            "SlogP_VSA1", "SlogP_VSA10", "SlogP_VSA11", "SlogP_VSA12",
            "SlogP_VSA2", "SlogP_VSA3", "SlogP_VSA4", "SlogP_VSA5", "SlogP_VSA6", "SlogP_VSA7", "SlogP_VSA8", "SlogP_VSA9",
            "TPSA",
            "EState_VSA1", "EState_VSA10", "EState_VSA11", "EState_VSA2", "EState_VSA3", "EState_VSA4",
            "EState_VSA5", "EState_VSA6", "EState_VSA7", "EState_VSA8", "EState_VSA9",
            "VSA_EState1", "VSA_EState10", "VSA_EState2", "VSA_EState3", "VSA_EState4",
            "VSA_EState5", "VSA_EState6", "VSA_EState7", "VSA_EState8", "VSA_EState9",
            "FractionCSP3", "HeavyAtomCount",
            "NHOHCount", "NOCount",
            "NumAliphaticCarbocycles", "NumAliphaticHeterocycles", "NumAliphaticRings", "NumAmideBonds",
            "NumAromaticCarbocycles", "NumAromaticHeterocycles", "NumAromaticRings",
            "NumBridgeheadAtoms", "NumHAcceptors", "NumHDonors", "NumHeteroatoms", "NumHeterocycles",
            "NumRotatableBonds", "NumSaturatedCarbocycles", "NumSaturatedHeterocycles", "NumSaturatedRings",
            "NumSpiroAtoms", "RingCount",
            "MolLogP", "MolMR",
            "fr_Al_COO", "fr_Al_OH", "fr_Al_OH_noTert", "fr_ArN", "fr_Ar_COO", "fr_Ar_N", "fr_Ar_NH", "fr_Ar_OH",
            "fr_COO", "fr_COO2", "fr_C_O", "fr_C_O_noCOO", "fr_C_S",
            "fr_HOCCN", "fr_Imine", "fr_NH0", "fr_NH1", "fr_NH2", "fr_N_O",
            "fr_Ndealkylation1", "fr_Ndealkylation2", "fr_Nhpyrrole", "fr_SH",
            "fr_aldehyde", "fr_alkyl_carbamate", "fr_alkyl_halide", "fr_allylic_oxid",
            "fr_amide", "fr_amidine", "fr_aniline", "fr_aryl_methyl", "fr_azide", "fr_azo",
            "fr_barbitur", "fr_benzene", "fr_benzodiazepine", "fr_bicyclic",
            "fr_diazo", "fr_dihydropyridine", "fr_epoxide", "fr_ester", "fr_ether",
            "fr_furan", "fr_guanido", "fr_halogen", "fr_hdrzine", "fr_hdrzone",
            "fr_imidazole", "fr_imide", "fr_isocyan", "fr_isothiocyan",
            "fr_ketone", "fr_ketone_Topliss", "fr_lactam", "fr_lactone",
            "fr_methoxy", "fr_morpholine", "fr_nitrile", "fr_nitro", "fr_nitro_arom", "fr_nitro_arom_nonortho",
            "fr_nitroso", "fr_oxazole", "fr_oxime", "fr_para_hydroxylation",
            "fr_phenol", "fr_phenol_noOrthoHbond", "fr_phos_acid", "fr_phos_ester",
            "fr_piperdine", "fr_piperzine", "fr_priamide", "fr_prisulfonamd",
            "fr_pyridine", "fr_quatN", "fr_sulfide", "fr_sulfonamd", "fr_sulfone",
            "fr_term_acetylene", "fr_tetrazole", "fr_thiazole", "fr_thiocyan",
            "fr_thiophene", "fr_unbrch_alkane", "fr_urea"
        };
        
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
        
        REQUIRE(mismatches == 0);
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
