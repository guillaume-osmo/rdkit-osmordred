#include "v_v2.h"
#include <RDKit/SmilesParse/SmilesParse.h>
#include <RDKit/SmilesParse/SmartsParse.h>
#include <RDKit/GraphMol/RingInfo.h>
#include <thread>
#include <mutex>
#include <algorithm>
#include <execution>
#include <cmath>
#include <limits>

namespace osmo {
namespace physchemprop {
namespace v2 {

// Feature mask (selected feature indices)
static const int FEATURE_MASK[500] = {322, 321, 325, 422, 300, 324, 193, 802, 1279, 2179, 828, 334, 1901, 318, 17, 533, 1896, 812, 526, 611, 581, 831, 66, 531, 335, 1366, 337, 381, 591, 1363, 573, 1731, 575, 2003, 232, 3868, 70, 912, 89, 317, 1372, 299, 1211, 95, 357, 631, 572, 974, 968, 1838, 701, 4041, 1986, 552, 396, 593, 3824, 1821, 589, 304, 147, 717, 86, 1924, 1842, 1359, 60, 629, 419, 867, 416, 67, 1805, 2001, 417, 1217, 763, 309, 682, 1825, 998, 3803, 294, 406, 1341, 116, 160, 1241, 3829, 133, 320, 654, 582, 467, 3813, 1831, 336, 413, 2000, 710, 628, 329, 643, 759, 1269, 1898, 636, 1320, 397, 1281, 2007, 1850, 79, 115, 375, 1991, 584, 1999, 722, 326, 183, 736, 333, 101, 1337, 1839, 689, 3845, 421, 1367, 328, 607, 3833, 697, 368, 651, 369, 715, 122, 1916, 810, 2148, 716, 1338, 1989, 1353, 1807, 698, 978, 858, 3800, 822, 1305, 848, 2166, 1811, 291, 863, 940, 600, 327, 673, 885, 311, 609, 2160, 214, 1977, 1377, 3769, 1115, 544, 1342, 1852, 1340, 1058, 653, 713, 1243, 649, 293, 1407, 1028, 687, 372, 1321, 1837, 695, 1165, 685, 1173, 1263, 298, 348, 233, 98, 158, 3852, 2082, 1259, 1335, 1062, 1812, 2002, 674, 295, 2184, 330, 45, 2211, 1885, 1827, 1351, 652, 545, 854, 182, 896, 2267, 1278, 220, 2214, 857, 1334, 1070, 718, 1959, 1198, 3818, 1295, 386, 159, 1150, 376, 559, 838, 1258, 1264, 323, 1510, 1046, 308, 411, 1038, 1332, 554, 3842, 595, 1045, 1636, 1995, 989, 1357, 2209, 1013, 777, 87, 1011, 1903, 1262, 2151, 384, 795, 626, 813, 1006, 297, 1271, 1378, 1061, 1149, 2359, 2271, 310, 618, 107, 1980, 944, 3844, 1611, 3839, 1020, 1304, 1170, 577, 1360, 586, 1434, 131, 459, 1155, 319, 797, 1834, 1963, 615, 871, 739, 768, 815, 1737, 578, 999, 1157, 840, 660, 1330, 1718, 2202, 927, 708, 1902, 412, 2247, 721, 859, 1943, 3786, 1829, 1855, 392, 219, 696, 1944, 3794, 588, 709, 655, 174, 784, 3867, 1822, 3203, 908, 1029, 1059, 1331, 556, 340, 948, 1956, 913, 2171, 1256, 313, 860, 2150, 2434, 580, 345, 394, 314, 1860, 2500, 385, 592, 177, 2162, 800, 1950, 1358, 0, 390, 3783, 365, 632, 1047, 675, 1200, 1942, 393, 2210, 1151, 774, 407, 779, 997, 2707, 344, 656, 1333, 943, 1887, 929, 549, 1993, 1843, 402, 1890, 769, 93, 370, 1990, 1817, 1373, 882, 1072, 624, 1329, 1849, 3765, 361, 1156, 627, 741, 3595, 1319, 919, 957, 1207, 1253, 2022, 3594, 3847, 1911, 869, 3785, 951, 3807, 599, 367, 782, 788, 642, 985, 724, 786, 377, 1163, 1895, 420, 302, 1067, 391, 1152, 1809, 1127, 331, 2240, 1865, 1055, 1806, 742, 1955, 558, 757, 808, 1143, 1054, 904, 1854, 1922, 383, 137, 982, 1708, 711, 3832, 1874, 2427, 1056, 312, 902, 765, 350, 546, 1313, 3850, 918, 316, 1856, 1154, 1826, 1553, 1350, 1863, 1097, 950, 1182, 398, 851, 1717, 1960, 1267, 1859, 564, 738, 735, 1804, 884, 1832, 100, 359, 2034, 1239, 834, 849, 292, 3801, 1012};

// SMARTS patterns (loaded from AbrahamSMARTS.cpp)
// TODO: Load from shared SMARTS header

// Golden ratios (50 pairs)
static const std::pair<int, int> GOLDEN_RATIOS[50] = {
    {89, 113}, {29, 127}, {91, 113}, {108, 127}, {29, 113},
    {39, 113}, {54, 113}, {0, 127}, {54, 127}, {34, 113},
    {18, 113}, {90, 127}, {90, 113}, {18, 127}, {89, 127},
    {55, 113}, {108, 113}, {2, 127}, {39, 127}, {91, 127},
    {51, 113}, {9, 127}, {107, 127}, {51, 127}, {13, 113},
    {93, 127}, {9, 113}, {55, 127}, {93, 113}, {11, 127},
    {52, 127}, {92, 127}, {92, 113}, {53, 113}, {2, 113},
    {53, 127}, {52, 113}, {56, 113}, {34, 127}, {56, 127},
    {59, 127}, {85, 127}, {85, 113}, {38, 113}, {110, 113},
    {110, 127}, {60, 113}, {59, 113}, {3, 113}, {84, 127},
};

VV2Predictor::VV2Predictor() {
    // Load SMARTS patterns (from shared header)
    // TODO: Load from AbrahamSMARTS.h
    
    // Initialize golden ratios
    golden_ratios_ = {
        {89, 113}, {29, 127}, {91, 113}, {108, 127}, {29, 113},
        {39, 113}, {54, 113}, {0, 127}, {54, 127}, {34, 113},
        {18, 113}, {90, 127}, {90, 113}, {18, 127}, {89, 127},
        {55, 113}, {108, 113}, {2, 127}, {39, 127}, {91, 127},
        {51, 113}, {9, 127}, {107, 127}, {51, 127}, {13, 113},
        {93, 127}, {9, 113}, {55, 127}, {93, 113}, {11, 127},
        {52, 127}, {92, 127}, {92, 113}, {53, 113}, {2, 113},
        {53, 127}, {52, 113}, {56, 113}, {34, 127}, {56, 127},
        {59, 127}, {85, 127}, {85, 113}, {38, 113}, {110, 113},
        {110, 127}, {60, 113}, {59, 113}, {3, 113}, {84, 127},
    };
    
    // Initialize RDKit descriptor calculator
    desc_calc_ = new RDKit::Descriptors::MoleculeDescriptorCalculator(
        RDKit::Descriptors::_descList
    );
}

VV2Predictor::~VV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool VV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
    if (!mol) {
        return true;
    }
    
    // Count rings (SSSR - Smallest Set of Smallest Rings)
    RDKit::RingInfo* ringInfo = mol->getRingInfo();
    int numRings = ringInfo->numRings();
    
    // Count heavy atoms (non-hydrogen)
    int numHeavyAtoms = mol->getNumHeavyAtoms();
    
    // Filter: max 10 rings AND max 200 heavy atoms
    // Molecules exceeding these limits can cause Osmordred to hang
    return (numRings > 10 || numHeavyAtoms > 200);
}

double VV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
    // Handle null molecule
    if (!mol) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    // Filter large molecules that cause processing to hang
    if (isMoleculeTooLarge(mol)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    // Extract unified features
    std::vector<double> all_features(4093, 0.0);
    extractUnifiedFeatures(mol, all_features);
    
    // Preprocess
    preprocessFeatures(all_features);
    
    // Select features using FEATURE_MASK
    std::vector<double> selected_features(500);
    for (int i = 0; i < 500; ++i) {
        int idx = FEATURE_MASK[i];
        if (idx >= 0 && idx < static_cast<int>(all_features.size())) {
            selected_features[i] = all_features[idx];
        } else {
            selected_features[i] = 0.0;
        }
    }
    
    // Predict with XGBoost (parallel if many trees)
    return predictParallel(selected_features, numThreads);
}

void VV2Predictor::predictBatch(
    const std::vector<const RDKit::ROMol*>& mols,
    std::vector<double>& results,
    int numThreads
) const {
    results.resize(mols.size());
    
    // Use RDKit threading
    unsigned int nThreads = getNumThreadsToUse(numThreads);
    
    // Parallel prediction using std::execution::par
    std::transform(
        std::execution::par_unseq,
        mols.begin(),
        mols.end(),
        results.begin(),
        [this, numThreads](const RDKit::ROMol* mol) {
            return this->predict(mol, numThreads);
        }
    );
}

double VV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
    // Parse SMILES (handles failures gracefully)
    RDKit::ROMol* mol = RDKit::SmilesToMol(smiles);
    if (!mol) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    // Check if molecule is too large before processing
    if (isMoleculeTooLarge(mol)) {
        delete mol;
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    double result = predict(mol, numThreads);
    delete mol;
    return result;
}

void VV2Predictor::predictFromSMILESBatch(
    const std::vector<std::string>& smiles_list,
    std::vector<double>& results,
    int numThreads
) const {
    results.resize(smiles_list.size());
    
    // Use RDKit threading
    unsigned int nThreads = getNumThreadsToUse(numThreads);
    
    // Parallel SMILES parsing and prediction
    // Each thread handles: parse SMILES -> extract features -> predict
    std::transform(
        std::execution::par_unseq,
        smiles_list.begin(),
        smiles_list.end(),
        results.begin(),
        [this, numThreads](const std::string& smiles) {
            return this->predictFromSMILES(smiles, numThreads);
        }
    );
}

void VV2Predictor::extractSMARTSFeatures(
    const RDKit::ROMol* mol,
    std::vector<double>& features
) const {
    if (!mol) return;
    
    // Extract 241 base SMARTS + 50 golden ratios = 291 features
    // Features start at index 0
    for (size_t i = 0; i < smarts_queries_.size(); ++i) {
        if (smarts_queries_[i]) {
            features[i] = static_cast<double>(mol->getSubstructMatches(*smarts_queries_[i]).size());
        }
    }
    
    // Golden ratios (50)
    for (size_t i = 0; i < 50; ++i) {
        int idx1 = GOLDEN_RATIOS[i].first;
        int idx2 = GOLDEN_RATIOS[i].second;
        if (idx1 < 241 && idx2 < 241 && features[idx1] > 0 && features[idx2] > 0) {
            features[241 + i] = features[idx1] / features[idx2];
        }
    }
}

void VV2Predictor::extractRDKitFeatures(
    const RDKit::ROMol* mol,
    std::vector<double>& features
) const {
    if (!mol) return;
    
    // Extract 217 RDKit descriptors
    // Features start at index 291
    std::vector<double> descs = desc_calc_->calcDescriptors(*mol);
    for (size_t i = 0; i < descs.size() && i < 217; ++i) {
        features[291 + i] = descs[i];
    }
}

void VV2Predictor::extractOsmordredFeatures(
    const RDKit::ROMol* mol,
    std::vector<double>& features
) const {
    if (!mol) return;
    
    // Extract 3585 Osmordred features
    // Features start at index 508
    try {
        RDKit::Osmordred::Calculate(*mol, features.data() + 508, 3585);
    } catch (...) {
        // Handle errors gracefully - features remain 0.0
    }
}

void VV2Predictor::extractUnifiedFeatures(
    const RDKit::ROMol* mol,
    std::vector<double>& features
) const {
    if (!mol) {
        features.assign(4093, 0.0);
        return;
    }
    
    features.assign(4093, 0.0);
    extractSMARTSFeatures(mol, features);
    extractRDKitFeatures(mol, features);
    extractOsmordredFeatures(mol, features);
}

void VV2Predictor::preprocessFeatures(std::vector<double>& features) const {
    // Convert inf to nan
    for (double& f : features) {
        if (!std::isfinite(f)) {
            f = std::numeric_limits<double>::quiet_NaN();
        }
    }
    
    // Apply arcsinh for |x| > 33
    for (double& f : features) {
        if (std::isfinite(f) && std::abs(f) > 33.0) {
            f = std::asinh(f);
        }
    }
}

void VV2Predictor::selectFeatures(
    const std::vector<double>& all_features,
    std::vector<double>& selected
) const {
    // Use FEATURE_MASK from header
    for (int i = 0; i < 500; ++i) {
        int idx = FEATURE_MASK[i];
        if (idx >= 0 && idx < static_cast<int>(all_features.size())) {
            selected[i] = all_features[idx];
        } else {
            selected[i] = 0.0;
        }
    }
}

} // namespace v2
} // namespace physchemprop
} // namespace osmo
