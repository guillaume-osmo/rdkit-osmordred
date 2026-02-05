#include "deltahf_v2.h"
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
static const int FEATURE_MASK[500] = {1320, 1916, 352, 1530, 316, 1498, 2155, 1844, 931, 411, 1429, 438, 1840, 481, 661, 1982, 1367, 376, 1204, 566, 813, 1296, 1414, 2211, 1349, 1258, 840, 4041, 1896, 564, 1845, 1173, 912, 561, 804, 535, 1948, 609, 95, 814, 1612, 1534, 1029, 3783, 563, 369, 1157, 359, 1972, 822, 1850, 805, 1125, 1828, 1295, 1926, 89, 375, 683, 2255, 660, 1151, 1219, 1163, 1159, 1205, 330, 1256, 1617, 647, 1892, 948, 823, 389, 686, 478, 334, 742, 1607, 357, 643, 664, 2060, 1922, 4044, 706, 777, 1818, 1206, 658, 3785, 530, 114, 1841, 2301, 1826, 922, 569, 3774, 570, 2009, 670, 1879, 1293, 1053, 1152, 1708, 399, 1809, 1891, 1372, 1310, 732, 1289, 543, 560, 548, 1824, 2210, 2208, 2152, 1927, 850, 297, 665, 1355, 534, 338, 1880, 957, 421, 1827, 403, 1180, 294, 1925, 545, 115, 2216, 298, 761, 1271, 1368, 652, 1936, 1307, 1171, 1021, 2231, 347, 361, 681, 1832, 674, 2251, 645, 1245, 1807, 309, 1117, 1872, 1350, 885, 641, 1881, 140, 1947, 2153, 1628, 2303, 626, 392, 2197, 904, 1167, 1520, 179, 2669, 332, 506, 2189, 344, 768, 223, 362, 138, 989, 433, 1983, 691, 1830, 1317, 1514, 368, 366, 567, 391, 1977, 855, 1288, 700, 1908, 1981, 1069, 565, 979, 1930, 741, 696, 849, 3781, 1155, 1109, 3399, 746, 672, 2405, 1162, 1038, 1232, 917, 1077, 1870, 3786, 1920, 651, 669, 1037, 913, 642, 315, 333, 370, 3765, 903, 1855, 1816, 343, 367, 3768, 2235, 2196, 1893, 1820, 1134, 899, 305, 393, 3796, 1829, 3769, 1087, 1825, 1127, 788, 1255, 1839, 336, 1813, 1885, 1102, 383, 1939, 546, 1093, 894, 733, 1149, 83, 655, 1854, 881, 1527, 302, 1984, 81, 1311, 1154, 860, 1897, 3329, 709, 1615, 308, 3780, 78, 514, 1625, 715, 1118, 573, 2174, 1867, 1303, 38, 2154, 1810, 1133, 1884, 1005, 1950, 215, 568, 1890, 1866, 841, 1598, 2275, 914, 1291, 1616, 2012, 770, 1059, 304, 1252, 409, 322, 724, 1302, 939, 1066, 86, 692, 687, 950, 1853, 3858, 1515, 314, 468, 858, 877, 1856, 303, 218, 832, 402, 773, 1736, 1193, 852, 1150, 2014, 318, 707, 1107, 714, 307, 1080, 1989, 20, 2198, 638, 388, 875, 815, 1822, 1371, 737, 386, 863, 750, 1588, 1088, 3798, 1712, 1146, 1516, 886, 1106, 547, 464, 1143, 1883, 1031, 1970, 869, 967, 3829, 1047, 1101, 1072, 970, 933, 1022, 1611, 1998, 321, 3502, 1158, 351, 601, 2228, 1104, 453, 1322, 876, 779, 806, 1966, 1836, 1168, 851, 329, 1120, 2164, 657, 959, 1014, 575, 2161, 919, 1837, 759, 603, 923, 107, 1282, 807, 690, 398, 842, 1013, 589, 1980, 1304, 883, 650, 2249, 723, 1172, 1055, 1928, 666, 1929, 1333, 1015, 331, 974, 1974, 2298, 1533, 2274, 1045, 2271, 291, 769, 1838, 1139, 941, 1869, 778, 1111, 799, 984, 820, 3770, 731, 1141, 3767, 656, 870, 79, 2261, 774, 326, 825, 1510, 1010, 663, 313, 182, 572, 550, 311, 1889, 1054, 340, 838, 574, 4072, 1166, 3777, 772, 1135, 2770, 2169, 355, 811, 1057, 744, 1335, 2260, 293, 833, 571, 1961, 1629, 1085};

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

deltaHfV2Predictor::deltaHfV2Predictor() {
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

deltaHfV2Predictor::~deltaHfV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool deltaHfV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double deltaHfV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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

void deltaHfV2Predictor::predictBatch(
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

double deltaHfV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void deltaHfV2Predictor::predictFromSMILESBatch(
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

void deltaHfV2Predictor::extractSMARTSFeatures(
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

void deltaHfV2Predictor::extractRDKitFeatures(
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

void deltaHfV2Predictor::extractOsmordredFeatures(
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

void deltaHfV2Predictor::extractUnifiedFeatures(
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

void deltaHfV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void deltaHfV2Predictor::selectFeatures(
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
