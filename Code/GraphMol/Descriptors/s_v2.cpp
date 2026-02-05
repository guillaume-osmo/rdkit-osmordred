#include "s_v2.h"
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
static const int FEATURE_MASK[700] = {1386, 1809, 552, 2196, 3768, 2160, 1278, 2158, 550, 1859, 1280, 1851, 1811, 319, 563, 650, 1987, 626, 547, 2669, 894, 885, 399, 4041, 759, 1151, 2179, 116, 1494, 598, 1146, 568, 213, 1988, 374, 3769, 1117, 1178, 554, 567, 409, 1340, 1391, 2194, 2153, 2152, 543, 584, 564, 1360, 1918, 1594, 683, 1908, 3034, 1115, 1961, 3509, 1101, 3810, 1327, 1986, 403, 1211, 1132, 1964, 2161, 1109, 926, 1545, 3785, 3831, 1298, 1153, 1826, 1478, 652, 591, 1962, 1611, 557, 566, 1534, 1998, 442, 1348, 856, 2205, 2519, 1635, 2007, 648, 1444, 131, 1048, 828, 1204, 1937, 1414, 46, 3202, 1953, 2530, 1926, 1879, 150, 2297, 1981, 1810, 971, 565, 183, 332, 300, 1598, 1813, 534, 3845, 1796, 549, 719, 1338, 717, 1319, 392, 217, 1124, 140, 940, 1345, 647, 326, 3781, 714, 3770, 2268, 179, 1837, 366, 1267, 2155, 1165, 1992, 1354, 1827, 768, 1320, 416, 3110, 1701, 715, 1007, 2359, 1854, 1855, 1155, 725, 869, 1696, 1927, 1125, 960, 1995, 1308, 1060, 341, 4058, 548, 896, 815, 1357, 343, 1584, 3958, 3979, 653, 744, 386, 1339, 1839, 1337, 1825, 581, 702, 1681, 1331, 2274, 1967, 3144, 2159, 689, 3042, 1906, 381, 1289, 70, 1051, 2046, 851, 524, 913, 546, 420, 1090, 3852, 1997, 2512, 1167, 522, 2162, 1920, 1714, 2173, 1497, 1910, 2252, 1925, 405, 649, 1426, 796, 302, 2409, 2310, 1411, 346, 1812, 180, 3805, 783, 1595, 2126, 377, 1119, 338, 1037, 890, 421, 513, 1001, 1316, 1254, 423, 317, 1325, 2891, 340, 474, 824, 733, 1999, 511, 959, 1014, 355, 1157, 1172, 86, 858, 398, 1903, 76, 660, 450, 1173, 1013, 818, 1031, 1625, 3766, 2001, 344, 813, 887, 1282, 1835, 1583, 2199, 921, 745, 716, 925, 1015, 805, 1636, 372, 1245, 753, 358, 2542, 1057, 1330, 1114, 706, 1512, 2197, 822, 1849, 3867, 923, 1948, 4011, 3074, 1086, 3848, 707, 1287, 544, 316, 2208, 1887, 1508, 705, 671, 1698, 182, 1797, 147, 1525, 1118, 825, 1695, 1302, 2501, 1518, 1977, 1347, 3207, 1256, 390, 1735, 434, 681, 2237, 966, 2472, 1107, 1093, 94, 2730, 106, 3765, 922, 701, 2027, 1897, 611, 497, 3594, 1966, 1304, 909, 1921, 1822, 1617, 350, 438, 229, 561, 3782, 651, 315, 852, 1022, 1983, 3250, 538, 1973, 2890, 353, 1170, 1025, 45, 1158, 636, 2002, 1717, 394, 2427, 643, 3829, 107, 2251, 832, 3566, 821, 736, 1120, 375, 811, 3767, 857, 1159, 3956, 1800, 1528, 2837, 400, 1055, 314, 3078, 3095, 486, 311, 1384, 1368, 1982, 3868, 2186, 1171, 2201, 3136, 1, 51, 595, 1154, 3849, 1596, 573, 1232, 419, 2270, 1917, 1329, 1523, 1219, 798, 356, 2163, 1305, 934, 410, 1366, 1830, 2298, 727, 1061, 2240, 518, 755, 376, 367, 3800, 3140, 585, 958, 948, 2241, 2428, 1697, 359, 2169, 3568, 2957, 771, 2275, 1310, 1313, 891, 1317, 1867, 779, 583, 586, 181, 1868, 1029, 1522, 115, 348, 2254, 1537, 3401, 897, 3811, 700, 1038, 912, 1985, 1685, 3607, 378, 2170, 3774, 1591, 3076, 1950, 1850, 2414, 1161, 1945, 1008, 1430, 1892, 2132, 1731, 1589, 1126, 997, 1808, 867, 3833, 1535, 1078, 2257, 1510, 814, 756, 362, 677, 2347, 2165, 899, 1258, 743, 1806, 3469, 724, 878, 347, 883, 1843, 385, 708, 1058, 1628, 545, 1929, 1095, 1527, 788, 3399, 757, 1324, 2211, 309, 932, 853, 1152, 1963, 670, 1821, 1334, 1976, 91, 742, 1147, 903, 876, 806, 2000, 1271, 1063, 1053, 1059, 1968, 2453, 787, 1169, 985, 841, 1160, 1116, 1065, 1786, 1139, 95, 983, 761, 1708, 866, 1076, 1164, 2306, 1030, 877, 391, 2230, 686, 1068, 1106, 401, 1547, 1377, 957, 1882, 1162, 839, 1074, 1299, 1017, 1149, 1516, 4089, 1924, 1144, 1142, 1958, 699, 2499, 2210, 697, 4060, 1845, 1782, 1000, 393, 305, 1582, 1071, 684, 773, 1201, 2743, 2255, 3888, 1309, 1829, 1607, 1143, 1984, 389, 3780, 1112, 849, 323, 904, 3036, 413, 918, 627, 778, 655, 965, 1587, 352, 704, 777, 1371, 3041, 2262, 1846, 66, 696, 2271, 772, 882, 769, 1045, 1860, 666, 3102, 28, 1099, 1353, 1971, 1082, 342, 1080, 408, 49, 1956, 820, 2038, 1341, 3255, 2188, 999, 1168, 294, 3441, 2100, 1056, 560, 1121, 741, 770, 1333, 1110, 3764, 318, 303, 1111, 1021, 1789, 840, 1446, 1694, 1581, 3135, 3104, 388, 919, 4069, 435, 947, 1323, 1991, 676, 1730};

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

SV2Predictor::SV2Predictor() {
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

SV2Predictor::~SV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool SV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double SV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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
    std::vector<double> selected_features(700);
    for (int i = 0; i < 700; ++i) {
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

void SV2Predictor::predictBatch(
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

double SV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void SV2Predictor::predictFromSMILESBatch(
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

void SV2Predictor::extractSMARTSFeatures(
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

void SV2Predictor::extractRDKitFeatures(
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

void SV2Predictor::extractOsmordredFeatures(
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

void SV2Predictor::extractUnifiedFeatures(
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

void SV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void SV2Predictor::selectFeatures(
    const std::vector<double>& all_features,
    std::vector<double>& selected
) const {
    // Use FEATURE_MASK from header
    for (int i = 0; i < 700; ++i) {
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
