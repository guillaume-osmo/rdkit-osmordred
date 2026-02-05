#include "bp_v2.h"
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
static const int FEATURE_MASK[800] = {1332, 1997, 3770, 1355, 1356, 1811, 3819, 1993, 2151, 565, 1263, 3813, 1375, 3778, 2152, 571, 1537, 374, 3767, 1804, 2153, 2182, 2303, 95, 399, 1977, 2669, 1230, 580, 624, 660, 1373, 1852, 1830, 876, 299, 1849, 1376, 777, 3768, 422, 3846, 642, 3766, 533, 319, 317, 581, 1243, 3801, 1861, 572, 398, 3837, 1165, 2061, 1434, 4012, 643, 1720, 410, 1258, 366, 1864, 545, 573, 587, 3783, 1611, 2310, 1824, 332, 1591, 1954, 3765, 1005, 726, 1857, 626, 3869, 1948, 1353, 1552, 2175, 389, 1535, 598, 1340, 1843, 2274, 1211, 1998, 1860, 421, 1986, 1374, 627, 409, 1966, 1348, 1149, 1378, 1321, 1534, 1342, 1856, 1055, 1313, 964, 795, 4011, 1372, 704, 1821, 1905, 969, 582, 1853, 3835, 450, 2520, 1854, 2297, 1956, 2026, 2013, 4049, 343, 1902, 2150, 550, 89, 1224, 3832, 1357, 1346, 4041, 561, 1991, 659, 1113, 3834, 1599, 3897, 1840, 3955, 1344, 298, 904, 552, 2208, 2261, 727, 715, 1497, 1835, 3820, 2211, 2198, 526, 823, 698, 1847, 3829, 323, 315, 1151, 1698, 3764, 395, 3769, 1310, 544, 445, 2254, 1283, 1516, 375, 675, 1072, 327, 1654, 362, 355, 1256, 346, 769, 966, 1204, 1506, 3781, 3780, 805, 3815, 892, 1129, 1832, 1288, 3034, 1813, 412, 2000, 2364, 2199, 669, 676, 2202, 1850, 297, 439, 1600, 1968, 1979, 1049, 895, 1036, 588, 1012, 1939, 534, 932, 4092, 353, 729, 1884, 3785, 52, 1837, 2181, 1058, 590, 750, 1914, 344, 1377, 1598, 352, 1957, 1594, 822, 799, 889, 604, 872, 2191, 1855, 1040, 701, 1607, 1015, 796, 1584, 569, 1365, 652, 920, 3407, 1354, 1037, 51, 2234, 586, 1060, 2205, 707, 1334, 2323, 846, 1029, 1892, 341, 1901, 148, 650, 600, 1447, 663, 1809, 2298, 348, 392, 584, 1172, 150, 236, 1114, 388, 1056, 498, 1885, 420, 974, 589, 531, 625, 1074, 1974, 180, 950, 685, 807, 385, 4058, 2169, 3805, 1906, 230, 90, 83, 1981, 3782, 2070, 434, 1109, 667, 662, 45, 2212, 994, 1021, 865, 1271, 1252, 2161, 3850, 680, 832, 2156, 733, 903, 1135, 1070, 702, 1161, 1807, 1367, 3843, 226, 176, 913, 984, 3848, 1301, 1312, 459, 936, 957, 2406, 1735, 1593, 679, 1635, 877, 1304, 2527, 1333, 101, 644, 1155, 1345, 868, 1878, 879, 767, 614, 1121, 886, 858, 329, 1117, 772, 1839, 1280, 1136, 1297, 1995, 1152, 1953, 1319, 1595, 993, 764, 564, 3904, 302, 1717, 152, 2253, 894, 304, 1891, 1014, 1011, 947, 811, 338, 1001, 1169, 2033, 2237, 1148, 1347, 1621, 1020, 17, 838, 1978, 646, 1818, 1289, 1846, 1019, 786, 1694, 1937, 596, 2512, 995, 1125, 2251, 1925, 979, 182, 961, 980, 818, 2771, 1950, 1692, 1967, 1013, 356, 2154, 303, 812, 2837, 1918, 725, 1079, 983, 743, 851, 636, 1156, 367, 797, 1137, 337, 516, 684, 1320, 391, 361, 1649, 761, 130, 1016, 2196, 1102, 978, 649, 1054, 1295, 1118, 1736, 706, 345, 1863, 1999, 378, 940, 3838, 294, 1085, 831, 480, 555, 474, 1349, 1101, 929, 350, 1171, 1890, 1008, 996, 316, 386, 1371, 1124, 1150, 1082, 1084, 2304, 1921, 1894, 1154, 780, 658, 305, 2414, 2166, 689, 692, 1873, 734, 1178, 1590, 1982, 1324, 845, 1050, 647, 1755, 372, 2427, 1817, 976, 3441, 1170, 3914, 311, 468, 1157, 2343, 783, 2710, 2379, 1110, 697, 4009, 973, 1041, 4037, 914, 1064, 2232, 682, 93, 560, 732, 804, 1339, 766, 1103, 2671, 423, 1116, 546, 2231, 1164, 3443, 575, 2003, 1073, 1282, 2362, 3898, 1111, 911, 778, 2006, 1926, 1062, 924, 2157, 518, 1936, 381, 170, 2792, 1293, 308, 1580, 1007, 415, 559, 438, 898, 543, 981, 854, 2268, 2163, 700, 859, 752, 1514, 713, 768, 991, 521, 1802, 1046, 2210, 1834, 1751, 747, 1094, 1078, 405, 3095, 958, 1959, 1366, 466, 2495, 4086, 2780, 653, 295, 1911, 688, 2346, 1281, 1888, 387, 223, 146, 921, 608, 1146, 850, 1106, 951, 1548, 413, 102, 1609, 1309, 347, 219, 887, 867, 1126, 1810, 2470, 787, 623, 1207, 1530, 788, 514, 1162, 905, 1893, 1159, 1153, 2347, 3842, 106, 798, 2230, 1722, 1127, 869, 3109, 998, 1631, 340, 1525, 1133, 896, 749, 4069, 3322, 946, 982, 771, 1872, 915, 673, 934, 1961, 1180, 843, 1000, 144, 668, 3607, 2348, 820, 3775, 891, 1859, 1786, 1322, 737, 1276, 1053, 2062, 475, 953, 828, 1632, 578, 852, 140, 1045, 296, 2250, 1976, 2216, 1308, 314, 394, 735, 871, 1990, 2490, 3399, 694, 134, 1077, 944, 2770, 3041, 670, 1167, 751, 1928, 1940, 109, 1145, 2149, 4006, 4067, 1612, 3797, 1030, 885, 1808, 1087, 1409, 1826, 857, 1750, 742, 1903, 1820, 1739, 756, 3074, 369, 2170, 1166, 3242, 931, 824, 977, 638, 1616, 4027, 1174, 333, 4072, 1801, 1827, 1626, 2273, 1168, 3847, 813, 1895, 1112, 1369, 1531, 570, 880, 1912, 1867, 3590, 2271, 3861, 641, 416, 393, 2249, 840, 987, 763, 893, 4071, 2001, 1032, 1806, 87, 3135, 2009, 2260, 3438, 1713, 791, 890, 2500, 603, 648, 1838, 1969, 3076};

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

BPV2Predictor::BPV2Predictor() {
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

BPV2Predictor::~BPV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool BPV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double BPV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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
    std::vector<double> selected_features(800);
    for (int i = 0; i < 800; ++i) {
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

void BPV2Predictor::predictBatch(
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

double BPV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void BPV2Predictor::predictFromSMILESBatch(
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

void BPV2Predictor::extractSMARTSFeatures(
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

void BPV2Predictor::extractRDKitFeatures(
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

void BPV2Predictor::extractOsmordredFeatures(
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

void BPV2Predictor::extractUnifiedFeatures(
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

void BPV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void BPV2Predictor::selectFeatures(
    const std::vector<double>& all_features,
    std::vector<double>& selected
) const {
    // Use FEATURE_MASK from header
    for (int i = 0; i < 800; ++i) {
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
