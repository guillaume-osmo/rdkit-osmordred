#include "logvp_v2.h"
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
static const int FEATURE_MASK[1000] = {3770, 333, 346, 1356, 323, 1340, 2151, 2234, 422, 1363, 533, 2301, 1902, 1916, 337, 1492, 2152, 374, 704, 101, 677, 1965, 324, 413, 535, 2213, 3034, 1326, 3765, 389, 298, 3820, 95, 1237, 3874, 1960, 2274, 3768, 2164, 398, 1162, 481, 330, 526, 1803, 410, 3407, 2187, 1891, 954, 2254, 1451, 3769, 4086, 1694, 3779, 730, 338, 703, 1865, 299, 343, 712, 1417, 731, 104, 686, 1377, 1889, 1821, 1506, 317, 1374, 558, 320, 1795, 2237, 1366, 894, 702, 1535, 1995, 397, 1953, 1843, 1591, 795, 394, 840, 300, 1735, 1005, 648, 352, 235, 386, 2005, 783, 4045, 1959, 319, 965, 1316, 631, 1996, 1948, 64, 3830, 228, 3777, 770, 663, 1901, 649, 82, 1966, 1375, 1791, 109, 550, 420, 359, 1596, 3794, 512, 3175, 1860, 399, 181, 1095, 304, 629, 1879, 396, 1944, 676, 3778, 1808, 1809, 1594, 331, 1873, 3849, 941, 579, 823, 1851, 1311, 551, 1413, 754, 353, 3075, 566, 89, 1615, 848, 3816, 1358, 1312, 332, 1354, 145, 903, 1695, 1153, 3819, 2240, 405, 1345, 685, 144, 2345, 1986, 406, 2210, 1341, 1957, 572, 846, 2798, 805, 1701, 557, 150, 421, 1630, 294, 621, 301, 315, 423, 334, 1151, 1855, 895, 931, 1514, 640, 1427, 520, 336, 3781, 148, 2211, 624, 1109, 2009, 1921, 438, 350, 1600, 2592, 3085, 3072, 841, 3567, 1320, 1161, 1831, 1159, 1166, 543, 1093, 3801, 642, 312, 950, 355, 1587, 325, 553, 1888, 1717, 886, 604, 3540, 790, 1347, 643, 797, 1904, 1635, 375, 510, 725, 1278, 3908, 1295, 833, 2355, 1852, 1968, 1021, 1737, 1905, 1838, 1980, 914, 409, 837, 714, 431, 1607, 83, 392, 788, 2007, 316, 589, 3888, 311, 722, 1729, 729, 1936, 1078, 1864, 1533, 1352, 789, 1810, 2153, 348, 1853, 75, 1126, 361, 516, 4066, 347, 2489, 3872, 376, 1826, 1894, 3916, 1878, 3897, 1055, 1146, 2163, 388, 2260, 1405, 1168, 402, 453, 857, 524, 1982, 63, 335, 650, 3946, 990, 1056, 1117, 2275, 1172, 1947, 1981, 2589, 2236, 1077, 814, 1357, 698, 2150, 630, 644, 562, 549, 716, 90, 1212, 412, 1169, 1198, 957, 2168, 477, 769, 425, 1789, 129, 706, 925, 439, 2177, 1304, 787, 1582, 1822, 2006, 2311, 1882, 1862, 1813, 1819, 385, 366, 627, 435, 3797, 2165, 2241, 1545, 1256, 4036, 3952, 2271, 544, 1133, 229, 2447, 180, 3780, 1149, 2394, 1861, 2771, 1688, 1825, 896, 401, 1029, 1624, 1232, 1991, 2710, 308, 1963, 1845, 1333, 297, 2409, 1710, 1534, 1154, 4074, 3846, 1224, 2423, 1185, 1518, 233, 580, 1167, 1322, 302, 391, 1054, 4011, 1850, 1645, 72, 3805, 737, 102, 1046, 309, 1598, 1175, 1962, 1045, 1141, 303, 1818, 3837, 806, 1342, 647, 2251, 1391, 1922, 1972, 3585, 328, 1087, 509, 1817, 1170, 33, 3833, 220, 1856, 1430, 15, 34, 740, 3785, 447, 613, 909, 3764, 16, 760, 2707, 1628, 327, 662, 1978, 18, 226, 3806, 236, 22, 3996, 1351, 429, 3399, 1925, 444, 3945, 1110, 329, 967, 779, 3842, 2711, 1926, 357, 610, 1288, 1070, 977, 1725, 1346, 3869, 1432, 362, 1835, 500, 1037, 1924, 1714, 834, 1147, 1067, 415, 1854, 2232, 339, 2175, 534, 608, 1287, 4062, 798, 43, 1025, 960, 341, 948, 1150, 2306, 511, 671, 1938, 693, 1961, 3053, 363, 2195, 1747, 758, 390, 1686, 1065, 1935, 3993, 1593, 1134, 111, 661, 1973, 1348, 1245, 1946, 966, 2230, 1934, 1097, 345, 1800, 1939, 38, 796, 695, 2166, 3919, 1122, 876, 1158, 3944, 2276, 3157, 1685, 525, 2257, 2205, 2407, 1806, 428, 1716, 1597, 1299, 1094, 470, 1977, 1063, 1619, 326, 313, 487, 318, 231, 1017, 989, 777, 314, 1353, 3271, 305, 968, 1796, 120, 161, 2445, 1367, 583, 3108, 1802, 3776, 734, 953, 1527, 1910, 937, 786, 1708, 2624, 684, 295, 855, 1317, 2189, 1160, 1955, 1548, 1903, 1086, 2272, 1309, 743, 3851, 106, 3887, 2495, 2155, 844, 2201, 407, 2346, 2058, 751, 1846, 739, 792, 2252, 1483, 996, 3210, 932, 419, 1319, 842, 291, 850, 1837, 2812, 459, 1884, 1917, 2298, 922, 683, 1720, 1131, 964, 2004, 545, 832, 3403, 924, 2256, 1484, 3141, 733, 310, 1998, 860, 1081, 427, 1993, 1824, 826, 395, 750, 849, 365, 1512, 2864, 2671, 949, 73, 1178, 1833, 39, 2270, 1355, 1841, 843, 1839, 1101, 2001, 1156, 3800, 2003, 1173, 3832, 3590, 1271, 1508, 2411, 878, 688, 827, 1171, 1157, 1064, 1335, 1255, 2776, 564, 2000, 3814, 943, 1001, 93, 2208, 224, 523, 3998, 762, 1285, 3954, 1096, 1914, 2359, 1443, 696, 4023, 1365, 824, 2674, 380, 825, 575, 141, 1857, 816, 464, 3767, 2250, 2194, 2159, 2171, 3848, 342, 474, 775, 1739, 86, 13, 3783, 4050, 1041, 728, 4092, 358, 1616, 1090, 904, 393, 1513, 746, 124, 1072, 2231, 1928, 1014, 372, 1289, 430, 1294, 1128, 675, 3850, 1525, 1291, 1617, 1292, 416, 874, 1816, 1815, 899, 1107, 1378, 1421, 1139, 1909, 1293, 322, 887, 2160, 1807, 1344, 4012, 1111, 132, 1152, 940, 1880, 588, 35, 1634, 654, 705, 1638, 838, 3142, 877, 3594, 354, 1636, 159, 1859, 1047, 869, 67, 2249, 1983, 1786, 933, 3473, 655, 3506, 2669, 3812, 774, 776, 1369, 780, 2405, 891, 3242, 3202, 3829, 2248, 138, 987, 1031, 1481, 2261, 1990, 489, 1130, 707, 680, 414, 2652, 1324, 822, 1872, 1497, 1073, 2269, 1085, 2176, 306, 1303, 3914, 1967, 1030, 2354, 2779, 369, 778, 1979, 1113, 1632, 2520, 865, 598, 1053, 921, 349, 1515, 820, 835, 1371, 935, 1118, 687, 1951, 1954, 1949, 2215, 3201, 1022, 1015, 723, 993, 732, 1919, 154, 4071, 839, 1507, 292, 3929, 652, 1974, 600, 2061, 1142, 1164, 378, 958, 3847, 1091, 3117, 1867, 1026, 2199, 862, 1016, 1023, 1155, 1135, 601, 3782, 1584, 1911, 561, 1531, 1310, 973, 1089, 0, 986, 1623, 726, 382, 3225, 1305, 1950, 1415, 851, 1999, 4033, 1407, 1208, 1829, 3876, 1057, 1334, 570, 1992, 530, 103, 1849, 2777, 4031, 868, 981, 3036, 656, 1336, 829, 2203, 892, 956, 923, 1952, 1898, 1165, 1522, 638, 913, 3109, 3774, 920, 660, 158, 563, 351, 1820, 915, 1551, 861, 784, 807, 2212, 1013, 881, 1581, 1811, 999, 3766, 871, 2599, 945, 1362, 859, 852, 218, 2149, 565, 1869, 21, 927, 1736, 1804, 1629, 1105, 1349};

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

logVPV2Predictor::logVPV2Predictor() {
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

logVPV2Predictor::~logVPV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool logVPV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double logVPV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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
    std::vector<double> selected_features(1000);
    for (int i = 0; i < 1000; ++i) {
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

void logVPV2Predictor::predictBatch(
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

double logVPV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void logVPV2Predictor::predictFromSMILESBatch(
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

void logVPV2Predictor::extractSMARTSFeatures(
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

void logVPV2Predictor::extractRDKitFeatures(
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

void logVPV2Predictor::extractOsmordredFeatures(
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

void logVPV2Predictor::extractUnifiedFeatures(
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

void logVPV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void logVPV2Predictor::selectFeatures(
    const std::vector<double>& all_features,
    std::vector<double>& selected
) const {
    // Use FEATURE_MASK from header
    for (int i = 0; i < 1000; ++i) {
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
