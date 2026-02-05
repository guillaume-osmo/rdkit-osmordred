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
static const int FEATURE_MASK[1000] = {333, 1375, 323, 337, 3770, 343, 411, 220, 1997, 374, 399, 353, 398, 649, 403, 228, 1374, 317, 1827, 95, 410, 1809, 392, 2005, 422, 3768, 363, 1285, 1956, 1377, 658, 82, 376, 2174, 511, 3780, 346, 319, 352, 331, 81, 406, 300, 856, 386, 324, 298, 1178, 299, 332, 1607, 389, 221, 1262, 1165, 1528, 1075, 1278, 621, 2502, 664, 624, 1005, 439, 1186, 3913, 1319, 2207, 3834, 1219, 1369, 3867, 1855, 814, 322, 301, 1129, 663, 905, 397, 1, 294, 625, 413, 2004, 1365, 2489, 2006, 3769, 960, 385, 325, 1446, 64, 557, 676, 1584, 1353, 2212, 1922, 315, 394, 2239, 1315, 3765, 2149, 467, 2345, 2022, 1976, 3783, 359, 1957, 544, 683, 179, 2170, 375, 2374, 1837, 320, 1252, 1260, 1628, 922, 3887, 229, 882, 1718, 823, 421, 3774, 998, 650, 330, 896, 724, 396, 750, 233, 1918, 643, 2490, 121, 1649, 1810, 1860, 877, 1882, 1987, 525, 506, 566, 1087, 657, 2312, 968, 329, 536, 2272, 1328, 534, 1737, 2025, 522, 759, 810, 1805, 1835, 3852, 1184, 1816, 521, 1324, 969, 1953, 3782, 438, 150, 365, 603, 302, 63, 929, 1512, 802, 83, 3136, 345, 760, 151, 933, 348, 1948, 181, 1170, 225, 101, 859, 1736, 1290, 1333, 2249, 1866, 388, 425, 1151, 126, 865, 571, 431, 1155, 355, 1914, 981, 4049, 1257, 362, 2046, 409, 728, 102, 936, 311, 1879, 468, 1157, 906, 1850, 1828, 987, 761, 1812, 1548, 1806, 850, 15, 835, 1944, 699, 1072, 1808, 138, 1049, 435, 366, 705, 90, 401, 805, 817, 1936, 1832, 1972, 2002, 1851, 1354, 1112, 738, 148, 1058, 972, 964, 2161, 1581, 2152, 3089, 775, 2342, 141, 1117, 826, 1427, 1095, 402, 1016, 953, 321, 1868, 2153, 2347, 1133, 769, 796, 1006, 297, 1950, 1885, 1646, 751, 546, 1949, 303, 547, 853, 3766, 76, 477, 328, 1047, 949, 295, 973, 628, 1853, 839, 1021, 709, 1854, 642, 1821, 370, 316, 1880, 583, 338, 2241, 1847, 371, 1958, 1017, 1306, 1684, 327, 954, 335, 1407, 1905, 1484, 372, 1979, 945, 1894, 1048, 803, 1348, 743, 1077, 997, 293, 889, 4041, 1998, 312, 2155, 822, 1232, 1152, 1037, 1318, 43, 1086, 821, 145, 3829, 1813, 1937, 1078, 551, 567, 1824, 956, 3399, 667, 1960, 1525, 1145, 39, 1135, 1973, 1527, 2268, 1982, 662, 632, 358, 2210, 1108, 935, 966, 1347, 369, 1281, 788, 1061, 918, 659, 845, 224, 350, 832, 419, 1534, 1934, 1820, 514, 1896, 1337, 869, 912, 786, 2405, 1906, 108, 340, 993, 1903, 133, 961, 1149, 1873, 314, 1054, 3946, 990, 860, 351, 886, 870, 1169, 1230, 86, 1350, 393, 975, 1688, 519, 1701, 576, 700, 1635, 1130, 2275, 1162, 1336, 1587, 978, 846, 387, 780, 309, 1109, 916, 1884, 1031, 2168, 361, 564, 444, 3906, 2164, 1163, 1099, 22, 582, 1271, 415, 1202, 3183, 3897, 808, 310, 296, 424, 1171, 768, 1978, 1940, 1060, 313, 799, 318, 1981, 4033, 1136, 885, 3053, 874, 305, 1079, 3074, 1892, 500, 1066, 1803, 510, 1523, 1939, 831, 341, 914, 652, 939, 2274, 412, 2216, 3838, 3848, 1091, 230, 1310, 1359, 635, 848, 25, 2412, 627, 1371, 2160, 1391, 588, 1046, 647, 977, 344, 1127, 423, 2159, 955, 395, 950, 1303, 622, 517, 665, 899, 1173, 1070, 1535, 1612, 2163, 524, 779, 122, 3877, 349, 1862, 1317, 193, 1925, 384, 1343, 849, 1022, 367, 994, 114, 1345, 304, 1057, 951, 1588, 1597, 291, 1124, 1013, 752, 903, 420, 2208, 836, 89, 306, 1355, 2169, 513, 770, 590, 3767, 1617, 1930, 861, 531, 2007, 2254, 1870, 1830, 2499, 1144, 2512, 723, 732, 416, 543, 339, 2344, 1116, 1097, 729, 578, 599, 619, 989, 1042, 1876, 34, 854, 887, 563, 767, 1971, 940, 1966, 1039, 910, 1138, 1153, 774, 986, 232, 2173, 2326, 555, 1871, 1841, 1920, 1082, 226, 1878, 1083, 1802, 3034, 655, 3804, 1911, 785, 1029, 907, 1105, 1111, 730, 1594, 2205, 894, 1912, 744, 391, 640, 1168, 13, 1838, 1339, 1886, 1692, 382, 3798, 1961, 800, 2455, 428, 1040, 1114, 1115, 1631, 211, 771, 140, 1358, 474, 1053, 3939, 932, 1732, 1717, 378, 1056, 773, 3847, 292, 3850, 383, 1913, 985, 1018, 697, 1785, 801, 180, 1030, 0, 1424, 552, 1253, 1340, 1506, 231, 772, 2165, 106, 1038, 648, 160, 2378, 1845, 828, 67, 827, 778, 891, 400, 1952, 948, 974, 965, 809, 550, 979, 680, 72, 923, 1128, 815, 1026, 944, 984, 604, 1010, 1055, 1341, 429, 3170, 4030, 2198, 3220, 137, 1909, 1405, 347, 1714, 1119, 354, 548, 2157, 3056, 1889, 307, 1258, 4074, 1834, 1293, 847, 227, 3868, 94, 958, 671, 1074, 2171, 867, 1910, 1019, 1141, 682, 390, 3073, 911, 1172, 2150, 1844, 782, 2592, 3504, 1690, 526, 766, 2541, 868, 3140, 158, 3108, 1893, 991, 561, 1351, 726, 651, 833, 812, 129, 645, 1624, 1332, 686, 1818, 1088, 1245, 132, 1902, 1977, 946, 742, 1545, 472, 1335, 1969, 68, 888, 755, 757, 2709, 813, 1858, 843, 1636, 523, 379, 1008, 449, 807, 1160, 1921, 368, 982, 1908, 2409, 781, 2255, 1122, 758, 1867, 418, 901, 3845, 754, 746, 580, 437, 1283, 2837, 646, 2414, 1370, 3818, 857, 864, 3805, 2167, 2001, 691, 957, 1600, 3246, 450, 336, 1857, 308, 2779, 1064, 1041, 863, 1020, 3888, 1126, 381, 113, 172, 1611, 2158, 1063, 749, 1044, 1967, 777, 1951, 1874, 1710, 1833, 1883, 1150, 783, 1012, 1654, 920, 80, 574, 3885, 1081, 1515, 731, 593, 1924, 154, 4044, 1067, 1916, 2303, 862, 2365, 1294, 1856, 3844, 2453, 3473, 1897, 1024, 1206, 73, 1143, 452, 765, 747, 2406, 3858, 1103, 737, 727, 1349, 898, 1807, 842, 3835, 1511, 879, 1859, 881, 1140, 1071, 3250, 852, 962, 357, 959, 1926, 342, 876, 161, 171, 762, 1025, 2302, 841, 753, 1304, 745, 1414, 1104, 787, 2774, 838, 952, 840, 1316, 1299, 2271, 1634, 3951, 2020, 1256, 2248, 1869, 3842, 1085, 1352, 3501, 1955, 1094, 834, 2009, 3764, 725, 1428, 1164, 4007, 1102, 855, 236, 795, 1305, 1932, 895, 2669, 890, 740, 1846, 1826, 1311, 572, 913, 851, 159, 688, 1943, 99, 592, 1593, 1089, 1863, 465, 1154, 804, 784, 937, 694, 3095, 2257, 1839, 1720, 2234, 1904, 569};

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
