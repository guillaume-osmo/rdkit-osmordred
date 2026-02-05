#include "logodv_ppmv_v2.h"
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
static const int FEATURE_MASK[1000] = {300, 1189, 2152, 337, 93, 1812, 413, 422, 320, 3810, 324, 1282, 939, 1906, 1813, 364, 3136, 997, 333, 2153, 3852, 3765, 334, 315, 346, 904, 1948, 229, 1892, 1916, 1297, 1147, 374, 3829, 1366, 733, 338, 4041, 822, 1292, 886, 1882, 2409, 1820, 1154, 1446, 3786, 389, 325, 876, 1029, 344, 995, 474, 1315, 1828, 3766, 2669, 398, 1349, 912, 1328, 227, 813, 348, 299, 1737, 643, 1912, 1635, 335, 323, 1695, 584, 1834, 1204, 1085, 2237, 2235, 867, 1232, 3102, 820, 3399, 1037, 1984, 409, 352, 3767, 1841, 1866, 2239, 805, 1888, 1005, 1206, 1150, 304, 1391, 555, 778, 403, 1161, 569, 1319, 128, 1944, 1863, 312, 1125, 892, 1978, 2269, 1913, 1152, 841, 916, 1946, 1904, 390, 890, 2148, 1323, 832, 1389, 230, 1279, 550, 1348, 1337, 598, 931, 419, 523, 907, 1036, 698, 968, 399, 353, 1607, 614, 1857, 879, 1629, 1335, 3834, 385, 1506, 1850, 965, 313, 966, 2162, 1025, 1252, 790, 2255, 1843, 878, 3819, 1006, 1103, 976, 450, 796, 3805, 1878, 768, 933, 3844, 298, 1067, 362, 1304, 1276, 949, 580, 3803, 831, 1825, 780, 1038, 752, 1168, 121, 150, 596, 3814, 1884, 411, 3850, 3500, 396, 1126, 2009, 2252, 1957, 1117, 640, 932, 871, 421, 654, 2212, 1535, 2155, 1372, 751, 893, 1013, 375, 825, 1530, 1051, 506, 1142, 1156, 1172, 1052, 3774, 732, 639, 3858, 1035, 587, 3833, 807, 845, 1908, 3768, 452, 996, 674, 941, 3401, 3841, 557, 787, 1891, 3839, 561, 990, 855, 1826, 1872, 1595, 1864, 811, 2209, 1510, 372, 3877, 779, 2248, 1961, 359, 590, 1171, 1141, 776, 341, 351, 1137, 521, 177, 1853, 1387, 1886, 1911, 3157, 991, 1065, 1046, 837, 723, 2272, 3811, 921, 1317, 764, 1924, 2344, 1096, 2149, 1330, 2254, 644, 326, 777, 369, 3849, 4087, 1030, 3843, 1533, 1118, 679, 953, 509, 1846, 3201, 1609, 1854, 1062, 402, 1313, 425, 756, 1020, 366, 3769, 1982, 1611, 1648, 3868, 1119, 717, 370, 903, 1040, 645, 1224, 331, 1238, 619, 1875, 750, 1095, 1068, 2271, 602, 3802, 3443, 295, 2215, 1805, 1511, 869, 706, 2268, 566, 849, 863, 3798, 570, 1855, 1867, 3135, 842, 330, 838, 546, 1696, 852, 1113, 1977, 1053, 1934, 1350, 1749, 3076, 791, 652, 174, 1958, 952, 1093, 3801, 1953, 1170, 775, 513, 2170, 664, 913, 3501, 3104, 1370, 648, 586, 340, 1288, 1947, 1959, 781, 663, 3566, 601, 1829, 851, 544, 1898, 812, 1112, 889, 1869, 2166, 1331, 735, 971, 910, 2412, 797, 3812, 3038, 1314, 1015, 1895, 319, 770, 2186, 1814, 1824, 1612, 1071, 318, 693, 985, 388, 310, 3437, 710, 1129, 937, 572, 908, 731, 3072, 1847, 176, 1840, 453, 1874, 760, 3041, 3073, 870, 989, 1107, 1160, 2414, 798, 3103, 843, 865, 874, 1109, 1008, 3861, 884, 734, 1023, 2201, 3078, 1937, 834, 1087, 1105, 983, 1018, 2346, 724, 630, 1176, 803, 1838, 3781, 2707, 327, 1735, 305, 1515, 1101, 872, 827, 1076, 1845, 671, 40, 1007, 564, 2430, 628, 850, 1104, 1091, 823, 814, 1918, 1859, 999, 624, 329, 686, 2005, 2471, 2050, 829, 547, 747, 511, 568, 1534, 982, 1822, 656, 1708, 687, 754, 689, 960, 2708, 1016, 979, 1110, 738, 1080, 1885, 1088, 3848, 2156, 1271, 1042, 861, 681, 984, 1508, 549, 1594, 593, 1009, 895, 1021, 950, 2167, 1926, 1086, 1950, 321, 806, 1809, 699, 667, 801, 1807, 1960, 836, 1044, 864, 1057, 1034, 959, 3144, 651, 757, 808, 1818, 1162, 675, 563, 934, 859, 3074, 894, 2241, 1139, 784, 3832, 1049, 1039, 2373, 860, 2026, 1108, 1316, 714, 848, 763, 1613, 800, 294, 2006, 381, 676, 1127, 1077, 1081, 2372, 576, 2022, 881, 1998, 1064, 1842, 1024, 3441, 1145, 1258, 964, 3773, 936, 1070, 696, 1178, 684, 339, 3776, 2099, 653, 1186, 1881, 817, 1434, 4048, 391, 2003, 729, 1289, 1928, 2001, 1083, 1131, 961, 317, 1165, 919, 332, 925, 692, 926, 303, 1322, 819, 1223, 718, 3845, 1256, 349, 940, 765, 923, 394, 378, 719, 844, 1694, 3036, 306, 1134, 818, 1050, 1054, 1169, 930, 309, 316, 1308, 944, 1055, 629, 975, 1714, 744, 15, 942, 3468, 350, 1909, 1157, 915, 665, 1512, 846, 970, 1001, 2150, 1353, 519, 682, 948, 459, 1063, 1136, 761, 1299, 293, 2274, 945, 2406, 1045, 314, 1079, 3764, 974, 862, 328, 1061, 1159, 1114, 981, 1073, 3800, 1720, 774, 762, 1717, 1868, 512, 1041, 1058, 358, 3467, 1887, 2143, 1352, 1910, 810, 1817, 824, 666, 1163, 2247, 828, 3788, 794, 918, 1980, 847, 1981, 1351, 662, 958, 816, 1883, 3777, 978, 2165, 3831, 773, 771, 1373, 1106, 3439, 711, 1979, 766, 1877, 565, 1121, 357, 1048, 769, 726, 659, 856, 1263, 857, 963, 1128, 737, 2251, 1718, 2405, 2174, 799, 767, 170, 1942, 1130, 1069, 947, 1354, 868, 1806, 885, 914, 1111, 1943, 854, 725, 749, 1929, 1278, 1358, 680, 909, 703, 3806, 2275, 1059, 1135, 1245, 1309, 708, 1336, 1140, 1032, 920, 2267, 311, 1082, 1056, 741, 877, 1033, 575, 146, 3770, 954, 1074, 1856, 728, 2779, 182, 746, 2210, 1356, 1939, 1155, 2298, 2307, 291, 2347, 3785, 1173, 508, 2261, 891, 3840, 1097, 967, 2359, 1099, 1920, 2172, 1345, 1254, 1585, 2000, 1871, 3089, 1894, 2216, 821, 1149, 830, 3594, 720, 1619, 1377, 1340, 1393, 655, 748, 883, 1031, 2177, 745, 1166, 2348, 1889, 1012, 1963, 449, 3835, 3996, 392, 1078, 3418, 833, 782, 297, 380, 3847, 2169, 875, 992, 715, 135, 1138, 880, 882, 1844, 986, 1320, 604, 972, 3040, 866, 446, 858, 1072, 2343, 1810, 2164, 2178, 673, 1014, 743, 1060, 951, 1983, 1339, 296, 516, 141, 957, 727, 955, 3403, 1269, 1355, 1022, 1712, 1090, 789, 3820, 292, 804, 788, 928, 677, 302, 1830, 2374, 906, 694, 356, 1816, 672, 496, 581, 987, 3407, 1839, 2342, 2168, 75, 574, 815, 740, 3110, 1893, 393, 1296, 1346, 2208, 1616, 379, 977, 1043, 1124, 988, 4018, 2191, 802, 603, 1897, 650, 395, 873, 962, 956, 1938, 2232, 3780, 3405, 368, 1952, 307, 1710, 336, 701, 1200, 1183, 2306, 649, 1305, 935, 697, 993, 1333, 685, 185, 1089, 367, 1608, 924, 510, 3142, 1870, 3108, 3836, 716, 905, 27, 695, 383, 758};

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

logODV_ppmvV2Predictor::logODV_ppmvV2Predictor() {
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

logODV_ppmvV2Predictor::~logODV_ppmvV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool logODV_ppmvV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double logODV_ppmvV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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

void logODV_ppmvV2Predictor::predictBatch(
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

double logODV_ppmvV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void logODV_ppmvV2Predictor::predictFromSMILESBatch(
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

void logODV_ppmvV2Predictor::extractSMARTSFeatures(
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

void logODV_ppmvV2Predictor::extractRDKitFeatures(
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

void logODV_ppmvV2Predictor::extractOsmordredFeatures(
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

void logODV_ppmvV2Predictor::extractUnifiedFeatures(
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

void logODV_ppmvV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void logODV_ppmvV2Predictor::selectFeatures(
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
