#include "deltahvap_v2.h"
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
static const int FEATURE_MASK[800] = {1296, 3782, 2707, 1979, 413, 2231, 637, 3820, 2151, 357, 598, 1204, 422, 1892, 2006, 3770, 369, 326, 3809, 549, 547, 3818, 520, 2251, 1995, 1151, 3844, 658, 67, 215, 1264, 3804, 1262, 3789, 1960, 3870, 1811, 2149, 2152, 611, 374, 410, 506, 608, 814, 1326, 550, 588, 1375, 1596, 1735, 562, 3819, 703, 3765, 3952, 3869, 1376, 1809, 2004, 913, 1191, 2254, 1373, 2304, 3864, 2274, 1804, 389, 1206, 328, 1805, 1374, 3764, 3799, 1100, 332, 576, 2669, 1822, 1591, 398, 931, 587, 1259, 685, 737, 563, 558, 551, 3780, 2303, 15, 325, 1816, 2208, 2237, 1948, 776, 660, 586, 1323, 1142, 1332, 1333, 813, 664, 3768, 3816, 421, 2153, 885, 3814, 1594, 1398, 1815, 1386, 713, 3834, 804, 624, 730, 632, 298, 686, 700, 95, 702, 649, 589, 218, 1110, 3839, 2197, 1956, 1861, 341, 948, 1067, 1114, 775, 711, 1533, 579, 3781, 1997, 919, 548, 667, 666, 1835, 673, 721, 508, 333, 1966, 882, 584, 771, 1161, 974, 1910, 805, 640, 628, 518, 1349, 1842, 585, 336, 366, 1813, 1377, 522, 2022, 3897, 1860, 704, 3766, 3803, 338, 3849, 631, 1710, 393, 695, 939, 822, 1855, 812, 438, 352, 1309, 731, 764, 1007, 1348, 317, 2297, 613, 709, 594, 1079, 1269, 2212, 650, 3609, 2165, 1192, 1055, 3034, 2302, 1171, 741, 1310, 924, 1111, 1005, 768, 759, 1430, 3767, 509, 806, 1975, 388, 1155, 1074, 722, 569, 3089, 1163, 1178, 674, 140, 1415, 826, 1201, 1531, 1803, 1288, 306, 677, 337, 916, 3797, 1737, 1506, 902, 90, 979, 343, 1092, 1140, 93, 1354, 718, 331, 580, 330, 1854, 3140, 572, 2005, 1812, 760, 1095, 1871, 1695, 1085, 353, 1903, 905, 1907, 2025, 1037, 795, 1371, 1038, 949, 831, 1636, 1162, 696, 610, 89, 940, 3769, 1980, 2201, 350, 394, 2154, 61, 1888, 1891, 1983, 3845, 1077, 2209, 2164, 1107, 1973, 904, 1066, 1080, 1981, 1338, 531, 751, 2183, 1275, 944, 155, 1102, 744, 302, 1843, 1921, 1146, 833, 1280, 1985, 582, 128, 1071, 643, 708, 3824, 772, 3202, 858, 348, 1058, 625, 40, 412, 987, 313, 1021, 950, 1130, 626, 1853, 2272, 1059, 543, 880, 439, 546, 1014, 777, 2174, 1134, 1057, 1629, 1887, 629, 647, 2298, 925, 1174, 735, 841, 470, 1113, 516, 838, 1256, 725, 1016, 830, 912, 3812, 739, 811, 865, 1895, 570, 1150, 34, 742, 1821, 1101, 1088, 1008, 2410, 1141, 557, 960, 135, 680, 1136, 1157, 312, 1156, 1720, 1046, 1347, 684, 1978, 386, 150, 733, 825, 1018, 1875, 701, 1118, 1828, 3504, 1133, 914, 683, 307, 2203, 575, 296, 727, 1087, 953, 951, 1041, 1858, 216, 3785, 3848, 1814, 1353, 778, 1126, 3846, 973, 959, 1806, 875, 1029, 1824, 1346, 1075, 918, 1048, 732, 1515, 1363, 1336, 1109, 779, 419, 556, 1867, 0, 657, 405, 1168, 788, 1881, 1022, 796, 1982, 561, 1355, 134, 1889, 1098, 1073, 1839, 1901, 1518, 327, 815, 1616, 1864, 898, 1534, 2070, 891, 2250, 638, 1112, 1918, 1064, 1634, 661, 1808, 743, 3838, 901, 1202, 816, 1845, 1056, 1514, 687, 2252, 3833, 1598, 878, 1061, 1285, 1977, 957, 2160, 1083, 1998, 1361, 786, 997, 1070, 146, 2409, 358, 1160, 862, 1863, 756, 675, 1006, 792, 769, 3441, 989, 560, 1957, 1340, 682, 1072, 1972, 874, 1819, 124, 975, 2271, 1337, 573, 692, 1857, 345, 415, 1013, 129, 1173, 2110, 401, 850, 861, 3077, 1965, 941, 1524, 1289, 1305, 892, 750, 309, 1880, 1877, 1054, 1017, 294, 1053, 1869, 2026, 1050, 1940, 868, 981, 734, 817, 308, 367, 807, 676, 824, 355, 1030, 396, 1324, 63, 344, 33, 597, 1911, 1615, 652, 392, 2182, 642, 391, 70, 1078, 2007, 1717, 295, 836, 932, 390, 74, 3829, 1138, 1625, 1329, 966, 1127, 363, 1372, 1069, 1872, 1104, 1295, 368, 689, 349, 318, 1708, 1040, 797, 126, 315, 3854, 100, 761, 399, 1984, 668, 859, 27, 969, 2406, 2186, 698, 339, 329, 3139, 2241, 903, 219, 1320, 1550, 181, 1089, 3144, 753, 566, 876, 887, 1328, 1893, 952, 961, 774, 808, 1291, 819, 1976, 321, 1105, 937, 840, 1873, 888, 1311, 965, 1958, 1914, 1164, 883, 1898, 1147, 2860, 2179, 636, 860, 2073, 983, 1019, 679, 1321, 340, 2210, 2150, 909, 789, 1950, 1045, 879, 3850, 935, 1128, 852, 316, 1894, 800, 1944, 1120, 347, 1941, 2155, 1697, 429, 1718, 1119, 1278, 934, 837, 1116, 1699, 2512, 848, 1158, 1327, 717, 428, 917, 766, 1818, 728, 1961, 1012, 889, 521, 991, 1149, 829, 513, 1039, 723, 2495, 654, 688, 784, 659, 1154, 1830, 1963, 1947, 545, 3841, 1876, 153, 1486, 2445, 749, 583, 2709, 1638, 2167, 845, 3072, 403, 870, 1117, 645, 601, 896, 922, 910, 107, 3136, 1153, 2473, 656, 943, 691, 4010, 3406, 1826, 2414, 1511, 199, 832, 719, 681, 1051, 2000, 787, 853, 801, 1810, 1167, 1654, 1096, 736, 593, 2168, 877, 1135, 147, 1002, 873, 1317, 1031, 1936, 214, 1820, 1010, 954, 762, 1065, 754};

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

deltaHvapV2Predictor::deltaHvapV2Predictor() {
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

deltaHvapV2Predictor::~deltaHvapV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool deltaHvapV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double deltaHvapV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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

void deltaHvapV2Predictor::predictBatch(
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

double deltaHvapV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void deltaHvapV2Predictor::predictFromSMILESBatch(
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

void deltaHvapV2Predictor::extractSMARTSFeatures(
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

void deltaHvapV2Predictor::extractRDKitFeatures(
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

void deltaHvapV2Predictor::extractOsmordredFeatures(
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

void deltaHvapV2Predictor::extractUnifiedFeatures(
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

void deltaHvapV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void deltaHvapV2Predictor::selectFeatures(
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
