#include "a_v2.h"
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
static const int FEATURE_MASK[700] = {2304, 2274, 338, 3765, 1936, 1903, 2237, 768, 1852, 398, 814, 410, 2152, 2779, 1843, 300, 823, 1288, 1916, 423, 1918, 389, 393, 181, 1167, 64, 140, 152, 1204, 1949, 818, 419, 3777, 809, 831, 1367, 432, 46, 742, 433, 217, 1478, 741, 1159, 1717, 1219, 1948, 805, 959, 1142, 1157, 1892, 1356, 1829, 2837, 45, 1983, 3095, 352, 3144, 3975, 452, 769, 876, 304, 687, 2302, 859, 940, 1155, 3960, 1071, 645, 1966, 2412, 1335, 150, 3767, 878, 374, 1690, 343, 1311, 376, 1078, 3074, 2212, 903, 973, 3460, 1483, 3842, 1515, 3829, 2241, 1823, 412, 2009, 650, 231, 1516, 908, 1015, 1879, 565, 690, 932, 721, 1714, 761, 1055, 3769, 983, 905, 776, 827, 180, 3817, 941, 303, 1891, 1590, 1206, 3919, 3950, 1281, 228, 2025, 226, 3993, 662, 1321, 2260, 660, 136, 826, 588, 397, 1256, 848, 4041, 948, 523, 2261, 106, 89, 851, 2132, 1251, 3781, 349, 1077, 833, 3868, 1580, 1925, 1854, 429, 1696, 1166, 2001, 2254, 151, 2149, 1443, 1351, 1847, 1086, 1340, 2303, 3974, 2050, 931, 1786, 3843, 1279, 358, 1163, 991, 924, 877, 810, 534, 395, 1898, 1154, 743, 1837, 392, 884, 770, 688, 1624, 345, 999, 882, 340, 18, 1125, 2256, 1046, 729, 986, 1106, 2165, 673, 308, 1980, 957, 1169, 1310, 1161, 2406, 879, 2407, 774, 1258, 1144, 355, 3835, 1866, 4011, 1853, 1143, 1061, 1944, 37, 738, 718, 3505, 758, 1818, 3918, 388, 1646, 1685, 674, 1150, 2174, 2471, 872, 1825, 1747, 1725, 1151, 684, 3766, 4092, 1800, 1596, 2561, 815, 695, 1016, 2512, 1534, 2208, 1145, 1970, 3034, 1352, 1718, 715, 1114, 1313, 601, 880, 2276, 979, 2006, 390, 1881, 3848, 2480, 1798, 1981, 1029, 1506, 692, 1822, 1064, 1021, 535, 591, 2472, 711, 1110, 646, 619, 1094, 2173, 4012, 1905, 3142, 1683, 975, 9, 1334, 3812, 982, 1808, 1135, 1850, 3939, 664, 942, 767, 1512, 1922, 1708, 804, 2355, 939, 1857, 1691, 700, 1091, 2255, 655, 907, 722, 1180, 486, 1014, 1893, 949, 1158, 510, 1322, 3785, 1867, 3972, 1117, 764, 1522, 3764, 1068, 341, 1978, 1834, 1621, 85, 739, 2416, 1584, 1882, 1595, 1085, 869, 1844, 497, 2501, 1828, 3399, 647, 544, 1062, 1038, 1350, 1294, 29, 1951, 1008, 958, 3596, 1096, 1926, 779, 667, 838, 2200, 182, 563, 965, 1484, 2216, 124, 158, 513, 1423, 2499, 2231, 701, 1112, 971, 2414, 1613, 753, 1875, 970, 947, 2409, 511, 1156, 786, 927, 1868, 1520, 956, 862, 864, 863, 1009, 881, 1006, 750, 1831, 1979, 960, 3945, 1688, 1083, 102, 832, 853, 762, 703, 752, 82, 1431, 1793, 2926, 865, 844, 709, 3959, 1869, 314, 112, 1160, 1134, 1044, 1309, 2347, 772, 969, 1152, 1118, 27, 984, 699, 2146, 1043, 1254, 843, 1050, 1732, 1820, 922, 887, 35, 105, 755, 972, 2430, 996, 811, 883, 693, 2240, 399, 1682, 677, 1536, 1054, 897, 1977, 26, 363, 33, 875, 1162, 316, 775, 951, 974, 886, 546, 719, 1371, 968, 861, 1897, 757, 1368, 1037, 1089, 1982, 1533, 894, 828, 20, 3930, 842, 584, 2000, 28, 825, 1129, 1026, 556, 858, 1153, 1354, 1984, 1005, 1088, 1735, 1109, 1697, 1878, 1544, 1480, 744, 1072, 3833, 574, 723, 1914, 1950, 871, 385, 1833, 2271, 394, 41, 2153, 976, 1095, 2269, 916, 3163, 796, 914, 1130, 808, 1113, 1746, 1121, 1048, 387, 611, 733, 93, 36, 1524, 131, 1087, 906, 1607, 663, 356, 849, 3322, 1080, 1073, 1007, 771, 1635, 730, 1333, 4007, 1634, 3845, 1876, 765, 751, 816, 2669, 1332, 689, 1018, 2160, 2298, 1020, 2190, 155, 1617, 3500, 1832, 359, 1873, 2310, 683, 16, 543, 889, 912, 1173, 1771, 1136, 2272, 137, 803, 1149, 307, 1954, 1320, 1855, 993, 1172, 1171, 3145, 1812, 3042, 522, 3595, 1816, 47, 1060, 366, 1066, 870, 788, 147, 1792, 1296, 3076, 312, 566, 1839, 3317, 1164, 332, 2013, 725, 3768, 2257, 25, 3078, 789, 1348, 2513, 1597, 710, 967, 1303, 2144, 1895, 1872, 1148, 1739, 375, 856, 15, 3072, 391, 354, 2168, 348, 3620, 1146, 1092, 236, 1057, 900, 2592, 1304, 1103, 653, 1645, 933, 1013, 342, 3249, 362, 487, 325, 4020, 854, 561, 950, 671, 860, 841, 4010, 911, 802, 1308, 559, 3140, 3824, 1985, 344, 913, 926, 1339, 918, 830, 1552, 735, 4088, 1042, 1583, 1710, 459, 3437, 302, 1998, 1104};

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

AV2Predictor::AV2Predictor() {
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

AV2Predictor::~AV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool AV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double AV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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

void AV2Predictor::predictBatch(
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

double AV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void AV2Predictor::predictFromSMILESBatch(
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

void AV2Predictor::extractSMARTSFeatures(
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

void AV2Predictor::extractRDKitFeatures(
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

void AV2Predictor::extractOsmordredFeatures(
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

void AV2Predictor::extractUnifiedFeatures(
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

void AV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void AV2Predictor::selectFeatures(
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
