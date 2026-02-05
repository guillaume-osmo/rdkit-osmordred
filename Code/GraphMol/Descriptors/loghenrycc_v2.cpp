#include "loghenrycc_v2.h"
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
static const int FEATURE_MASK[600] = {374, 2152, 338, 3767, 3766, 410, 1896, 3768, 1903, 562, 1852, 1948, 2003, 2303, 613, 1530, 1736, 1165, 2004, 2153, 412, 1535, 705, 560, 2002, 3870, 2254, 4044, 3769, 1534, 2209, 1865, 1968, 578, 1162, 3077, 422, 630, 1232, 2271, 389, 1289, 1726, 1329, 398, 1981, 1850, 3765, 1524, 520, 1527, 370, 580, 1093, 3442, 559, 1163, 349, 561, 2274, 1928, 384, 796, 1167, 588, 335, 1612, 1071, 416, 1288, 3210, 1943, 1386, 1635, 1287, 2158, 726, 2022, 921, 1063, 1808, 2275, 1809, 343, 1037, 1185, 518, 1175, 1339, 2200, 957, 1076, 1506, 1898, 1095, 569, 2319, 545, 814, 673, 2211, 741, 1314, 1371, 998, 1113, 1984, 1892, 592, 303, 1121, 903, 1206, 1946, 2182, 1771, 624, 1692, 342, 1835, 898, 3776, 1369, 1372, 1839, 3142, 1005, 2164, 530, 1999, 523, 1625, 1813, 727, 1158, 1982, 940, 704, 2541, 832, 534, 2242, 1843, 1855, 1966, 1377, 1983, 1812, 706, 1375, 1894, 18, 1977, 1638, 1085, 333, 3805, 806, 1346, 587, 357, 1305, 720, 363, 2255, 3880, 1384, 2161, 372, 1926, 3237, 366, 1030, 1989, 1974, 815, 1217, 2009, 150, 1533, 2542, 1309, 1713, 1009, 781, 1374, 983, 969, 548, 621, 1316, 1964, 904, 861, 831, 958, 1170, 927, 717, 650, 771, 1959, 1245, 1198, 840, 627, 1243, 1104, 327, 393, 226, 809, 1516, 801, 914, 350, 669, 1057, 648, 770, 3833, 736, 695, 697, 1081, 1065, 3847, 2507, 851, 3764, 600, 659, 3812, 1029, 915, 896, 231, 1129, 1087, 696, 959, 3807, 1171, 948, 2669, 905, 1151, 583, 1932, 817, 1696, 421, 3399, 2174, 985, 2262, 1140, 341, 707, 1697, 358, 521, 1312, 564, 304, 662, 547, 1178, 1918, 660, 1921, 1263, 1497, 376, 217, 776, 888, 883, 709, 1607, 657, 331, 1817, 774, 1014, 2163, 574, 1281, 2082, 846, 1867, 649, 1211, 1537, 687, 2208, 1880, 1811, 812, 1967, 158, 973, 519, 321, 1101, 3808, 753, 1731, 768, 698, 1956, 2162, 699, 807, 966, 1994, 842, 1053, 1335, 1088, 777, 1331, 336, 642, 1853, 1155, 652, 855, 975, 1310, 751, 213, 308, 1823, 718, 895, 1885, 299, 782, 386, 3796, 684, 516, 731, 1117, 1595, 685, 1109, 3548, 1202, 1957, 645, 667, 305, 332, 355, 723, 792, 344, 1797, 1897, 563, 2517, 787, 298, 1854, 769, 764, 513, 369, 73, 337, 1069, 956, 1068, 744, 1066, 3816, 3784, 625, 2150, 664, 2272, 1271, 1337, 601, 1149, 2171, 1323, 3868, 2983, 612, 816, 682, 1840, 661, 825, 1511, 239, 1060, 390, 3437, 602, 875, 1180, 2168, 328, 214, 1832, 1139, 1059, 860, 1252, 1822, 863, 2410, 2480, 1302, 1821, 1985, 70, 314, 2618, 1838, 94, 1789, 608, 1152, 1256, 986, 1826, 1347, 1727, 1695, 1138, 858, 300, 1945, 779, 1875, 799, 1600, 1147, 1837, 3783, 3782, 2201, 148, 1953, 954, 1350, 1630, 643, 388, 1615, 1352, 2170, 2210, 1594, 738, 1322, 1634, 1430, 1055, 772, 886, 1636, 757, 742, 854, 4086, 1844, 1064, 1150, 950, 558, 1349, 492, 1136, 676, 3602, 401, 1807, 1889, 3770, 3829, 990, 1654, 743, 1156, 663, 441, 3933, 1434, 791, 750, 715, 3819, 1831, 1820, 1026, 952, 1330, 415, 964, 1353, 1911, 316, 658, 1737, 1048, 1717, 989, 1998, 1033, 1699, 2864, 1882, 805, 1849, 977, 824, 319, 930, 1168, 1893, 330, 859, 1036, 1086, 631, 1070, 1017, 1134, 1015, 994, 884, 1084, 1357, 894, 788, 1973, 67, 1847, 1619, 963, 1145, 1073, 1804, 692, 3815, 297, 3575, 827, 306, 1348, 1006, 603, 3567, 1828, 1646, 907, 1058, 1127, 87, 1120, 876, 656, 1142, 1810, 1204, 356, 74, 3849, 351, 900, 862, 822, 765, 315, 1004, 712, 2257, 419, 679, 1799, 734, 737, 1351, 759, 1041, 38, 887, 879, 3244, 1868, 365, 215, 345, 1930};

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

logHenryccV2Predictor::logHenryccV2Predictor() {
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

logHenryccV2Predictor::~logHenryccV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool logHenryccV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double logHenryccV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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
    std::vector<double> selected_features(600);
    for (int i = 0; i < 600; ++i) {
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

void logHenryccV2Predictor::predictBatch(
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

double logHenryccV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void logHenryccV2Predictor::predictFromSMILESBatch(
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

void logHenryccV2Predictor::extractSMARTSFeatures(
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

void logHenryccV2Predictor::extractRDKitFeatures(
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

void logHenryccV2Predictor::extractOsmordredFeatures(
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

void logHenryccV2Predictor::extractUnifiedFeatures(
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

void logHenryccV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void logHenryccV2Predictor::selectFeatures(
    const std::vector<double>& all_features,
    std::vector<double>& selected
) const {
    // Use FEATURE_MASK from header
    for (int i = 0; i < 600; ++i) {
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
