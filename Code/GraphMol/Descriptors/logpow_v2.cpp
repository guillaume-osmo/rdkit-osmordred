#include "logpow_v2.h"
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
static const int FEATURE_MASK[600] = {2150, 421, 2491, 1680, 1781, 1998, 2518, 315, 1579, 562, 1368, 1374, 1908, 1921, 939, 1386, 1478, 2686, 1904, 1798, 2321, 833, 423, 1697, 348, 4079, 1545, 1278, 432, 2866, 1695, 228, 895, 997, 506, 1296, 1823, 682, 1997, 3768, 580, 3847, 1892, 660, 1021, 789, 1843, 2201, 3609, 3897, 3036, 132, 1698, 1548, 1926, 989, 4049, 571, 1907, 1486, 1809, 1692, 3766, 3625, 1854, 218, 3823, 497, 309, 3767, 886, 696, 496, 120, 3998, 1817, 464, 1431, 1906, 3782, 1591, 714, 670, 2771, 1890, 1095, 627, 2671, 1609, 1977, 470, 1619, 1943, 1847, 469, 2231, 1720, 1891, 934, 3868, 1980, 1171, 3815, 1793, 211, 81, 452, 94, 1799, 1649, 1796, 3437, 64, 1872, 1897, 90, 778, 673, 876, 614, 638, 2136, 4062, 159, 3852, 346, 3401, 1827, 2014, 1152, 823, 2471, 473, 2007, 2775, 3764, 1967, 832, 1331, 787, 869, 1530, 1930, 1795, 2002, 183, 2272, 1844, 3895, 1646, 940, 1974, 1587, 727, 2342, 511, 2006, 1927, 974, 585, 1029, 1168, 412, 2430, 3139, 1596, 921, 1324, 1156, 1690, 2252, 3829, 2410, 2344, 782, 409, 1039, 753, 2004, 2501, 826, 1158, 1155, 3914, 1356, 931, 1820, 1607, 434, 2212, 684, 1533, 1925, 1150, 344, 861, 3482, 1285, 523, 1600, 403, 1057, 3034, 400, 1797, 763, 1982, 537, 134, 1821, 1783, 107, 3886, 603, 1846, 1166, 304, 998, 2312, 2442, 216, 2774, 1013, 1819, 2299, 1932, 3613, 691, 922, 1170, 715, 1776, 831, 1316, 1684, 1484, 2236, 155, 1597, 39, 1825, 1133, 637, 316, 1611, 128, 2260, 1169, 694, 1586, 889, 3948, 2190, 1076, 2110, 2445, 1047, 3624, 990, 384, 967, 1786, 415, 860, 859, 3803, 1167, 805, 752, 374, 1141, 1911, 1160, 2232, 303, 1151, 1287, 1217, 1717, 574, 3256, 2199, 2153, 4058, 1822, 1873, 1161, 1086, 1162, 1910, 2529, 2000, 328, 2346, 13, 3348, 3102, 658, 1630, 1855, 1254, 693, 2061, 1841, 825, 151, 932, 1063, 1350, 2711, 4020, 4040, 1134, 372, 2163, 1528, 2003, 362, 1292, 769, 3567, 1328, 761, 3252, 4092, 1984, 628, 1922, 1616, 3140, 1754, 1415, 1213, 1142, 2306, 704, 1701, 2152, 1866, 521, 91, 1347, 55, 544, 3202, 966, 1064, 1788, 3955, 977, 350, 226, 1046, 545, 991, 2877, 1850, 1494, 733, 1840, 1136, 2149, 1583, 4023, 2522, 364, 1747, 1953, 816, 1968, 1372, 1176, 1973, 1022, 1412, 332, 1806, 2983, 2669, 1853, 367, 1289, 814, 572, 3888, 722, 664, 1323, 3845, 1093, 1172, 1097, 1299, 1889, 1879, 1288, 984, 1729, 79, 112, 386, 3940, 3509, 1818, 438, 427, 731, 721, 558, 1045, 2517, 311, 1828, 1628, 340, 554, 1696, 2274, 354, 2409, 3504, 937, 1060, 1691, 1787, 2387, 1978, 2512, 3441, 792, 3076, 973, 1715, 420, 1963, 2157, 3942, 2271, 1688, 1534, 1365, 1981, 43, 730, 1898, 1919, 232, 1087, 3977, 2427, 1313, 343, 230, 2090, 717, 1683, 368, 1791, 1106, 1753, 712, 1164, 1109, 948, 1839, 692, 1310, 2310, 1079, 2161, 2009, 347, 3210, 361, 1116, 1015, 1000, 1311, 1159, 735, 709, 1007, 212, 685, 388, 1962, 398, 1321, 2251, 949, 880, 1893, 1219, 3568, 1143, 725, 555, 1148, 771, 4015, 1006, 1157, 841, 589, 871, 965, 2837, 2298, 449, 960, 1367, 1950, 1725, 1206, 2026, 686, 1508, 1550, 1190, 695, 1800, 1074, 678, 734, 140, 2130, 768, 699, 1309, 522, 1895, 3330, 2013, 1085, 1960, 3072, 2169, 1071, 2297, 1510, 1483, 492, 2005, 770, 1946, 1632, 2197, 1993, 1650, 1517, 1585, 979, 877, 653, 635, 1687, 1634, 3983, 1551, 983, 899, 656, 1883, 1101, 1736, 395, 73, 1322, 982, 3780, 798, 1083, 3781, 751, 141, 3798, 4087, 655, 181, 1131, 3203, 1267, 1497, 624, 1880, 127, 102, 317, 772, 743, 887, 850, 1513, 2378, 1802, 942, 1291, 2882, 2165, 1653};

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

logPowV2Predictor::logPowV2Predictor() {
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

logPowV2Predictor::~logPowV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool logPowV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double logPowV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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

void logPowV2Predictor::predictBatch(
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

double logPowV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void logPowV2Predictor::predictFromSMILESBatch(
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

void logPowV2Predictor::extractSMARTSFeatures(
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

void logPowV2Predictor::extractRDKitFeatures(
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

void logPowV2Predictor::extractOsmordredFeatures(
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

void logPowV2Predictor::extractUnifiedFeatures(
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

void logPowV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void logPowV2Predictor::selectFeatures(
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
