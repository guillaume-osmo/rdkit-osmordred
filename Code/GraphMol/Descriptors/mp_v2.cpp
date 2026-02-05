#include "mp_v2.h"
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
static const int FEATURE_MASK[600] = {2152, 374, 1329, 2185, 3768, 2000, 653, 3776, 2153, 1852, 3765, 2274, 1989, 3769, 420, 410, 2007, 3819, 1334, 1992, 1925, 2001, 2002, 1998, 1999, 1331, 412, 571, 2196, 3815, 1995, 398, 2082, 3770, 2050, 2026, 1807, 1860, 1282, 4044, 2003, 1377, 2242, 1946, 3833, 1926, 534, 4023, 1810, 1981, 2669, 2136, 1170, 439, 3829, 2149, 1344, 427, 1021, 1497, 126, 511, 1483, 2275, 386, 2093, 2271, 232, 2231, 1335, 2148, 546, 1151, 2186, 1496, 1855, 1305, 388, 2257, 1782, 651, 573, 1994, 643, 438, 428, 151, 2472, 305, 362, 1986, 4040, 867, 1831, 1646, 1340, 1904, 1156, 689, 2983, 1972, 2302, 1692, 343, 413, 1837, 1320, 389, 1319, 1510, 329, 1515, 520, 1967, 1407, 302, 463, 633, 1918, 2262, 1621, 1815, 1348, 743, 1324, 1839, 1963, 1310, 1947, 1525, 851, 3798, 363, 3109, 1528, 1347, 652, 1078, 2240, 2189, 1029, 150, 3838, 1480, 1339, 1948, 1478, 2063, 358, 2228, 3821, 4049, 2298, 366, 1101, 2237, 2004, 1996, 2015, 1932, 2268, 3034, 1611, 3766, 2251, 1598, 518, 912, 2260, 357, 101, 587, 2154, 2184, 1593, 1893, 487, 1414, 233, 4092, 315, 1005, 581, 989, 316, 2187, 353, 423, 1846, 1165, 2259, 1370, 931, 1891, 3845, 1977, 3245, 2014, 1984, 813, 2234, 1372, 2770, 2253, 2379, 1150, 1894, 64, 1920, 1163, 1798, 176, 308, 2453, 3474, 649, 2480, 2235, 1813, 2273, 967, 1311, 1966, 563, 2198, 3973, 1785, 1133, 3877, 459, 1591, 859, 1580, 299, 470, 1982, 1960, 3713, 3500, 372, 1625, 307, 725, 1859, 2178, 582, 2693, 872, 2005, 1579, 3897, 3348, 345, 1596, 644, 1783, 1682, 3610, 2297, 1883, 2062, 1330, 1616, 1054, 886, 1858, 2864, 1877, 1584, 1512, 3865, 400, 786, 1535, 585, 434, 1069, 621, 2515, 662, 1524, 3959, 411, 924, 707, 2157, 4011, 1944, 2618, 1820, 2248, 1312, 726, 1735, 228, 2372, 1154, 182, 348, 2006, 1680, 642, 1968, 426, 424, 1993, 708, 2161, 2168, 3575, 671, 402, 2201, 2276, 1296, 1295, 1219, 1854, 489, 1171, 2195, 1545, 3767, 1976, 677, 631, 1110, 654, 2877, 824, 904, 1788, 905, 1871, 1908, 742, 1258, 2211, 1717, 83, 1895, 1921, 63, 1816, 712, 1142, 777, 705, 1896, 1943, 1695, 1648, 806, 107, 2965, 1322, 356, 1317, 355, 2600, 1876, 691, 1850, 3803, 753, 3108, 3994, 1864, 469, 1583, 1513, 2416, 1892, 1118, 778, 1157, 1970, 686, 1786, 2512, 761, 1718, 407, 1581, 1168, 815, 2216, 2880, 2272, 140, 1332, 1373, 1832, 2709, 1853, 555, 1690, 2711, 1840, 1006, 2013, 1056, 418, 1875, 1086, 4068, 390, 2025, 1959, 1843, 344, 2405, 698, 1726, 957, 106, 2122, 1289, 1683, 1737, 3439, 1863, 869, 1356, 718, 1337, 574, 1520, 446, 1104, 1685, 1371, 3596, 4079, 762, 958, 1094, 3231, 1143, 796, 2743, 717, 544, 347, 529, 1045, 1037, 1881, 226, 2308, 1699, 760, 1152, 1697, 825, 690, 421, 1412, 3919, 1874, 385, 1173, 1285, 2169, 764, 1117, 615, 795, 1781, 3796, 592, 894, 2542, 1867, 1838, 791, 2499, 1619, 539, 1636, 1594, 3941, 3502, 1733, 2387, 2495, 510, 1569, 868, 2155, 1613, 1180, 680, 1109, 1013, 4036, 409, 303, 1318, 3834, 1887, 1634, 1523, 1367, 1844, 1547, 2176, 1797, 1600, 1293, 1027, 408, 180, 124, 1725, 1061, 332, 1969, 1299, 670, 1292, 1368, 225, 1072, 914, 3277, 3210, 1527, 314, 2771, 2261, 1316, 3542, 1973, 3329, 4086, 2362, 998, 833, 2370, 903, 1534, 2599, 1681, 667, 735, 2239, 1306, 895, 639, 2174, 34, 2167, 734, 1951, 1516, 901, 640, 4025, 2867, 1077, 2348, 1829, 3694, 317, 553, 720, 722, 3183, 1974, 1434, 841, 217, 3041, 709, 2781, 1869, 2247, 1391, 338, 1149, 897, 3939, 3225, 1346, 1508, 1017, 1172, 1070, 3949, 3510, 1314, 992, 612, 2267, 3230, 3170, 349, 678, 736};

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

MPV2Predictor::MPV2Predictor() {
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

MPV2Predictor::~MPV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool MPV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double MPV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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

void MPV2Predictor::predictBatch(
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

double MPV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void MPV2Predictor::predictFromSMILESBatch(
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

void MPV2Predictor::extractSMARTSFeatures(
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

void MPV2Predictor::extractRDKitFeatures(
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

void MPV2Predictor::extractOsmordredFeatures(
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

void MPV2Predictor::extractUnifiedFeatures(
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

void MPV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void MPV2Predictor::selectFeatures(
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
