#include "e_v2.h"
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
static const int FEATURE_MASK[600] = {2006, 1285, 3769, 2001, 1278, 1998, 2000, 1597, 1397, 2779, 1206, 353, 1811, 1593, 1133, 2005, 1825, 1155, 3778, 114, 1537, 319, 544, 420, 2205, 1386, 1375, 3783, 1373, 550, 95, 399, 1330, 1372, 411, 732, 1977, 1376, 1847, 375, 1271, 1849, 405, 4044, 2412, 1846, 2198, 2153, 581, 2201, 832, 1170, 1219, 676, 1957, 1245, 3034, 571, 1165, 412, 1886, 2152, 1808, 1101, 1981, 1999, 1856, 126, 1842, 2022, 1850, 2669, 407, 2002, 4019, 2777, 1512, 2254, 580, 1037, 921, 1045, 34, 1983, 295, 2427, 1531, 1986, 726, 3785, 3768, 1805, 2046, 1859, 1698, 314, 366, 3807, 3824, 1619, 2007, 2771, 3842, 2162, 1180, 2233, 294, 299, 1964, 1888, 384, 3829, 3952, 3797, 346, 1384, 317, 389, 922, 410, 428, 563, 1599, 1877, 1800, 1733, 1987, 2297, 534, 3403, 184, 2271, 81, 1907, 1875, 1159, 518, 652, 3897, 3401, 2150, 115, 374, 1517, 2272, 1867, 1890, 2200, 2155, 1996, 672, 2074, 3770, 904, 1852, 1022, 424, 860, 1860, 1990, 1173, 572, 1310, 398, 561, 1857, 650, 372, 930, 1518, 2199, 388, 1851, 1953, 757, 1232, 2405, 966, 339, 421, 551, 654, 1844, 298, 1809, 1992, 545, 1186, 1309, 1911, 1843, 1813, 727, 355, 733, 1430, 1065, 409, 1134, 833, 3848, 464, 1167, 1694, 343, 1699, 1549, 352, 150, 2744, 63, 1357, 1947, 584, 2014, 1926, 649, 1163, 327, 3202, 570, 1905, 396, 864, 594, 1321, 3790, 1069, 336, 386, 2257, 1609, 362, 1127, 609, 1377, 691, 796, 393, 1322, 704, 1812, 664, 573, 1151, 1312, 842, 1861, 814, 4023, 999, 1508, 367, 524, 1654, 1161, 768, 1930, 1708, 1365, 2009, 1835, 438, 1814, 974, 2407, 841, 806, 1288, 940, 1973, 3796, 2185, 1919, 1946, 859, 1898, 219, 4092, 1991, 1533, 815, 122, 3802, 1076, 2274, 914, 3766, 1005, 1710, 1368, 2709, 624, 1838, 1068, 1146, 4000, 3808, 670, 1821, 1289, 1350, 1939, 2252, 2379, 3832, 1872, 776, 2026, 2197, 2275, 1014, 715, 1906, 780, 685, 394, 1613, 349, 2195, 1855, 3805, 1845, 392, 2202, 1331, 3780, 1960, 3776, 1473, 1510, 401, 452, 837, 957, 769, 338, 3142, 1854, 107, 1446, 1313, 678, 1092, 3256, 965, 1348, 1972, 1636, 1304, 523, 1944, 1853, 2157, 1589, 1299, 1516, 130, 671, 2838, 1282, 677, 1293, 1256, 3469, 797, 2671, 1713, 1128, 1160, 920, 1547, 2251, 579, 3798, 2268, 998, 1994, 3474, 1948, 2306, 970, 632, 1617, 648, 1093, 1985, 1976, 1334, 160, 725, 1892, 1091, 658, 1982, 2298, 1061, 543, 370, 2070, 1524, 1629, 2864, 379, 979, 1893, 1171, 144, 659, 2004, 805, 2512, 385, 1927, 967, 802, 874, 1057, 3908, 4041, 1997, 1052, 1009, 741, 3254, 2608, 500, 221, 656, 3040, 1415, 2262, 109, 3399, 296, 1164, 2249, 919, 1914, 1319, 840, 816, 1078, 2269, 3609, 990, 1305, 1169, 642, 1887, 2168, 870, 791, 956, 1995, 982, 872, 1287, 736, 302, 1108, 1193, 1141, 108, 1172, 695, 3104, 819, 1157, 2180, 2684, 2178, 681, 1594, 124, 867, 414, 936, 1051, 1031, 1730, 1918, 1074, 3835, 912, 1158, 858, 1868, 1724, 1625, 716, 582, 179, 3948, 1611, 1041, 2186, 1963, 475, 2234, 913, 514, 1011, 1347, 2208, 786, 1088, 2037, 2242, 2973, 1841, 1581, 1796, 2231, 933, 3812, 948, 3939, 3109, 756, 1535, 1978, 667, 788, 1144, 923, 1934, 730, 1485, 1536, 363, 365, 516, 879, 1142, 907, 368, 3036, 655, 926, 869, 1071, 1833, 2261, 703, 315, 1371, 881, 3765, 960, 1385, 1600, 1967, 893, 690, 895, 824, 2230, 2179, 1130, 1407, 1869, 1124, 1254, 2355, 1712, 1055, 831, 2303, 1102, 971, 1839, 762, 595, 773, 2216, 751, 369, 1335, 358, 1637, 1176, 159, 790, 1117, 1824, 3511, 753, 1332, 3594, 3868, 1337, 441, 1029, 823, 1828, 449, 1817, 1303, 3607};

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

EV2Predictor::EV2Predictor() {
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

EV2Predictor::~EV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool EV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double EV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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

void EV2Predictor::predictBatch(
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

double EV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void EV2Predictor::predictFromSMILESBatch(
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

void EV2Predictor::extractSMARTSFeatures(
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

void EV2Predictor::extractRDKitFeatures(
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

void EV2Predictor::extractOsmordredFeatures(
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

void EV2Predictor::extractUnifiedFeatures(
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

void EV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void EV2Predictor::selectFeatures(
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
