#include "logws_v2.h"
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
static const int FEATURE_MASK[500] = {2150, 421, 1853, 1278, 1998, 1997, 1996, 2007, 531, 319, 2003, 1331, 1377, 2185, 546, 1332, 730, 474, 1809, 1963, 332, 2195, 1819, 2002, 1989, 511, 691, 544, 646, 3769, 2198, 710, 586, 1288, 3952, 1159, 712, 806, 1918, 1093, 922, 682, 692, 3935, 3776, 1991, 631, 1307, 1308, 1021, 1125, 1157, 903, 2836, 652, 628, 2001, 369, 673, 694, 1729, 1061, 3797, 20, 1152, 1999, 2525, 1351, 1161, 1954, 1527, 795, 721, 1545, 1151, 3808, 2157, 1978, 355, 1585, 602, 1141, 2297, 3442, 1299, 1854, 3767, 1245, 1250, 2271, 2183, 1926, 1683, 2669, 353, 823, 422, 591, 508, 656, 395, 2006, 1628, 1384, 657, 226, 1202, 3183, 2151, 822, 1029, 2155, 1891, 1154, 914, 677, 1986, 2379, 4049, 650, 1117, 1177, 2257, 787, 1872, 1611, 1621, 690, 1994, 1700, 1510, 344, 1294, 1232, 1820, 2026, 1529, 1581, 3109, 729, 434, 1646, 1373, 1533, 658, 2210, 3820, 1884, 584, 676, 877, 1302, 1534, 654, 1062, 751, 1110, 931, 1626, 948, 1925, 400, 3765, 2156, 2176, 3996, 3139, 1892, 1346, 2268, 1635, 1356, 455, 1939, 1133, 3780, 613, 675, 399, 2254, 372, 1596, 769, 2471, 1966, 1077, 1636, 1269, 73, 374, 1738, 930, 1071, 2252, 2526, 2303, 683, 1484, 1347, 859, 1381, 116, 723, 1968, 1258, 1267, 705, 3870, 525, 1416, 727, 428, 581, 693, 939, 1201, 1037, 1846, 2251, 1129, 1162, 1995, 1984, 2618, 1722, 997, 1990, 3089, 1712, 2161, 3764, 1747, 3766, 1630, 860, 589, 1335, 1967, 1836, 1369, 1921, 4079, 1253, 704, 1150, 913, 134, 1932, 674, 2542, 2205, 3842, 2014, 1890, 998, 1844, 3802, 1118, 804, 1993, 1014, 949, 1171, 2298, 1158, 695, 3034, 731, 2000, 1085, 1552, 681, 596, 1812, 1796, 1293, 2274, 1357, 1054, 724, 912, 363, 1156, 689, 1946, 1808, 225, 2409, 734, 1838, 659, 685, 4011, 805, 2152, 813, 849, 607, 1298, 2272, 3399, 2153, 2004, 1981, 1138, 1843, 2099, 3831, 1957, 530, 1863, 688, 2009, 1494, 1517, 1697, 1845, 884, 547, 1825, 4092, 1977, 1271, 2256, 635, 2208, 1979, 368, 1730, 643, 2125, 409, 1587, 841, 1319, 735, 1629, 522, 604, 1350, 1649, 3824, 1397, 107, 1149, 707, 2234, 745, 878, 2232, 1548, 1525, 154, 3256, 966, 2089, 1072, 1079, 1313, 3786, 937, 2196, 305, 714, 3566, 1881, 1713, 988, 904, 2587, 687, 1444, 1176, 577, 2187, 315, 1682, 1616, 637, 235, 1180, 348, 1063, 588, 3038, 1992, 1754, 1982, 295, 1287, 3504, 1831, 1123, 832, 1511, 1084, 1363, 1589, 1866, 2472, 1615, 838, 1696, 2377, 366, 728, 3785, 1910, 1710, 1976, 944, 1372, 645, 1908, 1005, 4050, 1303, 1080, 647, 1699, 1354, 112, 1371, 1001, 882, 1883, 1943, 3242, 2133, 1594, 1311, 600, 1045, 1841, 2093, 1119, 1078, 1514, 1632, 521, 1172, 2190, 1842, 564, 1164, 1800, 427, 563, 1322, 583, 3796, 2489, 1949, 843, 524, 3857, 932, 350, 3803, 1789, 341, 1508, 1219, 3955, 1353, 1983, 1871, 940, 4020, 1653, 496, 1349, 1121, 1031, 1120, 3936, 3074, 3621, 3836, 2103, 890, 3439, 560, 651, 1920, 1082, 3781, 367, 2501, 3812, 1290, 362, 1348, 982, 2499, 1969, 327, 1169, 227, 1595};

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

logWSV2Predictor::logWSV2Predictor() {
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

logWSV2Predictor::~logWSV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool logWSV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double logWSV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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
    std::vector<double> selected_features(500);
    for (int i = 0; i < 500; ++i) {
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

void logWSV2Predictor::predictBatch(
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

double logWSV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void logWSV2Predictor::predictFromSMILESBatch(
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

void logWSV2Predictor::extractSMARTSFeatures(
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

void logWSV2Predictor::extractRDKitFeatures(
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

void logWSV2Predictor::extractOsmordredFeatures(
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

void logWSV2Predictor::extractUnifiedFeatures(
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

void logWSV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void logWSV2Predictor::selectFeatures(
    const std::vector<double>& all_features,
    std::vector<double>& selected
) const {
    // Use FEATURE_MASK from header
    for (int i = 0; i < 500; ++i) {
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
