#include "l_v2.h"
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
static const int FEATURE_MASK[500] = {515, 512, 1396, 1804, 3787, 3784, 1331, 1225, 3770, 1262, 3778, 1275, 1280, 527, 508, 3863, 1210, 1375, 1276, 1975, 397, 3789, 533, 1821, 1860, 1393, 1184, 3817, 1863, 1197, 517, 3783, 1211, 1249, 552, 2151, 2155, 1182, 548, 1811, 3777, 813, 1380, 411, 2213, 3869, 2003, 3814, 571, 584, 3844, 1357, 3780, 1259, 611, 1997, 337, 3819, 1995, 422, 1386, 1803, 578, 1263, 1861, 781, 1327, 2541, 89, 2530, 95, 3816, 1397, 557, 1372, 759, 323, 2851, 300, 1373, 298, 595, 1312, 580, 591, 550, 1858, 2191, 399, 520, 1374, 2312, 1825, 297, 2413, 3794, 70, 2160, 1316, 619, 795, 1718, 484, 2062, 1859, 2194, 1959, 658, 592, 1985, 1816, 78, 3839, 1320, 1973, 299, 799, 638, 398, 649, 913, 3779, 3769, 1890, 1161, 4018, 294, 1843, 1098, 1104, 3804, 587, 1384, 1809, 885, 1986, 1342, 1269, 804, 3768, 1600, 719, 374, 590, 621, 562, 805, 409, 220, 1228, 822, 2152, 1636, 352, 1348, 2274, 2007, 1377, 887, 2713, 3034, 551, 153, 145, 1328, 1534, 990, 669, 1383, 518, 1242, 2234, 806, 2561, 1527, 1972, 1827, 534, 586, 1154, 627, 1434, 628, 914, 2161, 1215, 1343, 1071, 1628, 2254, 365, 1936, 393, 3832, 3820, 1533, 1341, 315, 1146, 576, 1832, 701, 66, 1981, 2778, 385, 3249, 546, 1839, 1014, 126, 2261, 588, 1976, 1113, 1698, 1022, 2486, 1846, 1948, 1117, 1892, 1905, 3809, 1901, 1518, 1282, 727, 3764, 603, 1852, 356, 114, 1998, 684, 2153, 868, 1417, 1170, 714, 1017, 3825, 815, 1968, 582, 1813, 1919, 171, 86, 1040, 904, 509, 826, 1983, 1926, 3216, 636, 1525, 1065, 1289, 660, 1808, 671, 589, 2089, 389, 2298, 1864, 353, 969, 4004, 4025, 1285, 1515, 1356, 302, 923, 3765, 60, 499, 840, 931, 1208, 841, 2271, 3788, 710, 1545, 659, 1039, 1980, 706, 3836, 1030, 1918, 858, 1180, 788, 814, 585, 2208, 73, 2836, 912, 1694, 438, 1085, 1305, 384, 570, 624, 305, 1151, 2165, 1165, 930, 1896, 2471, 556, 2473, 707, 711, 1967, 2255, 886, 2004, 1355, 386, 1021, 421, 3798, 1635, 1426, 1033, 1005, 1189, 1594, 1535, 516, 1871, 1829, 2709, 1125, 1278, 1878, 1739, 343, 579, 2519, 735, 3850, 1335, 1147, 882, 2303, 1309, 1319, 939, 1044, 1982, 568, 326, 1061, 686, 1271, 74, 1332, 632, 3501, 667, 569, 1946, 622, 2272, 3868, 833, 1800, 72, 1063, 1131, 355, 1045, 375, 779, 1172, 863, 1363, 1129, 823, 1805, 629, 1836, 1616, 362, 593, 1166, 1865, 1156, 2864, 112, 349, 712, 1552, 674, 3805, 2237, 464, 583, 543, 1962, 1119, 2240, 1313, 3766, 3833, 222, 2500, 388, 1160, 859, 1153, 94, 1121, 637, 3219, 3846, 325, 725, 4023, 3074, 97, 685, 1089, 1508, 1150, 1853, 1979, 3848, 1334, 1617, 2164, 333, 786, 1953, 1057, 1066, 1067, 827, 1023, 1351, 958, 338, 1717, 860, 2173, 1011, 2302, 761, 1060, 3800, 2179, 2275, 3840, 1830, 932, 1841, 1597, 1072, 894, 3781, 3399, 1735, 956, 310, 903, 1059, 1163, 1019, 1126, 681, 1615, 832, 941, 3581, 791, 1287, 1497, 1158, 149, 2673, 1714, 634, 665, 1862, 1619, 1062, 2199, 3200, 974, 743};

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

LV2Predictor::LV2Predictor() {
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

LV2Predictor::~LV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool LV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double LV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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

void LV2Predictor::predictBatch(
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

double LV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void LV2Predictor::predictFromSMILESBatch(
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

void LV2Predictor::extractSMARTSFeatures(
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

void LV2Predictor::extractRDKitFeatures(
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

void LV2Predictor::extractOsmordredFeatures(
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

void LV2Predictor::extractUnifiedFeatures(
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

void LV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void LV2Predictor::selectFeatures(
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
