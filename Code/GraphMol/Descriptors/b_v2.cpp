#include "b_v2.h"
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
static const int FEATURE_MASK[600] = {535, 3767, 3766, 1861, 338, 399, 509, 1903, 1926, 563, 1834, 374, 1851, 366, 2234, 4041, 1288, 419, 1860, 411, 2251, 598, 1856, 885, 2215, 1594, 1977, 519, 1189, 1862, 1896, 526, 2152, 1159, 896, 616, 680, 1372, 1357, 562, 595, 750, 409, 1695, 678, 534, 759, 2212, 1905, 2208, 2153, 618, 2009, 1219, 551, 543, 650, 1925, 1286, 1289, 619, 1981, 768, 1986, 2210, 1897, 566, 571, 1206, 2775, 1634, 3768, 410, 531, 1796, 822, 560, 1264, 1245, 1258, 564, 823, 1216, 568, 1080, 839, 903, 328, 561, 1699, 335, 3831, 1380, 1323, 1621, 626, 1285, 1013, 3036, 395, 3789, 638, 124, 1368, 596, 2026, 1157, 3847, 1281, 1690, 1901, 314, 362, 3830, 2301, 1029, 2926, 2271, 3034, 2272, 2252, 832, 1365, 511, 68, 3823, 704, 1349, 1352, 1318, 1791, 1829, 550, 421, 698, 2093, 737, 1371, 345, 837, 1367, 687, 1209, 2050, 1735, 3501, 701, 603, 1269, 321, 604, 683, 1174, 332, 3824, 1927, 3140, 1164, 1921, 1598, 681, 637, 3770, 2254, 1720, 1805, 974, 341, 1271, 1093, 4033, 876, 966, 1990, 975, 1852, 1534, 1858, 3764, 1079, 1552, 2216, 1589, 336, 1708, 303, 387, 344, 786, 1062, 3805, 1849, 230, 1170, 340, 1948, 867, 1038, 733, 661, 453, 642, 1800, 2303, 584, 2240, 770, 567, 2161, 3858, 741, 573, 1886, 305, 1155, 1865, 3782, 1985, 547, 1635, 1968, 367, 1204, 1282, 294, 1582, 363, 1754, 370, 1910, 1126, 1141, 398, 1611, 1609, 1872, 1370, 353, 682, 597, 794, 1866, 220, 1058, 922, 1086, 412, 580, 2669, 1863, 1063, 934, 670, 1061, 1887, 3838, 1619, 1292, 440, 760, 692, 714, 717, 1736, 912, 1094, 1296, 738, 1444, 2274, 556, 1806, 1819, 1583, 672, 1719, 1350, 87, 740, 2002, 939, 1970, 2167, 1629, 805, 365, 1531, 1414, 1839, 2005, 1430, 1906, 1843, 820, 496, 1085, 3839, 1167, 1815, 874, 2410, 226, 1153, 304, 700, 653, 2162, 702, 342, 2342, 1149, 521, 2502, 2297, 694, 3139, 3883, 679, 1256, 3850, 1879, 1710, 3796, 707, 1876, 913, 1125, 1070, 1737, 1069, 1982, 886, 1825, 1974, 2255, 1835, 815, 391, 1075, 697, 1253, 824, 1810, 4050, 3852, 390, 748, 1994, 958, 1057, 2241, 1117, 736, 372, 636, 1118, 1139, 2173, 2261, 1812, 428, 1309, 925, 841, 706, 1143, 1067, 771, 923, 1293, 2237, 1305, 967, 1036, 1867, 877, 1162, 2156, 1979, 769, 795, 933, 1816, 26, 1580, 1786, 1898, 3800, 413, 1180, 621, 735, 742, 356, 1654, 2126, 1369, 919, 1042, 1171, 3136, 1511, 930, 831, 2154, 836, 376, 2196, 2160, 914, 1163, 378, 1332, 549, 389, 1046, 1836, 1178, 921, 2260, 1528, 291, 449, 3438, 388, 3829, 2209, 917, 1333, 2025, 1313, 1840, 1516, 3808, 157, 1976, 1535, 1003, 1353, 2232, 1483, 889, 1150, 973, 3993, 384, 381, 2770, 1147, 3399, 906, 350, 669, 1078, 860, 2149, 1254, 2006, 1828, 2007, 1081, 674, 1842, 1161, 840, 727, 983, 1520, 1515, 90, 510, 3898, 1377, 1868, 1966, 751, 849, 4089, 2137, 1154, 806, 548, 2774, 520, 232, 1890, 775, 1978, 646, 875, 3939, 302, 1850, 1127, 833, 703, 364, 325, 1294, 1617, 1068, 1045, 1625, 1717, 231, 1311, 1616, 887, 407, 3769, 982, 991, 1014, 150, 3565, 789, 871, 647, 1956, 478, 610, 1054, 761, 743, 1486, 2236, 347, 355, 711, 3567, 1135, 1997, 660, 2409, 3504, 1804, 3845, 358, 438, 1134, 1108, 3877, 3137, 910, 1837, 2275, 1892, 514, 643, 819, 1321, 574, 817, 1924, 337, 880, 296, 1165, 1958, 141, 1322, 2073, 1066, 1034, 732, 1827, 2257, 639, 394, 2512, 648, 1632, 1076, 393, 4014, 725, 2256, 825, 1110, 1533, 779, 665, 811, 1912, 401, 1895, 1095, 1718, 1782, 3919, 1156, 788, 1684, 83, 1310, 2298, 3868, 753, 1106, 1346};

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

BV2Predictor::BV2Predictor() {
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

BV2Predictor::~BV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool BV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double BV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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

void BV2Predictor::predictBatch(
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

double BV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void BV2Predictor::predictFromSMILESBatch(
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

void BV2Predictor::extractSMARTSFeatures(
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

void BV2Predictor::extractRDKitFeatures(
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

void BV2Predictor::extractOsmordredFeatures(
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

void BV2Predictor::extractUnifiedFeatures(
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

void BV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void BV2Predictor::selectFeatures(
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
