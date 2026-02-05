#include "flashpoint_v2.h"
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
static const int FEATURE_MASK[800] = {1386, 1997, 3770, 1373, 1328, 1263, 1811, 549, 2152, 551, 374, 569, 1357, 1804, 3768, 543, 1861, 410, 2153, 1375, 580, 2182, 1374, 2251, 398, 2163, 3819, 3820, 737, 2254, 1966, 2669, 603, 587, 1976, 935, 1331, 1492, 1018, 571, 2303, 701, 1157, 1860, 1803, 1389, 397, 519, 1852, 337, 565, 1369, 1993, 1512, 1409, 1988, 570, 3776, 1025, 518, 867, 411, 618, 1108, 3800, 1332, 389, 706, 1285, 2000, 823, 2274, 3765, 1253, 756, 3467, 2427, 1182, 846, 712, 1597, 1211, 2026, 2007, 1310, 534, 593, 997, 1259, 3229, 1612, 803, 813, 1994, 636, 562, 1363, 3769, 1881, 1893, 1849, 685, 853, 913, 759, 3840, 711, 3041, 946, 1809, 1529, 702, 1506, 2179, 714, 1138, 330, 1579, 624, 4086, 169, 1366, 1278, 1100, 965, 546, 2175, 298, 3786, 839, 626, 1005, 172, 512, 1101, 977, 3858, 341, 1844, 2150, 700, 1483, 1948, 456, 878, 3844, 1346, 3839, 2311, 1417, 564, 990, 1126, 709, 2415, 1619, 1349, 1909, 1969, 3442, 421, 838, 703, 1859, 591, 694, 3868, 1076, 1045, 941, 2249, 1873, 2380, 3454, 500, 226, 3094, 3774, 2004, 666, 423, 1334, 407, 857, 582, 331, 1866, 1982, 3089, 2234, 1371, 691, 434, 1161, 1636, 3040, 1312, 1833, 3826, 1914, 1172, 2228, 1020, 1073, 2162, 814, 363, 1981, 1059, 578, 4008, 2359, 3036, 1814, 798, 724, 907, 1048, 723, 224, 574, 4041, 1412, 322, 557, 1862, 2343, 338, 1961, 1038, 583, 741, 3829, 1864, 924, 1886, 1714, 185, 748, 357, 730, 905, 319, 1113, 670, 1055, 3777, 1897, 902, 1895, 976, 1189, 3405, 3407, 1843, 1853, 140, 911, 1863, 733, 1323, 1090, 3781, 316, 1019, 855, 1088, 1070, 940, 1949, 1508, 1888, 1891, 550, 2237, 658, 1951, 1180, 1903, 1939, 326, 1347, 923, 1111, 1064, 849, 3038, 665, 1830, 1139, 303, 3500, 1094, 3403, 1710, 1600, 959, 1174, 1520, 1333, 1840, 869, 1648, 1377, 1140, 4049, 1125, 954, 1826, 892, 1980, 975, 842, 1798, 1063, 3441, 1713, 1171, 563, 320, 1962, 559, 1013, 1977, 1117, 768, 1149, 980, 1060, 763, 347, 788, 1321, 3095, 390, 778, 2307, 1289, 1918, 304, 355, 1309, 422, 327, 716, 2194, 1089, 173, 1953, 1062, 1875, 3833, 183, 1224, 3767, 101, 1868, 526, 1979, 979, 1260, 2210, 679, 1037, 3406, 757, 94, 677, 1954, 585, 3766, 1614, 769, 520, 184, 949, 939, 3076, 2308, 1014, 1828, 213, 662, 3104, 1957, 753, 732, 817, 835, 1596, 810, 1327, 589, 568, 392, 4023, 856, 1593, 687, 1911, 1104, 3037, 314, 749, 802, 336, 1854, 734, 2009, 879, 707, 786, 219, 1634, 1887, 579, 649, 3042, 1694, 648, 1827, 2275, 1621, 3145, 742, 372, 1081, 515, 806, 801, 865, 1831, 796, 313, 323, 3817, 1717, 708, 1926, 1052, 2271, 610, 1818, 652, 718, 2143, 1842, 1012, 2433, 3102, 353, 2673, 1061, 2230, 812, 904, 1872, 2161, 522, 1091, 1824, 673, 684, 1252, 1114, 863, 1972, 1288, 1107, 1958, 807, 3801, 852, 1960, 870, 1807, 1370, 934, 356, 305, 301, 1822, 4044, 1837, 1812, 3401, 1889, 1152, 3053, 1074, 866, 1999, 1819, 1178, 513, 3271, 1855, 3594, 3952, 1109, 1580, 3522, 3135, 1156, 797, 3437, 1613, 1053, 1345, 1335, 725, 1007, 309, 419, 3848, 996, 1163, 1880, 1054, 1077, 885, 1058, 1829, 995, 3110, 1173, 596, 686, 908, 2272, 747, 1023, 3812, 1890, 302, 957, 894, 438, 3961, 340, 888, 659, 3870, 1286, 1882, 746, 660, 785, 1924, 771, 697, 1736, 2500, 1256, 4011, 386, 342, 1127, 296, 917, 750, 180, 545, 1046, 2616, 1150, 1894, 669, 829, 727, 1959, 1069, 792, 1547, 655, 2171, 770, 1407, 2430, 764, 418, 961, 1885, 1956, 2344, 399, 643, 841, 897, 728, 325, 1049, 1806, 1905, 2372, 1142, 1915, 1553, 962, 552, 1158, 317, 109, 704, 1680, 3832, 3764, 359, 2346, 921, 688, 2409, 300, 1291, 1339, 1787, 986, 1871, 1856, 2864, 779, 1688, 1535, 1066, 1591, 668, 3511, 830, 1793, 896, 877, 230, 385, 819, 674, 3072, 1968, 1072, 584, 2429, 4012, 950, 324, 150, 1518, 1119, 973, 315, 811, 744, 1701, 693, 2760, 1322, 2709, 535, 3922, 352, 507, 873, 887, 1147, 1497, 1867, 831, 2132, 1510, 1065, 366, 601, 960, 3136, 1154, 7, 233, 968, 1494, 2268, 1870, 816, 2250, 832, 294, 3230, 2174, 345, 1611, 232, 1095, 958, 470, 1607, 575, 1132, 1078, 1598, 657, 193, 1691, 533, 1934, 3469, 1134, 827, 805, 225, 933, 2216, 1153, 609, 909, 1735, 335, 3438, 926, 1105, 2005, 1021, 1834, 1970, 799, 1050, 646, 1415, 678, 1155, 1015, 384, 1071, 1696, 2157, 35, 868, 1133, 1057, 627, 820, 155, 1883, 1135, 361, 647, 3471, 765, 312, 3106, 953, 969, 378, 20, 83, 1434, 3402, 2501, 2877, 1820, 1513, 599, 1685, 699, 782, 972, 1162, 956, 745, 1001, 308, 217, 664, 1080, 211, 2377, 3202, 2248, 79, 952, 1771, 3107, 1515, 847, 1002, 2167, 182, 1737, 1391, 1168, 3797, 391, 1729, 1118, 561, 925, 4074, 898, 1115, 370, 729, 1791, 1876, 1281};

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

FlashpointV2Predictor::FlashpointV2Predictor() {
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

FlashpointV2Predictor::~FlashpointV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool FlashpointV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double FlashpointV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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

void FlashpointV2Predictor::predictBatch(
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

double FlashpointV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void FlashpointV2Predictor::predictFromSMILESBatch(
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

void FlashpointV2Predictor::extractSMARTSFeatures(
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

void FlashpointV2Predictor::extractRDKitFeatures(
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

void FlashpointV2Predictor::extractOsmordredFeatures(
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

void FlashpointV2Predictor::extractUnifiedFeatures(
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

void FlashpointV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void FlashpointV2Predictor::selectFeatures(
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
