#include "deltahc_v2.h"
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
static const int FEATURE_MASK[500] = {2008, 1363, 625, 831, 795, 4045, 533, 626, 217, 624, 397, 1979, 176, 369, 1853, 324, 233, 183, 232, 604, 67, 1193, 705, 376, 531, 710, 2150, 1842, 702, 1872, 1846, 365, 3776, 2406, 1180, 1377, 213, 697, 1805, 378, 524, 912, 1878, 1821, 352, 1836, 1813, 3785, 421, 706, 565, 422, 1865, 1900, 2001, 1980, 696, 876, 1866, 2006, 686, 2234, 1835, 854, 868, 1978, 885, 671, 894, 616, 1892, 298, 1313, 294, 1885, 1850, 1887, 1335, 3814, 1879, 997, 1369, 840, 588, 1806, 1960, 903, 922, 630, 3801, 3782, 1156, 2009, 574, 321, 1804, 913, 714, 1838, 850, 297, 348, 1359, 317, 78, 141, 1890, 660, 887, 921, 740, 346, 1271, 1891, 1930, 228, 1839, 931, 1616, 2179, 1368, 3869, 1355, 375, 3783, 965, 1347, 930, 1863, 631, 3501, 3868, 310, 562, 1896, 337, 1163, 628, 497, 526, 1493, 678, 1371, 1365, 614, 1021, 1520, 3835, 673, 570, 849, 741, 1884, 2164, 1844, 1241, 2203, 1029, 629, 700, 2210, 855, 704, 331, 632, 662, 1986, 783, 520, 535, 640, 1286, 79, 698, 1164, 3821, 1157, 1069, 862, 1232, 392, 751, 1889, 3822, 1166, 1961, 1321, 438, 732, 664, 3848, 3786, 787, 605, 1708, 1809, 343, 1966, 635, 805, 663, 1101, 1366, 765, 1984, 1340, 807, 716, 368, 3788, 853, 1393, 670, 2151, 1607, 1146, 1172, 882, 877, 579, 1854, 220, 516, 1141, 3764, 613, 1847, 1013, 300, 359, 652, 93, 556, 1158, 661, 1349, 1099, 814, 1843, 918, 15, 699, 35, 1145, 318, 95, 3866, 583, 895, 1346, 669, 568, 89, 2303, 722, 0, 332, 1226, 878, 1127, 1841, 329, 3781, 1214, 1717, 1888, 3780, 1350, 1333, 2779, 2169, 731, 330, 3828, 1997, 3790, 3847, 948, 179, 384, 559, 1945, 420, 595, 910, 1901, 844, 64, 3807, 1893, 555, 1824, 685, 1025, 1109, 325, 1064, 1943, 1291, 391, 1062, 770, 656, 966, 363, 593, 1548, 859, 309, 3906, 90, 86, 1087, 3829, 627, 601, 1389, 1149, 1267, 380, 2166, 728, 802, 3802, 742, 726, 1147, 690, 743, 2301, 729, 1160, 87, 41, 1343, 361, 1909, 1509, 122, 864, 599, 738, 2005, 939, 4055, 320, 1125, 760, 302, 973, 1279, 1162, 513, 370, 1364, 777, 718, 1159, 560, 347, 299, 3034, 1956, 1088, 3796, 549, 1144, 679, 147, 688, 295, 723, 667, 1129, 1834, 1154, 3805, 1351, 382, 947, 1182, 1186, 80, 1289, 140, 349, 1219, 1065, 1089, 763, 950, 694, 1515, 1957, 621, 603, 353, 315, 1136, 1169, 3800, 736, 884, 1867, 813, 644, 986, 3849, 796, 84, 322, 970, 612, 1120, 291, 434, 572, 581, 292, 724, 1615, 1203, 1110, 1299, 1117, 1311, 1514, 750, 303, 312, 768, 3850, 756, 715, 737, 328, 1050, 4004, 957, 747, 1897, 822, 1336, 2062, 730, 2201, 2490, 1281, 1118, 648, 1091, 395, 316, 319, 3824, 1046, 1832, 3867, 536, 915, 639, 3831, 389, 1547, 2248, 991, 1168, 88, 1178, 1078, 733, 636, 1053, 886, 734, 753, 12, 566, 1111, 1594, 2298, 24, 558, 1807, 684, 594, 771, 1898, 1908, 1152, 881, 2239, 1084, 778, 1148, 1372, 3858, 769, 1155, 1095};

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

deltaHcV2Predictor::deltaHcV2Predictor() {
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

deltaHcV2Predictor::~deltaHcV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool deltaHcV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double deltaHcV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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

void deltaHcV2Predictor::predictBatch(
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

double deltaHcV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void deltaHcV2Predictor::predictFromSMILESBatch(
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

void deltaHcV2Predictor::extractSMARTSFeatures(
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

void deltaHcV2Predictor::extractRDKitFeatures(
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

void deltaHcV2Predictor::extractOsmordredFeatures(
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

void deltaHcV2Predictor::extractUnifiedFeatures(
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

void deltaHcV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void deltaHcV2Predictor::selectFeatures(
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
