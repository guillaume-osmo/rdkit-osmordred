#include "dipolemoment_v2.h"
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
static const int FEATURE_MASK[400] = {2152, 1534, 1635, 1796, 374, 1842, 1854, 1525, 664, 1155, 4041, 1159, 354, 480, 1695, 885, 876, 363, 423, 2234, 1219, 1955, 2203, 2005, 768, 1171, 1981, 1548, 1163, 3768, 1810, 1838, 1232, 2153, 563, 903, 1912, 1365, 428, 1815, 352, 1157, 3769, 534, 531, 1866, 1619, 1911, 1178, 1164, 1807, 1982, 1718, 1334, 1126, 1717, 3833, 1977, 827, 2156, 894, 362, 1886, 3770, 339, 693, 922, 1872, 303, 1308, 1812, 1170, 1077, 90, 1626, 1345, 948, 435, 76, 1349, 2009, 1305, 1367, 3399, 1966, 1085, 1581, 721, 1817, 670, 1256, 1409, 1904, 773, 3812, 1843, 882, 353, 1093, 3774, 1895, 184, 635, 1323, 211, 802, 1884, 1070, 639, 496, 102, 750, 688, 814, 1206, 2248, 665, 1125, 741, 1960, 839, 3767, 1091, 823, 1061, 1856, 3800, 1332, 636, 226, 1150, 1391, 809, 304, 1162, 315, 2164, 390, 1333, 950, 1145, 1128, 1112, 1700, 682, 3845, 961, 999, 1897, 640, 383, 3829, 669, 1366, 1529, 1069, 1825, 351, 921, 3811, 3849, 671, 3835, 763, 702, 2216, 748, 780, 3471, 962, 1926, 913, 1035, 3858, 1119, 1607, 1599, 1110, 796, 701, 547, 2201, 1355, 643, 1137, 1055, 1850, 2137, 1943, 931, 1062, 1755, 1064, 690, 1079, 2251, 683, 1165, 384, 1063, 742, 1292, 1507, 986, 1986, 562, 1172, 878, 1109, 409, 877, 993, 850, 738, 366, 1143, 883, 1221, 1193, 1098, 1054, 1925, 955, 1881, 2231, 678, 711, 332, 699, 1889, 1828, 2171, 854, 1013, 1056, 1336, 566, 779, 1140, 716, 736, 764, 1893, 3766, 2376, 676, 708, 726, 1351, 1049, 689, 915, 3106, 350, 1508, 967, 2334, 376, 3782, 657, 1104, 308, 385, 810, 800, 2711, 712, 2003, 1871, 1271, 957, 987, 3441, 1873, 2237, 677, 1053, 2669, 1021, 379, 1708, 956, 307, 325, 3076, 128, 558, 1008, 1051, 1288, 1105, 3850, 1113, 1047, 305, 663, 1289, 1057, 1180, 743, 1636, 405, 925, 817, 2260, 886, 1152, 1826, 2359, 388, 1892, 1371, 150, 412, 1372, 745, 813, 943, 965, 774, 2303, 811, 3824, 306, 880, 1322, 799, 1078, 1875, 3765, 928, 847, 770, 1910, 3846, 434, 1022, 3764, 1258, 2255, 1007, 174, 544, 983, 1887, 1086, 680, 649, 1005, 797, 395, 973, 914, 926, 803, 1654, 861, 766, 645, 1048, 2166, 933, 837, 834, 1855, 977, 1835, 801, 1136, 782, 2307, 812, 651, 963, 3805, 841, 873, 990, 1111, 370, 819, 1019, 754, 421, 1515, 3838, 355, 2177, 2365, 510, 338, 747, 365, 858, 832, 1858, 740, 668, 951, 659};

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

DipoleMomentV2Predictor::DipoleMomentV2Predictor() {
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

DipoleMomentV2Predictor::~DipoleMomentV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool DipoleMomentV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double DipoleMomentV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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
    std::vector<double> selected_features(400);
    for (int i = 0; i < 400; ++i) {
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

void DipoleMomentV2Predictor::predictBatch(
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

double DipoleMomentV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void DipoleMomentV2Predictor::predictFromSMILESBatch(
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

void DipoleMomentV2Predictor::extractSMARTSFeatures(
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

void DipoleMomentV2Predictor::extractRDKitFeatures(
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

void DipoleMomentV2Predictor::extractOsmordredFeatures(
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

void DipoleMomentV2Predictor::extractUnifiedFeatures(
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

void DipoleMomentV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void DipoleMomentV2Predictor::selectFeatures(
    const std::vector<double>& all_features,
    std::vector<double>& selected
) const {
    // Use FEATURE_MASK from header
    for (int i = 0; i < 400; ++i) {
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
