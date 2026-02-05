#include "polarizability_v2.h"
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
static const int FEATURE_MASK[400] = {422, 2151, 1251, 1385, 1225, 323, 300, 397, 1804, 3837, 533, 3790, 1901, 523, 3786, 78, 351, 84, 3870, 588, 3807, 337, 2180, 293, 4045, 3839, 291, 317, 322, 318, 631, 3777, 215, 292, 3797, 3771, 1173, 1201, 1902, 3814, 1176, 3770, 3772, 2114, 624, 1342, 1200, 625, 2179, 596, 589, 1332, 626, 3849, 2008, 1517, 3822, 1202, 325, 1932, 3844, 3866, 577, 3819, 1204, 3852, 1341, 1240, 1829, 313, 210, 1357, 2006, 654, 1817, 1545, 1267, 1363, 80, 2001, 3776, 133, 658, 1629, 160, 927, 1239, 786, 2212, 1243, 996, 2214, 335, 1269, 1069, 1182, 1104, 1890, 630, 938, 159, 1987, 3796, 585, 570, 582, 914, 772, 2000, 2002, 595, 341, 1331, 723, 1305, 419, 627, 514, 3836, 1265, 831, 405, 741, 372, 1807, 3824, 297, 3136, 1373, 102, 1126, 3805, 2162, 413, 2771, 1998, 2005, 591, 334, 724, 3863, 1891, 1997, 336, 1352, 1963, 1535, 1976, 694, 594, 795, 1999, 695, 1874, 1266, 2406, 2216, 365, 3826, 2415, 815, 2003, 865, 667, 777, 1377, 1828, 1328, 1304, 388, 137, 1845, 756, 821, 1101, 2249, 1028, 1359, 2099, 1106, 1257, 1885, 978, 3144, 1351, 731, 2502, 858, 2269, 974, 181, 735, 1084, 616, 1896, 2414, 778, 1815, 1961, 35, 3795, 1371, 521, 561, 872, 429, 615, 926, 1372, 653, 1384, 860, 1747, 726, 1226, 1837, 1074, 1296, 2516, 866, 767, 3769, 875, 2499, 629, 571, 837, 774, 784, 824, 1124, 1322, 1065, 3501, 1211, 655, 1826, 1818, 1809, 1355, 3846, 1820, 1013, 874, 1141, 1144, 1217, 1930, 966, 73, 226, 1996, 587, 1029, 711, 1044, 389, 348, 2050, 359, 1720, 649, 315, 1128, 516, 548, 1149, 901, 1334, 1898, 1335, 917, 1059, 1172, 221, 2409, 1714, 19, 1416, 3851, 3834, 344, 1836, 1213, 3780, 869, 1010, 1053, 1350, 2229, 1133, 340, 3232, 752, 593, 2088, 1356, 2254, 1864, 1378, 1310, 761, 1694, 2297, 3612, 3798, 1979, 1152, 1230, 814, 3804, 2061, 2230, 513, 1827, 1156, 3775, 850, 1816, 947, 27, 576, 1043, 1330, 652, 979, 1070, 299, 1895, 572, 524, 370, 3135, 704, 633, 1967, 1321, 818, 3399, 1838, 1802, 2009, 688, 1516, 379, 1073, 1047, 298, 3821, 1994, 717, 1012, 599, 1135, 939, 551, 1025, 1808, 1297, 859, 995, 1094, 2156, 392, 2257, 1812, 728, 1850, 3504, 1799, 751, 1011, 760, 1153, 2867, 1289, 1160, 661, 1858, 378, 700, 920, 1143, 1186, 1844, 1288, 692, 1695, 1612, 1245, 1241, 581, 812, 1862, 1081, 522, 997, 2026, 1966, 1337};

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

PolarizabilityV2Predictor::PolarizabilityV2Predictor() {
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

PolarizabilityV2Predictor::~PolarizabilityV2Predictor() {
    delete desc_calc_;
    // Clean up SMARTS queries
    for (auto* mol : smarts_queries_) {
        delete mol;
    }
}

bool PolarizabilityV2Predictor::isMoleculeTooLarge(const RDKit::ROMol* mol) {
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

double PolarizabilityV2Predictor::predict(const RDKit::ROMol* mol, int numThreads) const {
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

void PolarizabilityV2Predictor::predictBatch(
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

double PolarizabilityV2Predictor::predictFromSMILES(const std::string& smiles, int numThreads) const {
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

void PolarizabilityV2Predictor::predictFromSMILESBatch(
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

void PolarizabilityV2Predictor::extractSMARTSFeatures(
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

void PolarizabilityV2Predictor::extractRDKitFeatures(
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

void PolarizabilityV2Predictor::extractOsmordredFeatures(
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

void PolarizabilityV2Predictor::extractUnifiedFeatures(
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

void PolarizabilityV2Predictor::preprocessFeatures(std::vector<double>& features) const {
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

void PolarizabilityV2Predictor::selectFeatures(
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
