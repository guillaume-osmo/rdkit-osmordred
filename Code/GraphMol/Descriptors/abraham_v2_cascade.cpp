// Abraham V2 Cascade Integration
// Uses CASCADE order: V → L → E → B → S → A
// Each model uses predictions from previous models as LSFER golden features

#include "Osmordred.h"
#include "AbrahamV2Cascade.h"
#include "AbrahamV2CascadeS.h"
#include "AbrahamV2CascadeA.h"
#include <GraphMol/ROMol.h>
#include <vector>
#include <cmath>
#include <algorithm>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// Calculate Abraham V2 parameters using cascade approach
// Returns: [A, B, S, E, L, V] in standard order
std::vector<double> calcAbrahamsV2Cascade(const ROMol& mol) {
    // Get 291 base features from CalcAbrahamsV2Features
    std::vector<double> base291 = calcAbrahamsV2Features(mol);
    
    if (base291.size() != 291) {
        return std::vector<double>(6, 0.0);
    }
    
    // CASCADE ORDER: V → L → E → B → S → A
    std::vector<double> preds;
    preds.reserve(6);
    
    // V (291 features only)
    double V = AbrahamV2Cascade::VModel::predict(base291);
    preds.push_back(V);
    
    // L (291 + LSFER features from V)
    auto lsfer_L = AbrahamV2Cascade::createLSFERFeatures(preds, 1);
    std::vector<double> feat_L = base291;
    feat_L.insert(feat_L.end(), lsfer_L.begin(), lsfer_L.end());
    double L = AbrahamV2Cascade::LModel::predict(feat_L);
    preds.push_back(L);
    
    // E (291 + LSFER features from V, L)
    auto lsfer_E = AbrahamV2Cascade::createLSFERFeatures(preds, 2);
    std::vector<double> feat_E = base291;
    feat_E.insert(feat_E.end(), lsfer_E.begin(), lsfer_E.end());
    double E = AbrahamV2Cascade::EModel::predict(feat_E);
    preds.push_back(E);
    
    // B (291 + LSFER features from V, L, E)
    auto lsfer_B = AbrahamV2Cascade::createLSFERFeatures(preds, 3);
    std::vector<double> feat_B = base291;
    feat_B.insert(feat_B.end(), lsfer_B.begin(), lsfer_B.end());
    double B = AbrahamV2Cascade::BModel::predict(feat_B);
    preds.push_back(B);
    
    // S (291 + LSFER features from V, L, E, B) - XGBoost
    auto lsfer_S = AbrahamV2Cascade::createLSFERFeatures(preds, 4);
    std::vector<double> feat_S = base291;
    feat_S.insert(feat_S.end(), lsfer_S.begin(), lsfer_S.end());
    double S = AbrahamSCascadeModel::predict(feat_S);
    preds.push_back(S);
    
    // A (291 + LSFER features from V, L, E, B, S) - XGBoost
    auto lsfer_A = AbrahamV2Cascade::createLSFERFeatures(preds, 5);
    std::vector<double> feat_A = base291;
    feat_A.insert(feat_A.end(), lsfer_A.begin(), lsfer_A.end());
    double A = AbrahamACascadeModel::predict(feat_A);
    
    // Return in standard order: [A, B, S, E, L, V]
    return {A, B, S, E, L, V};
}

} // namespace Osmordred
} // namespace Descriptors
} // namespace RDKit


