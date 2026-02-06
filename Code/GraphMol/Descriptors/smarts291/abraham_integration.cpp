// CLEAN Abraham Integration - Uses new query system
// Generated to match trained model EXACTLY
// EACH MODEL TYPE USES ITS OWN GOLDEN FEATURES!
//
// Note: The full calcAbrahams() function requires trained model headers
// (AbrahamGBAABRAHAM.h, etc.) which are generated separately.
// Define HAVE_ABRAHAM_MODELS to enable full Abraham parameter prediction.

#include "SMARTS291.h"
#include "abraham_queries.h"
#include <GraphMol/RWMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <vector>
#include <string>
#include <stdexcept>

#ifdef HAVE_ABRAHAM_MODELS
#include "AbrahamGBAABRAHAM.h"
#include "AbrahamGBSABRAHAM.h"
#include "AbrahamRidgeBABRAHAM.h"
#include "AbrahamRidgeEABRAHAM.h"
#include "AbrahamRidgeLABRAHAM.h"
#include "AbrahamRidgeVABRAHAM.h"
#endif

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// Extract 241 base features using the new query system
std::vector<double> extractAbrahamBaseFeatures(const RDKit::ROMol& mol) {
    std::vector<double> features;
    features.reserve(241);
    
    // Use the new query system that returns queries in EXACT model order
    auto queries = GetQueriesAbrahamBaseFeatures();
    
    // Safety check: if queries.size() is wrong, limit to 241
    size_t max_queries = (queries.size() > 241) ? 241 : queries.size();
    
    // Make a mutable RWMol copy for SubstructMatch (requires non-const reference)
    RDKit::RWMol mol_rw(mol);
    
    for (size_t i = 0; i < max_queries; ++i) {
        if (queries[i]) {
            std::vector<RDKit::MatchVectType> matches;
            RDKit::SubstructMatch(mol_rw, *queries[i], matches, true);
            features.push_back(static_cast<double>(matches.size()));
        } else {
            features.push_back(0.0);
        }
    }
    
    // Pad to 241 if we got fewer features
    while (features.size() < 241) {
        features.push_back(0.0);
    }

    return features;
}

std::vector<double> generateGoldenFeaturesA(const std::vector<double>& baseFeatures) {
    std::vector<double> goldenFeatures;
    goldenFeatures.reserve(50);

    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[89] / baseFeatures[113] : 0.0);  // abraham_BSEL_frag_3/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[29] / baseFeatures[127] : 0.0);  // abraham_A_frag_21/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[91] / baseFeatures[113] : 0.0);  // abraham_BSEL_frag_31/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[108] / baseFeatures[127] : 0.0);  // abraham_BSEL_frag_47/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[29] / baseFeatures[113] : 0.0);  // abraham_A_frag_21/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[39] / baseFeatures[113] : 0.0);  // abraham_A_frag_30/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[54] / baseFeatures[113] : 0.0);  // abraham_A_frag_44/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[0] / baseFeatures[127] : 0.0);  // abraham_A_additional_0/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[54] / baseFeatures[127] : 0.0);  // abraham_A_frag_44/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[34] / baseFeatures[113] : 0.0);  // abraham_A_frag_26/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[18] / baseFeatures[113] : 0.0);  // abraham_A_frag_11/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[90] / baseFeatures[127] : 0.0);  // abraham_BSEL_frag_30/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[90] / baseFeatures[113] : 0.0);  // abraham_BSEL_frag_30/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[18] / baseFeatures[127] : 0.0);  // abraham_A_frag_11/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[89] / baseFeatures[127] : 0.0);  // abraham_BSEL_frag_3/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[55] / baseFeatures[113] : 0.0);  // abraham_A_frag_45/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[108] / baseFeatures[113] : 0.0);  // abraham_BSEL_frag_47/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[2] / baseFeatures[127] : 0.0);  // abraham_A_additional_10/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[39] / baseFeatures[127] : 0.0);  // abraham_A_frag_30/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[91] / baseFeatures[127] : 0.0);  // abraham_BSEL_frag_31/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[51] / baseFeatures[113] : 0.0);  // abraham_A_frag_41/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[9] / baseFeatures[127] : 0.0);  // abraham_A_additional_4/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[107] / baseFeatures[127] : 0.0);  // abraham_BSEL_frag_46/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[51] / baseFeatures[127] : 0.0);  // abraham_A_frag_41/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[13] / baseFeatures[113] : 0.0);  // abraham_A_additional_8/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[93] / baseFeatures[127] : 0.0);  // abraham_BSEL_frag_33/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[9] / baseFeatures[113] : 0.0);  // abraham_A_additional_4/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[55] / baseFeatures[127] : 0.0);  // abraham_A_frag_45/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[93] / baseFeatures[113] : 0.0);  // abraham_BSEL_frag_33/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[11] / baseFeatures[127] : 0.0);  // abraham_A_additional_6/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[52] / baseFeatures[127] : 0.0);  // abraham_A_frag_42/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[92] / baseFeatures[127] : 0.0);  // abraham_BSEL_frag_32/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[92] / baseFeatures[113] : 0.0);  // abraham_BSEL_frag_32/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[53] / baseFeatures[113] : 0.0);  // abraham_A_frag_43/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[2] / baseFeatures[113] : 0.0);  // abraham_A_additional_10/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[53] / baseFeatures[127] : 0.0);  // abraham_A_frag_43/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[52] / baseFeatures[113] : 0.0);  // abraham_A_frag_42/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[56] / baseFeatures[113] : 0.0);  // abraham_A_frag_46/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[34] / baseFeatures[127] : 0.0);  // abraham_A_frag_26/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[56] / baseFeatures[127] : 0.0);  // abraham_A_frag_46/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[59] / baseFeatures[127] : 0.0);  // abraham_A_frag_49/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[85] / baseFeatures[127] : 0.0);  // abraham_BSEL_frag_26/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[85] / baseFeatures[113] : 0.0);  // abraham_BSEL_frag_26/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[38] / baseFeatures[113] : 0.0);  // abraham_A_frag_3/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[110] / baseFeatures[113] : 0.0);  // abraham_BSEL_frag_49/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[110] / baseFeatures[127] : 0.0);  // abraham_BSEL_frag_49/abraham_BSEL_frag_64
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[60] / baseFeatures[113] : 0.0);  // abraham_A_frag_5/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[59] / baseFeatures[113] : 0.0);  // abraham_A_frag_49/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[113] != 0.0 ? baseFeatures[3] / baseFeatures[113] : 0.0);  // abraham_A_additional_11/abraham_BSEL_frag_51
    goldenFeatures.push_back(baseFeatures[127] != 0.0 ? baseFeatures[84] / baseFeatures[127] : 0.0);  // abraham_BSEL_frag_25/abraham_BSEL_frag_64

    return goldenFeatures;
}

std::vector<double> generateGoldenFeaturesS(const std::vector<double>& baseFeatures) {
    std::vector<double> goldenFeatures;
    goldenFeatures.reserve(50);

    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[82] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_23/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[15] / baseFeatures[186] : 0.0);  // abraham_A_frag_0/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[66] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_0/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[89] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_3/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[175] / baseFeatures[186] : 0.0);  // abraham_new_[CX3H1]([CX4H3])=[CX3H2]/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[170] / baseFeatures[186] : 0.0);  // abraham_new_[CX3H0]([CX4H3])=[CX3H2]/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[12] / baseFeatures[186] : 0.0);  // abraham_A_additional_7/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[63] / baseFeatures[186] : 0.0);  // abraham_A_frag_7/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[4] / baseFeatures[186] : 0.0);  // abraham_A_additional_12/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[140] / baseFeatures[186] : 0.0);  // abraham_S_additional_0/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[8] / baseFeatures[186] : 0.0);  // abraham_A_additional_3/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[113] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_51/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[47] / baseFeatures[186] : 0.0);  // abraham_A_frag_38/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[60] / baseFeatures[186] : 0.0);  // abraham_A_frag_5/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[84] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_25/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[83] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_24/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[23] / baseFeatures[186] : 0.0);  // abraham_A_frag_16/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[90] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_30/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[151] / baseFeatures[186] : 0.0);  // abraham_S_additional_6/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[171] / baseFeatures[186] : 0.0);  // abraham_new_[CX3H0]=[CX3H0]/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[181] / baseFeatures[186] : 0.0);  // abraham_new_[CX3](=O)[OX2]/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[133] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_7/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[99] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_39/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[80] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_21/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[118] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_56/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[94] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_34/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[86] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_27/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[78] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_2/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[56] / baseFeatures[186] : 0.0);  // abraham_A_frag_46/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[182] / baseFeatures[186] : 0.0);  // abraham_new_[CX3]=[CX3]/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[108] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_47/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[124] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_61/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[169] / baseFeatures[186] : 0.0);  // abraham_new_[CX3H0]([CX4H3])=[CX3H1]/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[141] / baseFeatures[186] : 0.0);  // abraham_S_additional_1/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[98] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_38/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[64] / baseFeatures[186] : 0.0);  // abraham_A_frag_8/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[6] / baseFeatures[186] : 0.0);  // abraham_A_additional_14/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[150] / baseFeatures[186] : 0.0);  // abraham_S_additional_5/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[39] / baseFeatures[186] : 0.0);  // abraham_A_frag_30/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[167] / baseFeatures[186] : 0.0);  // abraham_new_[CX3H0]([CX4H3])([CX4H3])=[CX3H1]/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[184] / baseFeatures[186] : 0.0);  // abraham_new_[CX3]=[CX3][CX4]([CX4])([CX4])/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[72] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_14/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[58] / baseFeatures[186] : 0.0);  // abraham_A_frag_48/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[51] / baseFeatures[186] : 0.0);  // abraham_A_frag_41/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[41] / baseFeatures[186] : 0.0);  // abraham_A_frag_32/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[85] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_26/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[32] / baseFeatures[186] : 0.0);  // abraham_A_frag_24/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[92] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_32/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[14] / baseFeatures[186] : 0.0);  // abraham_A_additional_9/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[168] / baseFeatures[186] : 0.0);  // abraham_new_[CX3H0]([CX4H3])([CX4H3])=[CX3H2]/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]

    return goldenFeatures;
}

std::vector<double> generateGoldenFeaturesRidge(const std::vector<double>& baseFeatures) {
    std::vector<double> goldenFeatures;
    goldenFeatures.reserve(50);

    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[37] / baseFeatures[58] : 0.0);  // abraham_A_frag_29/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[82] / baseFeatures[186] : 0.0);  // abraham_BSEL_frag_23/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[186] != 0.0 ? baseFeatures[15] / baseFeatures[186] : 0.0);  // abraham_A_frag_0/abraham_new_[CX4H0]([CX4H2])([CX4H2])[CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[198] != 0.0 ? baseFeatures[147] / baseFeatures[198] : 0.0);  // abraham_S_additional_2/abraham_new_[CX4H1]([CX4H3])[CX4H2][CX4H2][CX4H3]
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[39] / baseFeatures[58] : 0.0);  // abraham_A_frag_30/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[50] / baseFeatures[58] : 0.0);  // abraham_A_frag_40/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[23] / baseFeatures[58] : 0.0);  // abraham_A_frag_16/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[9] / baseFeatures[58] : 0.0);  // abraham_A_additional_4/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[33] / baseFeatures[58] : 0.0);  // abraham_A_frag_25/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[46] / baseFeatures[58] : 0.0);  // abraham_A_frag_37/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[38] / baseFeatures[58] : 0.0);  // abraham_A_frag_3/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[48] / baseFeatures[58] : 0.0);  // abraham_A_frag_39/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[21] / baseFeatures[58] : 0.0);  // abraham_A_frag_14/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[29] / baseFeatures[58] : 0.0);  // abraham_A_frag_21/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[3] / baseFeatures[58] : 0.0);  // abraham_A_additional_11/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[34] / baseFeatures[58] : 0.0);  // abraham_A_frag_26/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[4] / baseFeatures[58] : 0.0);  // abraham_A_additional_12/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[18] / baseFeatures[58] : 0.0);  // abraham_A_frag_11/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[45] / baseFeatures[58] : 0.0);  // abraham_A_frag_36/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[54] / baseFeatures[58] : 0.0);  // abraham_A_frag_44/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[12] / baseFeatures[58] : 0.0);  // abraham_A_additional_7/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[35] / baseFeatures[58] : 0.0);  // abraham_A_frag_27/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[32] / baseFeatures[58] : 0.0);  // abraham_A_frag_24/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[40] / baseFeatures[58] : 0.0);  // abraham_A_frag_31/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[26] / baseFeatures[58] : 0.0);  // abraham_A_frag_19/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[43] / baseFeatures[58] : 0.0);  // abraham_A_frag_34/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[47] / baseFeatures[58] : 0.0);  // abraham_A_frag_38/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[56] / baseFeatures[58] : 0.0);  // abraham_A_frag_46/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[6] / baseFeatures[58] : 0.0);  // abraham_A_additional_14/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[27] / baseFeatures[58] : 0.0);  // abraham_A_frag_2/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[22] / baseFeatures[58] : 0.0);  // abraham_A_frag_15/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[13] / baseFeatures[58] : 0.0);  // abraham_A_additional_8/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[11] / baseFeatures[58] : 0.0);  // abraham_A_additional_6/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[52] / baseFeatures[58] : 0.0);  // abraham_A_frag_42/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[51] / baseFeatures[58] : 0.0);  // abraham_A_frag_41/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[8] / baseFeatures[58] : 0.0);  // abraham_A_additional_3/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[5] / baseFeatures[58] : 0.0);  // abraham_A_additional_13/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[55] / baseFeatures[58] : 0.0);  // abraham_A_frag_45/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[28] / baseFeatures[58] : 0.0);  // abraham_A_frag_20/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[2] / baseFeatures[58] : 0.0);  // abraham_A_additional_10/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[14] / baseFeatures[58] : 0.0);  // abraham_A_additional_9/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[0] / baseFeatures[58] : 0.0);  // abraham_A_additional_0/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[41] / baseFeatures[58] : 0.0);  // abraham_A_frag_32/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[19] / baseFeatures[58] : 0.0);  // abraham_A_frag_12/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[44] / baseFeatures[58] : 0.0);  // abraham_A_frag_35/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[17] / baseFeatures[58] : 0.0);  // abraham_A_frag_10/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[31] / baseFeatures[58] : 0.0);  // abraham_A_frag_23/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[53] / baseFeatures[58] : 0.0);  // abraham_A_frag_43/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[30] / baseFeatures[58] : 0.0);  // abraham_A_frag_22/abraham_A_frag_48
    goldenFeatures.push_back(baseFeatures[58] != 0.0 ? baseFeatures[7] / baseFeatures[58] : 0.0);  // abraham_A_additional_2/abraham_A_frag_48

    return goldenFeatures;
}


#ifdef HAVE_ABRAHAM_MODELS
std::vector<double> calcAbrahams(const RDKit::ROMol& mol) {
    // Extract 241 base features
    std::vector<double> baseFeatures = extractAbrahamBaseFeatures(mol);
    
    // Handle unexpected feature counts gracefully (e.g., 482 = 2x241, likely threading bug)
    if (baseFeatures.size() != 241) {
        if (baseFeatures.size() == 482) {
            // Likely duplication bug - take first 241 features
            baseFeatures.resize(241);
        } else if (baseFeatures.size() > 241) {
            // Take first 241 features
            baseFeatures.resize(241);
        } else {
            // Pad with zeros if too few features
            baseFeatures.resize(241, 0.0);
        }
    }
    
    // Generate golden features for each model type
    std::vector<double> goldenA = generateGoldenFeaturesA(baseFeatures);
    std::vector<double> goldenS = generateGoldenFeaturesS(baseFeatures);
    std::vector<double> goldenRidge = generateGoldenFeaturesRidge(baseFeatures);
    
    // Combine base + golden for each model
    std::vector<double> featuresA = baseFeatures;
    featuresA.insert(featuresA.end(), goldenA.begin(), goldenA.end());
    
    std::vector<double> featuresS = baseFeatures;
    featuresS.insert(featuresS.end(), goldenS.begin(), goldenS.end());
    
    std::vector<double> featuresRidge = baseFeatures;
    featuresRidge.insert(featuresRidge.end(), goldenRidge.begin(), goldenRidge.end());
    
    // Predict using each model with its correct features
    double A = GBModel_AABRAHAM::predict(featuresA);
    double B = RidgeModel_BABRAHAM::predict_unscaled(featuresRidge);
    double S = GBModel_SABRAHAM::predict(featuresS);
    double E = RidgeModel_EABRAHAM::predict_unscaled(featuresRidge);
    double L = RidgeModel_LABRAHAM::predict_unscaled(featuresRidge);
    double V = RidgeModel_VABRAHAM::predict_unscaled(featuresRidge);
    
    return {A, B, S, E, L, V};
}
#endif // HAVE_ABRAHAM_MODELS

std::vector<double> calcAbrahamsFeatures(const RDKit::ROMol& mol) {
    // Returns the 291 features for the A model (base + A's golden features)
    std::vector<double> baseFeatures = extractAbrahamBaseFeatures(mol);
    std::vector<double> goldenA = generateGoldenFeaturesA(baseFeatures);
    
    std::vector<double> allFeatures = baseFeatures;
    allFeatures.insert(allFeatures.end(), goldenA.begin(), goldenA.end());
    
    return allFeatures;
}

} // namespace Osmordred

namespace SMARTS291 {

bool hasSMARTS291Support() {
    return true;  // SMARTS291 support is always available when this code is compiled
}

} // namespace SMARTS291

} // namespace Descriptors
} // namespace RDKit
