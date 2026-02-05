// Auto-generated: XGBoost Cascade Meta37 Integration
// 13 production models with 2nd order golden features
#pragma once

#include "OsmordredXGBdD.h"
#include "OsmordredXGBdP.h"
#include "OsmordredXGBdH.h"
#include "OsmordredXGBDensity.h"
#include "OsmordredXGBRI.h"
#include "OsmordredXGBBP.h"
#include "OsmordredXGBlogPow.h"
#include "OsmordredXGBlogWS.h"
#include "OsmordredXGBlogVP.h"
#include "OsmordredXGBA_abraham.h"
#include "OsmordredXGBS_abraham.h"
#include "OsmordredXGBMP.h"
#include "OsmordredXGBlogHenrycc.h"

#include <vector>
#include <cmath>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// Create LSFER 2nd order features from Abraham parameters
inline std::vector<double> createLSFER(const std::vector<double>& abr) {
    // abr = [A, B, S, E, L, V]
    std::vector<double> lsfer;
    // Raw (6)
    for (double v : abr) lsfer.push_back(v);
    // Squares (6)
    for (double v : abr) lsfer.push_back(v * v);
    // Products (15)
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = i + 1; j < 6; ++j) {
            lsfer.push_back(abr[i] * abr[j]);
        }
    }
    return lsfer;  // 27 features
}

// Cascade prediction for all Meta37 properties
inline std::vector<double> calcMeta37Cascade(
    const std::vector<double>& osmordred,
    const std::vector<double>& abrahamv2
) {
    std::vector<double> results(13, 0.0);
    std::vector<double> lsfer = createLSFER(abrahamv2);
    
    // 0: dD (no cascade)
    results[0] = dDXGB::predict(osmordred, lsfer, {});
    
    // 1: dP (uses dD)
    results[1] = dPXGB::predict(osmordred, lsfer, {results[0]});
    
    // 2: dH (no cascade)
    results[2] = dHXGB::predict(osmordred, lsfer, {});
    
    // 3: Density (no cascade)
    results[3] = DensityXGB::predict(osmordred, lsfer, {});
    
    // 4: RI (uses Density, dD)
    results[4] = RIXGB::predict(osmordred, lsfer, {results[3], results[0]});
    
    // 5: BP (no cascade)
    results[5] = BPXGB::predict(osmordred, lsfer, {});
    
    // 6: logPow (no cascade)
    results[6] = logPowXGB::predict(osmordred, lsfer, {});
    
    // 7: logWS (no cascade)
    results[7] = logWSXGB::predict(osmordred, lsfer, {});
    
    // 8: logVP (no cascade)
    results[8] = logVPXGB::predict(osmordred, lsfer, {});
    
    // 9: A_abraham (no cascade)
    results[9] = A_abrahamXGB::predict(osmordred, lsfer, {});
    
    // 10: S_abraham (no cascade)
    results[10] = S_abrahamXGB::predict(osmordred, lsfer, {});
    
    // 11: MP (no cascade)
    results[11] = MPXGB::predict(osmordred, lsfer, {});
    
    // 12: logHenrycc (no cascade)
    results[12] = logHenryccXGB::predict(osmordred, lsfer, {});
    
    return results;
}

// Property indices
enum Meta37CascadeIndex {
    META37_DD = 0,
    META37_DP = 1,
    META37_DH = 2,
    META37_DENSITY = 3,
    META37_RI = 4,
    META37_BP = 5,
    META37_LOGPOW = 6,
    META37_LOGWS = 7,
    META37_LOGVP = 8,
    META37_A_ABRAHAM = 9,
    META37_S_ABRAHAM = 10,
    META37_MP = 11,
    META37_LOGHENRYCC = 12
};

}  // namespace Osmordred
}  // namespace Descriptors
}  // namespace RDKit
