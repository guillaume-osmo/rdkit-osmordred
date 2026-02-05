// Auto-generated master include for all Osmordred GBT models
// with parallel tree inference support using RDKit threading
#pragma once

#include "OsmordredGBTDDParallel.h"
#include "OsmordredGBTDPParallel.h"
#include "OsmordredGBTDHParallel.h"
#include "OsmordredGBTDensityParallel.h"
#include "OsmordredGBTRIParallel.h"
#include "OsmordredGBTBPParallel.h"
#include "OsmordredGBTLogPowParallel.h"
#include "OsmordredGBTLogWSParallel.h"
#include "OsmordredGBTLogVPParallel.h"
#include "OsmordredGBTDeltaHfParallel.h"
#include "OsmordredGBTDeltaHcParallel.h"
#include "OsmordredGBTPolarizabilityParallel.h"
#include "OsmordredGBTLogHenryccParallel.h"
#include "OsmordredGBTMPParallel.h"
#include "OsmordredGBTDipoleMomentParallel.h"
#include "OsmordredGBTFlashpointParallel.h"
#include "OsmordredGBTLogODTParallel.h"

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// Parallel Meta37 calculation structure
struct Meta37ParallelResults {
    double dD, dP, dH;
    double density, ri, bp;
    double logPow, logWS, logVP;
    double deltaHf, deltaHc;
    double polarizability, logHenrycc, mp;
    double dipoleMoment, flashpoint, logODT;
};

} // namespace Osmordred
} // namespace Descriptors
} // namespace RDKit
