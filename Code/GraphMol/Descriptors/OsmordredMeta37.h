// Meta37 Property Predictions - Osmordred Integration
// Auto-generated from production models trained on full data

#ifndef OSMORDRED_META37_H
#define OSMORDRED_META37_H

#include <vector>
#include <GraphMol/ROMol.h>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// Calculate all Meta37 properties at once
// Returns: [dD, dP, dH, Density, RI, BP, logPow, logWS, logVP, deltaHf, deltaHc, Polarizability]
std::vector<double> calcMeta37(const ROMol& mol);

// Individual property functions
double calcHansenDD(const ROMol& mol);
double calcHansenDP(const ROMol& mol);
double calcHansenDH(const ROMol& mol);
double calcDensity(const ROMol& mol);
double calcRefractiveIndex(const ROMol& mol);
double calcBoilingPoint(const ROMol& mol);
double calcLogPow(const ROMol& mol);
double calcLogWS(const ROMol& mol);
double calcLogVP(const ROMol& mol);
double calcDeltaHf(const ROMol& mol);
double calcDeltaHc(const ROMol& mol);
double calcMeta37Polarizability(const ROMol& mol);

}  // namespace Osmordred
}  // namespace Descriptors
}  // namespace RDKit

#endif // OSMORDRED_META37_H
