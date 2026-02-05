// Simple Meta37 Integration Header
#ifndef OSMORDRED_META37_SIMPLE_H
#define OSMORDRED_META37_SIMPLE_H

#include <RDGeneral/export.h>
#include <GraphMol/ROMol.h>
#include <vector>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// Calculate Meta37 properties using simple Osmordred + MW features
// Returns vector of 13 values:
//   0: dD, 1: dH, 2: dP, 3: BP, 4: logVP, 5: logPow, 6: logWS,
//   7: deltaHf, 8: deltaHc, 9: MP, 10: flashpoint, 11: logHenrycc, 12: dipolemoment
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcMeta37Simple(const ROMol& mol);

}  // namespace Osmordred
}  // namespace Descriptors
}  // namespace RDKit

#endif  // OSMORDRED_META37_SIMPLE_H

