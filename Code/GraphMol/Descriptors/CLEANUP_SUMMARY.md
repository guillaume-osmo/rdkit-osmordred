# PhysChemProp Code Cleanup Summary

## Overview
Cleaned up and organized PhysChemProp.cpp dependencies by:
1. Moving obsolete headers to `OLD/` folder
2. Organizing cascade model headers in `physchemprop/` folder
3. Updating includes and CMakeLists.txt

## File Organization

### ✅ physchemprop/ folder (Active Cascade Models)
Contains all active cascade model headers and related files:
- **Cascade exported headers** (16 files): All `*_cascade_exported.h` files
- **Core headers**:
  - `OsmordredCascadeIntegration.h`
  - `OsmordredCascadePolarizabilityV2.h`
  - `RDKit217Descriptors.h`
  - `logodt_top1000_features.h`
- **Source files**:
  - `RDKit217Descriptors.cpp`
  - `logodt_top1000_features.cpp`

### ✅ OLD/ folder (Obsolete Headers)
Contains 17 obsolete header files that are no longer used by PhysChemProp.cpp:
- `OsmordredMeta37Simple.h` (restored to main dir - still used by other files)
- `OsmordredCascadePolarizability.h` (replaced by V2)
- `OsmordredMeta37Logodt.h`
- `OsmordredCascadeMeta37BP.h` (replaced by bp_cascade_exported.h)
- `OsmordredCascadeMeta37LogVP.h` (replaced by logvp_cascade_exported.h)
- `OsmordredCascadeMeta37DeltaHf.h` (replaced by deltahf_cascade_exported.h)
- `OsmordredCascadeMeta37DeltaHc.h` (replaced by deltahc_cascade_exported.h)
- All other `OsmordredCascadeMeta37*.h` headers (replaced by cascade exported models)

## Changes Made

### PhysChemProp.cpp
- Removed obsolete includes
- Updated includes to use `physchemprop/` prefix for cascade models
- Cleaned up commented-out code references

### CMakeLists.txt
- Added `include_directories(Descriptors/physchemprop)` to include path
- Updated source file paths for moved .cpp files

### Wrap/rdMolDescriptors.cpp
- Updated `RDKit217Descriptors.h` include path to use `physchemprop/` prefix

## Notes
- `OsmordredMeta37Simple.h` was restored to main directory as it's still used by:
  - `OsmordredMeta37Simple.cpp`
  - `OsmordredMeta37Integration.cpp`
  - `OsmordredParallel.cpp`
  - `Wrap/rdMolDescriptors.cpp`

## Next Steps
1. Test compilation with full build
2. Verify all cascade models work correctly
3. Consider removing `OsmordredMeta37Simple.h` dependencies if those files become obsolete
