# Double Arcsinh Bug Fix Plan

## Bug Summary

**Location:** `PhysChemProp.cpp` (both original `calcPhysChemProp` and optimized `calcPhysChemPropOptimal`)

**Affected Properties:**
1. **logWS** - was copying LSFER/physics from logPow (already arcsinh'd)
2. **Flashpoint** - was copying osmordred/LSFER/physics from MP (already arcsinh'd)
3. **logViscosity** - depends on Flashpoint and logWS (indirectly affected)

**Root Cause:**
```cpp
// BUG: Copying from previous model's processed features
for (size_t i = 0; i < 54; ++i) {
    features_logws_full[lsfer_start + i] = features_logpow_full[lsfer_logpow_start + i];
    // ^^^ features_logpow_full already had arcsinh applied!
}
// Then applying arcsinh again:
for (int i = 0; i < N_ARCSINH_COLUMNS; ++i) {
    features_logws_full[col] = std::asinh(features_logws_full[col]);
    // ^^^ DOUBLE arcsinh on LSFER/physics columns!
}
```

## Training Code Analysis

**Good News:** The Python training scripts build features **correctly**:
```python
# Training code (71_train_logws_cascade.py)
osmordred = np.array(rdMolDescriptors.CalcOsmordred(mol))  # Fresh
lsfer_54 = create_lsfer_features(abraham_9)  # Fresh from cascade_9
physics_10 = create_physics_features(props, mw)  # Fresh
combined = np.concatenate([osmordred, abraham_9, cascade_6, lsfer_54, physics_10])
# Then single arcsinh applied once
```

**The models were trained correctly!** The bug is C++ inference only.

## Fix Status

### ✅ Already Fixed in `calcPhysChemPropOptimal`:
- logWS: Uses fresh `lsfer_54` and `physics_10` from cached cascade_9
- Flashpoint: Uses fresh `cache.osmordred`, `lsfer_54`, `physics_10`

### ⚠️ Needs Fix in Original `calcPhysChemProp`:
- logWS: Fixed (builds fresh LSFER/physics)
- Flashpoint: Fixed (builds fresh osmordred/LSFER/physics)

### Remaining Discrepancy
After fixes, `CalcPhysChemProp` and `CalcPhysChemPropOptimal` show small differences:
- logWS: ~0.03 difference
- Flashpoint: ~0.5-4°C difference

**Root Cause of Remaining Difference:**
The cascade values (logPow, logWS, etc.) used in downstream models come from earlier predictions. 
If earlier predictions differ slightly, this propagates through the cascade.

## Action Plan

### Phase 1: Verify Fix (No Rebuild Needed)
- [x] Fix `calcPhysChemPropOptimal` to use fresh features (DONE)
- [x] Fix original `calcPhysChemProp` to use fresh features (DONE)
- [ ] Rebuild and test

### Phase 2: Model Verification
After rebuild, verify that:
1. `CalcPhysChemPropOptimal` produces better predictions (closer to training labels)
2. Both functions produce identical results (or document why they differ)

### Phase 3: Add Timeout for Large Molecules
Add timeout protection for molecules with:
- nHeavyAtoms > 200
- nRings > 12
- Processing time > 10 seconds

## Timeout Implementation Plan

### Option A: Early Exit (Recommended)
```cpp
// At start of calcPhysChemProp / calcPhysChemPropOptimal:
unsigned int nHeavy = mol.getNumHeavyAtoms();
unsigned int nRings = RDKit::Descriptors::calcNumRings(mol);
if (nHeavy > 200 || nRings > 12) {
    // Return default/NaN values instead of hanging
    std::vector<double> results(NUM_PHYSCHEM_PROPS, std::nan(""));
    return results;
}
```

### Option B: Timeout with Thread (More Complex)
```cpp
#include <future>
#include <chrono>

std::vector<double> calcPhysChemPropWithTimeout(const ROMol& mol, int timeout_sec = 10) {
    auto future = std::async(std::launch::async, [&]() {
        return calcPhysChemProp(mol);
    });
    
    if (future.wait_for(std::chrono::seconds(timeout_sec)) == std::future_status::timeout) {
        // Return NaN values
        return std::vector<double>(NUM_PHYSCHEM_PROPS, std::nan(""));
    }
    return future.get();
}
```

### Recommended Approach
Use **Option A** (early exit) for simplicity and reliability. Timeout with threads can have issues with:
- Thread cleanup if computation is stuck in RDKit internals
- Memory leaks if thread is abandoned
- Platform-specific behavior

## Files to Modify

1. **PhysChemProp.cpp**
   - Add size checks at start of `calcPhysChemProp` and `calcPhysChemPropOptimal`
   - Already contains the double-arcsinh fixes

2. **Osmordred.h**
   - Add constants: `MAX_HEAVY_ATOMS = 200`, `MAX_RINGS = 12`

## Testing Checklist

- [ ] Test with small molecule (CCO) - should work normally
- [ ] Test with large molecule (>200 heavy atoms) - should return NaN quickly
- [ ] Test with complex ring system (>12 rings) - should return NaN quickly
- [ ] Compare CalcPhysChemProp vs CalcPhysChemPropOptimal - should match
- [ ] Verify predictions match Python training expectations

## Summary

| Issue | Status | Impact |
|-------|--------|--------|
| Double arcsinh in logWS | ✅ Fixed | Models get correct features |
| Double arcsinh in Flashpoint | ✅ Fixed | Models get correct features |
| Timeout for large molecules | ⏳ TODO | Prevents hanging |
| Small discrepancy between orig/opt | ⚠️ Investigate | Floating point? |
