# CalcPhysChemProp RDKit Installation Guide

This guide explains how to build and install RDKit with custom features for CalcPhysChemProp:
- **Osmordred**: 3585 molecular descriptors
- **RDKit217**: 217 RDKit descriptors (C++ implementation)
- **SMARTS291**: 291 Abraham features for property prediction

All three support **batch processing with parallel threading**.

## Prerequisites

### 1. Create Conda environment
```bash
conda create -n rdkit_osmordred_build python=3.11
conda activate rdkit_osmordred_build
conda install -c conda-forge boost libboost-python libboost-devel eigen cmake numpy pandas
```

### 2. Install LAPACK (for Osmordred)
```bash
# macOS
brew install lapack

# Linux
sudo apt-get install liblapack-dev liblapacke-dev
```

## Build Instructions

### 1. Clone the repository
```bash
git clone https://github.com/guillaume-osmo/rdkit-osmordred.git
cd rdkit-osmordred
git checkout master  # master has all features merged
```

### 2. Configure with CMake
```bash
mkdir -p build && cd build

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DRDK_BUILD_PYTHON_WRAPPERS=ON \
  -DRDK_BUILD_OSMORDRED_SUPPORT=ON \
  -DRDK_BUILD_DESCRIPTORS3D=ON \
  -DRDK_INSTALL_INTREE=ON \
  -DBoost_ROOT=$CONDA_PREFIX \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"
```

### 3. Build
```bash
make -j8
```

### 4. Install to conda environment
```bash
conda activate rdkit_osmordred_build

# Copy rdkit Python package to site-packages
SITE_PACKAGES=$(python -c "import site; print(site.getsitepackages()[0])")
cp -r rdkit $SITE_PACKAGES/

# Copy dylibs to conda lib (for rpath resolution)
cp lib/*.dylib $CONDA_PREFIX/lib/

echo "Installed to $SITE_PACKAGES"
```

## Verify Installation

**IMPORTANT**: Unset `DYLD_LIBRARY_PATH` if it points to another RDKit:
```bash
unset DYLD_LIBRARY_PATH
```

Then test:
```python
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

mol = Chem.MolFromSmiles('CCO')
smiles_list = ['CCO', 'CCC', 'c1ccccc1']

# SMARTS291 (291 Abraham features)
abraham = rdMolDescriptors.CalcAbrahamFeatures(mol)
print(f"CalcAbrahamFeatures: {len(abraham)} features")  # 291

# SMARTS291 Batch (parallel)
batch = rdMolDescriptors.CalcAbrahamFeaturesBatch(smiles_list, 4)
print(f"CalcAbrahamFeaturesBatch: {len(batch)} x {len(batch[0])} features")

# Osmordred (3585 features)
osmordred = rdMolDescriptors.CalcOsmordred(mol)
print(f"CalcOsmordred: {len(osmordred)} features")  # 3585

# Osmordred Batch (parallel)
batch = rdMolDescriptors.CalcOsmordredBatch(smiles_list, 4)
print(f"CalcOsmordredBatch: {len(batch)} x {len(batch[0])} features")

# RDKit217 (217 features)
rdkit217 = rdMolDescriptors.ExtractRDKitDescriptorsFromMolsBatch([mol], 1)
print(f"RDKit217: {len(rdkit217[0])} features")  # 217

# RDKit217 Batch from SMILES (parallel)
batch = rdMolDescriptors.ExtractRDKitDescriptorsBatch(smiles_list, 4)
print(f"ExtractRDKitDescriptorsBatch: {len(batch)} x {len(batch[0])} features")

print("\nALL FEATURES WORKING!")
```

## Available Functions

### Single Molecule
| Feature | Function | Output |
|---------|----------|--------|
| SMARTS291 | `CalcAbrahamFeatures(mol)` | 291 doubles |
| Osmordred | `CalcOsmordred(mol)` | 3585 doubles |
| Osmordred w/timeout | `CalcOsmordredWithTimeout(mol, 60)` | 3585 doubles (60s timeout) |

### Batch Processing (Parallel)
| Feature | Function | Output |
|---------|----------|--------|
| SMARTS291 | `CalcAbrahamFeaturesBatch(smiles_list, n_jobs)` | List of 291-d vectors |
| SMARTS291 | `CalcAbrahamFeaturesBatchFromMols(mols, n_jobs)` | List of 291-d vectors |
| Osmordred | `CalcOsmordredBatch(smiles_list, n_jobs)` | List of 3585-d vectors |
| Osmordred | `CalcOsmordredBatchFromMols(mols, n_jobs)` | List of 3585-d vectors |
| RDKit217 | `ExtractRDKitDescriptorsBatch(smiles_list, n_jobs)` | List of 217-d vectors |
| RDKit217 | `ExtractRDKitDescriptorsFromMolsBatch(mols, n_jobs)` | List of 217-d vectors |

### Utilities
| Function | Description |
|----------|-------------|
| `GetRDKit217DescriptorNames()` | 217 descriptor names |
| `GetOsmordredDescriptorNames()` | 3585 descriptor names |
| `HasOsmordredSupport()` | Check if Osmordred is available |

## Performance Notes

### Osmordred Timeouts
- Osmordred has a **60 second timeout** per molecule
- Very complex molecules (>10 rings or >200 heavy atoms) may hang
- Use batch functions which handle timeouts automatically

### Recommended Batch Sizes
- Process molecules in chunks of ~100 to avoid memory issues
- Use `n_jobs=0` for auto-detection of CPU cores

## Troubleshooting

### Import Error: Library not loaded
If you see `Library not loaded: @rpath/libRDKit...`:
```bash
# Make sure dylibs are in conda lib
cp /path/to/rdkit-osmordred/build/lib/*.dylib $CONDA_PREFIX/lib/
```

### Wrong RDKit version loaded
If features are missing, check `DYLD_LIBRARY_PATH`:
```bash
echo $DYLD_LIBRARY_PATH
unset DYLD_LIBRARY_PATH  # Clear if pointing to another RDKit
```

## Branch History

The `master` branch contains all features merged:
- Osmordred v2 (3585 features)
- RDKit217 (217 features)  
- SMARTS291 (291 Abraham features)
- All batch functions with parallel processing

## Running Tests

```bash
cd build
ctest -R "testOsmordred|testRDKit217|testSMARTS291" -V
```
