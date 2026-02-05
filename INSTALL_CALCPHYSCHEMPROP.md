# CalcPhysChemProp RDKit Installation Guide

This guide explains how to build and install RDKit with custom features for CalcPhysChemProp:
- **Osmordred**: 3585 molecular descriptors
- **RDKit217**: 217 RDKit descriptors (C++ implementation)
- **SMARTS291**: 291 Abraham features for property prediction

## Prerequisites

1. **Conda environment with Boost**:
```bash
conda create -n rdkit_build python=3.11
conda activate rdkit_build
conda install -c conda-forge boost libboost-python libboost-devel eigen cmake
```

2. **LAPACK** (for Osmordred):
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
git checkout calcphyschemprop-with-bindings
```

### 2. Configure with CMake
```bash
mkdir build && cd build

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DRDK_BUILD_PYTHON_WRAPPERS=ON \
  -DRDK_BUILD_OSMORDRED_SUPPORT=ON \
  -DRDK_BUILD_DESCRIPTORS3D=ON \
  -DRDK_INSTALL_INTREE=ON \
  -DBoost_ROOT=$CONDA_PREFIX \
  -DPYTHON_EXECUTABLE=$CONDA_PREFIX/bin/python \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX;/opt/homebrew/opt/lapack"
```

### 3. Build
```bash
cmake --build . --parallel 8
```

### 4. Install to conda environment
```bash
cmake --install .

# Copy to site-packages
SITE_PACKAGES=$(python -c "import site; print(site.getsitepackages()[0])")
cp -r ../rdkit $SITE_PACKAGES/
cp -r ../lib/*.dylib $SITE_PACKAGES/rdkit/  # macOS
# or: cp -r ../lib/*.so $SITE_PACKAGES/rdkit/  # Linux

# Fix rpaths (macOS only)
cd $SITE_PACKAGES/rdkit
for so in *.so Chem/*.so ML/*/*.so DataManip/*/*.so SimDivFilters/*.so; do
  [ -f "$so" ] && install_name_tool -add_rpath "@loader_path" "$so" 2>/dev/null
done
```

## Verify Installation

```python
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

mol = Chem.MolFromSmiles('CCO')

# SMARTS291 (291 Abraham features)
abraham = rdMolDescriptors.CalcAbrahamFeatures(mol)
print(f"CalcAbrahamFeatures: {len(abraham)} features")  # 291

# Osmordred (3585 features)
osmordred = rdMolDescriptors.CalcOsmordred(mol)
print(f"CalcOsmordred: {len(osmordred)} features")  # 3585

# RDKit217 (217 features)
rdkit217 = rdMolDescriptors.ExtractRDKitDescriptorsFromMolsBatch([mol], 1)
print(f"RDKit217: {len(rdkit217[0])} features")  # 217

# Get descriptor names
names = rdMolDescriptors.GetRDKit217DescriptorNames()
print(f"RDKit217 names: {len(names)} names")  # 217
```

## Branch Structure

```
calcphyschemprop-with-bindings (recommended)
├── osmordred-v2-release-2025.09.3 (3585 features)
├── rdkit217-release-2025.09.3 (217 features)
└── smarts291-release-2025.09.3 (291 features)
```

## Features

| Feature | Function | Output |
|---------|----------|--------|
| SMARTS291 | `CalcAbrahamFeatures(mol)` | 291 doubles |
| Osmordred | `CalcOsmordred(mol)` | 3585 doubles |
| RDKit217 | `ExtractRDKitDescriptorsFromMolsBatch([mols], n_jobs)` | 217 doubles per mol |
| RDKit217 Names | `GetRDKit217DescriptorNames()` | 217 strings |

## Golden Reference Tests

The branch includes C++ unit tests with golden reference data for validation:
- `testOsmordred`: Validates against NCI 1000 molecules
- `testRDKit217`: Validates against NCI 1000 molecules + name validation
- `testSMARTS291`: Validates against NCI 1000 molecules

Run tests:
```bash
cd build
ctest -R "testOsmordred|testRDKit217|testSMARTS291" -V
```
