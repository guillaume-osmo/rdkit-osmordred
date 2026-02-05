#!/usr/bin/env python
"""
Generate NCI 1000 Golden Reference for Osmordred v2.0

This script computes Osmordred descriptors for 1000 molecules from the NCI dataset
and stores them in a compressed format for comparison testing.

Usage:
    conda run -n osmo python generate_nci_1000_golden.py
"""

import gzip
import numpy as np
import os
import sys
import time

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

# Suppress RDKit warnings for cleaner output
RDLogger.DisableLog('rdApp.*')

def main():
    # Check for Osmordred support
    if not AllChem.HasOsmordredSupport():
        print("ERROR: RDKit was not built with Osmordred support!")
        print("Please use a build with RDK_BUILD_OSMORDRED_SUPPORT=ON")
        sys.exit(1)
    
    # Find NCI file
    rdbase = os.environ.get('RDBASE', '')
    
    # Try multiple locations
    possible_paths = [
        os.path.join(rdbase, 'Data', 'NCI', 'first_5K.smi') if rdbase else '',
        '/Users/guillaume-osmo/Github/rdkit-osmordred/Data/NCI/first_5K.smi',
        '/Users/guillaume-osmo/Github/rdkit/Data/NCI/first_5K.smi',
        'Data/NCI/first_5K.smi',
    ]
    
    nci_path = None
    for p in possible_paths:
        if p and os.path.exists(p):
            nci_path = p
            break
    
    if nci_path is None:
        print("ERROR: Could not find NCI first_5K.smi file")
        print("Searched locations:")
        for p in possible_paths:
            print(f"  - {p}")
        sys.exit(1)
    
    print(f"Using NCI file: {nci_path}")
    
    # Get descriptor names
    names = AllChem.GetOsmordredDescriptorNames()
    n_descriptors = len(names)
    print(f"Osmordred descriptor count: {n_descriptors}")
    
    # Process molecules
    smiles_list = []
    results_list = []
    max_mols = 1000
    
    print(f"\nProcessing up to {max_mols} molecules...")
    start_time = time.time()
    
    with open(nci_path, 'r') as f:
        for i, line in enumerate(f):
            if len(smiles_list) >= max_mols:
                break
            
            # Parse SMILES (tab-separated: smiles, name)
            parts = line.strip().split('\t')
            if not parts:
                continue
            smiles = parts[0].strip()
            if not smiles:
                continue
            
            # Parse molecule
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            
            # Skip very large molecules (>200 atoms or >10 rings)
            if mol.GetNumHeavyAtoms() > 200:
                continue
            try:
                ri = mol.GetRingInfo()
                if ri.NumRings() > 10:
                    continue
            except:
                pass
            
            # Calculate descriptors with timeout
            try:
                result = AllChem.CalcOsmordredWithTimeout(mol, 60.0)  # 60 second timeout
            except Exception as e:
                print(f"  Molecule {i}: Error - {e}")
                continue
            
            if len(result) != n_descriptors:
                print(f"  Molecule {i}: Wrong descriptor count ({len(result)} != {n_descriptors})")
                continue
            
            smiles_list.append(smiles)
            results_list.append(np.array(result, dtype=np.float64))
            
            if len(smiles_list) % 100 == 0:
                elapsed = time.time() - start_time
                print(f"  Processed {len(smiles_list)} molecules in {elapsed:.1f}s")
    
    elapsed = time.time() - start_time
    print(f"\nCompleted: {len(smiles_list)} molecules in {elapsed:.1f}s")
    print(f"Average: {elapsed / len(smiles_list):.3f}s per molecule")
    
    # Convert to numpy arrays
    smiles_array = np.array(smiles_list, dtype=object)
    results_array = np.array(results_list, dtype=np.float64)
    names_array = np.array(names, dtype=object)
    
    print(f"\nResults shape: {results_array.shape}")
    print(f"Non-NaN values: {np.sum(~np.isnan(results_array))}")
    print(f"NaN values: {np.sum(np.isnan(results_array))}")
    print(f"Inf values: {np.sum(np.isinf(results_array))}")
    
    # Save as compressed NPZ
    output_path = 'Code/GraphMol/Descriptors/test_data/nci_1000_osmordred_golden.npz'
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    np.savez_compressed(
        output_path,
        smiles=smiles_array,
        descriptors=results_array,
        descriptor_names=names_array,
        # Metadata
        n_molecules=len(smiles_list),
        n_descriptors=n_descriptors,
        version='osmordred-v2.0',
        tolerance=1e-5
    )
    
    file_size = os.path.getsize(output_path)
    print(f"\nSaved golden reference to: {output_path}")
    print(f"File size: {file_size / (1024*1024):.2f} MB")
    
    # Also save as gzipped CSV for the C++ test
    csv_path = 'Code/GraphMol/Descriptors/test_data/nci_1000_osmordred_golden.csv.gz'
    print(f"\nSaving CSV version to: {csv_path}")
    
    with gzip.open(csv_path, 'wt') as f:
        # Header
        f.write('smiles,' + ','.join(names) + '\n')
        # Data with full precision
        for smiles, result in zip(smiles_list, results_list):
            values_str = ','.join(f'{v:.15g}' for v in result)
            f.write(f'{smiles},{values_str}\n')
    
    csv_size = os.path.getsize(csv_path)
    print(f"CSV file size: {csv_size / (1024*1024):.2f} MB")
    
    # Verify the saved data
    print("\n--- Verification ---")
    loaded = np.load(output_path, allow_pickle=True)
    print(f"Loaded {loaded['n_molecules']} molecules")
    print(f"Loaded {loaded['n_descriptors']} descriptors")
    print(f"Version: {loaded['version']}")
    print(f"Tolerance: {loaded['tolerance']}")
    
    # Sample verification
    print("\nSample molecule 0:")
    print(f"  SMILES: {loaded['smiles'][0]}")
    print(f"  First 10 descriptors: {loaded['descriptors'][0][:10]}")
    
    print("\nGolden reference generation complete!")
    print(f"Files created:")
    print(f"  - {output_path} (NPZ format)")
    print(f"  - {csv_path} (CSV format)")

if __name__ == '__main__':
    main()
