"""
Test module for RDKit Osmordred functionality.

This module provides comprehensive tests for all Osmordred descriptor functions
available through the Python wrapper.
"""

import unittest
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD


def to_list(result):
    """Convert RDKit vector objects to Python lists."""
    if hasattr(result, '__iter__') and not isinstance(result, (list, tuple)):
        return list(result)
    return result


class TestOsmordred(unittest.TestCase):
    """Test class for Osmordred functionality."""

    def setUp(self):
        """Set up test molecules."""
        self.ethanol = Chem.MolFromSmiles("CCO")
        self.benzene = Chem.MolFromSmiles("c1ccccc1")
        self.complex_mol = Chem.MolFromSmiles(
            "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"
        )
        self.t_butyl_benzene = Chem.MolFromSmiles("CC(C)(C)c1ccc(C(C)(C)C)cc1")
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_basic_functionality(self):
        """Test basic Osmordred functions with simple molecules."""
        self.assertIsNotNone(self.ethanol)
        
        # Test atom counts (using version 2 as default)
        atom_counts = rdMD.CalcAtomCount(self.ethanol)
        atom_counts_list = to_list(atom_counts)
        self.assertGreater(len(atom_counts_list), 0)
        
        # Test bond counts
        bond_counts = rdMD.CalcBondCount(self.ethanol)
        bond_counts_list = to_list(bond_counts)
        self.assertGreater(len(bond_counts_list), 0)
        
        # Test molecular weight
        weight = rdMD.CalcWeight(self.ethanol)
        weight_list = to_list(weight)
        self.assertGreater(len(weight_list), 0)
        self.assertGreater(weight_list[0], 0.0)
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_aromatic_functionality(self):
        """Test aromatic-related functions."""
        self.assertIsNotNone(self.benzene)
        
        # Test aromatic counts
        aromatic_counts = rdMD.CalcAromatic(self.benzene)
        aromatic_counts_list = to_list(aromatic_counts)
        self.assertGreater(len(aromatic_counts_list), 0)
        
        # Test ring count
        ring_counts = rdMD.CalcRingCount(self.benzene)
        ring_counts_list = to_list(ring_counts)
        self.assertGreater(len(ring_counts_list), 0)
        
        # Test aromatic atom count
        aromatic_atoms = rdMD.CalcCountAromaticAtoms(self.benzene)
        self.assertIsInstance(aromatic_atoms, int)
        self.assertGreaterEqual(aromatic_atoms, 0)
        
        # Test aromatic bond count
        aromatic_bonds = rdMD.CalcCountAromaticBonds(self.benzene)
        self.assertIsInstance(aromatic_bonds, int)
        self.assertGreaterEqual(aromatic_bonds, 0)
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_matrix_descriptors(self):
        """Test matrix-based descriptors."""
        self.assertIsNotNone(self.ethanol)
        
        # Test distance matrix descriptors
        dist_descs = rdMD.CalcDistanceMatrix(self.ethanol)
        dist_descs_list = to_list(dist_descs)
        self.assertGreater(len(dist_descs_list), 0)
        
        dist_descs_eigen = rdMD.CalcDistanceMatrixEigen(self.ethanol)
        dist_descs_eigen_list = to_list(dist_descs_eigen)
        self.assertGreater(len(dist_descs_eigen_list), 0)
        
        # Test adjacency matrix descriptors
        adj_descs = rdMD.CalcAdjacencyMatrix(self.ethanol)
        adj_descs_list = to_list(adj_descs)
        self.assertGreater(len(adj_descs_list), 0)
        
        adj_descs_eigen = rdMD.CalcAdjacencyMatrixEigen(self.ethanol)
        adj_descs_eigen_list = to_list(adj_descs_eigen)
        self.assertGreater(len(adj_descs_eigen_list), 0)
        
        # Test detour matrix descriptors
        detour_descs = rdMD.CalcDetourMatrix(self.ethanol)
        detour_descs_list = to_list(detour_descs)
        self.assertGreater(len(detour_descs_list), 0)
        
        detour_descs_eigen = rdMD.CalcDetourMatrixEigen(self.ethanol)
        detour_descs_eigen_list = to_list(detour_descs_eigen)
        self.assertGreater(len(detour_descs_eigen_list), 0)
        
        # Test Barysz matrix descriptors
        barysz_descs = rdMD.CalcBaryszMatrix(self.ethanol)
        barysz_descs_list = to_list(barysz_descs)
        self.assertGreater(len(barysz_descs_list), 0)
        
        barysz_descs_eigen = rdMD.CalcBaryszMatrixEigen(self.ethanol)
        barysz_descs_eigen_list = to_list(barysz_descs_eigen)
        self.assertGreater(len(barysz_descs_eigen_list), 0)
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_physical_properties(self):
        """Test physical property calculations."""
        self.assertIsNotNone(self.ethanol)
        
        # Test molecular weight
        weight = rdMD.CalcWeight(self.ethanol)
        weight_list = to_list(weight)
        self.assertGreater(len(weight_list), 0)
        self.assertGreater(weight_list[0], 0.0)
        
        # Test van der Waals volume
        vdw_vol = rdMD.CalcVdwVolumeABC(self.ethanol)
        vdw_vol_list = to_list(vdw_vol)
        self.assertGreater(vdw_vol_list[0], 0.0)
        
        # Test McGowan volume
        mcgowan_vol = rdMD.CalcMcGowanVolume(self.ethanol)
        mcgowan_vol_list = to_list(mcgowan_vol)
        self.assertGreater(mcgowan_vol_list[0], 0.0)
        
        # Test polarizability
        polarizability = rdMD.CalcPolarizability(self.ethanol)
        polarizability_list = to_list(polarizability)
        self.assertGreater(polarizability_list[0], 0.0)
        
        # Test SLogP
        slogp = rdMD.CalcSLogP(self.ethanol)
        slogp_list = to_list(slogp)
        self.assertIsInstance(slogp_list, list)
        
        # Test LogS
        logs = rdMD.CalcLogS(self.ethanol)
        logs_list = to_list(logs)
        self.assertIsInstance(logs_list, list)
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_topological_descriptors(self):
        """Test topological descriptors."""
        self.assertIsNotNone(self.ethanol)
        
        # Test Wiener index
        wiener = rdMD.CalcWienerIndex(self.ethanol)
        wiener_list = to_list(wiener)
        self.assertGreater(len(wiener_list), 0)
        
        # Test Balaban J index
        balaban_j = rdMD.CalcBalabanJ(self.ethanol)
        balaban_j_list = to_list(balaban_j)
        self.assertGreater(len(balaban_j_list), 0)
        
        # Test Bertz CT index
        bertz_ct = rdMD.CalcBertzCT(self.ethanol)
        bertz_ct_list = to_list(bertz_ct)
        self.assertGreater(len(bertz_ct_list), 0)
        
        # Test Zagreb index
        zagreb = rdMD.CalcZagrebIndex(self.ethanol)
        zagreb_list = to_list(zagreb)
        self.assertGreater(len(zagreb_list), 0)
        
        # Test eccentric connectivity index
        eccentric = rdMD.CalcEccentricConnectivityIndex(self.ethanol)
        eccentric_list = to_list(eccentric)
        self.assertGreater(len(eccentric_list), 0)
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_chi_descriptors(self):
        """Test Chi descriptors."""
        self.assertIsNotNone(self.ethanol)
        
        # Test all Chi descriptors
        chi_descs = rdMD.CalcChi(self.ethanol)
        chi_descs_list = to_list(chi_descs)
        self.assertGreater(len(chi_descs_list), 0)
        
        # Test individual Chi descriptors
        chipath = rdMD.CalcChipath(self.ethanol)
        chipath_list = to_list(chipath)
        self.assertGreater(len(chipath_list), 0)
        
        chichain = rdMD.CalcChichain(self.ethanol)
        chichain_list = to_list(chichain)
        self.assertGreater(len(chichain_list), 0)
        
        chicluster = rdMD.CalcChicluster(self.ethanol)
        chicluster_list = to_list(chicluster)
        self.assertGreater(len(chicluster_list), 0)
        
        chipathcluster = rdMD.CalcChipathcluster(self.ethanol)
        chipathcluster_list = to_list(chipathcluster)
        self.assertGreater(len(chipathcluster_list), 0)
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_estate_descriptors(self):
        """Test EState descriptors."""
        self.assertIsNotNone(self.ethanol)
        
        # Test EState descriptors
        estate_descs = rdMD.CalcEState(self.ethanol, False)
        estate_descs_list = to_list(estate_descs)
        self.assertGreater(len(estate_descs_list), 0)
        
        # Test extended EState descriptors
        estate_descs_ext = rdMD.CalcEState(self.ethanol, True)
        estate_descs_ext_list = to_list(estate_descs_ext)
        self.assertGreater(len(estate_descs_ext_list), 0)
        
        # Test BEState descriptors
        bestate_descs = rdMD.CalcBEState(self.ethanol)
        bestate_descs_list = to_list(bestate_descs)
        self.assertGreater(len(bestate_descs_list), 0)
        
        # Test HEState descriptors
        hestate_descs = rdMD.CalcHEState(self.ethanol)
        hestate_descs_list = to_list(hestate_descs)
        self.assertGreater(len(hestate_descs_list), 0)
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_count_functions(self):
        """Test count-based functions."""
        self.assertIsNotNone(self.ethanol)
        
        # Test acidic group count
        acidic_count = rdMD.CalcAcidicGroupCount(self.ethanol)
        self.assertIsInstance(acidic_count, int)
        self.assertGreaterEqual(acidic_count, 0)
        
        # Test basic group count
        basic_count = rdMD.CalcBasicGroupCount(self.ethanol)
        self.assertIsInstance(basic_count, int)
        self.assertGreaterEqual(basic_count, 0)
        
        # Test aromatic atom count
        aromatic_atoms = rdMD.CalcCountAromaticAtoms(self.ethanol)
        self.assertIsInstance(aromatic_atoms, int)
        self.assertGreaterEqual(aromatic_atoms, 0)
        
        # Test aromatic bond count
        aromatic_bonds = rdMD.CalcCountAromaticBonds(self.ethanol)
        self.assertIsInstance(aromatic_bonds, int)
        self.assertGreaterEqual(aromatic_bonds, 0)
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_complex_molecule(self):
        """Test functions with complex molecules."""
        self.assertIsNotNone(self.complex_mol)
        
        # Test that all functions work with complex molecules
        # Functions that need version parameter
        functions_with_version = [
        ]
        
        # Functions that need no additional parameters
        functions_simple = [
            rdMD.CalcABCIndex,
            rdMD.CalcAcidBase,
            rdMD.CalcAromatic,
            rdMD.CalcBalabanJ,
            rdMD.CalcBertzCT,
            rdMD.CalcBondCount,
            rdMD.CalcVertexAdjacencyInformation,
            rdMD.CalcWeight,
            rdMD.CalcWienerIndex,
            rdMD.CalcVdwVolumeABC,
            rdMD.CalcTopoPSA,
            rdMD.CalcSLogP,
            rdMD.CalcHydrogenBond,
            rdMD.CalcLogS,
            rdMD.CalcLipinski,
            rdMD.CalcMcGowanVolume,
            rdMD.CalcPolarizability,
            rdMD.CalcRotatableBond,
            rdMD.CalcFragmentComplexity,
            rdMD.CalcConstitutional,
            rdMD.CalcTopologicalIndex,
        ]
        
        # Test functions with version parameter
        for func, version in functions_with_version:
            try:
                result = func(self.complex_mol, version)
                result_list = to_list(result)
                self.assertGreater(len(result_list), 0)
            except Exception as e:
                self.fail(f"Function {func.__name__} failed: {e}")
        
        # Test functions without additional parameters
        for func in functions_simple:
            try:
                result = func(self.complex_mol)
                result_list = to_list(result)
                self.assertGreater(len(result_list), 0)
            except Exception as e:
                self.fail(f"Function {func.__name__} failed: {e}")
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_edge_cases(self):
        """Test edge cases and error handling."""
        # Test with None molecule
        with self.assertRaises((TypeError, AttributeError)):
            rdMD.CalcAtomCount(None)
        
        # Test with empty SMILES
        empty_mol = Chem.MolFromSmiles("")
        if empty_mol is not None:
            atom_counts = rdMD.CalcAtomCount(empty_mol)
            atom_counts_list = to_list(atom_counts)
            self.assertIsInstance(atom_counts_list, list)
        
        # Test with single atom
        single_atom = Chem.MolFromSmiles("[He]")
        if single_atom is not None:
            atom_counts = rdMD.CalcAtomCount(single_atom)
            atom_counts_list = to_list(atom_counts)
            self.assertIsInstance(atom_counts_list, list)
            self.assertGreater(len(atom_counts_list), 0)
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_consistency(self):
        """Test that multiple calls return consistent results."""
        self.assertIsNotNone(self.ethanol)
        
        # Test atom counts consistency
        atom_counts1 = rdMD.CalcAtomCount(self.ethanol)
        atom_counts2 = rdMD.CalcAtomCount(self.ethanol)
        atom_counts1_list = to_list(atom_counts1)
        atom_counts2_list = to_list(atom_counts2)
        self.assertEqual(atom_counts1_list, atom_counts2_list)
        
        # Test bond counts consistency
        bond_counts1 = rdMD.CalcBondCount(self.ethanol)
        bond_counts2 = rdMD.CalcBondCount(self.ethanol)
        bond_counts1_list = to_list(bond_counts1)
        bond_counts2_list = to_list(bond_counts2)
        self.assertEqual(bond_counts1_list, bond_counts2_list)
        
        # Test molecular weight consistency
        weight1 = rdMD.CalcWeight(self.ethanol)
        weight2 = rdMD.CalcWeight(self.ethanol)
        weight1_list = to_list(weight1)
        weight2_list = to_list(weight2)
        self.assertEqual(weight1_list, weight2_list)
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_information_content(self):
        """Test information content descriptors."""
        self.assertIsNotNone(self.ethanol)
        
        # Test information content
        info_content = rdMD.CalcInformationContent(self.ethanol, 5)
        info_content_list = to_list(info_content)
        self.assertGreater(len(info_content_list), 0)
        
        # Test with different max radius
        info_content_r3 = rdMD.CalcInformationContent(self.ethanol, 3)
        info_content_r3_list = to_list(info_content_r3)
        self.assertGreater(len(info_content_r3_list), 0)
        
        info_content_r7 = rdMD.CalcInformationContent(self.ethanol, 7)
        info_content_r7_list = to_list(info_content_r7)
        self.assertGreater(len(info_content_r7_list), 0)
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_all_functions_exist(self):
        """Test that all expected functions are available."""
        expected_functions = [
            'CalcABCIndex', 'CalcAcidBase', 'CalcAromatic', 'CalcAtomCount',
            'CalcBalabanJ', 'CalcBertzCT', 'CalcBondCount', 'CalcVertexAdjacencyInformation',
            'CalcWeight', 'CalcWienerIndex', 'CalcVdwVolumeABC', 'CalcTopoPSA',
            'CalcSLogP', 'CalcHydrogenBond', 'CalcLogS', 'CalcLipinski',
            'CalcMcGowanVolume', 'CalcPolarizability', 'CalcRotatableBond',
            'CalcFragmentComplexity', 'CalcConstitutional', 'CalcTopologicalIndex',
            'CalcDetourMatrixEigen', 'CalcDetourMatrix', 'CalcDistanceMatrixEigen',
            'CalcDistanceMatrix', 'CalcAdjacencyMatrixEigen', 'CalcAdjacencyMatrix',
            'CalcCarbonTypes', 'CalcEccentricConnectivityIndex', 'CalcBaryszMatrix',
            'CalcBaryszMatrixEigen', 'CalcZagrebIndex', 'CalcMoeType',
            'CalcMolecularDistanceEdge', 'CalcEState', 'CalcWalkCount',
            'CalcTopologicalCharge', 'CalcChi', 'CalcPathCount', 'CalcKappaShapeIndex',
            'CalcRingCount', 'CalcMolecularId', 'CalcBCUT', 'CalcAutocorrelation',
            'CalcFramework', 'CalcExtendedTopochemicalAtom',
            'CalcChipath', 'CalcChichain', 'CalcChicluster', 'CalcChipathcluster',
            'CalcAcidicGroupCount', 'CalcBasicGroupCount', 'CalcCountAromaticAtoms',
            'CalcCountAromaticBonds', 'CalcBEState', 'CalcHEState',
            'CalcAlphaKappaShapeIndex', 'CalcAbrahams', 'CalcPol', 'CalcMR',
            'CalcFlexibility', 'CalcODT', 'CalcSchultz', 'CalcRNCGRPCG',
            'CalcAZV', 'CalcASV', 'CalcDSV', 'CalcAZS', 'CalcASZ', 'CalcDN2S',
            'CalcDN2I', 'CalcASI', 'CalcDSI', 'CalcASN', 'CalcDSN', 'CalcDN2N',
            'CalcANS', 'CalcANV', 'CalcAZN', 'CalcANZ', 'CalcANI', 'CalcDSZ',
            'CalcANN', 'CalcDN2Z', 'CalcANMat', 'CalcAZMat', 'CalcASMat',
            'CalcDSMat', 'CalcDN2Mat', 'CalcFrags', 'CalcAddFeatures',
            'CalcInformationContent',
            # SMARTS291 and RDKit217 functions
            'CalcAbrahamFeatures',
            'ExtractSMARTS291Batch',  # Multi-threaded batch SMARTS291 extraction
            'GetSMARTS291FeatureNames',  # Get 291 SMARTS feature names
            'ExtractRDKitDescriptorsFromMolsBatch',
            'GetRDKit217DescriptorNames'
        ]
        
        for func_name in expected_functions:
            self.assertTrue(hasattr(rdMD, func_name),
                          f"Function {func_name} not found in Osmordred module")


class TestSMARTS291(unittest.TestCase):
    """Test class for SMARTS291 (Abraham features) functionality."""

    def setUp(self):
        """Set up test molecules."""
        self.ethanol = Chem.MolFromSmiles("CCO")
        self.benzene = Chem.MolFromSmiles("c1ccccc1")
        self.acetic_acid = Chem.MolFromSmiles("CC(=O)O")
        self.amine = Chem.MolFromSmiles("CCN")
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_calc_abraham_features_exists(self):
        """Test that CalcAbrahamFeatures function exists."""
        self.assertTrue(hasattr(rdMD, 'CalcAbrahamFeatures'))
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_calc_abraham_features_basic(self):
        """Test CalcAbrahamFeatures returns 291 features."""
        features = rdMD.CalcAbrahamFeatures(self.ethanol)
        features_list = to_list(features)
        self.assertEqual(len(features_list), 291, 
                        f"Expected 291 features, got {len(features_list)}")
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_calc_abraham_features_consistency(self):
        """Test CalcAbrahamFeatures returns consistent results."""
        features1 = rdMD.CalcAbrahamFeatures(self.ethanol)
        features2 = rdMD.CalcAbrahamFeatures(self.ethanol)
        features1_list = to_list(features1)
        features2_list = to_list(features2)
        self.assertEqual(features1_list, features2_list)
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_calc_abraham_features_different_molecules(self):
        """Test CalcAbrahamFeatures for different molecules."""
        for mol in [self.ethanol, self.benzene, self.acetic_acid, self.amine]:
            features = rdMD.CalcAbrahamFeatures(mol)
            features_list = to_list(features)
            self.assertEqual(len(features_list), 291)
            # Check no NaN in base features (first 241)
            base_features = features_list[:241]
            nan_count = sum(1 for f in base_features if f != f)  # NaN check
            self.assertEqual(nan_count, 0, "Base features should not contain NaN")


class TestSMARTS291Batch(unittest.TestCase):
    """Test class for SMARTS291 batch extraction with multi-threading."""

    def setUp(self):
        """Set up test SMILES."""
        self.smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "CCN", "CC(C)(C)c1ccccc1"]
        self.single_smiles = ["CCO"]
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_extract_smarts291_batch_exists(self):
        """Test that ExtractSMARTS291Batch function exists."""
        self.assertTrue(hasattr(rdMD, 'ExtractSMARTS291Batch'))
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_get_smarts291_names_exists(self):
        """Test that GetSMARTS291FeatureNames function exists."""
        self.assertTrue(hasattr(rdMD, 'GetSMARTS291FeatureNames'))
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_get_smarts291_names_count(self):
        """Test GetSMARTS291FeatureNames returns 291 names."""
        names = rdMD.GetSMARTS291FeatureNames('A')
        names_list = to_list(names)
        self.assertEqual(len(names_list), 291,
                        f"Expected 291 names, got {len(names_list)}")
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_extract_smarts291_batch_single(self):
        """Test ExtractSMARTS291Batch with single molecule."""
        result = rdMD.ExtractSMARTS291Batch(self.single_smiles, 'A', 1)
        result_list = to_list(result)
        self.assertEqual(len(result_list), 1)
        self.assertEqual(len(to_list(result_list[0])), 291,
                        f"Expected 291 features, got {len(to_list(result_list[0]))}")
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_extract_smarts291_batch_multiple(self):
        """Test ExtractSMARTS291Batch with multiple molecules."""
        result = rdMD.ExtractSMARTS291Batch(self.smiles_list, 'A', 0)
        result_list = to_list(result)
        self.assertEqual(len(result_list), 5)
        for mol_features in result_list:
            mol_features_list = to_list(mol_features)
            self.assertEqual(len(mol_features_list), 291)
            
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_extract_smarts291_batch_threading(self):
        """Test ExtractSMARTS291Batch with different n_jobs values."""
        # Test with n_jobs=0 (auto)
        result_auto = rdMD.ExtractSMARTS291Batch(self.smiles_list, 'A', 0)
        result_auto_list = to_list(result_auto)
        self.assertEqual(len(result_auto_list), 5)
        
        # Test with n_jobs=1 (single thread)
        result_single = rdMD.ExtractSMARTS291Batch(self.smiles_list, 'A', 1)
        result_single_list = to_list(result_single)
        self.assertEqual(len(result_single_list), 5)
        
        # Test with n_jobs=2 (two threads)
        result_two = rdMD.ExtractSMARTS291Batch(self.smiles_list, 'A', 2)
        result_two_list = to_list(result_two)
        self.assertEqual(len(result_two_list), 5)
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_extract_smarts291_batch_consistency_with_single(self):
        """Test that batch results match single molecule results."""
        # Get batch results
        batch_result = rdMD.ExtractSMARTS291Batch(self.single_smiles, 'A', 1)
        batch_features = to_list(to_list(batch_result)[0])
        
        # Get single molecule result
        mol = Chem.MolFromSmiles("CCO")
        single_features = to_list(rdMD.CalcAbrahamFeatures(mol))
        
        # Compare
        self.assertEqual(len(batch_features), len(single_features))
        for i, (b, s) in enumerate(zip(batch_features, single_features)):
            if b == b and s == s:  # Both not NaN
                self.assertAlmostEqual(b, s, places=6,
                    msg=f"Feature {i} mismatch: batch={b}, single={s}")
                    
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_extract_smarts291_batch_empty_smiles(self):
        """Test ExtractSMARTS291Batch handles empty strings."""
        smiles_with_empty = ["CCO", "", "c1ccccc1"]
        result = rdMD.ExtractSMARTS291Batch(smiles_with_empty, 'A', 0)
        result_list = to_list(result)
        self.assertEqual(len(result_list), 3)


class TestRDKit217(unittest.TestCase):
    """Test class for RDKit217 descriptors functionality."""

    def setUp(self):
        """Set up test molecules."""
        self.ethanol = Chem.MolFromSmiles("CCO")
        self.benzene = Chem.MolFromSmiles("c1ccccc1")
        self.acetic_acid = Chem.MolFromSmiles("CC(=O)O")
        self.mols = [self.ethanol, self.benzene, self.acetic_acid]
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_get_rdkit217_names_exists(self):
        """Test that GetRDKit217DescriptorNames function exists."""
        self.assertTrue(hasattr(rdMD, 'GetRDKit217DescriptorNames'))
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_get_rdkit217_names_count(self):
        """Test GetRDKit217DescriptorNames returns 217 names."""
        names = rdMD.GetRDKit217DescriptorNames()
        names_list = to_list(names)
        self.assertEqual(len(names_list), 217,
                        f"Expected 217 names, got {len(names_list)}")
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_get_rdkit217_names_content(self):
        """Test GetRDKit217DescriptorNames contains expected descriptors."""
        names = rdMD.GetRDKit217DescriptorNames()
        names_list = to_list(names)
        expected_names = ['MolWt', 'HeavyAtomMolWt', 'MaxEStateIndex', 
                          'MinEStateIndex', 'MaxAbsEStateIndex']
        for name in expected_names:
            self.assertIn(name, names_list, f"Expected descriptor {name} not found")
            
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_extract_rdkit_descriptors_exists(self):
        """Test that ExtractRDKitDescriptorsFromMolsBatch function exists."""
        self.assertTrue(hasattr(rdMD, 'ExtractRDKitDescriptorsFromMolsBatch'))
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_extract_rdkit_descriptors_single(self):
        """Test ExtractRDKitDescriptorsFromMolsBatch with single molecule."""
        result = rdMD.ExtractRDKitDescriptorsFromMolsBatch([self.ethanol], 1)
        self.assertEqual(len(result), 1)
        self.assertEqual(len(result[0]), 217,
                        f"Expected 217 features, got {len(result[0])}")
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_extract_rdkit_descriptors_batch(self):
        """Test ExtractRDKitDescriptorsFromMolsBatch with multiple molecules."""
        result = rdMD.ExtractRDKitDescriptorsFromMolsBatch(self.mols, 0)
        self.assertEqual(len(result), 3)
        for mol_features in result:
            self.assertEqual(len(mol_features), 217)
            
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_extract_rdkit_descriptors_with_none(self):
        """Test ExtractRDKitDescriptorsFromMolsBatch handles None molecules."""
        mols_with_none = [self.ethanol, None, self.benzene]
        result = rdMD.ExtractRDKitDescriptorsFromMolsBatch(mols_with_none, 0)
        self.assertEqual(len(result), 3)
        # First and third should have valid features
        self.assertEqual(len(result[0]), 217)
        self.assertEqual(len(result[2]), 217)
        # Second (None) should have zeros
        self.assertEqual(len(result[1]), 217)
        
    @unittest.skipIf(rdMD.HasOsmordredSupport() == False, "No osmordred support")
    def test_extract_rdkit_descriptors_consistency(self):
        """Test ExtractRDKitDescriptorsFromMolsBatch returns consistent results."""
        result1 = rdMD.ExtractRDKitDescriptorsFromMolsBatch([self.ethanol], 1)
        result2 = rdMD.ExtractRDKitDescriptorsFromMolsBatch([self.ethanol], 1)
        # Compare feature values
        for f1, f2 in zip(result1[0], result2[0]):
            if f1 == f1 and f2 == f2:  # Both not NaN
                self.assertAlmostEqual(f1, f2, places=6)


if __name__ == '__main__':
    unittest.main() 
