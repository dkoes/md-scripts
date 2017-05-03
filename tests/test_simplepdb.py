#!/usr/bin/python
import unittest
import env
import simplepdb as pdb
import os
from plumbum.cmd import awk, head, tail

class IOTests(unittest.TestCase):
    '''
    Tests PDB file parsing and writing.
    '''
    def setUp(self):
        self.pdb = pdb.simplepdb('chignolin.pdb')
        self.maxDiff = None

    def test_natoms(self):
        self.assertEqual(self.pdb.natoms, 138)

    def test_ters(self):
        self.assertEqual(self.pdb.ters[0], 10)

    def test_recordname(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,0,6)}'])()
        #plumbum converts encoding to unicode and also adds an empty line at
        #the end...
        out = str(out)
        out = out.split('\n')
        self.assertEqual(self.pdb.mol_data['recordname'], out[:-1])

    def test_atomnum(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,7,5)}'])()
        out = out.split('\n')
        out = [int(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['atomnum'], out[:-1])

    def test_atomname(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,13,4)}'])()
        out = str(out)
        out = out.split('\n')
        self.assertEqual(self.pdb.mol_data['atomname'], out[:-1])

    def test_altloc(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,17,1)}'])()
        out = str(out)
        out = out.split('\n')
        self.assertEqual(self.pdb.mol_data['altloc'], out[:-1])

    def test_resname(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,18,3)}'])()
        out = str(out)
        out = out.split('\n')
        self.assertEqual(self.pdb.mol_data['resname'], out[:-1])

    def test_chainid(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,22,1)}'])()
        out = str(out)
        out = out.split('\n')
        self.assertEqual(self.pdb.mol_data['chainid'], out[:-1])

    def test_resnum(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,23,4)}'])()
        out = out.split('\n')
        out = [int(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['resnum'], out[:-1])

    def test_rescode(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,27,1)}'])()
        out = str(out)
        out = out.split('\n')
        self.assertEqual(self.pdb.mol_data['rescode'], out[:-1])

    def test_x(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,31,8)}'])()
        out = out.split('\n')
        out = [float(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['x'], out[:-1])

    def test_y(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,39,8)}'])()
        out = out.split('\n')
        out = [float(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['y'], out[:-1])

    def test_z(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,47,8)}'])()
        out = out.split('\n')
        out = [float(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['z'], out[:-1])

    def test_occupancy(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,55,6)}'])()
        out = out.split('\n')
        out = [float(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['occupancy'], out[:-1])

    def test_beta(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,61,6)}'])()
        out = out.split('\n')
        out = [float(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['beta'], out[:-1])

    def test_element(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,77,2)}'])()
        out = str(out)
        out = out.split('\n')
        self.assertEqual(self.pdb.mol_data['element'], out[:-1])

    def test_charge(self):
        out = (head['-n', 138, 'chignolin.pdb'] | awk['{print substr($0,79,2)}'])()
        out = out.split('\n')
        out = [int(elem) if elem.strip() else elem for elem in out]
        self.assertEqual(self.pdb.mol_data['charge'], out[:-1])

    def test_write(self):
        fname = 'tmp.pdb'
        self.pdb.writepdb(fname)
        new_pdb = pdb.simplepdb(fname)
        os.remove(fname)
        self.assertEqual(self.pdb, new_pdb)

class ActionTests(unittest.TestCase):
    '''
    Tests methods that manipulate molecular data.
    '''
    def setUp(self):
        self.receptor = pdb.simplepdb('receptor.pdb')
        self.ligand_h = pdb.simplepdb('LIG_h.pdb')
        self.ligand_noh = pdb.simplepdb('LIG_noh.pdb')
        self.complex = pdb.simplepdb('LIGreceptor_fixed.pdb')

    def test_group_by_residue(self):
        badres = pdb.simplepdb('LIGreceptor.pdb')
        badres.group_by_residue()
        self.assertEqual(badres.mol_data['resname'],
        self.complex.mol_data['resname'])

    def test_renumber_atoms(self):
        badnumber = pdb.simplepdb('LIGreceptor_number.pdb')
        badnumber.renumber_atoms()
        self.assertEqual(badnumber.mol_data['atomnum'],
        self.complex.mol_data['atomnum'])

    def test_renumber_residues(self):
        badnumber = pdb.simplepdb('LIGreceptor_number.pdb')
        badnumber.renumber_residues()
        self.assertEqual(badnumber.mol_data['resnum'],
        self.complex.mol_data['resnum'])

    def test_rename_atoms(self):
        ligand_h_copy = pdb.simplepdb(self.ligand_h)
        ligand_noh_copy = pdb.simplepdb(self.ligand_noh)
        ligand_h_copy.strip_hydrogen()
        ligand_noh_copy.rename_atoms()
        self.assertEqual(ligand_h_copy.mol_data['atomname'],
                ligand_noh_copy.mol_data['atomname'])

    def test_set_element(self):
        ligand_h_copy = pdb.simplepdb(self.ligand_h)
        ligand_noh_copy = pdb.simplepdb(self.ligand_noh)
        ligand_h_copy.strip_hydrogen()
        ligand_noh_copy.set_element()
        self.assertEqual(ligand_h_copy.mol_data['element'],
                ligand_noh_copy.mol_data['element'])

    def test_is_protein(self):
        self.assertTrue(self.receptor.is_protein())
        self.assertTrue(self.complex.is_protein())
        self.assertFalse(self.ligand_h.is_protein())

    def test_has_hydrogen(self):
        self.assertTrue(self.ligand_h.has_hydrogen())
        self.assertFalse(self.ligand_noh.has_hydrogen())

    def test_strip_hydrogen(self):
        ligand_copy = pdb.simplepdb(self.ligand_h)
        ligand_copy.strip_hydrogen()
        self.assertFalse(ligand_copy.has_hydrogen())

if __name__ == '__main__':
    unittest.main()
