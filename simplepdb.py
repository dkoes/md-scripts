#!/usr/bin/python
import pdb_util as util
from copy import deepcopy
import os, itertools
import collections

class simplepdb:
    '''
    Parses and writes PDB files, and exposes limited functionality for
    manipulating their contents with a particular focus on the kinds of
    manipulations required for setting up MD simulations. 

    Attributes:
        mol_data: Dictionary of PDB column names and their values.

        ters : Locations of breaks in the molecule, per the input PDB or
        resulting from simple operations such as merging molecules. Specified
        as a list of residues that appear immediately _before_ a break.

        connect: Connect records for the molecule. Format is a dictionary where
        the keys are the atoms that were found in the input connect record and
        the values are the list of atoms to which they were bonded (per that
        record). 

        natoms: Number of atoms in molecule(s).
    '''
    def __init__(self, other):
        '''
        Return a simplepdb object created by parsing an input PDB file or
        copying the contents of another object.
        Can't construct an object without such input because no utilities are 
        provided that could be used to construct a reasonable molecule.
        '''
        if isinstance(other, self.__class__):
            for k,v in other.__dict__.iteritems():
                setattr(self, k, deepcopy(v))
        else:
            assert os.path.isfile(other), 'simplepdb constructor requires \
            input PDB or object of the same type.\n'
            assert os.path.splitext(other)[-1] == '.pdb', 'Not a PDB file.\n'
            self.mol_data = self.parse_pdb(other)
            self.ters,self.connect = self.get_ters_and_connect(other)
            self.natoms = len(self.mol_data['atomnum'])

    def __eq__(self, other):
        '''
        Check equality of simplepdb objects based on the values of their fields
        '''
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return NotImplemented

    def parse_pdb(self, pdb):
        '''
        Return a dictionary of PDB column names and their values for ATOM and
        HETATM records in the provided PDB file
        '''
        #TODO: deal with multiple models
        parse = util.make_parser(util.pdb_fieldwidths)
        f = open(pdb, 'r')
        mol_data_list = [parse(line) for line in f if line.startswith('HETATM') or line.startswith('ATOM')]
        f.close()
        mol_data = {}
        for i,field in enumerate(util.pdb_fieldnames):
            if i in util.pdb_floatfields:
                fieldlist = [float(line[i]) if line[i].strip() else line[i] for
                        line in mol_data_list]
            elif i in util.pdb_intfields:
                fieldlist = [int(line[i]) if line[i].strip() else line[i] for
                        line in mol_data_list]
            else:
                fieldlist = [line[i] for line in mol_data_list]
            mol_data[field] = fieldlist
        if not mol_data['element']:
            mol_data = set_element(mol_data)
        return mol_data

    def get_ters_and_connect(self, pdb):
        '''
        Returns a list of breaks in a PDB file and a dictionary of CONECT records
        '''
        ters = []
        connect = collections.OrderedDict()
        last_line = ''
        parse = util.make_parser(util.pdb_connectfields)
        with open (pdb,'r') as f:
            for line in f:
                if line.startswith('TER'):
                    ter = line[23:27].strip()
                    if not ter:
                        ter = last_line[23:27].strip()
                    if ter: ters.append(int(ter)) 
                elif line.startswith('CONECT'):
                    contents = parse(line)
                    atom = int(contents[1])
                    bonds = [int(bond) for bond in contents[2:] if
                            bond.strip()]
                    if atom not in connect:
                        connect[atom] = bonds
                    else:
                        connect[atom] = connect[atom] + bonds
                elif line.startswith('ATOM') or line.startswith('HETATM'):
                    last_line = line
        ter = last_line[23:27].strip()
        if ter and not ter in ters:
            ters.append(ter)
        return ters,connect

    def group_by_residue(self):
        '''
        Rearrange atoms in a file so that atoms in the same residue are
        contiguous and orders residues monotonically by resnum
        '''
        unsorted_resmap = {}
        for old_idx in range(self.natoms):
            resnum = self.mol_data['resnum'][old_idx]
            if resnum not in unsorted_resmap:
                unsorted_resmap[resnum] = [old_idx]
            else:
                unsorted_resmap[resnum].append(old_idx)
        resmap = collections.OrderedDict(sorted(unsorted_resmap.items(), key=lambda t: t[0]))
        new_indices = list(itertools.chain.from_iterable(resmap.values()))
        new_mol_data = {}
        for key in self.mol_data:
            new_mol_data[key] = [self.mol_data[key][i] for i in new_indices]

        self.renumber_atoms()
        self.mol_data = new_mol_data
    
    def renumber_atoms(self, start_val=1):
        '''
        Renumber atoms so they start at start_val
        '''
        mapping = {}
        for i in range(self.natoms):
            old_val = self.mol_data['atomnum'][i]
            new_val = i + start_val
            self.mol_data['atomnum'][i] = new_val
            mapping[old_val] = new_val

        #TODO: ugly
        new_connect = collections.OrderedDict()
        for atom,bonds in self.connect.items():
            if atom in mapping:
                new_connect[mapping[atom]] = bonds
        
        for atom,bonds in new_connect.items():
            for bond in bonds:
                if bond in mapping:
                    idx = bonds.index(bond)
                    new_connect[atom][idx] = mapping[bond]

        self.connect = new_connect
    
    def renumber_residues(self, start_val=1):
        '''
        Renumber residues so they start at start_val in "first seen" order, desirable
        when there is a ligand at the end of data with an out-of-order resnum
        '''
        reslist = []
        for i,resnum in enumerate(self.mol_data['resnum']):
            name = self.mol_data['resname'][i] + str(resnum) 
            if name not in reslist:
                newidx = len(reslist)
                reslist.append(name)
            else:
                newidx = reslist.index(name)
            newnum = newidx + start_val
            self.mol_data['resnum'][i] = newnum
            if resnum in self.ters:
                ter_idx = self.ters.index(resnum)
                self.ters[ter_idx] = newnum 
    
    def rename_atoms(self):
        '''
        Generate unique atom names
        '''
        if self.has_unique_names():
            return
        for i,name in enumerate(self.mol_data['atomname']):
            self.mol_data['atomname'][i] = ''.join([char for char in
                    self.mol_data['element'][i]])
        
        occurrences = {}
        for i,atom in enumerate(self.mol_data['atomname']):
            if atom not in occurrences:
                occurrences[atom] = [i,1]
            else:
                occurrences[atom][1] += 1
            self.mol_data['atomname'][i] += str(occurrences[atom][1])
            self.mol_data['atomname'][i] = \
            '{:>{}s}'.format(self.mol_data['atomname'][i],
                    util.pdb_fieldwidths[3])
    
    def set_element(self):
        '''
        Set atom element based on atom name, but only if element not set.
        '''
        if not self.mol_data['element']:
            for i,name in enumerate(self.mol_data['atomname']):
                element = ''.join([char for char in name if char.isalpha()])
                self.mol_data['element'][i] = '{:>{}s}'.format(element,
                        util.pdb_fieldwidths[-2])

    def sanitize(self):
        '''
        Perform atom renumbering, residue renumbering, and regrouping atoms so
        residues are contiguous; if a small molecule, also uniquely names atoms
        and sets the element field
        '''
        self.group_by_residue()
        self.renumber_atoms()
        self.renumber_residues()
        if not self.is_protein():
            self.set_element()
            self.rename_atoms()

    def has_hydrogen(self):
        '''
        Returns true if hydrogens are present
        '''
        return 'H' in [elem.strip() for elem in self.mol_data['element']]
    
    def strip_hydrogen(self):
        '''
        Strip out all the hydrogens
        '''
        h_indices = [i for i,elem in enumerate(self.mol_data['element']) if elem.strip() ==
                'H']
        new_mol_data = {}
        for key in self.mol_data.keys():
            new_mol_data[key] = [self.mol_data[key][i] for i in
                    range(len(self.mol_data[key])) if i not in h_indices]
        self.mol_data = new_mol_data
    
    def is_protein(self, ff=''):
        '''
        Returns true if standard amino acid residues are present
        '''
        aa = util.get_available_res(ff).intersection(self.mol_data['resname'])
        return len(aa) > 0

    def has_unique_names(self):
        '''
        Returns true if atom names are unique
        '''
        #TODO: add to tests
        atom_ids = []
        for i in range(self.natoms):
            atom_ids.append(str(self.mol_data['resnum'][i]) +
                    self.mol_data['atomname'][i])
        counter = collections.Counter(atom_ids)
        if any(t > 1 for t in counter.values()):
            return False
        return True

    def set_recordname(self, newname, resnum=None):
        '''
        Set record name to ATOM or HETATM for residue number resnum or all
        resnums if no number is provided
        '''
        #TODO: add to tests
        assert newname=='ATOM' or newname=='HETATM', 'Record names must be one \
        of "ATOM" and "HETATM"'
        if not resnum:
            self.mol_data['recordname'] = [newname] * self.natoms
        else:
            self.mol_data['recordname'] = [newname for name in
            self.mol_data['recordname'] if self.mol_data['resnum'] == resnum]

    def set_resname(self, newname, oldname=''):
        '''
        Set resname to newname; if oldname is not specified then all resnames
        are updated to newname, otherwise just oldname is
        '''
        #TODO: add to tests
        if not oldname:
            self.mol_data['resname'] = [newname] * self.natoms
        else:
            self.mol_data['resname'] = [newname for name in
            self.mol_data['resname'] if name == oldname]

    def writepdb(self, fname, mols=[]):
        '''
        Write molecule data to a file; takes filename and a list of mols to be
        written. If no list is provided, just the calling molecule is written,
        otherwise only the molecules in the list are written. 
        '''
        if not mols:
            mols = [self]

        with open(fname, 'w') as f:
            start_atom = 1
            start_res = 1
            for mol in mols:
                mol.renumber_atoms(start_atom)
                mol.renumber_residues(start_res)
                for i in range(mol.natoms):
                    j = 0
                    for fieldwidth in util.pdb_fieldwidths:
                        if fieldwidth > 0:
                            fieldname = util.pdb_fieldnames[j]
                            if fieldname=='x' or fieldname=='y' or fieldname=='z':
                                output='{:.3f}'.format(mol.mol_data[fieldname][i])
                                f.write('{:>{}s}'.format(str(output),fieldwidth))
                            elif fieldname=='occupancy' or fieldname=='beta' and \
                            str(mol.mol_data[fieldname][i]).strip():
                                output='{:.2f}'.format(mol.mol_data[fieldname][i])
                                f.write('{:>{}s}'.format(str(output),fieldwidth))
                            else:
                                f.write('{:>{}s}'.format(str(mol.mol_data[fieldname][i]),fieldwidth))
                            j += 1
                        else:
                            f.write('{:>{}s}'.format('',abs(fieldwidth)))
                    f.write('\n')
                    if (i == mol.natoms-1):
                        f.write('TER\n')
                    elif (mol.mol_data['resnum'][i] in mol.ters and
                    mol.mol_data['resnum'][i] != mol.mol_data['resnum'][i+1]):
                        f.write('TER\n')
                start_atom = mol.mol_data['atomnum'][-1]+1
                start_res = mol.mol_data['resnum'][-1]+1

            for mol in mols:
                for atom,bonds in mol.connect.items():
                    bonds_seen = 0
                    for bond in bonds:
                        if not bonds_seen % 4:
                            if bonds_seen: f.write('\n')
                            f.write('CONECT')
                            f.write('{:>{}s}'.format(str(atom),util.pdb_connectfields[1]))
                        f.write('{:>{}s}'.format(str(bond),5)) 
                        bonds_seen = bonds_seen + 1
                    f.write('\n')

            f.write('END\n')

