#!/usr/bin/python
import os, re, math
from itertools import zip_longest

#Basic utilities used by simplepdb or prepareamber that are generally useful and
#not limited in their functionality to the simplepdb class

#The PDB file format is fixed width; the field widths, columns where each field
#starts (numbered from 0), and data assigned to the fields are as shown.
#Fieldwidth indices 2 6 10 16 are padding fields associated with 
#columns not assigned to anything by the spec
pdb_fieldwidths = (6, 5, -1, 4, 1, 3, -1, 1, 4, 1, -3, 8, 8, 8, 6, 6, -10, 2, 2)
pdb_fieldstarts = (0, 7, 13, 17, 18, 22, 23, 27, 31, 39, 47, 55, 61, 77, 79)
pdb_fieldnames = ('recordname', 'atomnum', 'atomname', 'altloc',
        'resname', 'chainid', 'resnum', 'rescode', 'x', 'y', 'z', 'occupancy',
        'beta', 'element', 'charge')
pdb_floatfields = (8, 9, 10, 11, 12)
pdb_intfields = (1, 6)
pdb_connectfields = (6, 5, 5, 5, 5, 5)

def get_dist(r1, r2):
    '''
    get distance between two atoms
    '''
    lens = []
    for val in [r1, r2]:
        try:
            lens.append(len(val))
        except TypeError as e:
            print(e)
            sys.exit()
    assert lens[0] == lens[1]
    total = 0
    for i in range(lens[0]):
        total += (r1[i] - r2[i]) * (r1[i] - r2[i])
    return math.sqrt(total)

def is_bonded(dist, thresh=2.0):
    '''
    bonded if dist < thresh, default from chodera lab's miniAmber which I am
    hoping is the LEAP criterion
    '''
    if dist < thresh:
        return True
    else:
        return False

def walk_to_end(atom, atomdict, mytype):
    '''
    when you walk the adjacency list in atomdict, do you make it to the end?
    '''
    at_end = True
    if atom[1] > 0:
        atom_index = atom[1]
        next_atom_type = atom[2]
        at_end = walk_to_end(atomdict[next_atom_type][atom_index], atomdict,
                next_atom_type)
    elif mytype == "O" or mytype == "OXT":
        return at_end
    else:
        at_end = False
    return at_end

def accumulate(iterable):
    '''
    Generate running sum from an iterator
    '''	
    total = next(iterable)
    yield total 
    for value in iterable:
    	total += value
    	yield total 

def make_parser(pdb_fieldwidths):
    '''
    Return a fixed-width file parser
    '''
    cuts = tuple(cut for cut in accumulate(abs(fw) for fw in pdb_fieldwidths))
    pads = tuple(fw < 0 for fw in pdb_fieldwidths)
    flds = tuple(zip_longest(pads, (0,)+cuts, cuts))[:-1]
    parse = lambda line: [line[i:j].rstrip('\n') for pad, i, j in flds if not pad]
    return parse

def get_base(fname):
    '''
    Strip path and extension from filename in order to generate other filenames
    '''
    return os.path.splitext(os.path.basename(fname))[0]

def get_fname(fname):
    '''
    Generate a filename that doesn't overwrite anything in the current
    directory
    '''
    i = 0
    name = fname.split('_')
    if len(name) > 1:
        base = '_'.join(name[:-1])
        ext = '_' + name[-1]
    else:
        base,ext = os.path.splitext(fname)
    while 1:
        if i == 0:
            fname = base + ext
        else:
            fname = base + str(i) + ext
        if os.path.isfile(fname):
            if i == 0:
                base = base + '_'
            i += 1
        else:
            return fname

def get_molname(molname):
    '''
    Generate a unique molname, important if there are multiple mols
    '''
    molname = molname[:3].upper()
    while os.path.isfile(molname + '.lib'):
        num = ''.join(m for m in molname if m.isdigit())
        num = int(num)+1 if num else 1
        numstr = str(num)
        molname = ''.join(m for m in molname if m.isalpha())
        molname = molname[:3-len(numstr)] + numstr
    return molname

def get_libs(ff):
    '''
    Get the libs that will be loaded by a leaprc; force field should be
    provided as a full path
    '''
    libs = []
    with open(ff, 'r') as f:
        if 'oldff' in os.path.realpath(f.name): #back up 3x
            path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(f.name))))+'/lib/'
        else: #backup 2x
            path = os.path.dirname(os.path.dirname(os.path.realpath(f.name)))+'/lib/'
        for line in f:
            if line.startswith('loadOff') or line.startswith('loadoff'):
                libs.append(path + line.split()[-1])

    return libs

def is_secret_peptide(mol_data):
    '''
    Checks whether something that wasn't identified as a peptide (i.e. didn't
    have properly named amino acid residues) might have a protein backbone and
    might need residues renamed (not by us, for now). Infers basic connectivity
    with atomic distances, returns True if there's evidence of even one peptide
    bond. Assumes chain is laid out in N->O order
    '''
    #relevant elements for identifying backbone
    nameslist = ["N", "C", "O"]
    #store tuple of locations of atoms for each name and whether they are
    #bonded to the correct next guy in the chain here
    atomdict = {}
    for name in nameslist:
        atomdict[name] = []
    #which ALREADY SEEN atom name(s) should be bonded to current atom if it's in a peptide
    last_atom = {}
    last_atom["C"] = ["C", "N", "O"]
    last_atom["N"] = ["C"]
    last_atom["O"] = ["C"]
    for i,atom in enumerate(mol_data.mol_data['element']):
        reduced_type = ''.join(char for char in atom if char.isalpha())
        if reduced_type in nameslist:
            me = (mol_data.mol_data['x'][i], mol_data.mol_data['y'][i],
                    mol_data.mol_data['z'][i])
            atomdict[reduced_type].append([me, -1, ''])
            be_my_neighbor = last_atom[reduced_type]
            for maybe_neighbor in be_my_neighbor:
                for j,seen in enumerate(atomdict[maybe_neighbor]):
                    if maybe_neighbor == reduced_type and j == len(atomdict[reduced_type])-1:
                        #we're at the current atom
                        break
                    elif is_bonded(get_dist(me, seen[0])) and atomdict[maybe_neighbor][j][1] < 0:
                        #update my relevant neighbor to indicate our adjacency
                        atomdict[maybe_neighbor][j][1] = len(atomdict[reduced_type]) - 1
                        atomdict[maybe_neighbor][j][2] = reduced_type
    result = False
    for atom in atomdict["N"]:
        result = walk_to_end(atom, atomdict, "N")
        if result:
            break
    return result

def get_units(lib):
    '''
    Get the names of the units defined by an AMBER library file
    '''
    units = []
    ext = os.path.splitext(lib)[-1]
    with open(lib,'r') as f:
        if ext[0:4] != 'prep':
            copy = False
            for line in f:
                if line.startswith('!!index array str'):
                    copy = True
                elif line.startswith('!'):
                    break
                elif copy:
                    units.append(line.strip().strip('"'))
        else:
            i = 0
            for line in f:
                if i == 4:
                    units.append(line.split()[0])
                    break
                i += 1
    return units

def get_charge(mol2):
    '''
    Get the net charge on a molecule from a structure file. Currently only
    works with mol2 format
    '''
    with open(mol2, 'r') as f:
        copy = False
        charge = 0.0
        for line in f:
            if line.startswith('@'):
                record = line.split('>')[-1]
                if record.strip() == 'ATOM':
                    copy = True 
                elif record != 'ATOM':
                    if copy == True:
                        break
                    else:
                        continue
            elif (copy == True):
                charge += float(line.split()[-1])
    return int(math.ceil(charge)) if charge > 0 else int(math.floor(charge))

def get_ions(libs):
    '''
    Find which ions are defined for chosen water model; returns a dict mapping
    the AMBER residue name to the element
    '''
    ions = {}
    for lib in libs:
        copy = False
        with open(lib, 'r') as f:
            for line in f:
                if line.startswith('MASS'):
                    copy = True
                elif line.startswith('NONBON'):
                    break
                elif copy:
                    contents = line.split()
                    if contents:
                        ambername = contents[0]
                        element = ''.join(c for c in contents[0] if c.isalpha())
                        ions[ambername] = element.upper()
    return ions

def get_available_res(ff=''):
    '''
    Return a set of the standard amino acid residues (plus water) defined by
    AMBER; specific to your system if you pass the force field you want to use,
    otherwise a generic set is returned based on ff14SB
    '''
    if ff and not os.path.isfile(ff):
        print('Force field not found.\n')
        ff = ''
    if ff:
        units = []
        libs = get_libs(ff)
        for lib in libs:
            units += get_units(lib)
        #include residues from atomic_ions.lib
        units += ['WAT', 'HOH', "AG", "AL", "Ag", "BA", "BR", "Be", "CA", "CD", 
        "CE", "CL", "CO", "CR", "CS", "CU", "CU1", "Ce", "Cl-", "Cr", "Dy", 
        "EU", "EU3", "Er", "F", "FE", "FE2", "GD3", "H3O+", "HE+", "HG", "HZ+", 
        "Hf", "IN", "IOD", "K+", "K", "LA", "LI", "LU", "MG", "MN", "NA", "NH4", 
        "NI", "Na+", "Nd", "PB", "PD", "PR", "PT", "Pu", "RB", "Ra", "SM", "SR", 
        "Sm", "Sn", "TB", "TL", "Th", "Tl", "Tm", "U4+", "V2+", "Y", "YB2", "ZN", "Zr"]
        units = set(units)
        with open(ff, 'r') as f:
            for line in f: # look for aliased residues
                m = re.search(r'(\S+)\s*=\s*(\S+)', line)
                if m and m.group(2) in units:
                    units.add(m.group(1))
        return units
    #amino acid residues that leap should recognize with a standard protein force
    #field, plus water
    return set(['ALA', 'GLY', 'SER', 'THR', 'LEU', 'ILE', 'VAL', 'ASN', 'GLN',
        'ARG', 'HID', 'HIE','HIS','HIP', 'TRP', 'PHE', 'TYR', 'GLU', 'ASP', 'LYS', 
        'LYN', 'PRO', 'CYS', 'CYX', 'MET', 'ASH', 'GLH', 'CYM', 'HYP', 'NALA', 
        'NGLY', 'NSER', 'NTHR', 'NLEU', 'NILE', 'NVAL', 'NASN', 'NGLN', 'NARG', 
        'NHID', 'NHIE', 'NHIP', 'NTRP', 'NPHE', 'NTYR', 'NGLU', 'NASP', 'NLYS', 
        'NPRO', 'NCYS', 'NCYX', 'NMET', 'NME', 'NHE', 'ACE', 'CALA', 'CGLY', 
        'CSER', 'CTHR', 'CLEU', 'CILE', 'CVAL', 'CASN', 'CGLN', 'CARG', 'CHID', 
        'CHIE', 'CHIP', 'CTRP', 'CPHE', 'CTYR', 'CGLU', 'CASP', 'CLYS', 'CPRO', 
        'CCYS', 'CCYX', 'CMET', 'CHYP']+ ['WAT', 'HOH', "AG", "AL", "Ag", "BA", 
        "BR", "Be", "CA", "CD", "CE", "CL", "CO", "CR", "CS", "CU", "CU1", "Ce", 
        "Cl-", "Cr", "Dy", "EU", "EU3", "Er", "F", "FE", "FE2", "GD3", "H3O+", 
        "HE+", "HG", "HZ+", "Hf", "IN", "IOD", "K+", "K", "LA", "LI", "LU", 
        "MG", "MN", "NA", "NH4", "NI", "Na+", "Nd", "PB", "PD", "PR", "PT", 
        "Pu", "RB", "Ra", "SM", "SR", "Sm", "Sn", "TB", "TL", "Th", "Tl", "Tm", 
        "U4+", "V2+", "Y", "YB2", "ZN", "Zr"])

#ions and cofactors that have parameters at the AMBER parameter database
#http://research.bmh.manchester.ac.uk/bryce/amber
#and a selection of other common molecules like SO4 that also have widely
#available parameters that aren't distributed with AMBER
common_ions_and_cofactors = set(['GDP', 'GTP', 'ADP', 'ATP', 'FMN', 'FAD',
    'NAD', 'NADH', 'NAH', 'NDP', 'NPD', 'NPH', 'ARP', 'HEM', 'HEME', 'MG', 
    'CO6', 'CAL', 'MNG', 'SO4', 'PO4'])

#modified amino acid residues that have parameters at the AMBER parameter
#database http://research.bmh.manchester.ac.uk/bryce/amber
modified_residues = set(['ZA', 'ZC', 'ZD', 'ZE', 'ZF', 'ZG', 'ZHD', 'ZHE',
    'ZHP', 'ZI', 'ZK', 'ZL', 'ZM', 'ZN', 'ZP', 'ZQ', 'ZR', 'ZS', 'ZT', 'ZV', 
    'ZW', 'ZY', 'M3L', 'CC4', 'CH4','K3M', 'TFL', 'HFL', 'NOR', 'ORN',
    'HEP', 'H1D', 'H2D', 'H1E', 'H2E', 'S1P', 'S2P', 'T1P', 'T2P', 'Y1P',
    'Y2P', 'SEP', 'THP', 'TYP'])

