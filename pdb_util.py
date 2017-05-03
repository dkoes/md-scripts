#!/usr/bin/python
import os
from itertools import izip_longest

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
    flds = tuple(izip_longest(pads, (0,)+cuts, cuts))[:-1]
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
    base,ext = fname.split('_')
    if ext: ext = '_' + ext
    while 1:
        if i == 0:
            fname = base + ext
        else:
            fname = base + str(i) + ext
        if os.path.isfile(fname):
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
        path = os.path.dirname(os.path.dirname(os.path.realpath(f.name)))+'/lib/'
        for line in f:
            if line.startswith('loadOff'):
                libs.append(path + line.split()[-1])

    return libs

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
        charge = 0
        for line in f:
            if line.startswith('@'):
                record = line.split('>')[-1]
                if record == 'ATOM':
                    copy = True 
                elif record != 'ATOM':
                    if copy == True:
                        break
                    else:
                        continue
            elif (copy == True):
                charge += float(line.split()[-1])
    return int(round(charge))


def get_available_res(ff=''):
    '''
    Return a set of the standard amino acid residues (plus water) defined by
    AMBER; specific to your system if you pass the force field you want to use,
    otherwise a generic set is returned based on ff14SB
    '''
    if ff and not os.path.isfile(ff):
        print 'Force field not found.\n'
        ff = ''
    if ff:
        units = []
        libs = get_libs(ff)
        for lib in libs:
            units += get_units(lib)
        units += ['WAT', 'HOH']
        return set(units)
    #amino acid residues that leap should recognize with a standard protein force
    #field, plus water
    return set(['ALA', 'GLY', 'SER', 'THR', 'LEU', 'ILE', 'VAL', 'ASN', 'GLN',
        'ARG', 'HID', 'HIE','HIP', 'TRP', 'PHE', 'TYR', 'GLU', 'ASP', 'LYS', 
        'LYN', 'PRO', 'CYS', 'CYX', 'MET', 'ASH', 'GLH', 'CYM', 'HYP', 'NALA', 
        'NGLY', 'NSER', 'NTHR', 'NLEU', 'NILE', 'NVAL', 'NASN', 'NGLN', 'NARG', 
        'NHID', 'NHIE', 'NHIP', 'NTRP', 'NPHE', 'NTYR', 'NGLU', 'NASP', 'NLYS', 
        'NPRO', 'NCYS', 'NCYX', 'NMET', 'NME', 'NHE', 'ACE', 'CALA', 'CGLY', 
        'CSER', 'CTHR', 'CLEU', 'CILE', 'CVAL', 'CASN', 'CGLN', 'CARG', 'CHID', 
        'CHIE', 'CHIP', 'CTRP', 'CPHE', 'CTYR', 'CGLU', 'CASP', 'CLYS', 'CPRO', 
        'CCYS', 'CCYX', 'CMET', 'CHYP', 'WAT', 'HOH'])

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

