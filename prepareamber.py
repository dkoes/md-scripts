#!/usr/bin/env python3

import argparse, re
import simplepdb as pdb
import pdb_util as util
from reparm_ligand import reparm  
import os, shutil, glob, sys, logging
from plumbum import FG, TEE
from plumbum.cmd import sed, grep, cut, uniq, wc
try:
    from plumbum.cmd import obabel
except ImportError:
    raise ImportError('Check that obabel is on your path')
try:
    from plumbum.cmd import antechamber, pdb4amber, tleap, pmemd_cuda
except ImportError as e:
    print(e)
    raise ImportError('Check that AMBER binaries are on your path')

try:
    from plumbum.cmd import parmchk
except ImportError:
    try:
        from plumbum.cmd import parmchk2 as parmchk
    except ImportError:
        raise ImportError('Check parmchk[2] on your path')

try:
    from plumbum.cmd import match_atomname
except ImportError:
    raise ImportError('Check match_atomname (one of the antechamber-related'
            ' AmberTools packages) has been built - starting with AmberTools18 it'
            ' is not built by default, but can be built manually from' 
            ' $AMBERHOME/AmberTools/src/antechamber/')

class Tee(object):
    '''For runfile, duplicate stdout to file'''
    def __init__(self, name=None):
        if name:
            self.file = open(name, 'w')
        else:
            self.file = None
    def writeln(self, data):
        print(data)
        if self.file:
            self.file.write(str(data)+'\n')
            self.file.flush()
        
runfile = Tee() #a global variable so I don't have to pass it around

def find_ff(amberhome, ffname):
    ff = ''
    if os.path.isfile(amberhome + '/dat/leap/cmd/' + ffname):
        ff = amberhome + '/dat/leap/cmd/' + ffname
    elif os.path.isfile(amberhome + '/dat/leap/cmd/leaprc.protein.' + ffname):
        ff = amberhome + '/dat/leap/cmd/leaprc.protein.' + ffname            
    elif os.path.isfile(amberhome + '/dat/leap/cmd/oldff/' + ffname):
        ff = amberhome + '/dat/leap/cmd/oldff/' + ffname
    elif os.path.isfile(amberhome + '/dat/leap/cmd/oldff/leaprc.' + ffname):
        ff = amberhome + '/dat/leap/cmd/oldff/leaprc.' + ffname            
    else:
        print("Warning: force field %s not found in %s! This is likely to cause \
problems later.\n"%(ffname,amberhome))
    return ff
        
def get_cmd(input_str):
    cmd_dict = {'.frcmod' : 'loadamberparams', '.lib' : 'loadoff', '.off' : 
            'loadoff', '.prep' : 'loadAmberPrep', '.mol2' : 'loadmol2', '.pdb' :
            'loadpdb'}
    ext = os.path.splitext(input_str)[-1]
    return cmd_dict[ext]

def get_waterbox(water_model):
    if water_model == 'leaprc.water.opc':
        return ' OPCBOX '
    elif water_model == 'leaprc.water.tip3p':
        return ' TIP3PBOX '
    elif water_model == 'leaprc.water.spce':
        return ' SPCBOX '
    elif water_model == 'leaprc.water.tip4pew':
        return ' TIP4PEWBOX '

def get_water_nickname(water_model):
    if water_model == 'leaprc.water.opc':
        return 'OPC'
    elif water_model == 'leaprc.water.tip3p':
        return 'TP3'
    elif water_model == 'leaprc.water.spce':
        return 'SPC'
    elif water_model == 'leaprc.water.tip4pew':
        return 'T4E'

def make_amber_parm(fname, base, ff, molname='', water_model = '', 
        wat_dist = 0, libs=[], frcmod = '', extra=None):
    '''
    Generate AMBER parameters with tleap
    '''
    inpcrd = base + '.inpcrd'
    prmtop = base + '.prmtop'
    if not molname: molname = fname[:3].upper()

    with open(base + '.tleap', 'w') as leap_input:
        for name in ff:
            leap_input.write('source %s\n' %name)
        leap_input.write('source leaprc.gaff\n')
                
        if water_model: # source before loading structures
            water_nickname = get_water_nickname(water_model)
            leap_input.write('WAT = %s\n' %water_nickname)
            leap_input.write('HOH = %s\n' %water_nickname)
            leap_input.write('source ' + water_model + '\n')
            
        for lib in libs:
            leap_input.write(get_cmd(lib) + ' ' + lib + '\n')
        leap_input.write(molname + '=' + get_cmd(fname) + ' ' + fname + '\n')
        
        if extra: 
            extracmds = open(extra).read()
            extracmds = extracmds.replace('MOLNAME',molname)
            print(extracmds)
            leap_input.write(extracmds)
        if water_model:
            leap_input.write('solvateoct '+ molname + get_waterbox(water_model) + 
                    str(wat_dist) + '\n' + 
                    'addions '+molname+' Na+ 0\n' + 
                    'addions '+molname+' Cl- 0\n')
        elif frcmod:
            leap_input.write('loadamberparams '+frcmod+'\n' + 
                    'saveoff '+molname+' '+base+'.lib\n')
        leap_input.write('saveamberparm ' + molname+ ' ' + prmtop + ' ' + inpcrd + '\n')
        leap_input.write('quit\n')
    command = tleap['-f', base + '.tleap'] 
    runfile.writeln(command)
    command & FG

def do_amber_min_constraint(fname, base):
    '''
    Do AMBER minimization with protein constraint
    '''
    try:
        numres = (grep['ATOM', fname] | cut['-b', '23-26'] | uniq | wc)()
        numres = numres.split()[0].strip()
    except Exception as e:
        if e[1] == 1:
            numres = 0
    with open(base + '_min1.in', 'w') as min_input:
        min_input.write(base + ': initial minimization solvent + ions\n' + 
            ' &cntrl\n' + 
            '  imin   = 1,\n' + 
            '  maxcyc = 1000,\n' +
            '  ncyc   = 500,\n' +
            '  ntb    = 1,\n' +
            '  ntr    = 1,\n' +
            '  cut    = 10.0\n' +
            ' /\n' +
            'Hold the protein fixed\n' +
            '500.0\n' +
            'FIND\n' +
            '* * S *\n' +
            '* * B *\n' +
            '* * 3 *\n' +
            '* * E *\n' +
            '* * M *\n' +
            'SEARCH\n' +
            'RES 1 ' + str(numres) + '\n' +
            'END\n' +
            'END\n')
    command = pmemd_cuda['-O', '-i', base+'_min1.in', '-o', base+'_min1.out', '-p',
            base+'.prmtop', '-c', base+'.inpcrd', '-r', base+'_min1.rst',
            '-ref', base+'.inpcrd'] 
    runfile.writeln(command)
    command & FG

def do_amber_min(base):
    '''
    Do unconstrained AMBER minimization 
    '''
    with open(base + '_min2.in', 'w') as min_input:
        min_input.write(base + ':  minimization whole system\n' + 
            ' &cntrl\n' + 
            '  imin   = 1,\n' + 
            '  maxcyc = 5000,\n' + 
            '  ncyc   = 2000,\n' + 
            '  ntb    = 1,\n' + 
            '  ntr    = 0,\n' + 
            '  cut    = 10.0\n' + 
            ' /\n')
    command = pmemd_cuda['-O', '-i', base+'_min2.in', '-o', base+'_min2.out', '-p',
            base+'.prmtop', '-c', base+'_min1.rst', '-r', base+'_min2.rst'] 
    runfile.writeln(command)
    command & FG

def do_amber_warmup(fname, base, temperature):
    '''
    Do AMBER MD to gradually increase system to target temp
    '''
    try:
        numres = (grep['ATOM', fname] | cut['-b', '23-26'] | uniq | wc)()
        numres = numres.split()[0].strip()
    except Exception as e:
        if e[1] == 1:
            numres = 0
    with open(base + '_md1.in','w') as md_input:
        md_input.write(' &cntrl\n' + 
              '  imin   = 0,\n' + 
              '  irest  = 0,\n' + 
              '  ntx    = 1,\n' + 
              '  ntb    = 1,\n' + 
              '  cut    = 10.0,\n' + 
              '  ntr    = 1,\n' + 
              '  ntc    = 2,\n' + 
              '  ntf    = 2,\n' + 
              '  tempi  = 0.0,\n' + 
              '  temp0  = ' + str(temperature) + '\n' + 
              '  ntt    = 3,\n' + 
              '  gamma_ln = 1.0,\n' + 
              '  nstlim = 50000, dt = 0.002, ntxo = 2,\n' + 
              '  ntpr = 1000, ntwx = 1000, ntwr = 10000,\n' + 
              '  ioutfm = 1\n' + 
              ' /\n' + 
            'Keep prot fixed with weak restraints\n' + 
            '10.0\n' + 
            'FIND\n' + 
            '* * S *\n' + 
            '* * B *\n' + 
            '* * 3 *\n'  
            '* * E *\n' + 
            '* * M *\n' + 
            'SEARCH\n' + 
            'RES 1 ' + str(numres) + '\n' + 
            'END\n' + 
            'END\n')
    command = pmemd_cuda['-O', '-i', base+'_md1.in', '-o', base+'_md1.out', '-p',
            base+'.prmtop', '-c', base+'_min2.rst', '-r', base+'_md1.rst',
            '-ref', base+'_min2.rst', '-x', base+'_md1.nc'] 
    runfile.writeln(command)
    command & FG

def do_amber_constant_pressure(base, temp):
    '''
    Do AMBER MD to equilibrate system at constant pressure
    '''
    with open(base + '_md2.in','w') as md_input:
        md_input.write(base + ': 100ps MD\n' + 
             ' &cntrl\n' + 
             '  imin = 0, irest = 1, ntx = 7,\n' + 
             '  ntb = 2, pres0 = 1.0, ntp = 1,\n' + 
             '  taup = 2.0,\n' + 
             '  cut = 10.0, ntr = 0,\n' + 
             '  ntc = 2, ntf = 2,\n' + 
             '  tempi = '+str(temp)+', temp0 = '+str(temp)+',\n' + 
             '  ntt = 3, gamma_ln = 1.0,\n' + 
             '  nstlim = 50000, dt = 0.002, ntxo = 2,\n' + 
             '  ntpr = 5000, ntwx = 5000, ntwr = 500000,\n' + 
             '  ioutfm = 1\n' + 
             ' /\n')
    command = pmemd_cuda['-O', '-i', base+'_md2.in', '-o', base+'_md2.out', '-p',
            base+'.prmtop', '-c', base+'_md1.rst', '-r', base+'_md2.rst',
            '-x', base+'_md2.nc'] 
    runfile.writeln(command)
    command & FG

def make_amber_production_input(base, args):
    '''
    Make input files for production run AMBER MD; length is in nanoseconds
    '''
    nstlim = int(float(args.prod_length) / .000002)
    irest = int(args.keep_velocities)
    ntx = 7 if args.keep_velocities else 1
    with open(base + '_md3.in','w') as md_input:
        md_input.write(base + ': ' + str(args.prod_length) + 'ns MD\n' + 
            ' &cntrl\n' + 
            '  imin = 0, irest = '+str(irest)+', ntx = '+str(ntx)+',\n' + 
            '  ntb = 2, pres0 = 1.0, ntp = 1,\n' + 
            '  taup = 2.0,\n' + 
            '  cut = 10.0, ntr = 0,\n' + 
            '  ntc = 2, ntf = 2,\n' + 
            '  tempi = {0}, temp0 = {0},\n'.format(args.temperature) + 
            '  ntt = 3, gamma_ln = 1.0,\n' + 
            '  nstlim = '+str(nstlim)+', dt = 0.002, ntxo = 2,\n' + 
            '  ntpr = 5000, ntwx = '+str(args.coord_dump_freq)+', ntwr = 500000,\n' + 
            '  ioutfm = 1\n' + 
             '/\n')

def do_amber_preproduction(fname, base, args, ff):
    '''
    Do minimization with constraints, minimization without constraints, initial
    MD as temperature is raised to target temp, second MD where system is
    equilibrated at constant pressure, and generate input files for production
    run MD but don't run it (becuz it's PREproduction, see?)
    '''
    do_amber_min_constraint(fname, base)
    do_amber_min(base)
    do_amber_warmup(fname, base, args.temperature)
    do_amber_constant_pressure(base, args.temperature)
    make_amber_production_input(base, args)

def do_amber_production(base, dorun):
    '''
    Does AMBER production run MD locally.  If dorun is false, only print command
    '''
    command = pmemd_cuda['-O', '-i', base+'_md3.in', '-o', base+'_md3.out', '-p',
            base+'.prmtop', '-c', base+'_md2.rst', '-r', base+'_md3.rst',
            '-x', base+'_md3.nc'] 
    runfile.writeln(command)
    if dorun:
        command & FG

def do_antechamber(fname, net_charge, ff, molname, base = ''):
    '''
    Run antechamber and get correctly named versions of the following: mol2
    with bcc charges, frcmod, lib, prmtop, inpcrd
    '''
    if not base: base = util.get_base(fname)
    ext = os.path.splitext(fname)[-1]
    ext = ext.lstrip('.')
    mol2 = base + '_amber.mol2'
    #TODO: known issues with phosphates (see PDB: 2PQC) when getting the net
    #charge from Gasteiger charges computed with Open Babel
    try:
        command = antechamber['-i', fname, '-fi', ext, '-o', mol2, '-fo', 'mol2', '-c',
                'bcc', '-nc', str(net_charge), '-s', '2']
        runfile.writeln(command)
        command()
    except Exception as e:
        passed = False
        charges = []
        if net_charge != 0:
            charges.append(0)
        if net_charge != -1:
            charges.append(-1)
        for charge in charges:
            try:
                command = antechamber['-i', fname, '-fi', ext, '-o', mol2, '-fo', 'mol2', '-c',
                        'bcc', '-nc', str(charge), '-s', '2']
                runfile.writeln(command)
                command()
                passed = True
                break
            except Exception as e:
                pass
        if not passed:
            print('Antechamber failed. Check {0} structure. Aborting...\n'.format(fname))
            sys.exit()

    frcmod = base + '.frcmod'
    parmchk['-i', mol2, '-f', 'mol2', '-o', frcmod]()
    make_amber_parm(mol2, base, ff, molname=molname, frcmod=frcmod)

def set_matches(fname, libs, reslist, orphaned_res, mol, force=False):
    '''
    Find whether any units defined by a lib are required; if they are, update
    the liblist to include that lib and remove the units it defines from
    orphaned_res
    '''
    units = util.get_units(fname)
    matches = set(units).intersection(reslist)
    #TODO: make this work for prep and check bonds as well as atoms.
    #require that atom names and connectivity match before adding lib
    libatoms = []
    ext = os.path.splitext(fname)[-1][0:4]
    if ext != 'prep':
        with open(fname,'r') as f:
            atomcopy = False
            for line in f:
                if line.startswith('!'):
                    if line.split()[0].split('.')[-1] == 'atoms':
                        atomcopy = True
                    else:
                        atomcopy = False
                elif atomcopy:
                    aname = line.split()[0].strip('"')
                    if not aname.startswith("H"): #ignore hydrogen
                        libatoms.append(aname)
    molatoms = set([name.strip() for name in mol.mol_data['atomname']])
    
    if not set(libatoms).issubset(molatoms) or ext == 'prep' and not force:            
        matches = set([])
        print("Unit(s) %s defines atoms that differ from undefined residue: %s"% (' '.join(units), ' '.join(set(libatoms) - molatoms)))
    if matches:
        #redefine a unit iff found in a user-provided lib, but warn the user about the
        #duplication. don't redefine if found locally. 
        if not matches.intersection(orphaned_res) and not force:
            print('Unit %s found in %s defined previously, \
            skipping to avoid redefinition\n' % (' '.join(match for match in
                matches), fname))                
        else:
            if not matches.intersection(orphaned_res):
                print('Unit %s found in %s defined previously, \
                adding user-provided lib but redefinition may cause problems\n' \
                % (' '.join(match for match in matches), fname))
            libs.add(fname)
            frcmod = util.get_base(fname) + '.frcmod'
            if os.path.isfile(frcmod):
                libs.add(frcmod)
            orphaned_res -= matches

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generates pre-production files for AMBER \
    MD. Can handle a receptor or ligand by themselves (checks for ligand library \
    files in the current directory and generates them if they don't exist) or \
    sets up the complex if given both a receptor and ligand.")

    parser.add_argument('-s', '--structures', nargs='+', required=True, help='Structures \
    for which you want to run a simulation. N.B. if more than one is provided \
    they will be simulated together.')

    parser.add_argument('-n', '--out_name', 
    help='Optionally provide a filename prefix for the output.  For multi-structure inputs, default is "complex."')

    parser.add_argument('-p', '--libs', nargs='+', required=False, help="Optionally specify \
    a prefix for nonstandard residue library files; this can include their path if they aren't \
    in the current directory. If the relevant files exist we assume you want \
    to use them, otherwise we assume this is where you want them to go and \
    derive the residue name accordingly.")

    parser.add_argument('-w', '--water_dist', default=12, help='Water box \
    distance; defaults to 12.')

    parser.add_argument('-wm', '--water_model', default='tip3p',
            help='Water model; OPC, SPCE, TIP4P, and TIP3P are available. \
            Defaults to TIP3P.')

    parser.add_argument('-ff', '--force_field', default='leaprc.protein.ff15ipq',
            help='Force field; defaults to leaprc.protein.ff15ipq.')

    parser.add_argument('-eff', '--extra_force_field', default='',
            help='Extra force fields (e.g. DNA, lipids); defaults to null')

    parser.add_argument('-t', '--temperature', default=300, help='Simulation \
    temperature; defaults to 300K.')

    parser.add_argument('-l', '--prod_length', default=100, help='Length of \
    production-run MD. This is used to generate the input files for that MD, \
    but note that by default this MD is not run. Units are nanoseconds.')

    parser.add_argument('-k', '--keep_velocities', default=False,
            action='store_true', help='Keep velocities from preproduction run \
                    when starting the production MD. Default is false.')

    parser.add_argument('-r', '--run_prod_md', default=False,
    action='store_true', help='Run production MD locally when all preprocessing \
    is finished. Default is false, because you might want to run it on a \
    cluster.')

    parser.add_argument('-noh', '--no_touch_hyd', dest='noh', default=False,
            action='store_true', help="Don't remove any hydrogens.")

    parser.add_argument('-nc', '--net_charge', help='Optionally specify a net \
            charge for small molecule parametrization with antechamber.')

    parser.add_argument('-parm', '--parm_only', action='store_true', default =
    False, help="Only generate the necessary ligand parameters, don't do the \
    preproduction MDs")

    parser.add_argument('-ui', '--uninteractive', action='store_true',
            default=False, help="Turn off interactive mode, which exists to let \
            you check the output of pdb4amber for potentially serious problems \
            with input structures.")

    parser.add_argument('-O', '--overwrite', action='store_false',
            default=True, help='Overwrite files in the current directory that \
            are generated by this script. Default is True.')

    parser.add_argument('-df', '--coord_dump_freq', default=5000,
    help='Frequency for dumping coordinates to traj_file. Defaults to 5000. \
    The old script referred to this as the "timestep."')
    
    parser.add_argument('--extra', help="File with additional leap commands to apply")
    
    parser.add_argument('-noff', '--no_openff', action='store_true', default=False,
    help='Do not use open force field\'s SMIRNOFF to reparameterize antechamber \
    parameterized ligands. Default is False')
    
    args = parser.parse_args()

    #Check whether AMBERHOME is set and the desired force field is available
    amberhome = os.environ['AMBERHOME']
    if not amberhome:
        print("Warning: AMBERHOME is not set! This is likely to cause problems \
        later.\n")
    else:
        ff = find_ff(amberhome, args.force_field)
        if args.extra_force_field:
            nff = find_ff(amberhome, args.extra_force_field)
        else:
            nff = ''

    #Find out which ions are defined with our water model
    ion_params = []
    if args.water_model != 'opc':
        ion_params.append(amberhome + '/dat/leap/parm/frcmod.ionsjc_' +
                args.water_model)
        ion_params.append(amberhome + '/dat/leap/parm/frcmod.ions234lm_126_' +
                args.water_model)
    else:
        ion_params.append(amberhome + '/dat/leap/parm/frcmod.ionsjc_tip4pew')
        ion_params.append(amberhome + '/dat/leap/parm/frcmod.ions234lm_126_tip4pew')
    ions = util.get_ions(ion_params)

    #Find out which water model we're using
    assert args.water_model in ['opc', 'tip3p', 'spce', 'tip4'], 'Unknown water \
model %s\n' %args.water_model
    if args.water_model == 'opc':
        args.water_model = 'leaprc.water.opc'
    elif args.water_model == 'tip3p':
        args.water_model = 'leaprc.water.tip3p'
    elif args.water_model == 'spce':
        args.water_model = 'leaprc.water.spce'
    elif args.water_model == 'tip4':
        args.water_model = 'leaprc.water.tip4pew'

    #do we have nonstandard residues?
    mol_data = {}
    nonstandard_res = {}
    standard_res = util.get_available_res(ff)
    ff = [ff]
    if nff:
        standard_res = standard_res.union(util.get_available_res(nff))
        ff = ff.append(nff)
    #pdb4amber seems to delete mercury (HG) along with hydrogens; for now my
    #hacky fix is to store the relevant atom info if mercury is present and add
    #the mercury back in after stripping...I'm preemptively doing this for
    #hafnium too
    metal_info = {}
   
    #if any structure was not provided in PDB format, we will attempt to create
    #one from what was provided using obabel, choosing a filename that will not
    #overwrite anything in the directory (optionally)
    for structure in args.structures:
        net_charge = None
        #the "structure" string in the args.structure list will be updated so
        #that it corresponds to the PDB we should use for subsequent steps
        assert os.path.isfile(structure),'%s does not exist\n' % structure
        #"base" is the base filename (no extension) from which others will be derived
        base = util.get_base(structure)
        ext = os.path.splitext(structure)[-1]
        if 'pdb' not in ext:
            #if it's a mol2, store the net_charge from the input because
            #conversion to a pdb and back to a mol2 with openbabel is not
            #guaranteed to result in the same partial charges
            if args.net_charge:
                net_charge = args.net_charge
            elif 'mol2' in ext:
                net_charge = util.get_charge(structure)
            outpdb = base + '.pdb'
            if not args.overwrite:
                outpdb = util.get_fname(outpdb)
            try:
                obabel[structure, '-O', outpdb, '-xn']()
            except Exception as e:
                print('Cannot create PDB from input, error {0} : {1}. Check \
{2}. Aborting...\n'.format(e.errno, e.strerror, outpdb))
                sys.exit()
            idx = args.structures.index(structure)
            args.structures[idx] = outpdb
            structure = outpdb
        mol_res = {}
        mol_data[structure] = pdb.simplepdb(structure)
        if not mol_data[structure].has_unique_names() and not mol_data[structure].is_protein():
            mol_data[structure].rename_atoms()
            print("Renaming atoms for",structure)
        mol_res[structure] = set(mol_data[structure].mol_data['resname'])
        ion_resnames = set(ions.keys())
        ions_present = set.intersection(mol_res[structure], ion_resnames)
        metal_info[structure] = {}
	#TODO: the pdb4amber problem this is trying to address still exists, but changes to this script's handling of ions mean this fails to solve the problem
        if 'HG' in ions_present:
            metal_info[structure]['HG'] = mol_data[structure].get_res_info({'resname':'HG'})
        if 'HF' in ions_present:
            metal_info[structure]['HF'] = mol_data[structure].get_res_info({'resname':'HF'})
        nonstandard_res[structure] = list(mol_res[structure] - standard_res -
                ion_resnames)

    #if nonstandard residues, do we have the necessary library files? 
    #check for prep, lib, and off; just add the frcmod if there is one
    libs = set([])
    # store the ligands that were parameterized by antechamber so they can be reparameterized by SMIRNOFF
    ante_lig = []
    if not args.libs:
        args.libs = []
    for struct,reslist in list(nonstandard_res.items()):
        #track which units you don't have libs for
        orphaned_res = set(reslist)
        base = util.get_base(struct)
        if orphaned_res:
        #try any user-provided locations first
            for user_lib in args.libs:
                for ext in ['.lib','.off','.prep']:
                    fname = user_lib + ext
                    if os.path.isfile(fname):
                        set_matches(fname, libs, reslist, orphaned_res,
                                mol_data[struct], True)
        #if residues are still undefined, check the current directory too
        if orphaned_res:
            local_libs = [name for name in glob.glob('*.lib') +
                    glob.glob('*.off') + glob.glob('*.prep')]
            for lib in local_libs:
                set_matches(lib, libs, reslist, orphaned_res, mol_data[struct])
   
        is_protein = mol_data[struct].is_protein()
        #for now, require that the ligand be provided separately from the protein -
        #that way we don't need to worry about differentiating between modified
        #residues (or other things we don't want to strip out of the protein) and
        #small molecules we can parametrize with antechamber
        assert not (is_protein and orphaned_res), \
        "Undefined units %s in protein - check for modified residues, ions, or \
cofactors\n" % ' '.join(orphaned_res)

        if is_protein and reslist:
            print("NOT RUNNING pdb4amber due to presence of modified residues.")
        elif is_protein and not args.noh: 
            fname = base + '_amber.pdb'
            command = pdb4amber['-y', '-i', struct, '-o', fname]
            runfile.writeln(command)
            code,stdout,stderr = command & TEE
            idx = args.structures.index(struct)
            args.structures[idx] = fname
            if not args.uninteractive:
                input('Read the above messages and then press any key to continue; note that prepareamber will insert any missing TER records but does not add ACE/NME caps...\n')
            #TODO: I mean, we _could_ add the ACE/NME caps...
            mol_data[fname] = pdb.simplepdb(fname)
            #if there were gaps, add appropriate TERs
            for line in stderr.splitlines():
                m = re.search(r'^gap .*between (\S+).(\d+)',line)
                if m:
                    gap_res = m.group(2)
                    mol_data[fname].add_ter(gap_res)
            for deleted_ion,data in metal_info[struct].items():
                current_residues = set(mol_data[fname].mol_data['resname'])
                if deleted_ion not in current_residues:
                    for res in data:
                        mol_data[fname].add_residue(res)

        assert len(orphaned_res)<2, "%s has multiple ligands; break them into \
separate files to process with antechamber\n" % struct

        #if we're handling a ligand and don't have library files, we will need at 
        #least the pdb-formatted data and a mol2 from which we can derive gasteiger 
        #charges for antechamber; make these with babel and find the net charge
        orphaned_res = list(orphaned_res)
        if orphaned_res:
            #"molname" will be the name of the unit for AMBER
            #TODO: check whether, if there are multiple ligands to be fit in
            #antechamber, the user has provided unit names for all of them or
            #they have distinct residue names 
            molname = orphaned_res[0]
            mol_data[struct].sanitize()
            mol_data[struct].set_resname(orphaned_res[0])
            tempname = base + '_temp.pdb'
            ligname = base + '_amber.pdb'
            if not args.overwrite:
                tempname = util.get_fname(tempname)
                ligname = util.get_fname(ligname)
            if os.path.isfile(ligname):
                os.remove(ligname)
            mol2 = base + '_amber.mol2'
            #create a PDB that has unique atom names, hydrogens, all HETATM
            #records, element names, and the correct residue name
            if not args.noh:
                mol_data[struct].writepdb(tempname)
                obabel[tempname, '-O', ligname, '-h','-xn']()
                os.remove(tempname)
                mol_data[struct] = pdb.simplepdb(ligname)
                mol_data[struct].sanitize()
                mol_data[struct].set_recordname('HETATM')
                os.remove(ligname)
                mol_data[struct].writepdb(ligname)
            else:
                mol_data[struct].writepdb(ligname)
            idx = args.structures.index(struct)
            args.structures[idx] = ligname
            mol_data[ligname] = mol_data[struct]
            #only compute if we didn't already get a value using an original mol2
            if net_charge is None:
                obabel[ligname, '-O', mol2]()
                net_charge = util.get_charge(mol2)
            #run antechamber
            print('Parametrizing unit %s with antechamber.\n' % ' '.join(orphaned_res))
            if util.is_secret_peptide(mol_data[ligname]):
                print('Warning: the ligand %s maybe actually be a peptide. If antechamber fails, check the residue names\n' %ligname)
            do_antechamber(ligname, net_charge, ff, molname, base)
            #add the libraries created in the last step to the libs list
            libs.add(base + '.lib')
            libs.add(base + '.frcmod')
            #Antechamber does not preserve the input atomnames. It comes packaged with
            #a program called match_atomname to cope with this. 
            command = match_atomname['-i', ligname, '-fi', 'pdb', '-r', mol2, '-fr',
                    'mol2', '-o', ligname, '-h', 1]
            runfile.writeln(command)
            command()
            mol_data[ligname] = pdb.simplepdb(ligname)
            ante_lig.append((ligname, mol2))

    #always add requeted frcmod files as they may apply to multiple ligands
    for lib in args.libs:
        if os.path.isfile(lib+'.frcmod'):
            libs.add(lib+'.frcmod')
            
    #ok, now we can be pretty sure we know what to do and that we are able to do it
    #create complex if there are multiple structures
    if args.out_name:
       complex_name = args.out_name + '.pdb'
    elif len(args.structures) > 1:
       complex_name = 'complex.pdb'
    else:
       complex_name = args.structures[0]

    runfile = Tee(os.path.splitext(complex_name)[0]+'.run')
    
    if len(args.structures) > 1:
        start_atom, start_res = 1,1
        if os.path.isfile(complex_name):
            os.remove(complex_name)
        final_mols = [mol_data[file] for file in args.structures]
        mol_data[args.structures[0]].writepdb(complex_name, final_mols)
    elif complex_name != args.structures[0]:
        mol_data[args.structures[0]].writepdb(complex_name)    
        
    base = util.get_base(complex_name)
    #make initial parameters files
    make_amber_parm(complex_name, base, ff, 'complex', args.water_model, args.water_dist, libs, extra=args.extra)
    if (len(ante_lig) > 0 and not args.no_openff):
        # Only do reparm if there are ligs to reparm    
        reparm(ante_lig, base)
    #run the two minimization and two pre-production  MDs
    if not args.parm_only: 
        do_amber_preproduction(complex_name, base, args, ff)
        #run the final production MD
        do_amber_production(base, args.run_prod_md)
