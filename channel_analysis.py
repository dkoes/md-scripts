#!/usr/bin/env python3

import numpy as np 
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader
import MDAnalysis.transformations as trans
from mdahole2.analysis import HoleAnalysis
from MDAnalysis.coordinates.memory import MemoryReader
from matplotlib.ticker import FuncFormatter

import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
import warnings
import re,glob,os,gzip,pickle,argparse,sys

default_bottom_selection = 'resname ASP and resid 238 and name CA'
default_top_selection = 'resname ALA and resid 263 and name CA'
default_ion_selection = 'resname CLA'
default_res_offset=44

def evaluate_box(U,bottom_selection = default_bottom_selection,
                 top_selection = default_top_selection):
    '''Return true if the box appears to be properly wrapped with the TMD 
    roughly in the center'''
    low = U.atoms.positions[:,2].min()
    high = U.atoms.positions[:,2].max()
    zdiff = high-low 
    if zdiff > 220:
        print(f'Unwrapped or unexpectedly large box. zdiff: {zdiff:.2f}')
        print('''Use cpptraj to wrap the trajectory with the membrane at the center:
parm md.psf
trajin md.nc
autoimage anchor :POPC
rms first @CA
strip :WAT outprefix strip
trajout strip.dcd dcd
''')
        return False
    top = U.select_atoms(top_selection).positions[:,2].mean()
    bottom = U.select_atoms(bottom_selection).positions[:,2].mean()
    mid = (top+bottom)/2
    frame = U.trajectory[0]
    print(f"TMD midpoint: {mid:.2f} Box midpoint: {zdiff/2+low:.2f}")
    if bottom-low < 30 or top > high-30:
       print("Not enough space at top or bottom")
       return False
    return True

def compute_ion_transitions(u,ion_selection = 'resname CLA',
                            bottom_selection = 'resname ASP and resid 238 and name CA',
                            top_selection = 'resname ALA and resid 263 and name CA'):
    '''Identify all the _full_ transitions through the transmembrane domain (as defined
    by the provided protein selections) for the specified ion.
    Returns the per-frame cummulative transition counts for both forward and backwards transitions.'''
    evaluate_box(u,bottom_selection,top_selection)

    top = u.select_atoms(top_selection)
    print(f"Selection defining top of TMD has {len(top)} atoms")
    bottom = u.select_atoms(bottom_selection)
    print(f"Selection defining bottom of TMD has {len(bottom)} atoms")

    cl = u.select_atoms(ion_selection)
    print(f'Number of ions ({ion_selection}): {len(cl)}')

    tmd_t = top.positions[:,2].mean()
    tmd_b = bottom.positions[:,2].mean()

    print(f'Top {tmd_t}')
    print(f'Bottom {tmd_b}')

    locs = np.array([(cl.positions[:,2] > tmd_t) + \
            ((cl.positions[:,2] <= tmd_t) & (cl.positions[:,2] > tmd_b))*2+ \
            (cl.positions[:,2] <= tmd_b)*4 for _ in u.trajectory]).T

    prevAB = np.zeros(locs.shape[0])
    prevBC = np.zeros(locs.shape[0])

    ABC = np.zeros(locs.shape[1])
    CBA = np.zeros(locs.shape[1])

    for i in range(1,locs.shape[1]):
        diff = locs[:,i]-locs[:,i-1]
        ABC[i] = ((diff == 2) & (prevAB == 1)).sum()
        prevAB[diff != 0] = 0
        prevAB[diff == 1] = 1    

        CBA[i] = ((diff == -1) & (prevBC == 1)).sum()
        prevBC[diff != 0] = 0
        prevBC[diff == -2] = 1
    
    return ABC.cumsum(),CBA.cumsum()

def plot_ion_density(u,outfile=None,ion_selection=default_ion_selection,
                     bottom_selection=default_bottom_selection,
                     top_selection=default_top_selection,
                     steps_per_ns=10):
    fig = plt.figure()
    cl = u.select_atoms(ion_selection)
    top = u.select_atoms(top_selection)
    bottom = u.select_atoms(bottom_selection)
    N = len(u.trajectory)
    tmd_t = [top.positions[:,2].mean() for _ in u.trajectory]
    tmd_b = [bottom.positions[:,2].mean() for _ in u.trajectory]
    positions = np.array([cl.positions[:,2] for _ in u.trajectory]).T
    
    for p in positions:
        plt.scatter(range(N),p,1,marker='.',alpha=0.5)
    plt.plot(tmd_t,color='k')
    plt.plot(tmd_b,color='k')
    plt.xlim(0,N); plt.xlabel('Time (ns)',fontsize=14); 
    plt.ylabel('Ion Z-coord',fontsize=14);
    
    formatter = FuncFormatter(lambda x,pos:  f'{x/steps_per_ns:g}')
    plt.gca().xaxis.set_major_formatter(formatter)
    
    if(outfile):
        plt.savefig(outfile,dpi=300,bbox_inches='tight')

def plot_transitions(u,minrads=None,outfile=None,ion_name='Cl',
                     ion_selection=default_ion_selection,
                     bottom_selection=default_bottom_selection,
                     top_selection=default_top_selection,
                     steps_per_ns=10):
    fig, ax1 = plt.subplots(figsize=(6, 4))
    ax2 = ax1.twinx()
    if minrads:
        ax1.axhline(y=2.3, color="k", linestyle="--")
        ax1.axhline(y=1.75, color="red", linestyle="--")
        ax1.plot(minrads,color='#1f77b4',zorder=10,linewidth=1)
        ax1.set_ylabel("Minimum Radius ($\mathrm{\AA}$)", color="#1f77b4", fontsize=14)
        ax1.set_ylim(0, 4)
        ax1.tick_params(axis="y", labelcolor="#1f77b4")

    ax1.set_xlabel('Time (ns)',fontsize=14);    

    ABC,CBA = compute_ion_transitions(u,ion_selection=ion_selection,bottom_selection=bottom_selection,top_selection=top_selection)
    plt.xlim(0,len(ABC))
    ax2.plot(ABC,color='#2ca02c',zorder=100,linewidth=3)
    ax2.plot(CBA,color='pink',zorder=99,linewidth=1)

    ax2.set_zorder(2)
    ax2.set_ylabel(f'Cummulative {ion_name} Transitions',fontsize=14,color='#2ca02c')
    ax2.tick_params(axis="y", labelcolor="#2ca02c")
    ax2.set_ylim(0, 60)

    formatter = FuncFormatter(lambda x,pos:  f'{x/steps_per_ns:g}')
    ax1.xaxis.set_major_formatter(formatter)
    
    if outfile:
        plt.savefig(outfile,dpi=300,bbox_inches='tight')

def process_hole_analysis(u,ha,bottom_selection=default_bottom_selection,
                        top_selection=default_top_selection,
                        res_offset=44):
    '''pull out radius profile between specified selections.  
    res_offset is a number to add to resids to get canonical numbering in labels
    Returns a dictionary with
    label_coords: z-coordinate of labels
    labels: residue labels (Ca)
    radii_zpos: position of radius bins
    radii: list of radii at each bin to compute summary statistics with
    min_radii: array indexed by frame of minimum radius in specified region at that time point
    '''
    m = re.search(r'resid (\d+)',top_selection)
    if not m:
        print(f'Missing resid in top_selection: {top_selection}')
    top_res = int(m.group(1))
    
    m = re.search(r'resid (\d+)',bottom_selection)
    if not m:
        print(f'Missing resid in bottom_selection: {bottom_selection}')
    bot_res = int(m.group(1))
    
    topz = u.select_atoms(f'protein and resid {top_res} and name CA').positions.mean(axis=0)[2]
    botz = u.select_atoms(f'protein and resid {bot_res} and name CA').positions.mean(axis=0)[2]
    
    #get coordinates of each residue
    zcoords = []
    labels = []
    for resi in range(bot_res,top_res+1):
        sel = u.select_atoms(f'protein and resid {resi} and name CA')
        z = sel.positions.mean(axis=0)[2]
        name = sel.resnames[0]
        zcoords.append(z)
        labels.append(f'{name}{resi+res_offset}')
    
    radii, edges = ha.bin_radii(range=(min(zcoords),max(zcoords)))
    # middle of each bin
    zpos = (edges+(edges[1]-edges[0])/2)[:-1]
    
    minrads = [ha.profiles[frame].radius[(botz < ha.profiles[frame].rxn_coord) & (ha.profiles[frame].rxn_coord < topz)].min() for frame in sorted(ha.profiles.keys())]

    return {'label_coords': zcoords,
            'labels': labels,
            'radii_zpos': zpos,
            'radii': radii,
            'min_radii': minrads}

def unique_name(u):
    '''Create a unique name for the universe based on its input files'''
    def sanitize(s):
        #return a more palitable string for a filename from a path
        return s.replace('.','/').lstrip('/').replace('/','_')
    topname = sanitize(u.filename)
    if hasattr(u.trajectory,'filenames'):
        trajname = '__'.join(map(sanitize,u.trajectory.filenames))
    else:
        trajname = sanitize(u.trajectory.filename)
    outprefix = f'{topname}-{trajname}'
    
    return outprefix
        
def align_traj(u, uref, bottom_selection=default_bottom_selection, 
                          top_selection=default_top_selection,inmem=True):
    '''Align trajectory to top/bottom.'''
                #make sure trajectory is aligned
    uref.trajectory[0]
    print("Aligning trajectory")
    selstr = f'({top_selection}) or ({bottom_selection})'
    sel = u.select_atoms(selstr)
    
    if len(sel) == 0:
        sys.stderr.write(f'Selection "{selstr}" is empty.')
        sys.exit(1)
        
    if inmem or isinstance(u.trajectory, MemoryReader):
        align.AlignTraj(u, uref, selstr, in_memory=True).run()
    else:
        align.AlignTraj(u, uref, selstr, filename=f'{outprefix}_aligned.dcd', in_memory=False).run()        
        u = mda.Universe(u.filename ,f'{outprefix}_aligned.dcd')
    
    return u
        
def perform_hole_analysis(u, outprefix=None, 
                          bottom_selection=default_bottom_selection, 
                          top_selection=default_top_selection,
                          res_offset=default_res_offset,
                          hole_exe='hole'):
    '''Run the hole analysis.  Will create a pickle file that it will reload if already present.
       Assumes aligned trajecotry
       Will also create a vmd script for visualization.
       u - MDAnalysis universe
       outprefix - prefix to use for created files
       bottom_selection - bottom of tmd
       top_selection - top of tmd
       res_offset - number to add to resid to get desired labels
    '''
    if outprefix == None:
        outprefix = unique_name(u)

    print(f"Using output prefix {outprefix}")
    try:
        ha = pickle.load(gzip.open(f'{outprefix}.pkl.gz'))
        print(f"Loaded cached result from {outprefix}.pkl.gz")
        return ha
    except:
        print("No cached result found.")
    

    selstr = f'({top_selection}) or ({bottom_selection})'
    sel = u.select_atoms(selstr)
        
    print("Performing hole analysis (this will take a long time)")
    pid = os.getpid() #can't use outprefix because hole will truncate long file names
    print(f"Temporary files will be created in current directory with prefix hole{pid}")
    center = (sel).positions.mean(axis=0)
    ha = HoleAnalysis(
        u, executable=hole_exe, 
        cpoint=center, 
        prefix=f'hole{pid}',
        cvect=[0, 0, 1],
        ignore_residues=['WAT','CLA','POT']
    ) 
    ha.run()
    ha.create_vmd_surface(filename=f"{outprefix}.vmd", dot_density=15)
    tmd = process_hole_analysis(u,ha,bottom_selection, top_selection, res_offset)
    ha.tmd = tmd
    with gzip.open(f'{outprefix}.pkl.gz','wb') as out:
        pickle.dump(ha,out)
    print(f'Created cached output {outprefix}.pkl.gz')
            
    ha.delete_temporary_files()
    return ha

def plot_profile(ha,outfile=None):
    '''Plot a profile of a single hole.'''
    radii = ha.tmd['radii']
    zpos = ha.tmd['radii_zpos']
    zcoords = ha.tmd['label_coords']
    labels = ha.tmd['labels']
    means = np.array([np.mean(r) for r in radii])
    stds = np.array([np.std(r) for r in radii])
    mins = np.array([np.min(r) for r in radii])
    maxs = np.array([np.max(r) for r in radii])
    los = np.array([np.percentile(r,5) for r in radii])
    his = np.array([np.percentile(r,95) for r in radii])
        
    fig = plt.figure(figsize=(4,6))
    plt.axvline(x=2.3, color="k", linestyle="--",linewidth=1)
    plt.axvline(x=1.75, color="red", linestyle="--",linewidth=1)    
    plt.plot(means,zpos,zorder=15)
    plt.fill_betweenx(
        zpos,
        los,
        his,
        facecolor="lightblue",
        alpha=0.2,zorder=10
    )
    plt.fill_betweenx(
        zpos,
        means - stds,
        means + stds,
        facecolor="lightblue",
        alpha=0.3,zorder=11
    )
    #plt.plot(mins,zpos,color='lightblue',linewidth=0.5)
    #plt.plot(maxs,zpos,color='lightblue',linewidth=0.5)

    plt.yticks(zcoords[1::2],labels[1::2]);
    plt.ylim(min(zpos),max(zpos))
    plt.xlim(0,10)
    plt.xlabel('Pore Radius ($\mathrm{\AA}$)',fontsize=14);
    
    if outfile:
        plt.savefig(outfile,dpi=300,bbox_inches='tight')
        
    return fig

def combine_tmd(hfiles):
    '''Merge the tmd attribute of all the provided ha.pkl.gz.
    radii are merged while the first ha file is used for the rest.'''
    from types import SimpleNamespace
    C = SimpleNamespace(tmd={})    
    C.tmd['allradii'] = []
    C.tmd['allradii_zpos'] = []
    C.tmd['alllabel_coords'] = []
    C.tmd['alllabels'] = []
    for hfile in hfiles:
        ha = pickle.load(gzip.open(hfile))
        C.tmd['allradii'].append(ha.tmd['radii'])        
        C.tmd['allradii_zpos'].append(ha.tmd['radii_zpos'])
        C.tmd['alllabel_coords'].append(ha.tmd['label_coords'])
        C.tmd['alllabels'].append(ha.tmd['labels'])
        
    C.tmd['radii'] = [np.hstack(r) for r in zip(*C.tmd['allradii'])]   
    C.tmd['radii_zpos'] = C.tmd['allradii_zpos'][0].copy()
    C.tmd['label_coords'] = C.tmd['alllabel_coords'][0].copy()
    C.tmd['labels'] = C.tmd['alllabels'][0].copy()
    return C
    
def plot_combined_tmd(C,outfile=None):
    '''Given a list of ha.pkl.gz files, plot the pore profiles.'''

    means = np.array([np.mean(r) for r in C.tmd['radii']])
    stds = np.array([np.std(r) for r in C.tmd['radii']])
    mins = np.array([np.min(r) for r in C.tmd['radii']])
    maxs = np.array([np.max(r) for r in C.tmd['radii']])
    los = np.array([np.percentile(r,5) for r in C.tmd['radii']])
    his = np.array([np.percentile(r,95) for r in C.tmd['radii']])

    zpos = C.tmd['radii_zpos']
    zcoords = C.tmd['label_coords']
    zlabels = C.tmd['labels']

    fig = plt.figure(figsize=(4,6))
    plt.axvline(x=2.3, color="k", linestyle="--",linewidth=1)
    plt.axvline(x=1.75, color="red", linestyle="--",linewidth=1)    
    plt.plot(means,zpos,zorder=15,color='#1f77b4',linewidth=3)
    plt.fill_betweenx(
        zpos,
        los,
        his,
        facecolor="lightblue",
        alpha=0.2,zorder=10
    )
    plt.fill_betweenx(
        zpos,
        means - stds,
        means + stds,
        facecolor="lightblue",
        alpha=0.3,zorder=11
    )

    for R in C.tmd['allradii']:
        m = np.array([np.mean(r) for r in R])
        plt.plot(m,zpos,zorder=5,linewidth=1,color='#1f77b4')

    plt.yticks(zcoords[1::2],zlabels[1::2]);
    plt.ylim(min(zpos),max(zpos))
    plt.xlim(0,10)
    plt.xlabel('Pore Radius ($\mathrm{\AA}$)',fontsize=14)
    
    if outfile:
        plt.savefig(outfile,dpi=300,bbox_inches='tight')
    return fig
    
    
def plot_separate_combined_tmd(C,outfile=None,labels=None):
    '''Given a list of ha.pkl.gz files, plot the pore profiles.'''

    zpos = C.tmd['radii_zpos']
    zcoords = C.tmd['label_coords']
    zlabels = C.tmd['labels']

    fig = plt.figure(figsize=(4,6))
    plt.axvline(x=2.3, color="k", linestyle="--",linewidth=1)
    plt.axvline(x=1.75, color="red", linestyle="--",linewidth=1)  

    for i,R in enumerate(C.tmd['allradii']):
        m = np.array([np.mean(r) for r in R])
        s = np.array([np.std(r) for r in R])
        line = plt.plot(m,zpos,zorder=15,label=labels[i] if labels else None)[0]
        c = line.get_color()
        plt.fill_betweenx(
            zpos,
            m - s,
            m + s,
            color=c,
            alpha=0.1,zorder=11
        )


    plt.yticks(zcoords[1::2],zlabels[1::2]);
    plt.ylim(min(zpos),max(zpos))
    plt.xlim(0,10)
    plt.xlabel('Pore Radius ($\mathrm{\AA}$)',fontsize=14)
    if labels:
        plt.legend()
    
    if outfile:
        plt.savefig(outfile,dpi=300,bbox_inches='tight')
    return fig        
    
def plot_combined_transitions(us,outfile=None,title=None,ion_selection='resname CLA',ion_name='Cl-',
                        bottom_selection=default_bottom_selection,
                        top_selection=default_top_selection,
                        steps_per_ns=10):
    '''Plot all ion transitions for list of universes provided.
       Trajectories should be aligned.'''
        
    ABCs = []
    CBAs = []
    radii = []
    for u in us:
        ABC,CBA = compute_ion_transitions(u,ion_selection=ion_selection,bottom_selection=bottom_selection,top_selection=top_selection)
        ABCs.append(ABC)
        CBAs.append(CBA)
    
    fig = plt.figure(figsize=(6,4))
    plt.xlim(0,len(ABCs[0]))
    plt.ylabel(f'Cummulative {ion_name} Transitions',fontsize=14)
    plt.ylim(0, 60)
    for ABC,CBA in zip(ABCs,CBAs):
        plt.plot(ABC,color='#2ca02c',alpha=0.5,zorder=10)
        plt.plot(CBA,color='pink',linewidth=1,alpha=0.1)
    
    plt.plot(np.array(ABCs).mean(axis=0),color='#2ca02c',zorder=10,linewidth=3)
    plt.plot(np.array(CBAs).mean(axis=0),color='pink',zorder=10,linewidth=3)
    if title:
        plt.title(title)
        
    formatter = FuncFormatter(lambda x,pos:  f'{x/steps_per_ns:g}')
    plt.gca().xaxis.set_major_formatter(formatter)
    plt.xlabel("Time (ns)",fontsize=14)
    
    tot = np.array(ABCs).mean(axis=0)[-1]
    ns = len(ABCs[0])/steps_per_ns
    rate = tot/ns
    plt.text(0.5,0.95, f'{rate:.2f} transitions/ns', horizontalalignment='center',verticalalignment='top',fontsize=14,transform=plt.gca().transAxes)
        
    if outfile != None:
        plt.savefig(outfile,dpi=300,bbox_inches='tight')
    return fig    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate analysis of membrane pore simulations.\nIMPORTANT: assumes a wrapped and aligned input')
    parser.add_argument("topology",help="Topology file")
    parser.add_argument("trajectory",nargs='+',help="Trajectory file(s)")
    parser.add_argument("--top_selection",default=default_top_selection,help="Selection string for top of hole.")
    parser.add_argument("--bottom_selection",default=default_bottom_selection,help="Selection string for bottom of hole.")
    parser.add_argument("--res_offset",default=44, type=int, help="Residue numbering adjustment")
    parser.add_argument("--hole_exe",default="hole",help="HOLE executable if not in path")
    parser.add_argument('--prefix',default=None,help="Prefix for output files")
    parser.add_argument('--steps_per_ns',default=10,help="Number of frames per a ns")
    parser.add_argument('--title',default=None,help="Plot title")

    args = parser.parse_args()

    #process trajectories individually
    for traj in args.trajectory:  #mdanalysis provides cryptic errors when traj doesn't exist so check ourselves    
        if not os.path.exists(traj):
            sys.stderr.write(f'{traj} does not exist\n')
            exit(-1)
               
    uref = mda.Universe(args.topology, args.trajectory[0])
    hfiles = []
    us = []
    for i,traj in enumerate(args.trajectory):
        u = mda.Universe(args.topology, traj)
        if len(u.select_atoms(args.bottom_selection)) == 0:
            print("Bottom selection did not select atoms:",args.bottom_selection)
            sys.exit(-1)
        if len(u.select_atoms(args.top_selection)) == 0:
            print("Top selection did not select atoms:",args.top_selection)
            sys.exit(-1)            
            
        u = align_traj(u, uref,bottom_selection=args.bottom_selection,
                     top_selection=args.top_selection)
        us.append(u)
        if args.prefix == None:
            prefix = unique_name(u)
        else:
            prefix = f'{args.prefix}_{i}'
        plot_ion_density(u,f'{prefix}_CLA_density.png', 'resname CLA',
                            bottom_selection=args.bottom_selection, 
                            top_selection=args.top_selection,
                            steps_per_ns=args.steps_per_ns)
        plot_ion_density(u,f'{prefix}_POT_density.png', 'resname POT',
                            bottom_selection=args.bottom_selection, 
                            top_selection=args.top_selection,
                            steps_per_ns=args.steps_per_ns)
        ha = perform_hole_analysis(u, outprefix=f'{prefix}_ha', 
                          bottom_selection=args.bottom_selection, 
                          top_selection=args.top_selection,
                          res_offset=args.res_offset,
                          hole_exe=args.hole_exe)
        hfiles.append(f'{prefix}_ha.pkl.gz')
        plot_profile(ha,f'{prefix}_profile.pdf')
        plot_transitions(u,minrads=ha.tmd['min_radii'],
                     outfile=f'{prefix}_CLA_transitions.pdf',
                     ion_name='Cl-',
                     ion_selection='resname CLA',
                     bottom_selection=args.bottom_selection,
                     top_selection=args.top_selection,
                     steps_per_ns=args.steps_per_ns)
        plot_transitions(u,minrads=ha.tmd['min_radii'],
                     outfile=f'{prefix}_POT_transitions.pdf',
                     ion_name='K+',
                     ion_selection='resname POT',
                     bottom_selection=args.bottom_selection,
                     top_selection=args.top_selection,
                     steps_per_ns=args.steps_per_ns)

    C = combine_tmd(hfiles)
    if args.prefix == None:
        prefix = unique_name(uref)
    else:
        prefix = args.prefix  
        
    if args.title == None:
        title = prefix
    else:
        title = args.title
    plot_combined_tmd(C,f'{prefix}_combined.pdf')
    plot_separate_combined_tmd(C,f'{prefix}_sep_combined.pdf')
    
    plot_combined_transitions(us,f'{prefix}_CLA_transitions_combined.pdf',title,
                ion_selection='resname CLA',ion_name='Cl-',
                bottom_selection=args.bottom_selection,
                top_selection=args.top_selection,
                steps_per_ns=args.steps_per_ns)
                
    plot_combined_transitions(us,f'{prefix}_POT_transitions_combined.pdf',title,
                ion_selection='resname POT',ion_name='K+',
                bottom_selection=args.bottom_selection,
                top_selection=args.top_selection,
                steps_per_ns=args.steps_per_ns)    
    
