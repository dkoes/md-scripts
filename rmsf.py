#!/usr/bin/env python3

'''Calculate RMSF of specified trajectory or trajectories (all with same topology).
   Outputs csv file and annotated structure.
'''
import MDAnalysis as mda
import argparse
from MDAnalysis.analysis import rms, align
import MDAnalysis.analysis
import os
import pandas as pd
from os.path import splitext, basename

parser = argparse.ArgumentParser(description='Identify most distinct frames of md trajectory\nIMPORTANT: assumes an aligned input')
parser.add_argument("topology",required=True,help='Topology file')
parser.add_argument("trajectory", nargs='+',required=True,'Trajectory file(s)')
parser.add_argument("--selection",default="backbone",required=False,'Selection to align to, default backbone')
parser.add_argument("--out",default=None,required=False,help='Output prefix for csv and pdb file')

args = parser.parse_args()

U = mda.Universe(args.topology, args.trajectory)

#strip water
notwat = U.select_atoms('not resname WAT')
u = mda.Merge(notwat).load_new(
         AnalysisFromFunction(lambda ag: ag.positions.copy(), notwat).run().results['timeseries'],
         format=MemoryReader)

#find average structure
average = align.AverageStructure(u, u, select=args.selection, ref_frame=0).run()
ref = average.results.universe

# align to this structure
aligner = align.AlignTraj(u, ref, select=args.selection, in_memory=True).run()

sel = u.select_atoms(args.selection)
R = rms.RMSF(sel).run()

if args.out == None:
    prefix = splitext(args.topology)[0]+'.rmsf'
else:
    prefix = args.out
    
res = zip(sel.indices,sel.resids,sel.resnames,R.results['rmsf'])
df = pd.DataFrame(res,columns=('index','resid','resname','rmsf'))
df.to_csv(f'{prefix}.csv')

u.add_TopologyAttr('tempfactors') # add empty attribute for all atoms
u.trajectory[0]
for residue, r_value in zip(sel.residues, R.results.rmsf):
    residue.atoms.tempfactors = r_value
    
u.atoms.write(f'{prefix}.pdb')
