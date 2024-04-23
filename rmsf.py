#!/usr/bin/env python3

'''Calculate RMSF of specified trajectory or trajectories (all with same topology).
   Outputs csv file and annotated structure.
'''
import MDAnalysis as mda
import argparse
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader
import MDAnalysis.analysis
import os
import pandas as pd
from os.path import splitext, basename
from MDAnalysis.coordinates.PDB import MultiPDBWriter

parser = argparse.ArgumentParser(description='Calculate RMSF')
parser.add_argument("topology",help='Topology file')
parser.add_argument("trajectory", nargs='+',help='Trajectory file(s)')
parser.add_argument("--selection",default="backbone",required=False,help='Selection to align to, default backbone')
parser.add_argument("--out",default=None,required=False,help='Output prefix for csv and pdb file')
args = parser.parse_args()

U = mda.Universe(args.topology, args.trajectory)

#strip water
notwat = U.select_atoms('not (resname WAT Cl- Na+)')
u = mda.Merge(notwat).load_new(
         AnalysisFromFunction(lambda ag: ag.positions.copy(), notwat).run().results['timeseries'],
         format=MemoryReader)

#find average structure
average = align.AverageStructure(u, u, select=args.selection, ref_frame=0).run()
ref = average.results.universe

# align to this structure
aligner = align.AlignTraj(u, ref, select=args.selection, in_memory=True).run()

sel = u.select_atoms('protein')
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

sel.tempfactors = R.results.rmsf

u.atoms.write(f'{prefix}.pdb')

