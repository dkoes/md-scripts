#!/usr/bin/env python3

import sys, MDAnalysis
import numpy as np
from os.path import splitext
from MDAnalysis.analysis.rms import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import argparse
import multiprocessing
from functools import partial

parser = argparse.ArgumentParser(description='Generate pairwise RMSD heatmap plot.\nIMPORTANT: assumes an aligned input')
parser.add_argument("topology")
parser.add_argument("trajectory")
parser.add_argument("--selection",default="backbone",help="MDAnalysis selection for computing RMSD",required=False)
parser.add_argument("--step",default=10, type=int, required=False,help="Frames to skip over")
parser.add_argument('--title',help="Graph title")
parser.add_argument('-o','--output',type=str,help="Output filename")
parser.add_argument('--max',type=float,help='Max RMSD value to consider',default=None)
args = parser.parse_args()

top = args.topology
traj = args.trajectory

base = splitext(top)[0]
if not args.title:
    args.title = base
    
if not args.output:
    args.output = base+'.png'
    
u1 = MDAnalysis.Universe(top,traj)
u2 = MDAnalysis.Universe(top,traj)

# arguments: topology, trajector, [selection], [graph title], [step]
# todo, switch to argparse
sel = args.selection

n = u1.trajectory.n_frames

div = args.step
downn = n//div 
if n % div != 0:
    downn += 1
rmat = np.zeros((downn,downn))

sel1 = u1.select_atoms(sel)
sel2 = u2.select_atoms(sel)

def rmsd_row(frame_index, u1, u2, sel1, sel2):
    row = []
    t1 = u1.trajectory[frame_index]
    for t2 in u2.trajectory[t1.frame::div]:
        row.append((t1.frame//div, t2.frame//div, rmsd(sel1.positions,sel2.positions)))
    return row
    
run_per_frame = partial(rmsd_row, u1=u1, u2=u2, sel1=sel1, sel2=sel2)

print('Processing')
pool = multiprocessing.Pool()
rows = pool.map(run_per_frame, range(0,len(u1.trajectory),div))
pool.close()
print('Done')

for row in rows:
    for i,j,r  in row:
        rmat[i,j] = rmat[j,i] = r
    

np.set_printoptions(threshold=np.inf,precision=2)


#find frame with most other frames under cutoff
cutoff = 2.0
cnts = np.array([len(row[row < cutoff]) for row in rmat])
pos = cnts.argmax()
print("Frame %d is within %.2f of %d frames" % (pos, cutoff, cnts[pos]))

import matplotlib.pylab as plt

plt.figure()
plt.title(args.title)

#multiples of 10, but no more than 6ish ticks
n = 10
while len(rmat)/n > 6:
    n += 10

if args.max:
    sns.heatmap(rmat,square=True,xticklabels=n,yticklabels=n,cmap='YlGnBu',cbar_kws={'label':'RMSD'},vmin=0,vmax=args.max)
else:
    sns.heatmap(rmat,square=True,xticklabels=n,yticklabels=n,cmap='YlGnBu',cbar_kws={'label':'RMSD'},vmin=0)
plt.xlabel("Frame #")
plt.ylabel("Frame #")
ax = plt.gca()
ax.tick_params(direction='out')
plt.tight_layout()

plt.savefig(args.output)
