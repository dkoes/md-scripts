#!/usr/bin/env python3

import sys, MDAnalysis
import numpy as np
from os.path import splitext
from MDAnalysis.analysis.rms import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import argparse

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
rmat = np.zeros((n//div,n//div))

sel1 = u1.select_atoms(sel)
sel2 = u2.select_atoms(sel)
#print len(sel1),len(sel2)
for t1 in u1.trajectory[::div]:
    for t2 in u2.trajectory[::div]:
        if t1.frame < t2.frame:
            rmat[t2.frame//div, t1.frame//div] = rmat[t1.frame//div, t2.frame//div] = rmsd(sel1.positions,sel2.positions)

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
