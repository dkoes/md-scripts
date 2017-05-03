#!/usr/bin/python

import sys, MDAnalysis
import numpy as np
from os.path import splitext
from MDAnalysis.analysis.rms import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
top = sys.argv[1]
traj = sys.argv[2]
u1 = MDAnalysis.Universe(top,traj)
u2 = MDAnalysis.Universe(top,traj)

# arguments: topology, trajector, [selection], [graph title], [step]
sel = "backbone"
if len(sys.argv) > 3:
    sel = sys.argv[3]

n = u1.trajectory.n_frames

if len(sys.argv) > 5:
    div = int(sys.argv[5])
else:
    div = 10
rmat = np.zeros((n/div,n/div))

sel1 = u1.select_atoms(sel)
sel2 = u2.select_atoms(sel)
for t1 in u1.trajectory[::div]:
    for t2 in u2.trajectory[::div]:
        if t1.frame < t2.frame:
            rmat[t2.frame/div, t1.frame/div] = rmat[t1.frame/div, t2.frame/div] = rmsd(sel1.coordinates(),sel2.coordinates())

np.set_printoptions(threshold=np.inf,precision=2)


#find frame with most other frames under cutoff
cutoff = 2.0
cnts = np.array([len(row[row < cutoff]) for row in rmat])
pos = cnts.argmax()
print "Frame %d is within %.2f of %d frames" % (pos, cutoff, cnts[pos])

import matplotlib.pylab as plt

plt.figure()

if len(sys.argv) > 4:
    name = sys.argv[4]
else:
    name = splitext(sys.argv[1])[0] 
plt.title(name)
plt.xlabel("Frame #")
plt.ylabel("Frame #")
ax = plt.gca()
im = ax.imshow(rmat)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im,cax=cax)
cbar.set_label('Backbone RMSD',rotation=-90,labelpad=15)

plt.savefig(name + '.png')
