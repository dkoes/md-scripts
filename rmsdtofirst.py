#!/usr/local/bin/python

import MDAnalysis, sys, argparse, numpy, math
import MDAnalysis.core.qcprot as qcp
from MDAnalysis.analysis.align import *
#take a topology, trajectory, and selction and output rmsd to the first frame
#throughout the trajectory
parser = argparse.ArgumentParser(description='Calculate RMSD to first frame (optionally perform alignment)')
parser.add_argument('topo',metavar='topology file')
parser.add_argument('traj',metavar='trajectory file')
parser.add_argument('sel',metavar='selection')
parser.add_argument('--minimize','-m', default=False, action='store_true', help="compute minimal rmsd for each frame")
args = parser.parse_args()

#read topo and trajectory
u = MDAnalysis.Universe(args.topo, args.traj)
refatoms = u.selectAtoms(args.sel);
trajatoms = u.selectAtoms(args.sel);

ref_com = refatoms.centerOfMass().astype(numpy.float32)

if args.minimize:
    ref_coordinates = refatoms.coordinates() - ref_com
else:
    ref_coordinates = refatoms.coordinates()
# allocate the array for selection atom coords
traj_coordinates = trajatoms.coordinates().copy()
natoms = trajatoms.numberOfAtoms()

# RMSD timeseries
frames = u.trajectory
nframes = len(frames)

# R: rotation matrix that aligns r-r_com, x~-x~com
#    (x~: selected coordinates, x: all coordinates)
# Final transformed traj coordinates: x' = (x-x~_com)*R + ref_com
rot = numpy.zeros(9,dtype=numpy.float64)      # allocate space for calculation
R = numpy.matrix(rot.reshape(3,3))

for k,ts in enumerate(frames):
    # shift coordinates for rotation fitting
    # selection is updated with the time frame

    if args.minimize:
        x_com = trajatoms.centerOfMass().astype(numpy.float32)
        traj_coordinates[:] = trajatoms.coordinates() - x_com
        r = rmsd(ref_coordinates, traj_coordinates)
    else:
        traj_coordinates[:] = trajatoms.coordinates() 
        diff = ref_coordinates - traj_coordinates
        r = math.sqrt(numpy.square(diff.flatten()).sum()/natoms)
    
    print k,r
    


