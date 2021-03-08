#!/usr/local/bin/python

import MDAnalysis
from MDAnalysis import *
from MDAnalysis.analysis.align import *
import sys

#align a trajectory to the backbone of a reference structure
# <topology> <trajectory> <reference> <output traj>

if len(sys.argv) != 5:
	print("Need  <topology> <trajectory> <reference> <output traj>")
	sys.exit(-1)
	
trj = Universe(sys.argv[1],sys.argv[2])

ref = Universe(sys.argv[1],sys.argv[3])
out = sys.argv[4];


rms_fit_trj(trj, ref, "backbone", out, "%s.rmsd" % out)
