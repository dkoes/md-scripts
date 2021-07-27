#!/usr/bin/env python

'''Split a trajectory into N equal pieces'''

import MDAnalysis, argparse
import os
from os.path import splitext, basename


parser = argparse.ArgumentParser(description='Split a trajectory into N equal pieces')
parser.add_argument("topology")
parser.add_argument("trajectory")
parser.add_argument("-n",type=int,default=10)

args = parser.parse_args()

format = None
if args.trajectory.endswith('.nc'):
    format = "ncdf" #not sure why .nc isn't auto recognized any more
model = MDAnalysis.Universe(args.topology,args.trajectory,format=format)

basename,ext = splitext(args.trajectory)
if ext == '.nc':
    ext = '.ncdf'
nframes = model.trajectory.n_frames
interval = int(nframes/args.n)
extra = nframes%args.n

for start in range(0,nframes,interval):
	end = start+interval
	if end+interval+extra > nframes:
		end = nframes
	fname = '%s_%d_%d%s'%(basename,start,end,ext)
	with MDAnalysis.Writer(fname, model.atoms.n_atoms) as W:
		for ts in model.trajectory[start:end]:
			W.write(model.atoms)
