#!/usr/bin/env python

import MDAnalysis, argparse
from MDAnalysis.analysis.rms import rmsd
import MDAnalysis.analysis
import os
from os.path import splitext, basename

#This finds the most distinct set of frames.  Does not perform alignment, just rmsd.

def compute_cluster_sizes(current, model, selection):
	"""Given a set of frames (current), assign each from in
	the trajectory to its closest current frame to get cluster sizes and radii"""
	cnts = [0]*len(current)
	maxrmsd = [0]*len(current)
	for ts in model.trajectory:
		frame = ts.frame
		minr = float('infinity')
		closest = 0
		for (pos,(refframe, refcoords)) in enumerate(current):
			r = rmsd(refcoords, selection.positions)			
			if r < minr:
				minr = r
				closest = pos
		cnts[closest] += 1
		if minr > maxrmsd[closest]:
			maxrmsd[closest] = minr
	return (cnts,maxrmsd)
		
def add_next_farthest(current, model, selection):
	"""Given a list of frames (frame number, coordinates) current
		identify the frame that has the maximum minimum distance
		from current and add it to current"""
	maxr = 0
	maxframe = 0
	maxcoords = None
	for ts in model.trajectory:
		frame = ts.frame
		minr = float('infinity')
		for (refframe, refcoords) in current:
			r = rmsd(refcoords, selection.positions)			
			if r < minr:
				minr = r
		if minr > maxr:
			maxframe = frame
			maxr = minr
			maxcoords = selection.positions
	
	current.append( (maxframe, maxcoords) )


parser = argparse.ArgumentParser(description='Identify most distinct frames of md trajectory\nIMPORTANT: assumes an aligned input')
parser.add_argument("topology")
parser.add_argument("trajectory")
parser.add_argument("size",type=int)
parser.add_argument("--selection",default="backbone",required=False)
parser.add_argument("--output_selection",default="not (resname WAT or resname HOH)",required=False)

args = parser.parse_args()

model = MDAnalysis.Universe(args.topology,args.trajectory)
selection = model.select_atoms(args.selection)

current = [(0, selection.positions)]

for i in range(args.size-1):
	add_next_farthest(current, model, selection)

name = splitext(basename(args.topology))[0]

for (frame,coords) in current:
	model.trajectory[frame] #has side effect of setting current frame
	fname = "%s_%d.pdb" % (name, frame)
	model.select_atoms(args.output_selection).write(fname)
	os.system("sed -i '/REMARK/d' %s" % fname) #remove remarks with binary characters

#print out frame numbers and their respective cnts
(cnts,rmsds) = compute_cluster_sizes(current,model,selection)
for i in range(len(current)):
	(frame,coords) = current[i]
	cnt = cnts[i]
	rmsd = rmsds[i]
	print('%d:\t%d\t%.3f' % (frame,cnt,rmsd))
