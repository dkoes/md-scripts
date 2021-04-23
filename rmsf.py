import sys
import MDAnalysis
import pandas as pd
from pandas import DataFrame
import numpy as np
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSF

def getTrajs():
	top = str(sys.argv[1])
	dcd_prefix = str(sys.argv[2])

	trajs = [];

	if (len(sys.argv) == 4):
		dcd_total = int(sys.argv[3])
		for i in range(dcd_total):
    			trajs.append(MDAnalysis.Universe(top, dcd_prefix + '{}.dcd'.format(i + 1)))
    		traj = MDAnalysis.Universe(top, [dcd_prefix + '{}.dcd'.format(i + 1) for i in range(10)])
	else:
		trajs.append(MDAnalysis.Universe(top, dcd_prefix))
		traj = trajs[0]

	return trajs, traj

def get_rmsf(t,label=None):
	sel = 'protein and not name H*'
	prot = t.select_atoms(sel)
    	average_coordinates = t.trajectory.timeseries(asel=prot).mean(axis=1)
    	# make a reference structure (need to reshape into a 1-frame "trajectory")
    	reference = MDAnalysis.Merge(prot).load_new(average_coordinates[:, None, :], order="afc")
    	aligner = align.AlignTraj(t, reference, select=sel, in_memory=True).run()
    	print('\tDone Aligning... calculating RMSF')
    	rmsfer = RMSF(prot).run()
    	ret = pd.DataFrame(zip(prot.resids,prot.resnames,rmsfer.rmsf),columns=('resid','resname','rmsf'))
    	if label != None:
        	ret['label'] = label
    	return ret


def rmsf_df(trajs, traj):
 	rmsfs = DataFrame()
 	protein = traj.select_atoms('protein')

 	rmsfs['resname'] = protein.residues.resnames
 	# rmsfs['res_num'] = protein.residues.resids

 	for i in range(len(trajs)):
		print('\tCalculating rmsf for traj {}'.format(i + 1))
    		rmsf = get_rmsf(trajs[i])
    		col_name = 'traj_%s' % str(i + 1)
    		rmsf = rmsf.groupby(['resid']).mean().reset_index().rename(columns={'rmsf': col_name})
    		rmsfs = pd.concat([rmsfs, rmsf], axis=1)
   	return rmsfs

if __name__ == '__main__' :
	print('Loading trajectories...')
	trajs, traj = getTrajs()
	print('Calculating rmsfs')
	rmsfs = rmsf_df(trajs, traj)
	print('Outputing csv...')
	rmsfs.to_csv('rmsfs.csv')
	print('Done')

