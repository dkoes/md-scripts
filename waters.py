import sys
import MDAnalysis
import pandas as pd
from pandas import DataFrame
import numpy as np
from MDAnalysis.analysis import distances
import multiprocessing


def getTraj():
	global dcd
	global curr_traj
	top = str(sys.argv[1])
	dcd = str(sys.argv[2])
	traj = MDAnalysis.Universe(top, dcd)
	curr_traj = traj	
	return traj


def get_res_waters(traj,cutoff=2):
	'''Given an mdanalysis trajectory and contact distance cutoff,
	count the number of water molecules within cutoff of each residue'''
	protein = traj.select_atoms('protein')
	resids = np.unique(protein.resids)
	ret = []
	for rid in resids:
		wat = traj.select_atoms('resname WAT and around %f resid %d'%(cutoff,rid))
		ret.append((rid,len(wat)))
	return pd.DataFrame(ret,columns=('resid','waters'))


def traj_frame_waters(i):
        curr_traj.trajectory[i]
        rc = get_res_waters(curr_traj)
        return rc.waters.to_numpy()


def waters_df():
 	waters = pd.DataFrame()
	protein = curr_traj.select_atoms('protein')
	waters['resid'] = protein.residues.resids
	waters['resname'] = protein.residues.resnames
 	print("Calculating waters for traj_%s" % dcd.split('_')[-1][0])
	pool = multiprocessing.Pool()
    	traj_waters = pool.map(traj_frame_waters,range(curr_traj.trajectory.n_frames))
    	traj_waters = np.array(traj_waters)
	waters['traj_{}'.format(dcd.split('_')[-1][0])] = (traj_waters.mean(axis=0))
   	return waters.set_index('resid')


if __name__ == '__main__' :
	print('Loading trajectories...')
	traj = getTraj()
	print('Calculating waters...')
	waters = waters_df()
	print('Outputing csv...')
	waters.to_csv(sys.argv[3])
	print('Done.')
