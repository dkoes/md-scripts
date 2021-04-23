import sys
import MDAnalysis
import pandas as pd
from pandas import DataFrame
import numpy as np
from MDAnalysis.analysis import distances
import multiprocessing
import os


def getTraj():
	global dcd
	global curr_traj
	top = str(sys.argv[1])
	dcd = str(sys.argv[2])
	traj = MDAnalysis.Universe(top, dcd)
	curr_traj = traj	
	return traj


def get_rescontacts(protein,cutoff=2):
	n = len(protein)
	self_distances = distances.self_distance_array(protein.positions)
	sq_dist_arr = np.zeros((n, n), dtype=np.float32)
	triu = np.triu_indices_from(sq_dist_arr, k=1)
	sq_dist_arr[triu] = self_distances
	sq_dist_arr.T[triu] = self_distances
	contacts = np.where(sq_dist_arr < cutoff, 1, 0)
	df = pd.DataFrame(zip(protein.resids,contacts.sum(axis=0)),columns=('resid','contacts'))
	rescontacts = df.groupby('resid').sum().reset_index()
	return rescontacts


def traj_contacts(i):
	print("Frame: %s" % i)
	protein = curr_traj.select_atoms('protein')
	curr_traj.trajectory[i]
	rc = get_rescontacts(protein)
	return rc.contacts.to_numpy()


def contacts_df():
 	contacts = pd.DataFrame()
	protein = curr_traj.select_atoms('protein')
	contacts['resid'] = protein.residues.resids
	contacts['resname'] = protein.residues.resnames
 	print("Calculating contacts for traj_%s" % dcd.split('_')[-1][0])
 	pool = multiprocessing.Pool(5)
 	traj_contacts_res = pool.map(traj_contacts, range(curr_traj.trajectory.n_frames))
 	traj_contacts_res = np.array(traj_contacts_res)
 	contacts['traj_{}'.format(dcd.split('_')[-1][0])] = (traj_contacts_res.mean(axis=0))
 	return contacts.set_index('resid')


if __name__ == '__main__' :
	print('Loading trajectories...')
	traj = getTraj()
	print('Calculating contacts...')
	contacts = contacts_df()
	print('Outputing csv...')
	contacts.to_csv(sys.argv[3])
	print('Done.')
