'''
Calulcates the average number of contacts with each protein atom across the trajectory.
For each frame, these interactions are grouped and summed by residue id to give give the
total number of atomic contacts for each residue. The contacts are then averaged acorss 
all frames to give the average total contacts for each residue.

tldr: Calculates the average total number of protein atomic contacts for each residue in the protein


Output: a single .csv file indexed by resid including the residue names and the average contacts
e.g.
resid | resname | contacts
    1 |     ASN |    96.53
...

Running the script:
python contacts.py MDAnalysis_supported_topolgy MDAnalysis_supported_trajectory csv_output_name.csv
e.g.
python contacts.py FXN_R165N_complex.prmtop FXN_R165N_complex_1.dcd contacts1.csv
'''

import sys
import MDAnalysis
import pandas as pd
from pandas import DataFrame
import numpy as np # numpy 1.16.6 
from MDAnalysis.analysis import distances
import multiprocessing
import os

'''
Loads the trajectory into curr_traj
'''
def getTraj():
	global dcd
	global curr_traj # must be glocal to work with multiprocessing
	top = str(sys.argv[1])
	dcd = str(sys.argv[2])
	traj = MDAnalysis.Universe(top, dcd)
	curr_traj = traj	


'''
Creates a distance matrix with all atoms then counts contacts if the distance
is less than the cutoff. Contacts are then grouped with and summed by resid 
and returned.
'''
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


'''
Calls the get_contacts funtion frame by frame
'''
def traj_contacts(i):
	print("Frame: %s" % i)
	protein = curr_traj.select_atoms('protein')
	curr_traj.trajectory[i]
	rc = get_rescontacts(protein)
	return rc.contacts.to_numpy()


'''
Creates and returns the final dataframe with contacts averaged for each 
residue over the entire trajectory.
'''
def contacts_df():
 	contacts = pd.DataFrame()
	protein = curr_traj.select_atoms('protein')
	contacts['resid'] = protein.residues.resids
	contacts['resname'] = protein.residues.resnames
 	print("Calculating contacts...")
 	pool = multiprocessing.Pool(5)
 	traj_contacts_res = pool.map(traj_contacts, range(curr_traj.trajectory.n_frames))
 	traj_contacts_res = np.array(traj_contacts_res)
 	contacts['contacts'] = (traj_contacts_res.mean(axis=0))
 	return contacts.set_index('resid')



if __name__ == '__main__' :
	print('Loading trajectories...')
	getTraj()
	print('Calculating contacts...')
	contacts = contacts_df()
	print('Outputing csv...')
	contacts.to_csv(sys.argv[3])
	print('Done.')
