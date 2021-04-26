'''
Calculates the average number of water contacts for each residue across the trajectory.
For each residue the waters are selected in a specified cutoff distance and then the counts are averaged
for that residue.

Output:
a .csv file with indexed by resid including residue names and the average water contacts
e.g.
resid | resname | waters
    1       LEU |    2.6 

Running the script
python waters.py MDAnalysis_supported_topolgy MDAnalysis_supported_trajectory csv_output_name.csv
e.g.
python contacts.py FXN_R165N_complex.prmtop FXN_R165N_complex1.dcd waters1.csv
'''
import sys
import MDAnalysis
import pandas as pd
from pandas import DataFrame
import numpy as np # numpy 1.16.6
from MDAnalysis.analysis import distances
import multiprocessing

'''
Load trajectory into curr_traj
'''
def getTraj():
	global dcd
	global curr_traj # must be global to work with multiprocessing
	top = str(sys.argv[1])
	dcd = str(sys.argv[2])
	traj = MDAnalysis.Universe(top, dcd)
	curr_traj = traj	


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


'''
calls get_res_contacts and returns numpy area of contacts for that frame
'''
def traj_frame_waters(i):
        curr_traj.trajectory[i]
        rc = get_res_waters(curr_traj)
        return rc.waters.to_numpy()


'''
Creates and returns the final dataframe with water contacts averaged for each residue over the
entire trajectory
'''
def waters_df():
 	waters = pd.DataFrame()
	protein = curr_traj.select_atoms('protein')
	waters['resid'] = protein.residues.resids
	waters['resname'] = protein.residues.resnames
 	print("Calculating waters...")
	pool = multiprocessing.Pool()
    	traj_waters = pool.map(traj_frame_waters,range(curr_traj.trajectory.n_frames))
    	traj_waters = np.array(traj_waters)
	waters['waters'] = (traj_waters.mean(axis=0))
   	return waters.set_index('resid')


if __name__ == '__main__' :
	print('Loading trajectories...')
	getTraj()
	print('Calculating waters...')
	waters = waters_df()
	print('Outputing csv...')
	waters.to_csv(sys.argv[3])
	print('Done.')
