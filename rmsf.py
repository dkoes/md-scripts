'''
Calculates the average RMSF for each residue in the trajectory using
MDAnalysis's RMSF method.

Output:
If called using shortcut for multiple residues, returned is a single csv called rmsfs.csv containg
the resname, resid, and trajectory number
e.g.
resname | resid | traj_1 | resid | traj_2 ...
LEU     |     1 |  1.272 |     1 |  1.410
...

If called using a single trajectory, returned is a csv called rmsfs.csv containing resid, resname
and RMSF

resid | resname | RMSF
    1 | LEU     | 1.27
...

Running the script:
Since calculating RMSF is not a very intense calculation, it may be useful to just
produce a single csv for all trajectories. If your trajetories are .dcd and have 
a naming scheme of prefix_number.dcd e.g. traj_1.dcd you can use the shortcut call.

python rmsf.py MDAnalysis_supported_topology trajectory_prefix total_number_of_trajs_starting_at_1
e.g. 
python rmsf.py FXN_R165N_complex.prmtop FXN_R165N_compex_ 10 

You can also call the script using a single trajectory:

python rmsfs.py MDAnalysis_supported_topology MDAnalysis_supported_trajectory
e.g.  
python rmsf.py FXN_R165N_complex.prmtop FXN_R165N_compex_1.dcd
'''
import sys
import MDAnalysis
import pandas as pd
from pandas import DataFrame
import numpy as np # numpy 1.16.6
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSF

'''
Loads trajectory(s) into trajs
'''
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


'''
Calculates the avergae RMSF for each residue in the trajectory
'''
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


'''
Creates and returns the final dataframe to become rmsfs.csv
'''
def rmsf_df(trajs, traj):
 	rmsfs = DataFrame()
 	protein = traj.select_atoms('protein')

 	rmsfs['resname'] = protein.residues.resnames

 	for i in range(len(trajs)):
		print('\tCalculating rmsf for traj {}'.format(i + 1))
    		rmsf = get_rmsf(trajs[i])
		if (len(trajs) == 1):
			col_name = 'RMSF'
		else:
    			col_name = 'traj_%s' % str(i + 1)
    		rmsf = rmsf.groupby(['resid']).mean().reset_index().rename(columns={'rmsf': col_name})
    		rmsfs = pd.concat([rmsfs, rmsf], axis=1)
   	if (len(trajs) == 1):
		return rmsfs.set_index('resid')
	else:
		return rmsfs


if __name__ == '__main__' :
	print('Loading trajectories...')
	trajs, traj = getTrajs()
	print('Calculating rmsfs')
	rmsfs = rmsf_df(trajs, traj)
	print('Outputing csv...')
	rmsfs.to_csv('rmsfs.csv')
	print('Done')

