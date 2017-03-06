#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:

# ./hbond.analysis.py pdb_file trajectory_location start_traj end_traj system_descriptor

# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import MDAnalysis
import MDAnalysis.analysis.hbonds as hbonds

# VARIABLE DECLARATION:
config_file = sys.argv[1]

necessary_parameters = ['pdb','traj_loc','start','end','selection_input_filename','hbond_distance_max','hbond_angle_min']
all_parameters = ['pdb','traj_loc','start','end','selection_input_filename','hbond_distance_max','hbond_angle_min']

flush = sys.stdout.flush

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
    print '%s' %(string)
    flush()

def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''

	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
#	parameters['alignment'] = 'protein'
#	parameters['solvent_resname'] = 'WAT'
#	parameters['write_summary'] = False
#	parameters['summary_filename'] = None

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
			sys.exit()

#def summary(nSteps):
#    sum_file = open('%s.hbond.summary' %(system),'w')
#    sum_file.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
#    sum_file.write('To recreate this analysis, run this line in terminal:\n .hbond.analysis.py pdb_file trajectory_location start_traj end_traj system_descriptor')
#    sum_file.write('\n\n')
#    sum_file.write('Progress output is written to:\n ')
#    sum_file.write('\nTotal number of steps analyzed: %d\n' %(nSteps))
#    sum_file.write('\nAtom Selections analyzed:n')
#    for i in range(nSel):
#        sum_file.write('%02d  %s %s %s\n' %(i,sel[i][1],sel[i][2],sel[i][3]))
#    sum_file.close()

# ----------------------------------------
# MAIN:
# ----------------------------------------
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)
from paramters['sel_list'] import *

nSel = len(sel)
hbond_distance_max = float(paramters['hbond_distance_max'])
hbond_angle_min = float(paramters['hbond_angle_min'])

# ----------------------------------------
# INITIALIZING UNIVERSE
u = MDAnalysis.Universe(parameters['pdb'])

h_sel_list = []
output_file_list = []
count = 0
for i in sel:
	h = hbonds.HydrogenBondAnalysis(u,selection1=i[1],selection2=i[2],selection1_type='both',update_selection1=False,update_selection2=False,detect_hydrogens='distance',distance=hbond_distance_max,angle=hbond_angle_min,donors=i[3],acceptors=i[4])
	h_list.append(h)
	temp = open('%02d.%s.num_hbonds.dat' %(count,i[0]),'w')
	output_file_list.append(temp)
	count += 1

# BEGINNING TO ANALYZE TRAJECTORIES
start = int(parameters['start'])
end = int(parameters['end'])
nSteps = 0
while start <= end:
	ffprint('Loading trajectory %s' %(start))
	u.load_new('%sproduction.%s/production.%s.dcd' %(parameters['traj_loc'],start,start))
	nSteps += len(u.trajectory)
	# Loop through all residue pair selections, calculate the Hydrogen Bond distances and angles for all pairs;  
	for i in range(nSel):
		h_list[i].run()
		# SAVE THE NUMBER OF HBONDS IN THE TIMESERIES
		np.savetxt(output_file_list[i],[int(j[1]) for j in h.count_by_time()])
		
		#htimeseries = h_list[i].timeseries     
		#for j in range(len(u.trajectory)):
		#	if len(htimeseries[j]) != 0:
		#		for k in range(len(htimeseries[j])):
		#			out1.write('%10d   %10d   %10d    %f   %f\n' %(temp,htimeseries[j][k][0],htimeseries[j][k][1],htimeseries[j][k][-2],htimeseries[j][k][-1]))
		#		temp += 1
		#	else:
		#		temp += 1
		#out1.close()
# Change directories back to parent directory, iterate to the next trajectory....REPEAT 
		#os.chdir('..')

	ffprint('Finished analyzing trajectory %02d\n' %(start_traj))
	start_traj += 1

ffprint('Analyzed %d steps.' %(nSteps))
#summary(nSteps)

