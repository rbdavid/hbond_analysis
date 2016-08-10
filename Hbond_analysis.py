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
from sel_list import *

# VARIABLE DECLARATION:

pdb = sys.argv[1]                   # point to a pdb or prmtop or psf file (untested for both prmtop and psf files)
traj_loc = sys.argv[2]              # point to the location of the trajectory files
start = int(sys.argv[3])
end = int(sys.argv[4])
system = sys.argv[5]

nSel = len(sel)

flush = sys.stdout.flush
mkdir = os.mkdir
chdir = os.chdir

solvent_resname = 'WAT'

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
    print '%s' %(string)
    flush()

def summary(nSteps):
    sum_file = open('%s.hbond.summary' %(system),'w')
    sum_file.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
    sum_file.write('To recreate this analysis, run this line in terminal:\n .hbond.analysis.py pdb_file trajectory_location start_traj end_traj system_descriptor')
    sum_file.write('\n\n')
    sum_file.write('Progress output is written to:\n ')
    sum_file.write('\nTotal number of steps analyzed: %d\n' %(nSteps))
    sum_file.write('\nAtom Selections analyzed:n')
    for i in range(nSel):
        sum_file.write('%02d  %s %s\n' %(i,sel[i][1],sel[i][2],sel[i][3]))
    sum_file.close()

# ----------------------------------------
# MAIN:
# ----------------------------------------
# INITIALIZING UNIVERSE
u = MDAnalysis.Universe(pdb)
h_list = []
for i in range(nSel):
	if sel[i][0] != sel[i-1][0] or i == 0:			# NOT GOOD BOOLEANS TO TEST
		mkdir('%s' %(sel[i][0]))
	
	if sel[i][3] == solvent_resname:
		h = hbonds.HydrogenBondAnalysis(u,selection1=sel[i][2],selection2=sel[i][3],selection1_type='both',update_selection1=False,update_selection2=True,detect_hydrogens='distance',distance=3.0,angle=120.0,donors=sel[i][4],acceptors=sel[i][5])
	else:
		h = hbonds.HydrogenBondAnalysis(u,selection1=sel[i][2],selection2=sel[i][3],selection1_type='both',update_selection1=False,update_selection2=False,detect_hydrogens='distance',distance=3.0,angle=120.0,donors=sel[i][4],acceptors=sel[i][5])

    h_list.append(h)

nSteps = 0
while start <= end:
	u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))
	nSteps += len(u.trajectory)
	start += 1

# BEGINNING TO ANALYZE TRAJECTORIES
start = int(sys.argv[3])
count = 0
while start <= end:
	ffprint('Loading trajectory %s' %(start))
	u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))

# Loop through all residue pair selections, calculate the Hydrogen Bond distances and angles for all pairs;  
	for i in range(nSel):
		chdir('%s' %(sel[i][0]))
		out1 = open('%s.%s.results.dat' %(sel[i][1]),system,'a')
		h_list[i].run()
		htimeseries = h_list[i].timeseries     
		temp = count
		for j in range(len(u.trajectory)):
			if len(htimeseries[j]) != 0:
				for k in range(len(htimeseries[j])):
					out1.write('%10d   %10d   %10d    %f   %f\n' %(temp,htimeseries[j][k][0],htimeseries[j][k][1],htimeseries[j][k][-2],htimeseries[j][k][-1]))
				temp += 1
			else:
				temp += 1
		out1.close()
# Change directories back to parent directory, iterate to the next trajectory....REPEAT 
		ffprint('Finished analyzing trajectory %02d\n' %(start_traj))
		os.chdir('..')
	count += len(u.trajectory)
	start_traj += 1

ffprint('Analyzed %d steps.' %(nSteps))
summary(nSteps)

