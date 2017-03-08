#!/bin/bash

START=$1
STOP=$2
TRAJ_LOC=$3
TOT_NUM_RESIDUES=$4
DIST_CUTOFF=$5
ANGLE_CUTOFF=$6

for ((i=$START;i<=$STOP;i++))
do
	printf -v x "%03d" $i
	mkdir $x.hbonding_analysis
	cd $x.hbonding_analysis
	traj_file=../$TRAJ_LOC/production.$i/production.$i.dcd
	sed -e s?AAA?$traj_file?g -e s/BBB/$TOT_NUM_RESIDUES/g -e s/CCC/$DIST_CUTOFF/g -e s/DDD/$ANGLE_CUTOFF/g < ../sample_cpptraj.in > $x.hbond.in
	time cpptraj -i $x.hbond.in > $x.hbond.log
	echo Done with Trajectory $i !
	cd ../
done

