#! /bim/bash

detected_spots=`ls *FP_0`
for detected_spot in $detected_spots
do
	echo $detected_spot
	grep -v "[a-z]" $detected_spot | sed 's/,//g' > tmp
	mv tmp  $detected_spot
done
matlab -nodesktop -nosplash < generate_transform.m > /dev/null
matlab -nodesktop -nosplash < transform.m > /dev/null

dir_name=`pwd | sed 's/\/.*\///'`
cp transform.mat ../$dir_name"_transform.mat"
