#!/bin/bash

folder_name=`ls | grep -x detected_beads_original`
if [ -z "$folder_name" ]; then 
	filesW1=`ls RFP_[1-9]`
	filesW2=`ls GFP_[1-9]`
	mkdir detected_beads_original
	
	grep ^[0-9] RFP_0 > tmpW1
	mv RFP_0 detected_beads_original/
	grep ^[0-9] GFP_0 > tmpW2
	mv GFP_0 detected_beads_original/
	
	for fileW1 in $filesW1
	do
		echo $fileW1
		grep ^[0-9] $fileW1 >> tmpW1
		mv $fileW1 detected_beads_original/
	done
	for fileW2 in $filesW2
	do
		echo $fileW2
		grep ^[0-9] $fileW2 >> tmpW2
		mv $fileW2 detected_beads_original/
	done

	mv tmpW2 GFP_0
	mv tmpW1 RFP_0
fi
