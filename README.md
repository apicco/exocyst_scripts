# exocyst_scripts
The collection of scripts used to compute distances between diffraction limited fluorescent spots of distinct fluorophores.
The distances were used to model the 3D architecture of the Exocyst and Cog complexes published in 

Andrea Picco, Ibai Irastorza-Azcarate, Tanja Specht, Dominik Böke, Irene Pazos, Anne-Sophie Rivier-Cordey, Damien P. Devos, Marko Kaksonen, Oriol Gallego, _The in vivo architecture of the exocyst provides structural basis for exocytosis_, Cell; DOI:_to be announced_

This is the repository to which the Cell _STAR methods_ points to.

# EXOCYST folder
This folder constains the scripts that compute the distances used to compute the molecular architecture of the exocyst complex. Note that this scripts run on old versions of R. 

# COG folder
This folder containts the script used to comput COG distances. The script has been adapted to newer versions of R and it also contains small differences to accommodate the difference in the image quality from the Exocsyt complex.

# Chromatic_aberration_correction 
This folder contains the script used to compute the chromatic aberration correction. It uses [Matlab](www.mathworks.com) and assumes that the red channel is the first frame (window 1, W1, in our convention), while the green channel is the second frame (W2). Obviously, this convention *must* hold true for all the images.
 
#Tracking
Tracking was done using [Particle Tracker](http://imagej.net/Particle_Tracker) from the [MOSAIC group](http://mosaic.mpi-cbg.de/MosaicToolboxSuite/ParticleTracker.html). The version of Particle Tracker used is the Version 1.5 from September 2006. This version was modified by us to do not prompt the trajectory GUI and to compute additional parameters for the spots, such as the central momenti of brigthness:

	//central moment of order 1,3,4,5
	this.particles[m].m1 += (float)(k + l) * c;
	this.particles[m].m3 += (float)(k * k * k + l * l * l) * c;
	this.particles[m].m4 += (float)(k * k * k *k + l * l * l * l) * c;
	this.particles[m].m5 += (float)(k * k * k * k *k + l * l * l * l * l) * c;
	//orientation and eccentricity of the spot
	//central momemts of order 11,02,20
	this.particles[m].mu11 += (float)(k * l * c);
	this.particles[m].mu20 += (float)(k * k * c);
	this.particles[m].mu02 += (float)(l * l * c); 

The software was run using a shell script. Here, we provide an example:

	#!/bin/bash
	#choose an option:
	#images: extract the centroids from the images and organize the datas, including the beads corrections
	#quant: run the beads correction and the R scripts for the quantification in each folder
	#all: is images + quant
	#then choose the dataset, either C or N
	#it is recommended to record the outputs of the the software on log, to check if some analysis failed:
	#bash au
	set -e
	set -u
	
	dobeads="FALSE"
	#dobeads="TRUE"
	#control_name="B9" #control name for the C terminus
	control_name="1p6" #control name for the N terminus
	
	#DIRECTORIES 
	current_directory=`pwd`
	software_directory="$current_directory/Analysis_SOFTWARE"
	
	if [[ $2 == "C" ]]
	then
		exocyst_image="$HOME/Desktop/Cterm/"
		exocyst_data="$HOME/Desktop/COG/C/"
		beads="$current_directory/BEADS/C/"
	elif [[ $2 == "N" ]]
	then
		exocyst_image="$HOME/Desktop/Nterm/"
		exocyst_data="$HOME/Desktop/COG/N/"
		beads="$current_directory/BEADS/N/"
	else
		echo "select whether is the C or N terminal"
	fi
	
	#list all the folders
	image_folders=`ls  "$exocyst_image" | grep -v BEADS`
	output=distances.data
	
	
	#sort the image folders
	if [[ $1 == "images" || (( $1 == "all" )) ]]
	then
		for  image_folder in $image_folders
		do
			cd "$exocyst_image$image_folder"
			image_files=`ls *_[0-9][0-9].tif`
			#run imageJ plugin
			for image_file in $image_files
			do
				echo $image_file
				#run imageJ
				#OLD java -jar /Applications/ImageJ/ij.jar $image_file -ijpath /Applications/ImageJ/ -macro /Applications/ImageJ/plugins/folder_BGS_MD
			   	java -jar /Applications/ImageJ/ImageJ.app/Contents/Resources/Java/ij.jar "$image_file" -ijpath /Applications/ImageJ/ -macro /Users/*/Google\ Drive/RESEARCH/COG_COMPLEX/ImageJ_PLUGINS/i_CRAI_Tracking.txt
			done
		
			#extract the code for the experiment (for example C1 in 111103_C1) so that experiments that had to be acquired on two days are the listed in the same folder. 
			#Note that each day has to be corrected with its corresponding transformation. 
			#This does not apply to $control_name, our internal control
		
			image_folder_name=`echo $image_folder | cut -c8-`
			image_folder_date=`echo $image_folder | cut -c-6`
		
			#if the folder $image_folder_name is not yet present
			if [ $image_folder_name == $control_name ]; then
				echo $exocyst_data$image_folder
				mkdir "$exocyst_data$image_folder"
				mkdir "$exocyst_data$image_folder/images"
				#load the desired software	
				cp "$software_directory/analyse_the_data.R" "$exocyst_data$image_folder"
				cp "$software_directory/transform.m"  "$exocyst_data$image_folder"
				cp "$software_directory/extract_W1W2.sh"  "$exocyst_data$image_folder"
				cp "$software_directory/parameters.data"  "$exocyst_data$image_folder"
				mv Traj_* "$exocyst_data$image_folder"
				mv *_BGN*tif "$exocyst_data$image_folder/images"
				cp "$beads$image_folder_date""_transform.mat" "$exocyst_data$image_folder"
			elif [ -z $(ls | grep -x "$exocyst_data$image_folder_name") ]; then
				echo "$exocyst_data$image_folder_name"
				mkdir "$exocyst_data$image_folder_name"
				mkdir "$exocyst_data$image_folder_name/images"
				#load the desired software	
				cp "$software_directory/analyse_the_data.R" "$exocyst_data$image_folder_name"
				cp "$software_directory/transform.m" "$exocyst_data$image_folder_name"
				cp "$software_directory/extract_W1W2.sh" "$exocyst_data$image_folder_name"
				cp "$software_directory/parameters.data"  "$exocyst_data$image_folder_name"
				mv Traj_* "$exocyst_data$image_folder_name"
				mv *_BGN*tif "$exocyst_data$image_folder_name/images"
				cp "$beads$image_folder_date""_transform.mat" "$exocyst_data$image_folder_name"
			else
				mv Traj_* "$exocyst_data$image_folder_name"
				mv *_BGN*tif "$exocyst_data$image_folder_name/images"
				cp "$beads$image_folder_date""_transform.mat" "$exocyst_data$image_folder_name"
			fi
			#rm *SATURATION.tif
		done
		echo "IMAGE ANALYSIS DONE"
	fi
	
	data_folders=`ls -d "$exocyst_data"*/`
	#if [[ $dobeads == "TRUE" || (( $1 == "all" )) ]]
	if [[ $dobeads == "TRUE" ]]
	then
		i=0
		echo "BEADS TRANSFORM"
		for data in $data_folders
		do
			i=`echo $i+1 | bc -l`
			if [ $i == 7 ] 
			then
				wait
				i=1
			fi
			cd "$data"
			echo $i" "$data
			bash extract_W1W2.sh &
		done
	fi
	wait
	if [[ $1 == "quant" || (( $1 == "all" )) ]]
	then
		i=0
		echo "R computation"
		for data in $data_folders
		do
			i=`echo $i+1 | bc -l`
			if [ $i == 7 ] 
			then
				wait
				i=1
			fi
			cd "$data"
			#if it's not all copy the files in
			if [[ $1 == "quant" ]];
			then
				cp "$software_directory/analyse_the_data.R" ./
				cp "$software_directory/transform.m" ./
				cp "$software_directory/extract_W1W2.sh" ./
			fi
	
			echo $i" "$data
			( { R batch --slave --quiet --vanilla <analyse_the_data.R >analysis.log; } || { echo $i" FAILED"; }
			echo $i" DONE" ) &
		done
		wait
		echo "folder dist est.se sigma est.se n" > $exocyst_data$output
		for data in $data_folders
		do
			cd "$data"
			data_folder_name=`pwd | sed 's/\/.*\///'`
			echo $data_folder_name
			if [ -a output.pdf ]; then mv output.pdf "$data_folder_name"".pdf"; fi
			#rm images/image*.tif
			if [ -f output.data ];
			then
				results=`cat output.data`
				cd "$exocyst_data"
				echo "$data_folder_name"" ""$results" >> "$output"
			else
				cd "$exocyst_data"
				echo $data_folder_name" FAILED" >> "$output"
			fi
		done
		cd "$exocyst_data"
		rm -f output.txt
		#parsing of the data
		cp "$software_directory/parser.sh" ./
		cp "$current_directory/Doc/Plate_design.csv" ./
		bash parser.sh
		rm parser.sh
		rm Plate_design.csv
		#finish
		echo "QUANTIFICATION DONE"
		
	fi
	
	if [[ $1 == "parse" ]]
	then
		for data in $data_folders
		do
			cd $data
			data_folder_name=`pwd | sed 's/\/.*\///'`
			echo $data_folder_name
			if [ -a output.pdf ]; then mv output.pdf $data_folder_name".pdf"; fi
			#rm images/image*.tif
			if [ -f output.data ];
			then
				results=`cat output.data`
				cd $exocyst_data
				echo $data_folder_name" "$results >> $output
			else
				cd $exocyst_data
				echo $data_folder_name" FAILED" >> $output
			fi
		done
	
		cd $exocyst_data
	#TO BE ADJUSTED# 	rm -f output.txt
	#TO BE ADJUSTED# 	#parsing of the data
	#TO BE ADJUSTED# 	cp "$software_directory/parser.sh" ./
	#TO BE ADJUSTED# 	cp "$current_directory/Doc/Plate_design.csv" ./
	#TO BE ADJUSTED# 	bash parser.sh
	#TO BE ADJUSTED# 	rm parser.sh
	#TO BE ADJUSTED# 	rm Plate_design.csv
		#finish
	fi
	
		if [[ $1 != "quant" &&  $1 != "all"  && $1 != "parse" && $1 != "images" ]];
	then
		echo "choose an option:\nimages\t extract the centroids from the images and organize the datas, including the beads corrections\nquant\t run the beads correction and the R scripts for the quantification in each folder\nall\t is images + quant\nparse does the parsing only"
	fi
	cd "$current_directory"
	
