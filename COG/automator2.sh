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
