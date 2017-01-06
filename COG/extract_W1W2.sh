#! /bin/bash

#echo 'RUNNING extract_W1W2.sh'

files=`ls Traj*`

images=`ls Traj* | sed s/Traj_// | sed s/_MD.txt/.tif/`
imagesMD=`ls Traj* | sed s/Traj_// | sed s/.txt/.tif/`
ctrl=1
echo "Warning: it is assumed that you had red (W1) channel as first frame in the stack and the green (W2) channel as second frame"
for file in $files
do
	file_date=`echo $file | cut -c6- | sed 's/_.*$//'`
	
	grep "^0" $file | sed 's/^0 //' | sed 's/,//g' > detected_spots_W1
	grep "^1" $file | sed 's/^1 //' | sed 's/,//g' > detected_spots_W2

	#add the correct date suffix to transform.mat
	sed "s/[0-9_]*transform.mat/$file_date"_transform.mat"/" transform.m > tmp.m
	mv tmp.m transform.m
	#execute transform.m	
	matlab -nosplash < transform.m > /dev/null
	if [ $ctrl -le 9 ]; then 
		file_name=`echo "detected_spot_0$ctrl"`
	else 
		file_name=`echo "detected_spot_$ctrl"`
	fi

	mv detected_spots_W1 $file_name"_W1"
	mv detected_spots_W1_warped $file_name"_W1_warped"
	mv detected_spots_W2 $file_name"_W2"
	ctrl=`echo $ctrl+1 | bc -l`
done


ctrl=1
for image in $images
do
	if [ $ctrl -le 9 ]; then 
		image_name=`echo "image_0$ctrl"`
	else 
		image_name=`echo "image_$ctrl"`
	fi
	cp images/$image images/$image_name".tif"
	ctrl=`echo $ctrl+1 | bc -l`
done

ctrl=1
for image in $imagesMD
do
	if [ $ctrl -le 9 ]; then 
		image_name=`echo "imageMD_0$ctrl"`
	else 
		image_name=`echo "imageMD_$ctrl"`
	fi
	cp images/$image images/$image_name".tif"
	ctrl=`echo $ctrl+1 | bc -l`
done

		
