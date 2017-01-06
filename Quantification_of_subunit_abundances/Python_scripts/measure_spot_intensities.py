import numpy as np
from os import listdir as ls
from skimage import filters #import filters
from skimage import morphology
from skimage.measure import label
from skimage.exposure import rescale_intensity
from skimage import img_as_float, img_as_uint
import skimage.external.tifffile as tiff
import skimage.io as io

#--------------------------------------------------
#	image correction
#--------------------------------------------------
def img_corr(img,radius):
	#define imgage objects
	img_corrected = np.zeros(shape=img.shape,dtype=img.dtype)#where the imgage corrected for photobleaching will be stored
	median_cytoplasmatic_bkg=np.zeros(shape=img.shape[0])#where to store the median value of the cytoplasmatic background of the cells
	for i in range(img.shape[0]):
		#compute the median of the cell
		img_median=filters.median(img[i,:,:],morphology.disk(radius))
		#compute the imgage corrected without cytoplasmatik background. All maths are done on signed float dtype and conferted in unsinged 16 bit format_CORRECTED
		img_corrected[i,:,:]=img_as_uint(
			(img_as_float(img[i,:,:])-img_as_float(img_median))
			)
	return img_corrected

def load_image(path,radius=17):
	im = tiff.imread(path)
	#images acquired with the CCD are on a 12 bit range and can be rescaled to 16 bit
	im_rescaled = rescale_intensity(im,in_range=(0,2**12-1),out_range='uint16')
	return img_corr(im_rescaled,radius)
	#return img_corr(im,radius)

#--------------------------------------------------
#	mask patches, used to identify
#	the good GFP patches to be 
#	quantified and to track the patches
#	through the stack.
#--------------------------------------------------
def dilation(input,iterations):
	iter=0
	while iter < iterations:
		input=morphology.dilation(input, morphology.ball(2)).astype(input.dtype)
		iter+=1
	return input

def threshold_erosion(input):
	i = 0;
	while i < input.shape[0]:
		if i == 12:
			io.imshow(input[i,:,:])
			io.show()
		disk_erosion = morphology.erosion(input[i,:,:],morphology.disk(1)).astype(input.dtype)
		square_erosion = morphology.erosion(input[i,:,:],morphology.square(2)).astype(input.dtype)
		frame_erosion = disk_erosion
		frame_erosion[ square_erosion > 0 ] = 2**16-1
		input[i,:,:] = frame_erosion
		i+=1
	return input

def erosion( image , n ) :

	if n == 0:
		n = 1
		print("n set to 1; eroding pixels with no neighbor makes no sense")
	if n > 26:
		n = 26
		print("n set to 26; number of neighbor pixels cannot exceed 26")

	brush = np.array([
			[#0
				[
					[ 1 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#1
				[
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#2
				[
					[ 0 , 0 , 1 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#3
				[
					[ 0 , 0 , 0 ],
					[ 1 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#4
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#5
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 1 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#6
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 1 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#7
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#8
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 1 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#9
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 1 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#10
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 1 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#11
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 1 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#12
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 1 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#13
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 1 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#14
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 1 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#15
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 1 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#16
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 1 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#17
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 1 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#18
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#19
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 1 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#20
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 1 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#21
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#22
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 1 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#23
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 1 , 0 , 0 ]] ,
				],
			[#24
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ]] ,
				],
			[#25
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 1 ]] ,
				]
			])
					
	image_output = np.zeros( shape = image.shape , dtype = image.dtype )
	
	eroded_images=np.zeros(shape=image.shape,dtype=image.dtype)
	
	for j in range(26):
	
		tmp_image = morphology.binary_erosion( image , brush[j,:,:,:] )
		eroded_images = eroded_images + tmp_image
		
		image_output[ eroded_images > 26 - n ] = 1
	
	return image_output


def mask_patches(image,threshold_value=None):
	#make a mask from the image image
	if threshold_value == None: 
		threshold_value=filters.threshold_yen( image )

	print 'Yen threshold: '+str(threshold_value)
	image_threshold=np.zeros(shape=image.shape,dtype=image.dtype)
	image_threshold[image > threshold_value ]=1

	image_eroded = erosion( image_threshold , 21 )
	
	return image_eroded

#--------------------------------------------------
#	Track the patches through the stack
#--------------------------------------------------
def measure_spot_intensities(image,image_mask,ctrl_mask):
	#Measure the spot intensities that are present in both a mask of the image and
	#a mask of a control image (the RFP images for the vacuole proteins).

	#label the mask of the image to be quantified and the mask of the RFP ctrl image.
	image_label=label(image_mask,connectivity=2)
	ctrl_label=label(ctrl_mask,connectivity=2)

	measurements=np.zeros(ctrl_label.max())
	for i in range(ctrl_label.max()):
		spot_at_the_edge=False #Define the variable that record if spots are at the edge of the z-stack
		for j in range(image_label.max()):
			match=image_label[np.all([ctrl_label==i+1,image_label==j+1],axis=0)]
			if match.size:
				#before proceeding with the quantification check that the 
				#spot to be quantified is not truncated at the beginning or 
				#at the end of the stack i.e.: that the label of the spot 
				#mask is not present either in any of the array entries of 
				#the first frame (0) or in any of those of the last frame 
				#(image_label.size[0]-1). The presence of a spot at the edge
				#is recorded as a variable, but the spot is not excluded yet
				#to allow the algorith to recongize whether there are multiple
				#spots within the patch of the ctrl_mask. If the spot selection
				#would be done at this point one could not exculde that an other
				#GFP patch, which colocalize within the same patch of the 
				#ctrl_mask (hence too close), would not start or end on the
				#edges of the stack.
				spot_at_the_edge=np.any([np.any(image_label[0,:,:]==j+1),np.any(image_label[image_label.shape[0]-1,:,:]==j+1)])
				#control it is not on the border of the frame
				if not spot_at_the_edge: 
					spot_at_the_edge=np.any([np.any(image_label[:,0,:]==j+1),np.any(image_label[:,image_label.shape[1]-1,:]==j+1)])
				if not spot_at_the_edge: 
					spot_at_the_edge=np.any([np.any(image_label[:,:,0]==j+1),np.any(image_label[:,:,image_label.shape[2]-1]==j+1)])

				#I want to avoid patches that are too close. Those patches
				#will have the same control_label, which is generated from
				#a dilated mask, but different image_labels. If close patches
				#exist then the value measurements[i] would be assigned more 
				#than once. If this happens, at the second assignment the 
				#value of the measurement will be non-zero. Therefore the 
				#measurement[i] will be deleted (assigned 0) and the loop
				#over range(image_labels.max()) will exit to proceed with the
				#next patch in the ctrl_mask.
				if not measurements[i]:
					measurements[i]=image[image_label==j+1].sum()
				else:
					measurements[i]=0
					break
		if spot_at_the_edge: 
			measurements[i]=0

	return measurements

#--------------------------------------------------
#	Analyse the images, with 
#	or without RFP mask
#--------------------------------------------------

def analysis(path_in,radius=17,GFP_pattern='GFP-FW',save_masks=True):
	#define the array in which all the measurements will be stored
	output_measurements=np.zeros(0)
	images=ls(path_in)
	
	GFP_images=[img for img in images if GFP_pattern in img]
	
	for i in range(len(GFP_images)):
		GFP_im = load_image(path_in+GFP_images[i],radius)
		print(path_in+GFP_images[i])
		
		mask = mask_patches(GFP_im)
		
		ctrl_mask =  dilation(mask,iterations = 2)
		#save the ctrl mask	
		if save_masks: tiff.imsave(path_in+GFP_images[i].replace(GFP_pattern,'_CtrlMask'),ctrl_mask)
		
		#comput the mask of the image and save it
		image_mask= dilation(mask,iterations=1)
		if save_masks: tiff.imsave(path_in+GFP_images[i].replace(GFP_pattern,'_ImageMask'),image_mask)
		output_measurements=np.concatenate((
			output_measurements,
			measure_spot_intensities(GFP_im,image_mask,ctrl_mask)
			))
	return output_measurements

def experiment(path,target_name,reference_name='Nuf2',GFP_pattern='GFP-FW',median_radius=17):
	reference=analysis(path+'/'+reference_name+'/',radius=median_radius,GFP_pattern=GFP_pattern)
	np.savetxt(path+'/'+reference_name+'_intensities.txt',reference)
	
	target=analysis(path+'/'+target_name+'/',GFP_pattern=GFP_pattern)
	np.savetxt(path+'/'+target_name+'_intensities.txt',target)
	
	return(reference,target)
	

