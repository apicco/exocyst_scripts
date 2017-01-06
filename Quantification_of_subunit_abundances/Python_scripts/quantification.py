#import measure_spot_intensities as msi
import measure_spot_intensities as msi
import numpy as np

#Sec5
images_path='/Volumes/MarkoKaksonenLab/Andrea/T1/Andrea/Original_Data/150617_2661_2920/CQ_out_of_CQ/'
msi.experiment(path=images_path,target_name='Sec5',reference_name='Nuf2',GFP_pattern='w1GFP',median_radius=9)

##Sec6
images_path='/Volumes/MarkoKaksonenLab/Andrea/UNIGE/MICROSCOPY_DATA/160804_2920_2662looped_out/Cells/'
msi.experiment(path=images_path,target_name='Sec6',reference_name='Nuf2',GFP_pattern='w1GFP',median_radius=9)

##Sec10
images_path='/Volumes/MarkoKaksonenLab/Andrea/T1/Andrea/Original_Data/150618_2664_2920/CQ_out_of_CQ/'
msi.experiment(path=images_path,target_name='Sec10',reference_name='Nuf2',GFP_pattern='w1GFP',median_radius=9)

#Sec15
images_path='/Volumes/MarkoKaksonenLab/Andrea/UNIGE/MICROSCOPY_DATA/160805_2920_2665looped_out/Cells/'
msi.experiment(path=images_path,target_name='Sec15',reference_name='Nuf2',GFP_pattern='w1GFP',median_radius=9)

#Exo84
images_path='/Volumes/MarkoKaksonenLab/Andrea/UNIGE/MICROSCOPY_DATA/160805_2920_2667looped_out/Cells/'
msi.experiment(path=images_path,target_name='Exo84',reference_name='Nuf2',GFP_pattern='w1GFP',median_radius=9)
