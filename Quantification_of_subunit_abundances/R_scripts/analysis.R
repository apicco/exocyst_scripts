#
source("Doc/quantification.R")

n_nuf2=c(280.6,16.1)

n_Sec5=quantify_number_of_molecules(
	path='/Volumes/MarkoKaksonenLab/Andrea/T1/Andrea/Original_Data/150617_2661_2920/CQ_out_of_CQ/',
	target='Sec5_intensities.txt',
	reference='Nuf2_intensities.txt',
	output='Sec5-GFP_Nuf2-GFP',
	n_reference=n_nuf2,
	target_name='Sec5',
	reference_name='Nuf2',
	log_reference=TRUE,
	log_target=TRUE,
	breaks='Freedman-Diaconis'
)


n_Sec6=quantify_number_of_molecules(
	path='/Volumes/MarkoKaksonenLab/Andrea/UNIGE/MICROSCOPY_DATA/160804_2920_2662looped_out/Cells/',
	target='Sec6_intensities.txt',
	reference='Nuf2_intensities.txt',
	output='Sec6-GFP_Nuf2-GFP',
	n_reference=n_nuf2,
	target_name='Sec6',
	reference_name='Nuf2',
	log_reference=TRUE,
	log_target=TRUE,
	breaks='Freedman-Diaconis'
)


n_Sec10=quantify_number_of_molecules(
	path='/Volumes/MarkoKaksonenLab/Andrea/T1/Andrea/Original_Data/150618_2664_2920/CQ_out_of_CQ/',
	target='Sec10_intensities.txt',
	reference='Nuf2_intensities.txt',
	output='Sec10-GFP_Nuf2-GFP',
	n_reference=n_nuf2,
	target_name='Sec10',
	reference_name='Nuf2',
	log_reference=TRUE,
	log_target=TRUE,
	breaks='Freedman-Diaconis'
)


n_Sec15=quantify_number_of_molecules(
	path='/Volumes/MarkoKaksonenLab/Andrea/UNIGE/MICROSCOPY_DATA/160805_2920_2665looped_out/Cells/',
	target='Sec15_intensities.txt',
	reference='Nuf2_intensities.txt',
	output='Sec15-GFP_Nuf2-GFP',
	n_reference=n_nuf2,
	target_name='Sec15',
	reference_name='Nuf2',
	log_reference=TRUE,
	log_target=TRUE,
	breaks='Freedman-Diaconis'
)

n_Exo84=quantify_number_of_molecules(
	path='/Volumes/MarkoKaksonenLab/Andrea/UNIGE/MICROSCOPY_DATA/160805_2920_2667looped_out/Cells/',
	target='Exo84_intensities.txt',
	reference='Nuf2_intensities.txt',
	output='Exo84-GFP_Nuf2-GFP',
	n_reference=n_nuf2,
	target_name='Exo84',
	reference_name='Nuf2',
	log_reference=TRUE,
	log_target=TRUE,
	breaks='Freedman-Diaconis'
)


print("--------------------------------")
print("QUANTIFICATIONS")
n_Sec5 = round(n_Sec5)
n_Sec6 = round(n_Sec6)
n_Sec10 = round(n_Sec10)
n_Sec15 = round(n_Sec15)
n_Exo84 = round(n_Exo84)

print(sprintf("Sec5: %d +- %d",n_Sec5[1],n_Sec5[2]))
print(sprintf("Sec6: %d +- %d",n_Sec6[1],n_Sec6[2]))
print(sprintf("Sec10: %d +- %d",n_Sec10[1],n_Sec10[2]))
print(sprintf("Sec15: %d +- %d",n_Sec15[1],n_Sec15[2]))
print(sprintf("Exo84: %d +- %d",n_Exo84[1],n_Exo84[2]))


pdf("quantifications.pdf")
par(omi=c(0,0,0,0))
par(cex.lab=2.5)
par(cex.axis=2)
par(cex.sub=1)
par(bty="n")
par(mar=c(9,10,2,1))
par(mgp=c(5,2,0))
par(xpd=TRUE)

br<-barplot(
	c( 
	n_Sec5[1],
	n_Sec6[1], 
	n_Sec10[1],
	n_Sec15[1],
	n_Exo84[1]
	),ylim=c(0,18),ylab="Number of molecules")
arrows(br, 
	c(
	n_Sec5[1]+n_Sec5[2],
	n_Sec6[1]+n_Sec6[2],
	n_Sec10[1]+n_Sec10[2],
	n_Sec15[1]+n_Sec15[2],
	n_Exo84[1]+n_Exo84[2]
	), br , 
	c(
	n_Sec5[1]-n_Sec5[2],
	n_Sec6[1]-n_Sec6[2],
	n_Sec10[1]-n_Sec10[2],
	n_Sec15[1]-n_Sec15[2],
	n_Exo84[1]-n_Exo84[2]
	),angle=90,cod=3)
axis(1,br,label=c( "Sec5" , "Sec6" , "Sec10" , "Sec15" , "Exo84" ) , las = 3)
dev.off()

