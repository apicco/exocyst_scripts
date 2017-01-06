#Quantifications and stats based on the quantification of the fluorescence intensities
quantify_number_of_molecules=function(path,target,reference,output,n_reference,target_name,reference_name,log_reference=TRUE,log_target=TRUE,...){

	t<-read.table(paste(path,target,sep="/"))[,1]
	r<-read.table(paste(path,reference,sep="/"))[,1]

	#remove the zeros, which mark wrong patches
	t<-t[which(t!=0)]
	r<-r[which(r!=0)]

	if (log_reference) {
		r<-log(r)
		r<-r[which(r>8)]
		xlab_reference=expression(log(Fi))
	} else 	xlab_reference=expression(Fi (a.u.))
	if (log_target) {
		t<-log(t)
		xlab_target=expression(log(Fi))
		xlab_ss=expression(log(Size))
	} else {	
		xlab_target=expression(Fi (a.u.))
		xlab_ss=expression(Size (px))
	}


	mt<-c(median(t),mad(t)/sqrt(length(t)))
	mr<-c(median(r),mad(r)/sqrt(length(r)))

	if (log_reference) reference=c(exp(-mr[1]),exp(-mr[1])*mr[2])
	else reference=c(1/mr[1],mr[2]/mr[1]^2)

	if (log_target) target=c(exp(mt[1]),exp(mt[1])*mt[2])
	else target=c(mt[1],mt[2])
	nr<-c(
		n_reference[1]*target[1]*reference[1],
		sqrt(
			(n_reference[2]*target[1]*reference[1])^2 +
			(n_reference[1]*target[2]*reference[1])^2 +
			(n_reference[1]*target[1]*reference[2])^2
		)
	)
			
	pdf(paste(path,paste(output,'hist.pdf',sep="_"),sep="/"),w=5,h=10)
	par(mfcol=c(2,1))
	par(cex=1.3)
	hist(t,xlab=xlab_target,main=sprintf('%s: %0.0f±%0.0f molec. (n=%d)',target_name,nr[1],nr[2],length(t)),col="#CDCDCD",...)
	if (log_reference) hist(r,xlab=xlab_reference,main=sprintf('%s log(Fi): %0.2f±%0.2f (n=%d)',reference_name,mr[1],mr[2],length(r)),col="#CDCDCD",...)
	else hist(r,xlab=xlab_reference,main=sprintf('%s Fi: %0.0f±%0.0f a.u. (n=%d)',reference_name,mr[1],mr[2],length(r)),col="#CDCDCD",...)

	dev.off()
	
	return(nr)
}
