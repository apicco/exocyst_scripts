#
library(bbmle)

eccentricity<-function(mu11,mu20,mu02){
	D<-(mu20-mu02)^2+4*mu11^2
	(mu20+mu02+sqrt(D))/(mu20+mu02-sqrt(D))
}
spot.plot<-function(x,y,image.data,imgID,crop.radius=5,input.mfrow=c(2,2),selected=TRUE){
	dims<-dim(image.data)	
	par(mfrow=input.mfrow)
	ctrl<-0	
	for (i in 1:length(x[,1])){
		
		#load image crop around x
		float.fx<-as.numeric(x[i,])[c(2,1)]+2
		fx<-as.numeric(floor(x[i,]))[c(2,1)]+2
		sel.x1<-max(c(1,fx[1]-crop.radius)):min(c(dims[1],fx[1]+crop.radius))
		sel.x2<-max(c(1,fx[2]-crop.radius)):min(c(dims[2],fx[2]+crop.radius))
		im.x<-image.data[sel.x1,sel.x2,1][,sort(seq(1:length(sel.x2)),decreasing=TRUE)]
		#load image crop around y
		float.fy<-as.numeric(y[i,])[c(2,1)]+2
		fy<-as.numeric(floor(y[i,]))[c(2,1)]+2
		sel.y1<-max(c(1,fy[1]-crop.radius)):min(c(dims[1],fy[1]+crop.radius))
		sel.y2<-max(c(1,fy[2]-crop.radius)):min(c(dims[2],fy[2]+crop.radius))
		im.y<-image.data[sel.y1,sel.y2,2][,sort(seq(1:length(sel.y2)),decreasing=TRUE)]

		if (selected){
			image(x=sel.x1,y=sel.x2,z=im.x,asp=1,main=sprintf("image %0.0f, sel W1 spot %0.0f",imgID,i),xlab="[pixel]",ylab="[pixel]")
			image(x=sel.y1,y=sel.y2,z=im.y,asp=1,main=sprintf("image %0.0f, sel W2 spot %0.0f",imgID,i),xlab="[pixel]",ylab="[pixel]")
		} else {
			image(x=sel.x1,y=sel.x2,z=im.x,asp=1,main=sprintf("image %0.0f, rej W1 spot %0.0f",imgID,i),xlab="[pixel]",ylab="[pixel]")
			image(x=sel.y1,y=sel.y2,z=im.y,asp=1,main=sprintf("image %0.0f, rej W2 spot %0.0f",imgID,i),xlab="[pixel]",ylab="[pixel]")
		}
		
		ctrl<-ctrl+2
		if (ctrl == input.mfrow[1]*input.mfrow[2]) ctrl<-0
	}
	if (ctrl > 0){
		for (i in 1:(input.mfrow[1]*input.mfrow[2]-ctrl)) plot.new()
	}
}

odca2<-function(x,mu=20,sigma=20,input.filter=1){

	x.old<-x
	fit<-c()
	scores<-c()
	fits<-list()
	it<-0
	
	LL<-function(mu,sigma,data) -sum(log(p(data,mu,sigma)))
	
	cutoff<-2*length(x.old)/5
	#if (length(x.old)/3 > 20) cutoff<-length(x.old)/3
	#else cutoff<-20
	#if (cutoff > length(x.old)) cutoff=floor(length(x.old)/2)

	if ((length(x.old)-cutoff)<=input.filter) input.filter=1

	print(sprintf("starting number of distances: %0.0f; cutoff: %0.0f",length(x.old),cutoff))
	while(length(x)>cutoff){
		l<-c() #clear l
		it<-it+1
		fits[[it]]<-mle2(LL,start=list(mu=mu,sigma=sigma),data=list(data=x))
		for (i in 1:length(x)) l<-c(l,LL(fits[[it]]@coef[[1]],fits[[it]]@coef[[2]],x[-i])) #log-likelihood without the i-th data given the mu and sigma of fits
		m<-min(l,na.rm=T)
		scores<-rbind(scores,c(which(x.old==x[which(l==m)]),LL(fits[[it]]@coef[[1]],fits[[it]]@coef[[2]],x),m,summary(fits[[it]])@coef[1,c(1,2)],summary(fits[[it]])@coef[2,c(1,2)])) #store the x.old that should be removed, the max likelihood without and with rejectig one possible outlier, and the mu and sigma 
		#scores<-rbind(scores,c(which(x.old==x[which(l==m)]),LL(fits[[it]]@coef[[1]],fits[[it]]@coef[[2]],x),m,fits[[it]]@coef[[1]],fits[[it]]@coef[[2]])) #store the x.old that should be removed, the max likelihood without and with rejectig one possible outlier, and the mu and sigma 
		x<-x[-which(l==m)]
	}
#	#METHOD #1
#	#I want to minimize the improvement in the m the outlier removing
#	#for each dataset X an element x is removed: X2 = X1 - x_i and the fit computed. deltad	is the difference in the fit between the two dataset: fit_2 - fit_1. If the removal of x_i was not so important the difference between the fit is minimal, therefore X1 was already a good estimate
#	delta.d<-cbind(scores[-length(scores[,1]),4],log(abs(scores[-length(scores[,1]),4]-scores[-1,4])))
#	#ONLY the data with fits greater than the fit with the maximal delta.d will be taken into account, in fact it can be that outliers give a min to delta.d too (if there are too many outliers for instance)
#	max.delta.d<-delta.d[which(delta.d[,2]==max(delta.d[,2],na.rm=T)),]
#	sel.delta.d<-delta.d[which(delta.d[,1]>=max.delta.d[1]),]
#	#select the min deltat and the sel out
#	min.delta.d<-sel.delta.d[which(sel.delta.d[,2]==min(sel.delta.d[,2],na.rm=T))]
#	sel.out<-which(scores[,4]==min.delta.d[1])
#	#the data will be those forming X1, which means all the x.old less those remuved up to X1 excluded (hence sel.out-1)
#	return(list(data=x.old[-scores[1:(sel.out-1),1]],LL=scores[sel.out,2],fit=summary(fits[[sel.out]]),delta.d=delta.d,sel.delta.d=sel.delta.d,sel=sel.out))
#	#METHOD #2
#	#I want to minimize the improvement in the m the outlier removing therefore I select the mu and sigma corresponding to the x[i] whose removal did improve the likelihood of the dataset the less.
#	#I compute the ratio between the likelihood with and without x[i]. this will be close greathere or = to one
#	Lratio<-as.vector(filter(scores[,2]/scores[,3],rep(1/input.filter,input.filter)))
#
#	#as removing a point is likely to decrease the likelihood: Lratio should be greater of equal to one
#
#	sel.out=which(Lratio==min(Lratio,na.rm=T))[1]
#
#	#the data will be those forming X1, which means all the x.old less those remuved up to X1 excluded (hence sel.out-1)
#	return(list(data=x.old[-scores[1:(sel.out-1),1]],LL=scores[sel.out,2],fit=summary(fits[[sel.out]]),delta.d=cbind(scores[,4],Lratio),sel.delta.d=1:length(Lratio),sel=sel.out))
#	#METHOD #3
#	#I want to select the dataset whose removal of each of his datapoints give the same likelihood. In other words all the points in the dataset are eligible to be there. The parameter to take into account is the standard deviation of the likelihoods compuded removing a point at the time. I want to minimize this likelihood.
#
#	sel.out=which(scores[,6]==min(scores[,6]))
#	return(list(data=x.old[-scores[1:(sel.out-1),1]],LL=scores[sel.out,2],fit=summary(fits[[sel.out]]),delta.d=cbind(scores[,4],scores[,6]),sel.delta.d=1:length(scores[,6]),sel=sel.out))
#	#METHOD #4
#	#When most of the outliers get rejected the distances estimates do not differ too much. First I identify the region in which distances accumulate, then I search for a local minima for the ratio between the likelihoods in this neighborhood
#	delta.d<-cbind(scores[-length(scores[,1]),4],log(abs(scores[-length(scores[,1]),4]-scores[-1,4])))
#	#ONLY the data with fits greater than the fit with the maximal delta.d will be taken into account, in fact it can be that outliers give a min to delta.d too (if there are too many outliers for instance)
#	max.delta.d<-delta.d[which(delta.d[,2]==max(delta.d[,2],na.rm=T)),]
#	sel.d<-which(delta.d[,1]>=max.delta.d[1])+1#the last distance was missing and the max was included with the >=, +1 exclde the max and include the last distance
#
#	h<-hist(scores[,4])
#	sel.break<-which(h$counts==max(h$counts,na.rm=T))+1#that's the starting distance vaule from which to look for the local minima in the ratio between likelihoods
#	look_for_the_minimum=TRUE
#	Lratio<-cbind(scores[,4],as.vector(filter(scores[,2]/scores[,3],rep(1/input.filter,input.filter))))
#	j=1
#	jj=0
#	while(look_for_the_minimum){
#		L<-Lratio[which(Lratio[,1]>h$breaks[sel.break-j] & Lratio[,1]<h$breaks[sel.break+jj]),]
#		L<-L[order(L[,1]),]
#		sel.min.L<-which(L[,2]==min(L[,2],na.rm=TRUE))
#		
#		#the minimum must not be in the boundaries of the selection	
#		if (sel.min.L > 1 & sel.min.L < length(L[,1])) {
#			look_for_the_minimum=FALSE
#		} else if (sel.min.L == length(L[,1])){
#			if ((sel.break+jj)<length(h$breaks)) jj=jj+1
#			else look_for_the_minimum=FALSE
#		} else if (sel.min.L ==1){
#			if ((sel.break-j)>1) j=j+1
#			else look_for_the_minimum=FALSE
#		}
#	}
#
#	sel.out=which(scores[,4]==L[sel.min.L,1])[1]
#	return(list(data=x.old[-scores[1:(sel.out-1),1]],LL=scores[sel.out,2],fit=summary(fits[[sel.out]]),delta.d=Lratio,sel.delta.d=1:length(scores[,1]),sel=sel.out))
	#METHOD 5

	pd.tmp<-1/abs(scores[-length(scores[,4]),4]-scores[-1,4])
	pd<-pd.tmp/sum(pd.tmp,na.rm=TRUE)
	
	pS.tmp<-1/abs(scores[-length(scores[,6]),6]-scores[-1,6])
	pS<-pS.tmp/sum(pS.tmp,na.rm=TRUE)

	S<--pd*log(pd)-pS*log(pS)
	#clear values of S that correspont to NA errors:
	scores[-length(scores[,1]),5]*scores[-length(scores[,1]),6]->nas
	nas[which(!is.na(nas))]<-1
	S<-S*nas
	if (length(scores[,1])>=10) delta.d<-cbind(scores[-length(scores[,1]),4],filter(log(abs(scores[-length(scores[,1]),4]-scores[-1,4])),rep(1/5,5)))
	else delta.d<-cbind(scores[-length(scores[,1]),4],filter(log(abs(scores[-length(scores[,1]),4]-scores[-1,4])),rep(1/floor(length(scores[,1])/2),floor(length(scores[,1])/2))))
	#ONLY the data with fits greater than the fit with the maximal delta.d will be taken into account, in fact it can be that outliers give a min to delta.d too (if there are too many outliers for instance)
	max.delta.d<-delta.d[which(delta.d[,2]==max(delta.d[,2],na.rm=T)),]
	input_min_d<-scan("parameters.data",comment.char="#")[1]
	if(is.na(input_min_d)) sel.d<-which(delta.d[,1]>=max.delta.d[1])
	else sel.d<-which(delta.d[,1]>=input_min_d)
	sel.out<-which(S==max(S[sel.d],na.rm=TRUE))[1]
#	save(scores,S,delta.d,max.delta.d,nas,sel.d,sel.out,x.old,cutoff,fits,file="tmp.Rdata")

	if (!is.na(sel.out) & sel.out>1) return(list(data=x.old[-scores[1:(sel.out-1),1]],LL=scores[sel.out,2],fit=summary(fits[[sel.out]]),delta.d=cbind(scores[-length(scores[,1]),4],S),sel.delta.d=1:length(S),sel=sel.out))
	else if (!is.na(sel.out) & sel.out==1) return(list(data=x.old,LL=scores[sel.out,2],fit=summary(fits[[sel.out]]),delta.d=cbind(scores[-length(scores[,1]),4],S),sel.delta.d=1:length(S),sel=sel.out))
	else return(list(data=x.old,LL=scores[,2],fit=summary(fits[[1]]),delta.d=cbind(scores[-length(scores[,1]),4],rep(0,length(S))),sel.delta.d=1:length(S),sel=sel.out))
	
}


iCrai<-function(directory="111103_F8/",image.directory="images/",plot.only=FALSE,sort.angle=TRUE,output.image=TRUE,output.spots=TRUE){
	percent.cutoff=0.5
	my.gridsize<-c(100,100)	
	#load the data, if already estrapolated from Traj_* otherwhise perform the extrapolation and clustering
	if (plot.only) {
		load(paste(directory,"output.Rdata",sep=""))
	} else {
		out<-iData(input.data.path=directory,image.path=image.directory)
		save(out,file=paste(directory,"output.Rdata",sep=""))
	}
	#shorter nomenclature	
	W1.data<-out$input.data$W1
	W2.data<-out$input.data$W2
	df.data<-out$df

	#select spots sharing most common features
	a.W1<-bkde2D(cbind(df.data[,"m2.W1"],df.data[,"ecc.W1"]),c(dpik(df.data[,"m2.W1"]),dpik(df.data[,"ecc.W1"])),gridsize=my.gridsize)
	a.W2<-bkde2D(cbind(df.data[,"m2.W2"],df.data[,"ecc.W2"]),c(dpik(df.data[,"m2.W2"]),dpik(df.data[,"ecc.W2"])),gridsize=my.gridsize)


	dx1<-a.W1$x1[2]-a.W1$x1[1]
	dx2<-a.W1$x2[2]-a.W1$x2[1]
	a.sort.W1<-sort(a.W1$fhat,decreasing=TRUE)
	a.sort.W2<-sort(a.W2$fhat,decreasing=TRUE)
	sum.sort.W1<-c(0)
	sum.sort.W2<-c(0)
	for ( i in 1:(length(a.sort.W1))) {
		sum.sort.W1<-sum.sort.W1+a.sort.W1[i]
		if (sum.sort.W1*dx1*dx2>=percent.cutoff) break
	}
	treshold.W1<-a.sort.W1[i-1]
	for ( i in 1:(length(a.sort.W2))) {
		sum.sort.W2<-sum.sort.W2+a.sort.W2[i]
		if (sum.sort.W2*dx1*dx2>=percent.cutoff) break
	}
	treshold.W2<-a.sort.W2[i-1]

	selected.region.W1<-cbind(a.W1$x1[which(a.W1$fhat>=treshold.W1,arr.ind=TRUE)[,1]],a.W1$x2[which(a.W1$fhat>=treshold.W1,arr.ind=TRUE)[,2]])
	no.selected.region.W1<-cbind(a.W1$x1[which(a.W1$fhat<treshold.W1,arr.ind=TRUE)[,1]],a.W1$x2[which(a.W1$fhat<treshold.W1,arr.ind=TRUE)[,2]])
	selected.region.W2<-cbind(a.W2$x1[which(a.W2$fhat>=treshold.W2,arr.ind=TRUE)[,1]],a.W2$x2[which(a.W2$fhat>=treshold.W2,arr.ind=TRUE)[,2]])
	no.selected.region.W2<-cbind(a.W2$x1[which(a.W2$fhat<treshold.W2,arr.ind=TRUE)[,1]],a.W2$x2[which(a.W2$fhat<treshold.W2,arr.ind=TRUE)[,2]])
	selected.df.data<-c()
	for (i in 1:length(df.data[,1])){
		tmp1.W1<-min(sqrt((df.data[i,"m2.W1"]-selected.region.W1[,1])^2+(df.data[i,"ecc.W1"]-selected.region.W1[,2])^2))
		tmp2.W1<-min(sqrt((df.data[i,"m2.W1"]-no.selected.region.W1[,1])^2+(df.data[i,"ecc.W1"]-no.selected.region.W1[,2])^2))
		tmp1.W2<-min(sqrt((df.data[i,"m2.W2"]-selected.region.W2[,1])^2+(df.data[i,"ecc.W2"]-selected.region.W2[,2])^2))
		tmp2.W2<-min(sqrt((df.data[i,"m2.W2"]-no.selected.region.W2[,1])^2+(df.data[i,"ecc.W2"]-no.selected.region.W2[,2])^2))
		if (tmp1.W1<tmp2.W1 & tmp1.W2<tmp2.W2) selected.df.data<-c(selected.df.data,i)
	}

	d<-dist(df.data,method="euclidean")
	cl<-hclust(d,method="ward")
	groups<-cutree(cl,2)

	#PAM HCLUSTER? or just KERNEL?
	
	x<-W1.data[selected.df.data,]
	y<-W2.data[selected.df.data,]
	
	#distance
	d<-x[,"dist"]*64.5
	
	#select non bad gaussian fitting (note that gaussian fitting could to be improved):
	p.cutoff<-0.35
	which(x[,"R^2"] > p.cutoff & y[,"R^2"] > p.cutoff)->p.sel
	
	#once selected the spots check angles with membrane and reject spots whose angle differ too much (this additional selection step might not be required)
	a<-cbind(abs(x[p.sel,15]/pi),abs(y[p.sel,15]/pi))
	a[which(a>1)]<-a[which(a>1)]-1
	
	a.fit<-lm(a[,2]~a[,1]+0)
	a.stats<-summary(a.fit)
	which(as.numeric(abs(a.stats$residuals))<0.05)->s
	#combine both angle sel and gaussian sel
	if (sort.angle) my.sel<-p.sel[s]
	else my.sel<-p.sel
	#--------------------------------------#
	#iterative optim (odca)  to remove outlayers#
	#--------------------------------------#
	selected.dist.with.outliers<-sort(d[my.sel])
	print("running odca")
	input_mu_and_sigma<-scan("parameters.data",comment.char="#")[c(2,3)]
	if (is.na(input_mu_and_sigma[1]) || is.na(input_mu_and_sigma[2])) o<-odca2(selected.dist.with.outliers)
	else o<-odca2(selected.dist.with.outliers,mu=input_mu_and_sigma[1],sigma=input_mu_and_sigma[2])
	print("odca done")
	selected.dist<-o$data
	selected.x<-x[my.sel,][pmatch(selected.dist,selected.dist.with.outliers),]
	selected.y<-y[my.sel,][pmatch(selected.dist,selected.dist.with.outliers),]
	rejected.x<-x[my.sel,][-pmatch(selected.dist,selected.dist.with.outliers),]
	rejected.y<-y[my.sel,][-pmatch(selected.dist,selected.dist.with.outliers),]
	
	#outputs
	output.mu<-c(o$fit@coef[[1,1]],o$fit@coef[[1,2]])
	output.sigma<-c(o$fit@coef[[2,1]],o$fit@coef[[2,2]])
	if (is.na(o$sel)) {
		output.mu*NA->output.mu
		output.sigma*NA->output.sigma
	}
	save(output.mu,output.sigma,selected.dist,selected.dist.with.outliers,file=paste(directory,"fit.Rdata",sep=""))
	pdf(paste(directory,"output.pdf",sep=""),useDingbats=FALSE)
	#output analysis
	par(mfcol=c(2,2))
	#pag1	
	plot(W1.data[,4],W1.data[,"W1.ecc"],xlab=quote(paste("W1 ",m[2])),ylab=" W1 ecc",pch=20)
	points(x[,4],x[,"W1.ecc"],pch=20,col="red")
	plot(W2.data[,4],W2.data[,"W2.ecc"],xlab=quote(paste("W2 ",m[2])),ylab="W2 ecc",pch=20)
	points(y[,4],y[,"W2.ecc"],pch=20,col="red")
	plot(W2.data[,"W2.ecc"],W1.data[,"W1.ecc"],xlab="W2 ecc",ylab="W1 ecc",pch=20)
	points(y[,"W2.ecc"],x[,"W1.ecc"],pch=20,col="red")
	plot(W2.data[,4],W1.data[,4],xlab=quote(paste("W2 ",m[2])),ylab=quote(paste("W1 ",m[2])),pch=20)
	points(y[,4],x[,4],pch=20,col="red")

	#pag2	
	plot(x[,"R^2"],y[,"R^2"],pch=20,xlab=quote(R[W1]^2),ylab=quote(R[W2]^2))
	points(x[p.sel,"R^2"],y[p.sel,"R^2"],pch=20,col="red")
	plot(a,xlab="W1 angle",ylab="W2 angle",pch=20)
	if (sort.angle) points(a[s,1],a[s,2],col="red",pch=20)
	else points(a[s,1],a[s,2],col="black",pch=20)
	plot.new()
	plot(seq(1,length(a.stats$residuals)),a.stats$residuals,ylab="angle diff.",xlab="spot pair",pch=20)
	if (sort.angle) points(seq(1,length(a.stats$residuals))[s],a.stats$residuals[s],pch=20,col="red")
	else points(seq(1,length(a.stats$residuals))[s],a.stats$residuals[s],pch=20,col="black")

	#pag3	
	plot(1:length(selected.dist.with.outliers),selected.dist.with.outliers,pch=16,xlab="",ylab="[nm]",col="grey",xaxt="n")
	points((1:length(selected.dist.with.outliers))[pmatch(selected.dist,selected.dist.with.outliers)],selected.dist.with.outliers[pmatch(selected.dist,selected.dist.with.outliers)],col="red")
	plot(o$delta.d,pch=16,ylab="S",xlab="d [nm]")
	points(o$sel.delta.d,col="red",lwd=2)
	points(o$delta.d[o$sel,][1],o$delta.d[o$sel,][2],col="green",pch=16,cex=2)
	the.distr<-p(sort(selected.dist),output.mu[1],output.sigma[1])
	hist(selected.dist,f=F,col="red",breaks="Freedman-Diaconis",xlim=c(0,max(hist(selected.dist.with.outliers,plot=F,breaks="Freedman-Diaconis")$breaks)),ylim=c(0,max(c(the.distr,hist(selected.dist.with.outliers,plot=F,breaks="Freedman-Diaconis")$density,hist(selected.dist,plot=F,breaks="Freedman-Diaconis")$density),na.rm=T)),xlab="dist. [nm]",main="distances")
	hist(selected.dist.with.outliers,add=T,f=F,col="grey",breaks="Freedman-Diaconis")
	if (!is.na(o$sel)) {
		hist(selected.dist,f=F,col="red",breaks="Freedman-Diaconis",xlab="dist. [nm]",main=bquote(paste(mu == .(round(output.mu[1],2)) %+-% .(round(output.mu[2],2))," nm",";  ",sigma == .(round(output.sigma[1],2)) %+-% .(round(output.sigma[2],2))," nm; ","n=",.(length(selected.dist)))),cex.main=1,ylim=c(0,max(c(the.distr,hist(selected.dist.with.outliers,plot=F,breaks="Freedman-Diaconis")$density,hist(selected.dist,plot=F,breaks="Freedman-Diaconis")$density),na.rm=T)))
		lines(sort(selected.dist),the.distr,col="black",lwd=3)
	} else plot.new()
	dev.off()
	#output image and spots (selected and rejected among selected)
	image.path<-paste(directory,image.directory,sep="")
	images<-paste(image.path,dir(image.path,pattern="imageMD_"),sep="")
	if (output.image){
		print("output image")
		for(imageID in 1:length(images)){
			print(imageID)
			png(paste(image.path,sprintf("image%02.0f.png",imageID),sep=""),1000,800)
			par(mfcol=c(1,1))
			image.data<-5*readImage(images[imageID])
			#plot image
			selected.in.img<-which(selected.y[,17]==imageID)
			rejected.in.img<-which(rejected.y[,17]==imageID)
			if (length(selected.in.img)) { 
				tmp.sel<-selected.y[selected.in.img,]
				tmp.selID<-seq(1,length(selected.in.img))
			} else {
				tmp.sel<-NULL
				tmp.selID<-NULL
			}
			if (length(rejected.in.img)) { 
				tmp.rej<-rejected.y[rejected.in.img,]
				tmp.rejID<-seq(1,length(rejected.in.img))
			} else {
				tmp.rej<-NULL
				tmp.rejID<-NULL
			}
			plot.image(image.data[,,2],selected.data=tmp.sel,selected.dataID=tmp.selID,rejected.data=tmp.rej,rejected.dataID=tmp.rejID,imageID=imageID)
			dev.off()
		}
	}
	if (output.spots){	
		#output spots
		print("output spots")
		pdf(paste(directory,sprintf("spots.pdf"),sep=""),useDingbats=FALSE)
		for(imageID in 1:length(images)){
			print(imageID)
			image.data<-5*readImage(images[imageID])
			#plot image
			selected.in.img<-which(selected.y[,17]==imageID)
			rejected.in.img<-which(rejected.y[,17]==imageID)
			if (length(selected.in.img)) spot.plot(selected.x[selected.in.img,],selected.y[selected.in.img,],image.data,imageID)
			if (length(rejected.in.img)) spot.plot(rejected.x[rejected.in.img,],rejected.y[rejected.in.img,],image.data,imageID,selected=FALSE)
		}
		dev.off()
	}
}
