#iCRAI/SLM detection software

library(EBImage)
library(cluster)
library(MASS)
library(gregmisc)
library(KernSmooth)

#check for the nearest spot to each spot and output the distance and the ID for the nearest spot
nearest.neigh<-function(input.data){
	nearests<-c()
	for (i in 1:length(input.data[,1])){
		dist<-sqrt((input.data[-i,1]-input.data[i,1])^2+(input.data[-i,2]-input.data[i,2])^2)
		#nearests<-rbind(nearests,c(min(dist),which(dist==min(dist))))
		nearests<-c(nearests,min(dist))
	}
	nearests
}


#contour image and compute distance from contour
otsuseg = function(x, range=c(0, 0.01)) {
  y = as.numeric(x)
  otsuf = function(t, y) {
    z = y>t
    mz = mean(z)
    if (mz==1) sd(y[z])
    else if (mz==0) sd(y[!z])
    else sd(y[z])*mz + sd(y[!z])*(1-mz)
  }
  optimize(otsuf, range, y=y)$minimum
}


get.contour<-function(image.data,input.range=c(0,0.01)){
	t = otsuseg(image.data, range=input.range)
	print(sprintf("ots segmentation: %f",t))
	seg = image.data>t
	
	mk = makeBrush(5,"diamond")
	mk2 = makeBrush(3,"box")
	seg=opening(seg,mk)
	seg=fillHull(seg)
	seg=bwlabel(seg)
	
	seg.erode<-erode(seg,mk2)
	contour<-seg-seg.erode
	#output<-which(contour>0,arr.ind=TRUE)
	#output[,c(2,1)]-1
	list(contour=contour,otsuseg=t)
}

#contour.dist<-function(data,image.data){
contour.dist<-function(input.data,image.contour,input.data2){
	dist<-c()
	angle<-c()
	angle.data<-c()
	tangent<-c()
	contour.coords<-which(image.contour$contour>0,arr.ind=TRUE)
	for (i in seq(1,length(input.data[,1]),by=1)){
		x.dists.tmp<-contour.coords[,2]-input.data[i,1]
		y.dists.tmp<-contour.coords[,1]-input.data[i,2]
		dists.tmp<-sqrt((x.dists.tmp)^2+(y.dists.tmp)^2)
		#store the minimal distance
		min.dist<-min(dists.tmp)
		if (length(min.dist)==1) dist<-c(dist,min.dist)
		else if (length(min.dist)==2 & abs(diff(min.dist))<0.05) dist<-c(dist,min.dist[1])
		else dist<-c(dist,NA)
		
		#compute the inclination (angle) of the tangent to the contour
		#get the closest point on the membrane
		sel.min<-which(dists.tmp==min.dist)
		#get the secon closed point on the membrane, together with the third closest point they will define the tangent at the membrane passing by the closest point	
		min2.dist<-min(dists.tmp[-sel.min])
		if (length(min2.dist)==1){
			sel.min2<-which(dists.tmp==min2.dist)
			min3.dist<-min(dists.tmp[-c(sel.min,sel.min2)])
			sel.min3<-which(dists.tmp==min3.dist)
		} else {
			sel.min2<-which(dists.tmp==min2.dist[1])
			sel.min3<-which(dists.tmp==min2.dist[2])
		}

		dx<-x.dists.tmp[sel.min2[1]]-x.dists.tmp[sel.min3[1]]
		dy<-y.dists.tmp[sel.min2[1]]-y.dists.tmp[sel.min3[1]]
		atan2(dy,dx)->angle.membrane.tmp

		angle<-c(angle,atan2(y.dists.tmp[sel.min[1]],x.dists.tmp[sel.min[1]])-angle.membrane.tmp)
		if(!is.null(input.data2)){
			dx2<-input.data[i,1]-input.data2[i,1]
			dy2<-input.data[i,2]-input.data2[i,2]
			angle.data<-c(angle.data,atan2(dy2,dx2)-angle.membrane.tmp)
		}
		#dist<-c(dist,min(sqrt((contour.coords[,2]-input.data[i,1])^2+(contour.coords[,1]-input.data[i,2])^2)))
	}
	if(!is.null(input.data2)) list(dist=dist,angle=angle,angle.data=angle.data,otsuseg=image.contour$otsuseg)
	else list(dist=dist,angle=angle,otsuseg=image.contour$otsuseg)
}

plot.image<-function(image.data,selected.data=NULL,selected.dataID=NULL,rejected.data=NULL,rejected.dataID=NULL,imageID){

	dims<-dim(image.data)
	image(x=seq(1,dims[1],length.out=nrow(image.data)),y=seq(1,dims[2],length.out=ncol(image.data)),z=image.data[,sort(seq(1,dims[2],length.out=ncol(image.data)),decreasing=TRUE)],col=topo.colors(20),asp=1,main=sprintf("image %02.0f",imageID),xlab="[pixel]",ylab="[pixel]")
	if (length(selected.data)) points(selected.data[,2],dims[2]-(selected.data[,1]),col="black",cex=2,lwd=2)
	if (length(selected.data) & length(selected.dataID)) text(selected.data[,2],dims[2]-(selected.data[,1]),selected.dataID,pos=4,col="black")
	if (length(rejected.data)) points(rejected.data[,2],dims[2]-(rejected.data[,1]),col="red",cex=2,lwd=2)
	if (length(rejected.data) & length(rejected.dataID)) text(rejected.data[,2],dims[2]-(rejected.data[,1]),rejected.dataID,pos=4,col="black")
}

gaussian.fit<-function(cm,image.data,crop.radius=5,weight=TRUE){
	output<-c()
	dims<-dim(image.data)
	for (i in 1:length(cm[,1])){
		fcm<-as.numeric(floor(cm[i,]))[c(2,1)]+2
		sel.x<-max(c(1,fcm[1]-crop.radius)):min(c(dims[1],fcm[1]+crop.radius))
		sel.y<-max(c(1,fcm[2]-crop.radius)):min(c(dims[2],fcm[2]+crop.radius))
		if (length(sel.x)<length(sel.y)){
	 		row.length<-min(length(sel.x),length(sel.y))
			col.length<-max(length(sel.x),length(sel.y))
		} else if (length(sel.x)>length(sel.y)){
			row.length<-max(length(sel.x),length(sel.y))
			col.length<-min(length(sel.x),length(sel.y))
		} else if (length(sel.x)==length(sel.y)){ 
			row.length<-length(sel.x)
			col.length<-length(sel.x)
		}

		coords<-cbind(as.vector(matrix(rep(sel.x,col.length),ncol=col.length,nrow=row.length,byrow=FALSE)),as.vector(matrix(rep(sel.y,row.length),ncol=col.length,nrow=row.length,byrow=TRUE)))
		im.coords<-as.numeric(image.data[coords])
		im.fit<-(im.coords-min(im.coords))/(max(im.coords)-min(im.coords))+0.001

		if (weight) gauss.fit<-lm(I(log(im.fit))~I((coords[,1]-fcm[1])^2)+I((coords[,2]-fcm[2])^2),weight=im.fit)
		else gauss.fit<-lm(I(log(im.fit))~I((coords[,1]-fcm[1])^2)+I((coords[,2]-fcm[2])^2))#,weight=im.fit)
		summary(gauss.fit)->g.summary
		
		output<-rbind(output,c(exp(gauss.fit$coef[1]),1/sqrt(-2*gauss.fit$coef[2]),1/sqrt(-2*gauss.fit$coef[3]),g.summary$r.squared,g.summary$coef[,4]))	
	}
	dimnames(output)<-list(c(),c("I","sigma_x","sigma_y","R^2","Pval I","Pval x","Pval y"))
	output
}


#angular.sel<-function(a,b){
#	a<-a/pi
#	b<-b/pi
#	
#	#focus on angles that have the same sign
#	ab<-cbind(sign(b)*a,sign(a)*b)
#	sel<-which(ab[,1]>0 & ab[,2]>0)
#	ab.QI<-ab[sel,]
#	
#	#select those that are around 3/2pi and bring them to pi/2
#	sel1<-which(ab.QI[,1]>1 & ab.QI[,2]>1)
#	ab.sel<-ab.QI
#	ab.sel[sel1,]<-ab.sel[sel1,]-1
#	
#	#removes the angles which are c(pi/2,3pi/2) or (3pi/2,pi/2) by computing the difference 
#	ab.diff<-abs(ab.sel[,1]-ab.sel[,2])
#	#this should be the x falling in the selected region
#	sel.diff<-ab.diff<1/6
#	
#}
#	
angular.sel<-function(a,b){
	a<-a/pi
	b<-b/pi
	
	#focus on angles that have the same sign
	ab<-cbind(sign(b)*a,sign(a)*b)
	sel<-which(ab[,1]>0 & ab[,2]>0)
	ab.QI<-ab[sel,]
	
	#select those that are around 3/2pi and bring them to pi/2
	sel1<-which(ab.QI[,1]>1 & ab.QI[,2]>1)
	ab.sel<-ab.QI
	ab.sel[sel1,]<-ab.sel[sel1,]-1
	
	#removes the angles which are c(pi/2,3pi/2) or (3pi/2,pi/2) by computing the difference 
	ab.diff<-abs(ab.sel[,1]-ab.sel[,2])
	#this should be the x falling in the selected region
	sel.diff<-ab.diff<1/6
	
}
	
#load the data collection and complete it with the eccentricity, orientation and distance features.
load.data.collection<-function(directory,image.directory=c()){
	collect.W1<-c()
	collect.W2<-c()
	W2<-dir(path=directory,pattern="_W2")
	W1<-dir(path=directory,pattern="_W1_warped")

	if (!is.null(image.directory)) {	
		image.path<-paste(directory,image.directory,sep="")
		images<-dir(image.path,pattern="image_")
		imagesMD<-dir(image.path,pattern="imageMD_")
	}
	for (i in 1:length(W2)){
		print(W1[i])
		print(W2[i])
		read.table(paste(directory,W2[i],sep=""))->tmp.data.W2
		read.table(paste(directory,W1[i],sep=""))->tmp.data.W1

		if (!is.null(image.directory)) {
			#load the image data for the image.id-th image
			image.data<-readImage(paste(image.path,images[i],sep=""))
			image.dataMD<-readImage(paste(image.path,imagesMD[i],sep=""))
			
			#compute the nearest neigh dist and the dist from contour for the image.id-th image	
			image.contour<-get.contour(image.data[,,2])
			nearest.W2<-nearest.neigh(tmp.data.W2)
			contour.W2<-contour.dist(tmp.data.W2,image.contour,tmp.data.W1)
			nearest.W1<-nearest.neigh(tmp.data.W1)
			contour.W1<-contour.dist(tmp.data.W1,image.contour,tmp.data.W2)#note that for the segmentation we use the image on channel W2
			
			gauss.W1<-gaussian.fit(tmp.data.W1,image.dataMD[,,1])
			gauss.W2<-gaussian.fit(tmp.data.W2,image.dataMD[,,2])

		}
			
		#append the data and complete them with the image id
		collect.W1<-rbind(collect.W1,cbind(tmp.data.W1,nearest.W1,contour.W1$dist,contour.W1$angle,-contour.W1$angle.data,rep(i,length(tmp.data.W1[,1])),gauss.W1))
		collect.W2<-rbind(collect.W2,cbind(tmp.data.W2,nearest.W2,contour.W2$dist,contour.W2$angle,contour.W2$angle.data,rep(i,length(tmp.data.W2[,1])),gauss.W2))
	}


	W1.u11 = collect.W1[,10]
	W1.u20 = collect.W1[,11]
	W1.u02 = collect.W1[,12]
	
	#compute eccentricity
	W1.ecc = (W1.u02 + W1.u20 + sqrt( (W1.u20 - W1.u02) * (W1.u20 - W1.u02) + 4 * W1.u11 * W1.u11))/(W1.u02 + W1.u20 - sqrt( (W1.u20 - W1.u02) * (W1.u20 - W1.u02) + 4 * W1.u11 * W1.u11))
	#compute orientation
	W1.ori = 2 * W1.u11 / (W1.u20 - W1.u02)

	
	W2.u11 = collect.W2[,10]
	W2.u20 = collect.W2[,11]
	W2.u02 = collect.W2[,12]

	#compute eccentricity
	W2.ecc = (W2.u02 + W2.u20 + sqrt( (W2.u20 - W2.u02) * (W2.u20 - W2.u02) + 4 * W2.u11 * W2.u11))/(W2.u02 + W2.u20 - sqrt( (W2.u20 - W2.u02) * (W2.u20 - W2.u02) + 4 * W2.u11 * W2.u11))
	#compute orientation
	W2.ori = 2 * W2.u11 / (W2.u20 - W2.u02)

	
	#compute dist
	dist = sqrt((collect.W1[,1]-collect.W2[,1])^2+(collect.W1[,2]-collect.W2[,2])^2)
	
	#compute angular selection to select spots with 90 deg from PM on both channels
	angular.sel(collect.W1[,15],collect.W2[,15])->membrane.angle
	
	collect.W1<-cbind(collect.W1,cbind(W1.ecc,W1.ori,dist,membrane.angle))
	collect.W2<-cbind(collect.W2,cbind(W2.ecc,W2.ori,dist,membrane.angle))

	
	list(W1=collect.W1,W2=collect.W2)
}


p<-function(x,mu,sigma) (x/sigma^2)*exp(-(mu^2+x^2)/(2*sigma^2))*besselI(x*mu/sigma^2,0)

iData<-function(input.data.path,image.path){
	#input.data.path<-"111103_F8/"
	#image.path<-"images/"
	
	contour.dist.cutoff<-22#18
	nearest.dist.cutoff<-10
	mu11<-10
	mu20<-11
	mu02<-12
	contour<-14
	nearest<-13
	eccl<-25

	input.data<-load.data.collection(input.data.path,image.path)
	#trivial sorting 
	#based on distance from membrane (contour)	
	sel.contour.dist<-which(input.data$W1[,contour] < contour.dist.cutoff & input.data$W2[,contour] < contour.dist.cutoff)
	input.data$W1<-input.data$W1[sel.contour.dist,]
	input.data$W2<-input.data$W2[sel.contour.dist,]

	#based on distance from neighbours 
	sel.nearest.dist<-which(input.data$W1[,nearest] > nearest.dist.cutoff & input.data$W2[,nearest] > nearest.dist.cutoff)
	input.data$W1<-input.data$W1[sel.nearest.dist,]
	input.data$W2<-input.data$W2[sel.nearest.dist,]

	data.frame(m2.W1=input.data$W1[,4],m2.W2=input.data$W2[,4],ecc.W1=input.data$W1[,eccl],ecc.W2=input.data$W2[,eccl])->my.df
	
	#d<-dist(my.df.data,method="euclidean")
	#cl<-hclust(d,method="ward")
	#groups<-cutree(cl,2)

	#list(input.data=input.data,group.1=which(groups==1),group.2=which(groups==2),df=df)
	list(input.data=input.data,df=my.df)
}	

gaussian.2D.p<-function(x,y,xp,yp,sigma.x,sigma.y,N){
	exp( - ( x - xp ) * ( x - xp ) / ( 2 * sigma.x ) - ( y - yp ) * ( y - yp ) / ( 2  * sigma.y ) ) / ( 2 * pi * sigma.x * sigma.y * N )
}

cost.function<-function(data){
	#compute a score for each points based on the eccentricity and second moment
}
