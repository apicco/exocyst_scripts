#
path = "/Users/andreapicco/Google Drive/RESEARCH/COG_COMPLEX"
folder = "Analysis_SOFTWARE"
source(paste(path,folder,"iCrai_1.0.R",sep="/"))
source(paste(path,folder,"analysis.R",sep="/"))

directory="./"
image.directory="images/"
plot.only=FALSE
#plot.only=TRUE
iCrai(directory,image.directory,plot.only,output.image=FALSE,output.spots=TRUE)

####
load("fit.Rdata")
write.table(matrix(c(round(output.mu[1],2),round(output.mu[2],2),round(output.sigma[1],2),round(output.sigma[2],2),length(selected.dist)),1,5),row.names=FALSE,col.names=FALSE,quote=FALSE,file="output.data")

