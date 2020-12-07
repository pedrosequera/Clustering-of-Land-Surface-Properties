##########READING DATA
#########################
##########################
require(tiff)

path="/Users/pedrosequera/Desktop/Pe/CCNY-CUNY/CUERG/Thesis/LCLU_Thesis/Urban_Subsets/Urban Only/"

####Reading the Data
albedo=readTIFF(paste(path,"albedo_1km_urbanonly.tif",sep=""))
emissivity=readTIFF(paste(path,"emissivity_1km_urbanonly.tif",sep=""))
LST=readTIFF(paste(path,"LST_1km_urbanonly.tif",sep=""))



LST=c(LST)
albedo=c(albedo)
emissivity=c(emissivity)

LST[LST==0]=NA
albedo[albedo==0]=NA
emissivity[emissivity==0]=NA





require(ggplot2)

downtown=data.frame(cbind(albedo,emissivity,LST))
####Eliminating NAs
downtown=downtown[complete.cases(downtown),]


#############################
############################# K-Medoids Clustering
library(fpc)
####Perform clustering without specifying number of groups
######pamk.result = pamk(downtown)
#####the dataset is too big for pamk. Let's use function "clara"


####For 6 groups
clara.result = clara(downtown[,1:2],k=7,samples=1000,sampsize=200)
clara.result$medoids
    # albedo emissivity
# 28992 0.17840475      0.967
# 2896  0.22813261      0.960
# 31630 0.15444693      0.933
# 45154 0.12466139      0.947
# 49800 0.15072021      0.971
# 40282 0.12388182      0.972
# 44932 0.09323987      0.968

####Plotting the 7 groups and silhouette
plot(albedo,emissivity,cex=0.1,col=clara.result$clustering)

png(file=paste(path,"Silhouette.png",sep=""))
plot(silhouette(clara.result),main="Silhouette plot for 7 classes")
dev.off()

boxplot(split(downtown$albedo,clara.result$clustering),ylab="albedo",xlab="categories",main="Categorization of Urban Albedos")
boxplot(split(downtown$emissivity,clara.result$clustering),ylab="emissivity",xlab="categories",main="Categorization of Urban Emissivities")









####Statistics of albedos by classes

clust.stats=matrix(NA,nrow=4,ncol=7)

for(i in 1:7){
	clust.stats[1,i]=min(downtown$albedo[clara.result$clustering==i])
	clust.stats[2,i]=mean(downtown$albedo[clara.result$clustering==i])
	clust.stats[3,i]=median(downtown$albedo[clara.result$clustering==i])
	clust.stats[4,i]=max(downtown$albedo[clara.result$clustering==i])
	rownames(clust.stats)=c("min","mean","median","max")
	colnames(clust.stats)=c("Class 1","Class 2","Class 3","Class 4","Class 5","Class 6","Class 7")
	stats=round(clust.stats,digits=2)
}

write.table(data.frame(stats),file=paste(path,"stats_albedo_clust.txt"),sep=" ")

clust.stats=matrix(NA,nrow=4,ncol=7)

for(i in 1:7){
	clust.stats[1,i]=min(downtown$emissivity[clara.result$clustering==i])
	clust.stats[2,i]=mean(downtown$emissivity[clara.result$clustering==i])
	clust.stats[3,i]=median(downtown$emissivity[clara.result$clustering==i])
	clust.stats[4,i]=max(downtown$emissivity[clara.result$clustering==i])
	rownames(clust.stats)=c("min","mean","median","max")
	colnames(clust.stats)=c("Class 1","Class 2","Class 3","Class 4","Class 5","Class 6","Class 7")
	stats=round(clust.stats,digits=2)
}

write.table(data.frame(stats),file=paste(path,"stats_emissivity_clust.txt"),sep=" ")



##########PLOTTING ALBEDO CLASSES vs EMISSIVITY and LST



png(file=paste(path,"Urban Albedo Classes vs Emissivity.png",sep=""))
plot(downtown$albedo,downtown$emissivity,cex=0.35,col=clara.result$clustering,main="Urban Classes Albedo vs Emissivity",xlab="Albedo",ylab="Emissivity",ylim=c(0.80,max(downtown$emissivity)))
legend(x="bottomright",legend=as.character(1:7),pch=rep(19,7),col=1:7,pt.bg=1:7,bg="gray",title="Classes")
dev.off()



png(file=paste(path,"Urban Albedo Classes vs LST.png",sep=""))
plot(downtown$albedo,downtown$LST,cex=0.35,col=clara.result$clustering,main="Urban Classes Albedo vs LST",xlab="Albedo",ylab=expression(paste("LST (",degree,"C)")))
legend(x="bottomright",legend=as.character(1:7),pch=rep(19,7),col=1:7,pt.bg=1:7,bg="gray",title="Classes")
dev.off()





############Adding the cluster result to the data frame
downtown=cbind(downtown,clara.result$clustering)

albedo=readTIFF(paste(path,"albedo_1km_urbanonly.tif",sep=""))
#####Data from the GeoTIFF file
lon=seq(-119.30906568500001,length.out=dim(albedo)[2],by=0.0083333300000)
lat=seq(34.86046837500000,length.out=dim(albedo)[1],by=-0.0083333300000)


######Create a raster with the Urban Classes
require(raster)
urban_classes=rep(0,dim(albedo)[1]*dim(albedo)[2])
urban_classes[as.numeric(rownames(downtown))]=downtown[,4]+(20-13)
urban_classes_raster=raster(matrix(urban_classes,nrow=dim(albedo)[1],ncol=dim(albedo)[2]),xmn=lon[1],xmx=lon[dim(albedo)[2]],ymn=lat[dim(albedo)[1]],ymx=lat[1])
writeRaster(urban_classes_raster,filename= "/Users/pedrosequera/Desktop/Pe/CCNY-CUNY/CUERG/Thesis/LCLU_Thesis/Urban_Subsets/Urban Only/urban_classes_raster.tif",format="GTiff")




######Create a raster with the Modis IGBP Land Class file in WRF
path="/Users/pedrosequera/Desktop/Pe/CCNY-CUNY/CUERG/Thesis/LCLU_Thesis/modis_landuse_30s_WRF/"

longitude=seq(-179.99583,by=0.0083333,length.out=43200)

latitude=seq(-89.99583,by=0.0083333,length.out=21600)
finfo=file.info(paste(path,"07201-08400.14401-15600",sep=""))
finfo$size

to.read=file(paste(path,"07201-08400.14401-15600",sep=""),"rb")

l=readBin(to.read,integer(),size=1,n=1200*1200,endian="big") #1 byte for each cell



#Create a grid and plot


lat=seq(latitude[14401],latitude[15600],length.out=1200)
lon=seq(longitude[7201],longitude[8400],length.out=1200)

modis=matrix(l,1200,1200,byrow=T)

readbinary.matrix=function(matrix){
	transformed.matrix=matrix(NA,dim(matrix)[1],dim(matrix)[2])
	for(i in 1:dim(matrix)[1]){
			transformed.matrix[i,]=matrix[dim(matrix)[1]-(i-1),]
	}
	return(transformed.matrix)
}



modis_raster=raster(readbinary.matrix(modis),xmn=lon[1],xmx=lon[1200],ymn=lat[1],ymx=lat[1200])
writeRaster(modis_raster,filename= "/Users/pedrosequera/Desktop/Pe/CCNY-CUNY/CUERG/Thesis/LCLU_Thesis/Urban_Subsets/Urban Only/modis_raster.tif",format="GTiff",)

########In ENVI. Sum urban_classes_raster and the modis layer in all_with_MODIS_1km.. Name ###the file newLCLU_subset 
######Stack inclusively modis_raster.tif and newLCLU_subset... Name it stack
#######Create mask from subset layer in stack... Name it mask_subset
######Apply mask_subset in modis_raster layer in stack... Name it masked_subset
####Add masked_subset and subset layer in stack... Name modis_modified







#####Reading ENVI
path="/Users/pedrosequera/Desktop/Pe/CCNY-CUNY/CUERG/Thesis/LCLU_Thesis/Urban_Subsets/Urban Only/"

library(caTools)

modis_modified = read.ENVI(paste(path,"modis_modified",sep=""),headerfile=paste(path,"modis_modified.hdr",sep=""))

vector = c(modis_modified)

vector=as.integer(vector)
storage.mode(vector) <- "integer"
typeof(vector)
#########Re-read the vector in the WRF Format
modis_modified_matrix = matrix(vector,1200,1200)

modis_modified_vector = modis_modified_matrix[1200,]
for(i in seq(1199,1,-1)){
	modis_modified_vector = append(modis_modified_vector,modis_modified_matrix[i,])
}

#####Writing file in binary form
to.write = file(paste(path,"07201-08400_14401-15600",sep=""),"wb")

m=writeBin(modis_modified_vector,to.write,size=1,endian="big",useBytes=F)

close(to.write)

