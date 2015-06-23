library(raster)
library(gstat)
library(spacetime)
library(reshape)
library(plot3D)
library(zoo)
library(rasterVis)
workDir <- "/Users/liziqi/Desktop/Nor_100X100"
setwd(workDir)
filenames <- list.files(pattern="*.tif", full.names=TRUE)

#create raster brick
r<-brick()
for(i in 1:120){
    r<-addLayer(r,raster(filenames[i]))
}

#get best summer and winter
summer = r$LST_2003_06
winter = r$LST_2008_01
worst = r$LST_2011_08
#interpolation
stfdf.winter = as.data.frame(winter,row.names=NULL)
stfdf.winter = melt(stfdf.winter)
stfdf.winter = stfdf.winter[,2]
stfdf.winter = data.frame(values = signif(stfdf.winter,6))
stfdf.winter = stfdf.winter-273.16
sp = SpatialPoints(coordinates(winter))

winter.df = SpatialPointsDataFrame(sp, stfdf.winter,proj4string=CRS(winter))
www = variogram(values~1, winter.df,width=1000)
www.fit = fit.variogram(www, vgm(40, "Exp", 50000, 1))
krig <- krige.cv(values~1, winter.df, vgm(18, "Exp", 50000, 1), nmax = 100, nfold=5)
idw <- krige.cv(values~1, winter.df, nmax = 100, nfold=5)#idw
fit = lm(krig$observed ~krig$var1.pred)
plot(x$observed~x$var1.pred)
summary(fit)

#summer
stfdf.summer = as.data.frame(summer,row.names=NULL)
stfdf.summer = melt(stfdf.summer)
stfdf.summer = stfdf.summer[,2]
stfdf.summer = data.frame(values = signif(stfdf.summer,6))
stfdf.summer = stfdf.summer-273.16
sp = SpatialPoints(coordinates(winter))
summer.df = SpatialPointsDataFrame(sp, stfdf.summer,proj4string=CRS(summer))

sss = variogram(values~1, summer.df,width=1000)
sss.fit = fit.variogram(sss, vgm(40, "Exp", 50000, 1))
x <- krige.cv(values~1, summer.df, vgm(40, "Exp", 50000, 1), nmax = 100, nfold=5)
x <- krige.cv(values~1, summer.df, nmax = 100, nfold=5)#idw
fit = lm(x$observed ~x$var1.pred)
summary(fit)

#plot
plot(idw$observed,idw$var1.pred,xlim = c(-50,-25),ylim=c(-50,-25),xlab="IDW_Predicted (°C)",ylab="Observed (°C)",main="IDW Cross Validation")
abline(a=0,b=1,col="red")
dev.copy(png,'myplot.png')
dev.off()

plot(krig$observed,krig$var1.pred,xlim = c(-50,-25),ylim=c(-50,-25),xlab="Krig_Predicted (°C)",ylab="Observed (°C)",main="Kriging Cross Validation")
abline(a=0,b=1,col="red")

plot(idw$observed,idw$residual,ylim = c(-5,5),xlim = c(-25,-50),xlab="IDW_Observed (°C)",ylab="Residuals (°C)",main="IDW Cross Validation")
plot(krig$observed,krig$residual,ylim = c(-5,5),xlim = c(-25,-50),xlab="Krig_Observed (°C)",ylab="Residuals (°C)",main="Kriging Cross Validation")
dev.copy(png,'myplot.png')
dev.off()

spplot(worst,at= seq(5,25,1),col.regions=colorRampPalette(c('blue', 'yellow','red')),ylab="Northing",xlab="Easting")

#get monthly annual average
monthly = ts(c(cellStats(r, 'mean')) ,start = c(2002,1),frequency = 12)
monthly.mean = tapply(monthly, cycle(air.ts), mean,na.rm=T)

#deseasonalized
for(i in 1:120){
    mon = i%%12
    if (mon==0) {mon = 12}
    r[[i]] = r[[i]] - monthly.mean[mon][[1]]
}


#variogram
mydata = as.data.frame(r,row.names=NULL)
mydata = melt(mydata)
mydata = mydata[,2]
mydata = data.frame(values = signif(mydata,6))
#mydata = mydata-273.16
YM <- as.yearmon(2002 + seq(0, 119)/12)
sp = SpatialPoints(coordinates(r))
stfdf = STFDF(sp, YM, mydata)
stplot(stfdf,col.regions=rev(heat.colors(60)))
ddd = variogramST(values~1,stfdf,width = 5000,cutoff=50000,tlags=0:5)

#get Map of NAs
na.count = is.na(subset(r,1))
for(i in 2:120){
    nalayer = is.na(subset(r,i))
    na.count = na.count + nalayer
}

na.count = as.data.frame(na.count,row.names=NULL)
na.count = melt(na.count)
na.count = na.count[,2]
na.count = data.frame(values = signif(na.count,6))
sp = SpatialPoints(coordinates(winter))

na.df = SpatialPointsDataFrame(sp, na.count,proj4string=CRS(winter))
na.df =as.data.frame(na.df)
scatterplot3d(na.df$x,na.df$y,na.df$values,type="h",pch="",zlab="Missing Value Counts",xlab="Easting",ylab="Northing",zlim=c(1,60),lwd = 2)

plot3D(na.count)

#
r1 = subset(r,1)
r1.df = as.data.frame(r1,row.names=NULL)
r1.df = melt(mydata)
r1.df = mydata[,2]
r1.df = data.frame(values = signif(r1.df,6))
r1.df = r1.df-273.16
YM <- as.yearmon("2002-01")
stfdf = STFDF(sp, YM, r1.df)



#2002 only
workDir <- "/Users/liziqi/Desktop/Nor_2002"
setwd(workDir)
filenames <- list.files(pattern="*.tif", full.names=TRUE)
r.2002<-brick()
for(i in 1:12){
    r.2002<-addLayer(r.2002,raster(filenames[i]))
}
for(i in 1:12){
    mon = i%%12
    if (mon==0) {mon = 12}
    r.2002[[i]] = r.2002[[i]] - monthly.mean[mon][[1]]
}

mydata = as.data.frame(r.2002,row.names=NULL)
mydata = melt(mydata)
mydata = mydata[,2]
mydata = data.frame(values = signif(mydata,6))

YM <- as.yearmon(2002 + seq(0, 11)/12)
sp = SpatialPoints(coordinates(r.2002))
stfdf = STFDF(sp, YM, mydata)
ddd = variogramST(values~1,stfdf,width = 5000,cutoff=50000,tlags=0:12)

plot(ddd,wireframe=T,col.regions=bpy.colors(1000),convertMonths=T,xlab=list("distance (km)", rot=30),plot.numbers=T,auto.key=T,all=T,scales=list(arrows=F, z = list(distance = 5)),zlim=c(0,1000),zlab="semivariance")



prodSumModel <- vgmST("productSum", space=vgm(16, "Exp", 40000), time= vgm(1000, "Per", 1),sill = 50,stAni=200000,nugget=0)



x=sample(1:100, 30, replace=F)
y=sample(1:100, 30, replace=F)
for (i in 1:30){
    cell<-c(r[x[i],y[i]])
    cell.up<-c(r[x[i],y[i]-1])
    cell.down<-c(r[x[i],y[i]+1])
    cell.left<-c(r[x[i]-1,y[i]])
    cell.right<-c(r[x[i]+1,y[i]])
    fit.up<-lm(formula = cell ~ cell.up, na.rm = TRUE)
    fit.down<-lm(formula = cell ~ cell.down, na.rm = TRUE)
    fit.left<-lm(formula = cell ~ cell.left, na.rm = TRUE)
    fit.right<-lm(formula = cell ~ cell.right, na.rm = TRUE)
    udlr<-c(summary(fit.up)$r.squared, summary(fit.down)$r.squared, summary(fit.left)$r.squared, summary(fit.right)$r.squared)
    print(udlr)
}

#Get station data
a = cbind(88.30,69.33)
xy = SpatialPoints(a,CRS("+init=epsg:4326"))
xy = spTransform(xy,crs("+init=epsg:3408"))
air = extract(r,xy)
m = colMeans(matrix(air, nrow=12))
m = m-273.16
acf(na.omit(air),120,xlab="time lag (month)",ylab="correlation",main="Temporal Autocorrelation",type= "correlation",ylim=c(-1,1))


