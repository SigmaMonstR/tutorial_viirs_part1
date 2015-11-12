##Process VIIRS Monthly GeoTIFFS 
##into summary files clipped by county-month units

library(raster)
library(sp)
library(rgdal)

start <- proc.time()[3]

periods <- c("2014-9","2014-11")

for(period in periods){
  print(paste(period,": Loading files",sep=""))
  setwd(paste("[SetDir]",period,sep=""))
  
  #Set projection
  wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  
  #Load Raster
  rast <- raster(paste("VIIRS-",period,".tif",sep=""))
  projection(rast) <- CRS(wgs84)
  
  #Load Shapefile
  us = shapefile("/Users/jeffreychen/Google Drive/DOC/034-VIIRS/cb_2014_us_county_20m/cb_2014_us_county_20m.shp")
  projection(us) <- CRS(wgs84)
  
  #Extract files
  print(paste(period,": Extracting raster statistics by county",sep=""))
  
  print(paste(period,": Extracting sum",sep=""))
  x <- extract(rast, us, fun = sum)
  
  print(paste(period,": Extracting raster mean",sep=""))
  y <- extract(rast, us, fun = mean)
  
  print(paste(period,": Extracting raster SD",sep=""))
  z <- extract(rast, us, fun = sd)
  
  print(paste(period,": Extracting raster quantiles",sep=""))
  d <- extract(rast, us, fun = function(x,...){quantile(x, probs = c(0,0.05,.1, .2, 0.3,0.4,.5, 0.6,0.7, 0.8, 0.9,0.95, 1),na.rm=TRUE)}) 
  
  result = cbind(us@data$GEOID,data.frame(x, y,z, d))
  colnames(result)<-c("GEOID","sum","mean","sd", "p0","p5","p10","p20","p30","p40","p50","p60","p70","p80","p90","p95","p100")
  result$GEOID <- as.character(result$GEOID)
  result$period <- period
  
  #Write
  print(paste(period,": Saving processed file",sep=""))
  write.csv(result,paste("TNL_",period,".csv",sep=""),row.names=F)

}

end <- proc.time()[3]
print(paste("Duration: ",(end-start)/60),sep="")
