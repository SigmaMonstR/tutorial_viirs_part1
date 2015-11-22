##Functions for manipulating VIIRS GeoTIFFs (ngdc.noaa.gov/eog/viirs/download_monthly.html)
##and Census TIGER county files (http://www2.census.gov/geo/tiger/GENZ2014/shp/cb_2014_us_county_20m.zip)

### MASQ: Function for clipping GeoTIFFs by shapefile boundaries
### masq(shp, rast, i)
# shp = polygon shapefile name
# rast = raster file name
# i = index of number of polygons in shapefile
### Note that function was written to handle a Census county shapefile

masq <- function(shp,rast,i){
  
  #Setup 
    polygon <- shp[i,] #extract one polygon
    extent <- extent(polygon) #single polygon extent
  
  #Raster extract
    outer <- crop(rast, extent) #extract raster by polygon extent
    inner <- mask(outer,polygon) #keeps values from raster extract that are within polygon
    
  #Convert raster to vector of radiances (measures in pixels)
    coords <- expand.grid(seq(extent@xmin,extent@xmax,(extent@xmax-extent@xmin)/(ncol(inner)-1)),
                          seq(extent@ymin,extent@ymax,(extent@ymax-extent@ymin)/(nrow(inner)-1)))
    data <- as.vector(inner)
    
  #Pull together County GEOID from shapefile, lat/lon, and radiance measure
    data <- cbind(paste(shp@data$STATEFP[i],shp@data$COUNTYFP[i],sep=""),coords, data) 
    colnames(data)<-c("GEOID","lon","lat","avg_rad") #note that 
    data <- data[!is.na(data$avg_rad),] #keep non-NA values only
  
  return(data)
}

### MASQ_QUANTILE: Function that extracts zonal statistics from rasters
# df = dataframe name
# id = field name of row id ("GEOID")
# col = field name for main variable ("avg_rad")
# int = interval of quantile calculation

masq_quantile <- function(df, id, col, int){
  
  #Calculate quantiles
    values <- as.data.frame(quantile(df[[col]], probs = seq(0,1,int),na.rm=TRUE))
    
  #Populate row identifiers (percentage, GEOID)
    values$pct <- row.names(values)
    colnames(values) <- c("avg_rad","pct")
    values$GEOID <- df[[id]][1]
  
  #Calculate Total Nighttime Light
    values$sum <- sum(df[[col]])

  return(values)
  
}
