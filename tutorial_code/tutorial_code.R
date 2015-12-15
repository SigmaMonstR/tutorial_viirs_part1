#################################
##VIIRS Day/Night Bands Data  ##
#################################

############
##OVERVIEW##
############

library(doParallel)
library(foreach)
library(raster)
library(sp)
library(rgdal)
library(ggmap)
library(plotly)

imagery = "/Users/jcc6/Google Drive/DOC/034-VIIRS/VIIRS_composites/imagery"
output = "/Users/jcc6/Google Drive/DOC/034-VIIRS/VIIRS_composites/output"
shp = "/Users/jcc6/Google Drive/DOC/034-VIIRS/VIIRS_composites/shp"

setwd(output)

##Grab list of TIF files
tifs = list.files(imagery,pattern = "\\.tif")

##Load in shapefile
wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

#Load in single test raster
rast <- raster(paste(imagery,"/",tifs[1],sep=""))
projection(rast) <- CRS(wgs84)

#Cities
cities <- c("New York, NY", "Los Angeles, CA","Chicago, IL", "Houston, TX",
            "Philadelphia, PA", "Phoenix, AZ", "San Antonio, TX", "San Diego, CA",     
            "Dallas, TX", "San Jose, CA", "Austin, TX", "Jacksonville, FL",
            "San Francisco, CA", "Indianapolis, IN", "Columbus, OH", "Fort Worth, TX",
            "Charlotte, NC", "Detroit, MI", "El Paso, TX", "Seattle, WA",
            "Denver, CO","Washington, DC", "Memphis, TN", "Boston, MA",
            "Nashville, TN", "Baltimore, MD", "Oklahoma City, OK", "Portland, OR",
            "Las Vegas, NV", "Louisville, KY","Milwaukee, WI","Albuquerque, NM",
            "Tucson, AZ","Fresno, CA","Sacramento, CA")

#MAP MONTAGE OF CITIES: Color intervals scaled to each respective city
par(mai=c(0,0,0,0),mfrow = c(7,5),bg='#001a4d', bty='n')
coords <- data.frame()
for(i in 1:length(cities)){
  
  ##Coords
  temp_coord <- geocode(cities[i], source = "google")
  coords <- rbind(coords,temp_coord)
  
  e <- extent(temp_coord$lon - 1, temp_coord$lon + 1,
              temp_coord$lat - 0.25, temp_coord$lat + 0.25)
  rc <- crop(rast, e)    
  
  ##Rescale brackets
  sampled <- as.vector(rc)
  clusters <- 15
  clust <- kmeans(sampled,clusters)$cluster
  combined <- as.data.frame(cbind(sampled,clust))
  brk <- sort(aggregate(combined[,1], list(combined[,2]), max)[,2])
  
  #Plots
  plot(rc, breaks=brk, col=colorRampPalette(c("#001a4d","#0066FF", "yellow"))(clusters), 
       legend=F,yaxt='n',xaxt='n',frame = F, asp=1.5)
  text(temp_coord$lon ,temp_coord$lat + 0.15,
       substr(cities[i],1,regexpr(",",cities[i])-1), 
       col="white", cex=1.25)
  
  rm(combined)
}

#MAP MONTAGE OF CITIES: Color intervals scaled to national sample
par(mai=c(0,0,0,0),mfrow = c(7,5),bg='#001a4d',bty='n')

set.seed(123)
sampled <- sample(rast, 20000)
clusters <- 15
clust <- kmeans(sampled,clusters)$cluster
combined <- as.data.frame(cbind(sampled,clust))
brk <- sort(aggregate(combined[,1], list(combined[,2]), max)[,2])

for(i in 1:length(cities)){
  
  ##Coords
  temp_coord <- coords[i,]
  e <- extent(temp_coord$lon - 1, temp_coord$lon + 1,
              temp_coord$lat - 0.25, temp_coord$lat + 0.25)
  rc <- crop(rast, e)    
  
  #Plots
  plot(rc, breaks=brk, col=colorRampPalette(c("#001a4d","#0066FF","yellow"))(clusters), 
       legend=F,yaxt='n',xaxt='n',frame = F, asp=1.5)
  text(temp_coord$lon ,temp_coord$lat + 0.15,
       substr(cities[i],1,regexpr(",",cities[i])-1), 
       col="white", cex=1.25)
  
  rm(combined)
}

##Radiance extraction function
#Function that extracts vector of radiances from GeoTIFF by shapefile polygon
##Written for Census bureau MSA shapefiles. Easily modified for others
masq <- function(shp,rast,i){
  
  #Setup 
  polygon <- shp[i,] #extract one polygon
  extent <- extent(polygon) #single polygon extent
  
  #Raster extract
  outer <- crop(rast, extent) #extract raster by polygon extent
  inner <- mask(outer,polygon) #keeps values from raster extract that are within polygon
  
  #Convert to matrix
  coords <- expand.grid(seq(extent@xmin,extent@xmax,(extent@xmax-extent@xmin)/(ncol(inner)-1)),
                        seq(extent@ymin,extent@ymax,(extent@ymax-extent@ymin)/(nrow(inner)-1)))
  data <- as.vector(inner)
  data <- cbind(as.character(shp@data$CBSAFP[i]),coords, data)  ##CBSAFP needs to be replaced with the Geo ID when using other shapefiles
  colnames(data)<-c("GEOID","lon","lat","avg_rad") #note that 
  data <- data[!is.na(data$avg_rad),] #keep non-NA values only
  
  return(data)
}

#MSA data from US Census Bureau
##Draw down MSA Shapefile
shape_direct <- function(url, shp) {
  library(rgdal)
  temp = tempfile()
  download.file(url, temp) ##download the URL taret to the temp file
  unzip(temp,exdir=getwd()) ##unzip that file
  return(readOGR(paste(shp,".shp",sep=""),shp))
}

msa <- shape_direct(url="http://www2.census.gov/geo/tiger/GENZ2014/shp/cb_2014_us_cbsa_20m.zip", 
                    shp= "cb_2014_us_cbsa_20m")
projection(msa) <- CRS(wgs84)

##Download and load MSA files from Census
msa_pop <- read.csv("http://www.census.gov/popest/data/metro/totals/2014/files/CBSA-EST2014-alldata.csv")
msa_pop <- msa_pop[msa_pop$LSAD=="Metropolitan Statistical Area",]
msa_pop <- msa_pop[order(msa_pop$POPESTIMATE2014),]
msa_pop$NAME <- as.character(msa_pop$NAME)

##Extract radiances by MSA
msa_list <- c(16180,19140,45820,42540,35620)

radiances <- data.frame()
for(i in msa_list){
  print(i)
  shp_temp <- msa[msa@data$GEOID==i,]
  if(regexpr("-",as.character(shp_temp@data$NAME)[1])[1]==-1){
    loc = as.character(substr(as.character(shp_temp@data$NAME)[1],1,regexpr(",",as.character(shp_temp@data$NAME)[1])-1))
  } else{
    loc = as.character(substr(as.character(shp_temp@data$NAME)[1],1,regexpr("-",as.character(shp_temp@data$NAME)[1])-1))
  }
  rad <- masq(shp_temp,rast,1)$avg_rad
  temp <- data.frame(loc = as.character(paste(loc,"(TNL = ",round(sum(rad),0),")",sep="")),avg_rad = rad)
  radiances <- rbind(radiances,temp)
}

radiances$log_avg_rad <- log(radiances$avg_rad)


##ggplot to prime the plotly function
ggplot(radiances, aes(x=log(avg_rad))) +
  geom_histogram(position="identity", alpha=0.4) +
  facet_grid(. ~ loc)


##Distribution comparison, set up styles
x <- list(
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
)
y <- list(
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
) 

#Plotly graph using ggplot functionality
ggplotly()  %>% layout(xaxis=x, yaxis=y)


#In direct comparison, TNL and population are positively correlated with a relatively high correlation coefficient of 0.78, showing promise as a proxy and can be tightened by incorporating the radiance distribution for approximating density. 


##Set up comparisons
registerDoParallel(cores=2)
extract <- foreach(i=1:nrow(msa@data), .combine=rbind) %dopar% {
  data <- masq(msa,rast,i)
  data.frame(GEOID = data$GEOID[1],sum = sum(data$avg_rad))
}
extract$GEOID <- as.numeric(as.character(extract$GEOID))

##Join in data
joined<-merge(extract, msa_pop[,c("CBSA","NAME","POPESTIMATE2014")],by.x="GEOID",by.y="CBSA")

colnames(joined) <- c("GEOID","TNL","MSA","Population")

plot_ly(joined, 
        x = log(TNL), 
        y = log(Population), 
        text = paste("MSA: ", MSA),
        mode = "markers", 
        color = TNL,colors="PuOr")  %>% 
  layout(title="Total Nighttime Light vs. Population", showlegend = F)

