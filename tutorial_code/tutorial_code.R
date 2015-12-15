#################################
##VIIRS Day/Night Bands Data  ##
#################################

##############
##OVERVIEW ##
#############

# #The processing steps are as follows:
# (1) Load a GeoTIFF (geo-referenced raster file), and assign a commonly used spatial projection that dictates the shape of a map canvas;
# (2) Visually inspect rasters by color coding raster values on a map using k-means clustering;
# (3) Crop raster data by geographic boundaries and extract associated radiance summary statistics;
# (4) Combine raster summary data with demographic data to uncover underlying patterns


##############
##Libraries ##
#############
#Start off by specifying the working directory as well as calling 7 libraries:
# - **doParallel**: Allows for parallel processing using multiple cores
# - **foreach**: Enables for-each statements to be used for parallel processing
# - **raster**: interface to raster (e.g. TIFF) data
# - **sp**: interface to spatial data processing
# - **rgdal**: interface for manipulating spatial data
# - **ggmap**: contains a Google geocoding wrapper 
# - **plotly**: API interface to plotly's graphing libraries providing interactive web visualizations

library(doParallel)
library(foreach)
library(raster)
library(sp)
library(rgdal)
library(ggmap)
library(plotly)

##############################
##Section 1: Obtain and Load##
##############################

###Obtain and Load VIIRS Day/Night Band: Satellite Data
# Satellite data is complex and processed in a number of ways to correct for a manifold of environmental conditions like stray light. 
# In the case of analyzing population demographics, we use VIIRS DNB monthly composites that omit records that have been effected by stray light. To obtain the data:
# 
# - Raster files can be obtained [here](http://ngdc.noaa.gov/eog/viirs/download_monthly.html);
# - Under each month's folder, select **Tile1_75N180W**, which contains data for North America;
# - Within this folder, download the file labeled **VCMCFG** containing data that excludes any stray light. 
# - There are two files. The file ending in **avg_rade9** contains average radiances; this should be set aside for use. The other file ending in "" contains the number of cloud free pixels included in each pixel's calculation.
# 
# Create a folder labeled "imagery" and save the **avg_rade9** .tif file into that directory. Below, we then specify the path for future use, open the first and only raster in the list of .tif files, then assign WGS84 as the spatial projection.

##Set directory path
  imagery = "your-path-here/imagery"

##Obtain a list of TIF files, load in the first file in list
  tifs = list.files(imagery,pattern = "\\.tif")
  rast <- raster(paste(imagery,"/",tifs[1],sep=""))

##Specify WGS84 as the projection of the raster file
  wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  projection(rast) <- CRS(wgs84)

###Metropolitan Statistical Areas (MSAs) data
  In addition, we use  MSA shapefiles and population data from the US Census Bureau. To draw down the shapefile, we write a function to download the shapefiles from the US Census Bureau website, then assign the same WGS84 projection to the shapefile. By assigning the same spatial projection, we can work across shapefiles and raster files.

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

  ##Download Population Estimates from Census
    msa_pop <- read.csv("http://www.census.gov/popest/data/metro/totals/2014/files/CBSA-EST2014-alldata.csv")
    msa_pop <- msa_pop[msa_pop$LSAD=="Metropolitan Statistical Area",]
    msa_pop <- msa_pop[order(msa_pop$POPESTIMATE2014),]
    msa_pop$NAME <- as.character(msa_pop$NAME) 


#######################
##Section 2: Mapping ##
#######################
  
##35 largest cities
cities <- c("New York, NY", "Los Angeles, CA","Chicago, IL", "Houston, TX",
            "Philadelphia, PA", "Phoenix, AZ", "San Antonio, TX", "San Diego, CA",     
            "Dallas, TX", "San Jose, CA", "Austin, TX", "Jacksonville, FL",
            "San Francisco, CA", "Indianapolis, IN", "Columbus, OH", "Fort Worth, TX",
            "Charlotte, NC", "Detroit, MI", "El Paso, TX", "Seattle, WA",
            "Denver, CO","Washington, DC", "Memphis, TN", "Boston, MA",
            "Nashville, TN", "Baltimore, MD", "Oklahoma City, OK", "Portland, OR",
            "Las Vegas, NV", "Louisville, KY","Milwaukee, WI","Albuquerque, NM",
            "Tucson, AZ","Fresno, CA","Sacramento, CA")

##Set graph layout:layout with no margins (mai), with 7 rows and 5 columns (mfrow), with a navy blue background (bg). 
  par(mai=c(0,0,0,0),mfrow = c(7,5),bg='#001a4d', bty='n')

#The process maps as follows:
# - First set a placeholder dataframe for coordinates
# - For each city, geocode the city name for the centroid using the geocode() function from ggmaps, append the coordinates to the placeholder. This coordinates file will be used again later.
# - Use the extent function to specify the spatial bounding box (the frame around the city) as +/-1 degree longitude and +/-0.25 degree latitude.
# - Crop the raster file (rast) by the bounding box.
# - Convert the cropped raster into a vector in order to get the radiances, then run a k-means clustering algorithm to find natural intervals within the radiance distribution. For each cluster, extract the maximum radiance.
# - Map the city with a navy to yellow color palette with intervals from the k-means clustering.

##Loop through data to create map with intervals calibrated to each city
    coords <- data.frame() ##placeholder

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

##Create map with light intervals based on national sample
  #Set layout
  par(mai=c(0,0,0,0),mfrow = c(7,5),bg='#001a4d', bty='n')
  
  #Run clustering
  set.seed(123) #set seed for reproducibility
  sampled <- sample(rast, 20000) #sample 20,000 pixels
  clusters <- 15 ##15 clusters
  clust <- kmeans(sampled,clusters)$cluster
  combined <- as.data.frame(cbind(sampled,clust))
  brk <- sort(aggregate(combined[,1], list(combined[,2]), max)[,2])
  
  ##Loop through each city
  for(i in 1:length(cities)){
    
    temp_coord <- coords[i,] ##re-use the coordinates 
    e <- extent(temp_coord$lon - 1, temp_coord$lon + 1,
                temp_coord$lat - 0.25, temp_coord$lat + 0.25)
    rc <- crop(rast, e)    
    
    #Plots
    plot(rc, breaks=brk, col=colorRampPalette(c("#001a4d","#0066FF", "yellow"))(clusters), 
         legend=F,yaxt='n',xaxt='n',frame = F, asp=1.5)
    text(temp_coord$lon ,temp_coord$lat + 0.15,
         substr(cities[i],1,regexpr(",",cities[i])-1), 
         col="white", cex=1.25)
    
    rm(combined)
  }

  
#######################
##Section 3: Mapping ##
#######################
  
# The function accepts three parameters:
# 
# - **shp**. A shapefile from the US Census Bureau
# - **rast**. A geotiff raster file
# - **i** is an index value that corresponds to an index position within a list of polygons in the shapefile.

  masq <- function(shp,rast,i){
  
    #Extract one polygon based on index value i
    polygon <- shp[i,] #extract one polygon
    extent <- extent(polygon) #extract the polygon extent 
    
    #Raster extract
    outer <- crop(rast, extent) #extract raster by polygon extent
    inner <- mask(outer,polygon) #keeps values from raster extract that are within polygon
    
    #Convert cropped raster into a vector
    #Specify coordinates
    coords <- expand.grid(seq(extent@xmin,extent@xmax,(extent@xmax-extent@xmin)/(ncol(inner)-1)),
    seq(extent@ymin,extent@ymax,(extent@ymax-extent@ymin)/(nrow(inner)-1)))
    #Convert raster into vector
    data <- as.vector(inner)
    
    #package data in neat dataframe
    data <- cbind(as.character(shp@data$CBSAFP[i]),coords, data) 
    colnames(data)<-c("GEOID","lon","lat","avg_rad") #note that 
    data <- data[!is.na(data$avg_rad),] #keep non-NA values only
    
    return(data)
  }

#####################################  
###Section 4: Maps Visualizations####
#####################################

  ##MSAs by GEOID
    msa_list <- c(16180,19140,45820,42540,35620)
  
  ##Placeholder
    radiances <- data.frame() 
  
  ##Loop MSA file
  for(i in msa_list){
    
    print(i)
    
    #Extract MSA i polygon
    shp_temp <- msa[msa@data$GEOID==i,]
    
    #Extract MSA abbreviated name
    if(regexpr("-",as.character(shp_temp@data$NAME)[1])[1]==-1){
      loc = as.character(substr(as.character(shp_temp@data$NAME)[1],1,regexpr(",",as.character(shp_temp@data$NAME)[1])-1))
    } else{
      loc = as.character(substr(as.character(shp_temp@data$NAME)[1],1,regexpr("-",as.character(shp_temp@data$NAME)[1])-1))
    }
    
    #Extract the radiances, append to radiances placeholder
    rad <- masq(shp_temp,rast,1)$avg_rad 
    temp <- data.frame(loc = as.character(paste(loc,"(TNL = ",round(sum(rad),0),")",sep="")), avg_rad = rad) 
    radiances <- rbind(radiances,temp)
  }

#Use ggplot to create histograms by MSA group. Preload.
  ggplot(radiances, aes(x=log(avg_rad))) +
    geom_histogram(position="identity", alpha=0.4) +
    facet_grid(. ~ loc)

#Remove all axes labels for style
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

#Initiate a plotly graph without axes
  ggplotly()  %>% layout(xaxis=x, yaxis=y)

####Scatterplot of TNL vs. Population####
  
##Extract and calculate TNL for each US MSA, parallelized
  registerDoParallel(cores=2)
  extract <- foreach(i=1:nrow(msa@data), .combine=rbind) %dopar% {
    data <- masq(msa,rast,i)
    data.frame(GEOID = data$GEOID[1],sum = sum(data$avg_rad))
  }
  extract$GEOID <- as.numeric(as.character(extract$GEOID))

##Join in data
  joined<-merge(extract, msa_pop[,c("CBSA","NAME","POPESTIMATE2014")],by.x="GEOID",by.y="CBSA")
  
  colnames(joined) <- c("GEOID","TNL","MSA","Population")

##Plot scatterplot
  plot_ly(joined, 
          x = log(TNL), 
          y = log(Population), 
          text = paste("MSA: ", MSA),
          mode = "markers", 
          color = TNL,colors="PuOr")  %>% 
    layout(title="Total Nighttime Light vs. Population", showlegend = F)
