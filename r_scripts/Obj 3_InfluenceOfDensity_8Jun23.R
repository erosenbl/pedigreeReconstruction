###Objective 3: Performance of estimator across varying local densities
##Goal: Delineate study area into 9 equal "study areas", with differing densities in each one to see how the pedigree estimator works at estimating densities across a range of possible values

#load required packages
library(adehabitatHR)
library(sf)

#Create storage for estimators
Nhat.calc.density <- data.frame(Simulation=NA,Area.sampled=NA, Inf.Grid.Cell=NA, Real.Samp.Int = NA,N.sampled.alive=NA,N.Sampled.alive.ad=NA, N.match = NA, N.match.alive = NA, N.unmatch.alive = NA,N.inferred.alive=NA, LCL=NA, Nhat=NA, UCL=NA, N=NA, N.ad=NA) 

###Step 1: Prepare spatial data and delineate 9 equal study areas####
#Bring in study area shape file and prep into proper projection
moose.study.area <- shapefile(x = "data/Study_Area/StudyDissolve.shp")  #Load study area
proj4string <- "+proj=utm +zone=19 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"  #Set projection of study area
moose.study.area <- spTransform(moose.study.area, proj4string, inverse=TRUE) #Transform shapefile projection

study_area_utm <- CRS("+init=epsg:32619") %>% 
  spTransform(moose.study.area, .) #ensure CRS is correct for study area

###Develop square sampling grid with 1 sq km grid cells, that can be used for spatial sampling later on in this simulation
grid <- raster(extent(study_area_utm), resolution = c(1000,1000), crs = proj4string(study_area_utm))#Create raster over study area
grid <- raster::extend(grid, c(1,1)) #Expand raster grid by 1 cell to fully encompass study area
gridPolygon <- rasterToPolygons(grid) #Convert raster cells to grid cells
gridPolygon$layer <- c(1:2925) #Give numeric valeus to each grid cell 
gridPoints <- rasterToPoints(grid, spatial=TRUE)
# sample_grid <- SpatialPixels(gridPoints, proj4string = proj4string(study_area_utm))

#Develop Inference Grid (9 cells of known density to sample from and estimate with pedigree reconstruction)
big_grid_dim_x <- (extent(study_area_utm)@xmax - extent(study_area_utm)@xmin)/3 #Determine width of inference grid cell (x dimension)
big_grid_dim_y <- (extent(study_area_utm)@ymax - extent(study_area_utm)@ymin)/3 #Determine height of inference grid cell (y dimension)

grid.big <- raster(extent(study_area_utm), resolution = c(big_grid_dim_x,big_grid_dim_y), crs = proj4string(study_area_utm)) #Create large grid cells
gridPolygon.big <- rasterToPolygons(grid.big)
gridPolygon.big$layer <- c(1,2,3,4,5,6,7,8,9)

#Assign sample grid cells to inference grid cells
pixel_in_cell <- rbind(over(gridPolygon, gridPolygon.big))
sample_grid_df <- SpatialPixelsDataFrame(points = gridPoints, data = data.frame(Pixel = gridPolygon$layer,  inf_cell = pixel_in_cell), proj4string = proj4string(study_area_utm))
names(sample_grid_df) <- c("Pixel","inf_cell")

#Read in population and point data from year 5 of objective 1 and 2
pop.obj3 <- readRDS("..data/popend.16Nov21.RDS")
points.obj3 <- readRDS("..data/pointsend.16Nov21.RDS")

#Filter moose that do not have summer points
points.obj3 <- points.obj3[which(points.obj3$Winter==TRUE),]
pop.obj3 <- pop.obj3[which(pop.obj3$Individual %in% points.obj3$id),]

#Calculate density of moose for each of 9 inference cells
mcp.inf <- mcp(points.obj3[,1], unin = "m", unout = "km2")
abund.inf <- over(gridPolygon.big, mcp.inf, returnList = TRUE) #Calculate animals in each grid cell

abund.count.inf <- sapply(abund.inf, nrow) #Saves # of animals in each grid cell
density.inf <- abund.count.inf/(raster::area(gridPolygon.big)/1000000)

###Step 2: Select grid cells to sample, and sample####
#Step 2.1 - calculate relative abundance for grid cells
#2.1.a. Extract number of home ranges in each grid cell
mcp.t <- mcp(points.obj3[which(points.obj3$Winter==TRUE),1], unin = "m", unout = "km2", )
abund <- over(sample_grid_df, mcp.t, returnList = TRUE) #Calculate animals in each grid cell
abund.count <- sapply(abund, nrow) #Saves # of animals in each grid cell
abund.count <- data.frame(Pixel = 1:length(abund.count), Abund = abund.count)

#2.1.b. Calculate relative abundance for each grid cell (0-1, scaled on maximum abundance)
abund.df <- data.frame(Pixel=(abund.count$Pixel), Moose = abund.count$Abund, Rel.Abund = abund.count$Abund/max(abund.count$Abund), inf_cell = sample_grid_df$inf_cell, Sampled = FALSE)
rel.abund.grid <- merge(sample_grid_df, abund.df, by = "Pixel",all.x = TRUE) #Combine into spatial polygon dataframe
raster.grid <- raster(rel.abund.grid, layer = 4, values = TRUE)  #rasterises grid cells, with values for relative abundance
poly.grid <- rasterToPolygons(raster.grid)  #convert abundance raster to polygon
poly.grid$Pixel <- rel.abund.grid$Pixel  #Label each pixel  
poly.grid$inf_cell <- sample_grid_df$inf_cell
row.names(poly.grid) <- as.character(poly.grid$Pixel)
#

#Step 2.2 - Identify which moose exist in each pixel
#Create list of moose that exist in each in each pixel
df <- suppressMessages(reshape2::melt(abund, level = 1))[,c("id","L1")] 
names(df) <- c("MooseID","Pixel")
df$Pixel <- as.integer(df$Pixel)

#Step 2.3 - Create home ranges for all moose at this time step
raster.list <- list(NA)

for(i in 1:nrow(pop.obj3)){
  bivn.kde <- kde2d(points.obj3@data[which(points.obj3@data$id==pop.obj3$Individual[i]),]$Easting, points.obj3@data[which(points.obj3@data$id==pop.obj3$Individual[i]),]$Northing, n = c(45,65), lims = c(262589.5,306589.5,4924607.3,4988607.3))
  raster.list[[i]] <- raster(bivn.kde)
  raster.list[[i]]@data@values <- (raster.list[[i]]@data@values/max(raster.list[[i]]@data@values))
  crs(raster.list[[i]]) <- "+proj=utm +zone=19 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  names(raster.list)[i] <- as.character(pop.obj3$Individual[i])}

#Step 2.4 - select grid cells for sampling
#Set simulation length and spatial sample
sims <- 100 #100 iterations
age.min <- 1 #minimum age for including in population count
spatial.sample <- round(seq(0.01,0.2,length.out=10)*298,0) #vector of spatial sampling effort
pop.temp <- pop.obj3 #assign population

Nhat.calc.density.complete <- NA #create storage for abundance estimates
running.df <- NA #create storage for moose UD values

#Extract all home ranges, across all pixels
for(i in 1:length(raster.list)){ #For each available moose
  temp.hr.list <- raster::extract(x = raster.list[[i]], y = sample_grid_df) #...calculate UD for each pixel
  temp.df <- data.frame(MooseID = names(raster.list)[i],Pixel = sample_grid_df$Pixel, UD = temp.hr.list) #... and save pixel x moose dataset
  running.df <- rbind(running.df,temp.df)
}

#Sample (This will take a while)
for(r in 1:9){ #Start inference cell loop (9 cells total)
for(s in 1:sims){ #Loop through iterations of simulation
  for(w in 1:length(spatial.sample))  {   #Loop through sampling intensity
    prob.sample <- sample(polygons(obj = poly.grid[which(poly.grid$inf_cell==r),]), size = spatial.sample[w], replace = FALSE, prob = poly.grid[which(poly.grid$inf_cell==r),]$Rel.Abund)#Survey cells with probabilities based on weighted by abundance
    
    sample.df <- data.frame(Pixel = poly.grid[which(poly.grid$Pixel %in% names(prob.sample)),"Pixel"],abund = poly.grid[which(poly.grid$Pixel %in% names(prob.sample)),"Rel.Abund"], poly.grid[which(poly.grid$Pixel %in% names(prob.sample)),"inf_cell"]) #extract cells based on those sampled in prob.sample.
    
    prob.sample <- SpatialPolygonsDataFrame(prob.sample, sample.df) #Make spatial
    
    ready.to.sample.df <- merge(x = running.df, y = prob.sample, by = "Pixel")
    ready.to.sample.df<- ready.to.sample.df[!is.na(ready.to.sample.df$UD),]
    
    for(i in 1:nrow(ready.to.sample.df)) {
      ready.to.sample.df$Sample[i] <- rbinom(1,1,prob = ready.to.sample.df$UD[i]) #Runs through moose in the sampled pixels, and sample them based on their UD value in a given pixel.
    }
    Sampled.DF <- subset(ready.to.sample.df, Sample ==1)
    
    for (i in 1:nrow(pop.temp)){    #For every animal in the population...
      if(pop.temp$Individual[i] %in% Sampled.DF$MooseID)  #If the animal was sampled
      {pop.temp$Spatial.Sample[i] <- 1}
      else  #If the animal was not sampled
      {pop.temp$Spatial.Sample[i] <- 0}
    }
    
    pop.temp$Spatial.Sample[which(pop.temp$Age.Class==0 & pop.temp$Spatial.Sample==0 & pop.temp$Mother %in% pop.temp$Individual[which(pop.temp$Spatial.Sample==1)])] <- 1 #If sampled individual is a mother of a dependent calf, sample calf as well.
    pop.temp$Spatial.Sample[which(pop.temp$Male==0 & pop.temp$Spatial.Sample==0 & pop.temp$Individual %in% pop.temp$Mother[which(pop.temp$Spatial.Sample==1 & pop.temp$Age.Class==0)])] <- 1 #If sampled individual is a calf dependent on its mother, sample mother as well.
    
    pop.spatial.samp <- subset(pop.temp, pop.temp$Spatial.Sample==1) #Save the animals that were sampled
    pop.spatial.samp$Inf.Mother<-rep(0,nrow(pop.spatial.samp)) #Fill in inferred individuals with placeholders
    pop.spatial.samp$Inf.Father<-rep(0,nrow(pop.spatial.samp)) #Fill in inferred individuals with placeholders
    pop.spatial.samp$Inf.Mother.Alive<-rep(0,nrow(pop.spatial.samp)) #Fill in inferred individuals with placeholders
    pop.spatial.samp$Inf.Father.Alive<-rep(0,nrow(pop.spatial.samp)) #Fill in inferred individuals with placeholders
    
    #Start to fill out out abundance estimate for a particular iteration
    Nhat.calc.density$Simulation <- s #iteration
    Nhat.calc.density$Area.sampled <-  spatial.sample[w] #Sample effort
    Nhat.calc.density$Inf.Grid.Cell <- r #which of the 9 areas are we estimating abundance for
    Nhat.calc.density$Real.Samp.Int <- nrow(pop.spatial.samp)/abund.count.inf[r] #what was the realized population sample intensity
    Nhat.calc.density$N.sampled.alive <- nrow(pop.spatial.samp) #Total number of living individuals sampled
    
    Nhat.calc.density$N.Sampled.alive.ad <- nrow(pop.spatial.samp[which(pop.spatial.samp$Age.Class>=1),]) #Total number of living adults sampled
   
     Nhat.calc.density$N.match <- (sum(unique(pop.spatial.samp$Individual)%in%unique(pop.spatial.samp$Mother)) + sum(unique(pop.spatial.samp$Individual)%in%unique(pop.spatial.samp$Father))) #How many parent-offspring pairs have been identified?
    
     Nhat.calc.density$N.match.alive <- (sum(unique(pop.spatial.samp$Individual) %in% unique(pop.spatial.samp$Mother) & unique(pop.spatial.samp$Individual) %in% unique(pop.temp$Individual[which(pop.temp$Age.Class>=age.min)])) + sum(unique(pop.spatial.samp$Individual)%in%unique(pop.spatial.samp$Father) & unique(pop.spatial.samp$Individual)%in%unique(pop.temp$Individual[which(pop.temp$Age.Class>=age.min)]))) #How many parent-offspring pairs have been identified with both individuals still alive?
    
     Nhat.calc.density$N.unmatch.alive <- Nhat.calc.density$N.Sampled.alive.ad - Nhat.calc.density$N.match.alive #How many sampled, living individuals were not matched with another sampled individual?
     
     #Inferring only live individuals
     Nhat.calc.density$Inf.Mother.Alive<-0 #placeholder
     Nhat.calc.density$Inf.Father.Alive<-0 #placeholder
     
     Nhat.calc.density$Inf.Mother.Alive <-  nrow(pop.temp[which(pop.temp$Spatial.Sample==1 & pop.temp$Mother %in% pop.temp$Individual[which(pop.temp$Age.Class>=1)] & !pop.temp$Mother %in% pop.temp$Individual[which(pop.temp$Spatial.Sample==1)] & pop.temp$Father %in% pop.temp$Individual[which(pop.temp$Spatial.Sample==1)]),]) #If mom is alive but isn't in the sample set and dad was sampled, we can infer mother
     
     Nhat.calc.density$Inf.Father.Alive <-  nrow(pop.temp[which(pop.temp$Spatial.Sample==1 & pop.temp$Father %in% pop.temp$Individual[which(pop.temp$Age.Class>=1)] & pop.temp$Mother %in% pop.temp$Individual[which(pop.temp$Spatial.Sample==1)] & !pop.temp$Father %in% pop.temp$Individual[which(pop.temp$Spatial.Sample==1)]),]) #If dad is alive but isn't in the sample set and mom was sampled, we can infer father
    
     Nhat.calc.density$N <- nrow(pop.temp[which(pop.temp$Individual %in% abund.inf[[r]]$id),]) #Calculate true abundance
     Nhat.calc.density$N.ad <- nrow(pop.temp[which(pop.temp$Individual %in% abund.inf[[r]]$id & pop.temp$Age.Class>=age.min),]) #Calculate true adult abundance
     
    Nhat.calc.density$TrueDensity <- density.inf[r] #record true density for 1 of the 9 inferred densities
     
     Nhat.calc.density.complete <- rbind(Nhat.calc.density.complete, Nhat.calc.density) #add calculations to storage object
     
  } #End sampling loop
  print(paste0("Simulation ", s, ", Iteration ", w))
} #End simulation loop
  print(r)
} #End inference cell loop (r)

#Some clean up needs to happen now...
Nhat.calc.density.complete <- Nhat.calc.density.complete[-1,] #drops first blank row

Nhat.calc.density.complete$N.inferred.alive <- Nhat.calc.density.complete$Inf.Mother.Alive + Nhat.calc.density.complete$Inf.Father.Alive #Calculate sum of living, inferred mothers and fathers.

Nhat.calc.density.complete$Inf.Grid.Cell <- as.factor(Nhat.calc.density.complete$Inf.Grid.Cell) #Turn inference grid cell (1-9) into a factor.

uninformative.prior <- 1 #Set uninformative prior to calculate population size.

### Step 3: Calculate population size####
joint.prob <- data.frame(Area.sampled = NA, Density = NA, det = NA, match = NA) #Create storage for storing joint probabilities to estimate abundance

for(i in 1:nrow(Nhat.calc.density.complete)){ 
  #Estimate shape parameters for probability of detection
  pdet.alpha <- uninformative.prior + Nhat.calc.density.complete$N.match.alive[i] 
  pdet.beta <- uninformative.prior + Nhat.calc.density.complete$N.inferred.alive[i]
  #Estimate shape parameters for probability of matching
  pmatch.alpha <- uninformative.prior + Nhat.calc.density.complete$N.match.alive[i]
  pmatch.beta <- uninformative.prior + (Nhat.calc.density.complete$N.Sampled.alive.ad[i] - Nhat.calc.density.complete$N.match.alive[i])
  
  det <- rbeta(n = 100,shape1=pdet.alpha, shape2 = pdet.beta)
  match  <-rbeta(n = 100,shape1=pmatch.alpha, shape2  = pmatch.beta)
  
  temp.joint.prob <- data.frame(Area.sampled = rep(Nhat.calc.density.complete$Area.sampled[i], length(det)), Density = rep(Nhat.calc.density.complete$TrueDensity[i], length(det)), det = det, match = match) #Temporary storage for known density and joint probabilities
  
  joint.prob <- rbind(joint.prob,temp.joint.prob)#Add temporary storage to final storage for all iterations in all grid cells.
  
  nhat.sa <- Nhat.calc.density.complete$N.match.alive[i]/(det*match) #Estimate abundance for a particular iteration
  
  Nhat.calc.density.complete$Nhat[i] <- (mean(nhat.sa) - Nhat.calc.density.complete$N.ad[i])/Nhat.calc.density.complete$N.ad[i] #Record average adult abundance estimate relative to known adult abundance N.ad
  Nhat.calc.density.complete$LCL[i] <- (quantile(x = nhat.sa, probs = c(0.025)) - Nhat.calc.density.complete$N.ad[i])/Nhat.calc.density.complete$N.ad[i] #Record 95% LCL
  Nhat.calc.density.complete$UCL[i] <-  (quantile(x = nhat.sa, probs = (.975)) - Nhat.calc.density.complete$N.ad[i])/Nhat.calc.density.complete$N.ad[i] #Record 95%UCL
  print(i)
}

joint.prob <- joint.prob[-1,]
joint.prob$Density <- round(joint.prob$Density,3)

Nhat.calc.density.complete #Return final results object
