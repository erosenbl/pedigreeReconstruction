#Function to sample proportion of study area (spatially-based sampling approach), with sampling restricted to the study area. As individuals are sampled, they receive a sample value, either indicating that they were never sampled (0), or the sampling year during which they were sampled. 

#The inputs include an indication if this the first year pixels are being surveyed (begin.sample), a SpatialPointsDataframe of animal locations (points), a SpatialPixelsDataFrame from which pixels are sampled, and logical arguments for figures to be returned (return.plot) and whether sampling pixels should be weighted by relative abundance (stratified).

#The output is a list of objects by simulation number/sampling intensity, listing individuals sampled, their coveriates, and their known parents for abundance estimation using "cre_est_updated_spatial.R"

spatial.samp.faster <- function(begin.sample, points = points.t, sample.grid = sampleGrid.df, return.plot = FALSE, stratified = FALSE){
  
  #Step 1 - calculate relative abundance for grid cells
  #1.a. Extract number of home ranges in each grid cell
  mcp.t <- mcp(points.t[which(points.t$Winter==TRUE),1], unin = "m", unout = "km2", ) #
  abund <- over(sampleGrid.df, mcp.t, returnList = TRUE) #Calculate animals in each grid cell
  abund.count <- sapply(abund, nrow) #Saves # of animals in each grid cell
  
  #1.b. Calculate relative abundance for each grid cell (0-1, scaled on maximum abundance)
  abund.df <- data.frame(Pixel=names(abund.count), Moose = abund.count, Rel.Abund = abund.count/max(abund.count), Sampled = FALSE) #Creates dataframe of relative abunce values for each pixel
  rel.abund.grid <- merge(sampleGrid.df, abund.df, by = "Pixel",all.x = TRUE) #Combine into spatial polygon dataframe
  rel.abund.grid.SA <- rel.abund.grid[which(rel.abund.grid$In.SA==1),] #clips to study area
  raster.grid <- raster(rel.abund.grid.SA, layer = 4, values = TRUE)  #rasterises grid cells, with values for relative abundance
  poly.grid <- rasterToPolygons(raster.grid)  #convert abundance raster to polygon
  poly.grid$Pixel <- rel.abund.grid.SA$Pixel  #Label each pixel  
  row.names(poly.grid) <- as.character(poly.grid$Pixel)
  
  #Step 2 - Identify which moose exist in each pixel
  #Create list of moose that exist in each in each pixel
  df <- suppressMessages(melt(abund, level = 1))[,c("id","L1")] 
  names(df) <- c("MooseID","Pixel")
  df$Pixel <- as.integer(df$Pixel)
  
  #Step 3 - Create home ranges for all moose at this time step
  raster.list <- list(NA)
  
  for(i in 1:nrow(pop.t)){
    bivn.kde <- kde2d(sigma.list.t[[which(names(sigma.list.t)==pop.t$Individual[i])]]$Easting, sigma.list.t[[which(names(sigma.list.t)==pop.t$Individual[i])]]$Northing, n = 100)
    raster.list[[i]] <- raster(bivn.kde)
    raster.list[[i]]@data@values <- (raster.list[[i]]@data@values/max(raster.list[[i]]@data@values))
    crs(raster.list[[i]]) <- "+proj=utm +zone=19 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    names(raster.list)[i] <- as.character(pop.t$Individual[i])}
  # raster.list[[i]]@data@names <- as.character(pop.t$Individual[i])}

  #Step 4 - select grid cells for sampling
  if(begin.sample==TRUE) #First year of sampling?
  {for(s in 1:sims){ #Loop through iterations of simulation
    for(w in 1:length(spatial.sample))  {   #Loop through sampling intensity
      #Survey cells with probabilities based on weighted by abundance
      if(stratified == TRUE) #If sampling is weighted by relative abundance...
      {prob.sample <- sample(polygons(obj = poly.grid), size = length(rel.abund.grid.SA)*spatial.sample[w], replace = FALSE, prob = rel.abund.grid.SA$Rel.Abund)} #Identify grid cells sampled at intensity w
      if(stratified == FALSE) #If sampling is completely random...
      {prob.sample <- polygons(obj = poly.grid)[names(polygons(obj = poly.grid)) %in% arrange(rel.abund.grid.SA@data, desc(Moose))[1:round(length(rel.abund.grid.SA)*spatial.sample[w]),"Pixel"]]}
      sample.df <- data.frame(abund = poly.grid[which(row.names(poly.grid) %in% names(prob.sample)),"Rel.Abund"])
      prob.sample <- SpatialPolygonsDataFrame(prob.sample, sample.df) #Identify grid cells sampled at intensity w
      
  #Step 5 - Identify moose sampled, using the UD as proxy for detection
      available.moose <- as.numeric(unique(df[which(df$Pixel %in% row.names(prob.sample)),"MooseID"])) #Develop list of individuals that overlap sampled pixels
      running.df <- NA #Create storage to accumulate individuals and their KUD value for each pixel in their home range
      
      for(i in 1:length(available.moose)){ #For each available moose...
        temp.poly <- rel.abund.grid.SA[which(rel.abund.grid.SA$Pixel %in% df$Pixel[which(df$MooseID==available.moose[i])]),] #filter pixels in the study area where this particular moose resides
        temp.hr.list <- raster::extract(x = raster.list[[which(names(raster.list)==available.moose[i])]], y = temp.poly) #Extract UD values
        temp.df <- data.frame(MooseID = available.moose[i],Pixel = temp.poly$Pixel, UD = temp.hr.list) #Record moose, pixels, and respective utilization values (UD)
        running.df <- rbind(running.df,temp.df) #Combine with running dataframe
      } #End moose loop
      
      ready.to.sample.df <- merge(x = df, y = running.df, by = c("MooseID", "Pixel"), all.x=TRUE) #Combine list of moose and their pixels (df) with moose UD values for each pixel
      ready.to.sample.df$UD[is.na(ready.to.sample.df$UD)] <- 0 #Set combinations with NA to 0 (not in use)
      
      for(i in 1:nrow(ready.to.sample.df)) { #Start loop to run through whether each moose in each sampled pixel, was sampled based on their UD probability
        ready.to.sample.df$Sample[i] <- rbinom(1,1,prob = ready.to.sample.df$UD[i]) 
      } #End moose sampling loop in sampled pixel loop
      Sampled.DF <- subset(ready.to.sample.df, Sample ==1) #Filter to only moose that were sampled in previous loop
      
      for (i in 1:nrow(pop.t)){    #For every animal in the population...
        if(pop.t$Individual[i] %in% Sampled.DF$MooseID)  #If the animal was sampled
        {pop.t$Spatial.Sample[i] <- 1} #Mark moose in pop.t that were sampled in spatially-based sampling scheme.
        else  #If the animal was not sampled
        {pop.t$Spatial.Sample[i] <- 0}
      }
        
        pop.t$Spatial.Sample[which(pop.t$Age.Class==0 & pop.t$Spatial.Sample==0 & pop.t$Mother %in% pop.t$Individual[which(pop.t$Spatial.Sample==1)])] <- 1 #If sampled individual is a mother of a dependent calf, sample calf as well.
        pop.t$Spatial.Sample[which(pop.t$Male==0 & pop.t$Spatial.Sample==0 & pop.t$Individual %in% pop.t$Mother[which(pop.t$Spatial.Sample==1 & pop.t$Age.Class==0)])] <- 1 #If sampled individual is a calf dependent on its mother, sample mother as well.
        
      pop.spatial.samp <- subset(pop.t, pop.t$Spatial.Sample==1) #Save the animals that were sampled
      pop.spatial.samp$Inf.Mother<-rep(0,nrow(pop.spatial.samp)) #Fill in inferred individuals with placeholders
      pop.spatial.samp$Inf.Father<-rep(0,nrow(pop.spatial.samp)) #Fill in inferred individuals with placeholders
      pop.spatial.samp$Inf.Mother.Alive<-rep(0,nrow(pop.spatial.samp)) #Fill in inferred individuals with placeholders
      pop.spatial.samp$Inf.Father.Alive<-rep(0,nrow(pop.spatial.samp)) #Fill in inferred individuals with placeholders
      
      pop.spatial.samp.list[[w]] <- pop.spatial.samp #Store sampled individuals in a list, to add to in additional years
      names(pop.spatial.samp.list)[[w]] <- as.character(paste0("Sample ",round(spatial.sample[w]*1650), " km")) #Give names for list storage indicating the area surveyed given a particular spatial sampling intensity
    } #End sampling loop
    sim.spatial.samp.list[[s]] <- pop.spatial.samp.list #add sampling results to running list
    names(sim.spatial.samp.list)[[s]] <- s #Name list item with simulation number
  } #End simulation iteration
  } #End begin.sample if statement
  
  else #If this is not the first year of sampling (cont. step)
  {for(s in 1:sims){ #Loop through iterations of simulation
    for(w in 1:length(spatial.sample)) {    #Loop through sampling intensity
      #Survey cells with probabilities based on weighted by abundance
      if(stratified == TRUE) #If sampling is weighted by relative abundance...
      {prob.sample <- sample(polygons(obj = poly.grid), size = length(rel.abund.grid.SA)*spatial.sample[w], replace = FALSE, prob = rel.abund.grid.SA$Rel.Abund)}
      if(stratified == FALSE) #If sampling is completely random...
      {prob.sample <- polygons(obj = poly.grid)[names(polygons(obj = poly.grid)) %in% arrange(rel.abund.grid.SA@data, desc(Moose))[1:round(length(rel.abund.grid.SA)*spatial.sample[w]),"Pixel"]]}
      
      sample.df <- data.frame(abund = poly.grid[which(row.names(poly.grid) %in% names(prob.sample)),"Rel.Abund"])
      prob.sample <- SpatialPolygonsDataFrame(prob.sample, sample.df) #Identify grid cells sampled at intensity w
      
      #Step 5 - Identify moose sampled, using the UD as proxy for detection
      available.moose <- as.numeric(unique(df[which(df$Pixel %in% row.names(prob.sample)),"MooseID"]))#Develop list of individuals that overlap sampled pixels 
      running.df <- NA #Create storage to accumulate individuals and their KUD value for each pixel in their home range
      
      for(i in 1:length(available.moose)){ #For each available moose...
        temp.poly <- rel.abund.grid.SA[which(rel.abund.grid.SA$Pixel %in% df$Pixel[which(df$MooseID==available.moose[i])]),] #filter pixels in the study area where this particular moose resides
        temp.hr.list <- raster::extract(x = raster.list[[which(names(raster.list)==available.moose[i])]], y = temp.poly)#Extract UD values
        temp.df <- data.frame(MooseID = available.moose[i],Pixel = temp.poly$Pixel, UD = temp.hr.list)#Record moose, pixels, and respective utilization values (UD)
        running.df <- rbind(running.df,temp.df)#Combine with running dataframe
      }
      
      ready.to.sample.df <- merge(x = df, y = running.df, by = c("MooseID", "Pixel"), all.x=TRUE)#Combine list of moose and their pixels (df) with moose UD values for each pixel
      ready.to.sample.df$UD[is.na(ready.to.sample.df$UD)] <- 0#Set combinations with NA to 0 (not in use)
      
      for(i in 1:nrow(ready.to.sample.df)) {#Start loop to run through whether each moose in each sampled pixel, was sampled based on their UD probability
        ready.to.sample.df$Sample[i] <- rbinom(1,1,prob = ready.to.sample.df$UD[i])
      }#End moose sampling loop in sampled pixel loop
      Sampled.DF <- subset(ready.to.sample.df, Sample == 1)#Filter to only moose that were sampled in loop
      
      for (i in 1:nrow(pop.t)){    #For every animal in the population...
        if (pop.t$Individual[i] %in% sim.spatial.samp.list[[s]][[w]]$Individual) {pop.t$Spatial.Sample[i] <- sim.spatial.samp.list[[s]][[w]]$Spatial.Sample[which(sim.spatial.samp.list[[s]][[w]]$Individual==pop.t$Individual[i])]} #Ensures that sampled animals remain recorded as sampled in the first year it was sampled
        if(pop.t$Individual[i] %in% Sampled.DF$MooseID & !pop.t$Individual[i] %in% sim.spatial.samp.list[[s]][[w]]$Individual)  #If the animal was sampled for the first time
        {pop.t$Spatial.Sample[i] <- z-samp.start+1} #Mark moose in pop.t that were sampled in spatially-based sampling scheme in year z-samp.start+1 (equals 2-5)
        if (!pop.t$Individual[i] %in% Sampled.DF$MooseID & !pop.t$Individual[i] %in% sim.spatial.samp.list[[s]][[w]]$Individual) #If an individual wasn't sampled this time around or previously, mark sample state 0
        {pop.t$Spatial.Sample[i] <- 0}
      } #End individual loop
      
      pop.t$Spatial.Sample[which(pop.t$Age.Class==0 & pop.t$Spatial.Sample==0 & pop.t$Mother %in% pop.t$Individual[which(pop.t$Spatial.Sample==z-samp.start+1)])] <- z-samp.start+1 #If sampled individual is a mother of a dependent calf, sample calf as well.
      pop.t$Spatial.Sample[which(pop.t$Male==0 & pop.t$Spatial.Sample==0 & pop.t$Individual %in% pop.t$Mother[which(pop.t$Spatial.Sample==z-samp.start+1 & pop.t$Age.Class==0)])] <- z-samp.start+1 #If sampled individual is a calf dependent on its mother, sample mother as well.
      # 
    pop.spatial.samp <- subset(pop.t, pop.t$Spatial.Sample==z-samp.start+1) #Save the animals that were sampled
    pop.spatial.samp$Inf.Mother<-rep(0,nrow(pop.spatial.samp)) #Fill in inferred individuals with placeholders
    pop.spatial.samp$Inf.Father<-rep(0,nrow(pop.spatial.samp))#Fill in inferred individuals with placeholders
    pop.spatial.samp$Inf.Mother.Alive<-rep(0,nrow(pop.spatial.samp))#Fill in inferred individuals with placeholders
    pop.spatial.samp$Inf.Father.Alive<-rep(0,nrow(pop.spatial.samp))#Fill in inferred individuals with placeholders
    
    sim.spatial.samp.list[[s]][[w]]$Age.Class <- sim.spatial.samp.list[[s]][[w]]$Age.Class + 1 #Increase the age of everyone who was sampled in previous years (this wasn't done in main code when everyone had their birthday)
    sim.spatial.samp.list[[s]][[w]]<-rbind(sim.spatial.samp.list[[s]][[w]],pop.spatial.samp) #Add results to running list
    } #End sampling loop
  } #End simulation iteration
} #End begin.sample if statement

  if(return.plot==TRUE) #Plots if requested
  {plot <- spplot(rel.abund.grid.SA, "Moose", col.regions=bpy.colors(20)) +
    latticeExtra::layer(sp.polygons(prob.sample, lwd = 4))
  return(plot)}

  return(sim.spatial.samp.list) #returns a list of objects by simulation number/sampling intensity, listing individuals sampled, their coveriates, and their known parents for abundance estimation using "cre_est_updated_spatial.R"
    
}