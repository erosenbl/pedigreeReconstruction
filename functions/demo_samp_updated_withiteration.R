#Function to sample proportion of population, with sampling restricted to the study area. As sampled individuals are sampled, they receive a sample value, either indicating that they were never sampled (0), or the sampling year during which they were sampled. 

#The inputs include an indication if this the first year animals are being surveyed (begin.sample), a running list of sampled individuals for each simulation:sample intensity combination (sim.samp.list), an integer z to track the year of the population projection, the number of simulations (sims), and whether to exclude animals with home ranges completely outside of the delineated study area (exclude.outliers).

#The output is  a list object tracking all living individuals in a population for a given iteration and sampling intensity, listing their identities, their covariates, when they are sampled, who their parents are, and whether we could infer their parents, living or dead, for abundance estimation using "cre_est_updated.R"

demo.samp.updated.withiteration <- function(begin.sample, sim.samp.list, z, sims, exclude.outliers,...) {
  if(exclude.outliers==TRUE) #If we want to exclude animals outside of our designated study area from sampling (note: animal maximum convex polygon (mcp) has to be completely outside of study area)
    {pop.mcp <- mcp(points.t[,1], unin = "m", unout = "km2") #Calculate MCPs for all individuals
    names(pop.mcp$id) <- pop.mcp@data$id #assign individual ids to mcps
    overlap.SA <- gIntersection(spgeom1 = pop.mcp, spgeom2 = study_area_utm, byid = TRUE,id = pop.mcp@data$id) #Identify animals with MCP that overlap the study area
  if(begin.sample==TRUE) #If this is the first year of sampling...
  {    
    pop.samp.list <- list() #Creates storage
    for(s in 1:sims){ #Start loop through iterations of simulation
    for(w in 1:length(sample))  {   #start nested loop through sampling intensity
    pop.t$Sample.P <- rep(0,nrow(pop.t))
    {for (i in 1:nrow(pop.t)){    #For every animal in the population...
      if(!pop.t$Individual[i] %in% names(overlap.SA))  #If the animal is not in the study area...
      {pop.t$Sample.P[i] <- NA #the animal has no sample probability
      pop.t$Sample[i] <- 0} #and the animal is not sampled
      else  #Otherwise, the animal is in the study area, and...
      {pop.t$Sample.P[i]<-runif(1,min=0,max=1) #we generate a random sample probability for the animal...
      pop.t$Sample[i]<-ifelse(pop.t$Sample.P[i]<=sample[w],1,0) #and the animal will be sampled if Sample.P is <= sampling intensity
      }
    } #End loop for every animal in the population
    pop.t$Sample[which(pop.t$Age.Class==0 & pop.t$Sample==0 & pop.t$Mother %in% pop.t$Individual[which(pop.t$Sample==1)])] <- 1 #If sampled individual is a mother of a dependent calf, sample calf as well.
    pop.t$Sample[which(pop.t$Male==0 & pop.t$Sample==0 & pop.t$Individual %in% pop.t$Mother[which(pop.t$Sample==1 & pop.t$Age.Class==0)])] <- 1 #If sampled individual is a calf dependent on its mother, sample mother as well.
    
    pop.samp <- subset(pop.t, pop.t$Sample==1)} #Save only animals that were sampled
    pop.samp$Inf.Mother<-rep(0,nrow(pop.samp)) #Fill in inferred individuals with placeholder
    pop.samp$Inf.Father<-rep(0,nrow(pop.samp)) #Fill in inferred individuals with placeholder
    pop.samp$Inf.Mother.Alive<-rep(0,nrow(pop.samp)) #Fill in inferred individuals with placeholder
    pop.samp$Inf.Father.Alive<-rep(0,nrow(pop.samp)) #Fill in inferred individuals with placeholder
    
    pop.samp.list[[w]] <- pop.samp #Store sampled individuals in a list, to add to in additional years of sampling
    names(pop.samp.list)[[w]] <- as.character(paste0("Sample ",w, "0%")) #Give names for list storage
  } #End sample intensity loop
    sim.samp.list[[s]] <- pop.samp.list #Add animals sampled in all sampling intensities to grand list
    names(sim.samp.list)[[s]] <- s #Add iteration as label
  } #End simulation iteration loop
  } #End begin.sample if statement loop
  
  else #If this is after the first year of sampling
  {for(s in 1:sims){ #Start loop through iterations of simulation
    for(w in 1:length(sample)) {    #start nested loop through sampling intensity
    for (i in 1:nrow(pop.t)){ #For every animal in the population...
      if (pop.t$Individual[i] %in% sim.samp.list[[s]][[w]]$Individual) {pop.t$Sample[i] <- sim.samp.list[[s]][[w]]$Sample[which(sim.samp.list[[s]][[w]]$Individual==pop.t$Individual[i])]} #Ensures that sampled animals remain sampled
      if (!pop.t$Individual[i] %in% sim.samp.list[[s]][[w]]$Individual & pop.t$Individual[i] %in% names(overlap.SA)) #If an individual has not been sampled AND overlaps the study area  
      { pop.t$Sample[i]==0 #Set their sample status to 0...
        pop.t$Sample.P[i]<-runif(1,min=0,max=1) #... we generate a random sample probability for the animal...
        pop.t$Sample[i]<-ifelse(pop.t$Sample.P[i]<=sample[w],(z-samp.start+1),0)} #and the animal will be sampled if Sample.P is <= sampling intensity, with a sample status matching the year of sampling.
    } #End individual loop
      pop.t$Sample[which(pop.t$Age.Class==0 & pop.t$Sample==0 & pop.t$Mother %in% pop.t$Individual[which(pop.t$Sample==z-samp.start+1)])] <- z-samp.start+1 #If sampled individual is a mother of a dependent calf, sample calf as well.
      pop.t$Sample[which(pop.t$Male==0 & pop.t$Sample==0 & pop.t$Individual %in% pop.t$Mother[which(pop.t$Sample==z-samp.start+1 & pop.t$Age.Class==0)])] <- z-samp.start+1 #If sampled individual is a calf dependent on its mother, sample mother as well.
    pop.sampled<-subset(pop.t, pop.t$Sample==(z-samp.start+1)) #Save only animals that were sampled
    pop.sampled$Inf.Mother<-rep(0,nrow(pop.sampled))#Fill in inferred individuals with placeholder
    pop.sampled$Inf.Father<-rep(0,nrow(pop.sampled))#Fill in inferred individuals with placeholder
    pop.sampled$Inf.Mother.Alive<-rep(0,nrow(pop.sampled))#Fill in inferred individuals with placeholder
    pop.sampled$Inf.Father.Alive<-rep(0,nrow(pop.sampled))#Fill in inferred individuals with placeholder
    sim.samp.list[[s]][[w]]$Age.Class <- sim.samp.list[[s]][[w]]$Age.Class + 1 #Increase the age of everyone who was sampled in previous years (this wasn't done in main code when everyone had their birthday)
    sim.samp.list[[s]][[w]]<-rbind(sim.samp.list[[s]][[w]],pop.sampled)} #End sample intensity loop
  } #End simulation iteration loop
  } #end else statement
}#End begin.sample else statement loop

  if(exclude.outliers==FALSE) #IF we don't want to exclude animals outside of the study area
  { if(begin.sample==TRUE) #If this is the first year of sampling...
  {for(s in 1:sims){ #Start loop through iterations of simulation
    for(w in 1:length(sample))  {   #start nested loop through sampling intensity
    pop.t$Sample.P <- rep(0,nrow(pop.t))
    for (i in 1:nrow(pop.t)){    #For every animal in the population...
      pop.t$Sample.P[i]<-runif(1,min=0,max=1) #generate random number for each animal...
      pop.t$Sample[i]<-ifelse(pop.t$Sample.P[i]<=sample[w],1,0) #and the animal will be sampled if Sample.P is <= sampling intensity
    }
    pop.t$Sample[which(pop.t$Age.Class==0 & pop.t$Sample==0 & pop.t$Mother %in% pop.t$Individual[which(pop.t$Sample==1)])] <- 1 #If sampled individual is a mother of a dependent calf, sample calf as well.
    pop.t$Sample[which(pop.t$Male==0 & pop.t$Sample==0 & pop.t$Individual %in% pop.t$Mother[which(pop.t$Sample==1 & pop.t$Age.Class==0)])] <- 1 #If sampled individual is a calf dependent on its mother, sample mother as well.
    
    pop.samp <- subset(pop.t, pop.t$Sample==1) #Save the animals that were sampled
    pop.samp$Inf.Mother<-rep(0,nrow(pop.samp)) #Fill in inferred individuals with placeholder
    pop.samp$Inf.Father<-rep(0,nrow(pop.samp)) #Fill in inferred individuals with placeholder
    pop.samp$Inf.Mother.Alive<-rep(0,nrow(pop.samp))#Fill in inferred individuals with placeholder
    pop.samp$Inf.Father.Alive<-rep(0,nrow(pop.samp))#Fill in inferred individuals with placeholder
    
    pop.samp.list[[w]] <- pop.samp #Store sampled individuals in a list, to add to in additional years
    names(pop.samp.list)[w] <- as.character(paste0("Sample ",w, "0%")) #Give names for list storage
  } #End sampling loop
    sim.samp.list[[s]] <- pop.samp.list
    names(sim.samp.list)[[s]] <- s
  } #End simulation iteration
  } #End begin.sample if statement
    
    else #If this is after the first year of sampling
    {for(s in 1:sims) {#Start loop through iterations of simulation
      for(w in 1:length(sample))  {   #start nested loop through sampling intensity
      for (i in 1:nrow(pop.t)){ #For every animal in the population...
        if (pop.t$Individual[i] %in% sim.samp.list[[s]][[w]]$Individual) {pop.t$Sample[i] <- sim.samp.list[[s]][[w]]$Sample[which(sim.samp.list[[s]][[w]]$Individual==pop.t$Individual[i])]} #Ensures that sampled animals remain sampled
        if (!pop.t$Individual[i] %in% sim.samp.list[[s]][[w]]$Individual)#If an individual has not been sampled AND overlaps the study area
        { pop.t$Sample[i]==0 #Set their sample status to 0...
          pop.t$Sample.P[i]<-runif(1,min=0,max=1)#... we generate a random sample probability for the animal...
          pop.t$Sample[i]<-ifelse(pop.t$Sample.P[i]<=sample[w],(z-samp.start+1),0)}#and the animal will be sampled if Sample.P is <= sampling intensity, with a sample status matching the year of sampling.
      } #End individual loop
        pop.t$Sample[which(pop.t$Age.Class==0 & pop.t$Sample==0 & pop.t$Mother %in% pop.t$Individual[which(pop.t$Sample==z-samp.start+1)])] <- z-samp.start+1 #If sampled individual is a mother of a dependent calf, sample calf as well.
        pop.t$Sample[which(pop.t$Male==0 & pop.t$Sample==0 & pop.t$Individual %in% pop.t$Mother[which(pop.t$Sample==z-samp.start+1 & pop.t$Age.Class==0)])] <- z-samp.start+1 #If sampled individual is a calf dependent on its mother, sample mother as well.
      pop.sampled<-subset(pop.t, pop.t$Sample==(z-samp.start+1))
      pop.sampled$Inf.Mother<-rep(0,nrow(pop.sampled))#Fill in inferred individuals with placeholder
      pop.sampled$Inf.Father<-rep(0,nrow(pop.sampled))#Fill in inferred individuals with placeholder
      pop.sampled$Inf.Mother.Alive<-rep(0,nrow(pop.sampled))#Fill in inferred individuals with placeholder
      pop.sampled$Inf.Father.Alive<-rep(0,nrow(pop.sampled))#Fill in inferred individuals with placeholder
      sim.samp.list[[s]][[w]]$Age.Class <- sim.samp.list[[s]][[w]]$Age.Class + 1 #Increase the age of everyone who was sampled in previous years (this wasn't done in main code when everyone had their birthday)
      sim.samp.list[[s]][[w]]<-rbind(sim.samp.list[[s]][[w]],pop.sampled)} #End sampling loop
    } # end simulation iteration
    } #end else statement
  } #end sampling excludes outlying animals
  
  return(sim.samp.list) #Returns a list object tracking all living individuals in a population for a given iteration and sampling intensity, listing their identities, their covariates, when they are sampled, who their parents are, and whether we could infer their parents, living or dead.
}
