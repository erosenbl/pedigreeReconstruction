#Function to summarize data for population estimation from spatial sampling

#The inputs include a dataframe Nhat.calc.spatial where total counts of sampled and inferred individuals are stored to estimate population size, a dataframe Nhat.calc.year.spatial to record annual counts sampled and inferred individuals for a particular year, an integer z to track the year of the population projection, the number of simulations (sims), a list sim.spatial.samp.list to track sampled individuals in each simulation/year/sampling intensity combination, and an age minimum to include in a population estimate.

#The output is Nhat.calc.spatial, a partially filled out dataframe that will be used to calculate joint probabilities and abundance measures.

cre.est.updated.spatial <- function(Nhat.calc.spatial, Nhat.calc.year.spatial, z, sims, sim.spatial.samp.list, study.area.est, age.min = 1, ...) {
  # First, we need to limit our summary to individuals inside of the study area
  pop.mcp <- mcp(points.t[,1], unin = "m", unout = "km2") #Calculate MCPs for all individuals
  names(pop.mcp$id) <- pop.mcp@data$id #assign individual ids to mcps
  overlap.SA <- gIntersection(spgeom1 = pop.mcp, spgeom2 = study_area_utm, byid = TRUE,id = pop.mcp@data$id)#Identify animals with MCP that overlap the study area
  pop.SA <- pop.t[which(pop.t$Individual %in% names(overlap.SA)),]#Limit inferences within study area
  
  for(s in 1:sims){ #Start loop through sampling iterations
    for(w in 1:length(spatial.sample)){ #Start nested loop for sampling intensities
      Nhat.calc.year.spatial$Year<-(z-samp.start+1)#Record year of sampling (Year 21 is year 1 od samplinfg, after 20 year burn-in)
      Nhat.calc.year.spatial$Simulation <- s #Record simulation
      Nhat.calc.year.spatial$Samp.Int <-spatial.sample[w] #Record target sampling intensity
      Nhat.calc.year.spatial$Real.Samp.Int <- length(unique(sim.spatial.samp.list[[s]][[w]]$Individual[which(sim.spatial.samp.list[[s]][[w]]$Spatial.Sample == (z-samp.start+1))]))/sum(!pop.t$Individual %in% sim.spatial.samp.list[[s]][[w]]$Individual[which(sim.spatial.samp.list[[s]][[w]]$Spatial.Sample <= (z-samp.start))]) #all animals sampled in a particular year and iteration (all alive)/unsampled individuals (excluding the current year of sampling). This will reveal the true sampling intensity.

      if(study.area.est == TRUE) {Nhat.calc.year.spatial$StudyAreaEst <-  TRUE} #Set placeholder for inferred mothers
      if(study.area.est == FALSE) {Nhat.calc.year.spatial$StudyAreaEst <-  FALSE} #Set placeholder for inferred fathers
      
      Nhat.calc.year.spatial$N.sampled.alive <-sum(unique(sim.spatial.samp.list[[s]][[w]]$Individual)%in%unique(pop.t$Individual))#Calculate the number of individuals that are sampled and currently alive
      Nhat.calc.year.spatial$N.Sampled.alive.ad <-sum(unique(pop.t[which(pop.t$Age.Class>=age.min),'Individual'])%in%(sim.spatial.samp.list[[s]][[w]]$Individual)) #Calculate the number of adults that are sampled and currently alive

      Nhat.calc.year.spatial$N.match <- (sum(unique(sim.spatial.samp.list[[s]][[w]]$Individual)%in%unique(sim.spatial.samp.list[[s]][[w]]$Mother)) + sum(unique(sim.spatial.samp.list[[s]][[w]]$Individual)%in%unique(sim.spatial.samp.list[[s]][[w]]$Father))) #Calculate the number of sampled individuals who are matched as parents
      Nhat.calc.year.spatial$N.match.alive <- (sum(unique(sim.spatial.samp.list[[s]][[w]]$Individual) %in% unique(sim.spatial.samp.list[[s]][[w]]$Mother) & unique(sim.spatial.samp.list[[s]][[w]]$Individual) %in% unique(pop.t$Individual[which(pop.t$Age.Class>=age.min)])) + sum(unique(sim.spatial.samp.list[[s]][[w]]$Individual)%in%unique(sim.spatial.samp.list[[s]][[w]]$Father) & unique(sim.spatial.samp.list[[s]][[w]]$Individual)%in%unique(pop.t$Individual[which(pop.t$Age.Class>=age.min)]))) #Calculate the number of sampled individuals who are matched as parents AND are currently alive
      Nhat.calc.year.spatial$N.unmatch.alive <- Nhat.calc.year.spatial$N.Sampled.alive.ad - Nhat.calc.year.spatial$N.match.alive #Calculate the number of unmatched adults as the number of sampled adults minus the number of living matched parents
      
      #Inferred individuals
      sim.spatial.samp.list[[s]][[w]]$Inf.Mother.Alive<-rep(0,nrow(sim.spatial.samp.list[[s]][[w]])) #Set placeholder for inferred mothers
      sim.spatial.samp.list[[s]][[w]]$Inf.Father.Alive<-rep(0,nrow(sim.spatial.samp.list[[s]][[w]])) #Set placeholder for inferred fathers
      
      sim.spatial.samp.list[[s]][[w]][which(sim.spatial.samp.list[[s]][[w]]$Mother %in% pop.t$Individual[which(pop.t$Age.Class>=age.min)] & !sim.spatial.samp.list[[s]][[w]]$Mother %in% sim.spatial.samp.list[[s]][[w]]$Individual & sim.spatial.samp.list[[s]][[w]]$Father %in% sim.spatial.samp.list[[s]][[w]]$Individual),'Inf.Mother.Alive']<-1 #If mom is alive but isn't in the sample set and dad was sampled, we can infer mother
    
      sim.spatial.samp.list[[s]][[w]][which(sim.spatial.samp.list[[s]][[w]]$Father %in% pop.t$Individual[which(pop.t$Age.Class>=age.min)] & sim.spatial.samp.list[[s]][[w]]$Mother %in% sim.spatial.samp.list[[s]][[w]]$Individual & !sim.spatial.samp.list[[s]][[w]]$Father %in% sim.spatial.samp.list[[s]][[w]]$Individual),'Inf.Father.Alive']<-1 #If dad is alive but isn't in the sample set and mom was sampled, we can infer father

      inf.mother.alive<-subset(sim.spatial.samp.list[[s]][[w]],sim.spatial.samp.list[[s]][[w]]$Inf.Mother.Alive==1) #Isolate all living, inferred mothers
      inf.mother.alive<-unique(inf.mother.alive$Mother) #Remove duplicate values (multiple offspring of a mother)
      inf.father.alive<-subset(sim.spatial.samp.list[[s]][[w]],sim.spatial.samp.list[[s]][[w]]$Inf.Father.Alive==1) #Isolate all living, inferred fathers
      inf.father.alive<-unique(inf.father.alive$Father) #Remove duplicate values (multiple offspring of a father)
      Nhat.calc.year.spatial$N.inferred.alive <- sum(length(inf.mother.alive)+length(inf.father.alive)) #Sum all inferred parents
      
      Nhat.calc.year.spatial$N <- nrow(pop.t) #Record total population size for year t
      Nhat.calc.year.spatial$N.ad <- length(unique(pop.t[which(pop.t$Age.Class>=age.min),1])) #record adult population size for year t
      Nhat.calc.year.spatial$N.studyarea <- nrow(pop.SA) #Record population size in study area for year t
      Nhat.calc.year.spatial$N.studyarea.ad <- length(unique(pop.SA[which(pop.SA$Age.Class>=age.min),1])) #Record adult population size in study area for year t.
  
      Nhat.calc.spatial<-rbind(Nhat.calc.spatial,Nhat.calc.year.spatial)#Add output to the running dataframe containing results from all simulations
      } #End sampling intensity
    } #End simulation iteration
  if(z==samp.start) {Nhat.calc.spatial <- Nhat.calc.spatial[-1,]} #Remove blank row as a result of coding approach
  return(Nhat.calc.spatial) #Return dataframe lisitng all simulation/sampling intensity/year combinations
  } #End function
