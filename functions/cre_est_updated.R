#Function to summarize population-based sampling data for population estimation using pedigree reconstruction. 

#The inputs include a dataframe Nhat.calc where total counts of sampled and inferred individuals are stored to estimate population size, a dataframe Nhat.calc.year to record annual counts sampled and inferred individuals for a particular year, an integer z to track the year of the population projection, the number of simulations (sims), a list sim.samp.list to track sampled individuals in each simulation/year/sampling intensity combination, and an age minimum to include in a population estimate.

#The output is Nhat.calc, a partially filled out dataframe that will be used to calculate joint probabilities and abundance measures.

cre.est.updated <- function(Nhat.calc, Nhat.calc.year, z, sims, sim.samp.list, study.area.est, age.min = 1, ...) {
  # First, we need to limit our summary to individuals inside of the study area
  pop.mcp <- mcp(points.t[,1], unin = "m", unout = "km2") #Calculate MCPs for all individuals
  names(pop.mcp$id) <- pop.mcp@data$id #assign individual ids to mcps
  overlap.SA <- gIntersection(spgeom1 = pop.mcp, spgeom2 = study_area_utm, byid = TRUE,id = pop.mcp@data$id)#Identify animals with MCP that overlap the study area
  pop.SA <- pop.t[which(pop.t$Individual %in% names(overlap.SA)),] #Limit inferences within study area
  
  for(s in 1:sims){ #Start loop through sampling iterations
    for(w in 1:length(sample)){ #Start nested loop for sampling intensities
      Nhat.calc.year$Year<-(z-samp.start+1) #Record year of sampling (Year 21 is year 1 od samplinfg, after 20 year burn-in)
      Nhat.calc.year$Simulation <- s #Record simulation
      Nhat.calc.year$Samp.Int <-sample[w] #Record target sampling intensity
      Nhat.calc.year$Real.Samp.Int <- length(unique(sim.samp.list[[s]][[w]]$Individual[which(sim.samp.list[[s]][[w]]$Sample == (z-samp.start+1))]))/sum(!pop.t$Individual %in% sim.samp.list[[s]][[w]]$Individual[which(sim.samp.list[[s]][[w]]$Sample <= (z-samp.start))]) #all animals sampled in a particular year and iteration (all alive)/unsampled individuals (excluding the current year of sampling). This will reveal the true sampling intensity.
      if(study.area.est == TRUE) {Nhat.calc.year$StudyAreaEst <-  TRUE} #Designate if estimate is restricted to study area
      if(study.area.est == FALSE) {Nhat.calc.year$StudyAreaEst <-  FALSE} #Designate if estimate is not restricted to study area
      Nhat.calc.year$N.sampled.alive <-sum(unique(sim.samp.list[[s]][[w]]$Individual)%in%unique(pop.t$Individual)) #Calculate the number of individuals that are sampled and currently alive
      Nhat.calc.year$N.Sampled.alive.ad <-sum(unique(pop.t[which(pop.t$Age.Class>=age.min),'Individual'])%in%(sim.samp.list[[s]][[w]]$Individual)) #Calculate the number of adults that are sampled and currently alive
      
      Nhat.calc.year$N.match <- (sum(unique(sim.samp.list[[s]][[w]]$Individual)%in%unique(sim.samp.list[[s]][[w]]$Mother)) + sum(unique(sim.samp.list[[s]][[w]]$Individual)%in%unique(sim.samp.list[[s]][[w]]$Father))) #Calculate the number of sampled individuals who are matched as parents
      Nhat.calc.year$N.match.alive <- (sum(unique(sim.samp.list[[s]][[w]]$Individual) %in% unique(sim.samp.list[[s]][[w]]$Mother) & unique(sim.samp.list[[s]][[w]]$Individual) %in% unique(pop.t$Individual[which(pop.t$Age.Class>=age.min)])) + sum(unique(sim.samp.list[[s]][[w]]$Individual)%in%unique(sim.samp.list[[s]][[w]]$Father) & unique(sim.samp.list[[s]][[w]]$Individual)%in%unique(pop.t$Individual[which(pop.t$Age.Class>=age.min)]))) #Calculate the number of sampled individuals who are matched as parents AND are currently alive
      Nhat.calc.year$N.unmatch.alive <- Nhat.calc.year$N.Sampled.alive.ad - Nhat.calc.year$N.match.alive #Calculate the number of unmatched adults as the number of sampled adults minus the number of living matched parents
      
      #Inferred individuals
      sim.samp.list[[s]][[w]]$Inf.Mother.Alive<-rep(0,nrow(sim.samp.list[[s]][[w]])) #Set placeholder for inferred mothers
      sim.samp.list[[s]][[w]]$Inf.Father.Alive<-rep(0,nrow(sim.samp.list[[s]][[w]])) #Set placeholder for inferred fathers
      
      sim.samp.list[[s]][[w]][which(sim.samp.list[[s]][[w]]$Mother %in% pop.t$Individual[which(pop.t$Age.Class>=age.min)] & !sim.samp.list[[s]][[w]]$Mother %in% sim.samp.list[[s]][[w]]$Individual & sim.samp.list[[s]][[w]]$Father %in% sim.samp.list[[s]][[w]]$Individual),'Inf.Mother.Alive']<-1 #If mom is alive but isn't in the sample set and dad was sampled, we can infer mother
    
      sim.samp.list[[s]][[w]][which(sim.samp.list[[s]][[w]]$Father %in% pop.t$Individual[which(pop.t$Age.Class>=age.min)] & sim.samp.list[[s]][[w]]$Mother %in% sim.samp.list[[s]][[w]]$Individual & !sim.samp.list[[s]][[w]]$Father %in% sim.samp.list[[s]][[w]]$Individual),'Inf.Father.Alive']<-1 #If dad is alive but isn't in the sample set and mom was sampled, we can infer father

      inf.mother.alive<-subset(sim.samp.list[[s]][[w]],sim.samp.list[[s]][[w]]$Inf.Mother.Alive==1) #Isolate all living, inferred mothers
      inf.mother.alive<-unique(inf.mother.alive$Mother) #Remove duplicate values (multiple offspring of a mother)
      inf.father.alive<-subset(sim.samp.list[[s]][[w]],sim.samp.list[[s]][[w]]$Inf.Father.Alive==1)#Isolate all living, inferred fathers
      inf.father.alive<-unique(inf.father.alive$Father) #Remove duplicate values (multiple offspring of a father)
      Nhat.calc.year$N.inferred.alive <- sum(length(inf.mother.alive)+length(inf.father.alive)) #Sum all inferred parents
      
      Nhat.calc.year$N <- nrow(pop.t) #Record total population size for year t
      Nhat.calc.year$N.ad <- length(unique(pop.t[which(pop.t$Age.Class>=age.min),1])) #record adult population size for year t
      Nhat.calc.year$N.studyarea <- nrow(pop.SA) #Record population size in study area for year t
      Nhat.calc.year$N.studyarea.ad <- length(unique(pop.SA[which(pop.SA$Age.Class>=age.min),1])) #Record adult population size in study area for year t.
  
      Nhat.calc<-rbind(Nhat.calc,Nhat.calc.year)#Add output to the running dataframe containing results from all simulations
      } #End sampling intensity
    } #End simulation iteration
  if(z==samp.start) {Nhat.calc <- Nhat.calc[-1,]} #Remove blank row as a result of coding approach
  return(Nhat.calc) #Return dataframe lisitng all simulation/sampling intensity/year combinations
  } #End function
