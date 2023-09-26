#Function to calculate detection and match probability shape parameters, and then calculate abundance using pedigree reconstruction under spatial sampling scheme. 

#The inputs include an object storing simulations from an individually-based population model with spatial sampling. The structure should follow that created by the simulation loop in SimulationCode_Rosenblatt_et_al.R (L201-403). It is assumed that this function will be applied to a 25 year simulation, with year 1-20 serving as a burn-in period, and years 21-25 serving as the years sampled

#The output is the same as the input, with additional columns with each simulations average abundance measure, and various confidence intervals and ranges reported.

nHat_est_spatial <- function(spatial.simulation.dataset) {
  
  #First, we create storage for posterior estimates of probability of detection (pdet.para.blank) and probability of matching individuals (pmatch.para.blank), using sample data from the population-based approach. These storage ojects store beta distribution shape parameters alpha and beta.
  pdet.spatial.para.blank <- list(data.frame("Samp.Int"  =  unique(spatial.simulation.dataset$Samp.Int),  "Year 1"  =  NA,  "Year 2"  =  NA,  "Year 3"  =  NA,  "Year 4"  =  NA,  "Year 5"  =  NA),  data.frame("Samp.Int"  =  unique(spatial.simulation.dataset$Samp.Int),  "Year 1"  =  NA,  "Year 2"  =  NA,  "Year 3"  =  NA,  "Year 4"  =  NA,  "Year 5"  =  NA))
  names(pdet.spatial.para.blank) <- c("alpha", "beta")
  
  pmatch.spatial.para.blank <- list(data.frame("Samp.Int"  =  unique(spatial.simulation.dataset$Samp.Int),  "Year 1"  =  NA,  "Year 2"  =  NA,  "Year 3"  =  NA,  "Year 4"  =  NA,  "Year 5"  =  NA),  data.frame("Samp.Int"  =  unique(spatial.simulation.dataset$Samp.Int),  "Year 1"  =  NA,  "Year 2"  =  NA,  "Year 3"  =  NA,  "Year 4"  =  NA,  "Year 5"  =  NA))
  names(pmatch.spatial.para.blank) <- c("alpha", "beta")
  
  #Create storage for beta distribution parameters - these lists will include copies of pdet.para and pmatch.para, each corresponding to a simulation.
  pdet.spatial.list <- list(NA)
  pmatch.spatial.list <- list(NA)
  
  uninformative.prior <- 1 #Set uninformative prior for all parameter estimation (results in a uniform beta distribution)
  
  for(i in 1:length(unique(spatial.simulation.dataset$Simulation))){#For every simulation...
    pdet.spatial.list[[i]] <- pdet.spatial.para.blank #Create a blank copy of probability of detection parameters in running list...
    names(pdet.spatial.list)[i] <- unique(spatial.simulation.dataset$Simulation)[i] #... and name the copy the simulation number
    
    pmatch.spatial.list[[i]] <- pmatch.spatial.para.blank #Create a blank copy of probability of mathing parameters in running list...
    names(pmatch.spatial.list)[i] <- unique(spatial.simulation.dataset$Simulation)[i] #... and name the copy the simulation number
    
    temp <- spatial.simulation.dataset[which(spatial.simulation.dataset$Simulation == unique(spatial.simulation.dataset$Simulation)[i]), ] #Isolate the summary values for this particular simulation
    for(j in 1:length(unique(spatial.simulation.dataset$Samp.Int))){ #Start nested loop through sampling intensities...
      temp2 <- temp[which(temp$Samp.Int == unique(temp$Samp.Int)[j]), ] #Isolate a sampling intensity, from the summary values for a particular simulation
      for(k in 1:length(unique(temp2$Year))){ #Start nested loop through year of simulation:sample intensity combination...
        {pdet.spatial.list[[i]]$alpha[j, 1+k] <- uninformative.prior + temp2$N.match.alive[k] #Calculate pdetection alpha shape parameter
        pdet.spatial.list[[i]]$beta[j, 1+k] <- uninformative.prior + temp2$N.inferred.alive[k] #Calculate pdetection beta shape parameter
        pmatch.spatial.list[[i]]$alpha[j, 1+k] <- uninformative.prior + temp2$N.match.alive[k]  #Calculate pmatch alpha shape parameter
        pmatch.spatial.list[[i]]$beta[j, 1+k] <- uninformative.prior + (temp2$N.Sampled.alive.ad[k] - temp2$N.match.alive[k])} #Calculate pmatch beta shape parameter
      } #Close year loop
    } #Close Sampling intensity loop
  } #close sampling iteration loop
  
  for(i in 1:nrow(spatial.simulation.dataset)){#Now, for each simulation:year:sample intensity combination...
    
    #...Isolate shape parameters for p detection and p match
    det.alpha <- pdet.spatial.list[[which(names(pdet.spatial.list) == spatial.simulation.dataset$Simulation[i])]]$alpha[which(pdet.spatial.list[[which(names(pdet.spatial.list) == spatial.simulation.dataset$Simulation[i])]]$alpha$Samp.Int == spatial.simulation.dataset$Samp.Int[i]),  as.numeric(spatial.simulation.dataset$Year[i])+1]
    det.beta <- pdet.spatial.list[[which(names(pdet.spatial.list) == spatial.simulation.dataset$Simulation[i])]]$beta[which(pdet.spatial.list[[which(names(pdet.spatial.list) == spatial.simulation.dataset$Simulation[i])]]$beta$Samp.Int == spatial.simulation.dataset$Samp.Int[i]),  as.numeric(spatial.simulation.dataset$Year[i])+1]
    match.alpha <- pmatch.spatial.list[[which(names(pmatch.spatial.list) == spatial.simulation.dataset$Simulation[i])]]$alpha[which(pmatch.spatial.list[[which(names(pmatch.spatial.list) == spatial.simulation.dataset$Simulation[i])]]$alpha$Samp.Int == spatial.simulation.dataset$Samp.Int[i]),  as.numeric(spatial.simulation.dataset$Year[i])+1]
    match.beta <- pmatch.spatial.list[[which(names(pmatch.spatial.list) == spatial.simulation.dataset$Simulation[i])]]$beta[which(pmatch.spatial.list[[which(names(pmatch.spatial.list) == spatial.simulation.dataset$Simulation[i])]]$beta$Samp.Int == spatial.simulation.dataset$Samp.Int[i]),  as.numeric(spatial.simulation.dataset$Year[i])+1]
    #Using these shape parameters, randomly sample from beta distributions to estimate detection and matching probabilities
    det <- rbeta(n  =  10000, shape1 = det.alpha,  shape2  =  det.beta)
    match  <-rbeta(n  =  10000, shape1 = match.alpha,  shape2   =  match.beta)
    
    #Calculate adult abundance by dividing the number of matched individuals by the joint probability of being detected(sampled) and matched with another individual.
    nhat.sa <- spatial.simulation.dataset$N.match.alive[i]/(det*match)
    
    #Calculate how far off this abundance estimate is from the known adult abundance in the study area (N.studyarea.ad) and for the entire population (N.ad)
    spatial.simulation.dataset$Nhat.sa[i] <- (mean(nhat.sa) - spatial.simulation.dataset$N.studyarea.ad[i])/spatial.simulation.dataset$N.studyarea.ad[i]
    spatial.simulation.dataset$Nhat[i] <- (mean(nhat.sa) - spatial.simulation.dataset$N.ad[i])/spatial.simulation.dataset$N.ad[i]
    
    #Calculate 95% and 99% bootstrap CIs, and min and max abundance estimates using pdetection and pmatch draws
    #Restricted to the study area (sa)
    spatial.simulation.dataset$LCL.sa.95[i] <- (quantile(x  =  nhat.sa,  probs  =  c(0.025)) - spatial.simulation.dataset$N.studyarea.ad[i])/spatial.simulation.dataset$N.studyarea.ad[i]
    spatial.simulation.dataset$UCL.sa.95[i] <-  (quantile(x  =  nhat.sa,  probs  =  (.975)) - spatial.simulation.dataset$N.studyarea.ad[i])/spatial.simulation.dataset$N.studyarea.ad[i]
    spatial.simulation.dataset$LCL.sa.99[i] <- (quantile(x  =  nhat.sa,  probs  =  c(0.005)) - spatial.simulation.dataset$N.studyarea.ad[i])/spatial.simulation.dataset$N.studyarea.ad[i]
    spatial.simulation.dataset$UCL.sa.99[i] <-  (quantile(x  =  nhat.sa,  probs  =  (.995)) - spatial.simulation.dataset$N.studyarea.ad[i])/spatial.simulation.dataset$N.studyarea.ad[i]
    spatial.simulation.dataset$min.sa[i] <- (min(nhat.sa) - spatial.simulation.dataset$N.studyarea.ad[i])/spatial.simulation.dataset$N.studyarea.ad[i]
    spatial.simulation.dataset$max.sa[i] <-  (max(nhat.sa) - spatial.simulation.dataset$N.studyarea.ad[i])/spatial.simulation.dataset$N.studyarea.ad[i]
    
    #Or for population-wide abundance estimates.
    spatial.simulation.dataset$LCL.95[i] <- (quantile(x  =  nhat.sa,  probs  =  c(0.025)) - spatial.simulation.dataset$N.ad[i])/spatial.simulation.dataset$N.ad[i]
    spatial.simulation.dataset$UCL.95[i] <-  (quantile(x  =  nhat.sa,  probs  =  (.975)) - spatial.simulation.dataset$N.ad[i])/spatial.simulation.dataset$N.ad[i]
    spatial.simulation.dataset$LCL.99[i] <- (quantile(x  =  nhat.sa,  probs  =  c(0.005)) - spatial.simulation.dataset$N.ad[i])/spatial.simulation.dataset$N.ad[i]
    spatial.simulation.dataset$UCL.99[i] <-  (quantile(x  =  nhat.sa,  probs  =  (.995)) - spatial.simulation.dataset$N.ad[i])/spatial.simulation.dataset$N.ad[i]
    spatial.simulation.dataset$min[i] <- (min(nhat.sa) - spatial.simulation.dataset$N.ad[i])/spatial.simulation.dataset$N.ad[i]
    spatial.simulation.dataset$max[i] <-  (max(nhat.sa) - spatial.simulation.dataset$N.ad[i])/spatial.simulation.dataset$N.ad[i]
    
  } 
  
  return(spatial.simulation.dataset)
} #End function