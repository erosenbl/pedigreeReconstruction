#Function to calculate detection and match probability shape parameters, and then calculate abundance using pedigree reconstruction.. 

#The inputs include an object storing simulations from an individually-based population model. The structure should follow that created by the simulation loop in SimulationCode_Rosenblatt_et_al.R (L201-403). It is assumed that this function will be applied to a 25 year simulation, with year 1-20 serving as a burn-in period, and years 21-25 serving as the years sampled

#The output is the same as the input, with additional columns with each simulations average abundance measure, and various confidence intervals and ranges reported.

nHat_est_population <- function(simulation.dataset) {

#First, we create storage for posterior estimates of probability of detection (pdet.para.blank) and probability of matching individuals (pmatch.para.blank), using sample data from the population-based approach. These storage ojects store beta distribution shape parameters alpha and beta.
pdet.para.blank <- list(data.frame("Samp.Int"  =  unique(simulation.dataset$Samp.Int),  
                                   "Year 1"  =  NA,  "Year 2"  =  NA,  
                                   "Year 3"  =  NA,  "Year 4"  =  NA,  
                                   "Year 5"  =  NA),  
                        data.frame("Samp.Int"  =  unique(simulation.dataset$Samp.Int),  
                                   "Year 1"  =  NA,  "Year 2"  =  NA,  
                                   "Year 3"  =  NA,  "Year 4"  =  NA,  
                                   "Year 5"  =  NA))
#Name the detection shape parameter storage
names(pdet.para.blank) <- c("alpha", "beta")


pmatch.para.blank <- list(data.frame("Samp.Int"  =  unique(simulation.dataset$Samp.Int), 
                                     "Year 1"  =  NA,  "Year 2"  =  NA,  
                                     "Year 3"  =  NA,  "Year 4"  =  NA,  
                                     "Year 5"  =  NA),  
                          data.frame("Samp.Int"  =  unique(simulation.dataset$Samp.Int), 
                                     "Year 1"  =  NA,  "Year 2"  =  NA,  
                                     "Year 3"  =  NA,  "Year 4"  =  NA,  
                                     "Year 5"  =  NA))

names(pmatch.para.blank) <- c("alpha", "beta")

#Create storage for beta distribution parameters - these lists will include copies of pdet.para and pmatch.para, each corresponding to a simulation.
pdet.list <- list(NA)
pmatch.list <- list(NA)

uninformative.prior <- 1 #Set uninformative prior for all parameter estimation (results in a uniform beta distribution)

for(i in 1:length(unique(simulation.dataset$Simulation))){ #For every simulation...
  pdet.list[[i]] <- pdet.para.blank #Create a blank copy of probability of detection parameters in running list...
  names(pdet.list)[i] <- unique(simulation.dataset$Simulation)[i] #... and name the copy the simulation number
  
  pmatch.list[[i]] <- pmatch.para.blank #Create a blank copy of probability of mathing parameters in running list...
  names(pmatch.list)[i] <- unique(simulation.dataset$Simulation)[i] #... and name the copy the simulation number
  
  temp <- simulation.dataset[which(simulation.dataset$Simulation == unique(simulation.dataset$Simulation)[i]), ] #Isolate the summary values for this particular simulation
  for(j in 1:length(unique(simulation.dataset$Samp.Int))){ #Start nested loop through sampling intensities...
    temp2 <- temp[which(temp$Samp.Int == unique(temp$Samp.Int)[j]), ] #Isolate a sampling intensity, from the summary values for a particular simulation
    for(k in 1:length(unique(temp2$Year))){ #Start nested loop through year of simulation:sample intensity combination...
      {pdet.list[[i]]$alpha[j, k+1] <- uninformative.prior + temp2$N.match.alive[k] #Calculate pdetection alpha shape parameter
      pdet.list[[i]]$beta[j, k+1] <- uninformative.prior + temp2$N.inferred.alive[k] #Calculate pdetection beta shape parameter
      pmatch.list[[i]]$alpha[j, k+1] <- uninformative.prior + temp2$N.match.alive[k]  #Calculate pmatch alpha shape parameter
      pmatch.list[[i]]$beta[j, k+1] <- uninformative.prior + (temp2$N.Sampled.alive.ad[k] - temp2$N.match.alive[k])} #Calculate pmatch beta shape parameter
    } #Close year loop
  } #Close Sampling intensity loop
} #close sampling iteration loop

for(i in 1:nrow(simulation.dataset)){ #Now, for each simulation:year:sample intensity combination...
  
  #...Isolate shape parameters for p detection and p match
  det.alpha <- pdet.list[[which(names(pdet.list) == simulation.dataset$Simulation[i])]]$alpha[which(pdet.list[[which(names(pdet.list) == simulation.dataset$Simulation[i])]]$alpha$Samp.Int == simulation.dataset$Samp.Int[i]),  as.numeric(simulation.dataset$Year[i])+1]
  det.beta <- pdet.list[[which(names(pdet.list) == simulation.dataset$Simulation[i])]]$beta[which(pdet.list[[which(names(pdet.list) == simulation.dataset$Simulation[i])]]$beta$Samp.Int == simulation.dataset$Samp.Int[i]),  as.numeric(simulation.dataset$Year[i])+1]
  match.alpha <- pmatch.list[[which(names(pmatch.list) == simulation.dataset$Simulation[i])]]$alpha[which(pmatch.list[[which(names(pmatch.list) == simulation.dataset$Simulation[i])]]$alpha$Samp.Int == simulation.dataset$Samp.Int[i]),  as.numeric(simulation.dataset$Year[i])+1]
  match.beta <- pmatch.list[[which(names(pmatch.list) == simulation.dataset$Simulation[i])]]$beta[which(pmatch.list[[which(names(pmatch.list) == simulation.dataset$Simulation[i])]]$beta$Samp.Int == simulation.dataset$Samp.Int[i]),  as.numeric(simulation.dataset$Year[i])+1]
  
  #Using these shape parameters, randomly sample from beta distributions to estimate detection and matching probabilities
  det <- rbeta(n  =  10000, shape1 = det.alpha,  shape2  =  det.beta)
  match  <-rbeta(n  =  10000, shape1 = match.alpha,  shape2   =  match.beta)
  
  #Calculate adult abundance by dividing the number of matched individuals by the joint probability of being detected(sampled) and matched with another individual.
  nhat.sa <- simulation.dataset$N.match.alive[i]/(det*match)
  
  #Calculate how far off this abundance estimate is from the known adult abundance in the study area (N.studyarea.ad) and for the entire population (N.ad)
  simulation.dataset$Nhat.sa[i] <- (mean(nhat.sa) - simulation.dataset$N.studyarea.ad[i])/simulation.dataset$N.studyarea.ad[i]
  simulation.dataset$Nhat[i] <- (mean(nhat.sa) - simulation.dataset$N.ad[i])/simulation.dataset$N.ad[i]
  
  #Calculate 95% and 99% bootstrap CIs, and min and max abundance estimates using pdetection and pmatch draws
  #Restricted to the study area (sa)
  simulation.dataset$LCL.sa.95[i] <- (quantile(x  =  nhat.sa,  probs  =  c(0.025)) - simulation.dataset$N.studyarea.ad[i])/simulation.dataset$N.studyarea.ad[i]
  simulation.dataset$UCL.sa.95[i] <-  (quantile(x  =  nhat.sa,  probs  =  (.975)) - simulation.dataset$N.studyarea.ad[i])/simulation.dataset$N.studyarea.ad[i]
  simulation.dataset$LCL.sa.99[i] <- (quantile(x  =  nhat.sa,  probs  =  c(0.005)) - simulation.dataset$N.studyarea.ad[i])/simulation.dataset$N.studyarea.ad[i]
  simulation.dataset$UCL.sa.99[i] <-  (quantile(x  =  nhat.sa,  probs  =  (.995)) - simulation.dataset$N.studyarea.ad[i])/simulation.dataset$N.studyarea.ad[i]
  simulation.dataset$min.sa[i] <- (min(nhat.sa) - simulation.dataset$N.studyarea.ad[i])/simulation.dataset$N.studyarea.ad[i]
  simulation.dataset$max.sa[i] <-  (max(nhat.sa) - simulation.dataset$N.studyarea.ad[i])/simulation.dataset$N.studyarea.ad[i]
  
  #Or for population-wide abundance estimates.
  simulation.dataset$LCL.95[i] <- (quantile(x  =  nhat.sa,  probs  =  c(0.025)) - simulation.dataset$N.ad[i])/simulation.dataset$N.ad[i]
  simulation.dataset$UCL.95[i] <-  (quantile(x  =  nhat.sa,  probs  =  (.975)) - simulation.dataset$N.ad[i])/simulation.dataset$N.ad[i]
  simulation.dataset$LCL.99[i] <- (quantile(x  =  nhat.sa,  probs  =  c(0.005)) - simulation.dataset$N.ad[i])/simulation.dataset$N.ad[i]
  simulation.dataset$UCL.99[i] <-  (quantile(x  =  nhat.sa,  probs  =  (.995)) - simulation.dataset$N.ad[i])/simulation.dataset$N.ad[i]
  simulation.dataset$min[i] <- (min(nhat.sa) - simulation.dataset$N.ad[i])/simulation.dataset$N.ad[i]
  simulation.dataset$max[i] <-  (max(nhat.sa) - simulation.dataset$N.ad[i])/simulation.dataset$N.ad[i]
} #End simulation:year:sample intensity combination loop

return(simulation.dataset)
} #End function