#Moose CRE Simulation####
#The following code simulated a spatially-explicit moose population,  and then samples individuals using two methods: Population based sampling,  or a proportion of unsampled individuals; or Spatially-based samples,  where a random proportion of grid cells are surveyed and animals using these cells are sampled within these selected grid cells. The script ends in saving various data files for precision,  bias,  and accuracy evaluations.

#Load Packages
library(plyr)
library(reshape2)
library(ggplot2)
library(popbio)
library(profvis)
library(raster)
library(rgdal)
library(dplyr)
library(rgeos)
library(MASS)
library(adehabitatHR)
library(REdaS)
library(lubridate)
library(rlist)
library(tidyverse)
library(ggplot2)
library(sp)

#load custom functions
source("Functions/gOverlap.R") #Calculates degree of overlap between potential mates
source("Functions/demo_samp_updated_withiteration.R") #Function to sample proportion of moose on landscape (not spatially-explicit)
source("Functions/cre_est_updated.R") #Function to calculate various pedigree reconstruction estimates
source("Functions/reproduce.R") #Function for moose reproduction
source("Functions/spatial.samp.faster.R")
source("Functions/cre_est_updated_spatial.R")

#set seed for simulation processes to be repeatable
set.seed(1234)

#Step 1: Build Simulation----
##Part A: Set Spatial Extent of the simulation####
###Pull in shapefile files for study area where simulation occurs
moose.study.area <- shapefile(x  =  "data/Study_Area/StudyDissolve.shp")  #Load study area
proj4string <- CRS("+proj=utm +zone=19 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")  #Set projection of study area
moose.study.area <- spTransform(moose.study.area,  CRSobj = proj4string,  inverse = TRUE) #Transform shapefile projection

study_area_utm <- sp::CRS("+init=epsg:32619") %>% 
  spTransform(moose.study.area,  .) #ensure CRS is correct for study area

###Develop square sampling grid with 1 sq km grid cells,  that can be used for spatial sampling later on in this simulation
grid <- raster(extent(study_area_utm),  resolution  =  c(1000, 1000),  crs  =  proj4string(study_area_utm))
grid <- raster::extend(grid,  c(1, 1))
gridPolygon <- rasterToPolygons(grid)
gridPoints <- rasterToPoints(grid,  spatial = TRUE)
sample_grid <- SpatialPixels(gridPoints,  proj4string  =  proj4string(study_area_utm))
pixel_in_SA <- sapply(over(sample_grid,  study_area_utm,  returnList  =  TRUE), nrow)
pixel_in_SA <- data.frame(Pixel  =  sample_grid@grid.index,  In.SA  =  pixel_in_SA) #Create pixel data

sampleGrid.df <- SpatialPixelsDataFrame(points  =  gridPoints,  data  =  pixel_in_SA,  proj4string  =  proj4string(study_area_utm)) #Combine spatial pixels and dataframe into a SpatialPixelsDataFrame

##Part B: Initial Population Demography####
Nstart <- 659 #Set up starting population

stable.age <- c(1.083016e-01, 7.178148e-02, 6.320842e-02, 5.565926e-02, 4.901173e-02, 4.315812e-02, 3.800363e-02,  3.346475e-02,  2.946796e-02,  2.594852e-02,  2.284942e-02,  1.874704e-02,  1.380596e-02,  8.652304e-03,  4.293857e-03,  2.011853e-03,  1.083016e-01,  6.800351e-02,  5.367999e-02,  4.237343e-02,  3.344836e-02,  2.640316e-02,  2.084189e-02,  1.645198e-02,  1.298672e-02,  1.025134e-02,  8.092109e-03,  5.621371e-03,  3.233268e-03,  1.420961e-03,  4.327576e-04,  9.169111e-05) #Stable age distribution for the vital rates used in this study (calculated using popbio package),  with proportions for 16 ages for each sex. 

N0 <- matrix(round(Nstart*stable.age),  ncol  =  1)  #Create an initial vector of population demography,  by multiplying the stable age distribution above,  by the initial population size.
N0[31:32] <- c(1, 1) #Adding two old males to match starting population size and correct for rounding

#Initial population matrix,  based on the number of initial individuals and sex ratio,  parentage unknown
pop.0 <- data.frame("Individual" = c(1:sum(N0)),  "Male" = c(rep(0, sum(N0[1:16])), rep(1, sum(N0[17:32]))),  "Age Class"  =  c(rep(0, N0[1]), rep(1, N0[2]), rep(2, N0[3]), rep(3, N0[4]), rep(4, N0[5]), rep(5, N0[6]),  rep(6, N0[7]),  rep(7, N0[8]), rep(8, N0[9]), rep(9, N0[10]), rep(10, N0[11]), rep(11, N0[12]),  rep(12, N0[13]), rep(13, N0[14]), rep(14, N0[15]), rep(15, N0[16]), rep(0, N0[17]), rep(1, N0[18]), rep(2, N0[19]), rep(3, N0[20]), rep(4, N0[21]), rep(5, N0[22]),  rep(6, N0[23]),  rep(7, N0[24]), rep(8, N0[25]), rep(9, N0[26]), rep(10, N0[27]), rep(11, N0[28]),  rep(12, N0[29]), rep(13, N0[30]), rep(14, N0[31]), rep(15, N0[32])),  "Mother" = NA,  "Father" = NA,  "Sample.P" = rep(NA, sum(N0)),  "Sample" = rep(0, sum(N0)),  "Spatial.Sample"  =  rep(0, sum(N0)))

pop.0$Sex <- "Female"
pop.0$Sex[which(pop.0$Male == 1)] <- "Male"
pop.0$Sex <- as.factor(pop.0$Sex)

pop.0$Age.Cat <- ifelse(pop.0$Age.Class == 0,  "Calf",  ifelse(pop.0$Age.Class == 1,  "Yearling",  "Adult"))

##Part C: Spatial data for each moose####
#Establish home ranges - centroids and distributions based on annual home ranges . No known mother calf relationships,  so no need to match them.
pop.0$mu <- coordinates(spsample(moose.study.area, n = nrow(pop.0), "random"))# #Randomly generate centroids for all animals

#Read in point data and home range data
  points <- readRDS("data/moose_gps_points.RDS") #Individual location points for building kernel utilization distribution from real, collared moose (cleaned to only include relevant data for this analysis)
  sigma.list <- readRDS("data/moose_HR.RDS") #Starting point values for moose in initial population
  rel.points <- readRDS("data/moose_relative_locations.RDS") #GPS coordinates for locations relative to home range centroid (to generate new locations and home ranges for individuals added to the population)
  
##Part D. Vital Rates####
#Survival Rates,  with senesence for moose 10 years old and up.
  s.summer.calf <- 0.7
  s.summer.adult.female<-0.95
  s.summer.adult.male<-0.9
  
  s.summer.10yo.female<-1-((((1-s.summer.adult.female)/s.summer.adult.female)*1.6)/(1+(((1-s.summer.adult.female)/s.summer.adult.female)*1.6)))
  s.summer.10yo.male<-1-((((1-s.summer.adult.male)/s.summer.adult.male)*1.6)/(1+(((1-s.summer.adult.male)/s.summer.adult.male)*1.6)))
  
  s.summer.11yo.female<-1-((((1-s.summer.10yo.female)/s.summer.10yo.female)*1.6)/(1+(((1-s.summer.10yo.female)/s.summer.10yo.female)*1.6)))
  s.summer.11yo.male<-1-((((1-s.summer.10yo.male)/s.summer.10yo.male)*1.6)/(1+(((1-s.summer.10yo.male)/s.summer.10yo.male)*1.6)))
  
  s.summer.12yo.female<-1-((((1-s.summer.11yo.female)/s.summer.11yo.female)*1.6)/(1+(((1-s.summer.11yo.female)/s.summer.11yo.female)*1.6)))
  s.summer.12yo.male<-1-((((1-s.summer.11yo.male)/s.summer.11yo.male)*1.6)/(1+(((1-s.summer.11yo.male)/s.summer.11yo.male)*1.6)))
  
  s.summer.13yo.female<-1-((((1-s.summer.12yo.female)/s.summer.12yo.female)*1.6)/(1+(((1-s.summer.12yo.female)/s.summer.12yo.female)*1.6)))
  s.summer.13yo.male<-1-((((1-s.summer.12yo.male)/s.summer.12yo.male)*1.6)/(1+(((1-s.summer.12yo.male)/s.summer.12yo.male)*1.6)))
  
  s.summer.14yo.female<-1-((((1-s.summer.13yo.female)/s.summer.13yo.female)*1.6)/(1+(((1-s.summer.13yo.female)/s.summer.13yo.female)*1.6)))
  s.summer.14yo.male<-1-((((1-s.summer.13yo.male)/s.summer.13yo.male)*1.6)/(1+(((1-s.summer.13yo.male)/s.summer.13yo.male)*1.6)))
  
  s.summer.15yo.female<-1-((((1-s.summer.14yo.female)/s.summer.14yo.female)*1.6)/(1+(((1-s.summer.14yo.female)/s.summer.14yo.female)*1.6)))
  s.summer.15yo.male<-1-((((1-s.summer.14yo.male)/s.summer.14yo.male)*1.6)/(1+(((1-s.summer.14yo.male)/s.summer.14yo.male)*1.6)))
  
  
  s.winter.calf <- 0.7
  s.winter.adult.female<-0.93
  s.winter.adult.male<-0.88
  
  s.winter.10yo.female<-1-((((1-s.winter.adult.female)/s.winter.adult.female)*1.6)/(1+(((1-s.winter.adult.female)/s.winter.adult.female)*1.6)))
  s.winter.10yo.male<-1-((((1-s.winter.adult.male)/s.winter.adult.male)*1.6)/(1+(((1-s.winter.adult.male)/s.winter.adult.male)*1.6)))
  
  s.winter.11yo.female<-1-((((1-s.winter.10yo.female)/s.winter.10yo.female)*1.6)/(1+(((1-s.winter.10yo.female)/s.winter.10yo.female)*1.6)))
  s.winter.11yo.male<-1-((((1-s.winter.10yo.male)/s.winter.10yo.male)*1.6)/(1+(((1-s.winter.10yo.male)/s.winter.10yo.male)*1.6)))
  
  
  s.winter.12yo.female<-1-((((1-s.winter.11yo.female)/s.winter.11yo.female)*1.6)/(1+(((1-s.winter.11yo.female)/s.winter.11yo.female)*1.6)))
  s.winter.12yo.male<-1-((((1-s.winter.11yo.male)/s.winter.11yo.male)*1.6)/(1+(((1-s.winter.11yo.male)/s.winter.11yo.male)*1.6)))
  
  
  s.winter.13yo.female<-1-((((1-s.winter.12yo.female)/s.winter.12yo.female)*1.6)/(1+(((1-s.winter.12yo.female)/s.winter.12yo.female)*1.6)))
  s.winter.13yo.male<-1-((((1-s.winter.12yo.male)/s.winter.12yo.male)*1.6)/(1+(((1-s.winter.12yo.male)/s.winter.12yo.male)*1.6)))
  
  s.winter.14yo.female<-1-((((1-s.winter.13yo.female)/s.winter.13yo.female)*1.6)/(1+(((1-s.winter.13yo.female)/s.winter.13yo.female)*1.6)))
  s.winter.14yo.male<-1-((((1-s.winter.13yo.male)/s.winter.13yo.male)*1.6)/(1+(((1-s.winter.13yo.male)/s.winter.13yo.male)*1.6)))
  
  s.winter.15yo.female<-1-((((1-s.winter.14yo.female)/s.winter.14yo.female)*1.6)/(1+(((1-s.winter.14yo.female)/s.winter.14yo.female)*1.6)))
  s.winter.15yo.male<-1-((((1-s.winter.14yo.male)/s.winter.14yo.male)*1.6)/(1+(((1-s.winter.14yo.male)/s.winter.14yo.male)*1.6)))
  

#Set up survival matrices for winter and summer,  with the number of rows matching the population matrix
S1<-matrix(c(s.winter.calf, s.winter.adult.female, s.winter.adult.female, s.winter.adult.female, s.winter.adult.female, s.winter.adult.female, s.winter.adult.female, s.winter.adult.female, s.winter.adult.female, s.winter.adult.female, s.winter.10yo.female,  s.winter.11yo.female, s.winter.12yo.female, s.winter.13yo.female, s.winter.14yo.female, s.winter.15yo.female, 
             s.winter.calf, s.winter.adult.male, s.winter.adult.male, s.winter.adult.male, s.winter.adult.male, s.winter.adult.male, s.winter.adult.male, s.winter.adult.male, s.winter.adult.male, s.winter.adult.male, s.winter.10yo.male, s.winter.11yo.male, s.winter.12yo.male, s.winter.13yo.male, s.winter.14yo.male, s.winter.15yo.male))

S2<-matrix(c(s.summer.calf, s.summer.adult.female, s.summer.adult.female, s.summer.adult.female, s.summer.adult.female, s.summer.adult.female, s.summer.adult.female, s.summer.adult.female, s.summer.adult.female, s.summer.adult.female, s.summer.10yo.female,  s.summer.11yo.female, s.summer.12yo.female, s.summer.13yo.female, s.summer.14yo.female, s.summer.15yo.female, 
             s.summer.calf, s.summer.adult.male, s.summer.adult.male, s.summer.adult.male, s.summer.adult.male, s.summer.adult.male, s.summer.adult.male, s.summer.adult.male, s.summer.adult.male, s.summer.adult.male, s.summer.10yo.male, s.summer.11yo.male, s.summer.12yo.male, s.summer.13yo.male, s.summer.14yo.male, s.summer.15yo.male))


##Reproduction---- Based on Stable Population
#Fecundity (Calves per mother)
f.calf <- 0
f.year <- 0.072
f.2year <- 0.765
f.adult <- 0.9
f.sen <- 0

sex.ratio<-0.5

#Set up reproductive matrix
F<-matrix(c(f.calf, f.year, f.2year, f.adult,  f.sen))

#Step 2: Run Simulation####
#There is an annual time step for this simulation,  with the following order of steps. The result will be a a number of objects saved that capture the ability of pedigree reconstruction to estimate abundance.

#Basic Stages:
#1. Identify candidate fathers from fall prior when breeding occurs
#2. Animals survive winter (Nov-April)
#3. Animals reproduce (May)
#4. Animals Age,  Disperse,  and Establish Summer HR
#5. Animals survive summer (May-October)
#6. Animals harvested (October)
#7. Estimate population size (Fall Sampling Only)
#8. Winter Home Ranges (HRs) Established (Starting condition for next iteration)

##Part A: Simulation Conditions ####
scope <- 25 #duration of simulation (time steps; years)
samp.start <- 21 #Time step (year) when sampling starts
sims <- 100 #number of simulated sampling events per year
sample <- seq(0.1, 0.9, 0.1) #Sampling Intensity (10-90% of unsampled individuals)
spatial.sample <- seq(0.0061, 0.122, length.out = 20) #set spatial sampling intensity (proportion of study area)
season <- "winter" #when sampling occurs
sample.SA.only <- TRUE #Should sampling be restricted to study area only?
age.min  =  1 #Age minimum for inference

pop.t<-pop.0 #set pop.t at initial population size and distribution
last.ID = tail(pop.0$Individual, 1)  #Track last used animal ID number for assigning new IDs to newborne calves throughout the simulation.
sigma.list.t <- sigma.list #Starting point values for moose in initial population

points.t <- points #Starting locations for all animals in the initial population

##Part B: Create Storage Objects for population size,  sampled individuals,  and population size calculation----
pop.traj<-matrix(nrow = 1, ncol = scope) #Track changes in abundance through th simulation

sim.samp.list <- list(NA) #List to store sampled individuals from each simulation,  using the population-based approach

sim.spatial.samp.list <- list(NA)#List to store sampled individuals from each simulation,  using the spatially-based approach.

Nhat.calc <- data.frame(Year = NA, Simulation = NA, Samp.Int = NA,  Real.Samp.Int  =  NA,  StudyAreaEst  =  NA, N.sampled.alive = NA, N.Sampled.alive.ad = NA,  N.match  =  NA,  N.match.alive  =  NA,  N.unmatch.alive  =  NA, N.inferred.alive = NA,  LCL.sa = NA,  Nhat.sa = NA,  UCL.sa = NA,  LCL.Ndet.sa = NA,  Nhat.Ndet.sa = NA,  UCL.Ndet.sa = NA,  LCL.Ndet.all.sa = NA,  Nhat.Ndet.all.sa = NA,  UCL.Ndet.all.sa = NA,  LCL = NA,  Nhat = NA,  UCL = NA,  LCL.Ndet = NA,  Nhat.Ndet = NA,  UCL.Ndet = NA,  LCL.Ndet.all = NA,  Nhat.Ndet.all = NA,  UCL.Ndet.all = NA,  N = NA,  N.ad = NA,  N.studyarea = NA,  N.studyarea.ad  =  NA) #Blank template for recording abundance estimates using population-based approach


Nhat.calc.spatial<-data.frame(Year = NA, Simulation = NA, Samp.Int = NA,  Real.Samp.Int  =  NA,  StudyAreaEst  =  NA, N.sampled.alive = NA, N.Sampled.alive.ad = NA,  N.match  =  NA,  N.match.alive  =  NA,  N.unmatch.alive  =  NA, N.inferred.alive = NA,  LCL.sa = NA,  Nhat.sa = NA,  UCL.sa = NA,  LCL.Ndet.sa = NA,  Nhat.Ndet.sa = NA,  UCL.Ndet.sa = NA,  LCL.Ndet.all.sa = NA,  Nhat.Ndet.all.sa = NA,  UCL.Ndet.all.sa = NA,  LCL = NA,  Nhat = NA,  UCL = NA,  LCL.Ndet = NA,  Nhat.Ndet = NA,  UCL.Ndet = NA,  LCL.Ndet.all = NA,  Nhat.Ndet.all = NA,  UCL.Ndet.all = NA,  N = NA,  N.ad = NA,  N.studyarea = NA,  N.studyarea.ad  =  NA) #Blank template for recording abundance estimates using spatial-based approach

#Save blank copies of these templates, as Nhat.calc and Nhat.calc.spatial will be modified.
Nhat.calc.blank <- Nhat.calc 
Nhat.calc.spatial.blank<-Nhat.calc.spatial

##Part C: Simulation####
ptm <- proc.time() #Start processing time (Warning: This simulation takes several hours to run!)
set.seed(12345) #Set seed for consistent results

for(z in 1:scope){ #Begin loop through duration of simulation (years)
  ###Part C.i: Record population size and identify candidate fathers-----
  pop.traj[z]<-nrow(pop.t) #Record abundance at beginning of annual time step
  cand.father.t <- pop.t[which(pop.t$Male == 1 & pop.t$Age.Class>= 2), c("Individual", "Age.Class")] #   Candidate fathers are listed prior to potential winter mortality 
  

  ##Part C.ii: Sample!-----
  if (season  ==  "winter"){ #If sampling occurs in the winter...
    if (z<samp.start) #.. AND if the time step is still in the burn in (earlier than the first sampling)...
    {Nhat.calc.year<-Nhat.calc.blank #keep population-based estimate df blank
    Nhat.calc.year.spatial<-Nhat.calc.spatial.blank} #and keep spatial based estimate df blank.
    
    if (z == samp.start) #If the time step is the first year of sampling, the following happens:
      
      #Population-based sampling...
    {sim.samp.list <- demo.samp.updated.withiteration(begin.sample  =  TRUE,  sim.samp.list  =  NULL,  z  =  z,  sims  =  sims,  exclude.outliers  =  sample.SA.only) #Sample through a vector of sampling efforts (sample object defined line 166), using the function demo.samp.updated.withiteration() to genetically a proportion of the population.
    Nhat.calc <- cre.est.updated(Nhat.calc  =  Nhat.calc,  Nhat.calc.year  =  Nhat.calc.year,  sim.samp.list  =  sim.samp.list,  z  =  z,  sims  =  sims,  study.area.est  =  sample.SA.only,  age.min  =  age.min) #Feed storage for abundance estimates (Nhat.calc) and sample inventory (sim.samp.list) into the cre.est.updated function, that will track sampled and matched individuals for each simulation.
    
      #Spatially-based sampling...
    sim.spatial.samp.list <- spatial.samp.faster(begin.sample  =  TRUE,  points  =  points.t,  sample.grid  =  sampleGrid.df,  return.plot  =  FALSE,  stratified  =  TRUE) #Sample through a vector of sampling efforts for area surveyed(spatial.sample object defined line 167), using the function spatial.samp.faster() to select grid cells to survey, and the resulting individuals sampled based on their utilization of surveyed grid cells.
    Nhat.calc.spatial <- cre.est.updated.spatial(Nhat.calc.spatial,  Nhat.calc.year.spatial,  z  =  z,  sims  =  sims,  sim.spatial.samp.list  =  sim.spatial.samp.list,  study.area.est  =  sample.SA.only,  age.min  =  1) #Feed storage for abundance estimates (Nhat.calc.spatial) and sample inventory (sim.spatial.samp.list) into the cre.est.updated function, that will track sampled and matched individuals for each simulation.
    } #End first year of sampling if statement
    
    if (z>samp.start) #If the time step is after the initial sampling year
    
      #Population-based sampling...
    {sim.samp.list <- demo.samp.updated.withiteration(begin.sample  =  FALSE,  sim.samp.list  =  sim.samp.list,  z  =  z,  sims  =  sims,  exclude.outliers  =  sample.SA.only)#Sample through a vector of sampling efforts (sample object defined line 166), using the function demo.samp.updated.withiteration() to genetically a proportion of the population, accounting for previous years that were sampled.
    Nhat.calc <- cre.est.updated(Nhat.calc  =  Nhat.calc,  Nhat.calc.year  =  Nhat.calc.year,  sim.samp.list  =  sim.samp.list,  z  =  z,  sims  =  sims,  study.area.est  =  sample.SA.only,  age.min  =  age.min)#Feed storage for abundance estimates (Nhat.calc) and sample inventory (sim.samp.list) into the cre.est.updated function, that will track sampled and matched individuals for each simulation, adding to those from previous years of sampling
    
    #Spatially-based sampling...
    sim.spatial.samp.list <- spatial.samp.faster(begin.sample  =  FALSE,  points  =  points.t,  sample.grid  =  sampleGrid.df,  return.plot  =  FALSE,  stratified  =  TRUE)#Sample through a vector of sampling efforts for area surveyed(spatial.sample object defined line 167), using the function spatial.samp.faster() to select grid cells to survey, and the resulting individuals sampled based on their utilization of surveyed grid cells.
    Nhat.calc.spatial <- cre.est.updated.spatial(Nhat.calc.spatial,  Nhat.calc.year.spatial,  z  =  z,  sims  =  sims,  sim.spatial.samp.list  =  sim.spatial.samp.list,  study.area.est  =  sample.SA.only,  age.min  =  1)#Feed storage for abundance estimates (Nhat.calc.spatial) and sample inventory (sim.spatial.samp.list) into the cre.est.updated function, that will track sampled and matched individuals for each simulation.
    } #End if statement specifying that this is not the first year of sampling.
    
    if(z == scope){ #if this is the last year of sampling...
      pop.final <- pop.t #record population as the final population
      points.final <- points.t #record individuals' locations as the final spatial distribution of the population.
    }
  } #end of winter sampling if statement  
  
##Part C.iii: Winter survival----
#Females experience senescence after 13 YO,  males at 10 YO
pop.t$Survive.Winter<-rep(NA, nrow(pop.t))

#Female winter survival fate by age, drawn from survival rates stored in S1
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 0), 'Survive.Winter'] <- rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 0)), 1, S1[1])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 1), 'Survive.Winter'] <- rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 1)), 1, S1[2])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 2), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 2)), 1, S1[3])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 3), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 3)), 1, S1[4])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 4), 'Survive.Winter'] <- rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 4)), 1, S1[5])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 5), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 5)), 1, S1[6])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 6), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 6)), 1, S1[7])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 7), 'Survive.Winter'] <- rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 7)), 1, S1[8])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 8), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 8)), 1, S1[9])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 9), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 9)), 1, S1[10])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 10), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 10)), 1, S1[11])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 11), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 11)), 1, S1[12])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 12), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 12)), 1, S1[13])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 13), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 13)), 1, S1[14])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 14), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 14)), 1, S1[15])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class >= 15), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class >= 15)), 1, S1[16])

#Male winter survival fate by age, drawn from survival rates stored in S1
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 0), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 0)), 1, S1[17])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 1), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 1)), 1, S1[18])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 2), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 2)), 1, S1[19])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 3), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 3)), 1, S1[20])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 4), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 4)), 1, S1[21])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 5), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 5)), 1, S1[22])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 6), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 6)), 1, S1[23])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 7), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 7)), 1, S1[24])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 8), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 8)), 1, S1[25])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 9), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 9)), 1, S1[26])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 10), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 10)), 1, S1[27])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 11), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 11)), 1, S1[28])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 12), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 12)), 1, S1[29])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 13), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 13)), 1, S1[30])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 14), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 14)), 1, S1[31])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class >= 15), 'Survive.Winter']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class >= 15)), 1, S1[32])

pop.t<-subset(pop.t,  pop.t$Survive.Winter == 1) #Keep only surviving individuals in the population
pop.t$Survive.Winter<-NULL #Remove winter survival field from the population inventory (to be rewritten every time step)

##Part C.iv: Spring recruitment----
#Use the reproduce() function to take the population (pop.t), a dataframe of potential fathers (cand.father.t), a SpatialPointsDataFrame of locations for all living animals, and fecundity information (F), simulate reporduction and add new individuals to the population, along with the identity of their parents. All other individuals celebrate their birthday.
reproduction <- reproduce(population  =  pop.t,  potential.fathers  =  cand.father.t,  locations  =  points.t,  year  =  z,  fecundity  =  F,  max.age  =  14) 
pop.t <- reproduction$population #Update the population to include the newborne animals

offspring.t <- reproduction$neonates #Record the new animals seperately as well
mates <- reproduction$mates #Record mates
last.ID = tail(pop.t$Individual, 1) #and record the last animal ID used for next year's cohort of newbornes.

##Part C.v: Spring Dispersal and updated locations----
#First, drop spatial data for animals that perished. We've kept the data for fathers that died after mating because they still need to be linked to a particular mate and offspirng.
sigma.list.t[which(!names(sigma.list.t) %in% pop.t$Individual)] <- NULL

##Second, yearlings disperse:
disperse <- pop.t[which(pop.t$Age.Class  ==  1), ] #Subset yearlings
disperse$bearing <- runif(n  =  nrow(disperse), 0, 359) #randomly draw a dispersal direction
disperse$dist <- NA #Create space for dispersal distance

#Generate sex-specific dispersal distances...
disperse$dist[which(disperse$Male  ==  0)] <- rnorm(n  =  length(which(disperse$Male  ==  0)),  mean  =  2300,  sd  =  400) #for females
disperse$dist[which(disperse$Male  ==  1)] <- rnorm(n  =  length(which(disperse$Male  ==  1)),  mean  =  9100,  sd  =  3100) #and for males

disperse[, 'mu'] <- cbind(x = disperse$mu[, 'x']+sin(deg2rad(disperse[, 'bearing']))*disperse[, 'dist'], y = disperse$mu[, 'y']+cos(deg2rad(disperse[, 'bearing']))*disperse[, 'dist']) #Calculate new centroid using dispersal bearing and distance

disperse$bearing <- NULL  #Clear bearing to match pop.t columns
disperse$dist <- NULL #Clear distance to match pop.t columns
pop.t[which(pop.t$Age.Class  ==  1), ] <- disperse #Update dispersed individuals in pop.t


#Third, assign calves with their mother's HR centroids and HRs
for(i in 1:nrow(pop.t)){  #Going by individual
  if(pop.t$Age.Class[i] == 0) #If a calf...
  {pop.t[i, 'mu'] <- pop.t[which(pop.t[, 'Individual'] == pop.t[i, 'Mother']), 'mu']} #... assign mother's centroid
  if(pop.t$Age.Class[i]>0) next #If not a calf,  skip.
}
  
  list.ref<- length(sigma.list.t) #Count of individuals with spatial data (yet to add spatial data for calves)
for (i in 1:nrow(offspring.t)){ #Add spatial data for calves in a loop...
  sigma.list.t[[list.ref+i]] <-  sigma.list.t[[which(names(sigma.list.t) == offspring.t$Mother[i])]] #... by assigning mom's location data to appropriate calf and adding said calf as a new object in sigma.list.t
} #End loop adding spatial data for calves
  
  names(sigma.list.t) <- as.character(pop.t$Individual) #Reassign names for the spatial list from pop.t, which will have the same order as the list of locations (sigma.list.t)

#Finally, update locations for individuals that have entered new age classes with the reproduction() function by drawing relative point data matrices from VT collar data.

for(i in 1:nrow(pop.t)){ #Working through all the individuals in the population
  
  if(pop.t$Age.Class[i] == 0 | pop.t$Age.Class[i]>3) next #Skip newbornes (just assigned spatial data) and established adults (who stay in their established ranges)
  
  if(pop.t$Age.Class[i] == 1) #If the individual has aged into the yearling class...
  {moose <- eval(parse(text  =  paste0("rel.points[['", pop.t$Sex[i], "']]$'Yearling'")))[[sample(1:length(eval(parse(text  =  paste0("rel.points[['", pop.t$Sex[i], "']]$'Yearling'")))),  size  =  1,  replace  =  T)]] #... randomly sample a yearling moose in the rel.points object from actual collar data...
  Sigma <- sample(x  =  moose,  size = 1,  replace = T)} #... and assign these relative points to the new yearling.
  
  if(pop.t$Age.Class[i] == 2) #If the individual has aged into the young adult class...
  {moose <- eval(parse(text  =  paste0("rel.points[['", pop.t$Sex[i], "']]$'Adult'")))[[sample(1:length(eval(parse(text  =  paste0("rel.points[['", pop.t$Sex[i], "']]$'Adult'")))),  size  =  1,  replace  =  T)]]#... randomly sample a young adult moose in the rel.points object from actual collar data...
  Sigma <- sample(x  =  moose,  size = 1,  replace = T)}#... and assign these relative points to the new young adult.
  
  if(pop.t$Age.Class[i] == 3) #If the individual has aged into the adult class...
  {moose <- eval(parse(text  =  paste0("rel.points[['", pop.t$Sex[i], "']]$'Adult'")))[[sample(1:length(eval(parse(text  =  paste0("rel.points[['", pop.t$Sex[i], "']]$'Adult'")))),  size  =  1,  replace  =  T)]]#... randomly sample an adult moose in the rel.points object from actual collar data...
  Sigma <- sample(x  =  moose,  size = 1,  replace = T)}#... and assign these relative points to the new young adult.
  
  #Take the newly assigned relative collar locations and add it to the individuals location data set in sigma.list.t.
  sigma.list.t[[i]] <- Sigma[[1]]
  
  #Calculate coordinates using these relative locations and the individual's centroid
  sigma.list.t[[i]]$Easting <- sigma.list.t[[i]]$Easting + pop.t$mu[i, 'x']
  sigma.list.t[[i]]$Northing <- sigma.list.t[[i]]$Northing + pop.t$mu[i, 'y']
}
  
  comb.points <- bind_rows(sigma.list.t,  .id  =  'id')  #Replace location data in sigma.list.t for these animals that have changed age classes
  points.t <- sp::SpatialPointsDataFrame(coords  =  comb.points[c('Easting',  'Northing')], data  =  comb.points,  proj4string  =  CRS("+init=epsg:32619")) #Update points.t (alternative formatted object, but contains the same information as sigma.list.t)
  
##Part C.vi: Summer Mortality----
#Here individuals survive based on their age-class mortality rates (object)
pop.t$Survive.Summer<-rep(NA, nrow(pop.t))
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 0), 'Survive.Summer'] <- rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 0)), 1, S2[1])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 1), 'Survive.Summer'] <- rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 1)), 1, S2[2])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 2), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 2)), 1, S2[3])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 3), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 3)), 1, S2[4])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 4), 'Survive.Summer'] <- rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 4)), 1, S2[5])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 5), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 5)), 1, S2[6])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 6), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 6)), 1, S2[7])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 7), 'Survive.Summer'] <- rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 7)), 1, S2[8])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 8), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 8)), 1, S2[9])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 9), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 9)), 1, S2[10])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 10), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 10)), 1, S2[11])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 11), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 11)), 1, S2[12])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 12), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 12)), 1, S2[13])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 13), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 13)), 1, S2[14])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 14), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 14)), 1, S2[15])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class >= 15), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class >= 15)), 1, S2[16])


#Male summer survival rates
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 0), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 0)), 1, S2[17])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 1), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 1)), 1, S2[18])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 2), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 2)), 1, S2[19])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 3), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 3)), 1, S2[20])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 4), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 4)), 1, S2[21])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 5), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 5)), 1, S2[22])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 6), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 6)), 1, S2[23])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 7), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 7)), 1, S2[24])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 8), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 8)), 1, S2[25])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 9), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 9)), 1, S2[26])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 10), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 10)), 1, S2[27])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 11), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 11)), 1, S2[28])
pop.t[which(pop.t$Male == 1 & pop.t$Age.Class == 12), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 1 & pop.t$Age.Class == 12)), 1, S2[29])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 13), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 13)), 1, S2[30])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class == 14), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class == 14)), 1, S2[31])
pop.t[which(pop.t$Male == 0 & pop.t$Age.Class >= 15), 'Survive.Summer']<-rbinom(length(which(pop.t$Male == 0 & pop.t$Age.Class >= 15)), 1, S2[32])

pop.t<-subset(pop.t,  pop.t$Survive.Summer == 1) #Filter out animals that perished during the summer
pop.t$Survive.Summer<-NULL #Remove field for summer survival fate
sigma.list.t[which(!names(sigma.list.t) %in% pop.t$Individual)] <- NULL #Remove home ranges for animals that perished breeding occurs (start of next year of simulation).
comb.points <- bind_rows(sigma.list.t,  .id  =  'id') #Update point data...
points.t <- sp::SpatialPointsDataFrame(coords  =  comb.points[c('Easting',  'Northing')], data  =  comb.points,  proj4string  =  CRS("+init=epsg:32619")) #... and overwrite points.t object that contains the SpatialPointsDataFrame.

print(paste0("Year ", z, " complete")) #Flash the year completed on the console
} #End of year loop
proc.time() - ptm #Report elapsed time.

#Step 3: Abundance estimation from sampled animals and detected and inferred relationships----

#Run the nHat_est_population() function to calculate abundance and various CIs for population-based sampling
Nhat.calc <- nHat_est_population(Nhat.calc)

#Run the nHat_est_spatial() function to calculate abundance and various CIs for spatially-based sampling
Nhat.calc.spatial <- nHat_est_population(Nhat.calc.spatial)


#Step 4: Calculate CV, SME, and SRMSE----
#Population-based sampling
Nhat.calc %>% 
  mutate(., Nhat.ad = (1+Nhat) * N.ad) %>% 
  dplyr::select(., Year, Real.Samp.Int, Samp.Int, Nhat.ad, N.ad, N.studyarea.ad) %>% 
  group_by(., Samp.Int, Year) %>% 
  summarize(StD.ad = sd(Nhat.ad),
            CV.ad = StD.ad/mean(Nhat.ad),
            SME.ad.sa = mean((Nhat.ad) - N.studyarea.ad)/mean(N.studyarea.ad),
            SRMSE.sa.ad = (1/mean(N.studyarea.ad))*sqrt(sum((Nhat.ad-mean(N.studyarea.ad))^2)/length(Samp.Int))) %>% 
  select(., -StD.ad)-> results.summary

#Spatially-based sampling
Nhat.calc.spatial %>% 
  mutate(., Nhat.ad = (1+Nhat) * N.ad) %>% 
  mutate(Samp.Int = as.factor(round(Samp.Int*1650))) %>% 
  filter(., Samp.Int %in% c("10","30","50","70", "91","111","131", "151","171","191","201")) %>% 
  dplyr::select(., Year, Real.Samp.Int, Samp.Int, Nhat.ad, N.ad, N.studyarea.ad) %>% 
  group_by(., Samp.Int, Year) %>% 
  summarize(StD.ad = sd(Nhat.ad),
            CV.ad = StD.ad/mean(Nhat.ad),
            SME.ad = mean((Nhat.ad) - N.ad)/mean(N.ad),
            SME.ad.sa = mean((Nhat.ad) - N.studyarea.ad)/mean(N.studyarea.ad),
            SRMSE.ad = (1/mean(N.ad))*sqrt(sum((Nhat.ad-mean(N.ad))^2)/length(Samp.Int)),
            SRMSE.sa.ad = (1/mean(N.studyarea.ad))*sqrt(sum((Nhat.ad-mean(N.studyarea.ad))^2)/length(Samp.Int))) -> results.summary.spatial

