#Function to simulate spatially explicit reproduction

#The inputs include a dataframe listing all individuals currently living in the population (population), a dataframe of potential fathers and their age (potential.fathers), a SpatialPointsDataFrame listing all individuals locations (locations), year of the simulation (year; same as z in other functions), a reproductive matrix listing female fecundity corresponding to age classes (fecundity; see function code), and the maximum age of reproduction (max.age).

#The outputs include: a updated population list that includes the new offspring (population), a list of newborn animals with covariate data and the identity of their parents (neonates), and a matrix listing the degree of home range overlap for males and females for assigning a male mate to a female that gives birth to a calf.

reproduce <- function(population, potential.fathers, locations, year, fecundity, max.age) {
  cand.mother.t <- data.frame("Individual"=population[which(population$Male==0 & population$Age.Class>= 1 & population$Age.Class<=max.age),1], "Age.Class"=population[which(population$Male==0 & population$Age.Class>=1 &  population$Age.Class<=max.age),3], "Mate"=NA, "Offspring"=NA)   #Creates list of all potential mothers in a population, ages 1 - max.age
  
  cand.mother.t[which(cand.mother.t$Age.Class == 1), 'Offspring'] <- rbinom(n = length(which(cand.mother.t$Age.Class == 1)), size = 1, prob = fecundity[2]) #Randomly draws 0/1 calving fate for all yearlings, based on yearling fecundity value.
  cand.mother.t[which(cand.mother.t$Age.Class == 2), 'Offspring'] <- rbinom(n = length(which(cand.mother.t$Age.Class == 2)), size = 1, prob = fecundity[3])   #Randomly draws 0/1 calving fate for all young adults (2 years old), based on young adult fecundity value.
  cand.mother.t[which(cand.mother.t$Age.Class >= 3), 'Offspring'] <- rbinom(n = length(which(cand.mother.t$Age.Class >= 3)), size = 1, prob = fecundity[4])    #Randomly draws 0/1 calving fate for all adults (>= 3 years old), based on adult fecundity value.
  
  cand.mother.t<-subset(cand.mother.t,cand.mother.t$Offspring>0) #Drop subset potential mothers to those who had offspring
  
  #Assign fathers that spatially overlap with mothers during the previous fall breeding season
  mate.df <- c(as.character(cand.mother.t$Individual), as.character(potential.fathers$Individual)) #Identify all mates
  mate.points <- locations[which(locations$id %in% mate.df),] #Limit point data for mates only
  mates.mcp <- mcp(mate.points[,1], unin = "m", unout = "km2")  #Draw mcp for all mates
  
  if(year==1) {      #If this is the first year of the simulation...
    mates <- matrix(NA, ncol = length(potential.fathers$Individual), nrow = length(cand.mother.t$Individual), dimnames = list(as.character(cand.mother.t$Individual), as.character(potential.fathers$Individual))) #Create blank mate storage matrix (will store whom mates with whom)
    for(i in 1:nrow(cand.mother.t)) { #Begin loop to calculate overlap between males and females
      t1 <- mates.mcp[which(mates.mcp$id==cand.mother.t$Individual[i]),] #Select polygon for mother
      for(j in 1:length(potential.fathers$Individual)) {  #Cycle through fathers
        t2 <- mates.mcp[which(mates.mcp$id==potential.fathers$Individual[j]),] #Select candidate father polygon
        mates[i,j] <- gOverlap(t1,t2) #calculate overlap
      }
      if(sum(mates[i,])==0) {cand.mother.t$Mate[i] <- NA} else {cand.mother.t$Mate[i] <- as.numeric(sample(names(mates[i,]),1,replace=TRUE, prob = mates[i,]))}}}   #If no mates were identified, mother does not reproduce. Otherwise, a male is selected with weights based on the proportion of overlap. 
  
  if(year>1){ #If this is not the first year of the simulation
    mates.last <- mates   #Store last year's mate matrix to save computation time
    mates <- matrix(NA, ncol = length(potential.fathers$Individual), nrow = length(cand.mother.t$Individual), dimnames = list(as.character(cand.mother.t$Individual), as.character(potential.fathers$Individual)))
    for(i in 1:nrow(cand.mother.t)) {
      t1 <- mates.mcp[which(mates.mcp$id==cand.mother.t$Individual[i]),] #Select polygon for mother
      for(j in 1:length(potential.fathers$Individual)) {  #Cycle through fathers
        if(potential.fathers$Age.Class[j]>=4 & cand.mother.t$Age.Class[i]>=4 & rownames(mates)[i] %in% rownames(mates.last) & colnames(mates)[j] %in% colnames(mates.last)) {mates[i,j] <- mates.last[which(row.names(mates.last)==row.names(mates)[i]), which(colnames(mates.last)==colnames(mates)[j])]  #Use overlap already calculated for any adult (4+ YO) pair in the previous year (saves computation time)
        } else {  #Otherwise, take the MCP and calculate the overlap
          t2 <- mates.mcp[which(mates.mcp$id==potential.fathers$Individual[j]),]
          mates[i,j] <- gOverlap(t1,t2)
        }}
      if(sum(mates[i,])==0) {cand.mother.t$Mate[i] <- NA} else {cand.mother.t$Mate[i] <- as.numeric(sample(names(mates[i,]),1,replace=TRUE, prob = mates[i,]))  #If no mates were identified, mother does not reproduce. Otherwise, a male is selected with weights based on the proportion of overlap.
      }}
  }
  
  mother.t<-subset(cand.mother.t,!is.na(cand.mother.t$Mate)) #Drop anyone who does not have a mate
  
  #Create storage for offspring that match format of population dataframe
  offspring.t <- data.frame("Individual"=seq(last.ID+1, last.ID+sum(mother.t$Offspring), 1), "Male"=rep(NA, sum(mother.t$Offspring)), 
                            "Age Class"=rep(NA, sum(mother.t$Offspring)), "Mother"=rep(NA, sum(mother.t$Offspring)),
                            "Father"=rep(NA, sum(mother.t$Offspring)),
                            "Sample.P"=rep(NA,sum(mother.t$Offspring)), "Sample"=rep(0,sum(mother.t$Offspring)), "Spatial.Sample"=rep(0,sum(mother.t$Offspring)), "Sex" = rep("Female", sum(mother.t$Offspring)), "Age.Cat"=rep("Calf",sum(mother.t$Offspring)),"mu"=rep(NA, sum(mother.t$Offspring)))
  
  #Create offspring from mothers dataframe
  counter<-1 #Keeps track of row of offspring dataframe
  for(i in 1:nrow(mother.t)){    #loop through each mother
    for(j in 1:mother.t[i,'Offspring']){    #loop through each offspring
      offspring.t[counter,'Mother']<-mother.t[i,'Individual'] #Assign mother
      offspring.t[counter,'Father']<-mother.t[i,'Mate']       #Assign father
      counter<-counter+1                                      #moves on to next offspring, same mother or new mother
    }
  }
  
  offspring.t[,'Male']<-rbinom(nrow(offspring.t),1,sex.ratio) #Assign sex to offspring based on 1:1 sex ratio at birth
  offspring.t[,'Age.Class']<-0  #Assign calf age
  offspring.t[which(offspring.t$Male==1),'Sex'] <- "Male"
  
  #At the birth pulse, every animal ages before the next time step
  population[,'Age.Class']<-population[,'Age.Class']+1
  #population[,'Age.Cat'] <- ifelse(pop.0$Age.Class==1, "Yearling", ifelse(pop.0$Age.Class==2, "TwoYO", "Adult"))
  
  #Merge offspring and newly aged animals
  population<-rbind(population,offspring.t) #Merge Offspring and Adults
  return(list("population" = population, "neonates" = offspring.t, "mates" = mates))
}