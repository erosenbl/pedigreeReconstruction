#Calculates degree of overlap between potential mates.

#The inputs are the minimum convex polygons for an individual of interest (hr1) and a potential mate (hr2).

#The output is proportion of the individul's home range shared by the potential mate.

gOverlap <- function(hr1, hr2, ...) { #Calculates degree of overlap between potential mates
  a <- gIntersection(hr1, hr2, ...)
  if (is.null(a)) {
    return(0)
  }
  gArea(a, byid=TRUE) / gArea(hr1, byid=TRUE) 
}
