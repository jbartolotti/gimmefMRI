
applyShorten <- function(vector,lookup){
  #input a vector and lookup table.
  #For anything in the vector thats in the lookup table, replace
  #with its short form.
  sn <- unlist(lapply(vector,
                      function(x){
                        val<-lookup$short[lookup$long==x]
                        if(length(val)==0){val<-x}
                        return(val)
                      }))
  return(sn)
}

filenameDateTime <- function()
{
  datetime <- Sys.time()
  datetime <- gsub('-','',datetime)
  datetime <- gsub(':','',datetime)
  datetime <- gsub(' ','_',datetime)
  return(datetime)
}

useCores <- function(maxcores){
  #if maxcores is >= 1, use that many, but not more than actual available.
  # if maxcores < 1, treat as percent of number of cores, rounded down, at least one.
  # if maxcores is 0 (or neg), use one core
  if(maxcores > 0){
    numcores <- parallel::detectCores()
    if(maxcores >= 1){usecores <- min(maxcores, numcores)} else
    {usecores <- max(1, floor(numcores*maxcores))}
  } else{usecores <- 1}
  return(usecores)
}
