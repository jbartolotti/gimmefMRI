
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
