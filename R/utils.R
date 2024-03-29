
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

UTILS.getFilename <- function(input_filename)
{
  if(input_filename=='gui')
  {
    if(Sys.info()[['sysname']] == 'Linux')
    {
      myfile <- my.file.browse()
    } else {
      myfile <- file.choose()
    }
  } else{
    myfile <- input_filename
  }
  return(myfile)
}


#https://stackoverflow.com/questions/9122600/r-command-line-file-dialog-similar-to-file-choose

#' Text-based interactive file selection.
#'@param root the root directory to explore
#'             (default current working directory)
#'@param multiple boolean specifying whether to allow
#'                 multiple files to be selected
#'@return character vector of selected files.
#'@examples
#'fileList <- my.file.browse()
my.file.browse <- function (root=getwd(), multiple=F) {
  # .. and list.files(root)
  x <- c( dirname(normalizePath(root)), list.files(root,full.names=T) )
  isdir <- file.info(x)$isdir
  obj <- sort(isdir,index.return=T,decreasing=T)
  isdir <- obj$x
  x <- x[obj$ix]
  lbls <- sprintf('%s%s',basename(x),ifelse(isdir,'/',''))
  lbls[1] <- sprintf('../ (%s)', basename(x[1]))

  files <- c()
  sel = -1
  while ( TRUE ) {
    sel <- menu(lbls,title=sprintf('Select file(s) (0 to quit) in folder %s:',root))
    if (sel == 0 )
      break
    if (isdir[sel]) {
      # directory, browse further
      files <- c(files, my.file.browse( x[sel], multiple ))
      break
    } else {
      # file, add to list
      files <- c(files,x[sel])
      if ( !multiple )
        break
      # remove selected file from choices
      lbls <- lbls[-sel]
      x <- x[-sel]
      isdir <- isdir[-sel]
    }
  }
  return(files)
}
