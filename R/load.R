timecoursesToRdat <- function(timecourse_locations){
  #timecourse_locations can be
    # - a csv that contains absolute paths for each timecourse
    # - a vector that contains absolute paths
    # - a list with any number of basedirs that contain fields for basedir, and then the filenames within there
  #

}

readModelInput <- function(file_timecourses, file_modelspec, file_lists, file_shortnames){
  mm <- list()
  tc <- read.csv(file_timecourses)
  mspec <- read.csv(file_modelspec)
  mylists <- read.csv(file_lists)
  sn <- read.csv(file_shortnames)

  mm$timecourses <- tc
  mm$model_spec <- mspec
  mm$lists <- list()
  for(listname in colnames(mylists)){
    mm$lists[[listname]] <- mylists[!(is.na(mylists[,listname]) | mylists[,listname]==''), listname]
    }
  mm$shortnames <- sn
  return(mm)
  }
