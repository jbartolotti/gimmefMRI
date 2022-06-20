timecoursesToRdat <- function(timecourse_locations){
  #timecourse_locations can be
    # - a csv that contains absolute paths for each timecourse
    # - a vector that contains absolute paths
    # - a list with any number of basedirs that contain fields for basedir, and then the filenames within there
  #

}

readXLSXinput <- function(xlsx_file){
  cntrl <- as.data.frame(readxl::read_excel(xlsx_file, sheet = 'control', col_names = c('help','name','value')))
  tc <- as.data.frame(readxl::read_excel(xlsx_file, sheet = 'timecourses',col_names = TRUE))
  mspec_wide <- as.data.frame(readxl::read_excel(xlsx_file, sheet = 'models', col_names = FALSE))
  mylists <- as.data.frame(readxl::read_excel(xlsx_file, sheet = 'lists',col_names = TRUE))
  sn <- as.data.frame(readxl::read_excel(xlsx_file, sheet = 'abbreviations',col_names = TRUE))
  figs <- as.data.frame(readxl::read_excel(xlsx_file, sheet = 'figures',col_names = TRUE))


  mspec <- t(mspec_wide[,3:(dim(mspec_wide)[2])])
  colnames(mspec) <- mspec_wide[,2]
  rownames(mspec) <- 1:(dim(mspec)[1])
  mspec <- as.data.frame(mspec)
  mspec$group_thresh <- as.numeric(mspec$group_thresh)
  mspec$subgroup_thresh <- as.numeric(mspec$subgroup_thresh)
  mspec$subgroups <- as.logical(mspec$subgroups)
  mspec$apriori_subgroups <- as.logical(mspec$apriori_subgroups)

  return(readModelInput(cntrl,tc,mspec,mylists,sn))
}

readCSVInput <- function(datadir, file_timecourses = NA, file_modelspec = NA, file_lists = NA, file_shortnames = NA){
  if (is.na(file_timecourses)) {file_timecourses = file.path(datadir,'timecourses.csv')}
  if (is.na(file_modelspec)) {file_modelspec = file.path(datadir,'model_spec.csv')}
  if (is.na(file_lists)) {file_lists = file.path(datadir,'info_lists.csv')}
  if (is.na(file_shortnames)) {file_shortnames = file.path(datadir,'shortnames.csv')}
  tc <- read.csv(file_timecourses)
  mspec <- read.csv(file_modelspec)
  mylists <- read.csv(file_lists)
  sn <- read.csv(file_shortnames)

  return(readModelInput(NA,tc,mspec,mylists,sn))
}

readModelInput <- function(cntrl,tc,mspec,mylists,sn){
  mm <- list()
  mm$cntrl <- cntrl
  mm$timecourses <- tc
  mm$model_spec <- mspec
  mm$lists <- list()
  for(listname in colnames(mylists)){
    mm$lists[[listname]] <- mylists[!(is.na(mylists[,listname]) | mylists[,listname]==''), listname]
    }
  mm$shortnames <- sn
  return(mm)
}
