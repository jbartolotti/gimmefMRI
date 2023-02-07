timecoursesToRdat <- function(timecourse_locations){
  #timecourse_locations can be
    # - a csv that contains absolute paths for each timecourse
    # - a vector that contains absolute paths
    # - a list with any number of basedirs that contain fields for basedir, and then the filenames within there
  #

}

readXLSXinput <- function(xlsx_file,savedir){
  cntrl_wide <- as.data.frame(readxl::read_excel(xlsx_file, sheet = 'control', col_names = c('help','name','value')))
  tc <- as.data.frame(readxl::read_excel(xlsx_file, sheet = 'timecourses',col_names = TRUE))
  mspec_wide <- as.data.frame(readxl::read_excel(xlsx_file, sheet = 'models', col_names = TRUE))
  mylists <- as.data.frame(readxl::read_excel(xlsx_file, sheet = 'lists',col_names = TRUE))
  sn <- as.data.frame(readxl::read_excel(xlsx_file, sheet = 'abbreviations',col_names = TRUE))
  figs <- as.data.frame(readxl::read_excel(xlsx_file, sheet = 'figures',col_names = TRUE))

  cntrl <- list()
  for(i in 1:dim(cntrl_wide)[1]){
    cntrl[[cntrl_wide$name[i]]] <- cntrl_wide$value[i]
  }
  cntrl$generate_model_code <- as.logical(cntrl$generate_model_code)
  cntrl$run_model_code <- as.logical(cntrl$run_model_code)
  cntrl$generate_figure_code <- as.logical(cntrl$generate_figure_code)
  cntrl$run_figure_code <- as.logical(cntrl$run_figure_code)

  cntrl$script_save_directory <- sub('{PWD}',savedir,cntrl$script_save_directory, fixed = TRUE)
  cntrl$model_save_directory <- sub('{PWD}',savedir,cntrl$model_save_directory, fixed = TRUE)
  cntrl$figure_save_directory <- sub('{PWD}',savedir,cntrl$figure_save_directory, fixed = TRUE)
  #  sub('PWD',savedir,'PWD/models')

  mspec <- t(mspec_wide[,3:(dim(mspec_wide)[2])])
  colnames(mspec) <- mspec_wide[,2]
  rownames(mspec) <- 1:(dim(mspec)[1])
  mspec <- as.data.frame(mspec)
  mspec$group_thresh <- as.numeric(mspec$group_thresh)
  mspec$subgroup_thresh <- as.numeric(mspec$subgroup_thresh)
  mspec$subgroups <- as.logical(mspec$subgroups)
  mspec$apriori_subgroups <- as.logical(mspec$apriori_subgroups)

  figs2 <- t(figs[,3:(dim(figs)[2])])
  colnames(figs2) <- figs[,2]
  rownames(figs2) <- 1:(dim(figs2)[1])
  figs2 <- as.data.frame(figs2)

  return(readModelInput(cntrl,tc,mspec,mylists,sn,figs2))
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

readModelInput <- function(cntrl,tc,mspec,mylists,sn,figs){
  mm <- list()
  mm$cntrl <- cntrl
  mm$timecourses <- tc
  mm$model_spec <- mspec
  mm$lists <- list()
  for(listname in colnames(mylists)){
    mm$lists[[listname]] <- mylists[!(is.na(mylists[,listname]) | mylists[,listname]==''), listname]
    }
  mm$shortnames <- sn
  mm$figures <- figs
  return(mm)
}


applyCustomSettings <- function(mm, run, models){

  #Update Model List in config file
  if (length(models) == 1 && (models == 'all' || models == 'use_config')){
    model_list <- mm$model_spec$model_name
  } else{
    model_list <- mm$model_spec$model_name
    if(any(!(models %in% model_list))){
      warning(sprintf('WARNING: The following models are not present in the config file and will be skipped:\n%s',paste(models[!(models %in% model_list)], collapse = '\n')))
      model_list <- model_list[model_list %in% models]
    }
  }
  mm$model_spec <- mm$model_spec[mm$model_spec$model_name %in% model_list]

  #Update Runlist in config file
  run_options <- c('generate_model_code','run_model_code','generate_figure_code','run_figure_code')
  if( !( length(run) == 1 && (run == 'use_config')) ){
    if(length(run) == 1 && run == 'all'){
      mm$cntrl$generate_model_code = TRUE
      mm$cntrl$run_model_code = TRUE
      mm$cntrl$generate_figure_code = TRUE
      mm$cntrl$run_figure_code = TRUE
    } else{
      if(any(!(run %in% run_options))){
        warning(sprintf("WARNING: Allowable run options are '%s'. The following run options were ignored:\n%s",
                        paste(run_options, collapse = "', '"),
                        paste(run[!(run %in% run_options)], collapse = '\n')))
      }

      mm$cntrl$generate_model_code = 'generate_model_code' %in% run
      mm$cntrl$run_model_code = 'run_model_code' %in% run
      mm$cntrl$generate_figure_code = 'generate_figure_code' %in% run
      mm$cntrl$run_figure_code = 'run_figure_code' %in% run
    }
  }


  return(mm)
}

