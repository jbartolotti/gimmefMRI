
#input:
#model_spec csv or rds
  #each row a model that will be run. Columns are:
    #group threshold
    #subgroup threshold
    #network name (name found in network_list)
    #conditions to censor

#network_list csv or rds
  #each column header is the network name.
  #rows are the predictors to include. Mark exogenous predictors somehow?

#alldata csv or rds
  #columns: each ROI, slicenum, time, condition, censor, subnum, subgroup, run
  #rows: TR*sub*run

#roinames
  #columns: longroi, shortroi

initializeGimmeFolders <- function(mm){

  savedir <- mm$cntrl$model_save_directory
  dir.create(savedir, showWarnings = FALSE)
  for (modname in mm$model_spec$model_name ){
    thismod <- mm$model_spec[mm$model_spec$model_name == modname,]
    #make target directory
    dir.create(file.path(savedir,thismod$model_name),showWarnings = FALSE)
    unlink(file.path(savedir,thismod$model_name,'input_files'),recursive = TRUE)
    dir.create(file.path(savedir,thismod$model_name,'input_files'),showWarnings = FALSE)

    #fill target directory with subject .csv files

    #subset timecourse data for this run, and rois in this network
    roicols <- mm$lists[[thismod$network_name]]
    datacols <- c('slicenum', 'time', 'condition', 'censor', 'subnum', 'subgroup', 'run')
    colkeep <- c(roicols, datacols)
    rowkeep <- mm$timecourses$run == thismod$run
    dd <- mm$timecourses[rowkeep,colkeep]

    #for each subject, write timecourses to file, after censoring
    for(s in mm$lists[[thismod$subselect]])
    {
      thissub <- subset(dd, subnum == s)
      thisgroup <- unique(thissub$subgroup)

      #remove censored TRs and conditions
      thissub[thissub$censor == 0, roicols] = NA
      if (!is.na(thismod$censor_conditions)){
        for(censor_cond in mm$lists[[thismod$censor_conditions]]){
          thissub[thissub$condition == censor_cond, roicols] = NA
        }
      }
      short_roicols <- applyShorten(roicols,mm$shortnames)
      colnames(thissub) <- applyShorten(colnames(thissub),mm$shortnames)

      if (dim(thissub)[1] == 0 ){warning(sprintf('Empty dataset for %s. File skipped.',s))}
      else{
        write.csv(thissub[,short_roicols],
                file.path(savedir,thismod$model_name,'input_files',
                          sprintf('%s_%s.csv',thisgroup,s)), row.names = FALSE)
      }
    }

  }

}

doreplace <- function(line, pattern_list, target){
  for(pattern in pattern_list){
    if(!is.na(line[pattern])){
      message(line)
      line[target] <- gsub(paste0('{',pattern,'}'),line[pattern],line[target], fixed = TRUE)
    }
  }
  return(line)
}

#' @export
getTimecourse <- function(write_file = 'extract_timecourses.sh', config_file = 'gui', sub_list = 'all', roi_list = 'all'){

  if(config_file=='gui'){
    if(Sys.info()[['sysname']] == 'Linux'){
      myfile <- my.file.browse()
    } else {
      myfile <- file.choose()
    }
  } else{
    myfile <- config_file
  }
  # Read config file and replace search patterns e.g. {BASE_DIR} with appropriate columns
  config <- readr::read_csv(myfile, col_types = readr::cols(.default = "c"))
  for(linenum in 1:dim(config)[1]){
    line <- config[linenum,]
    line <- doreplace(line,c('BASE_DIR','GROUP','ID','RUN'),'DATA_DIR')
    line <- doreplace(line,c('DATA_DIR','GROUP','ID','RUN'),'ROI_FILENAME')
    line <- doreplace(line,c('DATA_DIR','GROUP','ID','RUN'),'FUNCTIONAL_FILENAME')
    line <- doreplace(line,c('DATA_DIR','GROUP','ID','RUN'),'MASK_FILENAME')
    line <- doreplace(line,c('DATA_DIR','GROUP','ID','RUN'),'CENSOR_FILENAME')
    line <- doreplace(line,c('DATA_DIR','GROUP','ID','RUN'),'SUB_ROI_DIR')
    line <- doreplace(line,c('DATA_DIR','GROUP','ID','RUN'),'TIMECOURSE_DIR')

    config[linenum,] <- line
  }

  # process sub_list and roi_list to ensure requested entries are present in the config file
  config_sub_list <- config$ID[config$TYPE == 'sub']
  config_roi_list <- config$ID[config$TYPE == 'roi']

  if(length(sub_list) == 1 && sub_list == 'all'){sub_list <- config_sub_list}
  if(length(roi_list) == 1 && roi_list == 'all'){roi_list <- config_roi_list}

  if(any(! sub_list %in% config_sub_list)){
    message(sprintf('WARNING: These requested subjects were not found in the config file: %s', paste(sub_list[! sub_list %in% config_sub_list],collapse=', ')))
  }
  if(any(! roi_list %in% config_roi_list)){
    message(sprintf('WARNING: These requested ROIs were not found in the config file: %s', paste(roi_list[! roi_list %in% config_roi_list],collapse=', ')))
  }


  #write big shell file that will extract timecourses from all subjects from all rois to individual tc files
  if(is.na(write_file)){
    xtc_fileConn <- NA}
  else {
      xtc_fileConn <- file(write_file,'wb')
  }

  filelocs <- list()
  for(s in sub_list){
    line <- config[config$ID == s,]
      filelocs[[s]] <- getTimecourse_OneSub(xtc_fileConn, s, config[config$ID %in% roi_list,],
                                            line$DATA_DIR,
                                            line$TIMECOURSE_DIR,
                                            line$FUNCTIONAL_FILENAME,
                                            line$MASK_FILENAME,
                                            submaskdir = line$DATA_DIR,
                                            subroi_dir = line$SUB_ROI_DIR,
                                            remove_subrois = FALSE
                                            )
  }
  if(!is.na(xtc_fileConn)){close(xtc_fileConn)}
  return(list(filelocs = filelocs, config = config))
}

genTimecoursesCSV_config <- function(tcfilename, config = 'gui'){
  if(config_file=='gui'){
    myfile <- file.choose()
  } else{
    myfile <- config_file
  }
  # Read config file and replace search patterns e.g. {BASE_DIR} with appropriate columns
  config <- read.csv(myfile)
  for(linenum in 1:dim(config)[1]){
    line <- config[linenum,]
    line <- doreplace(line,c('BASE_DIR','GROUP','ID','RUN'),'DATA_DIR')
    line <- doreplace(line,c('DATA_DIR','GROUP','ID','RUN'),'ROI_FILENAME')
    line <- doreplace(line,c('DATA_DIR','GROUP','ID','RUN'),'FUNCTIONAL_FILENAME')
    line <- doreplace(line,c('DATA_DIR','GROUP','ID','RUN'),'MASK_FILENAME')
    line <- doreplace(line,c('DATA_DIR','GROUP','ID','RUN'),'CENSOR_FILENAME')
    line <- doreplace(line,c('DATA_DIR','GROUP','ID','RUN'),'SUB_ROI_DIR')
    line <- doreplace(line,c('DATA_DIR','GROUP','ID','RUN'),'TIMECOURSE_DIR')

    config[linenum,] <- line
  }


}

genTimecoursesCSV <- function(tcfilename, filelocs, config = NA){
    allrois <- as.character()
    for(s in names(filelocs)){
      allrois <- c(allrois,names(filelocs[[s]]))
    }
    allrois <- unique(allrois)

    eachdf <- list()

    for(s in names(filelocs)){
      if( !(typeof(config)=='logical' && is.na(config)) ){
        group <- unlist(config[config$ID == s,'GROUP'])
        thisrun <- unlist(config[config$ID == s,'RUN'])
        thiscensor <- unlist(config[config$ID == s,'CENSOR_FILENAME'])
      } else{
        group <- NA
        thisrun <- NA
        thiscensor <- NA
      }
      for(r in names(filelocs[[s]])){
        thistc_filename <- filelocs[[s]][[r]]$timecourse_loc
        thistc <- read.csv(thistc_filename,header = FALSE)

        # initialize this subject's timecourse df if it doesn't exist
        if(! s %in% names(eachdf)){
          tclength <- dim(thistc)[1]
          eachdf[[s]]  <- data.frame(slicenum = 1:tclength, time = rep(NA,1,tclength), condition = rep(NA,1,tclength), censor = rep(0,1,tclength), subnum = rep(s,1,tclength), subgroup = rep(group,1,tclength), run = rep(thisrun,1,tclength), stringsAsFactors = FALSE)
          eachdf[[s]][,allrois] <- NA
          if(!is.na(thiscensor)){
            cens <- read.csv(thiscensor, header = FALSE)
            eachdf[[s]]$censor <- cens[1:tclength,1]
          }
        }
        eachdf[[s]][,r] <- thistc
      }
    }
    alldf <- do.call('rbind',eachdf)
    write.csv(alldf, file = tcfilename,row.names = FALSE)
  }

getTimecourse_OneSub <- function(fileConn, subname, roi_df,
                                 subdir,
                                 timecoursedir,
                                 sub_functional_filename,
                                 sub_mask_filename,
                                 submaskdir = subdir,
                                 subroi_dir = file.path(subdir,'subject_rois'),
                                 remove_subrois = FALSE
                                 ){
  filelocs = list()
  for(roiname in roi_df$ID){
    line <- roi_df[roi_df$ID == roiname,]
    filelocs[[roiname]] <- getTimecourse_OneSub_OneROI(fileConn, subname, roiname,
                                                   subdir,
                                                   line$DATA_DIR,
                                                   timecoursedir,
                                                   sub_functional_filename,
                                                   sub_mask_filename,
                                                   line$ROI_FILENAME,
                                                   submaskdir = submaskdir, subroi_dir = subroi_dir, remove_subrois = remove_subrois)
  }
  return(filelocs)
}

getTimecourse_OneSub_OneROI <- function(fileConn, subname, roiname,
                                        subdir,
                                        roidir,
                                        timecoursedir,
                                        subfunc_loc, sub_mask_loc, roi_loc,
                                        submaskdir = subdir, subroi_dir = file.path(subdir,'subject_rois'), remove_subrois = FALSE){
  #resample the roi to the subjects functional
#  subfunc_loc <-  file.path(subdir,sub_functional_filename)
  roisub_loc <- file.path(subroi_dir,sprintf('SUB_%s_ROI_%s+tlrc',subname,roiname))
  dir.create(file.path(subroi_dir), showWarnings = FALSE)
  dir.create(file.path(timecoursedir), showWarnings = FALSE)
#  roi_loc <- file.path(roidir, roi_filename)
  line_3dresample <- sprintf('3dresample -master %s -prefix %s -inset %s -overwrite',subfunc_loc, roisub_loc, roi_loc)

  #combine subroi with the subject mask
#  sub_mask_loc <- file.path(submaskdir, sub_mask_filename)
  line_3dcalc <- sprintf("3dcalc -a %s -b %s -expr \'a*b\' -nscale -overwrite -prefix %s", roisub_loc, sub_mask_loc, roisub_loc)

  #save timecourse in this mask to text
  timecourse_loc <- file.path(timecoursedir,  sprintf('SUB_%s_ROI_%s_timecourse.txt',subname,roiname))
  line_3dmaskave <- sprintf("3dmaskave -mask %s -quiet %s > %s",roisub_loc, subfunc_loc, timecourse_loc)


  if(!is.na(fileConn)){write(c(line_3dresample, line_3dcalc, line_3dmaskave), fileConn, append = TRUE)}
  return(list(subfunc_loc = subfunc_loc, roisub_loc = roisub_loc, roi_loc = roi_loc, timecourse_loc = timecourse_loc))
}
