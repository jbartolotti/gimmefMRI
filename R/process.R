
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

initializeGimmeFolders <- function(savedir,mm){
  for (i in 1:dim(mm$model_spec)[1] ){
    thismod <- mm$model_spec[i,]
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


getTimecourse <- function(sublist, subdf, roi_list){

  #write big shell file that will extract timecourses from all subjects from all rois to individual tc files
  extract_timecourse_filename <- file.path('...')
  xtc_fileConn <- file(extract_timecourse_filename,'wb')

  filelocs <- list()
  for(s in sublist){
      filelocs[[s]] <- getTimecourse_OneSub(xtc_fileConn, s, roi_list,
                                            #subdir, roidir, timecoursedir,
                                            ...)
  }

  #run the shell file

  #read in all those tc files and save them to the big timecourses.csv to use as input to gimmefMRI()

}

getTimecourse_OneSub <- function(fileConn, subname, roi_list,
                                 subdir, roidir, timecoursedir,
                                 sub_functional_filename, sub_mask_filename,
                                 submaskdir = subdir,
                                 subroi_dirname = 'subject_rois',
                                 remove_subrois = FALSE
                                 ){
  filelocs = list()
  for(roiname in roi_list$name){
    roi_filename <- roi_list$filename[roi_list$name == roiname]
    filelocs[[roi]] <- getTimecourse_OneSub_OneROI(fileConn, subname, roiname,
                                                   subdir, roidir, timecoursedir,
                                                   sub_functional_filename, sub_mask_filename, roi_filename,
                                                   submaskdir = submaskdir, subroi_dirname = subroi_dirname, remove_subrois = remove_subrois)
  }
  return(filelocs)
}

getTimecourse_OneSub_OneROI <- function(fileConn, subname, roiname,
                                        subdir, roidir, timecoursedir,
                                        sub_functional_filename, sub_mask_filename, roi_filename,
                                        submaskdir = subdir, subroi_dirname = 'subject_rois', remove_subrois = FALSE){
  #resample the roi to the subjects functional
  subfunc_loc <-  file.path(subdir,sub_functional_filename)
  roisub_loc <- file.path(subdir,subroi_dirname,sprintf('SUB_%s_ROI_%s',subname,roiname))
  dir.create(file.path(subdir, subroi_dirname), showWarnings = FALSE)
  roi_loc <- file.path(roidir, roi_filename)
  line_resample <- sprintf('3dresample -master %s -prefix %s -inset %s -overwrite',subfunc_loc, roisub_loc, roi_loc)

  #combine subroi with the subject mask
  sub_mask_loc <- file.path(submaskdir, sub_mask_filename)
  line_3dcalc <- sprintf("3dcalc -a %s -b %s -expr \'a*b\' -nscale -overwrite -prefix %s", roisub_loc, sub_mask_loc, roisub_loc)

  #save timecourse in this mask to text
  timecourse_loc <- file.path(timecoursedir,  sprintf('SUB_%s_ROI_%s_timecourse.txt',subname,roiname))
  line_3dmaskave <- sprintf("3dmaskave -mask %s -quiet %s > %s",roisub_loc, subfunc_loc, timecourse_loc)

  write(c(line_3dresample, line_3dcalc, line_3dmaskave), fileConn, append = TRUE)
  return(list(subfunc_loc = subfunc_loc, roisub_loc = roisub_loc, roi_loc = roi_loc))
}
