
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
      filelocs[[s]] <- getTimecourse_OneSub(...)
  }

  #run the shell file

  #read in all those tc files and save them to the big timecourses.csv to use as input to gimmefMRI()

}

getTimecourse_OneSub <- function(... ,
                                 sub_functional, sub_mask, roi_list){
  filelocs = list()
  for(roi in roi_list){
    filelocs[[roi]] <- getTimecourse_OneSub_OneROI(...)
  }
  return(filelocs)
}

getTimecourse_OneSub_OneROI <- function(fileConn, subname, roiname,
                                        subdir, roidir,
                                        sub_functional_filename, sub_mask_filename, roi_filename,
                                        submaskdir = subdir, subroi_dirname = 'subject_rois', remove_subrois = FALSE){
  #resample the roi to the subjects functional
  subfunc_loc <-  file.path(subdir,sub_functional_filename)
  roisub_loc <- file.path(subdir,subroi_dirname,sprintf('SUB_%s_ROI_%s'))
  roi_loc <- file.path(roidir, roi_filename)
  line_resample <- sprintf('3dresample -master %s -prefix %s -inset %s -overwrite',subfunc_loc, roisub_loc, roi_loc)
  #combine subroi with the subject mask
     #3dcalc -a ROIA_SUB9999 -b SUB9999_MASK -exrp 'a*b' -nscale -overwrite -prefix ROIA_SUB9999
  #save timecourse in this mask to text
     #3dmaskave -mask ROIA_SUB9999 -quiet SUB9999_FUNC > ROIA_SUB999_TIMECOURSE.txt
  write(c(line_3dresample, line_3dcalc, line_3dmaskave), fileConn, append = TRUE)
  return(list(subfunc_loc = subfunc_loc, roisub_loc = roisub_loc, roi_loc = roi_loc))
}
