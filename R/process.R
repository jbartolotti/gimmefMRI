
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

writeGimmeCode <- function(savedir, mm, maxcores = .5){
  numcores <- useCores(maxcores)
  nummodels <- dim(mm$model_spec)[1]
  if (numcores == 1){
    fortext <- as.character()
    #  fortext <- sprintf('for(i in 1:%s){',nummodels)
  }else{
    fortext <- c(sprintf('registerDoParallel(%s)',numcores),sprintf('foreach(i=1:%s) %%dopar%% {',nummodels))
  }

  runmodel_fileConn <- file(file.path(savedir,'runmodels.R'),'wb')
  write(fortext, runmodel_fileConn,append = TRUE)

  for(mi in 1:nummodels){
    thismod <- mm$model_spec[mi,]
    input_basedir = file.path(savedir,thismod$model_name,'input_files')
    output_basedir = file.path(savedir,thismod$model_name)
    write(c(
        sprintf('  gimme(data = %s,',input_basedir),
        sprintf('      out = %s,',output_basedir),
                '      sep = \',\'',
                '      ar = TRUE,',
                '      plot = TRUE,',
        sprintf('      groupcutoff = %s,',thismod$group_thresh)),
      runmodel_fileConn, append = TRUE)

    #subgroups
    if(!is.na(thismod$subgroup_thresh) && length(thismod$subgroup_thresh)>0){
      write(c(
                '      subgroup = TRUE,',
        sprintf('      subcutoff = %s,',thismod$subgroup_thresh)),
        runmodel_fileConn, append = TRUE)
    } else{
      write(c(
        '      subgroup = FALSE,'),
        runmodel_fileConn, append = TRUE)
    }
    #apriori subgroups
    if(!is.na(thismod$subgroup_thresh) && length(thismod$subgroup_thresh)>0 &&
       !is.na(thismod$subgroup_names) && length(thismod$subgroup_names)>0){
      write(c(
        sprintf('      confirm_subgroup = %s,','SUBGROUPS')),
        runmodel_fileConn, append = TRUE)
    } else {
      write(c(
                '      confirm_subgroup = NULL,'),
        runmodel_fileConn, append = TRUE)
    }
    #apriori connections
    if(FALSE){
    } else{
    write(c(
      '      paths = NULL,'),
      runmodel_fileConn, append = TRUE)
    }
    #exogenous factors
    if(FALSE){
    } else{
    write(c(
     '      exogenous = NULL)'),
     runmodel_fileConn, append = TRUE)
    }
   #write this gimme code to file. Will either run sequentially or parallelized

  }
  #close parallel loop if necessary
  if (numcores > 1){
    write(c(
      '}'),
      runmodel_fileConn, append = TRUE)
  }

}

