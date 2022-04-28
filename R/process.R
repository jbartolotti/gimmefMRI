
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

writeGimmeCode <- function(savedir, mm, maxcores = .5){

  runmodel_fileConn <- file(file.path(savedir,'runmodels.R'),'wb')
  nummodels <- dim(mm$model_spec)[1]

  write('modelspecs <- list()',runmodel_fileConn, append = TRUE)
  write('cs_subgroups <- list()',runmodel_fileConn, append = TRUE)

  for(mi in 1:nummodels){
    thismod <- mm$model_spec[mi,]
    input_basedir = file.path(savedir,thismod$model_name,'input_files')
    output_basedir = file.path(savedir,thismod$model_name)

    dosubgroup <- !is.na(thismod$subgroup_thresh) && length(thismod$subgroup_thresh)>0
    if(dosubgroup){subcutoff <- thismod$subgroup_thresh}else{subcutoff <- 'NULL'}

    #write code to generate cs_subgroup dataframe from filelist to the .R file
    if(dosubgroup &&!is.na(thismod$subgroup_names) && length(thismod$subgroup_names)>0){

      grouplist <- mm$lists[[thismod$subgroup_names]]
      grouplist_string <- paste('grouplist <- c(',
                                paste(unlist(lapply(grouplist,function(x){return(sprintf("\'%s\'",x))})),collapse=','),
                                ')',collapse='')
      write(c(
        grouplist_string,
        "groupnum_lookup <- list()",
        "index <- 0",
        "for(g in grouplist){index <- index+1; groupnum_lookup[g] = index}",
        sprintf("cs_subgroups[[%s]] <- data.frame(filename = dir(\'%s\'), groupnum = 0, stringsAsFactors = FALSE)",mi,input_basedir),
        sprintf("cs_subgroups[[%s]]$filename <- unlist(lapply(cs_subgroups[[%s]]$filename, function(x){gsub(\'.csv\',\'\',x)}))",mi,mi),
        sprintf("cs_subgroups[[%s]]$groupnum <- unlist(lapply(cs_subgroups[[%s]]$filename, function(x){groupnum_lookup[[strsplit(x,\'_\')[[1]][1]]]}))",mi,mi)),
        runmodel_fileConn, append = TRUE
      )
      cssubgroup_write <- "\'cs_subgroups[[i]]\'"
    } else{cssubgroup_write <- 'NULL'}

    write(c(
      sprintf("modelspecs[[%s]] <- list(",mi),
      sprintf("      data = \'%s\'",input_basedir),
      sprintf("      out = \'%s\',",output_basedir),
      "      sep = \',\',",
      "      ar = TRUE,",
      "      plot = TRUE,",
      sprintf("      groupcutoff = %s,",thismod$group_thresh),
      sprintf("      subgroup = %s,",dosubgroup),
      sprintf("      subcutoff = %s,",subcutoff),
      sprintf("      confirm_subgroup = %s",cssubgroup_write),
      "      paths = NULL",
      "      exogenous = NULL",
      ")"),
      runmodel_fileConn,append = TRUE)
  }

  numcores <- useCores(maxcores)
  if (numcores == 1){
    fortext <- 'i <- 1'
    #  fortext <- sprintf('for(i in 1:%s){',nummodels)
  }else{
    fortext <- c(sprintf('registerDoParallel(%s)',numcores),sprintf('foreach(i=1:%s) %%dopar%% {',nummodels))
  }
  write(fortext, runmodel_fileConn,append = TRUE)


    write(c(
      "  gimme(data = modelspecs[[i]]$data,",
      "      out = modelspecs[[i]]$out,",
      "      sep = modelspecs[[i]]$sep,",
      "      ar = modelspecs[[i]]$ar,",
      "      plot = modelspecs[[i]]$plot,",
      "      groupcutoff = modelspecs[[i]]$groupcutoff,",
      "      subgroup = modelspecs[[i]]$subgroup,",
      "      subcutoff = modelspecs[[i]]$subcutoff,",
      "      confirm_subgroup = modelspecs[[i]]$confirm_subgroup,",
      "      paths = modelspecs[[i]]$paths,",
      "      exogenous = modelspecs[[i]]$exogenous",
      "  )"),
      runmodel_fileConn, append = TRUE)



  #close parallel loop if necessary
  if (numcores > 1){
    write(c(
      '}'),
      runmodel_fileConn, append = TRUE)
  }
  close(runmodel_fileConn)
}


oldfun <- function(){
  for(mi in 1:nummodels){
    thismod <- mm$model_spec[mi,]
    input_basedir = file.path(savedir,thismod$model_name,'input_files')
    output_basedir = file.path(savedir,thismod$model_name)

    #write code to generate cs_subgroup dataframe from filelist to the .R file
    if(!is.na(thismod$subgroup_thresh) && length(thismod$subgroup_thresh)>0 &&
       !is.na(thismod$subgroup_names) && length(thismod$subgroup_names)>0){

      grouplist <- mm$lists[[thismod$subgroup_names]]
      grouplist_string <- paste('grouplist <- c(',
            paste(unlist(lapply(grouplist,function(x){return(sprintf("\'%s\'",x))})),collapse=','),
            ')',collapse='')
      write(c(
        grouplist_string,
        "groupnum_lookup <- list()",
        "index <- 0",
        "for(g in grouplist){index <- index+1; groupnum_lookup[g] = index}",
        sprintf("cs_subgroups = data.frame(filename = dir(\'%s\'), groupnum = 0, stringsAsFactors = FALSE)",input_basedir),
        "cs_subgroups$filename = unlist(lapply(cs_subgroups$filename, function(x){gsub(\'.csv\',\'\',x)}))",
        "cs_subgroups$groupnum = unlist(lapply(cs_subgroups$filename, function(x){groupnum_lookup[[strsplit(x,\'_\')[[1]][1]]]}))"),
        runmodel_fileConn, append = TRUE
        )
    }

    write(c(
        sprintf("  gimme(data = \'%s\',",input_basedir),
        sprintf("      out = \'%s\',",output_basedir),
                "      sep = \',\',",
                "      ar = TRUE,",
                "      plot = TRUE,",
        sprintf("      groupcutoff = %s,",thismod$group_thresh)),
      runmodel_fileConn, append = TRUE)

    #subgroups
    if(!is.na(thismod$subgroup_thresh) && length(thismod$subgroup_thresh)>0){

      write(c(
                "      subgroup = TRUE,",
        sprintf("      subcutoff = %s,",thismod$subgroup_thresh)),
        runmodel_fileConn, append = TRUE)
    } else{
      write(c(
        "      subgroup = FALSE,"),
        runmodel_fileConn, append = TRUE)
    }
    #apriori subgroups
    if(!is.na(thismod$subgroup_thresh) && length(thismod$subgroup_thresh)>0 &&
       !is.na(thismod$subgroup_names) && length(thismod$subgroup_names)>0){
      write(c(
        sprintf("      confirm_subgroup = %s,",'SUBGROUPS')),
        runmodel_fileConn, append = TRUE)
    } else {
      write(c(
                "      confirm_subgroup = NULL,"),
        runmodel_fileConn, append = TRUE)
    }
    #apriori connections
    if(FALSE){
    } else{
    write(c(
      "      paths = NULL,"),
      runmodel_fileConn, append = TRUE)
    }
    #exogenous factors
    if(FALSE){
    } else{
    write(c(
     "      exogenous = NULL)"),
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

