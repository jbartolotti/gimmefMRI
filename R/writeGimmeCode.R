
writeGimmeCode <- function(mm, maxcores = 1){

  #Open file and write header
  dir.create(mm$cntrl$script_save_directory, showWarnings = FALSE)
  model_filename <- file.path(mm$cntrl$script_save_directory,'run_models.R')
  runmodel_fileConn <- file(model_filename,'wb')
  WGC_header(runmodel_fileConn, maxcores)

  #Write model specifications for each model
  nummodels <- dim(mm$model_spec)[1]
  for(mi in 1:nummodels){
    thismod <- mm$model_spec[mi,]
    input_basedir = file.path(mm$cntrl$model_save_directory, thismod$model_name, 'input_files')
    output_basedir = file.path(mm$cntrl$model_save_directory, thismod$model_name)

    #Set subgroup variables, and write confirmatory_subgroup definitions to file, if applicable
    if(thismod$subgroups){
      subcutoff <- thismod$subgroup_thresh
      if(!(is.na(thismod$subgroup_names) || thismod$subgroup_names %in% c('NA','na','n/a','N/A'))){
        cssubgroup_write <- WGC_subgroup(runmodel_fileConn, mi, thismod, mm, input_basedir)
      } else{
        cssubgroup_write <- 'NULL'
      }
    } else{
      subcutoff <- 'NULL'
      cssubgroup_write <- 'NULL'
    }
    #get apriori paths
    paths_write <- WGC_paths(thismod,mm)

    #get exogenous predictors
    exogenous_write <- WGC_exogenous(thismod,mm)

    #get modulatory effects. i.e., multiply one var with another, like task with each node
    modu_write <- WGC_modulatory(thismod,mm)

    #Standardize predictors, i.e. mean 0 and SD 1.
    standardize_write <- WGC_standardize(thismod,mm)

    #Write model specifications for this model to file
    WGC_modelspecs(runmodel_fileConn, mi,input_basedir,output_basedir,
                   thismod, subcutoff,
                   cssubgroup_write,paths_write,exogenous_write,modu_write,standardize_write)
  }

  #Write parallelization header, if applicable
  WGC_parallel_header(runmodel_fileConn, maxcores, nummodels)

  #Write gimme code to file, to makeuse of modelspecs within parallelization loop
  write(c(
    "  gimme::gimme(data = modelspecs[[i]]$data,",
    "      out = modelspecs[[i]]$out,",
    "      sep = modelspecs[[i]]$sep,",
    "      ar = modelspecs[[i]]$ar,",
    "      header = modelspecs[[i]]$header,",
    "      plot = modelspecs[[i]]$plot,",
    "      groupcutoff = modelspecs[[i]]$groupcutoff,",
    "      subgroup = modelspecs[[i]]$subgroup,",
    "      subcutoff = modelspecs[[i]]$subcutoff,",
    "      confirm_subgroup = modelspecs[[i]]$confirm_subgroup,",
    "      paths = modelspecs[[i]]$paths,",
    "      exogenous = modelspecs[[i]]$exogenous,",
    "      standardize = modelspecs[[i]]$standardize",
    "  )"),
    runmodel_fileConn, append = TRUE)

  #close parallel loop if necessary
  WGC_parallel_close(runmodel_fileConn, maxcores)

  close(runmodel_fileConn)
  return(model_filename)
}

WGC_modulatory <- function(thismod,mm){
  modu <- thismod$modulatory_predictors
  if (!(is.na(modu) || modu %in% c('NA','na','n/a','N/A')) && length(modu)>0){
    all_mods <- as.character()
    modulist <- strsplit(gsub('[(]','',modu),')')[[1]]
    for (amod in modulist){
      #if either side is a listname, expand it to the list. Apply shorten to all nodes.
      sides <- strsplit(amod,'[*]')[[1]]
      if(mm$lists[sides[1]]=='NULL'){Lside <- sides[1]}else{Lside <- unlist(mm$lists[sides[1]])}
      if(mm$lists[sides[2]]=='NULL'){Rside <- sides[2]}else{Rside <- unlist(mm$lists[sides[2]])}
      Lside <- applyShorten(Lside,mm$shortnames)
      Rside <- applyShorten(Rside,mm$shortnames)

      #combine each lefthand term with all righthand terms
      for(LS in Lside){
        all_mods <- c(all_mods,sprintf("\'%s*%s\'",LS,Rside))
      }
      #remove entries where the same node is on left and right sides
      duplicateLR <- unlist(lapply(all_mods, function(x){spl <- gsub("\'",'',strsplit(x,'[*]')[[1]]); return(spl[1]==spl[2])}))
      all_mods <- all_mods[!duplicateLR]
    }
    all_mods <- unique(all_mods)
    modu_write <- sprintf('c(%s)',paste(all_mods,collapse = ','))
  } else{modu_write <- 'NULL'}
  return(modu_write)
}

WGC_exogenous <- function(thismod,mm){
  exo <- thismod$exogenous_predictors
  if (!(is.na(exo) || exo %in% c('NA','na','n/a','N/A'))  && length(exo)>0){
    exolist <- strsplit(gsub('[(]','',exo),')')[[1]]
    shorten_exolist <- applyShorten(exolist,mm$shortnames)
    exo_write <- sprintf('c(%s)',paste(sprintf("\'%s\'",shorten_exolist),collapse = ','))
  } else{exo_write <- 'NULL'}
  return(exo_write)
}

WGC_standardize <- function(thismod,mm){
  stan <- thismod$standardize_predictors
  if (stan || stan %in% c('TRUE','true','yes')){
    stan_write <- 'TRUE'
  } else{stan_write <- 'FALSE'}
  return(stan_write)
}



WGC_paths <- function(thismod, mm){
  if (!(is.na(thismod$apriori_paths) || thismod$apriori_paths %in% c('NA','na','n/a','N/A')) && length(thismod$apriori_paths)>0){
    paths <- as.character()
    pathlist <- strsplit(gsub('[(]','',thismod$apriori_paths),')')[[1]]
    L2R <- grep('>',pathlist)
    R2L <- grep('<',pathlist)
    TILDE <- grep('~',pathlist)

    for(i in L2R){
        nodes <- strsplit(pathlist[i],'>')[[1]]
        paths <- c(paths,sprintf("\'%s ~ %s\'",
                                             applyShorten(nodes[2],mm$shortnames),
                                             applyShorten(nodes[1],mm$shortnames)))
    }
    for(i in R2L){
      nodes <- strsplit(pathlist[i],'<')[[1]]
      paths <- c(paths,sprintf("\'%s ~ %s\'",
                                           applyShorten(nodes[1],mm$shortnames),
                                           applyShorten(nodes[2],mm$shortnames)))
    }
    for(i in TILDE){
      nodes <- strsplit(pathlist[i],'~')[[1]]
      paths <- c(paths,sprintf("\'%s ~ %s\'",
                               applyShorten(nodes[1],mm$shortnames),
                               applyShorten(nodes[2],mm$shortnames)))
    }
    paths_write <- sprintf('c(%s)',paste(paths, collapse = ','))
  }else{paths_write <- 'NULL'}
  return(paths_write)
}

WGC_header <- function(runmodel_fileConn, maxcores){
  if(useCores(maxcores)>1){
    write(c(
    'library(doParallel)'),
    runmodel_fileConn, append = TRUE)
  }
  write(c(
    'library(gimme)',
    'modelspecs <- list()',
    'cs_subgroups <- list()'),
    runmodel_fileConn, append = TRUE)
}

WGC_parallel_header <- function(runmodel_fileConn, maxcores, nummodels) {
  numcores <- useCores(maxcores)
  if (numcores == 1){
    fortext <- sprintf('for(i in 1:%s){',nummodels)
  }else{
    fortext <- c(sprintf('registerDoParallel(%s)',numcores),sprintf('foreach(i=1:%s) %%dopar%% {',nummodels))
  }
  write(fortext, runmodel_fileConn,append = TRUE)
}

WGC_parallel_close <- function(runmodel_fileConn, maxcores){
  if (useCores(maxcores) > 1){
    write(c(
      '}'),
      runmodel_fileConn, append = TRUE)
  } else if(useCores(maxcores) == 1){
    write(c(
      '}'),
      runmodel_fileConn, append = TRUE)
  }
}

WGC_subgroup <- function(runmodel_fileConn, mi, thismod, mm, input_basedir){
  #write code to generate cs_subgroup dataframe from filelist to the .R file
  if(!is.na(thismod$subgroup_names) && length(thismod$subgroup_names)>0){

#    grouplist <- mm$lists[[thismod$subgroup_names]]
#    grouplist_string <- paste('grouplist <- c(',
#                              paste(unlist(lapply(grouplist,function(x){return(sprintf("\'%s\'",x))})),collapse=','),
#                              ')',collapse='')
    grouplist <- strsplit(gsub('[(]','',thismod$subgroup_names),')')[[1]]
    write("groupnum_lookup <- list()",runmodel_fileConn, append = TRUE)
    index <- 0
    for(g in grouplist){
      index <- index+1
      write(sprintf('groupnum_lookup$%s <- %s',g,index),runmodel_fileConn, append = TRUE)
    }
#      "index <- 0",
#      "for(g in grouplist){index <- index+1; groupnum_lookup[g] = index}",

    write(c(
      sprintf("cs_subgroups[[%s]] <- data.frame(filename = dir(\'%s\'), groupnum = 0, stringsAsFactors = FALSE)",mi,input_basedir),
      sprintf("cs_subgroups[[%s]]$filename <- unlist(lapply(cs_subgroups[[%s]]$filename, function(x){gsub(\'.csv\',\'\',x)}))",mi,mi),
      sprintf("cs_subgroups[[%s]]$groupnum <- unlist(lapply(cs_subgroups[[%s]]$filename, function(x){groupnum_lookup[[strsplit(x,\'_\')[[1]][1]]]}))",mi,mi)),
      runmodel_fileConn, append = TRUE
    )
    cssubgroup_write <- sprintf("cs_subgroups[[%s]]",mi)
  } else{cssubgroup_write <- 'NULL'}
  return(cssubgroup_write)
}

WGC_modelspecs <- function(runmodel_fileConn,
                           mi,input_basedir,output_basedir,
                           thismod, subcutoff,
                           cssubgroup_write,paths_write,exogenous_write,modu_write,standardize_write){
  write(c(
    sprintf("modelspecs[[%s]] <- list(",mi),
    sprintf("      data = \'%s\',",input_basedir),
    sprintf("      out = \'%s\',",output_basedir),
    "      sep = \',\',",
    "      ar = TRUE,",
    "      plot = TRUE,",
    "      header = TRUE,",
    sprintf("      groupcutoff = %s,",thismod$group_thresh),
    sprintf("      subgroup = %s,",thismod$subgroups),
    sprintf("      subcutoff = %s,",subcutoff),
    sprintf("      confirm_subgroup = %s,",cssubgroup_write),
    sprintf("      paths = %s,",paths_write),
    sprintf("      exogenous = %s,",exogenous_write),
    sprintf("      mult_vars = %s,",modu_write),
    sprintf("      standardize = %s",standardize_write),
    ")"),
    runmodel_fileConn,append = TRUE)
}

