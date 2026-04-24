
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

initializeGimmeFolders <- function(mm, savedir = NA){

  if(is.na(savedir)){
    savedir <- mm$cntrl$model_save_directory
  }
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
    if(toupper(thismod$run) == 'ALL') {
      #all
      rowkeep <- rep(TRUE,length(mm$timecourses$run))
    } else{
      if(class(thismod$run) == 'character'){
        #character. search for each entry within parentheses
      spl <- strsplit(thismod$run,'[()]')[[1]]
      run_name_lengths <- unlist(lapply(spl,nchar))
      spl <- spl[run_name_lengths > 0]
      eachrun <- spl[! spl %in% c('(',')')]
      rowkeep <- rep(FALSE,length(mm$timecourses$run))
      for(thisrun in eachrun){
        rowkeep <- rowkeep | mm$timecourses$run == thisrun
        }
      } else {
        #numeric
        rowkeep <- mm$timecourses$run == thismod$run
        }
      }
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
        if(thismod$apriori_subgroups){
          prefix <- sprintf('%s_',thisgroup)
        } else{
          prefix <- ''
        }
        outfilename <- file.path(savedir,thismod$model_name,'input_files',sprintf('%s%s.csv',prefix,s))


        write.csv(thissub[,short_roicols], outfilename, row.names = FALSE)
      }
    }

  }

}

doreplace <- function(line, pattern_list, target){
  for(pattern in pattern_list){
    if(!is.na(line[pattern])){
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
  config_sub_list <- unique(config$ID[config$TYPE == 'sub'])
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
    allruns <- config[config$ID == s,]
    for(ll in 1:dim(allruns)[1]){
    line <- allruns[ll,]
    uid <- paste(line$ID,line$RUN,sep = '_')
      filelocs[[uid]] <- getTimecourse_OneSub(xtc_fileConn, s, line$RUN, config[config$ID %in% roi_list,],
                                            line$DATA_DIR,
                                            line$TIMECOURSE_DIR,
                                            line$FUNCTIONAL_FILENAME,
                                            line$MASK_FILENAME,
                                            submaskdir = line$DATA_DIR,
                                            subroi_dir = line$SUB_ROI_DIR,
                                            remove_subrois = FALSE
                                            )
    }
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

genTimecoursesCSV <- function(tcfilename, filelocs, config = NA, verbose = FALSE){
  if(!(typeof(config)=='logical' && is.na(config))){
    config$UID <- paste(config$ID,config$RUN,sep = '_')
  }

  allrois <- as.character()
    for(s in names(filelocs)){
      allrois <- c(allrois,names(filelocs[[s]]))
    }
    allrois <- unique(allrois)

    eachdf <- list()

    for(s in names(filelocs)){
      if( !(typeof(config)=='logical' && is.na(config)) ){
        thissub <- unlist(config[config$UID == s,'ID'])
        group <- unlist(config[config$UID == s,'GROUP'])
        thisrun <- unlist(config[config$UID == s,'RUN'])
        thiscensor <- unlist(config[config$UID == s,'CENSOR_FILENAME'])
      } else{
        thissub <- NA
        group <- NA
        thisrun <- NA
        thiscensor <- NA
      }
      for(r in names(filelocs[[s]])){
        thistc_filename <- filelocs[[s]][[r]]$timecourse_loc
        if(!is.na(file.info(thistc_filename)$size) && file.info(thistc_filename)$size > 0){
          if(verbose){message(sprintf('Reading %s',thistc_filename))}
          didreadtc <- FALSE
          tryCatch(
            {thistc <- read.csv(thistc_filename,header = FALSE)
            didreadtc <- TRUE
          }, error = function(e){message('error reading ',thistc_filename)})
          if(didreadtc){

          # initialize this subject's timecourse df if it doesn't exist
          if(! s %in% names(eachdf)){
            tclength <- dim(thistc)[1]
            eachdf[[s]]  <- data.frame(slicenum = 1:tclength, time = rep(NA,1,tclength), condition = rep(NA,1,tclength), censor = rep(0,1,tclength), subnum = rep(thissub,1,tclength), subgroup = rep(group,1,tclength), run = rep(thisrun,1,tclength), stringsAsFactors = FALSE)
            eachdf[[s]][,allrois] <- NA
            if(!is.na(thiscensor)){
              intext <- readLines(thiscensor)
              if(startsWith(intext[1],'-CENSORTR')){
               cens <- rep(1,tclength)
               censtext <- gsub('-CENSORTR ','',intext)
               # Split the input text by space and comma to separate chunks
               chunks <- strsplit(censtext, "[, ]+")[[1]]
               run_prefix <- sprintf('%s:',thisrun)
               chunks_thisrun <- chunks[grepl(run_prefix,chunks)]
               chunks_thisrun <- gsub(run_prefix,'',chunks_thisrun)

               for(index in chunks_thisrun){
                 if (grepl("\\.\\.", index)) {
                   # If it's a range, expand it
                   range_parts <- as.integer(unlist(strsplit(index, "\\.\\.")))
                   cens[(range_parts[1]+1):(range_parts[2]+1)] <- 0
                 } else {
                   # If it's a single index, mark it
                   cens[as.integer(index)+1] <- 0
                 }
               }
               eachdf[[s]]$censor <- cens


              } else {

              if(verbose){message(sprintf('Reading %s',thiscensor))}
                if(file.size(thiscensor) > 1){
                  cens <- read.csv(thiscensor, header = FALSE)
                } else {
                  #if the censor file is empty, assume that nothing is censored
                  warning(sprintf('Censor File %s is empty. Applying no censoring.',thiscensor))
                  cens <- as.matrix(rep(1,tclength))
                }
              eachdf[[s]]$censor <- cens[1:tclength,1]
              }

            }
          }
          eachdf[[s]][,r] <- thistc
        }
        }
      }
    }
    alldf <- do.call('rbind',eachdf)
    write.csv(alldf, file = tcfilename,row.names = FALSE)
  }

getTimecourse_OneSub <- function(fileConn, subname, runname, roi_df,
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
    filelocs[[roiname]] <- getTimecourse_OneSub_OneROI(fileConn, subname, runname, roiname,
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

getTimecourse_OneSub_OneROI <- function(fileConn, subname, runname, roiname,
                                        subdir,
                                        roidir,
                                        timecoursedir,
                                        subfunc_loc, sub_mask_loc, roi_loc,
                                        submaskdir = subdir, subroi_dir = file.path(subdir,'subject_rois'), remove_subrois = FALSE){
  #resample the roi to the subjects functional
#  subfunc_loc <-  file.path(subdir,sub_functional_filename)
  roisub_loc <- file.path(subroi_dir,sprintf('SUB_%s_RUN_%s_ROI_%s+tlrc',subname, runname, roiname))
  dir.create(file.path(subroi_dir), showWarnings = FALSE)
  dir.create(file.path(timecoursedir), showWarnings = FALSE)
#  roi_loc <- file.path(roidir, roi_filename)
  line_3dresample <- sprintf('3dresample -master %s -prefix %s -inset %s -overwrite',subfunc_loc, roisub_loc, roi_loc)

  #combine subroi with the subject mask
#  sub_mask_loc <- file.path(submaskdir, sub_mask_filename)
  line_3dcalc <- sprintf("3dcalc -a %s -b %s -expr \'a*b\' -nscale -overwrite -prefix %s", roisub_loc, sub_mask_loc, roisub_loc)

  #save timecourse in this mask to text
  timecourse_loc <- file.path(timecoursedir,  sprintf('SUB_%s_RUN_%s_ROI_%s_timecourse.txt',subname, runname, roiname))
  line_3dmaskave <- sprintf("3dmaskave -mask %s -quiet %s > %s",roisub_loc, subfunc_loc, timecourse_loc)


  if(!is.na(fileConn)){write(c(line_3dresample, line_3dcalc, line_3dmaskave), fileConn, append = TRUE)}
  return(list(subfunc_loc = subfunc_loc, roisub_loc = roisub_loc, roi_loc = roi_loc, timecourse_loc = timecourse_loc))
}


#' Compute Network Metrics from GIMME Model Output (Internal Implementation)
#'
#' Computes node-level and network-level graph metrics for individual networks
#' based on GIMME model output. Saves individual subject files and a combined summary.
#' This is the internal implementation. Use the wrapper in interface.R for exported function.
#'
#' @param model_dir Path to GIMME model output folder (should contain indivPathEstimates.csv)
#' @param ignore_lags Logical. If TRUE (default), excludes lagged connections from analysis
#' @param save_individual Logical. If TRUE (default), saves individual network metrics to individual/ folder
#' @param save_summary Logical. If TRUE (default), saves combined summary file to model_dir
#' @param verbose Logical. If TRUE (default), prints progress messages
#'
#' @return A list containing:
#'   \item{node_metrics}{Data frame with node-level metrics for all subjects}
#'   \item{network_metrics}{Data frame with network-level metrics for all subjects}
computeNetworkMetrics_internal <- function(model_dir,
                                   ignore_lags = TRUE,
                                   save_individual = TRUE,
                                   save_summary = TRUE,
                                   verbose = TRUE) {

  # Check if igraph is available
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required but not installed. Please install it with: install.packages('igraph')")
  }

  # Read indivPathEstimates.csv
  path_file <- file.path(model_dir, "indivPathEstimates.csv")
  if (!file.exists(path_file)) {
    stop(sprintf("File not found: %s", path_file))
  }

  if (verbose) message(sprintf("Reading path estimates from: %s", path_file))
  path_data <- read.csv(path_file, stringsAsFactors = FALSE)

  # Filter out lagged connections if requested
  if (ignore_lags) {
    # Check if lhs or rhs contains "lag"
    original_rows <- nrow(path_data)
    path_data <- path_data[!grepl("lag", path_data$lhs, ignore.case = TRUE) &
                           !grepl("lag", path_data$rhs, ignore.case = TRUE), ]
    if (verbose) message(sprintf("Filtered out %d lagged connections", original_rows - nrow(path_data)))
  }

  # Get unique subject-condition combinations
  subjects <- unique(path_data$file)
  if (verbose) message(sprintf("Found %d unique subject-condition combinations", length(subjects)))

  # Initialize storage for results
  all_node_metrics <- list()
  all_network_metrics <- list()

  # Process each subject-condition
  for (subj in subjects) {
    if (verbose) message(sprintf("Processing: %s", subj))

    # Extract edges for this subject
    subj_edges <- path_data[path_data$file == subj, ]

    # Parse subject and condition
    parsed <- parseSubjectCondition(subj)

    # Build network and compute metrics
    metrics <- computeSubjectNetworkMetrics(subj_edges, parsed$file, parsed$subject, parsed$condition)

    # Save individual file if requested
    if (save_individual && !is.null(metrics$node_metrics)) {
      individual_dir <- file.path(model_dir, "individual")
      if (!dir.exists(individual_dir)) {
        dir.create(individual_dir, showWarnings = FALSE)
      }
      individual_file <- file.path(individual_dir, sprintf("%sNetworkMetrics.csv", parsed$file))
      write.csv(metrics$node_metrics, individual_file, row.names = FALSE)
      if (verbose) message(sprintf("  Saved individual metrics to: %s", individual_file))
    }

    # Store results
    all_node_metrics[[subj]] <- metrics$node_metrics
    all_network_metrics[[subj]] <- metrics$network_metrics
  }

  # Combine all results
  combined_node_metrics <- do.call(rbind, all_node_metrics)
  combined_network_metrics <- do.call(rbind, all_network_metrics)

  # Save summary file if requested
  if (save_summary) {
    summary_file <- file.path(model_dir, "networkMetrics_summary.csv")
    write.csv(combined_network_metrics, summary_file, row.names = FALSE)
    if (verbose) message(sprintf("Saved summary metrics to: %s", summary_file))
  }

  return(list(
    node_metrics = combined_node_metrics,
    network_metrics = combined_network_metrics
  ))
}


#' Parse Subject and Condition from File Identifier
#'
#' @param file_id Character string like "2003_A" or "2003"
#' @return List with file, subject, and condition
parseSubjectCondition <- function(file_id) {
  parts <- strsplit(file_id, "_")[[1]]

  if (length(parts) > 1) {
    # Has condition
    subject <- parts[1]
    condition <- parts[2]
  } else {
    # No condition
    subject <- parts[1]
    condition <- NA
  }

  return(list(
    file = file_id,
    subject = subject,
    condition = condition
  ))
}


#' Compute Network Metrics for a Single Subject
#'
#' @param edges Data frame with columns: lhs, rhs, beta
#' @param file_id Full file identifier (e.g., "2003_A")
#' @param subject Subject ID
#' @param condition Condition (or NA)
#' @return List with node_metrics and network_metrics data frames
computeSubjectNetworkMetrics <- function(edges, file_id, subject, condition) {

  # Check if there are any edges
  if (nrow(edges) == 0) {
    warning(sprintf("No edges found for %s. Returning NA metrics.", file_id))
    return(list(
      node_metrics = data.frame(
        file = file_id,
        subject = subject,
        condition = condition,
        node = NA,
        in_degree = NA,
        out_degree = NA,
        in_strength = NA,
        out_strength = NA,
        betweenness = NA,
        clustering = NA
      ),
      network_metrics = data.frame(
        file = file_id,
        subject = subject,
        condition = condition,
        n_edges = 0,
        density = NA,
        mean_strength = NA,
        total_strength = NA,
        global_efficiency = NA,
        mean_clustering = NA,
        modularity = NA
      )
    ))
  }

  # Get all unique nodes
  all_nodes <- unique(c(edges$lhs, edges$rhs))

  # Build directed weighted graph
  g <- igraph::graph_from_data_frame(
    d = edges[, c("rhs", "lhs", "beta")],  # Note: rhs -> lhs based on your description
    directed = TRUE,
    vertices = all_nodes
  )

  # Rename edge attribute to 'weight'
  igraph::E(g)$weight <- abs(igraph::E(g)$beta)  # Use absolute values for weights
  g <- igraph::delete_edge_attr(g, "beta")

  # Compute node-level metrics
  node_metrics <- data.frame(
    file = file_id,
    subject = subject,
    condition = condition,
    node = igraph::V(g)$name,
    in_degree = igraph::degree(g, mode = "in"),
    out_degree = igraph::degree(g, mode = "out"),
    in_strength = igraph::strength(g, mode = "in"),
    out_strength = igraph::strength(g, mode = "out"),
    betweenness = tryCatch(igraph::betweenness(g, directed = TRUE, weights = NA),
                           error = function(e) rep(NA, igraph::vcount(g))),
    clustering = tryCatch(igraph::transitivity(g, type = "local", isolates = "zero"),
                          error = function(e) rep(NA, igraph::vcount(g)))
  )

  # Compute network-level metrics
  n_edges <- igraph::ecount(g)
  n_nodes <- igraph::vcount(g)
  density <- igraph::edge_density(g)
  total_strength <- sum(igraph::E(g)$weight)
  mean_strength <- mean(igraph::E(g)$weight)

  # Global efficiency
  global_efficiency <- tryCatch({
    # Compute shortest paths with weights (using 1/weight as distance)
    if (n_edges > 0 && all(igraph::E(g)$weight > 0)) {
      distances <- igraph::distances(g, weights = 1/igraph::E(g)$weight, mode = "out")
      distances[is.infinite(distances)] <- 0
      if (n_nodes > 1) {
        sum(1/distances[distances > 0]) / (n_nodes * (n_nodes - 1))
      } else {
        NA
      }
    } else {
      NA
    }
  }, error = function(e) NA)

  # Mean clustering coefficient
  mean_clustering <- tryCatch({
    local_cc <- igraph::transitivity(g, type = "local", isolates = "zero")
    mean(local_cc, na.rm = TRUE)
  }, error = function(e) NA)

  # Modularity (using Louvain algorithm)
  modularity <- tryCatch({
    # Convert to undirected for modularity
    g_undir <- igraph::as.undirected(g, mode = "collapse")
    clusters <- igraph::cluster_louvain(g_undir, weights = igraph::E(g_undir)$weight)
    igraph::modularity(clusters)
  }, error = function(e) NA)

  network_metrics <- data.frame(
    file = file_id,
    subject = subject,
    condition = condition,
    n_edges = n_edges,
    density = density,
    mean_strength = mean_strength,
    total_strength = total_strength,
    global_efficiency = global_efficiency,
    mean_clustering = mean_clustering,
    modularity = modularity
  )

  return(list(
    node_metrics = node_metrics,
    network_metrics = network_metrics
  ))
}


#' Compare Network Conditions (Internal Implementation)
#'
#' Computes within-subject network comparison metrics between two conditions
#' (e.g., ON vs OFF medication). Handles both single-model and dual-model cases.
#'
#' @param model_dir_A Path to first condition model (or single model with both)
#' @param model_dir_B Path to second condition (NULL if single model)
#' @param condition_A Condition identifier for first network (default "A")
#' @param condition_B Condition identifier for second network (default "B")
#' @param ignore_lags Logical. If TRUE (default), excludes lagged connections
#' @param save_output Logical. If TRUE (default), saves output files
#' @param output_dir Where to save results (NULL for auto-generated path)
#' @param verbose Logical. If TRUE (default), prints progress messages
#'
#' @return A list containing comparison_summary and individual_comparisons data frames
compareNetworkConditions_internal <- function(model_dir_A,
                                              model_dir_B = NULL,
                                              condition_A = "A",
                                              condition_B = "B",
                                              ignore_lags = TRUE,
                                              save_output = TRUE,
                                              output_dir = NULL,
                                              verbose = TRUE) {
  
  # Determine output directory
  if (is.null(output_dir)) {
    if (!is.null(model_dir_B)) {
      # Dual model case: create comparison folder in parent of model_dir_A
      parent_dir <- dirname(model_dir_A)
      dir_A_name <- basename(model_dir_A)
      dir_B_name <- basename(model_dir_B)
      output_dir <- file.path(parent_dir, sprintf("comparison_%s_vs_%s", dir_A_name, dir_B_name))
    } else {
      # Single model case: create comparison folder in parent
      parent_dir <- dirname(model_dir_A)
      dir_name <- basename(model_dir_A)
      output_dir <- file.path(parent_dir, sprintf("comparison_%s_%s_vs_%s", dir_name, condition_A, condition_B))
    }
  }
  
  if (verbose) message(sprintf("Output directory: %s", output_dir))
  
  # Read path estimates from both conditions
  if (!is.null(model_dir_B)) {
    # Dual model case
    if (verbose) message("Reading from dual model folders...")
    path_file_A <- file.path(model_dir_A, "indivPathEstimates.csv")
    path_file_B <- file.path(model_dir_B, "indivPathEstimates.csv")
    
    if (!file.exists(path_file_A)) stop(sprintf("File not found: %s", path_file_A))
    if (!file.exists(path_file_B)) stop(sprintf("File not found: %s", path_file_B))
    
    data_A <- read.csv(path_file_A, stringsAsFactors = FALSE)
    data_B <- read.csv(path_file_B, stringsAsFactors = FALSE)
    
  } else {
    # Single model case
    if (verbose) message("Reading from single model folder...")
    path_file <- file.path(model_dir_A, "indivPathEstimates.csv")
    
    if (!file.exists(path_file)) stop(sprintf("File not found: %s", path_file))
    
    all_data <- read.csv(path_file, stringsAsFactors = FALSE)
    
    # Split by condition
    data_A <- all_data[grepl(sprintf("_%s$", condition_A), all_data$file), ]
    data_B <- all_data[grepl(sprintf("_%s$", condition_B), all_data$file), ]
  }
  
  # Filter out lagged connections if requested
  if (ignore_lags) {
    original_A <- nrow(data_A)
    original_B <- nrow(data_B)
    data_A <- data_A[!grepl("lag", data_A$lhs, ignore.case = TRUE) & 
                     !grepl("lag", data_A$rhs, ignore.case = TRUE), ]
    data_B <- data_B[!grepl("lag", data_B$lhs, ignore.case = TRUE) & 
                     !grepl("lag", data_B$rhs, ignore.case = TRUE), ]
    if (verbose) message(sprintf("Filtered %d lagged connections from condition A, %d from condition B",
                                 original_A - nrow(data_A), original_B - nrow(data_B)))
  }
  
  # Extract subject IDs from both conditions
  subjects_A <- unique(data_A$file)
  subjects_B <- unique(data_B$file)
  
  # Parse to get base subject IDs
  parse_subjects <- function(file_ids, condition) {
    base_ids <- sapply(file_ids, function(x) {
      parsed <- parseSubjectCondition(x)
      return(parsed$subject)
    })
    return(data.frame(file = file_ids, subject = base_ids, stringsAsFactors = FALSE))
  }
  
  subj_map_A <- parse_subjects(subjects_A, condition_A)
  subj_map_B <- parse_subjects(subjects_B, condition_B)
  
  # Get all unique subject IDs
  all_subjects <- unique(c(subj_map_A$subject, subj_map_B$subject))
  
  if (verbose) message(sprintf("Found %d subjects in condition A, %d in condition B, %d total unique",
                               nrow(subj_map_A), nrow(subj_map_B), length(all_subjects)))
  
  # Initialize storage
  comparison_results <- list()
  individual_comparisons <- list()
  
  # Process each subject
  for (subj_id in all_subjects) {
    if (verbose) message(sprintf("Processing subject: %s", subj_id))
    
    # Get file IDs for this subject in each condition
    file_A <- subj_map_A$file[subj_map_A$subject == subj_id]
    file_B <- subj_map_B$file[subj_map_B$subject == subj_id]
    
    # Check if subject exists in both conditions
    has_A <- length(file_A) > 0
    has_B <- length(file_B) > 0
    
    if (!has_A && !has_B) next  # Skip if somehow neither (shouldn't happen)
    
    # Compute comparison metrics
    comp_result <- compareSubjectNetworks(
      subj_id = subj_id,
      edges_A = if (has_A) data_A[data_A$file == file_A[1], ] else NULL,
      edges_B = if (has_B) data_B[data_B$file == file_B[1], ] else NULL,
      condition_A = condition_A,
      condition_B = condition_B
    )
    
    comparison_results[[subj_id]] <- comp_result$summary
    individual_comparisons[[subj_id]] <- comp_result$details
    
    # Save individual comparison file if requested
    if (save_output && !is.null(comp_result$details)) {
      individual_dir <- file.path(output_dir, "individual")
      dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)
      individual_file <- file.path(individual_dir, sprintf("%s_comparison.csv", subj_id))
      write.csv(comp_result$details, individual_file, row.names = FALSE)
      if (verbose) message(sprintf("  Saved individual comparison to: %s", individual_file))
    }
  }
  
  # Combine all results
  comparison_summary <- do.call(rbind, comparison_results)
  
  # Save summary file if requested
  if (save_output) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    summary_file <- file.path(output_dir, "comparison_summary.csv")
    write.csv(comparison_summary, summary_file, row.names = FALSE)
    if (verbose) message(sprintf("Saved comparison summary to: %s", summary_file))
  }
  
  return(list(
    comparison_summary = comparison_summary,
    individual_comparisons = individual_comparisons
  ))
}


#' Compare Networks for a Single Subject Across Conditions
#'
#' @param subj_id Subject ID
#' @param edges_A Data frame with edges for condition A (or NULL if missing)
#' @param edges_B Data frame with edges for condition B (or NULL if missing)
#' @param condition_A Name of condition A
#' @param condition_B Name of condition B
#' @return List with summary and details data frames
compareSubjectNetworks <- function(subj_id, edges_A, edges_B, condition_A, condition_B) {
  
  # Handle missing conditions
  has_A <- !is.null(edges_A) && nrow(edges_A) > 0
  has_B <- !is.null(edges_B) && nrow(edges_B) > 0
  
  if (!has_A && !has_B) {
    # Neither condition exists
    return(list(
      summary = data.frame(
        subject = subj_id,
        has_condition_A = FALSE,
        has_condition_B = FALSE,
        jaccard_similarity = NA,
        edge_overlap_count = NA,
        edge_overlap_pct = NA,
        edges_only_A = NA,
        edges_only_B = NA,
        strength_correlation = NA,
        mean_strength_diff = NA,
        mean_strength_A = NA,
        mean_strength_B = NA,
        sd_strength_A = NA,
        sd_strength_B = NA
      ),
      details = NULL
    ))
  }
  
  # Create edge identifiers (from -> to)
  if (has_A) {
    edges_A$edge_id <- paste(edges_A$rhs, edges_A$lhs, sep = "->")
    edge_set_A <- unique(edges_A$edge_id)
    weights_A <- setNames(abs(edges_A$beta), edges_A$edge_id)
  } else {
    edge_set_A <- character(0)
    weights_A <- numeric(0)
  }
  
  if (has_B) {
    edges_B$edge_id <- paste(edges_B$rhs, edges_B$lhs, sep = "->")
    edge_set_B <- unique(edges_B$edge_id)
    weights_B <- setNames(abs(edges_B$beta), edges_B$edge_id)
  } else {
    edge_set_B <- character(0)
    weights_B <- numeric(0)
  }
  
  # Compute set operations
  all_edges <- unique(c(edge_set_A, edge_set_B))
  shared_edges <- intersect(edge_set_A, edge_set_B)
  only_A <- setdiff(edge_set_A, edge_set_B)
  only_B <- setdiff(edge_set_B, edge_set_A)
  
  # Compute metrics
  n_shared <- length(shared_edges)
  n_union <- length(all_edges)
  n_only_A <- length(only_A)
  n_only_B <- length(only_B)
  
  jaccard_similarity <- if (n_union > 0) n_shared / n_union else NA
  edge_overlap_pct <- if (n_union > 0) n_shared / n_union * 100 else NA
  
  # Strength statistics
  mean_strength_A <- if (has_A) mean(abs(edges_A$beta), na.rm = TRUE) else NA
  mean_strength_B <- if (has_B) mean(abs(edges_B$beta), na.rm = TRUE) else NA
  sd_strength_A <- if (has_A) sd(abs(edges_A$beta), na.rm = TRUE) else NA
  sd_strength_B <- if (has_B) sd(abs(edges_B$beta), na.rm = TRUE) else NA
  
  # Strength correlation and difference (for shared edges only)
  if (n_shared > 0) {
    shared_weights_A <- weights_A[shared_edges]
    shared_weights_B <- weights_B[shared_edges]
    strength_correlation <- cor(shared_weights_A, shared_weights_B, use = "complete.obs")
    mean_strength_diff <- mean(abs(shared_weights_A - shared_weights_B), na.rm = TRUE)
  } else {
    strength_correlation <- NA
    mean_strength_diff <- NA
  }
  
  # Create summary
  summary <- data.frame(
    subject = subj_id,
    has_condition_A = has_A,
    has_condition_B = has_B,
    jaccard_similarity = jaccard_similarity,
    edge_overlap_count = n_shared,
    edge_overlap_pct = edge_overlap_pct,
    edges_only_A = n_only_A,
    edges_only_B = n_only_B,
    strength_correlation = strength_correlation,
    mean_strength_diff = mean_strength_diff,
    mean_strength_A = mean_strength_A,
    mean_strength_B = mean_strength_B,
    sd_strength_A = sd_strength_A,
    sd_strength_B = sd_strength_B
  )
  
  # Create detailed edge-by-edge comparison
  if (length(all_edges) > 0) {
    details <- data.frame(
      subject = subj_id,
      edge = all_edges,
      edge_from = sapply(strsplit(all_edges, "->"), `[`, 1),
      edge_to = sapply(strsplit(all_edges, "->"), `[`, 2),
      in_A = all_edges %in% edge_set_A,
      in_B = all_edges %in% edge_set_B,
      weight_A = ifelse(all_edges %in% edge_set_A, weights_A[all_edges], NA),
      weight_B = ifelse(all_edges %in% edge_set_B, weights_B[all_edges], NA),
      weight_diff = NA,
      stringsAsFactors = FALSE
    )
    
    # Compute weight difference for shared edges
    shared_idx <- details$in_A & details$in_B
    details$weight_diff[shared_idx] <- abs(details$weight_A[shared_idx] - details$weight_B[shared_idx])
  } else {
    details <- NULL
  }
  
  return(list(
    summary = summary,
    details = details
  ))
}


#' Correlate GIMME Edge Weights with a Behavioral Outcome (Internal)
#'
#' For each group-level or subgroup-level non-lag edge, correlates individual
#' path beta weights with a behavioral outcome, separately per subgroup defined
#' in the participants file.
#'
#' @param model_dir Path to GIMME model output folder containing indivPathEstimates.csv
#' @param participants_file Path to participants tsv/csv (NULL = look for
#'   participants.tsv in model_dir)
#' @param behavioral_col Column in participants file to use as outcome. NULL =
#'   auto-detect.
#' @param save_output Logical. Save correlation summary CSV to model_dir.
#' @param verbose Logical. Print progress messages.
#'
#' @return A list with:
#'   \item{correlations}{Data frame of correlation results, one row per edge x subgroup}
#'   \item{edge_data}{Wide data frame used for correlations (subjects x edges + behavioral)}
#'   \item{behavioral_col}{The behavioral column used}
correlateBehavior_internal <- function(model_dir,
                                       participants_file = NULL,
                                       behavioral_col = NULL,
                                       save_output = TRUE,
                                       verbose = TRUE) {

  # Load path estimates
  path_file <- file.path(model_dir, "indivPathEstimates.csv")
  if (!file.exists(path_file)) stop(sprintf("File not found: %s", path_file))
  paths <- utils::read.csv(path_file, stringsAsFactors = FALSE)

  # Keep only group and subgroup-level non-lag edges
  paths <- paths[!grepl("lag", paths$rhs, ignore.case = TRUE) &
                 !grepl("lag", paths$lhs, ignore.case = TRUE), ]
  paths <- paths[paths$level %in% c("group", "subgroup1", "subgroup2"), ]

  if (nrow(paths) == 0) stop("No group/subgroup non-lag edges found in indivPathEstimates.csv")

  # Identify group/subgroup edges (those that appear for >1 subject at that level)
  paths$edge <- paste(paths$rhs, paths$lhs, sep = "->")

  # Wide format: one row per subject, one column per edge
  edge_wide <- reshape(
    paths[, c("file", "edge", "beta")],
    idvar = "file",
    timevar = "edge",
    direction = "wide",
    v.names = "beta"
  )
  names(edge_wide) <- sub("^beta\\.", "", names(edge_wide))

  # Load participants file
  ptcp_info <- loadParticipantsFile(
    participants_file = participants_file,
    gimme_dir = model_dir,
    behavioral_col = behavioral_col
  )
  ptcp <- ptcp_info$data
  behavioral_col <- ptcp_info$behavioral_col

  if (verbose) message(sprintf("Using behavioral column: '%s'", behavioral_col))

  # Match participant IDs between GIMME output and participants file
  # GIMME file column may be e.g. "101" while participant_id is "sub-101"
  # Try stripping "sub-" prefix for matching
  ptcp$match_id <- sub("^sub-", "", ptcp$participant_id)
  edge_wide$match_id <- sub("^sub-", "", edge_wide$file)

  merged <- merge(edge_wide, ptcp[, c("match_id", "group", behavioral_col)],
                  by = "match_id", all.x = FALSE)

  if (nrow(merged) == 0) {
    stop(paste(
      "No subjects matched between indivPathEstimates.csv and participants file.",
      sprintf("GIMME subjects (first 5): %s", paste(head(edge_wide$file, 5), collapse = ", ")),
      sprintf("Participants file IDs (first 5): %s", paste(head(ptcp$participant_id, 5), collapse = ", ")),
      sep = "\n"
    ))
  }

  if (verbose) message(sprintf("Matched %d subjects", nrow(merged)))

  # Identify edge columns (everything except file, match_id, group, behavioral_col)
  non_edge_cols <- c("file", "match_id", "group", behavioral_col)
  edge_cols <- setdiff(names(merged), non_edge_cols)

  # Run correlations per subgroup x edge
  subgroups <- sort(unique(merged$group))
  results <- list()

  for (sg in subgroups) {
    sg_data <- merged[merged$group == sg, ]
    for (ec in edge_cols) {
      vals <- sg_data[[ec]]
      beh  <- sg_data[[behavioral_col]]
      # Only non-NA pairs
      ok <- !is.na(vals) & !is.na(beh)
      n_valid <- sum(ok)
      if (n_valid < 3) {
        results[[paste(sg, ec, sep = "|||")]] <- data.frame(
          subgroup = sg, edge = ec, n = n_valid,
          r = NA_real_, r2 = NA_real_, pval = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }
      ct <- stats::cor.test(vals[ok], beh[ok], method = "pearson")
      results[[paste(sg, ec, sep = "|||")]] <- data.frame(
        subgroup = sg, edge = ec, n = n_valid,
        r = ct$estimate, r2 = ct$estimate^2, pval = ct$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  correlations <- do.call(rbind, results)
  rownames(correlations) <- NULL

  if (save_output) {
    out_file <- file.path(model_dir, sprintf("behaviorCorrelations_%s.csv", behavioral_col))
    utils::write.csv(correlations, out_file, row.names = FALSE)
    if (verbose) message(sprintf("Saved correlation results to: %s", out_file))
  }

  return(list(
    correlations = correlations,
    edge_data = merged,
    behavioral_col = behavioral_col
  ))
}
