#Todo:
# allowing multiple runs e.g. (1)(2) isn't working
# allow standardization of inputs. Detect if inputs are too different and prompt to standardize.
# convert a condition variable that's in text to a factor
# create additional figures, e.g. weights scatter plot


#' @export
gimmefMRI_templates <- function(writedir = getwd()){
  file.copy(
    system.file('extdata','get_timecourses.csv', package = 'gimmefMRI'),
    file.path(writedir,'get_timecourses.csv'))
  file.copy(
    system.file('extdata','DemoGIMME.xlsx', package = 'gimmefMRI'),
    file.path(writedir,'DemoGIMME.xlsx'))
  }

#' @export
getTC <- function(
    write_file = 'extract_timecourses.sh',
    config_file = 'gui',
    sub_list = 'all',
    roi_list = 'all',
    timecourse_filename = 'timecourses.csv',
    run = c('generate_timecourse_code','run_timecourse_code','run_timecourse_csv','write_timecourse_xlsx'),
    filelocs = NA,
    config = NA,
    verbose = FALSE
    ){

  if('generate_timecourse_code' %in% run){
    #write afni code to extract timecourses for each subject from each roi
    returndat <- getTimecourse(write_file = write_file, config_file = config_file, sub_list = sub_list, roi_list = roi_list)
  }
  if('run_timecourse_code' %in% run){
    #run the shell file
    if(any(grep('/',write_file))){
        system(sprintf('%s',write_file))
      } else {
        system(sprintf('./%s',write_file))
      }
  }
  if('run_timecourse_csv' %in% run){
    #read in all those tc files and save them to the big timecourses.csv to use as input to gimmefMRI()
    if(!('generate_timecourse_code' %in% run)){
      returndat <- getTimecourse(write_file = NA, config_file = config_file, sub_list = sub_list, roi_list = roi_list)
    }
    genTimecoursesCSV(timecourse_filename, returndat$filelocs, returndat$config,verbose)
  }
  if('write_timecourse_xlsx' %in% run){
    #copy XLSX_FILENAME to sprintf('backup_%s',XLSX_FILENAME)
    #my_timecourses <- read.csv(CSV_FILE)
    #all_sheets <- list()
    #sheetnames <- readxl::excel_sheets(XLSX_FILE)
    #for sheet in sheetnames
    # all_sheets[[sheet]] <- readxl::readxl(XLSX_FILENAME, sheet = sheet)
    #end
    #all_sheets$timecourses <- my_timecourses
    #writexl::write_xlsx(all_sheets, XLSX_FILENAMEF)

  }

  return(returndat)
}

#' @export
TCCrossCorr <- function(
  timecourse_file = 'gui',
  cormat_filename = 'timecourses_cormat.csv',
  sub_list = 'all',
  roi_list = 'all',
  censor_values = 0
  ){
  #Show warnings during execution. Revert to user setting on exit.
  oldwarn <- getOption("warn")
  on.exit(options(warn = oldwarn))
  options(warn = 1)

  #Get filename and load timecourse file
  myfile <- UTILS.getFilename(timecourse_file)
  tc <- read.csv(myfile)
  tc$uid <- paste(tc$subnum, tc$run, sep = '_')


  #Process roi_list, warn if requested rois not present.
  standard_cols <- c('slicenum','time','condition','censor','subnum','subgroup','run')
  file_rois <- names(tc)[! names(tc) %in% standard_cols]
  if(roi_list == 'all'){
    roi_list <- file_rois
  } else {
    if(any(!roi_list %in% file_rois)){
      warning(sprintf('These requested ROIs were not found in %s and will be skipped: %s', myfile, paste(roi_list[!roi_list %in% file_rois]  , collapse = ', ')))
    }
    roi_list <- roi_list[roi_list %in% file_rois]
  }

  #Process sub_list, warn if requested participants not present.
  file_ptcp <- unique(tc$subnum)
  if(sub_list == 'all'){
    sub_list <- file_pctp
  } else {
    if(any(!sub_list %in% file_pctp)){
      warning(sprintf('These requested Participants were not found in %s and will be skipped: %s', myfile, paste(sub_list[!sub_list %in% file_ptcp]  , collapse = ', ')))
    }
    sub_list <- sub_list[sub_list %in% file_ptcp]
  }

  #Censor requested values from timecourses
  if(censor){
    tc[tc$censor %in% censor_values, roi_list] <- NA
  }

  each_cormat <- list()
  for(u in unique(tc$uid)){
    thistc <- subset(tc, uid == u)
    this_cormat <- cor(thistc[,roi_list], use = 'pairwise.complete.obs')
    ut <- upper.tri(this_cormat, diag = FALSE)

    cordf <- as.data.frame(cbind(roi1 = rownames(this_cormat)[row(this_cormat)[ut]],
                                 roi2 = colnames(this_cormat)[col(this_cormat)[ut]],
                                 pearson_r = c(this_cormat[ut])))
    cordf$connection <- paste(cordf$roi1,cordf$roi2, sep = '~')

    cordf <- cbind(cordf, condition = thistc$condition[1], subnum = thistc$subnum[1], subgroup = thistc$subgroup[1], run = thistc$run[1], uid = thistc$uid[1])
    cordf$z <- atanh(as.numeric(cordf$pearson_r)) #Fisher's Z transformation: Inverse hyperbolic tangent function
    each_cormat[[u]] <- cordf[,c('uid','subnum','subgroup','condition','run','roi1','roi2','connection','pearson_r','z')]
  }
  all_cormat <- do.call('rbind',each_cormat)

  write.csv(all_cormat, cormat_filename, row.names = FALSE)

}

#' @export
gimmefMRI <- function(mode = 'interactive', run = 'use_config', models = 'use_config', filename = NULL, filepath = NULL){
  sinfo <- Sys.info()
  if(mode == 'example' || mode == 'demo'){
      myfile <- system.file('extdata','DemoGIMME.xlsx', package = 'gimmefMRI')
      myfile_dir <- getwd()
  } else if(mode == 'interactive') {
    if(sinfo[['sysname']] == 'Linux'){
      myfile <- my.file.browse()
    } else {
      myfile <- file.choose()
    }
    myfile <- paste(strsplit(myfile,'\\',fixed=TRUE)[[1]],collapse='/')

    myfile_dir <- paste(strsplit(myfile,'/')[[1]][1:(length(strsplit(myfile,'/')[[1]])-1)],collapse = '/')
  } else if(mode == 'filename'){
    myfile <- filename
    myfile_dir <- filepath
  }

  mm_raw <- readXLSXinput(myfile,myfile_dir)

  mm <- applyCustomSettings(mm_raw,run,models)

  mm <- runGimmeSteps(mm,mm_raw)

}

runGimmeSteps <- function(mm,mm_raw){
  runmodel_filename = as.character()
  runfigure_filename = as.character()

  if(mm$cntrl$generate_model_code){
    initializeGimmeFolders(mm)
    runmodel_filename <- writeGimmeCode(mm)
    message(sprintf('GIMME R code written to %s.',
                    file.path(mm$cntrl$model_script_directory,runmodel_filename)))
  }

  if(mm$cntrl$run_model_code){
    if(length(runmodel_filename) == 0){
      runmodel_filename <- file.path(mm$cntrl$script_save_directory, 'run_models.R')
    }
    if(file.exists(runmodel_filename)){
      message('Running GIMME model code now.')
      source(runmodel_filename, print.eval = TRUE)
      #system(sprintf('Rscript %s',runmodel_filename))
    } else{
      warning('WARNING: file %n not found. Skipped running models.', runmodel_filename)
    }

  }

  if(mm$cntrl$generate_figure_code){
    message('generate figure code not implemented yet. ')

  }
  if(mm$cntrl$run_figure_code){
    gimmefMRI_figures(mm_raw)

  }


  }

#' @export
gimmefMRIcustom <- function(
      datadir = 'C:/Users/j186b025/Documents/GitHub/jbartolotti/gimmefMRI/demodat',
      savedir = '//kumc.edu/data/Research/Hoglund/Bartolotti_J/gimme_toolkit/models',
      run_now = TRUE,
      timecourse_file = NA,
      model_spec_file = NA,
      info_lists_file = NA,
      shortnames_file = NA
      ){

  mm <- readCSVInput(datadir, timecourse_file, model_spec_file, info_lists_file, shortnames_file)
                 #timecourse = file.path(datadir,'timecourses.csv'),
                 #model_spec = file.path(datadir,'model_spec.csv'),
                 #info_lists = file.path(datadir,'info_lists.csv'),
                 #shortnames = file.path(datadir,'shortnames.csv'))

  #Creates model folders and fills input_files folders with
  #subject timecourses from mm$timecourses based on
  #model parameters in mm$model_spec
  initializeGimmeFolders(mm, savedir = savedir)

  #Creates an .R file containing code to run all models in parallel
  runmodel_filename <- writeGimmeCode(savedir, mm)
  message(sprintf('GIMME R code written to %s.',runmodel_filename))
  if(run_now){
    message('Running GIMME model code now.')
    #source(runmodel_filename)
    system(sprintf('Rscript %s',runmodel_filename))
  }

}

gimmefMRI_figures_old <- function(modeldir = '//kumc.edu/data/Research/Hoglund/Bartolotti_J/gimme_toolkit/models/BGR_g50_sg50',
                              savedir = '//kumc.edu/data/Research/Hoglund/Bartolotti_J/gimme_toolkit/models/BGR_g50_sg50/figures',
                              modelname = 'BGR_g50_sg50',

                              datadir = 'C:/Users/j186b025/Documents/GitHub/jbartolotti/gimmefMRI/demodat',
                              timecourse_file = NA,
                              model_spec_file = NA,
                              info_lists_file = NA,
                              shortnames_file = NA
                              ){

  mm <- readCSVInput(datadir, timecourse_file, model_spec_file, info_lists_file, shortnames_file)
  network_name <- mm$model_spec[mm$model_spec$model_name == modelname,'network_name']
  nodes <- mm$lists[[network_name]]
  shortnodes <- applyShorten(nodes,mm$shortnames)

  SHOW_LAG <- FALSE
  SHOW_UNCONNECTED <- TRUE
  SHOW_INDIVIDUAL <- TRUE
  ALL_COLOR <- 'black'
  INDIVIDUAL_COLOR <- 'grey'
  SUBGROUP1_COLOR <- '#009900'
  SUBGROUP2_COLOR <- '#00CCFF'
  SHAPE = 'circle' #circle, square, triangle, diamond, heart
  LAYOUT = 'circle' #circle, spring
  EDGE_CURVE <- 1
  NODE_WIDTH <- 2
  NODE_COLOR <- 'black'
  ARROW_SIZE <- 4
  FIGWIDTH <- 1000
  FIGHEIGHT <- 1000


  counts <- read.table(file.path(modeldir,'summaryPathCounts.csv'), sep = ',', stringsAsFactors = FALSE, header = TRUE)
  empty_nodes <- counts[0,]
  for(n in shortnodes){empty_nodes[dim(empty_nodes)[1]+1,] = c(n,'~',n,rep(0,sum(grepl('count',names(empty_nodes)))))}
  counts <- rbind(counts,empty_nodes)
  counts$lagged <- grepl('lag',counts$rhs)
  for(col in grep('count',names(counts))){counts[,col] <- as.numeric(counts[,col])}

  ccedge <- counts
  if(!SHOW_LAG){
    #remove rows where the predictor includes "lag"
    ccedge <- subset(counts, !lagged)
  }
  if(!SHOW_INDIVIDUAL){
    #remove the column with individual counts
    ccedge <- ccedge[,!(names(ccedge) %in% 'count.ind')]
  }
  if(!SHOW_UNCONNECTED){
    #remove rows where the sum of counts across columns is 0
    totalConn <- rowSums(ccedge[,grepl('count',names(ccedge))]) #total connection count across groups/subgroups/individuals
    ccedge <- ccedge[totalConn != 0,]
  }

  left_nodes <- unique(ccedge$lhs)
  right_nodes <- unique(ccedge$rhs)

  all_nodes <- unique(c(left_nodes,right_nodes))
  ccmat <- matrix(0,length(all_nodes),length(all_nodes))
  colnames(ccmat) <- all_nodes
  rownames(ccmat) <- all_nodes
  colormat <- matrix('#000000',length(all_nodes),length(all_nodes))
  colnames(colormat) <- all_nodes
  rownames(colormat) <- all_nodes

  for(left in left_nodes){
    for(right in right_nodes){
      thisconn <- ccedge[ccedge$lhs == left & ccedge$rhs == right,]
      thiscolor <- '#000000'
      thisweight <- 0
      if(dim(thisconn)[1] > 0){
        #individual first, so that it can get overwritten if needed. This occurs if a connection appears in one subgroup, and also in some individuals in the other subgroup.
        if('count.ind' %in% colnames(thisconn) && thisconn$count.ind > 0){
          thiscolor <- INDIVIDUAL_COLOR
          thisweight <- thisconn$count.ind
        }
        if('count.group' %in% colnames(thisconn) && thisconn$count.group > 0){
          thiscolor <- ALL_COLOR
          thisweight <- thisconn$count.group
        }
        if('count.subgroup1' %in% colnames(thisconn) && thisconn$count.subgroup1 > 0){
          thiscolor <- SUBGROUP1_COLOR
          thisweight <- thisconn$count.subgroup1
        }
        if('count.subgroup2' %in% colnames(thisconn) && thisconn$count.subgroup2 > 0){
          thiscolor <- SUBGROUP2_COLOR
          thisweight <- thisconn$count.subgroup2
        }

      }
      ccmat[left,right] <- thisweight
      colormat[left,right] <- thiscolor
    }
  }
  node_order <- all_nodes

  dir.create(savedir,showWarnings = FALSE)
  png(file.path(savedir, 'figure.png'), width = FIGWIDTH, height = FIGHEIGHT)
  qgraph::qgraph(ccmat[node_order,node_order],
         layout = LAYOUT,
         shape = SHAPE,
         border.width = NODE_WIDTH,
         border.color = NODE_COLOR,
         labels = node_order,

         edge.color = colormat[node_order,node_order],
         fade = FALSE,
         curve = EDGE_CURVE,
         asize = ARROW_SIZE
         )
  dev.off()

#  cc2 = ccedge
#  cc2$weight = 0
#  cc2$color = ''
#  cc2$weight[cc2$count.ind > 0] = cc2$count.ind[cc2$count.ind > 0]
#  cc2$color[cc2$count.ind > 0] = 'grey'
#  cc2$weight[cc2$count.group > 0] = cc2$count.group[cc2$count.group > 0]
#  cc2$color[cc2$count.group > 0] = 'black'
#  cc2$weight[cc2$count.subgroup1 > 0] = cc2$count.subgroup1[cc2$count.subgroup1 > 0]
#  cc2$color[cc2$count.subgroup1 > 0] = 'blue'
#  cc2$weight[cc2$count.subgroup2 > 0] = cc2$count.subgroup2[cc2$count.subgroup2 > 0]
#  cc2$color[cc2$count.subgroup2 > 0] = 'green'
#qgraph(cc2[,c('lhs','rhs','weight')],
#       layout = LAYOUT,
#       edge.color = cc2$color,
#       shape = SHAPE,
#       fade = FALSE
#)
}
#gimmefMRI(datadir = '~/R-Drive/Bartolotti_J/gimme_toolkit/demodat', savedir = '~/R-Drive/Bartolotti_J/gimme_toolkit/models')


#' Compute Network Metrics from GIMME Model Output
#'
#' Computes node-level and network-level graph metrics for individual networks
#' based on GIMME model output. Saves individual subject files and a combined summary.
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
#'
#' @examples
#' \dontrun{
#' # Compute metrics for a GIMME model
#' results <- computeNetworkMetrics("path/to/model_dir")
#'
#' # View network-level summary
#' head(results$network_metrics)
#'
#' # View node-level details
#' head(results$node_metrics)
#' }
#'
#' @export
computeNetworkMetrics <- function(model_dir,
                                   ignore_lags = TRUE,
                                   save_individual = TRUE,
                                   save_summary = TRUE,
                                   verbose = TRUE) {
  computeNetworkMetrics_internal(model_dir, ignore_lags, save_individual, save_summary, verbose)
}


#' Compare Network Conditions for Within-Subject Analysis
#'
#' Computes within-subject network comparison metrics between two conditions
#' (e.g., ON vs OFF medication). Supports both single-model (conditions in same folder)
#' and dual-model (conditions in separate folders) designs.
#'
#' @param model_dir_A Path to first condition model folder, or single model folder containing both conditions
#' @param model_dir_B Path to second condition model folder (NULL if using single model folder)
#' @param condition_A Condition identifier for first network (default "A")
#' @param condition_B Condition identifier for second network (default "B")
#' @param ignore_lags Logical. If TRUE (default), excludes lagged connections from analysis
#' @param save_output Logical. If TRUE (default), saves comparison files
#' @param output_dir Custom output directory path (NULL for auto-generated path in parent of model_dir_A)
#' @param verbose Logical. If TRUE (default), prints progress messages
#'
#' @return A list containing:
#'   \item{comparison_summary}{Data frame with comparison metrics for all subjects}
#'   \item{individual_comparisons}{List of detailed edge-by-edge comparisons per subject}
#'
#' @details
#' The function computes the following within-subject metrics:
#' \itemize{
#'   \item Jaccard similarity - proportion of edges shared between conditions
#'   \item Edge overlap - count and percentage of shared edges
#'   \item Edges unique to each condition
#'   \item Strength correlation - correlation of edge weights for shared connections
#'   \item Mean strength difference - average absolute difference in edge weights
#'   \item Strength distribution statistics (mean, SD) for each condition
#' }
#'
#' Output files are saved to a comparison folder with individual files and summary:
#' \itemize{
#'   \item comparison_[modelA]_vs_[modelB]/individual/[subject]_comparison.csv
#'   \item comparison_[modelA]_vs_[modelB]/comparison_summary.csv
#' }
#'
#' @examples
#' \dontrun{
#' # Dual model case (conditions in separate folders)
#' results <- compareNetworkConditions(
#'   model_dir_A = "models/4mod_AT1T2",
#'   model_dir_B = "models/4mod_BT1T2"
#' )
#'
#' # Single model case (conditions in same folder)
#' results <- compareNetworkConditions(
#'   model_dir_A = "models/combined_model",
#'   model_dir_B = NULL,
#'   condition_A = "ON",
#'   condition_B = "OFF"
#' )
#'
#' # View summary
#' head(results$comparison_summary)
#' }
#'
#' @export
compareNetworkConditions <- function(model_dir_A,
                                     model_dir_B = NULL,
                                     condition_A = "A",
                                     condition_B = "B",
                                     ignore_lags = TRUE,
                                     save_output = TRUE,
                                     output_dir = NULL,
                                     verbose = TRUE) {
  compareNetworkConditions_internal(model_dir_A, model_dir_B, condition_A, condition_B,
                                   ignore_lags, save_output, output_dir, verbose)
}


#' Plot Network Metrics Across Conditions
#'
#' Generates visualizations of network metrics showing distributions across conditions
#' and comparison metrics using beeswarm plots. Creates one figure per metric.
#'
#' @param model_dir_A Path to first condition model folder (or NULL if only plotting comparisons)
#' @param model_dir_B Path to second condition model folder (or NULL)
#' @param comparison_dir Path to comparison results folder (or NULL)
#' @param condition_A_label Label for condition A in plots (default "A")
#' @param condition_B_label Label for condition B in plots (default "B")
#' @param group_file Path to CSV file with subject_id and group columns (or NULL)
#' @param output_dir Custom directory to save figures (NULL for auto-detection)
#' @param save_figures Logical. If TRUE (default), saves figures as PNG files
#' @param verbose Logical. If TRUE (default), prints progress messages
#'
#' @return Invisible list of ggplot objects
#'
#' @details
#' This function generates publication-ready figures for network metrics:
#' \itemize{
#'   \item Network-level metrics: Creates plots comparing conditions A and B (if both provided)
#'   \item Comparison metrics: Creates plots for Jaccard similarity, edge overlap, strength correlation, etc.
#'   \item Uses beeswarm plots to show individual data points with mean and 95% CI
#' }
#'
#' If group_file is provided, plots will be stratified by group, showing condition_group
#' combinations on network metric plots and group categories on comparison metric plots.
#'
#' Network metric plots include:
#' \itemize{
#'   \item Number of edges
#'   \item Network density
#'   \item Mean/total edge strength
#'   \item Global efficiency
#'   \item Mean clustering coefficient
#'   \item Modularity
#' }
#'
#' Comparison plots include (when comparison_dir provided):
#' \itemize{
#'   \item Jaccard similarity
#'   \item Edge overlap percentage
#'   \item Strength correlation
#'   \item Mean strength difference
#'   \item Unique edges per condition
#' }
#'
#' @examples
#' \dontrun{
#' # Plot both conditions and comparisons
#' plotNetworkMetrics(
#'   model_dir_A = "models/4mod_AT1T2",
#'   model_dir_B = "models/4mod_BT1T2",
#'   comparison_dir = "models/comparison_4mod_AT1T2_vs_4mod_BT1T2",
#'   condition_A_label = "ON medication",
#'   condition_B_label = "OFF medication"
#' )
#'
#' # Plot only comparison metrics
#' plotNetworkMetrics(
#'   comparison_dir = "models/comparison_4mod_AT1T2_vs_4mod_BT1T2"
#' )
#' }
#'
#' @export
plotNetworkMetrics <- function(model_dir_A = NULL,
                               model_dir_B = NULL,
                               comparison_dir = NULL,
                               condition_A_label = "A",
                               condition_B_label = "B",
                               group_file = NULL,
                               output_dir = NULL,
                               save_figures = TRUE,
                               verbose = TRUE) {
  plotNetworkMetrics_internal(model_dir_A, model_dir_B, comparison_dir,
                              condition_A_label, condition_B_label, group_file,
                              output_dir, save_figures, verbose)
}
