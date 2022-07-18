
#' @export
getTC <- function(
    write_file = 'extract_timecourses.sh',
    config_file = 'gui',
    sub_list = 'all',
    roi_list = 'all',
    write_sh = TRUE,
    run_now = TRUE,
    write_csv = TRUE,
    filelocs = NA,
    config = NA
    ){

  if(write_sh){
    #write afni code to extract timecourses for each subject from each roi
    returndat <- getTimecourse(write_file = write_file, config_file = config_file, sub_list = sub_list, roi_list = roi_list)
  }
  if(run_now){
    #run the shell file
    system(sprintf('Rscript %s',write_file))
  }
  if(write_csv){
    #read in all those tc files and save them to the big timecourses.csv to use as input to gimmefMRI()
    genTimecoursesCSV('timecourses.csv', returndat$filelocs, returndat$config)
  }
}


#' @export
gimmefMRI <- function(){
  myfile <- file.choose()
  mm <- readXLSXinput(myfile)

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
  initializeGimmeFolders(savedir, mm)

  #Creates an .R file containing code to run all models in parallel
  runmodel_filename <- writeGimmeCode(savedir, mm)
  message(sprintf('GIMME R code written to %s.',runmodel_filename))
  if(run_now){
    message('Running GIMME model code now.')
    #source(runmodel_filename)
    system(sprintf('Rscript %s',runmodel_filename))
  }

}

gimmefMRI_figures <- function(modeldir = '//kumc.edu/data/Research/Hoglund/Bartolotti_J/gimme_toolkit/models/BGR_g50_sg50',
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
