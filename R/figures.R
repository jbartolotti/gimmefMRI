
gimmefMRI_figures <- function(mm){
  for(i in 1:dim(mm$figures)[1]){
    thisfig <- mm$figures[i,]
  }
  network_name <- mm$model_spec[mm$model_spec$model_name %in% thisfig$model_name,'network_name']
  nodes <- mm$lists[[network_name]]
  shortnodes <- applyShorten(nodes,mm$shortnames)

  SHOW_LAG <- as.logical(thisfig$show_lag_connections)
  if(is.na(SHOW_LAG)){SHOW_LAG <- FALSE}
  SHOW_UNCONNECTED <- as.logical(thisfig$show_unconnected_nodes)
  if(is.na(SHOW_UNCONNECTED)){SHOW_UNCONNECTED <- TRUE}
  SHOW_INDIVIDUAL <- as.logical(thisfig$show_individual_connections)
  if(is.na(SHOW_INDIVIDUAL)){SHOW_INDIVIDUAL <- FALSE}

  ALL_COLOR <- thisfig$color_all_groups
  INDIVIDUAL_COLOR <- thisfig$color_individuals
  SUBGROUP1_COLOR <- thisfig$color_subgroup1
  SUBGROUP2_COLOR <- thisfig$color_subgroup2
  SHAPE = thisfig$node_shape #circle, square, triangle, diamond, heart
  LAYOUT = thisfig$network_style #circle, spring
  EDGE_CURVE <- as.numeric(thisfig$edge_curvedness)
  NODE_WIDTH <- as.numeric(thisfig$node_thickness)
  NODE_COLOR <- thisfig$node_color
  ARROW_SIZE <- as.numeric(thisfig$edge_arrow_size)
  FIGWIDTH <- suppressWarnings(as.numeric(thisfig$width_pixels))
  FIGHEIGHT <- suppressWarnings(as.numeric(thisfig$height_pixels))

  if(is.na(FIGWIDTH)){ FIGWIDTH <- as.numeric(thisfig$width_inches)*as.numeric(thisfig$DPI)}
  if(is.na(FIGHEIGHT)){ FIGHEIGHT <- as.numeric(thisfig$height_inches)*as.numeric(thisfig$DPI)}

  modeldir <- file.path(mm$cntrl$model_save_directory,thisfig$model_name)

  savedir <- mm$cntrl$figure_save_directory
  savedir <- sub('{MODEL}',thisfig$model_name, savedir, fixed = TRUE)

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
