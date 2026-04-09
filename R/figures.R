
gimmefMRI_figures <- function(mm){
  # for compatibility with old versions. If there is no figure enable column, then enable them all.
  if(! 'enable' %in% colnames(mm$figures)){
    mm$figures$enable = TRUE
  }
  mm$figures$enable <- as.logical(mm$figures$enable)
  for(i in 1:dim(mm$figures)[1]){
    thisfig <- mm$figures[i,]
    if(thisfig$enable){

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
    }
  }
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

gimme_dotplot <- function(){
  basedir <- "P:/PHI_HL_CMH2228_NN-RCT/CM2228_NN-RCT/MR_Data/GIMME"
  modeldir <- "models/T1T2"
  filename <- 'indivPathEstimates.csv'
  dat <- read.csv(file.path(basedir, modeldir, filename))

  dat$group <- unlist(lapply(dat$file, function(x){strsplit(x, '_')[[1]][1]}))
  dat$color <- dat$group
  dat$color[dat$color == 'A'] <- 'blue'
  dat$color[dat$color == 'B'] <- 'red'
  dat$color <- factor(dat$color, levels = c('blue','red'))


  dat$pid <- unlist(lapply(dat$file, function(x){strsplit(x, '_')[[1]][2]}))
  dat$time <- unlist(lapply(dat$file, function(x){strsplit(x, '_')[[1]][3]}))

  dat$conn <- paste0(dat$lhs,dat$op,dat$rhs)
  dat$conn_clean <- dat$conn
  dat$conn_clean[dat$conn_clean == 'L_lPFC~R_lPFC'] <- 'R_lPFC -> L_lPFC'
  dat$conn_clean[dat$conn_clean == 'L_Put~L_GP'] <- 'L_GP -> L_Putamen'
  dat$conn_clean[dat$conn_clean == 'R_Put~R_GP'] <- 'R_GP -> L_Putamen'
  dat$conn_clean[dat$conn_clean == 'rACC~vmPFC'] <- 'vmPFC -> rACC'


  dat$islag <- grepl('lag',dat$rhs)
  dat$conclean_group <- paste(dat$conn_clean, dat$color, sep = ' ')

  library(Hmisc)
  ggplot2::ggplot(subset(dat, !islag & level == 'group'), ggplot2::aes(x = conclean_group, y = beta, color = color)) + ggbeeswarm::geom_quasirandom() +
    ggplot2::stat_summary(fun='mean', geom = 'point', color = 'black') +
    ggplot2::stat_summary(fun.data= Hmisc::mean_cl_normal, geom = 'pointrange', color = 'black') +
    scale_color_manual(values= c('blue','red')) +
    theme_bw() +
    labs(x = '') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave('T1T2_dot.png', width = 8, height = 6)

  library(lme4)


  dat2 <- subset(dat, conn == 'rACC~vmPFC')

  dat2$group <- factor(dat2$group)
  contrasts(dat2$group) <- c(-.5, .5)
  a <- lmer(beta ~ group +  + (1 | pid), REML = FALSE, data = dat2)
  a.base <- lmer(beta ~ + (1 | pid), REML = FALSE, data = dat2)

  dat3 <- subset(dat, conn == 'L_lPFC~R_lPFC')

  dat3$group <- factor(dat3$group)
  contrasts(dat3$group) <- c(-.5, .5)
  b <- lmer(beta ~ group +  + (1 | pid), REML = FALSE, data = dat3)
  b.base <- lmer(beta ~ + (1 | pid), REML = FALSE, data = dat3)

  print('vmPFC -> rACC')
  anova(a.base, a)
  summary(a)
  print('R_lPFC -> L_lPFC')
  anova(b.base, b)
  summary(b)



  library(tidyR)
  # Use the reshape function to reshape the dataframe
  dat2_wide <- reshape(dat2,
                     timevar = "group",
                     idvar = "pid",
                     direction = "wide")

  dat2_wide$delta <- dat2_wide$beta.B - dat2_wide$beta.A

  ggplot2::ggplot(dat2_wide, ggplot2::aes(x = 1, y = delta)) + ggbeeswarm::geom_quasirandom() + theme_bw() + ggplot2::stat_summary(fun = 'mean', geom = 'point', color = 'black') +
    ggplot2::stat_summary(fun.data = Hmisc::mean_cl_normal, geom = 'pointrange', color = 'black') + labs(x = 'vmPFC -> rACC') + coord_cartesian(x = c(0,2)) +
    ggsave('delta_dot_vmpfc_racc.png', width = 4, height = 4)


  dat3_wide <- reshape(dat3,
                       timevar = "group",
                       idvar = "pid",
                       direction = "wide")

  dat3_wide$delta <- dat3_wide$beta.B - dat3_wide$beta.A

  ggplot2::ggplot(dat3_wide, ggplot2::aes(x = 1, y = delta)) + ggbeeswarm::geom_quasirandom() + theme_bw() + ggplot2::stat_summary(fun = 'mean', geom = 'point', color = 'black') +
    ggplot2::stat_summary(fun.data = Hmisc::mean_cl_normal, geom = 'pointrange', color = 'black') + labs(x = 'R_lPFC -> L_lPFC') + coord_cartesian(x = c(0,2))
    ggsave('delta_dot_RL_lPFC.png', width = 4, height = 4)

  }


# Detect node/coordinate columns in a simple node coordinates CSV.
summaryPathCounts_detectCoordinateColumns <- function(node_coordinates) {
  lower_names <- tolower(names(node_coordinates))

  x_candidates <- c("x", "xcoord", "x_coord", "x.coordinate", "x_coordinate")
  y_candidates <- c("y", "ycoord", "y_coord", "y.coordinate", "y_coordinate")
  label_candidates <- c("node", "label", "roi", "shortroi", "name", "region", "id")

  x_index <- match(TRUE, lower_names %in% x_candidates)
  y_index <- match(TRUE, lower_names %in% y_candidates)
  label_index <- match(TRUE, lower_names %in% label_candidates)

  if (is.na(x_index) || is.na(y_index)) {
    stop("node_coordinates_file must contain x and y columns")
  }

  if (is.na(label_index)) {
    non_numeric <- which(!vapply(node_coordinates, is.numeric, logical(1)))
    if (length(non_numeric) == 0) {
      stop("node_coordinates_file must contain a node label column")
    }
    label_index <- non_numeric[1]
  }

  return(list(
    label = names(node_coordinates)[label_index],
    x = names(node_coordinates)[x_index],
    y = names(node_coordinates)[y_index]
  ))
}

summaryPathCounts_buildEdgeData <- function(counts) {
  edge_specs <- list(
    count.group = list(color = "grey20", label = "Group"),
    count.subgroup1 = list(color = "blue", label = "Subgroup 1"),
    count.subgroup2 = list(color = "red", label = "Subgroup 2")
  )

  edge_rows <- list()
  row_index <- 0
  for (i in seq_len(nrow(counts))) {
    for (count_name in names(edge_specs)) {
      if (count_name %in% names(counts) && !is.na(counts[i, count_name]) && counts[i, count_name] > 0) {
        row_index <- row_index + 1
        edge_rows[[row_index]] <- data.frame(
          from = counts$rhs[i],
          to = counts$lhs[i],
          category = edge_specs[[count_name]]$label,
          color = edge_specs[[count_name]]$color,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(edge_rows) == 0) {
    return(data.frame(from = character(), to = character(), category = character(), color = character(), stringsAsFactors = FALSE))
  }

  edges <- do.call(rbind, edge_rows)
  pair_id <- paste(edges$from, edges$to, sep = "->")
  edges$offset_index <- 0

  for (this_pair in unique(pair_id)) {
    pair_rows <- which(pair_id == this_pair)
    if (length(pair_rows) > 1) {
      edges$offset_index[pair_rows] <- seq(-(length(pair_rows) - 1) / 2,
                                           (length(pair_rows) - 1) / 2,
                                           length.out = length(pair_rows))
    }
  }

  return(edges)
}

summaryPathCounts_drawArrow <- function(from_x,
                                        from_y,
                                        to_x,
                                        to_y,
                                        offset_index,
                                        offset_step,
                                        node_radius,
                                        color,
                                        arrow_lwd,
                                        arrow_length) {
  dx <- to_x - from_x
  dy <- to_y - from_y
  edge_length <- sqrt(dx^2 + dy^2)
  if (edge_length == 0) {
    return(invisible(NULL))
  }

  ux <- dx / edge_length
  uy <- dy / edge_length
  px <- -uy
  py <- ux
  offset_x <- px * offset_index * offset_step
  offset_y <- py * offset_index * offset_step

  start_x <- from_x + offset_x + ux * node_radius
  start_y <- from_y + offset_y + uy * node_radius
  end_x <- to_x + offset_x - ux * node_radius
  end_y <- to_y + offset_y - uy * node_radius

  graphics::arrows(start_x, start_y, end_x, end_y,
                   col = color,
                   lwd = arrow_lwd,
                   length = arrow_length)
}

#' Plot directed connections from a summaryPathCounts file (internal implementation)
#'
#' Uses a finished `summaryPathCounts.csv` file and a node coordinate CSV to
#' create a directed network figure with nodes fixed at user-specified x-y
#' coordinates.
#'
#' @param summary_counts_file Path to `summaryPathCounts.csv`
#' @param node_coordinates_file Path to node coordinates CSV
#' @param output_file Optional PNG output path. If `NULL`, plots to active device.
#' @param width Width in inches for PNG output
#' @param height Height in inches for PNG output
#' @param point_cex Point size for nodes
#' @param label_cex Text size for labels
#' @param arrow_lwd Line width for arrows
#' @param arrow_length Arrowhead size
#' @param verbose Logical. If `TRUE`, prints progress messages
#'
#' @return Invisibly returns a list containing the filtered nodes and edges.
plotNetworkCoords_internal <- function(summary_counts_file,
                                       node_coordinates_file,
                                       output_file = NULL,
                                       width = 7,
                                       height = 7,
                                       point_cex = 2,
                                       label_cex = 1,
                                       arrow_lwd = 2,
                                       arrow_length = 0.12,
                                       show_axes = FALSE,
                                       verbose = TRUE) {
  if (!file.exists(summary_counts_file)) {
    stop(sprintf("summary_counts_file not found: %s", summary_counts_file))
  }
  if (!file.exists(node_coordinates_file)) {
    stop(sprintf("node_coordinates_file not found: %s", node_coordinates_file))
  }

  counts <- read.csv(summary_counts_file, stringsAsFactors = FALSE)
  node_coordinates <- read.csv(node_coordinates_file, stringsAsFactors = FALSE)

  required_cols <- c("lhs", "rhs")
  missing_cols <- setdiff(required_cols, names(counts))
  if (length(missing_cols) > 0) {
    stop(sprintf("summary_counts_file is missing required columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  coordinate_cols <- summaryPathCounts_detectCoordinateColumns(node_coordinates)
  coords <- node_coordinates[, c(coordinate_cols$label, coordinate_cols$x, coordinate_cols$y)]
  names(coords) <- c("node", "x", "y")
  coords$node <- as.character(coords$node)
  coords$x <- suppressWarnings(as.numeric(coords$x))
  coords$y <- suppressWarnings(as.numeric(coords$y))
  coords <- coords[!duplicated(coords$node), , drop = FALSE]

  if (any(is.na(coords$x)) || any(is.na(coords$y))) {
    stop("node_coordinates_file contains non-numeric x or y values")
  }

  present_count_cols <- intersect(c("count.group", "count.subgroup1", "count.subgroup2"), names(counts))
  if (length(present_count_cols) == 0) {
    stop("summary_counts_file must contain at least one of count.group, count.subgroup1, or count.subgroup2")
  }

  counts$lhs <- as.character(counts$lhs)
  counts$rhs <- as.character(counts$rhs)
  counts <- counts[!grepl("lag", counts$rhs, ignore.case = TRUE), , drop = FALSE]

  for (count_name in present_count_cols) {
    counts[[count_name]] <- suppressWarnings(as.numeric(counts[[count_name]]))
  }

  keep_rows <- rep(FALSE, nrow(counts))
  for (count_name in present_count_cols) {
    keep_rows <- keep_rows | (!is.na(counts[[count_name]]) & counts[[count_name]] > 0)
  }
  counts <- counts[keep_rows, , drop = FALSE]

  if (nrow(counts) > 0) {
    missing_nodes <- setdiff(unique(c(counts$lhs, counts$rhs)), coords$node)
    if (length(missing_nodes) > 0) {
      stop(sprintf("Missing coordinates for nodes: %s", paste(missing_nodes, collapse = ", ")))
    }
  }

  edges <- summaryPathCounts_buildEdgeData(counts)

  x_span <- diff(range(coords$x))
  y_span <- diff(range(coords$y))
  plot_span <- max(c(x_span, y_span, 1))
  node_radius <- plot_span * 0.03
  offset_step <- plot_span * 0.03

  # Use the same half-range on both axes so 1 unit = 1 unit in both directions
  axis_pad <- plot_span * 0.15
  half_range <- plot_span / 2 + axis_pad
  center_x <- mean(range(coords$x))
  center_y <- mean(range(coords$y))

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(output_file, width = width, height = height, units = "in", res = 300)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  if (!show_axes) {
    graphics::par(mar = c(0.2, 0.2, 0.2, 0.2))
  }

  graphics::plot(coords$x, coords$y,
                 type = "n",
                 axes = show_axes,
                 xlab = if (show_axes) "x" else "",
                 ylab = if (show_axes) "y" else "",
                 asp = 1,
                 xlim = c(center_x - half_range, center_x + half_range),
                 ylim = c(center_y - half_range, center_y + half_range),
                 bty = if (show_axes) "o" else "n")

  graphics::points(coords$x, coords$y, pch = 16, cex = point_cex)

  # Label to the right by default; to the left for nodes near the leftmost edge
  x_range <- diff(range(coords$x))
  left_threshold <- min(coords$x) + x_range * 0.2
  label_pos <- ifelse(coords$x <= left_threshold, 2L, 4L)
  graphics::text(coords$x, coords$y, labels = coords$node, pos = label_pos, cex = label_cex)

  if (nrow(edges) > 0) {
    for (i in seq_len(nrow(edges))) {
      from_row <- coords[coords$node == edges$from[i], , drop = FALSE]
      to_row <- coords[coords$node == edges$to[i], , drop = FALSE]
      summaryPathCounts_drawArrow(from_row$x, from_row$y,
                                  to_row$x, to_row$y,
                                  edges$offset_index[i],
                                  offset_step,
                                  node_radius,
                                  edges$color[i],
                                  arrow_lwd,
                                  arrow_length)
    }

    present_categories <- unique(edges[, c("category", "color")])
    graphics::legend("topright",
                     legend = present_categories$category,
                     col = present_categories$color,
                     lwd = arrow_lwd,
                     bty = "n")
  } else {
    warning("No non-lag paths with positive group or subgroup counts were found")
  }

  if (verbose) {
    message(sprintf("Plotted %d nodes and %d arrows", nrow(coords), nrow(edges)))
    if (!is.null(output_file)) {
      message(sprintf("Saved figure to: %s", output_file))
    }
  }

  return(invisible(list(nodes = coords, edges = edges, filtered_counts = counts)))
}

#' Plot Network Metrics Across Conditions (Internal Implementation)
#'
#' Generates ggplot2 visualizations of network metrics, showing distributions
#' across conditions and comparison metrics using beeswarm plots.
#'
#' @param model_dir_A Path to first condition model folder (or NULL)
#' @param model_dir_B Path to second condition model folder (or NULL)
#' @param comparison_dir Path to comparison results folder (or NULL)
#' @param condition_A_label Label for condition A (default "A")
#' @param condition_B_label Label for condition B (default "B")
#' @param group_file Path to CSV with subject_id and group columns (or NULL)
#' @param output_dir Directory to save figures (NULL for auto-detection)
#' @param save_figures Logical. If TRUE (default), saves figures to output_dir
#' @param verbose Logical. If TRUE (default), prints progress messages
#'
#' @return List of ggplot objects
plotNetworkMetrics_internal <- function(model_dir_A = NULL,
                                        model_dir_B = NULL,
                                        comparison_dir = NULL,
                                        condition_A_label = "A",
                                        condition_B_label = "B",
                                        group_file = NULL,
                                        output_dir = NULL,
                                        save_figures = TRUE,
                                        verbose = TRUE) {
  
  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
    stop("Package 'ggbeeswarm' is required. Install with: install.packages('ggbeeswarm')")
  }
  
  # Load group assignments if provided
  group_data <- NULL
  if (!is.null(group_file)) {
    if (!file.exists(group_file)) {
      warning(sprintf("Group file not found: %s", group_file))
    } else {
      group_data <- read.csv(group_file, stringsAsFactors = FALSE)
      if (verbose) message(sprintf("Loaded group assignments for %d subjects", nrow(group_data)))
    }
  }
  
  # Determine output directory
  if (is.null(output_dir)) {
    if (!is.null(comparison_dir)) {
      output_dir <- file.path(comparison_dir, "figures")
    } else if (!is.null(model_dir_A)) {
      output_dir <- file.path(dirname(model_dir_A), "network_metric_figures")
    } else if (!is.null(model_dir_B)) {
      output_dir <- file.path(dirname(model_dir_B), "network_metric_figures")
    } else {
      stop("At least one of model_dir_A, model_dir_B, or comparison_dir must be provided")
    }
  }
  
  if (save_figures) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (verbose) message(sprintf("Saving figures to: %s", output_dir))
  }
  
  plot_list <- list()
  
  # Load network metrics data from both conditions
  network_data <- NULL
  
  if (!is.null(model_dir_A)) {
    metrics_file_A <- file.path(model_dir_A, "networkMetrics_summary.csv")
    if (file.exists(metrics_file_A)) {
      data_A <- read.csv(metrics_file_A, stringsAsFactors = FALSE)
      data_A$condition_label <- condition_A_label
      network_data <- data_A
      if (verbose) message(sprintf("Loaded %d subjects from condition A", nrow(data_A)))
    } else {
      warning(sprintf("File not found: %s", metrics_file_A))
    }
  }
  
  if (!is.null(model_dir_B)) {
    metrics_file_B <- file.path(model_dir_B, "networkMetrics_summary.csv")
    if (file.exists(metrics_file_B)) {
      data_B <- read.csv(metrics_file_B, stringsAsFactors = FALSE)
      data_B$condition_label <- condition_B_label
      if (is.null(network_data)) {
        network_data <- data_B
      } else {
        network_data <- rbind(network_data, data_B)
      }
      if (verbose) message(sprintf("Loaded %d subjects from condition B", nrow(data_B)))
    } else {
      warning(sprintf("File not found: %s", metrics_file_B))
    }
  }
  
  # Merge with group data if provided
  if (!is.null(network_data) && !is.null(group_data)) {
    # Extract subject ID from 'subject' column
    network_data$subject_id <- as.numeric(network_data$subject)
    group_data$subject_id <- as.numeric(group_data$subject_id)
    network_data <- merge(network_data, group_data[, c("subject_id", "group")], 
                         by = "subject_id", all.x = TRUE)
    
    # Create combined condition_group factor for x-axis
    network_data$condition_group <- paste(network_data$condition_label, 
                                          network_data$group, sep = "_")
    if (verbose) message("Merged group assignments with network data")
  }
  
  # Generate plots for network-level metrics
  if (!is.null(network_data)) {
    network_metrics <- c("n_edges", "density", "mean_strength", "total_strength",
                        "global_efficiency", "mean_clustering", "modularity")
    
    metric_labels <- c(
      n_edges = "Number of Edges",
      density = "Network Density",
      mean_strength = "Mean Edge Strength",
      total_strength = "Total Network Strength",
      global_efficiency = "Global Efficiency",
      mean_clustering = "Mean Clustering Coefficient",
      modularity = "Modularity"
    )
    
    # Determine x-axis variable based on whether groups are available
    x_var <- if (!is.null(group_data) && "condition_group" %in% colnames(network_data)) {
      "condition_group"
    } else {
      "condition_label"
    }
    x_label <- if (x_var == "condition_group") "Condition & Group" else "Condition"
    
    for (metric in network_metrics) {
      if (metric %in% colnames(network_data)) {
        if (verbose) message(sprintf("Generating plot for: %s", metric))
        
        p <- ggplot2::ggplot(network_data, ggplot2::aes(x = .data[[x_var]], y = .data[[metric]])) +
          ggbeeswarm::geom_quasirandom(size = 2, alpha = 0.6) +
          ggplot2::stat_summary(fun.data = Hmisc::mean_cl_normal, geom = "pointrange",
                                color = "red", size = 0.8, linewidth = 1) +
          ggplot2::labs(
            title = metric_labels[metric],
            x = x_label,
            y = metric_labels[metric]
          ) +
          ggplot2::theme_classic(base_size = 12) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            axis.text = ggplot2::element_text(color = "black"),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
          )
        
        plot_list[[metric]] <- p
        
        if (save_figures) {
          filename <- file.path(output_dir, sprintf("network_%s.png", metric))
          ggplot2::ggsave(filename, p, width = 6, height = 5, dpi = 300)
        }
      }
    }
  }
  
  # Load and plot comparison metrics
  if (!is.null(comparison_dir)) {
    comparison_file <- file.path(comparison_dir, "comparison_summary.csv")
    if (file.exists(comparison_file)) {
      comp_data <- read.csv(comparison_file, stringsAsFactors = FALSE)
      if (verbose) message(sprintf("Loaded comparison data for %d subjects", nrow(comp_data)))
      
      # Filter to subjects with both conditions
      comp_data_both <- comp_data[comp_data$has_condition_A & comp_data$has_condition_B, ]
      
      if (nrow(comp_data_both) > 0) {
        # Merge with group data if provided
        if (!is.null(group_data)) {
          comp_data_both$subject_id <- as.numeric(comp_data_both$subject)
          group_data$subject_id <- as.numeric(group_data$subject_id)
          comp_data_both <- merge(comp_data_both, group_data[, c("subject_id", "group")], 
                                 by = "subject_id", all.x = TRUE)
        }
        
        comparison_metrics <- c("jaccard_similarity", "edge_overlap_pct",
                               "strength_correlation", "mean_strength_diff")
        
        comparison_labels <- c(
          jaccard_similarity = "Jaccard Similarity",
          edge_overlap_pct = "Edge Overlap (%)",
          strength_correlation = "Strength Correlation",
          mean_strength_diff = "Mean Strength Difference"
        )
        
        # Determine x-axis variable based on whether groups are available
        x_var_comp <- if (!is.null(group_data) && "group" %in% colnames(comp_data_both)) {
          "group"
        } else {
          NA
        }
        
        for (metric in comparison_metrics) {
          if (metric %in% colnames(comp_data_both)) {
            if (verbose) message(sprintf("Generating plot for: %s", metric))
            
            if (!is.na(x_var_comp)) {
              # Plot with groups on x-axis
              p <- ggplot2::ggplot(comp_data_both, ggplot2::aes(x = .data[[x_var_comp]], y = .data[[metric]])) +
                ggbeeswarm::geom_quasirandom(size = 2, alpha = 0.6) +
                ggplot2::stat_summary(fun.data = Hmisc::mean_cl_normal, geom = "pointrange",
                                      color = "red", size = 0.8, linewidth = 1) +
                ggplot2::labs(
                  title = comparison_labels[metric],
                  x = "Group",
                  y = comparison_labels[metric]
                ) +
                ggplot2::theme_classic(base_size = 12) +
                ggplot2::theme(
                  plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                  axis.text = ggplot2::element_text(color = "black")
                )
            } else {
              # Plot without groups (all subjects together)
              p <- ggplot2::ggplot(comp_data_both, ggplot2::aes(x = 1, y = .data[[metric]])) +
                ggbeeswarm::geom_quasirandom(size = 2, alpha = 0.6) +
                ggplot2::stat_summary(fun.data = Hmisc::mean_cl_normal, geom = "pointrange",
                                      color = "red", size = 0.8, linewidth = 1) +
                ggplot2::labs(
                  title = comparison_labels[metric],
                  x = sprintf("%s vs %s", condition_A_label, condition_B_label),
                  y = comparison_labels[metric]
                ) +
                ggplot2::theme_classic(base_size = 12) +
                ggplot2::theme(
                  plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                  axis.text = ggplot2::element_text(color = "black"),
                  axis.text.x = ggplot2::element_blank(),
                  axis.ticks.x = ggplot2::element_blank()
                )
            }
            
            plot_list[[paste0("comparison_", metric)]] <- p
            
            if (save_figures) {
              filename <- file.path(output_dir, sprintf("comparison_%s.png", metric))
              ggplot2::ggsave(filename, p, width = 5, height = 5, dpi = 300)
            }
          }
        }
        
        # Additional comparison plot: edges only in A vs only in B
        if ("edges_only_A" %in% colnames(comp_data_both) && 
            "edges_only_B" %in% colnames(comp_data_both)) {
          
          # Reshape data for plotting, incorporating groups if available
          if (!is.na(x_var_comp) && "group" %in% colnames(comp_data_both)) {
            # Create condition_group combinations for x-axis
            edge_diff_data <- data.frame(
              subject = rep(comp_data_both$subject, 2),
              condition_group = rep(NA, nrow(comp_data_both) * 2),
              n_edges = c(comp_data_both$edges_only_A, comp_data_both$edges_only_B)
            )
            
            # Fill condition_group for first half (Only A)
            edge_diff_data$condition_group[1:nrow(comp_data_both)] <- 
              paste(sprintf("Only %s", condition_A_label), 
                    comp_data_both$group, sep = "_")
            
            # Fill condition_group for second half (Only B)
            edge_diff_data$condition_group[(nrow(comp_data_both) + 1):nrow(edge_diff_data)] <- 
              paste(sprintf("Only %s", condition_B_label), 
                    comp_data_both$group, sep = "_")
            
            x_axis_var <- "condition_group"
            x_axis_label <- sprintf("%s & Group", "Condition")
          } else {
            # Without groups, just condition
            edge_diff_data <- data.frame(
              subject = rep(comp_data_both$subject, 2),
              condition = rep(c(sprintf("Only %s", condition_A_label), 
                               sprintf("Only %s", condition_B_label)), 
                             each = nrow(comp_data_both)),
              n_edges = c(comp_data_both$edges_only_A, comp_data_both$edges_only_B)
            )
            x_axis_var <- "condition"
            x_axis_label <- "Condition"
          }
          
          p <- ggplot2::ggplot(edge_diff_data, ggplot2::aes(x = .data[[x_axis_var]], y = n_edges)) +
            ggbeeswarm::geom_quasirandom(size = 2, alpha = 0.6) +
            ggplot2::stat_summary(fun.data = Hmisc::mean_cl_normal, geom = "pointrange",
                                  color = "red", size = 0.8, linewidth = 1) +
            ggplot2::labs(
              title = "Unique Edges per Condition",
              x = x_axis_label,
              y = "Number of Unique Edges"
            ) +
            ggplot2::theme_classic(base_size = 12) +
            ggplot2::theme(
              plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
              axis.text = ggplot2::element_text(color = "black"),
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
            )
          
          plot_list[["comparison_unique_edges"]] <- p
          
          if (save_figures) {
            filename <- file.path(output_dir, "comparison_unique_edges.png")
            ggplot2::ggsave(filename, p, width = 6, height = 5, dpi = 300)
          }
        }
      } else {
        warning("No subjects with both conditions found in comparison data")
      }
    } else {
      warning(sprintf("Comparison file not found: %s", comparison_file))
    }
  }
  
  if (verbose) message(sprintf("Generated %d plots", length(plot_list)))
  
  return(invisible(plot_list))
}
