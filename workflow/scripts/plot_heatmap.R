# functions for plotting deeptools-style heatmaps ================================================
# function for quickly computing a coverage matrix from a bigwig file. Core code taken from seqPlots package
coverage_matrix <- function(file,regions, upstream = 1000, downstream = 1000, region_width = 1, bin_size =10) {
  require(rtracklayer)
  require(GenomicRanges)
  
  # check arguments
  if (!is.character(file)) stop("file should be a character vector")
  
  if (!file.exists(file)) stop("one or more input files not found")
  
  # get regions to include in heatmap
  gr <- GenomicRanges::resize(regions, region_width, fix='center')
  gr <- suppressWarnings( GenomicRanges::promoters(gr, upstream, downstream ) )
  gr <- trim(gr)
  hm_bins  <- seq(-upstream, downstream, by=bin_size )
  
  # extract_matrix function taken from seqPlots package
  extract_matrix <- function(track, gr, size, ignore_strand) {
    sum <- suppressWarnings(.Call(
      'BWGFile_summary', path.expand(path(track)),
      as.character(seqnames(gr)), ranges(gr), 
      S4Vectors::recycleIntegerArg(size, "size", length(gr)), 
      "mean", as.numeric(NA_real_), PACKAGE='rtracklayer'
    ))
    M <- do.call( rbind, sum )
    if (!ignore_strand) 
      M[as.character(strand(gr))=='-',] <- M[
        as.character(strand(gr))=='-', ncol(M):1]
    return(M)
  }
  
  # read bigwig file to matrix
  bwf <- BigWigFile(file)
  coverage_matrix <- extract_matrix(bwf, gr, length(hm_bins), ignore_strand = TRUE)
  
  return(coverage_matrix)
}


plot_heatmap <- function(bw, regions, 
                         upstream = 1000, downstream = 1000, region_width = 1,
                         row_split = NULL,
                         order_by_samples = NULL ,
                         colors = c("#440154FF", "#21908CFF", "#FDE725FF"),
                         scale_group = NULL,
                         individual_scales = FALSE, 
                         return_heatmap_list = FALSE,
                         min_percentile = 0.01,
                         max_percentile = 0.98,
                         row_order = NULL,
                         plot_average = TRUE,
                         average_gp = gpar(lwd = 3, col = n_groups),
                         ...) {
  require(rtracklayer)
  require(GenomicRanges)
  require(tidyverse)
  require(EnrichedHeatmap)
  require(RColorBrewer)
  
  # check arguments
  if (!is.character(bw)) stop("bw should be a character vector")
  
  if (!all(file.exists(bw))) stop("one or more input files not found")
  
  if (any(!str_detect(bw, ".bw$"))) stop("bw file extension missing from one or more files, check input file formats")
  
  
  if (!is.null(scale_group)) {
    if (!is.numeric(scale_group)) stop("scale_group should be a numeric vector")
    if (length(scale_group) != 1 & length(scale_group) != length(bw)) stop("scale_group should be a single number or a vector of numbers of the same length as bw")
  }
  
  # get sample names
  if (is.null(names(bw))) {
    names(bw) <- str_replace(basename(bw), ".bw", "")
  }
  
  # read in coverage matrices
  message("reading coverage matrix from bigWig files")
  mat_list <-bw %>% 
    map(coverage_matrix, regions = regions, upstream = upstream, downstream = downstream, region_width = region_width) 
  
  # check for missing values 
  any_NA <- mat_list %>%
    map(is.na) %>%
    map(any) %>%
    unlist()
  
  if (any(any_NA)) {
    warning("NAs detected in coverage matrices. Replacing NAs with 0")
    
    replace_mat_na <- function(mat) {
      mat[is.na(mat)] <- 0
      return(mat)
    }
    
    mat_list <- mat_list %>%
      map(replace_mat_na)
  }
  
  
  # get row order
  if (is.null(row_order)) {
    if(is.null(order_by_samples)) {
      order <- mat_list %>% 
        map(as.data.frame) %>% 
        bind_cols %>% 
        rowMeans() %>% 
        order(decreasing = TRUE)
    } else if (length(order_by_samples) == 1) {
      order <- mat_list[[order_by_samples]] %>% 
        rowMeans() %>% 
        order(decreasing = TRUE)
    } else {
      order <- mat_list[order_by_samples] %>% 
        map(as.data.frame) %>%
        bind_cols() %>%
        rowMeans() %>% 
        order(decreasing = TRUE)
    }
  } else {
    order  <- row_order
  }
  
  # # apply gaussian smoothing
  # if (smooth_heatmap) {
  # message("applying gaussian smoothing to heatmap")
  # mat_list <- mat_list %>%
  #   map(gaussianFilter)
  # }
  #   
  message("converting matrix to enrichedHeatmap normalizedMatrix object")
  mat_list <- mat_list %>%
    map(as.normalizedMatrix, k_upstream = upstream/10, k_downstream = downstream/10, k_target = region_width, extend = c(upstream,downstream),
        smooth = FALSE,
        keep = c(min_percentile, max_percentile))
  
  # get color scale for heatmaps
  if (!individual_scales) {
    col_min <- mat_list %>%
      map(na.exclude) %>%
      map(min) %>%
      unlist() %>%
      min()
    
    col_max <- mat_list %>%
      map(na.exclude) %>%
      map(max) %>%
      unlist() %>%
      max()
    
    col_fun <- circlize::colorRamp2(seq(col_min, col_max, length.out = length(colors)), colors)
    
  } 
  
  # construct heatmap_list object
  hm_list <- NULL
  
  # set legend
  if (!is.null(row_split)) {
    n_groups <- seq_along(unique(row_split))
    lgd = Legend(at = sort(unique(row_split)), title = "Group", 
                 type = "lines", legend_gp = average_gp)
  } else {
    n_groups <- 1
    lgd <- NULL
  }
  
  for (i in seq_along(mat_list)) {
    # in individual_scales and/or scale_group is set, determine separate scale for individual heatmaps or groups of heatmaps
    if (individual_scales) {
      col_fun <- circlize::colorRamp2(seq(min(na.exclude(mat_list[[i]])), max(na.exclude(mat_list[[i]])), length.out = length(colors)), colors)
    } 
    
    if (!is.null(scale_group)) {
      if (individual_scales) {
        warning("scale_groups and individual_scales are both set. scale_groups will override individual_scales for determining color scales")
      }
      
      if (length(scale_group == length(mat_list))) {
        col_list <- list()
        
        group_i <- which(scale_group == scale_group[i])
        
        col_min <- mat_list[group_i] %>%
          map(na.exclude) %>%
          map(min) %>%
          unlist() %>%
          min()
        
        col_max <- mat_list[group_i] %>%
          map(na.exclude) %>%
          map(max) %>%
          unlist() %>%
          max()
        
        
        
        col_fun <-  circlize::colorRamp2( seq(col_min, col_max, length.out = length(colors)), colors)
        
        
      } else {
        col_fun <- circlize::colorRamp2(seq(min(na.exclude(mat_list[[scale_group]])), max(na.exclude(mat_list[[scale_group]])), length.out = length(colors)), colors)
      }
      
    } 
    
    # get range for metaplot
    hm_range <- range(attr(col_fun, "breaks"))
    metaplot_range <- c(hm_range[1], ceiling(hm_range[2] * 0.75))
    
    # add heatmap to list
    message("creating heatmap for sample ", i)
    hm_list <- hm_list + EnrichedHeatmap(mat_list[[i]],  col = col_fun,
                                         row_split = row_split,
                                         row_title_rot = 0,
                                         pos_line = FALSE, 
                                         row_order = order,
                                         column_title = names(bw)[i], 
                                         name = names(bw)[i],
                                         top_annotation = if (plot_average) HeatmapAnnotation(lines = anno_enriched(ylim = metaplot_range, gp = average_gp))  else {NULL},
                                         ...)
    
  }
  
  if (!return_heatmap_list) {
    if (is.null(lgd)) {
      draw(hm_list, merge_legends = TRUE)
    } else {
      draw(hm_list,annotation_legend_list = list(lgd), merge_legend = TRUE, heatmap_legend_side = "bottom", 
           annotation_legend_side = "bottom")
    }
  } else {
    return(hm_list)
  }
  
}

plot_heatmap_minimal <- function(bw, regions, 
                                 upstream = 1000, downstream = 1000, region_width = 1,
                                 row_split = NULL,
                                 order_by_samples = NULL ,
                                 colors = c("#440154FF", "#21908CFF", "#FDE725FF"),
                                 scale_group = NULL,
                                 individual_scales = FALSE, 
                                 return_heatmap_list = FALSE,
                                 min_percentile = 0.01,
                                 max_percentile = 0.98,
                                 row_order = NULL,
                                 ...) {
  require(rtracklayer)
  require(GenomicRanges)
  require(tidyverse)
  require(ComplexHeatmap)
  require(RColorBrewer)
  
  # check arguments
  if (!is.character(bw)) stop("bw should be a character vector")
  
  if (!all(file.exists(bw))) stop("one or more input files not found")
  
  if (any(!str_detect(bw, ".bw$"))) stop("bw file extension missing from one or more files, check input file formats")
  
  
  if (!is.null(scale_group)) {
    if (!is.numeric(scale_group)) stop("scale_group should be a numeric vector")
    if (length(scale_group) != 1 & length(scale_group) != length(bw)) stop("scale_group should be a single number or a vector of numbers of the same length as bw")
  }
  
  # get sample names
  if (is.null(names(bw))) {
    names(bw) <- str_replace(basename(bw), ".bw", "")
  }
  
  # read in coverage matrices
  message("reading coverage matrix from bigWig files")
  mat_list <-bw %>% 
    map(coverage_matrix, regions = regions, upstream = upstream, downstream = downstream, region_width = region_width) 
  
  # check for missing values 
  any_NA <- mat_list %>%
    map(is.na) %>%
    map(any) %>%
    unlist()
  
  if (any(any_NA)) {
    warning("NAs detected in coverage matrices. Replacing NAs with 0")
    
    replace_mat_na <- function(mat) {
      mat[is.na(mat)] <- 0
      return(mat)
    }
    
    mat_list <- mat_list %>%
      map(replace_mat_na)
  }
  
  
  # get row order
  if (is.null(row_order)) {
    if(is.null(order_by_samples)) {
      order <- mat_list %>% 
        map(as.data.frame) %>% 
        bind_cols %>% 
        rowMeans() %>% 
        base::order(decreasing = TRUE)
    } else if (length(order_by_samples) == 1) {
      order <- mat_list[[order_by_samples]] %>% 
        rowMeans() %>% 
        base::order(decreasing = TRUE)
    } else {
      order <- mat_list[order_by_samples] %>% 
        map(as.data.frame) %>%
        bind_cols() %>%
        rowMeans() %>% 
        base::order(decreasing = TRUE)
    }
  } else {
    order  <- row_order
  }
  
  # replace outlier values based on quantiles
  replace_outliers <- function(mat, min, max) {
    q1 <-  quantile(mat, min, na.rm = TRUE)
    q2 <-  quantile(mat, max, na.rm = TRUE)
    mat[mat <= q1] <-  q1
    mat[mat >= q2] <-  q2
    return(mat)
  }
  
  mat_list <- mat_list %>%
    map(replace_outliers, min = min_percentile, max = max_percentile)
  
  
  
  # get color scale for heatmaps
  if (!individual_scales) {
    col_min <- mat_list %>%
      map(na.exclude) %>%
      map(min) %>%
      unlist() %>%
      min()
    
    col_max <- mat_list %>%
      map(na.exclude) %>%
      map(max) %>%
      unlist() %>%
      max()
    
    col_fun <- circlize::colorRamp2(seq(col_min, col_max, length.out = length(colors)), colors)
    
  } 
  
  # construct heatmap_list object
  hm_list <- NULL
  
  for (i in seq_along(mat_list)) {
    # in individual_scales and/or scale_group is set, determine separate scale for individual heatmaps or groups of heatmaps
    if (individual_scales) {
      col_fun <- circlize::colorRamp2(seq(min(na.exclude(mat_list[[i]])), max(na.exclude(mat_list[[i]])), length.out = length(colors)), colors)
    } 
    
    if (!is.null(scale_group)) {
      if (individual_scales) {
        warning("scale_groups and individual_scales are both set. scale_groups will override individual_scales for determining color scales")
      }
      
      if (length(scale_group == length(mat_list))) {
        col_list <- list()
        
        group_i <- which(scale_group == scale_group[i])
        
        col_min <- mat_list[group_i] %>%
          map(na.exclude) %>%
          map(min) %>%
          unlist() %>%
          min()
        
        col_max <- mat_list[group_i] %>%
          map(na.exclude) %>%
          map(max) %>%
          unlist() %>%
          max()
        
        
        
        col_fun <-  circlize::colorRamp2( seq(col_min, col_max, length.out = length(colors)), colors)
        
        
      } else {
        col_fun <- circlize::colorRamp2(seq(min(na.exclude(mat_list[[scale_group]])), max(na.exclude(mat_list[[scale_group]])), length.out = length(colors)), colors)
      }
      
    } 
    
    
    # add heatmap to list
    
    
    message("creating heatmap for sample ", i)
    hm_list <- hm_list + ComplexHeatmap::Heatmap(mat_list[[i]],  col = col_fun,
                                                 row_split = row_split,
                                                 row_order = order,
                                                 cluster_rows = FALSE, 
                                                 cluster_columns = FALSE,
                                                 row_title = NULL,
                                                 ...)
    
  }
  if (!return_heatmap_list) {
    
    draw(hm_list)
  } else {
    return(hm_list)
  }
  
}

# function for generating metaplots/average plots ==============================
plot_average <- function(
    bw, 
    regions, 
    upstream = 1000, downstream = 1000, region_width = 1,
    row_split = NULL,
    colors = NULL,
    summary_fun = mean,
    line_width = 1,
    center_label = "center"
) {
  
  require(rtracklayer)
  require(GenomicRanges)
  require(tidyverse)
  
  # check arguments
  if (!is.character(bw)) stop("bw should be a character vector")
  
  if (!all(file.exists(bw))) stop("input file not found")
  
  if (any(!str_detect(bw, ".bw$"))) stop("bw file extension missing from input file, check input file format")
  
  # get sample names
  if (is.null(names(bw))) {
    names(bw) <- str_replace(basename(bw), ".bw", "")
  }
  
  # read in coverage matrices
  message("reading coverage matrix from bigWig file")
  mat_list <-bw %>% 
    map(coverage_matrix, regions = regions, upstream = upstream, downstream = downstream, region_width = region_width) 
  
  # check for missing values 
  any_NA <- mat_list %>%
    map(is.na) %>%
    map(any) %>%
    unlist()
  
  if (any(any_NA)) {
    warning("NAs detected in coverage matrices. Replacing NAs with 0")
    
    replace_mat_na <- function(mat) {
      mat[is.na(mat)] <- 0
      return(mat)
    }
    
    mat_list <- mat_list %>%
      map(replace_mat_na)
  }
  
  # convert coverage matrix to tibble
  coverage_tb <- mat_list[[1]] %>% 
    as.data.frame() %>% 
    tibble()
  
  if (!is.null(row_split)) {
  # add column with row groupings
  coverage_tb$group <- row_split
  
  # compute average signal within groups
  coverage_tb <- coverage_tb %>% 
    group_by(group) %>% 
    summarise(across(.fns = summary_fun))
  
  # set column names
  colnames(coverage_tb)[2:ncol(coverage_tb)] <- 1:(ncol(coverage_tb) - 1)
  
  
  # reformat tibble for ggplot
  coverage_tb <- coverage_tb %>% 
    pivot_longer(2:ncol(coverage_tb), names_to = "position", values_to = "signal") %>% 
    mutate(position = as.integer(position))
  
  # initialise plot
  p <- coverage_tb %>% 
    ggplot(aes(x = position, y = signal, color = group)) + geom_line(lwd = line_width * 2)
  
  } else {
    # compute average signal within groups
    coverage_tb <- coverage_tb %>% 
      summarise(across(.fns = summary_fun))
    
    # set column names
    colnames(coverage_tb) <- 1:(ncol(coverage_tb))
    
    
    # reformat tibble for ggplot
    coverage_tb <- coverage_tb %>% 
      pivot_longer(1:ncol(coverage_tb), names_to = "position", values_to = "signal") %>% 
      mutate(position = as.integer(position))
    
    # initialise plot
    p <- coverage_tb %>% 
      ggplot(aes(x = position, y = signal)) + geom_line(lwd = line_width * 2)
  }
  # get position start and stop for axis labels
  start <- min(coverage_tb$position)
  end <- max(coverage_tb$position)
  
  # generate plot
  p <- p +
    coord_cartesian(clip = "off") +
    # viridis::scale_color_viridis(discrete=TRUE) +
    theme(text = element_text(size = 16),
          panel.border = element_rect(color = "black", fill=NA, size=line_width),
          axis.ticks = element_line(colour = "black", size = line_width),
          axis.ticks.length=unit(line_width / 5, "cm"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank()
    ) + 
    scale_x_continuous(
      
      breaks = c(start,end/2,end), 
      labels = c(paste0("-", upstream / 1000, "KB"),center_label, paste0("+", upstream / 1000, "KB")),
      expand = c(0,0,0,0)
      ) 
  
  # return plot
  return(p)
  
}

