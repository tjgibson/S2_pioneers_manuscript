# function to extract peak summits from narrowPeak file ==================================================================
# peak start and end values will be replaced with the summit location
extract_summits <- function(gr, extend_width = 1L) {
  # verify input is a GRanges object
  if (!inherits(gr, "GRanges") ) {
    stop("x must be a GRanges object")
  }
  # verify that gr object is in narrowPeak format
  if ( !all(names(mcols(gr))  %in% c("name", "score", "signalValue", "pValue",  "qValue",  "peak"))) {
    stop(strwrap("GRanges object does not appear to be in narrowPeak format. Object should contain the following metadata columns: name, score, signalValue, pValue, qValue and peak"))
  }
  
  if (all(gr$peak == -1)) {
    stop("All values for 'peak' column == -1, indicating that summits were not called. Verify that your peak caller called peak summits")
  }
  
  # replace peak start and end values with summit location
  summit <- start(gr) + gr$peak - 1
  start(gr) <- summit
  end(gr) <- summit
  
  # resize peaks to desired width
  gr <- resize(gr, width = extend_width, fix = "center")
  
  # remove now meaningless peak column
  gr$peak <- NULL
  
  return(gr)
}

# function to merge multiple sets of ChIP peaks into a non-overlapping master set =========================================================
# supports multiple methods
# peaks should be a GRangesList
# method should be either "overlap" or "merge": 
# overlap method will determine overlaps among all peak sets and combine non-overlapping peaks
# merge method will combine overlapping peaks into a single peak
combine_peaks <- function(peaks, method = "overlap", min_overlap = 1L) {
  require(GenomicRanges)
  
  # check arguments
  if (length(peaks) < 2) {
    stop("Only one set of peaks provided. For merging peaks, provide multiple peak sets.")
  }
  
  if (!inherits(peaks, "GRangesList") ) {
    stop("one or more peak set is not a Granges object. peaks must be a GRangesList")
  }
  
  if (!any(method == c("overlap", "merge"))) {
    stop("Invalid method provided. method should be 'overlap' or 'merge' ")
  }
  
  if (!is_integer(min_overlap)) {
    stop("min_overlap must be an integer")
  }
  
  
  # merge peaks using "overlap" method
  if (method == "overlap") {
    combined_peaks <- peaks[[1]]
    for (i in 2:length(peaks)) {
      i_set <- peaks[[i]]
      combined_peaks <- c(combined_peaks,
                          subsetByOverlaps(i_set, combined_peaks, minoverlap = min_overlap, invert = TRUE))
    }
    
  }
  # merge peaks using "merge" method
  if (method == "merge") {
    combined_peaks <- GenomicRanges::reduce(unlist(peaks))
    
  }
  mcols(combined_peaks) <- NULL
  return(combined_peaks)
  
}

# function to build an overlap table =========================================================
# takes a Granges list object
# combines all nonoverlapping peaks from all samples, then builds a logical table indicating which peaks from each sample overlap each peak on the master list
peak_overlap_table <- function(peaks, method = "overlap", min_overlap = 1L, combine_peaks = TRUE) {
  
  # check arguments
  if (length(peaks) < 2) {
    stop("Only one set of peaks provided. For merging peaks, provide multiple peak sets.")
  }
  
  if (!inherits(peaks, "GRangesList") ) {
    stop("one or more peak set is not a Granges object. peaks must be a GRangesList")
  }
  
  if (!any(method == c("overlap", "merge"))) {
    stop("Invalid method provided. method should be 'overlap' or 'merge' ")
  }
  
  if (!is_integer(min_overlap)) {
    stop("min_overlap must be an integer")
  }
  
  
  # set names
  if(is.null(names(peaks))) {
    names(peaks) <- paste0("sample_", seq((peaks)))
  }
  
  
  if (combine_peaks) {
    # get a master set of all nonoverlapping peaks from all peak sets
    all_peaks_gr <- combine_peaks(peaks, method = method, min_overlap = min_overlap)
    
  } else {
    all_peaks_gr <- peaks[[1]]
  }
  all_peaks_df <- as.data.frame(all_peaks_gr) %>%
    dplyr::select(1:5)
  
  # build a logical overlap table indicating which peaks were detected in which samples
  for (i in 1:length(peaks)) {
    sample_name <- names(peaks)[i]
    i_peaks <- peaks[[i]]
    
    all_peaks_df[,sample_name] <- FALSE
    overlaps <- findOverlaps(all_peaks_gr, i_peaks, minoverlap = min_overlap)@from
    all_peaks_df[overlaps,sample_name] <- TRUE
  }
  return(all_peaks_df)
}

# function to check if a sequence is a palindrome ===============================
is_palindrome <- function(x) {
  require(Biostrings)
  return(ifelse(x == reverseComplement(x), TRUE, FALSE) )
}
