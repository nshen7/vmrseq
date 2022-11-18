#' Pool single-cell file together into an HDF5-based SummarizedExperiment
#' object (with or w/o sparse matrix representation).
#'
#' @description In each cell, sites with 0 < meth_read/total_read < 1 are removed.
#'
#' @param cellFiles Vector of character strings indicating single-cell file
#' paths you wish to pool into SE object(s). Cell files should be in BED-like
#' format, where the first 5 columns in each file must be: <chr>, <pos>, <strand>,
#' <meth_read>, <total_read>, in strict order. Each row contains info of one CpG.
#' Strand is assumed to be properly filpped. Cell files can be zipped ones in .gz
#' format.
#' @param writeDir A single character string indicating a folder directory where
#' you wish to store the processed SE object(s).
#' @param chrNames Single character string or a vector of character strings.
#' Only chromosomes listed in \code{selectChrs} will be processed.
#' @param sparseNAdrop Logical value indicating whether or not to use NA-dropped
#' sparseMatrix representation. 'NA-dropped' means replacing NAs as zeros and
#' then represent 0 as a very small positive value close to 0 so that the data
#' matrix is converted to a traditional sparseMatrix form with 0 as default
#' entry value. (see \code{sparseNAMatrix-class} from package \code{recommenderlab}).
#' Using NA-dropped sparseMatrix representation helps to save RAM during data
#' processing, as well as storage space saving to disk. Default value is TRUE.
#'
#' @importFrom recommenderlab dropNA
#' @importFrom data.table fread
#' @importFrom dplyr filter mutate select
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom HDF5Array saveHDF5SummarizedExperiment
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom DelayedArray rowSums

#'
#' @return
#' @export
#'
#' @examples
#'

data.pool <- function(cellFiles,
                      writeDir,
                      chrNames,
                      # sepChrs = TRUE,
                      sparseNAdrop = TRUE) {

  # TODO: making checks on input data format

  # if (!sepChrs) stop("Sorry, `sepChrs=FALSE` has not been implemented yet!")
  if (!sparseNAdrop) stop("Sorry, `sparseNAdrop=FALSE` has not been implemented yet!")

  message("Start processing chromosome-by-chromosome...")

  for (chr in chrNames) {
    message(paste0("...", chr, ": "), appendLF = FALSE)

    # Gather genomic coordinates
    pos_full <- extractCoord(cellFiles[1], chr)
    for (i in 2:length(cellFiles)) {
      pos_temp <- extractCoord(cellFiles[i], chr)
      pos_new <- pos_temp[which(!pos_temp %in% pos_full)]
      pos_full <- c(pos_full, pos_new)
      # cat(i)
    }
    pos_full <- sort(unique(pos_full))
    message("coordinates gathered; ", appendLF = FALSE)

    # Process cell info and assemble into matrix form and summarized info to a GRanges obj
    if (sparseNAdrop) {
      M_mat <- NULL
      for (i in 1:length(cellFiles)) {
        M_mat <- cbind(M_mat,
                       fillNA(cellFiles[i], chr, pos_full) %>% as.matrix() %>% recommenderlab::dropNA())
        cat(i, " ")
      }
      gr <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = pos_full, end = pos_full))
      gr$meth <- DelayedArray::rowSums(M_mat)
      gr$total <- DelayedArray::rowSums(M_mat > 0)
    } else {
      break # TODO
    }
    message("cells processed; ", appendLF = FALSE)

    # Write out processed data for current chromosome
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(M_mat = M_mat), rowRanges = gr)
    se <- se[granges(se)$total > 0, ]

    n <- nchar(writeDir); if(substring(writeDir, n, n)=='/') writeDir <- substring(writeDir,1, n-1)
    HDF5Array::saveHDF5SummarizedExperiment(
      x = se,
      dir = paste0(writeDir, "/", chr),
      replace = TRUE
    )
    message("processed file written out.")
  }
}


#' Extract genomic coordinates of a particular chromosome in a cell file.
#'
#' @param file
#' @param chr
#'
#' @return
#' @export
#'
#' @examples
extractCoord <- function(file, chr) {
  data.table::fread(cmd = paste0("zcat ", file, " | awk '$1==\"", chr, "\"' | awk '{print $2}'"))$V1
}

#' Extract and process methylation info of a particular chromosome from a cell file.
#'
#' @param file
#' @param chr
#'
#' @return
#' @export
#'
#' @examples
extractInfo <- function(file, chr) {
  df <- data.table::fread(cmd = paste0("zcat ", file, " | awk '$1==\"", chr, "\"' | awk '{print $2\"\t\"$4\"\t\"$5}'"))
  colnames(df) <- c("pos", "meth", "total")
  df <- df %>%
    dplyr::filter(meth/total == 1 | meth/total == 0) %>%
    dplyr::filter(!duplicated(pos)) %>%
    dplyr::mutate(bool = round(meth/total)) %>%
    dplyr::select(pos, bool)
  return(df)
}

#' Fill NA's in missing cites per cell
#'
#' @param file
#' @param chr
#' @param pos_full
#'
#' @return
#' @export
#'
#' @examples
fillNA <- function(file, chr, pos_full) {
  # Read in cell info
  df <- extractInfo(file, chr)

  # Fill in NA for uncovered position
  filled <- rep(NA, length(pos_full))
  ind <- which(pos_full %in% df$pos)
  if(length(ind) != nrow(df)) stop("'pos_full' does not contain all possible positions")
  filled[ind] <- df$bool

  return(filled)
}

