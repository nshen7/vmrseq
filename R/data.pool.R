#' Pool single-cell file together into an HDF5-based SummarizedExperiment
#' object with sparse matrix representation.
#'
#' @description This function pools individual-cell CpG read files into a
#' SummarizedExperiment object so that it can be input to the \code{vmrseq.smooth}
#' function. Note that in each cell, sites with hemimethylation or intermediate
#' methylation levels (i.e., 0 < meth_read/total_read < 1) will be removed.
#'
#' @param cellFiles Vector of character strings indicating single-cell file
#' paths you wish to pool into SE object(s). Cell files should be in BED-like
#' format, where the first 5 columns in each file must be: <chr>, <pos>, <strand>,
#' <meth_read>, <total_read>, in strict order. Each row contains info of one CpG.
#' Strand is assumed to be properly filpped. Cell files can be zipped ones in .gz
#' format.
#' @param sep the field separator character. Values on each line of cell file are
#' separated by this character.
#' @param writeDir A single character string indicating a folder directory where
#' you wish to store the processed SE object(s). The SE will be stored in HDF5 format.
#' @param chrNames Single or a vector of character strings representing chromosome
#' names in cell files. Only chromosomes listed in \code{selectChrs} will be
#' processed.
#' @param colData A DataFrame or data.frame object containing colData for the SE
#' object (applied to all chromosomes).
#' @param sparseNAdrop Logical value indicating whether or not to use NA-dropped
#' sparseMatrix representation. 'NA-dropped' means replacing NAs as zeros and
#' then represent 0 as a very small positive value close to 0 so that the data
#' matrix is converted to a traditional sparseMatrix form with 0 as default
#' entry value. (see \code{sparseNAMatrix-class} from package \code{recommenderlab}).
#' Using NA-dropped sparseMatrix representation helps to save RAM during data
#' processing, as well as storage space saving to disk. Currently only supports
#' TRUE!!
#'
#' @importFrom recommenderlab dropNA
#' @importFrom data.table fread
#' @importFrom dplyr filter mutate select
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom HDF5Array saveHDF5SummarizedExperiment
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom DelayedArray rowSums
#' @importFrom S4Vectors DataFrame

#'
#' @return Directly write out to the `writeDir` and does not return anything.
#' @export
#'

data.pool <- function(cellFiles,
                      sep,
                      writeDir,
                      chrNames,
                      colData = NULL,
                      sparseNAdrop = TRUE) {

  # TODO: making checks on input data format
  chrNames <- as.character(chrNames)

  if (!all(chrNames %in% unlist(fread(cellFiles[1], select = 1, colClasses = 'character'))))
    stop('Chromosomes not all found in first cell file (`cellFiles[1]`)!')

  if (!file.exists(writeDir)) dir.create(writeDir)

  if (!sparseNAdrop) stop("Sorry, `sparseNAdrop=FALSE` has not been implemented yet!")

  if (!is.null(colData)) {
    if (nrow(colData) != length(cellFiles)) stop('Number of rows of colData should equal to length of cellFiles!')
  }

  message("Start processing chromosome-by-chromosome...")


  for (chr in chrNames) {
    message(paste0("...", chr, ": "), appendLF = FALSE)

    # Gather genomic coordinates
    pos_full <- extractCoord(cellFiles[1], chr, sep)
    for (i in 2:length(cellFiles)) {
      pos_temp <- extractCoord(cellFiles[i], chr, sep)
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
        M_mat <- cbind(
          M_mat,
          fillNA(cellFiles[i], chr, pos_full, sep) %>%
            as.matrix() %>%
            recommenderlab::dropNA()
        )
        # cat(i, " ")
      }
      gr <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = pos_full, end = pos_full))
      gr$meth <- as.integer(round(DelayedArray::rowSums(M_mat)))
      gr$total <- as.integer(DelayedArray::rowSums(M_mat > 0))
    } else {
      break # TODO
    }
    message("cells processed; ", appendLF = FALSE)

    # Write out processed data for current chromosome
    if (is.null(colData)) {
      se <- SummarizedExperiment::SummarizedExperiment(assays = list(M_mat = M_mat),
                                                       rowRanges = gr)
    } else {
      se <- SummarizedExperiment::SummarizedExperiment(assays = list(M_mat = M_mat),
                                                       rowRanges = gr,
                                                       colData = S4Vectors::DataFrame(colData))
    }
    se <- se[granges(se)$total > 0, ]

    n <- nchar(writeDir); if(substring(writeDir, n, n)=='/') writeDir <- substring(writeDir,1, n-1)
    HDF5Array::saveHDF5SummarizedExperiment(
      x = se,
      dir = paste0(writeDir, "/", chr),
      replace = TRUE
    )
    message("processed file written out.")
    if (sparseNAdrop) message("Important note: the `M_mat` assay is stored in NA-dropped sparse matrix format!")
  }
}

