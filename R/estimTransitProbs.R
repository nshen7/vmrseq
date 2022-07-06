#' Estimate transition probilities conditioning on CpG-CpG distance
#'
#' @description This function first computes transition probilities
#' conditioning on CpG-CpG distance in each cell. Then the probs from
#' individual cells are smoothed over CpG-CpG distance using `loess`
#' with inverse-variance fitting.
#'
#' @param list_cells a list of data.frame objects. Each data.frame
#' contains information of 1 cell and should have 3 columns in strict
#' order of: (chr), (pos), (binary methyl value). Column names are not
#' necessary.
#' @param max_dist_bp positive integer value indicating the maximum
#' CpG-CpG distance in base pairs before the transition probabilities
#' reach constant value. Default value is 2000.
#' @param buffer_bp length of buffer in base pairs used to fit smoothing
#' curve near maximum distance and plot diagnostics. Default is 3000.
#' @param BPPARAM a \code{BiocParallelParam} object to specify the parallel
#' backend. The default option is \code{BiocParallel::bpparam()} which will
#' automatically creates a cluster appropriate for the operating system.
#' @param ... additional arguments passed into the `loess` function.
#'
#' @return a 'transitProbs' object. Postfixes rule in the output variables:
#' P(0|0) => '00'; P(0|1) => '01'; P(1|0) => '10'; P(1|1) => '11'.
#'
#' @import dplyr
#' @import stats
#' @importFrom BiocParallel bplapply register MulticoreParam bpparam
#'
#' @export
#'
#' @examples
#'
estimTransitProbs <- function(list_cells,
                              max_dist_bp = 2000,
                              buffer_bp = 3000,
                              BPPARAM = bpparam(),
                              ...) {

  # Register the parallel backend
  BiocParallel::register(BPPARAM)
  backend <- paste0("BiocParallel:", class(bpparam())[1])

  if (bpparam()$workers == 1) {
    mes <- "Using a single core (backend: %s)."
    message(sprintf(mes, backend))
    parallel <- FALSE
  } else {
    mes <- paste0("Parallelizing using %s workers/cores ",
                  "(backend: %s).")
    message(sprintf(mes, bpparam()$workers, backend))
    parallel <- TRUE
  }

  message("Computing transition probs within cells.")
  # Proportions of sites in categories 00, 01, 10, 11 conditioning on
  # CpG distance
  smr_cells <- do.call(rbind,
                       list_cells %>% BiocParallel::bplapply(.computeProb1Cell,
                                               max_dist_bp = max_dist_bp,
                                               buffer_bp = buffer_bp)
  )

  return(.estimTransitProbsFromSummary(
    smr_cells = smr_cells,
    max_dist_bp = max_dist_bp, buffer_bp = buffer_bp,
    parallel = parallel, ...
    )
  )
}


.computeProb1Cell <- function(df, max_dist_bp, buffer_bp){

  df <- df %>% na.omit()
  if(!all(df[,3] %in% 0:1)) stop("Methylation value should be either 0 or 1.")

  colnames(df) <- c("chr", "pos", "MF")
  smr_df <- df %>%
    group_by(chr) %>%
    mutate(MF_lag1 = lag(MF), dist_bp = c(NA, diff(pos))) %>%
    filter(dist_bp <= max_dist_bp + buffer_bp) %>%
    group_by(dist_bp) %>%
    summarise(N_00 = sum(MF == 0 & MF_lag1==0), # number of sites from 0 to 0
              N_10 = sum(MF == 1 & MF_lag1==0), # number of sites from 0 to 1
              N_01 = sum(MF == 0 & MF_lag1==1), # number of sites from 1 to 0
              N_11 = sum(MF == 1 & MF_lag1==1)  # number of sites from 1 to 1
    ) %>%
    mutate(p_00 = N_00 / (N_00 + N_10), # P(X_i = 0 | X_{i-1} = 0)
           p_01 = N_01 / (N_01 + N_11) # P(X_i = 0 | X_{i-1} = 1)
    ) %>%
    mutate(p_10 = 1 - p_00, # P(X_i = 1 | X_{i-1} = 0)
           p_11 = 1 - p_01 # P(X_i = 1 | X_{i-1} = 1)
    ) %>%
    dplyr::select(-c(N_00, N_10, N_01, N_11))
  return(smr_df)
}

.estimTransitProbsFromSummary <- function(smr_cells, max_dist_bp, buffer_bp, ...){

  stopifnot("max(smr_cells$dist_bp) not equal to max_dist_bp + buffer_bp." =
              max(smr_cells$dist_bp) == max_dist_bp + buffer_bp)

  message("Computing mean and var of the probs across cells.")
  train_data <- smr_cells %>%
    group_by(dist_bp) %>%
    summarise(pbar_00 = mean(p_00, na.rm = T),
              pbar_01 = mean(p_01, na.rm = T),
              pbar_10 = mean(p_10, na.rm = T),
              pbar_11 = mean(p_11, na.rm = T),
              var_00 = var(p_00, na.rm = T),
              var_01 = var(p_01, na.rm = T),
              var_10 = var(p_10, na.rm = T),
              var_11 = var(p_11, na.rm = T)) %>%
    add_row(dist_bp = 1,
            pbar_00 = NA, pbar_01 = NA, pbar_10 = NA, pbar_11 = NA,
            var_00 = NA, var_01 = NA, var_10 = NA, var_11 = NA,
            .before = 1)

  message("Loess smoothing over probs.")
  x <- 1:(max_dist_bp+buffer_bp) # starting from 1 so that row number equal to distance
  smoothed_probs <- with(train_data %>% na.omit() %>% filter(var_00!=0, var_01!=0, var_10!=0, var_11!=0),
                         data.frame(phat_00 = loess(pbar_00 ~ log(dist_bp), weights = 1/var_00, ...) %>% predict(newdata = log(x)),
                                    phat_01 = loess(pbar_01 ~ log(dist_bp), weights = 1/var_01, ...) %>% predict(newdata = log(x)),
                                    phat_10 = loess(pbar_10 ~ log(dist_bp), weights = 1/var_10, ...) %>% predict(newdata = log(x)),
                                    phat_11 = loess(pbar_11 ~ log(dist_bp), weights = 1/var_11, ...) %>% predict(newdata = log(x))
                         )
  )
  tp <- new("transitProbs",
            max_dist_bp = max_dist_bp, buffer_bp = buffer_bp,
            transit_probs = smoothed_probs[1:max_dist_bp,],
            buffer_probs = smoothed_probs[(1:buffer_bp)+max_dist_bp,],
            train = train_data %>% dplyr::select(-dist_bp))
  return(tp)
}

