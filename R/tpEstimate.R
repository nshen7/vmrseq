#' Estimate transition probilities conditioning on CpG-CpG distance
#'
#' @description This function first computes transition probilities
#' conditioning on CpG-CpG distance in each cell. Then the probs from
#' individual cells are smoothed over CpG-CpG distance using `loess`
#' with inverse-variance fitting.
#'
#' @param list a list of data.frame objects. Each data.frame
#' contains information of 1 unit of training data (can be a cell or a
#' subtype) and should have 3 columns in strict rder of: (chr), (pos),
#' (binary methyl value). Column names are not necessary.
#' @param max_dist_bp positive integer value indicating the maximum
#' CpG-CpG distance in base pairs before the transition probabilities
#' reach constant value. Default value is 2000.
#' @param buffer_bp length of buffer in base pairs used to fit smoothing
#' curve near maximum distance and plot diagnostics. Default is 3000.
#' @param lags a vector indicating possible number of lagged CpGs for
#' estimation of transition probabilities.
#' @param BPPARAM a \code{BiocParallelParam} object to specify the parallel
#' backend. The default option is \code{BiocParallel::bpparam()} which will
#' automatically creates a cluster appropriate for the operating system.
#' @param degree 'degree' argument for `loess` function
#' @param span 'span' argument for `loess` function
#' @param ... additional arguments passed into the `loess` function.
#'
#' @return a 'transitProbs' object. Postfixes rule in the output variables:
#' P(0|0) => '00'; P(0|1) => '01'; P(1|0) => '10'; P(1|1) => '11'.
#'
#' @importFrom dplyr %>% mutate lag ungroup select filter group_by summarise right_join arrange
#' @importFrom methods new
#' @importFrom stats loess
#' @importFrom BiocParallel bplapply register MulticoreParam bpparam
#'
#' @export
#' 
#' @examples
#' 
#' # load example data
#' data(cell_1, cell_2, cell_3)
#' 
#' # process the data to align with input format
#' library(dplyr)
#' df_1 <- cell_1 %>% mutate(meth = round(mc_count / total)) %>% select(chr, pos, meth)
#' df_2 <- cell_2 %>% mutate(meth = round(mc_count / total)) %>% select(chr, pos, meth)
#' df_3 <- cell_3 %>% mutate(meth = round(mc_count / total)) %>% select(chr, pos, meth)
#' 
#' # run tpEstimate
#' list <- list(df_1, df_2, df_3)
#' tpEstimate(list)
#'
tpEstimate <- function(list,
                        max_dist_bp = 2000,
                        buffer_bp = 3000,
                        lags = 1:10,
                        BPPARAM = bpparam(),
                        degree = 2, span = 0.02,
                        ...) {

  # Register the parallel backend
  BiocParallel::register(BPPARAM)
  backend <- paste0("BiocParallel:", class(bpparam())[1])

  for (i in length(list)) {
    list[[i]] <- list[[i]] %>% na.omit()
    if(!all(unlist(list[[i]][,3]) %in% 0:1))
      stop("Methylation value should be either integer 0 or 1.")
  }

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

  message("Computing transition probs within training unit.")
  # Proportions of sites in categories 00, 01, 10, 11 conditioning on
  # CpG distance
  if (parallel) {
    smr_units <- do.call(
      rbind,
      list %>% BiocParallel::bplapply(.computeProb1Unit,
                                      max_dist_bp = max_dist_bp,
                                      buffer_bp = buffer_bp,
                                      lags = lags)
    )
  } else {
    smr_units <- do.call(
      rbind,
      list %>% lapply(.computeProb1Unit,
                      max_dist_bp = max_dist_bp,
                      buffer_bp = buffer_bp,
                      lags = lags)
    )
  }

  return(
    .estimTransitProbsFromSummary(
      smr_units = smr_units,
      max_dist_bp = max_dist_bp, buffer_bp = buffer_bp,
      degree = degree, span = span,
      parallel = parallel, ...
    )
  )
}

# compute probs in one training unit
.computeProb1Unit <- function(df, max_dist_bp, buffer_bp, lags){

  colnames(df) <- c("chr", "pos", "state")

  .computeLag <- function(lag) {
    lag_df <- df %>%
      group_by(chr) %>%
      dplyr::mutate(state_lag = dplyr::lag(state, n = lag),
                    dist_bp = c(rep(NA, lag), diff(pos, lag = lag))) %>%
      dplyr::ungroup() %>%
      dplyr::select(state, state_lag, dist_bp) %>%
      dplyr::filter(dist_bp <= max_dist_bp + buffer_bp) %>%
      group_by(dist_bp) %>%
      summarise(N_00 = sum(state == 0 & state_lag==0), # number of sites from 0 to 0
                N_10 = sum(state == 1 & state_lag==0), # number of sites from 0 to 1
                N_01 = sum(state == 0 & state_lag==1), # number of sites from 1 to 0
                N_11 = sum(state == 1 & state_lag==1)  # number of sites from 1 to 1
      )
    return(lag_df)
  }

  smr_df <- do.call(rbind, lapply(lags, .computeLag)) %>%
    dplyr::group_by(dist_bp) %>%
    dplyr::summarise(N_00 = sum(N_00), # number of sites from 0 to 0
                     N_10 = sum(N_10), # number of sites from 0 to 1
                     N_01 = sum(N_01), # number of sites from 1 to 0
                     N_11 = sum(N_11)  # number of sites from 1 to 1
    ) %>%
    dplyr::mutate(p_00 = N_00 / (N_00 + N_10), # P(X_i = 0 | X_{i-1} = 0)
                  p_01 = N_01 / (N_01 + N_11) # P(X_i = 0 | X_{i-1} = 1)
    ) %>%
    dplyr::mutate(p_10 = 1 - p_00, # P(X_i = 1 | X_{i-1} = 0)
                  p_11 = 1 - p_01 # P(X_i = 1 | X_{i-1} = 1)
    ) %>%
    # dplyr::select(-c(N_00, N_10, N_01, N_11)) %>%
    dplyr::right_join(data.frame(dist_bp = 1:(max_dist_bp+buffer_bp)), by = "dist_bp") %>%
    dplyr::arrange(dist_bp)

  return(smr_df)
}

# loess fitting from cell-level summary data
.estimTransitProbsFromSummary <- function(smr_units, max_dist_bp, buffer_bp, degree, span, ...){

  stopifnot("max(smr_units$dist_bp) not equal to max_dist_bp + buffer_bp." =
              max(smr_units$dist_bp) == max_dist_bp + buffer_bp)

  message("Computing mean and var of the probs across cells.")
  train_data <- smr_units %>%
    dplyr::filter(dist_bp > 0) %>%
    dplyr::group_by(dist_bp) %>%
    dplyr::summarise(pbar_00 = mean(p_00, na.rm = TRUE),
              pbar_01 = mean(p_01, na.rm = TRUE),
              pbar_10 = mean(p_10, na.rm = TRUE),
              pbar_11 = mean(p_11, na.rm = TRUE),
              var_00 = var(p_00, na.rm = TRUE),
              var_01 = var(p_01, na.rm = TRUE),
              var_10 = var(p_10, na.rm = TRUE),
              var_11 = var(p_11, na.rm = TRUE)) %>%
    # dplyr::add_row(dist_bp = 1,
    #         pbar_00 = NA, pbar_01 = NA, pbar_10 = NA, pbar_11 = NA,
    #         var_00 = NA, var_01 = NA, var_10 = NA, var_11 = NA,
    #         .before = 1) %>%
    dplyr::right_join(data.frame(dist_bp = 1:(max_dist_bp+buffer_bp)), by = "dist_bp")

  message("Loess smoothing over probs.")

  .selectCols <- function(subscript) {
    # 'subscript' should be one of the strings '00', '01', '10', '11'
    train_cols <- train_data %>%
      dplyr::select(dist_bp,
                    paste0("pbar_", subscript),
                    paste0("var_", subscript)) %>%
      dplyr::filter(eval(parse(text = paste0("var_", subscript))) > 0)
    return(train_cols)
  }

  x <- 1:(max_dist_bp+buffer_bp) # starting from 1 so that row number equal to distance
  smoothed_probs <- data.frame(
    phat_00 = with(
      .selectCols('00'),
      loess(pbar_00 ~ log(dist_bp), weights = 1/var_00, degree = degree, span = span, ...) %>%
        predict(newdata = log(x))
    ),
    phat_01 = with(
      .selectCols('01'),
      loess(pbar_01 ~ log(dist_bp), weights = 1/var_01, degree = degree, span = span, ...) %>%
        predict(newdata = log(x))
    ),
    phat_10 = with(
      .selectCols('10'),
      loess(pbar_10 ~ log(dist_bp), weights = 1/var_10, degree = degree, span = span, ...) %>%
        predict(newdata = log(x))
    ),
    phat_11 = with(
      .selectCols('11'),
      loess(pbar_11 ~ log(dist_bp), weights = 1/var_11, degree = degree, span = span, ...) %>%
        predict(newdata = log(x))
    )
  )

  tp <- new("transitProbs",
            max_dist_bp = max_dist_bp, buffer_bp = buffer_bp,
            transit_probs = smoothed_probs[1:max_dist_bp,],
            buffer_probs = smoothed_probs[(1:buffer_bp)+max_dist_bp,],
            train = train_data %>% dplyr::select(-dist_bp))
  return(tp)
}

