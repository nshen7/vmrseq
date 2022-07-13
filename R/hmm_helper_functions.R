# For computing REFARRAY, i.e., array of log factorial of length `max_size` -- for the purpose of easy computation
.calRefArray <- function(max_size) cumsum(log(1:max_size))

# Outputs the log factorial with Stirling's approximation when int >= max_size
.logFactorial <- function(int, REFARRAY){
  # `int` is a non-negative integer
  if(int == 0) return(0)
  else if(int <= length(REFARRAY)) return(REFARRAY[int]) # = cumsum(log(1:int))
  else return(int*log(int) - int + 0.5*log(2*pi*int))
}

.logChoose <- function(n, k, REFARRAY){
  if (k > n) return(-Inf)
  else return(.logFactorial(n, REFARRAY) - .logFactorial(k, REFARRAY) - .logFactorial(n - k, REFARRAY))
}

# Computes the log-choose(n,k) matrix
.calChoiceArray <- function(REFARRAY){
  # To store the values when n or k equals 0, we add one dimension for both row and column
  # thus the choose(n, k) is CHOICEARRAY[n+1, k+1]
  max_size <- length(REFARRAY)
  CHOICEARRAY <- matrix(0, max_size+1, max_size+1)

  for(n in 1:(max_size+1))
    for(k in 1:n)
      CHOICEARRAY[n,k] <- .logChoose(n-1, k-1, REFARRAY)

  rownames(CHOICEARRAY) <- colnames(CHOICEARRAY) <- 0:max_size
  return(CHOICEARRAY)
}

# mod0_u <- readRDS("code/deprecated/estim_emiBetaPrior_ZIBBregression/model_unmeth_ZIBBregression.rds")
# mod0_m <- readRDS("code/deprecated/estim_emiBetaPrior_ZIBBregression/model_meth_BBregression.rds")
# data0_u <- readRDS("data/interim/deprecated/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_unmethClust.rds")
# data0_m <- readRDS("data/interim/deprecated/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_methClust.rds")

.priorParams <- function(med_cov, type,
                         mod_u = mod0_u, mod_m = mod0_m,
                         data_u = data0_u, data_m = data0_m){
  # Estimates parameters in beta-mixture priors based on number of cells.
  # `med_cov` is the median across-cell coverage in the input dataset.
  # `type` should be 'u' or 'm' indicating grouping type.

  stopifnot("`type` should be either 'u' or 'm'." = type %in% c('u','m'))
  stopifnot("`med_cov` should be positive." = med_cov > 0)

  if(type == 'u'){
    nu <- mod_u %>% predict(what = "nu", type = "response",
                            newdata = data.frame(med_cov = med_cov),
                            data = data_u)
    mu <- mod_u %>% predict(what = "mu", type = "response",
                            newdata = data.frame(med_cov = med_cov),
                            data = data_u)
    sigma <- mod_u %>% predict(what = "sigma", type = "response",
                               newdata = data.frame(med_cov = med_cov),
                               data = data_u)
    return(c(nu = nu, mu = mu, sigma = sigma))
  } else {
    mu <- mod_m %>% predict(what = "mu", type = "response",
                            newdata = data.frame(med_cov = med_cov),
                            data = data_m)
    sigma <- mod_m %>% predict(what = "sigma", type = "response",
                               newdata = data.frame(med_cov = med_cov),
                               data = data_m)
    return(c(mu = mu, sigma = sigma))
  }
}

# .priorParams <- function(med_cov, type){
#   # Estimates parameters in beta-mixture priors based on number of cells.
#   # `med_cov` is the median across-cell coverage in the input dataset.
#   # `type` should be 'u' or 'm' indicating grouping type.
#
#   stopifnot("`type` should be either 'u' or 'm'." = type %in% c('u','m'))
#   stopifnot("`med_cov` should be positive." = med_cov > 0)
#   smr <- fread(here::here("data/interim/estim_emiBetaPrior_inflated/emiBetaPrior_subtype_summary.csv"))
#
#   if(type == 'u'){
#     model_w <- lm(w_u ~ log(med_cov), data = smr)
#     w <- model_w %>% predict(newdata = data.frame(med_cov = med_cov)) %>% max(0) %>% min(1)
#
#     alpha <- lm(alpha_u ~ I(1/med_cov), data = smr) %>% predict(newdata = data.frame(med_cov = med_cov)) #%>% max(0.5)
#     beta <- lm(beta_u ~ I(1/med_cov), data = smr) %>% predict(newdata = data.frame(med_cov = med_cov)) #%>% max(3)
#   } else {
#     model_w <- lm(w_m ~ log(med_cov), data = smr)
#     w <- model_w %>% predict(newdata = data.frame(med_cov = med_cov)) %>% max(0) %>% min(1)
#
#     alpha <- lm(alpha_m ~ I(1/med_cov), data = smr) %>% predict(newdata = data.frame(med_cov = med_cov)) %>% max(1)
#     beta <- lm(beta_m ~ I(1/med_cov), data = smr) %>% predict(newdata = data.frame(med_cov = med_cov)) %>% max(0.01)
#   }
#
#   return(c(w, alpha, beta) %>% unname())
# }
#

# .calMethArray <- function(par_u, par_m, REFARRAY){
.calMethArray <- function(par_u, par_m, max_size){
  # To store the values when n or k equals 0, we add one dimension for both row and column
  # thus the betaBinom(n, k) - i.e. k methylated from n reads - is METHARRAY[n+1, k+1] or UNMETHARRAY[n+1, k+1]
  # `par_u` should be 3 parameters in beta-mixture prior (unmethylated): c(w_1, beta_1, beta_2), where the mixture is w_1*dbeta(1, beta_1) + (1-w_1)*dbeta(1, beta_2)
  # `par_m` should be 3 parameters in beta-mixture prior (methylated): c(w_1, alpha_1, alpha_2), where the mixture is w_1*dbeta(alpha_1, 1) + (1-w_1)*dbeta(alpha_2, 1)

  # max_size <- length(REFARRAY)
  # print(par_u); print(par_m)
  METHARRAY <- matrix(data = 0, nrow = max_size+1, ncol = max_size+1)
  UNMETHARRAY <- matrix(data = 0, nrow = max_size+1, ncol = max_size+1)
  for (n in 1:(max_size+1)) {
    for (k in 1:n) {
      UNMETHARRAY[n,k] <- dZIBB(x = k-1, mu = par_u['mu'], sigma = par_u['sigma'],
                                nu = par_u['nu'], bd = n-1)
      METHARRAY[n,k] <- dBB(x = k-1, mu = par_m['mu'], sigma = par_m['sigma'],
                            bd = n-1)
      # UNMETHARRAY[n,k] <- ifelse(k == 1,
      #                            yes = par_u[1] + (1-par_u[1])*dbb(0, n-1, par_u[2], par_u[3]),
      #                            no = (1-par_u[1])*dbb(k-1, n-1, par_u[2], par_u[3]))
      # METHARRAY[n,k] <- ifelse(k == n,
      #                          yes = par_m[1] + (1-par_m[1])*dbb(n-1, n-1, par_m[2], par_m[3]),
      #                          no = (1-par_m[1])*dbb(k-1, n-1, par_m[2], par_m[3]))
    }
  }
  rownames(METHARRAY) <- colnames(METHARRAY) <- 0:max_size
  rownames(UNMETHARRAY) <- colnames(UNMETHARRAY) <- 0:max_size

  return(list(METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY))
}

# ==== Functions for computing emission probability ====


## Function for computing emission probability (1-grouping)
.calEmissionProb1Grp <- function(state_1g, total_read, meth_read, METHARRAY, UNMETHARRAY)
  ifelse(state_1g, METHARRAY[total_read+1, meth_read+1], UNMETHARRAY[total_read+1, meth_read+1])
# Compute emission probability (in 1-grouping case) of observing `meth_read` methylated read in `total_read` total reads
# `state_1g` is the binary code of underlying methylation states of the grouping: 1 - methylated; 0 - unmethylated.
# Tests:
# .calEmissionProb1Grp(state_1g = 1, total_read = 20, meth_read = 0, METHARRAY, UNMETHARRAY)
# .calEmissionProb1Grp(state_1g = 1, total_read = 20, meth_read = 10, METHARRAY, UNMETHARRAY)
# .calEmissionProb1Grp(state_1g = 1, total_read = 20, meth_read = 20, METHARRAY, UNMETHARRAY)
# sum(map_dbl(0:20, ~.calEmissionProb1Grp(state_1g = 0, total_read = 20, meth_read = .x, METHARRAY, UNMETHARRAY)))
# sum(map_dbl(0:20, ~.calEmissionProb1Grp(state_1g = 1, total_read = 20, meth_read = .x, METHARRAY, UNMETHARRAY)))


.translateState2Grp <- function(state_2g){
  # Translates state_2g 0, 1 or 2 into binary grouping state vectors
  # Three states are encoded as: `0`->(0,0), `1`->(1,0), `2`->(1,1)
  state_2g <- as.integer(state_2g)
  # stopifnot("`state_2g` should be 0, 1 or 2." = state_2g %in% 0:2)
  group_state_vector <- c(0,0)
  if(state_2g %in% 1:2) group_state_vector[1] <- 1
  if(state_2g %in% 2) group_state_vector[2] <- 1
  return(group_state_vector)
}

.translateMethFrac2Grp <- function(state_2g, pi1){
  # Compute expected methylated fraction, e.g. if pi1 = 0.7, (0,1) is translated to 0.3
  group_state_vector <- .translateState2Grp(state_2g)
  return(sum(group_state_vector*c(pi1, 1-pi1)))
}

.calEmissionProb2Grp <- function(state_2g, total_read, meth_read, pi1,
                                 REFARRAY, CHOICEARRAY, METHARRAY, UNMETHARRAY){
  # Compute emission probability (in 2-grouping case) of observing `meth_read` methylated read in `total_read` total reads
  # `state_2g` is the integer code of underlying methylation states of the 2 groupings (`0` = (0,0), `1` = (1,0), `2` = (1,1)).
  # `pi1` is the prevalence of the methylated grouping

  max_size <- length(REFARRAY)
  # stopifnot("state_2g` is either 0, 1 or 2." = state_2g %in% 0:2)
  # stopifnot("Number of methylated reads succeeds coverage." = meth_read <= total_read)
  # stopifnot("`total_read` is greater than maximum expected coverage." = total_read <= max_size)

  max_size <- length(REFARRAY)
  p <- .translateMethFrac2Grp(state_2g, pi1)
  prob <- 0
  for (i in 0:total_read) {
    for (j in 0:meth_read) {
      if (j <= i & meth_read-j <= total_read-i) {
        if (p < 0.01 & i == 0) log_lik_1 <- 0
        else if (p > 0.99 & i == total_read) log_lik_1 <- 0
        else log_lik_1 <- CHOICEARRAY[total_read+1, i+1] + i*log(p) + (total_read-i)*log(1-p)
        log_lik_2 <- log(METHARRAY[i+1, j+1]) + log(UNMETHARRAY[total_read-i+1, meth_read-j+1])
        log_lik_i_j <- log_lik_1 + log_lik_2
        prob <- prob + exp(log_lik_i_j)
      }
    }
  }
  return(prob)
}
# Tests:
# .calEmissionProb2Grp(state_2g = 0, total_read = 20, meth_read = 0, pi1, REFARRAY, CHOICEARRAY, METHARRAY, UNMETHARRAY)
# .calEmissionProb2Grp(state_2g = 0, total_read = 20, meth_read = 10, pi1, REFARRAY, CHOICEARRAY, METHARRAY, UNMETHARRAY)
# .calEmissionProb2Grp(state_2g = 0, total_read = 20, meth_read = 20, pi1, REFARRAY, CHOICEARRAY, METHARRAY, UNMETHARRAY)
# sum(map_dbl(0:20, ~.calEmissionProb2Grp(state_2g = 0, total_read = 20, meth_read = .x, pi1, REFARRAY, CHOICEARRAY, METHARRAY, UNMETHARRAY)))
# sum(map_dbl(0:20, ~.calEmissionProb2Grp(state_2g = 1, total_read = 20, meth_read = .x, pi1, REFARRAY, CHOICEARRAY, METHARRAY, UNMETHARRAY)))
# sum(map_dbl(0:20, ~.calEmissionProb2Grp(state_2g = 2, total_read = 20, meth_read = .x, pi1, REFARRAY, CHOICEARRAY, METHARRAY, UNMETHARRAY)))




# ==== Function for loading transition probabilities ====

# tp0 <- read_rds(here::here("code/package_functions/transitProbs_27subtypes_1350cells_Luo2017&Liu2021.rds"))

# `probs` should be a data.frame of 4 columns [P(0|0),P(0|1),P(1|0),P(1|1)], while row number represents CpG-CpG distance
# `pos` is an atomic vector of genomic positions which shall be used to compute CpG-CpG distances and load transition probs
.loadTransitProbs <- function(pos, all_probs = tp0@transit_probs)
  all_probs[pmin(nrow(all_probs), diff(pos)), ]



# ==== Viterbi algorithm for 1- and 2-grouping cases ====

.Viterbi1Grp <- function(pos, total_reads, meth_reads, tp = NULL,
                         METHARRAY, UNMETHARRAY){
  if (any(c(length(total_reads), length(meth_reads)) != length(pos)))
    stop("Wrong dimensions of input data.")

  # compute transition probability
  if(is.null(tp)) trans_probs <- .loadTransitProbs(pos = pos)
  else trans_probs <- .loadTransitProbs(pos = pos, all_probs = tp0@transit_probs)

  # Two states are encoded as: `0`->unmethylated, `1`->methylated
  state_nums <- 0:1
  num_CpG <- length(total_reads)

  # Memory used during Viterbi computation and for storing results
  V <- matrix(-Inf, nrow = num_CpG, ncol = 2) # for storing log likelihood
  traceback <- matrix(0, nrow = num_CpG, ncol = 2) # for storing tracebback pointers

  # Distribution of initial state is set to uniform (proportional to 1)
  for (i in 1:num_CpG) {
    for (j in state_nums+1) {
      log_em_prob <- log(.calEmissionProb1Grp(state_1g = j-1,
                                              total_read = total_reads[i], meth_read = meth_reads[i],
                                              METHARRAY, UNMETHARRAY))
      if (i == 1) { ## likelihood for the first CpG
        V[i, j] <- log_em_prob
      } else {
        best_prior_state <- 0
        for (j_prior in state_nums+1) {
          log_tr_prob <- log(trans_probs[i-1, j_prior+2*j-2])
          oldVal <- V[i,j]
          newVal <- log_em_prob + log_tr_prob + V[i-1, j_prior]
          if (oldVal < newVal){
            V[i,j] <- newVal
            best_prior_state <- j_prior - 1
          }else if (oldVal == newVal){ ## tie-breaking in Viterbi
            if (sample(0:1,1) == 0) best_prior_state <- j_prior - 1
          }
        }
        traceback[i,j] <- best_prior_state
      }
    }
  }
  # print(V); print(traceback)

  # Viterbi traceback
  state_path <- rep(0, num_CpG)
  loglik_path <- rep(-Inf, num_CpG)

  state_path[num_CpG] <- which.max(V[num_CpG,])-1
  loglik_path[num_CpG] <- max(V[num_CpG,])

  for (i in (num_CpG-1):1) {
    state_path[i] <- traceback[i+1, state_path[i+1]+1]
    loglik_path[i] <- V[i, state_path[i]+1]
  }
  return(data.frame(state_path, loglik_path))
}
# Tests:
# pos = cumsum(sample(100:2000, 10))
# total_reads = sample(30:100, 10)
# .Viterbi1Grp(pos, total_reads, meth_reads = round(total_reads*0.1), tp = NULL, METHARRAY, UNMETHARRAY)
# .Viterbi1Grp(pos, total_reads, meth_reads = round(total_reads*0.9), tp = NULL, METHARRAY, UNMETHARRAY)



.Viterbi2Grp <- function(total_reads, meth_reads, trans_probs, pi1,
                         REFARRAY, CHOICEARRAY, METHARRAY, UNMETHARRAY){
  # should be applied on 1 region
  # `trans_probs` should have four columns:  [P(0|0), P(0|1), P(1|0), P(1|1)]

  if (any(c(length(total_reads), length(meth_reads)) != nrow(trans_probs)+1))
    stop("Wrong dimensions of input data.")

  # Three states are encoded as: `0`->(0,0), `1`->(1,0), `2`->(1,1)
  state_nums <- 0:2 # only allow three states
  num_group <- 2 # assuming 2 groupings
  num_CpG <- length(total_reads)

  # Memory used during Viterbi computation and for storing results
  V <- matrix(-Inf, nrow = num_CpG, ncol = 3) # for storing log likelihood
  traceback <- matrix(0, nrow = num_CpG, ncol = 3) # for storing tracebback pointers

  # Distribution of initial state is set to uniform (proportional to 1)
  for (i in 1:num_CpG) {
    for (j in state_nums+1) {
      log_em_prob <- log(.calEmissionProb2Grp(state_2g = j-1, total_read = total_reads[i], meth_read = meth_reads[i],
                                              pi1, REFARRAY, CHOICEARRAY, METHARRAY, UNMETHARRAY))
      if (i == 1) { ## likelihood for the first CpG
        V[i, j] <- log_em_prob
      } else {
        best_prior_state <- 0
        for (j_prior in state_nums+1) {
          log_tr_prob = 0
          for (k in 1:num_group) {
            j_prior_bin <- .translateState2Grp(j_prior-1)[k]
            j_bin <- .translateState2Grp(j-1)[k]
            log_tr_prob <- log_tr_prob + log(trans_probs[i-1, j_prior_bin+2*j_bin+1])
            # cat("k = ", k, "; log_tr_prob = ", log_tr_prob, "\n")
          }
          # cat(i, j, j_prior, log_em_prob, log_tr_prob, V[i-1, j_prior], "\n") # for debug
          oldVal <- V[i,j]
          newVal <- log_em_prob + log_tr_prob + V[i-1, j_prior]
          if (oldVal < newVal){
            V[i,j] <- newVal
            best_prior_state <- j_prior - 1
          }else if (oldVal == newVal){
            ## tie-breaking in Viterbi
            if (sample(0:1,1) == 0) best_prior_state <- j_prior - 1
          }
        }
        traceback[i,j] <- best_prior_state
        # cat(i,j,best_prior_state,"\n")
      }
    }
  }
  # print(V)
  # print(traceback)

  # Viterbi traceback
  state_path <- rep(0, num_CpG)
  bin_state_path <- matrix(0, nrow = num_CpG, ncol = 2)
  loglik_path <- rep(-Inf, num_CpG)

  state_path[num_CpG] <- which.max(V[num_CpG,])-1
  bin_state_path[num_CpG,] <- .translateState2Grp(state_path[num_CpG])
  loglik_path[num_CpG] <- max(V[num_CpG,])

  for (i in (num_CpG-1):1) {
    state_path[i] <- traceback[i+1, state_path[i+1]+1]
    bin_state_path[i, ] <- .translateState2Grp(state_path[i])
    loglik_path[i] <- V[i, state_path[i]+1]
  }

  colnames(bin_state_path) <- c("methState_group1", "methState_group2")
  return(data.frame(state_path, bin_state_path, loglik_path))
}
# Tests:
# pi1 = 0.3
# pos = cumsum(sample(100:2000, 10))
# total_reads = sample(30:100, 10)
# trans_probs = .loadTransitProbs(pos = pos)
#
# meth_reads = round(total_reads*0.2)
# .Viterbi2Grp(total_reads, meth_reads, trans_probs,
#          pi1, REFARRAY, CHOICEARRAY, METHARRAY, UNMETHARRAY)
#






# ==== Functions for prevalence optimization in 2-grouping case ====

## Functions for interim calculation in exponentiated gradient update

.gradntFunc <- function(x, total_read, meth_read, CHOICEARRAY, METHARRAY, UNMETHARRAY) {
  # function g_k(x) in the derivative of likelihood
  numer <- 0
  denom <- 0
  for (i in 0:total_read) {
    for (j in 0:meth_read) {
      if (j <= i & meth_read-j <= total_read-i) {
        log_c_ijk <- CHOICEARRAY[total_read+1, i+1] + log(METHARRAY[i+1, j+1]) + log(UNMETHARRAY[total_read-i+1, meth_read-j+1])
        numer <- numer + exp(log_c_ijk + (i-1)*log(x) + (total_read-i-1)*log(1-x)) * (i-x*total_read)
        denom <- denom + exp(log_c_ijk + i*log(x) + (total_read-i)*log(1-x))
      }
    }
  }
  val <- - numer / denom
  return(val)
}

.loglikGrad <- function(pi1, state_path, total_reads, meth_reads, CHOICEARRAY, METHARRAY, UNMETHARRAY) {
  # Gradient of the log-likelihood with respect to (\pi_1, \pi_2)
  if (any(c(length(state_path), length(total_reads)) != length(meth_reads)))
    stop("Wrong dimensions of input data.")

  grad1 <- 0#; grad2 <- 0 ## grad2 is always 0
  for (k in 1:length(state_path)) {
    state_2g <- state_path[k]
    total_read <- total_reads[k]
    meth_read <- meth_reads[k]

    if (state_2g == 1) { ## Three states are encoded as: `0`->(0,0), `1`->(1,0), `2`->(1,1)
      p <- .translateMethFrac2Grp(state_2g, pi1)
      grad1 <- grad1 + .gradntFunc(p, total_read, meth_read, CHOICEARRAY, METHARRAY, UNMETHARRAY)
    }
  }
  return(grad1)
}

.prevOptimSnglInit <- function(pos, total_reads, meth_reads, pi1_init, tp = NULL,
                               epsilon = 1e-4, backtrack = T, eta = ifelse(backtrack, 0.05, 0.005), max_iter = 200,
                               CHOICEARRAY, METHARRAY, UNMETHARRAY){
  # `tp` is an transitProbs object, which stores the transition probability distribution
  # `pi1_init` is the set initial value of \pi_1
  # `backtrack` is logical value indicating whether to use backtracking line search
  # `eta` is the learning rate (default: ifelse(backtrack, 0.05, 0.005))
  # `epsilon` is the convergence upper bound of \pi_1
  # `max_iter` is the maximum number of iteration
  # `eta` is the learning rate
  max_size <- nrow(CHOICEARRAY) - 1

  if(is.null(tp)) trans_probs <- .loadTransitProbs(pos = pos)
  else trans_probs <- .loadTransitProbs(pos = pos, all_probs = tp0@transit_probs)

  pi_1 <- pi1_init; pi_2 <- 1 - pi_1
  init_vit <- .Viterbi2Grp(total_reads, meth_reads, trans_probs, pi_1, REFARRAY, CHOICEARRAY, METHARRAY, UNMETHARRAY)
  state_path <- init_vit$state_path
  loglik <- init_vit$loglik_path[nrow(init_vit)]

  for (t in 1:max_iter) {
    old_pi_1 <- pi_1

    ## EG update (Experimented with unit test: slower when momentum is added)
    pi_2 <- 1 - pi_1
    grad <- .loglikGrad(pi_1, state_path, total_reads, meth_reads, CHOICEARRAY, METHARRAY, UNMETHARRAY)
    pi_1 <- pi_1 * exp(-eta*grad)
    pi_1 <- pi_1 / (pi_1 + pi_2)

    # viterbi update
    vit <- .Viterbi2Grp(total_reads, meth_reads, trans_probs, pi_1, REFARRAY, CHOICEARRAY, METHARRAY, UNMETHARRAY)
    old_loglik <- loglik
    loglik <- vit$loglik_path[nrow(init_vit)]

    # backtracking line search
    if(backtrack) {
      alpha <- 0.5; beta <- 0.7
      if(loglik > old_loglik - alpha*eta*grad^2) eta <- beta*eta
    }

    if(t == max_iter) print("Max iteration achieved.")
    if(abs(loglik - old_loglik) < epsilon & abs(old_pi_1 - pi_1) < epsilon) break
  }

  return(list(pi1_init = pi1_init, optim_pi_1 = pi_1, vit_path = vit[, 2:3], loglik = loglik, n_iter = t))
}
# Tests:
# pos = cumsum(sample(100:2000, 10))
# total_reads = sample(30:100, 10)
# meth_reads = round(total_reads*runif(length(pos), 0.5, 0.6))
# .prevOptimSnglInit(pos, total_reads, meth_reads, tp = NULL, pi1_init = 0.75,
#                    epsilon = 1e-4, backtrack = F, eta = 0.005, max_iter = 100,
#                    CHOICEARRAY = CHOICEARRAY, METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY)
# .prevOptimSnglInit(pos, total_reads, meth_reads, tp = NULL, pi1_init = 0.75,
#                    epsilon = 1e-4, backtrack = T, eta = 0.05, max_iter = 100,
#                    CHOICEARRAY = CHOICEARRAY, METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY)
# .prevOptimSnglInit(pos, total_reads, meth_reads, tp = NULL, pi1_init = 0.25,
#                    epsilon = 1e-4, backtrack = F, eta = 0.005, max_iter = 100,
#                    CHOICEARRAY = CHOICEARRAY, METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY)
# .prevOptimSnglInit(pos, total_reads, meth_reads, tp = NULL, pi1_init = 0.25,
#                    epsilon = 1e-4, backtrack = T, eta = 0.05, max_iter = 100,
#                    CHOICEARRAY = CHOICEARRAY, METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY)




.prevOptimMultiInit <- function(pos, total_reads, meth_reads, inits, tp = NULL,
                                epsilon = 1e-4, backtrack = T, eta = ifelse(backtrack, 0.05, 0.005), max_iter = 200,
                                CHOICEARRAY, METHARRAY, UNMETHARRAY){
  loglik <- -Inf; optim_pi_1 <- -1
  for (i in 1:length(inits)) {

    # ## Skip the next initial value if the current optimized prevalence is larger than it. e.g., skip init=0.5 if optim_pi_1=0.55 from init=0.25
    # if (i < length(inits) & optim_pi_1 >= inits[i]) next

    pi1_init <- inits[i]
    res_temp <- .prevOptimSnglInit(pos, total_reads, meth_reads, pi1_init, tp,
                                   epsilon, backtrack, eta, max_iter,
                                   CHOICEARRAY, METHARRAY, UNMETHARRAY)
    if (res_temp$loglik > loglik) {
      loglik <- res_temp$loglik
      res <- res_temp
      optim_pi_1 <- res$optim_pi_1
    }
  }
  # if (is.null(findHVR(res$vit_path, min_nCpG_inHVR))) res$optim_pi_1 <- NA
  return(res)
}
