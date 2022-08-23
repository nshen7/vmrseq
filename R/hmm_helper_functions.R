###### utils =====
# For computing REFARRAY, i.e., array of log factorial of length `max_cov` -- for the purpose of easy computation
.calRefArray <- function(max_cov) cumsum(log(1:max_cov))

# Outputs the log factorial with Stirling's approximation when int >= max_cov
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
  max_cov <- length(REFARRAY)
  CHOICEARRAY <- matrix(0, max_cov+1, max_cov+1)

  for(n in 1:(max_cov+1))
    for(k in 1:n)
      CHOICEARRAY[n,k] <- .logChoose(n-1, k-1, REFARRAY)

  rownames(CHOICEARRAY) <- colnames(CHOICEARRAY) <- 0:max_cov
  return(CHOICEARRAY)
}

# Estimates parameters in beta-mixture priors based on number of cells.
.priorParams <- function(med_cov, type){
  # `med_cov` is the median across-cell coverage in the input dataset.
  # `type` should be 'u' or 'm' indicating grouping type.

  stopifnot("`type` should be either 'u' or 'm'." = type %in% c('u','m'))
  # stopifnot("`med_cov` should be a positive integer." = med_cov > 0 & round(med_cov)==med_cov)
  med_cov <- round(med_cov) %>% max(1)

  if (type == 'u') {
    pars <- params_u
    max_mc <- max(pars$med_cov)
    if (med_cov > max_mc) {
      return(c(nu = pars$nu[max_mc], mu = pars$mu[max_mc], sigma = pars$sigma[max_mc]))
    } else {
      return(c(nu = pars$nu[med_cov], mu = pars$mu[med_cov], sigma = pars$sigma[med_cov]))
    }
  } else {
    pars <- params_m
    return(c(mu = pars$mu[1], sigma = pars$sigma[1]))
  }
}

.calMethArray <- function(par_u, par_m, max_cov){
  # To store the values when n or k equals 0, we add one dimension for both row and column
  # thus the betaBinom(n, k) - i.e. k methylated from n reads - is METHARRAY[n+1, k+1] or UNMETHARRAY[n+1, k+1]
  # `par_u` should be 3 parameters in beta-mixture prior (unmethylated): c(w_1, beta_1, beta_2), where the mixture is w_1*dbeta(1, beta_1) + (1-w_1)*dbeta(1, beta_2)
  # `par_m` should be 3 parameters in beta-mixture prior (methylated): c(w_1, alpha_1, alpha_2), where the mixture is w_1*dbeta(alpha_1, 1) + (1-w_1)*dbeta(alpha_2, 1)

  METHARRAY <- matrix(data = 0, nrow = max_cov+1, ncol = max_cov+1)
  UNMETHARRAY <- matrix(data = 0, nrow = max_cov+1, ncol = max_cov+1)
  for (n in 1:(max_cov+1)) {
    for (k in 1:n) {
      UNMETHARRAY[n,k] <- gamlss.dist::dZIBB(x = k-1, mu = par_u['mu'], sigma = par_u['sigma'],
                                nu = par_u['nu'], bd = n-1)
      METHARRAY[n,k] <- gamlss.dist::dBB(x = k-1, mu = par_m['mu'], sigma = par_m['sigma'],
                            bd = n-1)
    }
  }
  rownames(METHARRAY) <- colnames(METHARRAY) <- 0:max_cov
  rownames(UNMETHARRAY) <- colnames(UNMETHARRAY) <- 0:max_cov

  return(list(METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY))
}

###### functions for computing emission probability ====


# Function for computing emission probability (1-grouping)
.calEmissionProb1Grp <- function(state_1g, total, meth, METHARRAY, UNMETHARRAY)
  ifelse(state_1g, METHARRAY[total+1, meth+1], UNMETHARRAY[total+1, meth+1])
# Compute emission probability (in 1-grouping case) of observing `meth` methylated read in `total` total reads
# `state_1g` is the binary code of underlying methylation states of the grouping: 1 - methylated; 0 - unmethylated.

# Tests:
# .calEmissionProb1Grp(state_1g = 1, total = 20, meth = 0, METHARRAY, UNMETHARRAY)
# .calEmissionProb1Grp(state_1g = 1, total = 20, meth = 10, METHARRAY, UNMETHARRAY)
# .calEmissionProb1Grp(state_1g = 1, total = 20, meth = 20, METHARRAY, UNMETHARRAY)
# sum(map_dbl(0:20, ~.calEmissionProb1Grp(state_1g = 0, total = 20, meth = .x, METHARRAY, UNMETHARRAY)))
# sum(map_dbl(0:20, ~.calEmissionProb1Grp(state_1g = 1, total = 20, meth = .x, METHARRAY, UNMETHARRAY)))


# Translates state_2g 0, 1 or 2 into binary grouping state vectors
.translateState2Grp <- function(state_2g){
  # Three states are encoded as: `0`->(0,0), `1`->(1,0), `2`->(1,1)
  state_2g <- as.integer(state_2g)
  # stopifnot("`state_2g` should be 0, 1 or 2." = state_2g %in% 0:2)
  group_state_vector <- c(0,0)
  if(state_2g %in% 1:2) group_state_vector[1] <- 1
  if(state_2g %in% 2) group_state_vector[2] <- 1
  return(group_state_vector)
}

# Compute expected methylated fraction, e.g. if pi1 = 0.7, (0,1) is translated to 0.3
.translateMethFrac2Grp <- function(state_2g, pi1){
  group_state_vector <- .translateState2Grp(state_2g)
  return(sum(group_state_vector*c(pi1, 1-pi1)))
}

.calEmissionProb2Grp <- function(state_2g, total, meth, pi1,
                                 CHOICEARRAY, METHARRAY, UNMETHARRAY){
  # Compute emission probability (in 2-grouping case) of observing `meth` methylated read in `total` total reads
  # `state_2g` is the integer code of underlying methylation states of the 2 groupings (`0` = (0,0), `1` = (1,0), `2` = (1,1)).
  # `pi1` is the prevalence of the methylated grouping

  max_cov <- nrow(CHOICEARRAY) - 1
  p <- .translateMethFrac2Grp(state_2g, pi1)
  # cat("state =", state_2g, "pi1 =", pi1, "p =", p, "\n") ## DEBUG
  prob <- 0
  for (i in 0:total) {
    for (j in 0:meth) {
      if (j <= i & meth-j <= total-i) {
        if (p < 0.0001 & i == 0) log_lik_1 <- 0
        else if (p > 0.9999 & i == total) log_lik_1 <- 0
        else log_lik_1 <- CHOICEARRAY[total+1, i+1] + i*log(p) + (total-i)*log(1-p)
        log_lik_2 <- log(METHARRAY[i+1, j+1]) + log(UNMETHARRAY[total-i+1, meth-j+1])
        log_lik_i_j <- log_lik_1 + log_lik_2
        prob <- prob + exp(log_lik_i_j)
      }
    }
  }
  return(prob)
}
# Tests:
# .calEmissionProb2Grp(state_2g = 0, total = 20, meth = 0, pi1, CHOICEARRAY, METHARRAY, UNMETHARRAY)
# .calEmissionProb2Grp(state_2g = 0, total = 20, meth = 10, pi1, CHOICEARRAY, METHARRAY, UNMETHARRAY)
# .calEmissionProb2Grp(state_2g = 0, total = 20, meth = 20, pi1, CHOICEARRAY, METHARRAY, UNMETHARRAY)
# sum(map_dbl(0:20, ~.calEmissionProb2Grp(state_2g = 0, total = 20, meth = .x, pi1, CHOICEARRAY, METHARRAY, UNMETHARRAY)))
# sum(map_dbl(0:20, ~.calEmissionProb2Grp(state_2g = 1, total = 20, meth = .x, pi1, CHOICEARRAY, METHARRAY, UNMETHARRAY)))
# sum(map_dbl(0:20, ~.calEmissionProb2Grp(state_2g = 2, total = 20, meth = .x, pi1, CHOICEARRAY, METHARRAY, UNMETHARRAY)))




# Function for loading transition probabilities from tp
.loadTransitProbs <- function(pos, all_probs = tp0@transit_probs)
  # `probs` should be a data.frame of 4 columns [P(0|0),P(0|1),P(1|0),P(1|1)], while row number represents CpG-CpG distance
  # `pos` is an atomic vector of genomic positions which shall be used to compute CpG-CpG distances and load transition probs
  all_probs[pmin(nrow(all_probs), diff(pos)), ]



###### Viterbi algorithm for 1- and 2-grouping cases =====
.Viterbi1Grp <- function(totals, meths, trans_probs, METHARRAY, UNMETHARRAY){
  if (any(c(length(totals), length(meths)) != nrow(trans_probs)+1))
    stop("Wrong dimensions of input data.")

  # Two states are encoded as: `0`->unmethylated, `1`->methylated
  state_nums <- 0:1
  num_cpg <- length(totals)

  # Memory used during Viterbi computation and for storing results
  V <- matrix(-Inf, nrow = num_cpg, ncol = 2) # for storing total log likelihood
  P <- matrix(-Inf, nrow = num_cpg, ncol = 2) # for storing transition and emission probs
  colnames(P) <- c("emission", "transition")
  traceback <- matrix(0, nrow = num_cpg, ncol = 2) # for storing tracebback pointers

  # Distribution of initial state is set to uniform (proportional to 1)
  for (i in 1:num_cpg) {
    for (j in state_nums+1) {
      log_em_prob <- log(.calEmissionProb1Grp(state_1g = j-1,
                                              total = totals[i], meth = meths[i],
                                              METHARRAY, UNMETHARRAY))
      if (i == 1) { ## likelihood for the first CpG
        V[i, j] <- log_em_prob
        P[i, 'emission'] <- log_em_prob
        P[i, 'transition'] <- NA
      } else {
        best_prior_state <- 0
        for (j_prior in state_nums+1) {
          log_tr_prob <- log(trans_probs[i-1, j_prior+2*j-2])
          oldVal <- V[i,j]
          newVal <- log_em_prob + log_tr_prob + V[i-1, j_prior]
          if (oldVal < newVal) {
            V[i,j] <- newVal
            P[i, 'emission'] <- log_em_prob
            P[i, 'transition'] <- log_tr_prob
            best_prior_state <- j_prior - 1
          } else if (oldVal == newVal) { ## tie-breaking in Viterbi
            if (sample(0:1,1) == 0) best_prior_state <- j_prior - 1
          }
        }
        traceback[i,j] <- best_prior_state
      }
    }
  }
  # print(V); print(traceback)

  # Viterbi traceback
  state_path <- rep(0, num_cpg)
  loglik_path <- rep(-Inf, num_cpg)

  state_path[num_cpg] <- which.max(V[num_cpg,])-1
  loglik_path[num_cpg] <- max(V[num_cpg,])

  for (i in (num_cpg-1):1) {
    state_path[i] <- traceback[i+1, state_path[i+1]+1]
    loglik_path[i] <- V[i, state_path[i]+1]
  }
  return(data.frame(state_path, loglik_path, P))
}


.calTransProb2Grp <- function(i, trans_probs) {
  # Compute 2-group transition probability matrix for one CpG pair
  # All possible states in 2-group: (0,0), (1,0), (1,1)
  # Hence row/cols of prob mat for first grouping: 0, 1, 1; for second grouping: 0, 0, 1

  state1g_mat_grp1 <- matrix(c(0, 1, 1), nrow = 3, ncol = 3)
  state2g_mat_grp1 <- state1g_mat_grp1 + 2*t(state1g_mat_grp1) + 1

  state1g_mat_grp2 <- matrix(c(0, 0, 1), nrow = 3, ncol = 3)
  state2g_mat_grp2 <- state1g_mat_grp2 + 2*t(state1g_mat_grp2) + 1

  # Transition prob mat (rows = 'from', cols = 'to')
  prob_mat <- apply(state2g_mat_grp1, 1:2, function(j) trans_probs[i, j]) *
    apply(state2g_mat_grp2, 1:2, function(j) trans_probs[i, j])
  prob_mat <- prob_mat / rowSums(prob_mat) # Scale to normalized probs (each row adds up to 1)

  return(prob_mat)
}


.Viterbi2Grp <- function(totals, meths, trans_probs, pi1,
                         CHOICEARRAY, METHARRAY, UNMETHARRAY){
  # should be applied on 1 region
  # `trans_probs` should have four columns:  [P(0|0), P(0|1), P(1|0), P(1|1)]

  if (any(c(length(totals), length(meths)) != nrow(trans_probs)+1))
    stop("Wrong dimensions of input data.")

  # Three states are encoded as: `0`->(0,0), `1`->(1,0), `2`->(1,1)
  state_nums <- 0:2 # only allow three states
  num_group <- 2 # assuming 2 groupings
  num_cpg <- length(totals)

  # Memory used during Viterbi computation and for storing results
  V <- matrix(-Inf, nrow = num_cpg, ncol = 3) # for storing log likelihood
  P <- matrix(-Inf, nrow = num_cpg, ncol = 2) # for storing transition and emission probs
  colnames(P) <- c("emission", "transition")
  traceback <- matrix(0, nrow = num_cpg, ncol = 3) # for storing tracebback pointers

  # Distribution of initial state is set to uniform (proportional to 1)
  for (i in 1:num_cpg) {

    # Compute transition prob matrix for 3-state model
    if(i > 1) tp_i_mat <- .calTransProb2Grp(i-1, trans_probs)

    for (j in state_nums+1) {
      log_em_prob <- log(.calEmissionProb2Grp(state_2g = j-1, total = totals[i], meth = meths[i],
                                              pi1, CHOICEARRAY, METHARRAY, UNMETHARRAY))
      if (i == 1) { ## likelihood for the first CpG
        V[i, j] <- log_em_prob
        P[i, 'emission'] <- log_em_prob
        P[i, 'transition'] <- NA
      } else {
        best_prior_state <- 0
        for (j_prior in state_nums+1) {

          # log_tr_prob <- 0
          # for (k in 1:num_group) {
          #   j_prior_bin <- .translateState2Grp(j_prior-1)[k]
          #   j_bin <- .translateState2Grp(j-1)[k]
          #   log_tr_prob <- log_tr_prob + log(trans_probs[i-1, j_prior_bin+2*j_bin+1])
          #   # cat("k = ", k, "; log_tr_prob = ", log_tr_prob, "\n")
          # }

          log_tr_prob <- log(tp_i_mat[j_prior, j])
          # cat(i, j, j_prior, exp(log_tr_prob), "\n") # DEBUG

          oldVal <- V[i,j]
          newVal <- log_em_prob + log_tr_prob + V[i-1, j_prior]
          if (oldVal < newVal){
            V[i,j] <- newVal
            P[i, 'emission'] <- log_em_prob
            P[i, 'transition'] <- log_tr_prob
            best_prior_state <- j_prior - 1
          }else if (oldVal == newVal){
            ## tie-breaking in Viterbi
            if (sample(0:1,1) == 0) best_prior_state <- j_prior - 1
          }
        }
        traceback[i,j] <- best_prior_state
        # cat(i, j, best_prior_state,"\n") # DEBUG
      }
    }
  }
  # print(V)
  # print(traceback)

  # Viterbi traceback
  state_path <- rep(0, num_cpg)
  bin_state_path <- matrix(0, nrow = num_cpg, ncol = 2)
  loglik_path <- rep(-Inf, num_cpg)

  state_path[num_cpg] <- which.max(V[num_cpg,])-1
  bin_state_path[num_cpg,] <- .translateState2Grp(state_path[num_cpg])
  loglik_path[num_cpg] <- max(V[num_cpg,])

  for (i in (num_cpg-1):1) {
    state_path[i] <- traceback[i+1, state_path[i+1]+1]
    bin_state_path[i, ] <- .translateState2Grp(state_path[i])
    loglik_path[i] <- V[i, state_path[i]+1]
  }

  colnames(bin_state_path) <- c("methState_group1", "methState_group2")
  return(data.frame(state_path, bin_state_path, loglik_path, P))
}
# Tests:
# pi1 = 0.3; N = 10
# pos = cumsum(sample(100:2000, 10))
# totals = sample(1:N, 10)
# trans_probs = .loadTransitProbs(pos = pos, all_probs = vmrseq:::tp0@transit_probs)
# meths = round(totals*0.3)
# REFARRAY <- .calRefArray(max_cov = N)
# CHOICEARRAY <- .calChoiceArray(REFARRAY)
# list <- .calMethArray(par_u = .priorParams(med_cov = round(median(totals)), type = "u"),
#                       par_m = .priorParams(med_cov = round(median(totals)), type = "m"),
#                       max_cov = N)
# METHARRAY <- list$METHARRAY; UNMETHARRAY <- list$UNMETHARRAY
# .Viterbi2Grp(totals, meths, trans_probs,
#              pi1, CHOICEARRAY, METHARRAY, UNMETHARRAY)
# .Viterbi1Grp(totals, meths, trans_probs, METHARRAY, UNMETHARRAY)



###### functions for prevalence optimization in 2-grouping case ====

# Function g_k(x) in the derivative of likelihood
.gradntFunc <- function(x, total, meth, CHOICEARRAY, METHARRAY, UNMETHARRAY) {
  numer <- 0
  denom <- 0
  for (i in 0:total) {
    for (j in 0:meth) {
      if (j <= i & meth-j <= total-i) {
        log_c_ijk <- CHOICEARRAY[total+1, i+1] + log(METHARRAY[i+1, j+1]) + log(UNMETHARRAY[total-i+1, meth-j+1])
        numer <- numer + exp(log_c_ijk + (i-1)*log(x) + (total-i-1)*log(1-x)) * (i-x*total)
        denom <- denom + exp(log_c_ijk + i*log(x) + (total-i)*log(1-x))
      }
    }
  }
  val <- - numer / denom
  return(val)
}

# Gradient of the log-likelihood with respect to (\pi_1, \pi_2)
.loglikGrad <- function(pi1, state_path, totals, meths, CHOICEARRAY, METHARRAY, UNMETHARRAY) {
  if (any(c(length(state_path), length(totals)) != length(meths)))
    stop("Wrong dimensions of input data.")

  grad1 <- 0#; grad2 <- 0 ## grad2 is always 0
  for (k in 1:length(state_path)) {
    state_2g <- state_path[k]
    total <- totals[k]
    meth <- meths[k]

    if (state_2g == 1) { ## Three states are encoded as: `0`->(0,0), `1`->(1,0), `2`->(1,1)
      p <- .translateMethFrac2Grp(state_2g, pi1)
      grad1 <- grad1 + .gradntFunc(p, total, meth, CHOICEARRAY, METHARRAY, UNMETHARRAY)
    }
  }
  return(grad1)
}

.prevOptimSnglInit <- function(pos, totals, meths, pi1_init, trans_probs,
                               epsilon = 1e-3, backtrack = T, eta = ifelse(backtrack, 0.05, 0.005), max_iter = 200,
                               CHOICEARRAY, METHARRAY, UNMETHARRAY){
  # `tp` is an transitProbs object, which stores the transition probability distribution
  # `pi1_init` is the set initial value of \pi_1
  # `backtrack` is logical value indicating whether to use backtracking line search
  # `eta` is the learning rate (default: ifelse(backtrack, 0.05, 0.005))
  # `epsilon` is the convergence upper bound of \pi_1
  # `max_iter` is the maximum number of iteration
  # `eta` is the learning rate
  max_cov <- nrow(CHOICEARRAY) - 1

  pi_1 <- pi1_init; pi_2 <- 1 - pi_1
  init_vit <- .Viterbi2Grp(totals, meths, trans_probs, pi_1, CHOICEARRAY, METHARRAY, UNMETHARRAY)
  state_path <- init_vit$state_path
  old_loglik <- loglik <- init_vit$loglik_path[nrow(init_vit)]

  for (t in 1:max_iter) {
    # print(t) # DEBUG
    old_pi_1 <- pi_1

    ## EG update (Experimented with unit test: slower when momentum is added)
    pi_2 <- 1 - pi_1
    grad <- .loglikGrad(pi_1, state_path, totals, meths, CHOICEARRAY, METHARRAY, UNMETHARRAY)

    # Decrease learning rate (eta) if updated pi_1 goes to 0 or 1
    try_pi_1 <- pi_1 * exp(-eta*grad)
    try_pi_1 <- try_pi_1 / (try_pi_1 + pi_2)
    # print(try_pi_1) # DEBUG
    while (try_pi_1 < 0.01 | try_pi_1 > 0.99) {
      eta <- 0.1 * eta
      try_pi_1 <- pi_1 * exp(-eta*grad)
      try_pi_1 <- try_pi_1 / (try_pi_1 + pi_2)
      # print(try_pi_1) # DEBUG
    }

    # Viterbi update
    vit <- .Viterbi2Grp(totals, meths, trans_probs, try_pi_1, CHOICEARRAY, METHARRAY, UNMETHARRAY)
    if (vit$loglik_path[nrow(init_vit)] <= loglik) { # decrease learning rate if likelihood does NOT increase
      eta <- 0.5 * eta
      try_pi_1 <- pi_1 * exp(-eta*grad)
      try_pi_1 <- try_pi_1 / (try_pi_1 + pi_2)
      pi_1 <- try_pi_1
    } else { # only adopt updated pi if likelihood increases
      old_loglik <- loglik
      loglik <- vit$loglik_path[nrow(init_vit)]
      state_path <- vit$state_path
      pi_1 <- try_pi_1
      if(backtrack) { # backtracking line search
        alpha <- 0.5; beta <- 0.7
        if(loglik > old_loglik - alpha*eta*grad^2) eta <- beta*eta
      }
    }

    # cat(pi_1, loglik, "\n") # DEBUG
    # if(t == max_iter) print("Max iteration achieved.")
    if(abs(loglik - old_loglik) < epsilon & abs(old_pi_1 - pi_1) < epsilon) break
  }

  return(list(pi1_init = pi1_init, optim_pi_1 = pi_1, vit_path = vit[, -1], loglik = loglik, n_iter = t))
}



###### wrapper for 1- and 2-group optimization ======

.solve1Grp <- function(pos, totals, meths, tp, METHARRAY, UNMETHARRAY) {
  trans_probs <- .loadTransitProbs(pos = pos, all_probs = tp@transit_probs)
  vit <- .Viterbi1Grp(totals, meths, trans_probs, METHARRAY, UNMETHARRAY)
  return(list(vit_path = vit, loglik = vit[nrow(vit), 'loglik_path']))
}

.solve2Grp <- function(gradient,
                       pos, totals, meths, tp,
                       inits, epsilon, backtrack,
                       eta, max_iter,
                       CHOICEARRAY, METHARRAY, UNMETHARRAY){

  trans_probs <- .loadTransitProbs(pos = pos, all_probs = tp@transit_probs)

  if (gradient) {
    loglik <- -Inf; optim_pi_1 <- -1
    for (i in 1:length(inits)) {
      # cat(i, "\n") # DEBUG
      pi1_init <- inits[i]
      res_temp <- .prevOptimSnglInit(pos, totals, meths, pi1_init, trans_probs,
                                     epsilon, backtrack, eta, max_iter,
                                     CHOICEARRAY, METHARRAY, UNMETHARRAY)
      if (res_temp$loglik > loglik) {
        loglik <- res_temp$loglik
        res <- res_temp
        optim_pi_1 <- res$optim_pi_1
      }
    }
  } else {
    loglik <- -Inf
    for (i in 1:length(inits)) {
      pi1 <- inits[i]
      vit <- .Viterbi2Grp(totals, meths, trans_probs, pi1,
                          CHOICEARRAY, METHARRAY, UNMETHARRAY)
      if (vit$loglik_path[nrow(vit)] > loglik) {
        loglik <- vit$loglik_path[nrow(vit)]
        res <- list("pi1_init" = pi1,
                    "optim_pi_1" = pi1,
                    "vit_path" = vit[, -1],
                    "loglik" = vit$loglik_path[nrow(vit)],
                    "n_iter" = 0)
      }
    }
  }
  return(res)
}

# Tests:
# pos = cumsum(sample(2:500, 100))
# totals = sample(5:20, 100, replace = T)
# meths = round(totals*runif(length(pos), 0.1, 0.9))
# REFARRAY <- .calRefArray(max_cov = max(totals))
# CHOICEARRAY <- .calChoiceArray(REFARRAY)
# list <- .calMethArray(par_u = .priorParams(med_cov = round(median(totals)), type = "u"),
#                       par_m = .priorParams(med_cov = round(median(totals)), type = "m"),
#                       max_cov = N)
# METHARRAY <- list$METHARRAY; UNMETHARRAY <- list$UNMETHARRAY
# .solve1Grp(pos, totals, meths, tp0, METHARRAY, UNMETHARRAY)
# .solve2Grp(pos, totals, meths, inits = c(.2,.5,.8), tp0,
#            epsilon = 1e-3, backtrack = T,
#            eta = ifelse(backtrack, 0.05, 0.005), max_iter = 200,
#            CHOICEARRAY, METHARRAY, UNMETHARRAY)


# functions for detect VMRs in predicted state sequence
.callVMR <- function(state_seq_2g, min_n, max_n_merge){

  state_seq_2g <- as.data.frame(state_seq_2g) # Formatting

  is_vml <- as.logical(abs(state_seq_2g[[1]] - state_seq_2g[[2]]))
  len <- length(is_vml)

  i <- j <- 1; start_ind <- end_ind <- NULL
  while (j <= len) {
    if (i==j) { # either i,j are non-VML or start of a VMR
      if (is_vml[i]) {
        start_ind <- c(start_ind, i)
        if (j == len) end_ind <- c(end_ind, j)
        j <- j + 1
      } else {
        i <- j <- j + 1
      }
    } else { # i at start of VMR, j goes to end of VMR
      if (is_vml[j]) {
        if (j == len) end_ind <- c(end_ind, j)
        j <- j + 1
      } else {
        if (j + max_n_merge <= len & is_vml[j + max_n_merge]) {
          j <- j + max_n_merge
        } else {
          end_ind <- c(end_ind, j-1)
          i <- j <- j + max_n_merge + 1
        }
      }
      # if (j == len) end_ind <- c(end_ind, j)
    }
  }

  if (length(start_ind) != length(end_ind))
    stop("Unmatched length of start_ind and end_ind.")
  inds <- data.frame(start_ind = start_ind, end_ind = end_ind)
  if (nrow(inds) > 0) inds <- inds %>% filter(end_ind - start_ind + 1 >= min_n)
  if (nrow(inds) == 0) inds <- NULL
  return(inds)
}
