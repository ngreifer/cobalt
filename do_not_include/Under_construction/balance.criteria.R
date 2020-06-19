#Criterion
#Criteria to use in methods that involve a tuning parameter. Also called "stop.method" in twang.
#Separate functions for different treatment types.

#init functions
init_es <- function(covs, treat, s.weights = NULL, estimand = "ATE") {
    needs.splitting <- FALSE
    if (!is.matrix(covs)) {
        if (is.data.frame(covs)) {
            if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else covs <- as.matrix(covs)
        }
        else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
        else stop("covs must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(covs)) stop("covs must be a data.frame or numeric matrix.")
    
    bin.vars <- process.bin.vars(mat = covs)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        covs <- do.call("splitfactor", list(covs, drop.first ="if2",
                                            split.with = bin.vars))
        bin.vars <- attr(covs, "split.with")[[1]]
    }
    
    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector or factor.")
    
    possibly.supplied <- c("covs", "treat", "s.weights")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(covs))
    
    if (!is_binary(treat)) stop("treat must be a binary variable.")
    
    s.d.denom <- get.s.d.denom(estimand = estimand, treat = treat, quietly = TRUE)
    
    denoms <- compute_s.d.denom(covs, treat = treat, 
                                s.d.denom = s.d.denom, s.weights = s.weights, 
                                bin.vars = bin.vars, to.sd = rep(TRUE, ncol(covs)))
    
    list(treat = treat,
         covs = covs,
         bin.vars = bin.vars,
         s.weights = s.weights,
         s.d.denom = denoms)
}
init_ks <- function(covs, treat, s.weights = NULL) {
    needs.splitting <- FALSE
    if (!is.matrix(covs)) {
        if (is.data.frame(covs)) {
            if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else covs <- as.matrix(covs)
        }
        else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
        else stop("covs must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(covs)) stop("covs must be a data.frame or numeric matrix.")
    
    bin.vars <- process.bin.vars(mat = covs)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        covs <- do.call("splitfactor", list(covs, drop.first ="if2",
                                            split.with = bin.vars))
        bin.vars <- attr(covs, "split.with")[[1]]
    }
    
    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector or factor.")
    
    possibly.supplied <- c("covs", "treat", "s.weights")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(covs))
    
    if (!is_binary(treat)) stop("treat must be a binary variable.")
    
    list(treat = treat,
         covs = covs,
         bin.vars = bin.vars,
         s.weights = s.weights)
}
init_mahalanobis <- function(covs, treat, s.weights = NULL, estimand = "ATE") {
    needs.splitting <- FALSE
    if (!is.matrix(covs)) {
        if (is.data.frame(covs)) {
            if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else covs <- as.matrix(covs)
        }
        else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
        else stop("covs must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(covs)) stop("covs must be a data.frame or numeric matrix.")
    
    bin.vars <- process.bin.vars(mat = covs)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        covs <- do.call("splitfactor", list(covs, drop.first ="if2",
                                            split.with = bin.vars))
        bin.vars <- attr(covs, "split.with")[[1]]
    }
    
    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector or factor.")
    
    possibly.supplied <- c("covs", "treat", "s.weights")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(covs))
    
    if (!is_binary(treat)) stop("treat must be a binary variable.")
    
    s.d.denom <- get.s.d.denom(estimand = estimand, treat = treat, quietly = TRUE)
    
    if (s.d.denom %in% as.character(treat)) {
        sigma <- cov.wt(covs[treat == s.d.denom,,drop = FALSE], s.weights[treat == s.d.denom])$cov
        if (any(zeros <- diag(sigma) == 0)) {
            sigma_all <- cov.wt(covs, s.weights)$cov
            sigma[zeros,] <- sigma_all[zeros,]
            sigma[,zeros] <- sigma_all[,zeros]
            if (any(zeros <- diag(sigma) == 0)) diag(sigma)[zeros]  <- 1
        }
    }
    else {
        sigma <- cov.wt(covs, s.weights)$cov
        if (any(zeros <- diag(sigma) == 0)) diag(sigma)[zeros]  <- 1
    }
    sigma_inv <- chol2inv(chol(sigma))
    
    list(treat = treat,
         covs = covs,
         bin.vars = bin.vars,
         s.weights = s.weights,
         s.d.denom = s.d.denom,
         sigma_inv = sigma_inv)
}
init_energy.dist <- function(covs, treat, s.weights = NULL, estimand = "ATE", focal = NULL, improved = TRUE) {
    needs.splitting <- FALSE
    if (!is.matrix(covs)) {
        if (is.data.frame(covs)) {
            if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else covs <- as.matrix(covs)
        }
        else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
        else stop("covs must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(covs)) stop("covs must be a data.frame or numeric matrix.")
    
    if (needs.splitting) {
        covs <- do.call("splitfactor", list(covs, drop.first ="if2"))
    }
    
    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector or factor.")
    
    possibly.supplied <- c("covs", "treat", "s.weights")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(covs))
    
    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)
    
    if (treat.type == "continuous") stop("treat must be a binary or multinomial variable.")
    
    f.e.r <- process.focal.and.estimand(focal, estimand, treat, treat.type)
    focal <- f.e.r[["focal"]]
    estimand <- f.e.r[["estimand"]]

    covs <- mat_div(center(covs, at = col.w.m(covs, s.weights)),
                    sqrt(col.w.v(covs, s.weights)))
    
    d <- as.matrix(dist(covs))
    
    n <- length(treat)
    levels_treat <- unique(treat)
    diagn <- diag(n)
    
    for (t in levels_treat) s.weights[treat == t] <- s.weights[treat == t]/mean_fast(s.weights[treat == t])
    
    tmat <- vapply(levels_treat, function(t) treat == t, logical(n))
    nt <- colSums(tmat)
    
    J <- setNames(lapply(levels_treat, function(t) s.weights*tmat[,t]/nt[t]), levels_treat)
    
    if (is_null(focal)) {
        J0 <- as.matrix(s.weights/n)
        
        M2_array <- vapply(levels_treat, function(t) -2 * tcrossprod(J[[t]]) * d, diagn)
        M1_array <- vapply(levels_treat, function(t) 2 * J[[t]] * d %*% J0, J0)
        
        M2 <- rowSums(M2_array, dims = 2)
        M1 <- rowSums(M1_array)
        
        if (improved) {
            all_pairs <- combn(levels_treat, 2, simplify = FALSE)
            M2_pairs_array <- vapply(all_pairs, function(p) -tcrossprod(J[[p[1]]]-J[[p[2]]]) * d, diagn)
            M2 <- M2 + rowSums(M2_pairs_array, dims = 2)
        }
    }
    else {
        J0_focal <- as.matrix(J[[focal]])
        clevs <- levels_treat[levels_treat != focal]
        
        M2_array <- vapply(clevs, function(t) -2 * tcrossprod(J[[t]]) * d, diagn)
        M1_array <- vapply(clevs, function(t) 2 * J[[t]] * d %*% J0_focal, J0_focal)
        
        M2 <- rowSums(M2_array, dims = 2)
        M1 <- rowSums(M1_array)
    }
    
    list(M1 = M1,
         M2 = M2,
         s.weights = s.weights)
}
init_p <- function(covs, treat, s.weights = NULL) {
    needs.splitting <- FALSE
    if (!is.matrix(covs)) {
        if (is.data.frame(covs)) {
            if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else covs <- as.matrix(covs)
        }
        else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
        else stop("covs must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(covs)) stop("covs must be a data.frame or numeric matrix.")
    
    bin.vars <- process.bin.vars(mat = covs)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        covs <- do.call("splitfactor", list(covs, drop.first ="if2",
                                            split.with = bin.vars))
        bin.vars <- attr(covs, "split.with")[[1]]
    }
    
    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector or factor.")
    
    possibly.supplied <- c("covs", "treat", "s.weights")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(covs))
    
    if (get.treat.type(assign.treat.type(treat)) != "continuous") stop("treat must be a numeric non-binary variable.")
    
    s.d.denom <- get.s.d.denom.cont(quietly = TRUE)
    
    denoms <- compute_s.d.denom(covs, treat = treat, 
                                s.d.denom = s.d.denom, s.weights = s.weights, 
                                bin.vars = bin.vars, to.sd = rep(TRUE, ncol(covs)))

    list(treat = treat,
         covs = covs,
         bin.vars = bin.vars,
         s.weights = s.weights,
         s.d.denom = denoms)
}
init_s <- function(covs, treat, s.weights = NULL) {
    needs.splitting <- FALSE
    if (!is.matrix(covs)) {
        if (is.data.frame(covs)) {
            if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else covs <- as.matrix(covs)
        }
        else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
        else stop("covs must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(covs)) stop("covs must be a data.frame or numeric matrix.")
    
    bin.vars <- process.bin.vars(mat = covs)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        covs <- do.call("splitfactor", list(covs, drop.first ="if2",
                                            split.with = bin.vars))
        bin.vars <- attr(covs, "split.with")[[1]]
    }
    
    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector or factor.")
    
    possibly.supplied <- c("covs", "treat", "s.weights")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(covs))
    
    if (get.treat.type(assign.treat.type(treat)) != "continuous") stop("treat must be a numeric non-binary variable.")
    
    for (i in 1:ncol(covs)) if (!bin.vars[i]) covs[,i] <- rank(covs[,i], na.last = "keep")
    treat <- rank(treat, na.last = "keep")
    
    s.d.denom <- get.s.d.denom.cont(quietly = TRUE)
    
    denoms <- compute_s.d.denom(covs, treat = treat, 
                                s.d.denom = s.d.denom, s.weights = s.weights, 
                                bin.vars = bin.vars, to.sd = rep(TRUE, ncol(covs)))
    
    list(treat = treat,
         covs = covs,
         bin.vars = bin.vars,
         s.weights = s.weights,
         s.d.denom = denoms)
}

bal_criterion <- function(treat.type = "binary", criterion, list = FALSE) {
    treat.type <- match_arg(treat.type, c("binary", "multinomial", "continuous"))
    criteria <- switch(treat.type,
           binary = list(
               es.mean = list(
                   fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", init = NULL, ...) {
                       if (is_null(init)) {
                           init <- init_es(covs, treat, s.weights, estimand)
                       }
                       else {
                           if (any(c("covs", "treat", "s.weights", "bin.vars", "s.d.denom") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       }
                       
                       es <- col_w_smd(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights, 
                                               bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE)
                       mean_fast(es)
                   },
                   init = init_es
               ),
               es.max = list(
                   fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", init = NULL, ...) {
                       if (is_null(init)) {
                           init <- init_es(covs, treat, s.weights, estimand)
                       }
                       else {
                           if (any(c("covs", "treat", "s.weights", "bin.vars", "s.d.denom") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       }
                       
                       es <- col_w_smd(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights, 
                                               bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE)
                       max(es)
                   },
                   init = init_es
               ),
               es.rms = list(
                   fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", init = NULL, ...) {
                       if (is_null(init)) {
                           init <- init_es(covs, treat, s.weights, estimand)
                       }
                       else {
                           if (any(c("covs", "treat", "s.weights", "bin.vars", "s.d.denom") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       }
                       
                       es <- col_w_smd(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights, 
                                               bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE)
                       rms(es)
                   },
                   init = init_es
               ),
               ks.mean = list(
                   fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                       if (is_null(init)) {
                           init <- init_ks(covs, treat, s.weights)
                       }
                       else {
                           if (any(c("covs", "treat", "s.weights", "bin.vars") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       }
                       
                       ks <- col_w_ks(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights, 
                                              bin.vars = init$bin.vars)
                       mean_fast(ks)
                   },
                   init = init_ks
               ),
               ks.max = list(
                   fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                       if (is_null(init)) {
                           init <- init_ks(covs, treat, s.weights)
                       }
                       else {
                           if (any(c("covs", "treat", "s.weights", "bin.vars") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       }
                       
                       ks <- col_w_ks(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
                                              bin.vars = init$bin.vars)
                       max(ks)
                   },
                   init = init_ks
               ),
               ks.rms = list(
                   fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                       if (is_null(init)) {
                           init <- init_ks(covs, treat, s.weights)
                       }
                       else {
                           if (any(c("covs", "treat", "s.weights", "bin.vars") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       }
                       
                       ks <- col_w_ks(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
                                              bin.vars = init$bin.vars)
                       rms(ks)
                   },
                   init = init_ks
               ),
               mahalanobis = list(
                   fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", init = NULL, ...) {
                       if (is_null(init)) {
                           init <- init_mahalanobis(covs, treat, s.weights, estimand)
                       }
                       else {
                           if (any(c("covs", "treat", "s.weights", "bin.vars", "sigma_inv") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       } 
                       
                       mean.diffs <- col_w_smd(init$covs, init$treat, weights, s.weights = init$s.weights, 
                                                       bin.vars = init$bin.vars, std = FALSE)
                       drop(t(mean.diffs) %*% init$sigma_inv %*% mean.diffs)
                   },
                   init = init_mahalanobis
               ),
               energy.dist = list(
                   fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = "ATE", focal = NULL, improved = TRUE, init = NULL, ...) {
                       
                       if (is_null(init)) {
                           init <- init_energy.dist(covs, treat, s.weights, estimand, focal, improved)
                       }
                       else {
                           if (any(c("M1", "M2", "s.weights") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       } 
                       
                       M1 <- init[["M1"]]
                       M2 <- init[["M2"]]
                       
                       if (is_null(weights)) weights <- rep(1, nrow(M2))
                       
                       weights <- weights * init[["s.weights"]]
                       
                       return(drop(.5 * t(weights) %*% M2 %*% weights + t(M1) %*% weights))
                   },
                   init = init_energy.dist
               )
           ),
           multinomial =  list(
               es.mean = list(
                   fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", pairwise = TRUE, init = NULL, ...) {
                       if (is_null(init)) {
                           init <- init_es(covs, treat, s.weights, estimand)
                       }
                       else {
                           if (any(c("covs", "treat", "s.weights", "bin.vars", "s.d.denom") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       }
                       
                       if (pairwise) 
                       es <- col_w_smd(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights, 
                                       bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE)
                       mean_fast(es)
                   },
                   init = init_es
               ),
               energy.dist = list(
                   fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = "ATE", focal = NULL, improved = TRUE, init = NULL, ...) {
                       
                       if (is_null(init)) {
                           init <- init_energy.dist(covs, treat, s.weights, estimand, focal, improved)
                       } 
                       else {
                           if (any(c("M1", "M2", "s.weights") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       } 
                       
                       M1 <- init[["M1"]]
                       M2 <- init[["M2"]]
                       
                       if (is_null(weights)) weights <- rep(1, nrow(M2))
                       
                       weights <- weights * init[["s.weights"]]
                       
                       return(drop(.5 * t(weights) %*% M2 %*% weights + t(M1) %*% weights))
                   },
                   init = init_energy.dist
               )
           ),
           
           continuous = list(
               p.mean = list(
                   fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                       if (is_null(init)) {
                           init <- init_p(covs, treat, s.weights)
                       }
                       else {
                           if (any(c("covs", "treat", "s.weights", "bin.vars", "s.d.denom") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       }
                       
                       p <- col_w_cov(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights, 
                                      bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE,
                                      std = TRUE)
                       mean_fast(p)
                   },
                   init = init_p
               ),
               p.max = list(
                   fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                       if (is_null(init)) {
                           init <- init_p(covs, treat, s.weights)
                       }
                       else {
                           if (any(c("covs", "treat", "s.weights", "bin.vars", "s.d.denom") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       }
                       
                       p <- col_w_cov(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights, 
                                      bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE,
                                      std = TRUE)
                       max(p)
                   },
                   init = init_p
               ),
               p.rms = list(
                   fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                       if (is_null(init)) {
                           init <- init_p(covs, treat, s.weights)
                       }
                       else {
                           if (any(c("covs", "treat", "s.weights", "bin.vars", "s.d.denom") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       }
                       
                       p <- col_w_cov(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights, 
                                      bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE,
                                      std = TRUE)
                       rms(p)
                   },
                   init = init_p
               ),
               s.mean = list(
                   fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                       if (is_null(init)) {
                           init <- init_s(covs, treat, s.weights)
                       }
                       else {
                           if (any(c("covs", "treat", "s.weights", "bin.vars", "s.d.denom") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       }
                       
                       s <- col_w_cov(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights, 
                                      bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE,
                                      std = TRUE)
                       mean_fast(s)
                   },
                   init = init_s
               ),
               s.max = list(
                   fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                       if (is_null(init)) {
                           init <- init_s(covs, treat, s.weights)
                       }
                       else {
                           if (any(c("covs", "treat", "s.weights", "bin.vars", "s.d.denom") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       }
                       
                       s <- col_w_cov(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights, 
                                      bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE,
                                      std = TRUE)
                       max(s)
                   },
                   init = init_s
               ),
               s.rms = list(
                   fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                       if (is_null(init)) {
                           init <- init_s(covs, treat, s.weights)
                       }
                       else {
                           if (any(c("covs", "treat", "s.weights", "bin.vars", "s.d.denom") %nin% names(init))) {
                               stop("init was not correctly supplied.")
                           }
                       }
                       
                       s <- col_w_cov(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights, 
                                      bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE,
                                      std = TRUE)
                       rms(s)
                   },
                   init = init_s
               )
           )
    )
    if (isTRUE(list)) return(names(criteria))
    if (missing(criterion)) stop("'criterion' must be provided.", call. = FALSE)
    criterion <- match_arg(criterion, names(criteria))
    
    out <- criteria[[criterion]]
    attr(out, "treat.type") <- treat.type
    attr(out, "criterion") <- criterion
    class(out) <- "bal_criterion"
    out
}

print.bal_criterion <- function(x, ...) {
    cat("A bal_criterion object\n")
    cat(paste0("  treatment type: ", attr(x, "treat.type"), "\n"))
    cat(paste0("  criterion: ", attr(x, "criterion"), " (", bal_criterion.to.phrase(attr(x, "criterion")), ")\n"))
    cat("With components:\n - $init (for initialization)\n - $fun (to compute the criterion value)\n")
    invisible(x)
}
bal_criterion.to.phrase <- function(criterion) {
    
    switch(criterion,
           "es.mean" = "average absolute standardized mean difference",
           "es.max" = "maximum absolute standardized mean difference",
           "es.rms" = "root-mean-square absolute standardized mean difference",
           "ks.mean" = "average Kolmogorov–Smirnov statistic",
           "ks.max" = "maximum Kolmogorov–Smirnov statistic",
           "ks.rms" = "root-mean-square Kolmogorov-Smirnov statistic",
           "mahalanobis" = "sample mahalanobis distance",
           "energy.dist" = "energy distance",
           stop(paste0("\"", criterion, "\" is not an allowed criterion."))
    )

}
