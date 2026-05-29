skew.diff <- function(x, group, weights, var.type) {
    #Calculates weighted skew and skew difference for groups; uses formulas at http://www.nematrian.com/WeightedMomentsAndCumulants
    skew.cols <- data.frame(Skew.C = NA, Skew.T = NA, Skew.Diff = NA)
    if (is_null(weights)) weights <- rep(1, length(x))
    if (var.type=="Contin.") {
        # skew for each group
        x0 <- x[group==0];       x1 <- x[group==1];
        w0 <- weights[group==0]; w1 <- weights[group==1]
        m0 <- w.m(x0, w0);       m1 <- w.m(x1, w1)
        wvp <- function(x, w, m) return(sum(w*(x-m)^2, na.rm=TRUE)/sum(w, na.rm=TRUE)) #weighted variance population
        Fv <- function(w) return((sum(w0, na.rm=TRUE)^2-sum(w0^2, na.rm=TRUE))/(sum(w0, na.rm=TRUE)^2)) #adjustment for sample variance
        wskp <- function(x, w, m, sp) return(sum(w*((x-m)/sp)^3, na.rm=TRUE)/sum(w, na.rm=TRUE)) #weighted skew population
        Fsk <- function(w, sp, s) return(((sum(w0, na.rm=TRUE)^3 - 3*sum(w0^2, na.rm=TRUE)*sum(w0, na.rm=TRUE) + 2*sum(w0^3, na.rm=TRUE))/(sum(w0, na.rm=TRUE)^3))*(sp/s)^3 ) #adjustment for sample skew
        #group==0
        sp0 <- sqrt(wvp(x0, w0, m0))
        s0 <- sqrt(wvp(x0, w0, m0)/Fv(w0))
        skp0 <- wskp(x0, w0, m0, sp0)
        sk0 <- skp0/Fsk(w0, sp0, s0)
        #group==1
        sp1 <- sqrt(wvp(x1, w1, m1))
        s1 <- sqrt(wvp(x1, w1, m1)/Fv(w1))
        skp1 <- wskp(x1, w1, m1, sp1)
        sk1 <- skp1/Fsk(w1, sp1, s1)
        
        skew.cols[, ] <- c(sk0, sk1, sk1-sk0)
    }
    return(skew.cols)
}

col_w_edist <- function(mat, treat, weights = NULL, s.weights = NULL, bin.vars, subset = NULL, na.rm = TRUE, ...) {
    needs.splitting <- FALSE
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(to.split <- vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("'mat' must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("'mat' must be a data.frame or numeric matrix.")
    
    bin.vars <- process.bin.vars(bin.vars, mat)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                   names(A) %nin% c("data", "var.name", "drop.first",
                                    "drop.level", "split.with")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2",
                                           split.with = bin.vars),
                                      A))
        bin.vars <- attr(mat, "split.with")[[1]]
    }
    
    if (!is.atomic(treat) || !is_binary(treat)) stop("treat must be a binary variable.")
    
    possibly.supplied <- c("mat", "treat", "weights", "s.weights", "subset")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
    
    if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
    else if (anyNA(as.logical(subset))) stop("'subset' must be a logical vector.")
    
    weights <- weights * s.weights
    
    weights <- weights[subset]
    treat <- treat[subset]
    mat <- mat[subset, , drop = FALSE]
    
    tval1 <- treat[1]
    edist <- rep(NA_real_, NCOL(mat))
    
    if (any(!bin.vars)) {
        mat[,!bin.vars] <- mat[,!bin.vars, drop = FALSE]/sqrt(col.w.v(mat[,!bin.vars]))
        
        t1 <- treat == tval1
        
        weights[t1] <- weights[t1]/mean(weights[t1])
        weights[!t1] <- weights[!t1]/mean(weights[!t1])
        
        J <- diag(t1/sum(t1) - (!t1)/sum(!t1))
        
        edist[!bin.vars] <- vapply(which(!bin.vars), function(i) {
            x <- mat[,i]
            d <- as.matrix(dist(x))
            sqrt(drop(-t(weights) %*% J %*% d %*% J %*% weights)/2)
        }, numeric(1L))
    }
    if (any(bin.vars)) {
        edist[bin.vars] <- (abs(col.w.m(mat[treat == tval1, bin.vars, drop = FALSE], weights[treat == tval1], na.rm = na.rm) - 
                                    col.w.m(mat[treat != tval1, bin.vars, drop = FALSE], weights[treat != tval1], na.rm = na.rm)))
    }
    
    setNames(edist, colnames(mat))
    
}
col_pair_diff <- function(mat, treat, strata = NULL, std = TRUE, s.d.denom = "pooled", bin.vars, subset = NULL, ...) {
    needs.splitting <- FALSE
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(to.split <- vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else mat <- as.matrix(mat)
        }
        else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
        else stop("mat must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
    
    bin.vars <- process.bin.vars(bin.vars, mat)
    
    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                   names(A) %nin% c("data", "var.name", "drop.first",
                                    "drop.level", "split.with")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2",
                                           split.with = bin.vars),
                                      A))
        bin.vars <- attr(mat, "split.with")[[1]]
    }
    
    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector or factor.")
    if (!is.atomic(std) || anyNA(as.logical(std)) ||
        length(std) %nin% c(1L, NCOL(mat))) {
        stop("std must be a logical vector with length equal to 1 or the number of columns of mat.")
    }
    if (is_not_null(subset)) stop("subset cannot be non-NULL.", call. = FALSE)
    
    possibly.supplied <- c("mat", "treat", "subset")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }
    
    if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
    else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
    
    if (!is_binary(treat[subset])) stop("treat must be a binary variable.")
    
    if (length(std) == 1L) std <- rep(std, NCOL(mat))
    
    tval1_0 <- treat[1]
    
    denoms <- compute_s.d.denom(mat, treat = treat, 
                                s.d.denom = s.d.denom, 
                                bin.vars = bin.vars, to.sd = std)
    std_mat <- mat_div(mat, denoms)
    
    if (is_not_null(strata)) {
        out <- apply(std_mat, 2, function(x) {
            diffs <- sapply(unique(strata), function(s) {
                mean_fast(x[strata == s & treat == tval1_0]) - mean_fast(x[strata == s & treat != tval1_0])
            })
            return(sqrt(mean_fast(diffs^2)))
        })
    }
    else {
        out <- sqrt((col_w_sd(mat, bin.vars = bin.vars, subset = treat == tval1_0)/denoms)^2 +
                        (col_w_sd(mat, bin.vars = bin.vars, subset = treat != tval1_0)/denoms)^2 +
                        (col_w_smd(mat, treat, std = FALSE, bin.vars = bin.vars)/denoms)^2)
    }
    out
    
}

STATS[["pair.diffs"]] <- {list(
    type = "bin",
    threshold = "pair.threshold",
    Threshold = "Pair.Threshold",
    disp_stat = "disp.pair.diff",
    bin_only = FALSE,
    abs = function(x) abs_(x),
    bal.tab_column_prefix = "Pair.Diff", #Also which.stat in love.plot
    threshold_range = c(0, Inf),
    balance_tally_for = "within-pair mean differences",
    variable_with_the_greatest = "within-pair mean difference", #also which.stat2 in love.plot
    love.plot_xlab = function(...) {
        A <- list(...)
        binary <- A$binary #attr(x, "print.options")$binary
        continuous <- A$continuous #attr(x, "print.options")$continuous
        abs <- A$abs
        var_type <- A$var_type #B[["type"]]
        stars <- A$stars
        
        #All std, no std, some std
        if ((binary == "std" || sum(var_type == "Binary") == 0) && 
            (continuous == "std" || sum(var_type != "Binary") == 0)) {
            xlab.diff <- "Standardized Mean Differences"
        } 
        else if ((binary == "raw" || sum(var_type == "Binary") == 0) && 
                 (continuous == "raw" || sum(var_type != "Binary") == 0)) {
            xlab.diff <- "Mean Differences"
        }
        else {
            stars <- match_arg(stars, c("none", "std", "raw"))
            if (stars == "none") {
                xlab.diff <- "Mean Differences"
            }
            else if (stars == "std") {
                xlab.diff <- "Mean Differences"
            }
            else if (stars == "raw") {
                xlab.diff <- "Standardized Mean Differences"
            }
        }
        
        xlab <- if (abs) paste("Absolute", xlab.diff) else xlab.diff
        return(xlab)
    },
    love.plot_add_stars = function(SS.var, variable.names, ...) {
        A <- list(...)
        binary <- A$binary #attr(x, "print.options")$binary
        continuous <- A$continuous #attr(x, "print.options")$continuous
        var_type <- A$var_type #B[["Type"]]
        stars <- A$stars
        star_char = A$star_char #args$star_char
        
        #All std, no std, some std
        if (!((binary == "std" || sum(var_type == "Binary") == 0) && 
              (continuous == "std" || sum(var_type != "Binary") == 0)) 
            &&
            !((binary == "raw" || sum(var_type == "Binary") == 0) && 
              (continuous == "raw" || sum(var_type != "Binary") == 0))) {
            
            stars <- match_arg(stars, c("none", "std", "raw"))
            if (stars == "none") {
                warning("Standardized mean differences and raw mean differences are present in the same plot. \nUse the 'stars' argument to distinguish between them and appropriately label the x-axis.", call. = FALSE)
            }
            else {
                if (length(star_char) != 1 || !is.character(star_char)) star_char <- "*"
                
                vars_to_star <- setNames(rep(FALSE, length(variable.names)), variable.names)
                if (stars == "std") {
                    if (binary == "std") vars_to_star[variable.names[var_type == "Binary"]] <- TRUE
                    if (continuous == "std") vars_to_star[variable.names[var_type != "Binary"]] <- TRUE
                }
                else if (stars == "raw") {
                    if (binary == "raw") vars_to_star[variable.names[var_type == "Binary"]] <- TRUE
                    if (continuous == "raw") vars_to_star[variable.names[var_type != "Binary"]] <- TRUE
                }
                new.variable.names <- setNames(variable.names, variable.names)
                names(new.variable.names)[vars_to_star[variable.names]] <- paste0(variable.names[vars_to_star[variable.names]], star_char)
                SS.var <- do.call(f.recode, c(list(SS.var), new.variable.names))
            }
        }
        
        return(SS.var)
    },
    baseline.xintercept = 0,
    threshold.xintercepts = function(threshold, abs) {
        if (abs) c(lower = base::abs(threshold))
        else c(lower = -base::abs(threshold), upper = base::abs(threshold))
    },
    love.plot_axis_scale = ggplot2::scale_x_continuous,
    fun = function(C, treat, weights, std, s.d.denom, abs, s.weights, bin.vars, weighted.weights = weights, subset = NULL, ...) {
        if (is_null(strata <- attr(weights, "match.strata"))) {
            stop("Within-pair differences cannot be used unless matching strata were specified.", call. = FALSE)
        }
        
        
        col_w_smd(C, treat = treat, weights = weights, 
                  std = std, s.d.denom = s.d.denom,
                  abs = abs, s.weights = s.weights, bin.vars = bin.vars,
                  weighted.weights = weighted.weights,
                  subset = NULL)
    }
)}