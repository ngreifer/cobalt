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

reconstruct_factors <- function(df) {
    bin <- apply(df, 2, is_binary)
    
    factor.groups <- rep(NA_integer_, sum(bin))
    for (i in 1:sum(bin)) {
        if (i == 1) factor.groups <- 1L
        else {
            found <- FALSE
            k <- 1
            while (k < i && !found) {
                if (all(check_if_zero(df[,which(bin)[i]] * df[,which(bin)[k]]))) {
                    found <- TRUE
                    factor.groups[i] <- factor.groups[k]
                }
                else {
                    k <- k + 1
                }
            }
            if (!found) factor.groups[i] <- factor.groups[i-1] + 1
        }
    }
    
    f.list <- lapply(unique(factor.groups), function(f) {
        g.names <- colnames(df)[which(bin)[factor.groups==f]]
        #Find common name
        sn <- strsplit(g.names, "")
        
        stop <- FALSE
        k <- 1
        while (k <= min(nchar(g.names)) && !stop) {
            if (nunique.gt(sapply(sn, function(x) paste0(x[1:k], collapse = "")), 1)) {
                end.same <- k - 1
                stop <- TRUE
            }
            else k <- k + 1
        }
        if (end.same == -1) {
            #No sameness in names
        }
        else {
            symbols <- c("_", ".", " ", ":")
            if (any(symbols %in% sn[[1]][1:end.same])) {
                sep <- symbols[symbols == sn[[1]][1:end.same][last(which(sn[[1]][1:end.same] %in% symbols))]]
                last_sep_index <- last(which(sn[[1]][1:end.same] == sep))
                f.name <- paste0(sn[[1]][1:(last_sep_index - 1)], collapse = "")
                f.levels <- sapply(g.names, function(x) substr(x, last_sep_index + 1, nchar(x)))
            }
            else {
                sep <- ""
                f.name <- paste0(sn[[1]][1:end.same], collapse = "")
                f.levels <- sapply(g.names, function(x) substr(x, nchar(f.name) + 1, nchar(x)))
            }
            if (any(check_if_zero(rowSums(df[,g.names, drop = FALSE])))) {
                other <- "other"
                while (other %in% f.levels) other <- paste0(other, "_")
            }
            else other <- NULL
        }
        return(c(f.name = f.name, sep = sep, other = other))
    })
    for (f in f.list) df <- unsplitfactor(df, var.name = f["f.name"], sep = f["sep"], dropped.level = f["other"])
    return(df)
}
