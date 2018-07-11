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

