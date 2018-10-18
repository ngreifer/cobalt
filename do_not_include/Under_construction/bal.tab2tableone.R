bal.tab2tableone <- function(b, un = TRUE) {
    co.names <- b$print.options$co.names
    
    vars <- unique(sapply(co.names, function(x) x$component[1]))
    varNumerics <- vars[sapply(vars, function(x) (x %in% rownames(b$Balance)) && (b$Balance[x, "Type"] != "Binary"))]
    varFactors <- setdiff(vars, varNumerics)
    
    logiFactors <- vars %in% cat.vars
    
    missing.indicators <- rownames(b$Balance)[sapply(co.names, function(co) "na" %in% co$type)]
    vars.w.missing <- sapply(co.names[missing.indicators], function(co) paste0(co$component[co$type != "na"], collapse = ""))
    names(missing.indicators) <- vars.w.missing
    
    varLabels <- setNames(vector("list", length(vars)), vars)
    
    missing.count <- setNames(rep(NA, length(vars)), vars)
    
    ##Create ContTable
    result <- lapply(names(b$print.options$treat.names), function(x) {
        mat <- matrix(NA_real_, nrow = length(varNumerics), ncol = 12,
                      dimnames = list(varNumerics,
                                      c("n", "miss", "p.miss", "mean", "sd", "median", "p25", "p75", 
                                        "min", "max", "skew", "kurt")))
        for (i in varNumerics) {
            mat[i, "n"] <- b$Observations["All", b$print.options$treat.names[x]]
            mat[i, "p.miss"] <- 100*ifelse(i %in% names(missing.indicators), b$Balance[missing.indicators[i], paste0("M.", x, ".Un")], 0)
            mat[i, "miss"] <- round(mat[i, "p.miss"]*mat[i, "n"]/100)
            mat[i, "mean"] <- b$Balance[i, paste0("M.", x, ".Un")]
            
        }
        mat
    })
    names(result) <- b$print.options$treat.names
    
    if (length(result) > 1 ) {
        ## Add an attribute for the stratifying variable name
        attributes(result) <- c(attributes(result),
                                list(strataVarName = "treat"))
    }
    
    ## list of smds
    smds <- matrix(abs(b$Balance[varNumerics, "Diff.Un"]), ncol = 1, 
                   dimnames = list(varNumerics, paste0(names(b$print.options$treat.names), collapse = " vs ")))
    
    ## Return object
    ## Give an S3 class
    class(result) <- c("ContTable", class(result))
    
    ## Give additional attributes
    attributes(result) <- c(attributes(result),
                            list(smd     = smds))
    
    ContTable <- result
    
    
    ##Create CatTable
    result <- lapply(names(b$print.options$treat.names), function(x) {
        out <- setNames(lapply(varFactors, function(v) {
            co.names.ind.with.v <- sapply(co.names, function(co) {
                v == co$component[1] && (length(co$component) == 1 || (co$type[2] == "fsep" && co$type[3] == "level"))
            })
            num.levels <- sum(co.names.ind.with.v)
            if (num.levels > 1) {
                v.levels <- sapply(co.names[co.names.ind.with.v],
                                   function(co) ifelse(length(co$component) == 1, "1", co$component[co$type == "level"][1]))
                d <- data.frame(n = rep(b$Observations["All", b$print.options$treat.names[x]], length(v.levels)),
                                miss = NA,
                                p.miss = 100*rep(ifelse(v %in% names(missing.indicators), b$Balance[missing.indicators[v], paste0("M.", x, ".Un")], 0), length(v.levels)),
                                level = v.levels,
                                freq = NA,
                                percent = 100*sapply(v.levels, function(lev) {
                                    b$Balance[paste0(v, attr(co.names, "sep")["factor"], lev, sep = ""), paste0("M.", x, ".Un")]
                                }),
                                cum.percent = NA,
                                row.names = seq_along(v.levels))
                d[["miss"]] <- round(d[["p.miss"]]*d[["n"]]/100)
                d[["freq"]] <- round((d[["percent"]])*(d[["n"]]-d[["miss"]])/100)
                d[["cum.percent"]] <- cumsum(d[["percent"]])
            }
            else {
                has.level <- length(co.names[co.names.ind.with.v][[1]]$component) > 1
                v.levels <- c("0", ifelse(has.level, co.names[co.names.ind.with.v][[1]]$component[3], ""))
                d <- data.frame(n = rep(b$Observations["All", b$print.options$treat.names[x]], 2),
                                miss = NA,
                                p.miss = 100*rep(ifelse(v %in% names(missing.indicators), b$Balance[missing.indicators[v], paste0("M.", x, ".Un")], 0), 2),
                                level = ifelse(v.levels == "", "1", v.levels),
                                freq = NA,
                                percent = c(NA, 100*b$Balance[paste0(c(v, v.levels[2][v.levels[2] != ""]), collapse = attr(co.names, "sep")["factor"]), paste0("M.", x, ".Un")]),
                                cum.percent = NA,
                                row.names = seq_along(v.levels))
                d[["miss"]] <- round(d[["p.miss"]]*d[["n"]]/100)
                d[["percent"]][1] <- 100 - d[["percent"]][2]
                d[["freq"]] <- round((d[["percent"]])*(d[["n"]]-d[["miss"]])/100)
                d[["cum.percent"]] <- cumsum(d[["percent"]])
            }
            d
        }), varFactors)
    })
    names(result) <- b$print.options$treat.names
    
    if (length(result) > 1 ) {
        ## Add an attribute for the stratifying variable name
        attributes(result) <- c(attributes(result),
                                list(strataVarName = "treat"))
    }
    
    ## list of smds
    smds <- matrix(abs(b$Balance[varFactors, "Diff.Un"]), ncol = 1, 
                   dimnames = list(varFactors, paste0(names(b$print.options$treat.names), collapse = " vs ")))
    
    class(result) <- c("CatTable", class(result))
    
    ## Give additional attributes
    attributes(result) <- c(attributes(result),
                            list(smd     = smds))
    
    CatTable <- result
    
    ## Compute total missingness
    
    missing.count[varNumerics] <- sapply(varNumerics, function(v) {
        sum(sapply(ContTable, function(t) t[v, "miss"]))
    })
    missing.count[varFactors] <- sapply(varFactors, function(v) {
        sum(sapply(CatTable, function(t) t[[v]][1, "miss"]))
    })
    percentMissing <- 100*missing.count/sum(b$Observations["All", ])
    
    ## Create a list for output
    ## Either one of the two tables may be NULL
    TableOneObject <- list(ContTable = ContTable,
                           CatTable  = CatTable,
                           MetaData  = list(vars        = vars,
                                            ## describes which pos is vars is factor
                                            logiFactors = logiFactors,
                                            ## names of vars of each type
                                            varFactors  = varFactors,
                                            varNumerics = varNumerics,
                                            ## Missing data percentage for each variable (no strata).
                                            percentMissing = percentMissing,
                                            ## Variable labels
                                            varLabels = varLabels))
    
    ## Give a class
    class(TableOneObject) <- "TableOne"
    
    ## Return the object
    return(TableOneObject)
}
