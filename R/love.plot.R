love.plot <- function(x, stats, abs, agg.fun = NULL, 
                      var.order = NULL, drop.missing = TRUE, drop.distance = FALSE, 
                      thresholds = NULL, line = FALSE, stars = "none", grid = FALSE, 
                      limits = NULL, colors = NULL, shapes = NULL, alpha = 1, size = 3, 
                      wrap = 30, var.names = NULL, title, sample.names, labels = FALSE,
                      position = "right", themes = NULL, ...) {
    
    #Replace .all and .none with NULL and NA respectively
    .call <- match.call()
    .alls <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.all)), logical(1L))
    .nones <- vapply(seq_along(.call), function(z) identical(.call[[z]], quote(.none)), logical(1L))
    if (any(c(.alls, .nones))) {
        .call[.alls] <- expression(NULL)
        .call[.nones] <- expression(NA)
        return(eval.parent(.call))
    }
    
    if (missing(stats)) stats <- NULL
    
    #Re-call bal.tab with disp.v.ratio or disp.ks if stats = "v" or "k".
    if (typeof(.call[["x"]]) == "language") { #if x is not an object (i.e., is a function call)
        
        replace.args <- function(m) {
            #m is bal.tab call or list (for do.call)
            m[["un"]] <- TRUE
            if (is_not_null(stats)) m[["stats"]] <- stats
            
            if (any(names(m) == "agg.fun")) m[["agg.fun"]] <- NULL
            
            if (any(names(m) %pin% "abs")) m[["abs"]] <- abs
            
            if (any(names(m) %pin% "thresholds")) m["thresholds"] <- list(NULL)
            
            return(m)
        }
        
        if (deparse1(.call[["x"]][[1]]) %in% c("bal.tab", methods("bal.tab"))) { #if x i bal.tab call
            .call[["x"]] <- replace.args(.call[["x"]])
            x <- eval.parent(.call[["x"]])
            
        }
        else if (deparse1(.call[["x"]][[1]]) == "do.call") { #if x is do.call
            d <- match.call(eval(.call[["x"]][[1]]), .call[["x"]])
            if (deparse1(d[["what"]]) %in% c("bal.tab", methods("bal.tab"))) {
                d[["args"]] <- replace.args(d[["args"]])
                x <- eval.parent(d)
            }
        }
    }
    
    tryCatch(force(x), error = function(e) stop(conditionMessage(e), call. = FALSE))
    
    if ("bal.tab" %nin% class(x)) {
        #Use bal.tab on inputs first, then love.plot on that
        .call2 <- .call
        .call2[[1]] <- quote(bal.tab)
        .call2[["x"]] <- x
        
        .call2["thresholds"] <- list(NULL)
        
        .call[["x"]] <- .call2
        
        return(eval.parent(.call))
    }
    
    args <- list(...)
    
    #shape (deprecated)
    #un.color (deprecated)
    #adj.color (deprecated)
    #cluster.fun (deprecated)
    #star_char
    
    p.ops <- c("which.cluster", "which.imp", "which.treat", "which.time", "disp.subclass")
    for (i in p.ops) {
        if (rlang::has_name(args, i)) attr(x, "print.options")[[i]] <- args[[i]]
    }
    
    #Using old argument names
    if (is_not_null(args$cluster.fun) && is_null(agg.fun)) agg.fun <- args$cluster.fun
    if (is_not_null(args$no.missing)) drop.missing <- args$no.missing
    
    Agg.Fun <- NULL
    subtitle <- NULL
    
    #Process abs
    if (missing(abs)) {
        abs <- if_null_then(attr(x, "print.options")[["abs"]], TRUE)
    }
    
    #Process stats
    if (is_null(stats)) stats <- attr(x, "print.options")$stats
    stats <- match_arg(stats, all_STATS(attr(x, "print.options")$type), several.ok = TRUE)
    
    #Get B and config
    if ("bal.tab.subclass" %in% class(x)) {
        subclass.names <- names(x[["Subclass.Balance"]])
        sub.B <- do.call("cbind", lapply(subclass.names, function(s) {
            sub <- x[["Subclass.Balance"]][[s]]
            sub.B0 <- setNames(sub[endsWith(names(sub), ".Adj")],
                               gsub(".Adj", paste0(".Subclass ", s), names(sub)[endsWith(names(sub), ".Adj")]))
            return(sub.B0) }))
        B <- cbind(x[["Balance.Across.Subclass"]], sub.B, variable.names = row.names(x[["Balance.Across.Subclass"]]))
        if (attr(x, "print.options")$disp.subclass) attr(x, "print.options")$weight.names <- c("Adj", paste("Subclass", subclass.names))
        else attr(x, "print.options")$weight.names <- "Adj"
        subtitle <- "Across Subclasses"
        config <- "agg.none"
        facet <- NULL
    }
    else {
        B_list <- unpack_bal.tab(x)
        namesep <- attr(B_list, "namesep")
        class_sequence <- attr(B_list, "class_sequence")
        if (is_not_null(class_sequence)) {
            #Multiple layers present
            facet_mat <- as.matrix(do.call(rbind, strsplit(names(B_list), namesep, fixed = TRUE)))
            facet <- unname(vapply(class_sequence, switch, character(1L),
                                   bal.tab.cluster = "cluster",
                                   bal.tab.msm = "time",
                                   bal.tab.multi = "treat",
                                   bal.tab.imp = "imp", NULL))
            dimnames(facet_mat) <- list(names(B_list), facet)
            
            for (b in seq_along(B_list)) {
                B_list[[b]][["variable.names"]] <- rownames(B_list[[b]])
                for (i in facet) {
                    if (i == "imp") B_list[[b]][[i]] <- factor(paste("Imputation:", facet_mat[b, i]),
                                                               levels = paste("Imputation:", sort(unique(as.numeric(facet_mat[b, i])))))
                    else B_list[[b]][[i]] <- facet_mat[b, i]
                }
            }
            
            #Process which. so that B_list can be shortened
            agg.over <- character(0)
            for (i in facet) {
                which. <- attr(x, "print.options")[[paste0("which.", i)]]
                if (is_null(which.)) {
                    #All levels; facet_mat stays the same.
                }
                else if (anyNA(which.)) {
                    agg.over <- c(agg.over, i)
                }
                else {
                    if (i == "treat") {
                        treat_levels <- attr(x, "print.options")$treat_names_multi
                        if (is.numeric(which.)) which. <- treat_levels[which.]
                        if (any(which. %nin% treat_levels)) stop("All values in 'which.treat' must be names or indices of treatment levels.", call. = FALSE)
                        if (attr(x, "print.options")$pairwise) {
                            vs.combs <- cbind(vs.tmp <- as.matrix(expand.grid(treat_levels, treat_levels, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)), 
                                              apply(vs.tmp, 1, paste, collapse = " vs. "))
                            vs.combs <- vs.combs[vs.combs[,3] %in% facet_mat[, i],]
                            if (length(which.) == 1) facet_mat <- facet_mat[facet_mat[,i] %in% vs.combs[,3][vs.combs[,1] == which. | vs.combs[,2] == which.], , drop = FALSE]
                            else facet_mat <- facet_mat[facet_mat[,i] %in% vs.combs[,3][vs.combs[,1] %in% which. & vs.combs[,2] %in% which.], , drop = FALSE]
                        }
                        else {
                            vs.combs <- cbind(vs.tmp <- as.matrix(data.frame("Others", treat_levels, stringsAsFactors = FALSE)), 
                                              apply(vs.tmp, 1, paste, collapse = " vs. "))
                            vs.combs <- vs.combs[vs.combs[,3] %in% facet_mat[, i],]
                            facet_mat <- facet_mat[facet_mat[,i] %in% vs.combs[,3][vs.combs[,2] %in% which.], , drop = FALSE]
                        }
                    }
                    else {
                        if (is.numeric(which.) && max(which.) <= nunique(facet_mat[,i])) {
                            if (i == "imp") facet_mat <- facet_mat[facet_mat[,i] %in% as.character(which.), ,drop = FALSE]
                            facet_mat <- facet_mat[facet_mat[,i] %in% sort(unique(facet_mat[,i]))[which.], ,drop = FALSE]
                        }
                        else if (is.character(which.) && all(which. %in% unique(facet_mat[,i]))) {
                            facet_mat <- facet_mat[facet_mat[,i] %in% which., ,drop = FALSE]
                        }
                        else stop(paste0("The argument to 'which.", i, "' must be .none, .all, or the desired levels or indices of ", switch(i, time = "time points", i), "."), call. = FALSE)
                    }
                }
                
            }
            B_list <- B_list[rownames(facet_mat)]
            
            stat.cols <- expand.grid_string(vapply(stats, function(s) STATS[[s]]$bal.tab_column_prefix, character(1L)),
                                            c("Un", attr(x, "print.options")[["weight.names"]]),
                                            collapse = ".")
            cols.to.keep <- c("variable.names", "Type", facet, stat.cols)
            
            for (b in seq_along(B_list)) {
                B_list[[b]] <- B_list[[b]][cols.to.keep]
            }
            
            B_stack <- do.call(rbind, c(B_list, list(make.row.names = FALSE)))
            
            if (is_not_null(agg.over)) {
                if (is_null(agg.fun)) {
                    if (any(c("treat", "time") %in% agg.over)) agg.fun <- "max"
                    else agg.fun <- "range"
                }
                agg.fun <- tolower(agg.fun)
                Agg.Fun <- firstup(agg.fun <- match_arg(agg.fun, c("range", "max", "mean")))
                if (agg.fun == "max") abs <- TRUE
                
                if (abs) B_stack[stat.cols] <- lapply(stat.cols, function(sc) abs_(B_stack[[sc]], ratio = startsWith(sc, "V.Ratio")))
                
                facet <- facet[facet %nin% agg.over]
                
                aggregate_B <- function(FUN, B) {
                    B_agged <- aggregate(B[stat.cols], 
                                         by = B[c("variable.names", "Type", facet)], 
                                         FUN = FUN)
                    names(B_agged)[names(B_agged) %in% stat.cols] <- paste.(firstup(FUN), names(B_agged)[names(B_agged) %in% stat.cols])
                    return(B_agged)
                }
                
                if (agg.fun == "range") {
                    B <- Reduce(function(x, y) merge(x, y, by = c("variable.names", "Type", facet), 
                                                     sort = FALSE),
                                lapply(c("min", "mean", "max"), aggregate_B, B_stack))
                }
                else {
                    B <- aggregate_B(agg.fun, B_stack)
                }
                
                subtitle1 <- paste0(Agg.Fun, " across ", word_list(vapply(agg.over, switch, character(1L),
                                                                          "cluster" = "clusters",
                                                                          "time" = "time points",
                                                                          "treat" = "treatment pairs",
                                                                          "imp" = "imputations")))
                config <- paste.("agg", agg.over)
            }
            else {
                B <- B_stack
                subtitle1 <- NULL
                config <- "agg.none"
            }
            
            one.level.facet <- facet[vapply(B[facet], all_the_same, logical(1L))]
            if (is_not_null(one.level.facet)) {
                subtitle2 <- paste(vapply(one.level.facet, function(olf) {
                    paste(firstup(olf), B[1,olf], sep = ": ")
                }, character(1L)), collapse = ", ")
            }
            else subtitle2 <- NULL
            
            B[names(B) %in% one.level.facet] <- NULL
            
            if (sum(facet %nin% one.level.facet) > 1) stop(paste("At least one of", word_list(paste.("which", facet), "or", quotes = 1), "must be .none or of length 1."), call. = FALSE)
            
            facet <- facet[facet %nin% one.level.facet]
            
            subtitle <- paste(c(subtitle1, subtitle2), collapse = "\n")
            
            #one.level.facet - go in subtitle
            #facet - go in facet
            #agg.over - aggregated (e.g., averaged) over
            
        }
        else {
            #Single-layer bal.tab
            B <- data.frame(B_list, variable.names = rownames(B_list))
            
            facet <- one.level.facet <- agg.over <- NULL
            
            stat.cols <- expand.grid_string(vapply(stats, function(s) STATS[[s]]$bal.tab_column_prefix, character(1L)),
                                            c("Un", attr(x, "print.options")[["weight.names"]]),
                                            collapse = ".")
            cols.to.keep <- c("variable.names", "Type", stat.cols)
            B <- B[cols.to.keep]
            
            config <- "agg.none"
            
            subtitle <- NULL
        }
        subclass.names <- NULL
    }
    
    if (is_not_null(facet) && length(stats) > 1) {
        stop("'stats' can only have a length of 1 when faceting by other dimension (e.g., cluster, treatment).", call. = FALSE)
    }
    
    if (is_not_null(agg.fun) && config == "agg.none") {
        warning("No aggregation will take place, so 'agg.fun' will be ignored. Remember to set 'which.<ARG> = .none' to aggregate across <ARG>.", call. = FALSE)
    }
    
    #Process variable names
    if (is_not_null(var.names)) {
        if (is.data.frame(var.names)) {
            if (ncol(var.names)==1) {
                if (is_not_null(row.names(var.names))) {
                    new.labels <- setNames(unlist(as.character(var.names[,1])), rownames(var.names))
                }
                else warning("'var.names' is a data frame, but its rows are unnamed.", call. = FALSE)
            }
            else {
                if (all(c("old", "new") %in% names(var.names))) {
                    new.labels <- setNames(unlist(as.character(var.names[,"new"])), var.names[,"old"])
                }
                else {
                    if (ncol(var.names)>2) warning("Only using first 2 columns of 'var.names'.", call. = FALSE)
                    new.labels <- setNames(unlist(as.character(var.names[,2])), var.names[,1])
                }
            } 
        }
        else if (is_(var.names, "atomic")) {
            if (is_not_null(names(var.names))) {
                new.labels <- setNames(as.character(var.names), names(var.names))
            }
            else warning("'var.names' is a vector, but its values are unnamed.", call. = FALSE)
        }
        else if (is_(var.names, "list")) {
            if (all(sapply(var.names, function(x) is_(x, c("character", "factor"))))) {
                if (is_not_null(names(var.names))) {
                    new.labels <- unlist(var.names) #already a list
                }
                else warning("'var.names' is a list, but its values are unnamed.", call. = FALSE)
            }
            else warning("'var.names' is a list, but its values are not the new names of the variables.", call. = FALSE)
        }
        else warning("Argument to 'var.names' is not one of the accepted structures and will be ignored.\n  See help(love.plot) for details.", immediate.=TRUE, call. = FALSE)
        
        co.names <- attr(x, "print.options")[["co.names"]]
        seps <- attr(co.names, "seps")
        for (i in names(co.names)) {
            comp <- co.names[[i]][["component"]]
            type <- co.names[[i]][["type"]]
            
            if (i %in% names(new.labels) && !is.na(new.labels[i])) {
                co.names[[i]][["component"]] <- new.labels[i]
                co.names[[i]][["type"]] <- "base"
            }
            else {
                if ("isep" %in% type) {
                    named.vars <- character(sum(type == "isep") + 1)
                    sep.inds <- c(which(type == "isep"), length(comp) + 1)
                    named.vars <- lapply(seq_along(sep.inds), function(k) {
                        inds <- (if (k == 1) seq(1, sep.inds[k] - 1) 
                                 else seq(sep.inds[k-1] + 1, sep.inds[k] - 1))
                        var <- comp[inds]
                        var.is.base <- type[inds] == "base"
                        pasted.var <- paste(var, collapse = "")
                        if (pasted.var %in% names(new.labels)) return(new.labels[pasted.var])
                        else return(paste(ifelse(var.is.base & var %in% names(new.labels) & !is.na(new.labels[var]), new.labels[var], var), collapse = ""))
                    })
                    co.names[[i]][["component"]] <- do.call("paste", c(unname(named.vars), list(sep = seps["int"])))
                }
                else co.names[[i]][["component"]] <- ifelse(type == "base" & comp %in% names(new.labels) & !is.na(new.labels[comp]), new.labels[comp], comp)
            }
        }
        
        recode.labels <- setNames(names(co.names), 
                                  vapply(co.names, function(x) paste0(x[["component"]], collapse = ""), character(1L)))
        
        B[["variable.names"]] <- do.call(f.recode, c(list(B[["variable.names"]]), recode.labels))
    }
    
    distance.names <- as.character(unique(B[["variable.names"]][B[["Type"]] == "Distance"], nmax = sum(B[["Type"]] == "Distance")))
    if (drop.distance) {
        B <- B[B[["variable.names"]] %nin% distance.names, , drop = FALSE]
    }
    
    #Process variable order
    if (is_not_null(var.order) && "love.plot" %nin% class(var.order)) {
        if (is_null(attr(x, "print.options")$nweights) ||
            attr(x, "print.options")$nweights == 0) {
            ua <- c("Unadjusted", "Alphabetical")
            names(ua) <- c("unadjusted", "alphabetical")
        }
        else if (attr(x, "print.options")$nweights == 1) {
            ua <- c("Adjusted", "Unadjusted", "Alphabetical")
            names(ua) <- c("adjusted", "unadjusted", "alphabetical")
        }
        else {
            ua <- c("Unadjusted", attr(x, "print.options")$weight.names, "Alphabetical")
            names(ua) <- c("unadjusted", attr(x, "print.options")$weight.names, "alphabetical")
        }
        var.order <- ua[match_arg(var.order, tolower(ua))]
    }
    
    #Process sample names
    ntypes <- if (is_null(subclass.names)) length(attr(x, "print.options")$weight.names) + 1 else 2
    if (!missing(sample.names)) {
        if (!is.vector(sample.names, "character")) {
            warning("The argument to 'sample.names' must be a character vector. Ignoring 'sample.names'.", call. = FALSE)
            sample.names <- NULL
        }
        else if (length(sample.names) %nin% c(ntypes, ntypes - 1)) {
            warning("The argument to 'sample.names' must contain as many names as there are sample types, or one fewer. Ignoring 'sample.names'.", call. = FALSE)
            sample.names <- NULL
        }
    }
    else sample.names <- NULL
    
    #Process limits
    if (is_not_null(limits)) {
        if (!is_(limits, "list")) {
            limits <- list(limits)
        }
        if (any(vapply(limits, 
                       function(l) !is_(l, "numeric") || length(l) %nin% c(0L, 2L), 
                       logical(1L)))) {
            warning("'limits' must be a list of numeric vectors of legnth 2. Ignoring 'limits'.", call. = FALSE)
            limits <- NULL
        }
        
        if (is_not_null(names(limits))) {
            names(limits) <- stats[pmatch(names(limits), stats, duplicates.ok = TRUE)]
            limits <- limits[!is.na(names(limits))]
        }
        else {
            names(limits) <- stats[1:length(limits)]
        }
    }
    
    #Setting up appearance
    
    #Alpha (transparency)
    if (is.numeric(alpha[1]) && 
        !anyNA(alpha[1]) && 
        between(alpha[1], c(0,1))) alpha <- alpha[1]
    else {
        warning("The argument to 'alpha' must be a number between 0 and 1. Using 1 instead.", call. = FALSE)
        alpha <- 1
    }
    
    #Color
    if (is_not_null(args[["colours"]])) colors <- args[["colours"]]
    
    if (is_null(colors)) {
        if (shapes.ok(shapes, ntypes) && length(shapes) > 1 && length(shapes) == ntypes) {
            colors <- rep("black", ntypes)
        }
        else colors <- gg_color_hue(ntypes)
    }
    else {
        if (length(colors) == 1) colors <- rep(colors, ntypes)
        else if (length(colors) > ntypes) {
            colors <- colors[seq_len(ntypes)]
            warning(paste("Only using first", ntypes, "value", if (ntypes > 1) "s " else " ", "in 'colors'."), call. = FALSE)
        }
        else if (length(colors) < ntypes) {
            warning("Not enough colors were specified. Using default colors instead.", call. = FALSE)
            colors <- gg_color_hue(ntypes)
        }
        
        if (!all(vapply(colors, isColor, logical(1L)))) {
            warning("The argument to 'colors' contains at least one value that is not a recognized color. Using default colors instead.", call. = FALSE)
            colors <- gg_color_hue(ntypes)
        }
        
    }
    # colors[] <- vapply(colors, col_plus_alpha, character(1L), alpha = alpha)
    fill <- colors
    
    #Shapes
    if (is_null(shapes)) {
        shapes <- assign.shapes(colors)
    }
    else {
        #check shapes
        if (!shapes.ok(shapes, ntypes)) {
            warning(paste("The argument to shape must be", ntypes, "valid shape", if (ntypes > 1) "s." else ".", " See ?love.plot for more information.\nUsing default shapes instead."), call. = FALSE)
            shapes <- assign.shapes(colors)
        }
        else if (length(shapes) == 1) shapes <- rep(shapes, ntypes)
    }
    
    #Size
    if (is.numeric(size)) size <- size[1]
    else {
        warning("The argument to size must be a number. Using 3 instead.", call. = FALSE)
        size <- 3
    }
    
    stroke <- rep(0, ntypes)
    size0 <- size <- rep(size, ntypes)
    
    shapes.with.fill <- grepl("filled", shapes, fixed = TRUE)
    stroke[shapes.with.fill] <- size[shapes.with.fill]/3
    size[shapes.with.fill] <- size[shapes.with.fill]* .58
    
    # stroke <- .8*size
    
    if (is_not_null(facet)) {
        if (is_not_null(var.order) && "love.plot" %nin% class(var.order) && tolower(var.order) != "alphabetical") {
            warning("'var.order' cannot be set with faceted plots (unless \"alphabetical\"). Ignoring 'var.order'.", call. = FALSE)
            var.order <- NULL
        }
    }
    
    agg.range <- isTRUE(Agg.Fun == "Range")
    
    #Process thresholds
    thresholds <- if_null_then(attr(x, "print.options")$thresholds[stats], 
                               process_thresholds(thresholds, stats))
    
    #Title
    if (missing(title)) title <- "Covariate Balance"
    else title <- as.character(title)
    # if (missing(subtitle)) subtitle <- as.character(subtitle)
    
    #Process themes
    if (is_not_null(themes)) {
        if (!is.vector(themes, "list")) {
            themes <- list(themes)
        }
        if (any(vapply(themes, 
                       function(t) !all(c("theme", "gg") %in% class(t)), 
                       logical(1L)))) {
            warning("'themes' must be a list of \"theme\" objects. Ignoring 'themes'.", call. = FALSE)
            themes <- NULL
        }
        
        if (is_not_null(names(themes))) {
            names(themes) <- stats[pmatch(names(themes), stats, duplicates.ok = TRUE)]
            themes <- themes[!is.na(names(themes))]
        }
        else {
            names(themes) <- stats[1:length(themes)]
        }
    }
    
    plot.list <- make_list(stats)
    for (s in stats) {
        variable.names <- as.character(B[["variable.names"]])
        
        #Get SS
        if (agg.range) {
            SS <- do.call("rbind", 
                          lapply(c("Un", attr(x, "print.options")$weight.names),
                                 function(w) data.frame(var = variable.names,
                                                        type = B[["Type"]],
                                                        min.stat = B[[paste.("Min", STATS[[s]]$bal.tab_column_prefix, w)]],
                                                        max.stat = B[[paste.("Max", STATS[[s]]$bal.tab_column_prefix, w)]],
                                                        mean.stat = B[[paste.("Mean", STATS[[s]]$bal.tab_column_prefix, w)]],
                                                        Sample = switch(w, "Un"= "Unadjusted", 
                                                                        "Adj" = "Adjusted", w),
                                                        B[facet],
                                                        row.names = NULL)))
            
            if (all(sapply(SS[c("min.stat", "max.stat", "mean.stat")], is.na))) 
                stop(paste("No balance statistics to display. This can occur when", 
                           STATS[[s]]$disp_stat, 
                           "= FALSE and quick = TRUE in the original call to bal.tab()."), call. = FALSE)
            
            missing.stat <- all(is.na(SS[["mean.stat"]]))
            if (missing.stat) stop(paste0(word_list(firstup(STATS[[s]]$balance_tally_for)), 
                                          " cannot be displayed. This can occur when ", 
                                          word_list(STATS[[s]]$disp_stat, and.or = "and", is.are = TRUE), 
                                          " FALSE and quick = TRUE in the original call to bal.tab()."), call. = FALSE)
            
            gone <- character(0)
            for (i in levels(SS$Sample)) {
                if (all(sapply(SS[SS$Sample == i, c("min.stat", "max.stat", "mean.stat")], is.na))) {
                    gone <- c(gone, i)
                    if (i == "Unadjusted") warning("Unadjusted values are missing. This can occur when un = FALSE and quick = TRUE in the original call to bal.tab().", call. = FALSE, immediate. = TRUE)
                    SS <- SS[SS[["Sample"]] != i,]
                }
            }
            
            dec <- FALSE
            
            if (is_not_null(plot.list[[1]])) var.order <- plot.list[[1]]
            
            if (is_not_null(var.order)) {
                if ("love.plot" %in% class(var.order)) {
                    old.vars <- levels(var.order$data$var)
                    old.vars[endsWith(old.vars, "*")] <- substr(old.vars[endsWith(old.vars, "*")], 1, nchar(old.vars[endsWith(old.vars, "*")])-1)
                    if (any(SS[["var"]] %nin% old.vars)) {
                        warning("The love.plot object in 'var.order' doesn't have the same variables as the current input. Ignoring 'var.order'.", call. = FALSE)
                        var.order <- NULL
                    }
                    else {
                        SS[["var"]] <- factor(SS[["var"]], levels = old.vars[old.vars %in% SS[["var"]]])
                    }
                }
                else if (tolower(var.order) == "alphabetical") {
                    if ("time" %in% facet) {
                        covnames0 <- make_list(length(unique(SS[["time"]])))
                        for (i in seq_along(covnames0)) {
                            if (i == 1) {
                                covnames0[[i]] <- sort(levels(SS[["var"]][SS[["time"]] == i]))
                            }
                            else {
                                covnames0[[i]] <- sort(setdiff(levels(SS[["var"]][SS[["time"]] == i]), unlist(covnames0[seq_along(covnames0) < i])))
                            }
                        }
                        covnames <- unlist(covnames0)
                    }
                    else covnames <- sort(levels(SS[["var"]]))
                    SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[!covnames %in% distance.names]), sort(distance.names, decreasing = TRUE)))
                    
                }
                else if (var.order %in% ua) {
                    if (var.order %in% gone) {
                        warning(paste0("'var.order' was set to \"", tolower(var.order), "\", but no ", tolower(var.order), " ", STATS[[s]]$balance_tally_for, " were calculated. Ignoring 'var.order'."), call. = FALSE, immediate. = TRUE)
                        var.order <- NULL
                    }
                    else {
                        v <- as.character(SS[["var"]][order(SS[["mean.stat"]][SS[["Sample"]]==var.order], decreasing = dec, na.last = FALSE)])
                        
                        SS[["var"]] <- factor(SS[["var"]], 
                                              levels=c(v[v %nin% distance.names], 
                                                       sort(distance.names, decreasing = TRUE)))
                    }
                }
                
            }
            if (is_null(var.order)) {
                covnames <- as.character(unique(SS[["var"]]))
                SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[!covnames %in% distance.names]), sort(distance.names, decreasing = TRUE)))
            }
            # SS[, "Sample"] <- factor(SS[, "Sample"], levels = c("Adjusted", "Unadjusted"))
            SS[["Sample"]] <- factor(SS[["Sample"]])
            if (s == "mean.diffs" && any(base::abs(SS[["max.stat"]]) > 5, na.rm = TRUE)) warning("Large mean differences detected; you may not be using standardized mean differences for continuous variables.", call.=FALSE)
            if (length(stats) == 1 && drop.missing) SS <- SS[!is.na(SS[["min.stat"]]),]
            SS[["stat"]] <- SS[["mean.stat"]]
        }
        else {
            
            SS <- do.call("rbind", 
                          lapply(c("Un", attr(x, "print.options")$weight.names),
                                 function(w) data.frame(var = variable.names,
                                                        type = B[["Type"]],
                                                        stat = B[[ifelse(is_null(Agg.Fun), paste.(STATS[[s]]$bal.tab_column_prefix, w),
                                                                         paste.(Agg.Fun, STATS[[s]]$bal.tab_column_prefix, w))]],
                                                        Sample = switch(w, "Un"= "Unadjusted",
                                                                        "Adj" = "Adjusted", w),
                                                        B[facet],
                                                        row.names = NULL)
                          ))
            
            missing.stat <- all(is.na(SS[["stat"]]))
            if (missing.stat) stop(paste0(word_list(firstup(STATS[[s]]$balance_tally_for)), 
                                          " cannot be displayed. This can occur when ", 
                                          word_list(STATS[[s]]$disp_stat, and.or = "and", is.are = TRUE), 
                                          " FALSE and quick = TRUE in the original call to bal.tab()."), call. = FALSE)
            
            gone <- character(0)
            for (i in levels(SS$Sample)) {
                if (all(is.na(SS[["stat"]][SS$Sample==i]))) {
                    gone <- c(gone, i)
                    if (i == "Unadjusted") warning("Unadjusted values are missing. This can occur when un = FALSE and quick = TRUE in the original call to bal.tab().", call. = FALSE, immediate. = TRUE)
                    SS <- SS[SS[["Sample"]]!=i,]
                }
            }
            
            if (abs) {
                SS[["stat"]] <- abs_(SS[["stat"]], ratio = s == "variance.ratios")
            }
            dec <- FALSE
            
            if (is_not_null(plot.list[[1]])) var.order <- plot.list[[1]]
            
            if (is_not_null(var.order)) {
                if ("love.plot" %in% class(var.order)) {
                    old.vars <- levels(var.order$data$var)
                    old.vars[endsWith(old.vars, "*")] <- substr(old.vars[endsWith(old.vars, "*")], 1, nchar(old.vars[endsWith(old.vars, "*")])-1)
                    if (any(SS[["var"]] %nin% old.vars)) {
                        warning("The love.plot object in 'var.order' doesn't have the same variables as the current input. Ignoring 'var.order'.", call. = FALSE)
                        var.order <- NULL
                    }
                    else {
                        SS[["var"]] <- factor(SS[["var"]], levels = old.vars[old.vars %in% SS[["var"]]])
                    }
                }
                else if (tolower(var.order) == "alphabetical") {
                    if ("time" %in% facet) {
                        covnames0 <- make_list(length(unique(SS[["time"]])))
                        for (i in seq_along(covnames0)) {
                            if (i == 1) {
                                covnames0[[i]] <- sort(levels(SS[["var"]][SS[["time"]] == i]))
                            }
                            else {
                                covnames0[[i]] <- sort(setdiff(levels(SS[["var"]][SS[["time"]] == i]), unlist(covnames0[seq_along(covnames0) < i])))
                            }
                        }
                        covnames <- unlist(covnames0)
                    }
                    else covnames <- sort(levels(SS[["var"]]))
                    SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[!covnames %in% distance.names]), sort(distance.names, decreasing = TRUE)))
                    
                }
                else if (var.order %in% ua) {
                    if (var.order %in% gone) {
                        warning(paste0("var.order was set to \"", tolower(var.order), "\", but no ", tolower(var.order), " ", STATS[[s]]$balance_tally_for, " were calculated. Ignoring var.order."), call. = FALSE, immediate. = TRUE)
                        var.order <- NULL
                    }
                    else {
                        v <- as.character(SS[["var"]][order(SS[["stat"]][SS[["Sample"]]==var.order], decreasing = dec, na.last = FALSE)])
                        
                        SS[["var"]] <- factor(SS[["var"]], 
                                              levels=c(v[v %nin% distance.names], 
                                                       sort(distance.names, decreasing = TRUE)))
                    }
                }
                
            }
            if (is_null(var.order)) {
                covnames <- as.character(unique(SS[["var"]])) #Don't use levels here to preserve original order
                SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[covnames %nin% distance.names]), sort(distance.names, decreasing = TRUE)))
            }
            SS[["Sample"]] <- factor(SS[["Sample"]])
            if (s == "mean.diffs" && any(base::abs(SS[["stat"]]) > 5, na.rm = TRUE)) warning("Large mean differences detected; you may not be using standardized mean differences for continuous variables.", call.=FALSE)
            if (length(stats) == 1 && drop.missing) SS <- SS[!is.na(SS[["stat"]]),]
        }
        
        SS <- SS[order(SS[["var"]], na.last = FALSE),]
        SS[["var"]] <- factor(SS[["var"]])
        
        #Make the plot
        #library(ggplot2)
        
        baseline.xintercept <- STATS[[s]]$baseline.xintercept
        if (is_not_null(thresholds[[s]])) threshold.xintercepts <- STATS[[s]]$threshold.xintercepts(thresholds[[s]], abs)
        else threshold.xintercepts <- NULL
        xlab <- STATS[[s]]$love.plot_xlab(abs = abs, binary = attr(x, "print.options")$binary,
                                          continuous = attr(x, "print.options")$continuous,
                                          var_type = B[["Type"]],
                                          stars = stars)
        SS[["var"]] <- STATS[[s]]$love.plot_add_stars(SS[["var"]], 
                                                      variable.names = variable.names,
                                                      binary = attr(x, "print.options")$binary,
                                                      continuous = attr(x, "print.options")$continuous,
                                                      var_type = B[["Type"]],
                                                      stars = stars,
                                                      star_char = args$star_char)
        scale_Statistics <- STATS[[s]]$love.plot_axis_scale
        
        apply.limits <- FALSE
        SS[["on.border"]] <- FALSE
        if (is_not_null(limits[[s]])) {
            if (limits[[s]][2] < limits[[s]][1]) {
                limits[[s]] <- c(limits[[s]][2], limits[[s]][1])
            }
            
            if (limits[[s]][1] >= baseline.xintercept) limits[[s]][1] <- baseline.xintercept - .05*limits[[s]][2]
            if (limits[[s]][2] <= baseline.xintercept) limits[[s]][2] <- baseline.xintercept - .05*limits[[s]][1]
            
            if (identical(scale_Statistics, ggplot2::scale_x_log10)) limits[[s]][limits[[s]] <= 1e-2] <- 1e-2
            
            if (agg.range) {
                
                if (any(SS[["mean.stat"]] < limits[[s]][1], na.rm = TRUE)) {
                    SS[["on.border"]][SS[["mean.stat"]] < limits[[s]][1]] <- TRUE
                    SS[["mean.stat"]][SS[["mean.stat"]] < limits[[s]][1]] <- limits[[s]][1]
                    SS[["max.stat"]][SS[["max.stat"]] < limits[[s]][1]] <- limits[[s]][1]
                    SS[["min.stat"]][SS[["min.stat"]] < limits[[s]][1]] <- limits[[s]][1]
                }
                if (any(SS[["mean.stat"]] > limits[[s]][2], na.rm = TRUE)) {
                    SS[["on.border"]][SS[["mean.stat"]] > limits[[s]][2]] <- TRUE
                    SS[["mean.stat"]][SS[["mean.stat"]] > limits[[s]][2]] <- limits[[s]][2]
                    SS[["max.stat"]][SS[["max.stat"]] > limits[[s]][2]] <- limits[[s]][2]
                    SS[["min.stat"]][SS[["min.stat"]] > limits[[s]][2]] <- limits[[s]][2]
                    # warning("Some points will be removed from the plot by the limits.", call. = FALSE)
                }
                # warning("Some points will be removed from the plot by the limits.", call. = FALSE)
            }
            else {
                if (any(SS[["stat"]] < limits[[s]][1], na.rm = TRUE)) {
                    SS[["on.border"]][SS[["stat"]] < limits[[s]][1]] <- TRUE
                    SS[["stat"]][SS[["stat"]] < limits[[s]][1]] <- limits[[s]][1]
                }
                if (any(SS[["stat"]] > limits[[s]][2], na.rm = TRUE)) {
                    SS[["on.border"]][SS[["stat"]] > limits[[s]][2]] <- TRUE
                    SS[["stat"]][SS[["stat"]] > limits[[s]][2]] <- limits[[s]][2]
                    # warning("Some points will be removed from the plot by the limits.", call. = FALSE)
                }
            }
            
            apply.limits <- TRUE
        }
        
        if (is_not_null(sample.names)) {
            if (length(sample.names) == ntypes - 1) {
                levels(SS$Sample)[-1] <- sample.names
            }
            else if (length(sample.names) == ntypes) {
                levels(SS$Sample) <- sample.names
            }
        }
        
        lp <- ggplot2::ggplot(aes(y = .data$var, x = .data$stat, group = .data$Sample), data = SS) +
            ggplot2::theme(panel.background = element_rect(fill = "white"),
                           axis.text.x = element_text(color = "black"),
                           axis.text.y = element_text(color = "black"),
                           panel.border = element_rect(fill = NA, color = "black"),
                           plot.background = element_blank(),
                           legend.background = element_blank(),
                           legend.key = element_blank()
            ) +
            ggplot2::scale_shape_manual(values = shapes) +
            ggplot2::scale_size_manual(values = size) +
            ggplot2::scale_discrete_manual(aesthetics = "stroke", values = stroke) +
            ggplot2::scale_fill_manual(values = fill) +
            ggplot2::scale_color_manual(values = colors) +
            ggplot2::labs(y = NULL, x = wrap(xlab, wrap))
        
        lp <- lp + ggplot2::geom_vline(xintercept = baseline.xintercept,
                                       linetype = 1, color = "gray5")
        
        if (is_not_null(threshold.xintercepts)) {
            lp <- lp + ggplot2::geom_vline(xintercept = threshold.xintercepts,
                                           linetype = 2, color = "gray8")
        }
        
        if (agg.range) {
            position.dodge <- ggplot2::position_dodge(.5*(size0[1]/3))
            if (line == TRUE) { #Add line except to distance
                f <- function(q) {q[["stat"]][q$type == "Distance"] <- NA; q}
                lp <- lp + ggplot2::layer(geom = "path", data = f, 
                                          position = position.dodge, 
                                          stat = "identity",
                                          mapping = aes(x = .data$mean.stat, color = .data$Sample), 
                                          params = list(size = size0[1]*.8/3, na.rm = TRUE,
                                                        alpha = alpha))
            }
            
            lp <- lp +
                ggplot2::geom_linerange(aes(y = .data$var, xmin = .data$min.stat, xmax = .data$max.stat,
                                            color = .data$Sample), position = position.dodge,
                                        size = size0[1]*.8/3,
                                        alpha = alpha, 
                                        orientation = "y",
                                        show.legend = FALSE,
                                        na.rm = TRUE) +
                ggplot2::geom_point(aes(y = .data$var, 
                                        x = .data$mean.stat, 
                                        shape = .data$Sample,
                                        size = .data$Sample,
                                        stroke = .data$Sample,
                                        color = .data$Sample),
                                    fill = "white", na.rm = TRUE,
                                    alpha = alpha,
                                    position = position.dodge)
            
        }
        else {
            if (is_null(subclass.names) || !attr(x, "print.options")$disp.subclass) {
                if (isTRUE(line)) { #Add line except to distance
                    f <- function(q) {q[["stat"]][q$type == "Distance"] <- NA; q}
                    lp <- lp + ggplot2::layer(geom = "path", data = f(SS),
                                              position = "identity", stat = "identity",
                                              mapping = aes(color = .data$Sample),
                                              params = list(size = size0[1]*.8/3,
                                                            na.rm = TRUE,
                                                            alpha = alpha))
                }
                lp <- lp + ggplot2::geom_point(data = SS, aes(shape = .data$Sample,
                                                              size = .data$Sample,
                                                              stroke = .data$Sample,
                                                              color = .data$Sample),
                                               fill = "white", 
                                               na.rm = TRUE,
                                               alpha = alpha)
                
            }
            else {
                SS.u.a <- SS[SS$Sample %in% c("Unadjusted", "Adjusted"),]
                SS.u.a$Sample <- factor(SS.u.a$Sample)
                if (line == TRUE) { #Add line except to distance
                    f <- function(q) {q[["stat"]][q$type == "Distance"] <- NA; q}
                    lp <- lp + ggplot2::layer(geom = "path", data = f(SS.u.a),
                                              position = "identity", stat = "identity",
                                              mapping = aes(color = .data$Sample),
                                              params = list(size = size*.8,
                                                            na.rm = TRUE,
                                                            alpha = alpha))
                }
                lp <- lp + ggplot2::geom_point(data = SS.u.a,
                                               aes(shape = .data$Sample,
                                                   size = .data$Sample,
                                                   stroke = .data$Sample,
                                                   color = .data$Sample),
                                               fill = "white",
                                               na.rm = TRUE)
                lp <- lp + ggplot2::geom_text(data = SS[SS$Sample %nin% c("Unadjusted", "Adjusted"),],
                                              mapping = aes(label = gsub("Subclass ", "", .data$Sample)),
                                              size = 2.5*size0[1]/3, na.rm = TRUE)
            }
            
            
        }
        
        if (!drop.distance && is_not_null(distance.names)) {
            lp <- lp + ggplot2::geom_hline(linetype = 1, color = "black",
                                           yintercept = nunique(SS[["var"]]) - length(distance.names) + .5)
        }
        if (apply.limits) {
            lp <- lp + scale_Statistics(limits = limits[[s]], expand = c(0, 0))
        }
        else {
            lp <- lp + scale_Statistics()
        }
        
        if (isFALSE(grid)) {
            lp <- lp + ggplot2::theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank())
        }
        else {
            lp <- lp + ggplot2::theme(panel.grid.major = element_line(color = "gray87"),
                                      panel.grid.minor = element_line(color = "gray90"))
        }
        
        if (is_not_null(facet)) {
            lp <- lp + ggplot2::facet_grid(f.build(".", facet), drop = FALSE) + ggplot2::labs(x = xlab)
        }
        
        class(lp) <- c(class(lp), "love.plot")
        plot.list[[s]] <- lp
    }
    
    if (length(stats) > 1 || isTRUE(args$use.grid)) {
        
        if (!rlang::is_string(position)) {
            position <- NA_character_
        }
        else position <- match_arg(position, 
                                   c("right", "left", "top", "bottom", "none"))
        
        #Process labels
        if (isTRUE(labels)) labels <- LETTERS[seq_along(plot.list)]
        else if (is_null(labels) || isFALSE(labels)) labels <- NULL
        else if (!is.atomic(labels) || length(labels) != length(plot.list)) {
            warning("'labels' must be TRUE or a string with the same length as 'stats'. Ignoring 'labels'.", call. = FALSE)
            labels <- NULL
        }
        else labels <- as.character(labels)
        
        # p <- ggpubr::ggarrange(plotlist = plot.list, common.legend = TRUE, legend = position, 
        #                align = "hv", nrow = 1)
        # if (is_not_null(subtitle)) {
        #     p <- ggpubr::annotate_figure(p, top = ggpubr::text_grob(subtitle, size = 11))
        # }
        # p <- ggpubr::annotate_figure(p, top = ggpubr::text_grob(title, size = 13.2))
        # 
        # P <- attr(P, "plots")
        
        plots.to.combine <- plot.list
        for (i in seq_along(plots.to.combine)) {
            if (i > 1) {
                plots.to.combine[[i]] <- plots.to.combine[[i]] + 
                    ggplot2::theme(axis.text.y=element_blank(),
                                   axis.ticks.y=element_blank(),
                                   legend.position = "none")
            }
            else {
                plots.to.combine[[i]] <- plots.to.combine[[i]] + ggplot2::theme(legend.position = "none")
            }
            
            if (is_not_null(labels)) {
                plots.to.combine[[i]] <- plots.to.combine[[i]] + ggplot2::labs(title = labels[i])
            }
            
            if (is_not_null(themes[[stats[i]]])) {
                plots.to.combine[[i]] <- plots.to.combine[[i]] + themes[[stats[i]]]
            }
        }
        
        g <- ggarrange_simple(plots = plots.to.combine, nrow = 1)
        title.grob <- grid::textGrob(title, gp = grid::gpar(fontsize=13.2))
        subtitle.grob <- grid::textGrob(subtitle, gp = grid::gpar(fontsize=13.2))
        
        if (position == "none") {
            p <- gridExtra::arrangeGrob(grobs = list(g), nrow = 1)
        }
        else {
            legg <- ggplot2::ggplotGrob(plots.to.combine[[1]] + ggplot2::theme(legend.position = position))
            leg <- legg$grobs[[which(legg$layout$name == "guide-box")]]
            
            if (position == "left") {
                p <- gridExtra::arrangeGrob(grobs = list(leg, g), nrow = 1, 
                                            widths = grid::unit.c(sum(leg$widths), grid::unit(1, "npc") - sum(leg$widths)))
            }
            else if (position == "right") {
                p <- gridExtra::arrangeGrob(grobs = list(g, leg), nrow = 1, 
                                            widths = grid::unit.c(grid::unit(1, "npc") - sum(leg$widths), sum(leg$widths)))
            }
            else if (position == "top") {
                p <- gridExtra::arrangeGrob(grobs = list(leg, g), nrow = 2,
                                            heights = grid::unit.c(sum(leg$heights), grid::unit(1, "npc") - sum(leg$heights)))
            }
            else if (position == "bottom") {
                p <- gridExtra::arrangeGrob(grobs = list(g, leg), nrow = 2,
                                            heights = grid::unit.c(grid::unit(1, "npc") - sum(leg$heights), sum(leg$heights)))
            }
        }
        
        if (is_not_null(subtitle)) {
            p <- gridExtra::arrangeGrob(p, top = subtitle.grob)
        }
        p <- gridExtra::arrangeGrob(p, top = title.grob)
        
        grid::grid.newpage()
        grid::grid.draw(p)
        
        attr(p, "plots") <- plot.list
        class(p) <- c(class(p), "love.plot")
        
        return(invisible(p))
    }
    else {
        
        p <- plot.list[[1]] + 
            ggplot2::labs(title = title, subtitle = subtitle) +
            ggplot2::theme(plot.title = element_text(hjust = 0.5),
                           plot.subtitle = element_text(hjust = 0.5),
                           legend.position = position)
        
        if (is_not_null(themes[[1]])) {
            p <- p + themes[[1]]
        }
        
        return(p)
        
    }
    
}
autoplot.bal.tab <- love.plot
plot.bal.tab <- autoplot.bal.tab