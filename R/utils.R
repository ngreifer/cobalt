#Strings
word_list <- function(word.list = NULL, and.or = "and", is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately
  
  word.list <- setdiff(word.list, c(NA_character_, ""))
  
  if (is_null(word.list)) {
    out <- ""
    attr(out, "plural") <- FALSE
    return(out)
  }
  
  word.list <- add_quotes(word.list, quotes)
  
  L <- length(word.list)
  
  if (L == 1L) {
    out <- word.list
    if (is.are) out <- paste(out, "is")
    attr(out, "plural") <- FALSE
    return(out)
  }
  
  if (is_null(and.or) || isFALSE(and.or)) {
    out <- toString(word.list)
  }
  else {
    and.or <- match_arg(and.or, c("and", "or"))
    
    if (L == 2L) {
      out <- sprintf("%s %s %s",
                     word.list[1L],
                     and.or,
                     word.list[2L])
    }
    else {
      out <- sprintf("%s, %s %s",
                     toString(word.list[-L]),
                     and.or,
                     word.list[L])
    }
  }
  
  if (is.are) out <- sprintf("%s are", out)
  
  attr(out, "plural") <- TRUE
  
  out
}
add_quotes <- function(x, quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }
  
  if (isTRUE(quotes)) {
    quotes <- '"'
  }
  
  if (chk::vld_string(quotes)) {
    return(paste0(quotes, x, str_rev(quotes)))
  }
  
  if (!chk::vld_count(quotes) || quotes > 2L) {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }
  
  if (quotes == 0L) {
    return(x)
  }
  
  x <- {
    if (quotes == 1L) sprintf("'%s'", x)
    else sprintf('"%s"', x)
  }
  
  x
}
expand_grid_string <- function(..., collapse = "") {
  do.call("paste", c(expand.grid(...), sep = collapse))
}
num_to_superscript <- function(x) {
  nums <- setNames(c("\u2070",
                     "\u00B9",
                     "\u00B2",
                     "\u00B3",
                     "\u2074",
                     "\u2075",
                     "\u2076",
                     "\u2077",
                     "\u2078",
                     "\u2079"),
                   as.character(0:9))
  x <- as.character(x)
  splitx <- strsplit(x, "", fixed = TRUE)
  
  vapply(splitx, function(y) paste(nums[y], collapse = ""), character(1L))
}
ordinal <- function(x) {
  if (is_null(x) || !is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }
  
  if (length(x) > 1L) {
    out <- setNames(vapply(x, ordinal, character(1L)), names(x))
    return(out)
  }
  
  x0 <- abs(x)
  out <- paste0(x0, switch(substring(x0, nchar(x0), nchar(x0)),
                           "1" = "st",
                           "2" = "nd",
                           "3" = "rd",
                           "th"))
  if (x < 0) out <- sprintf("-%s", out)
  
  setNames(out, names(x))
}
firstup <- function(x) {
  #Capitalize first letter
  substr(x, 1L, 1L) <- toupper(substr(x, 1L, 1L))
  x
}
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
  if (NROW(df) == 0L || NCOL(df) == 0L) {
    return(df)
  }
  
  if (!is.data.frame(df)) {
    df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
  }
  
  rn <- rownames(df)
  cn <- colnames(df)
  
  infs <- o.negs <- array(FALSE, dim = dim(df))
  nas <- is.na(df)
  nums <- vapply(df, is.numeric, logical(1L))
  for (i in which(nums)) {
    infs[, i] <- !nas[, i] & !is.finite(df[[i]])
  }
  
  for (i in which(!nums)) {
    if (can_str2num(df[[i]])) {
      df[[i]] <- str2num(df[[i]])
      nums[i] <- TRUE
    }
  }
  
  o.negs[, nums] <- !nas[, nums] & df[nums] < 0 & round(df[nums], digits) == 0
  df[nums] <- round(df[nums], digits = digits)
  
  for (i in which(nums)) {
    df[[i]] <- format(df[[i]], scientific = FALSE, justify = "none", trim = TRUE,
                      drop0trailing = !identical(as.character(pad), "0"))
    
    if (!identical(as.character(pad), "0") && any(grepl(".", df[[i]], fixed = TRUE))) {
      s <- strsplit(df[[i]], ".", fixed = TRUE)
      s_lengths <- lengths(s)
      digits.r.of.. <- rep.int(0, NROW(df))
      digits.r.of..[s_lengths > 1L] <- nchar(vapply(s[s_lengths > 1L], `[[`, character(1L), 2L))
      max.dig <- max(digits.r.of..)
      
      dots <- ifelse(s_lengths > 1L, "", if (nzchar(as.character(pad))) "." else pad)
      pads <- vapply(max.dig - digits.r.of.., function(n) strrep(pad, n), character(1L))
      
      df[[i]] <- paste0(df[[i]], dots, pads)
    }
  }
  
  df[o.negs] <- paste0("-", df[o.negs])
  
  # Insert NA placeholders
  df[nas] <- na_vals
  df[infs] <- "N/A"
  
  if (is_not_null(rn)) rownames(df) <- rn
  if (is_not_null(cn)) names(df) <- cn
  
  df
}
text_box_plot <- function(range.list, width = 12L) {
  full.range <- range(unlist(range.list))
  if (all_the_same(full.range)) {
    for (i in seq_along(range.list)) {
      range.list[[i]][1L] <- range.list[[i]][1L] - 1e-6
      range.list[[i]][2L] <- range.list[[i]][2L] + 1e-6
    }
    full.range <- range(unlist(range.list))
  }
  ratio <- diff(full.range) / (width + 1)
  rescaled.range.list <- lapply(range.list, function(x) round(x / ratio))
  rescaled.full.range <- round(full.range / ratio)
  d <- make_df(c("Min", space(width + 1L), "Max"),
               names(range.list),
               "character")
  d[["Min"]] <- vapply(range.list, function(x) x[1L], numeric(1L))
  d[["Max"]] <- vapply(range.list, function(x) x[2L], numeric(1L))
  for (i in seq_row(d)) {
    spaces1 <- rescaled.range.list[[i]][1L] - rescaled.full.range[1L]
    dashes <- max(c(0L, diff(rescaled.range.list[[i]]) - 2L))
    spaces2 <- max(c(0L, diff(rescaled.full.range) - (spaces1 + 1L + dashes + 1L)))
    
    d[i, 2L] <- sprintf("%s|%s|%s",
                        space(spaces1),
                        strrep("-", dashes),
                        space(spaces2))
  }
  
  d
}

equivalent.factors2 <- function(f1, f2) {
  if (!chk::vld_character_or_factor(f1) ||
      !chk::vld_character_or_factor(f2)) {
    return(FALSE)
  }
  
  f1 <- as.factor(f1)
  f2 <- as.factor(f2)
  ll1 <- levels(f1)
  ll2 <- levels(f2)
  f1 <- as.integer(f1)
  f2 <- as.integer(f2)
  nl1 <- length(ll1)
  nl2 <- length(ll2)
  
  dims <- c(nl1, nl2)
  dn <- list(ll1, ll2)
  
  bin <- f1 + nl1 * (f2 - 1L)
  pd <- nl1 * nl2
  
  tab_ <- array(tabulate(bin, pd), dims, dimnames = dn)
  
  all(colSums(tab_ != 0) %in% 0:1) && all(rowSums(tab_ != 0) %in% 0:1)
}

paste. <- function(..., collapse = NULL) {
  #Like paste0 but with sep = ".'
  paste(..., sep = ".", collapse = collapse)
}
wrap <- function(s, nchar, ...) {
  vapply(s, function(s_) {
    x <- strwrap(s_, width = nchar, ...)
    paste(x, collapse = "\n")
  }, character(1L))
}
strsplits <- function(x, splits, fixed = TRUE, ...) {
  #Link strsplit but takes multiple split values.
  #Only works for one string at a time (in x).
  for (split in splits) {
    x <- unlist(strsplit(x, split, fixed = TRUE, ...))
  }
  
  x[nzchar(x)] # Remove empty values
}
#' @exportS3Method NULL
c.factor <- function(..., recursive = TRUE) {
  #c() for factors
  unlist(list(...), recursive = recursive)
}
can_str2num <- function(x) {
  if (is.numeric(x) || is.logical(x)) {
    return(TRUE)
  }
  
  nas <- is.na(x)
  x_num <- suppressWarnings(as.numeric(as.character(x[!nas])))
  
  !anyNA(x_num)
}
str2num <- function(x) {
  nas <- is.na(x)
  if (!is.numeric(x) && !is.logical(x)) x <- as.character(x)
  x_num <- suppressWarnings(as.numeric(x))
  is.na(x_num)[nas] <- TRUE
  x_num
}
trim_string <- function(x, char = " ", symmetrical = TRUE, recursive = TRUE) {
  if (is_null(x)) {
    return(x)
  }
  
  sw <- startsWith(x, char)
  ew <- endsWith(x, char)
  
  if (symmetrical) {
    while (any(sw & ew)) {
      x[sw & ew] <- gsub("^.|.$", "", x[sw & ew])
      
      if (!recursive) {
        break
      }
      
      sw <- startsWith(x, char)
      ew <- endsWith(x, char)
    }
  }
  else {
    asw <- any(sw)
    aew <- any(ew)
    
    while (asw || aew) {
      if (asw) {
        x[sw] <- gsub("^.", "", x[sw])
      }
      
      if (aew) {
        x[ew] <- gsub(".$", "", x[ew])
      }
      
      if (!recursive) {
        break
      }
      
      if (asw) {
        sw <- startsWith(x, char)
        asw <- any(sw)
      }
      
      if (aew) {
        ew <- endsWith(x, char)
        aew <- any(ew)
      }
    }
  }
  
  x
}
space <- function(n) {
  strrep(" ", n)
}
str_rev <- function(x) {
  vapply(lapply(strsplit(x, NULL), rev), paste, character(1L), collapse = "")
}

#Numbers
check_if_zero <- function(x) {
  # this is the default tolerance used in all.equal
  tolerance <- .Machine$double.eps^0.5
  abs(x) < tolerance
}
between <- function(x, range, inclusive = TRUE, na.action = FALSE) {
  if (!all(is.numeric(x))) {
    stop("'x' must be a numeric vector.", call. = FALSE)
  }
  
  if (length(range) != 2L) {
    stop("'range' must be of length 2.", call. = FALSE)
  }
  
  if (anyNA(range) || !is.numeric(range)) {
    stop("'range' must contain numeric entries only.", call. = FALSE)
  }
  
  if (range[2L] < range[1L]) {
    range <- c(range[2L], range[1L])
  }
  
  if (anyNA(x) && (length(na.action) != 1L || !is.logical(na.action))) {
    stop("'na.action' must be a logical vector of length 1.", call. = FALSE)
  }
  
  out <- {
    if (inclusive) ifelse(is.na(x), na.action, x >= range[1L] & x <= range[2L])
    else ifelse(is.na(x), na.action, x > range[1L] & x < range[2L])
  }
  
  out
}

check_if_int <- function(x) {
  #Checks if integer-like
  if (is.integer(x)) rep.int(TRUE, length(x))
  else if (is.numeric(x)) check_if_zero(x - round(x))
  else rep.int(FALSE, length(x))
}

#Statistics
binarize <- function(variable, zero = NULL, one = NULL) {
  var.name <- deparse1(substitute(variable))
  if (is.character(variable) || is.factor(variable)) {
    variable <- factor(variable, nmax = if (is.factor(variable)) nlevels(variable) else NA)
    unique.vals <- levels(variable)
  }
  else {
    unique.vals <- unique(variable)
  }
  
  if (length(unique.vals) == 1L) {
    return(setNames(rep.int(1L, length(variable)), names(variable)))
  }
  
  if (length(unique.vals) != 2L) {
    .err(sprintf("cannot binarize %s: more than two levels", var.name))
  }
  
  if (is_not_null(zero)) {
    if (!any(unique.vals == zero)) {
      .err(sprintf("the argument to `zero` is not the name of a level of %s", var.name))
    }
    
    return(setNames(as.integer(variable != zero), names(variable)))
  }
  
  if (is_not_null(one)) {
    if (!any(unique.vals == one)) {
      .err(sprintf("the argument to `one` is not the name of a level of %s", var.name))
    }
    
    return(setNames(as.integer(variable == one), names(variable)))
  }
  
  variable.numeric <- {
    if (can_str2num(unique.vals)) str2num(variable)
    else as.numeric(as.factor(variable))
  }
  
  zero <- {
    if (0 %in% variable.numeric) 0
    else min(variable.numeric, na.rm = TRUE)
  }
  
  setNames(as.integer(variable.numeric != zero), names(variable))
}
ESS <- function(w) {
  sum(w)^2 / sum(w^2)
}
center <- function(x, at = NULL, na.rm = TRUE) {
  if (is.data.frame(x)) {
    x <- as.matrix.data.frame(x)
    type <- "df"
  }
  
  if (!is.numeric(x)) {
    stop("'x' must be numeric.")
  }
  
  if (is.array(x) && length(dim(x)) > 2L) {
    stop("'x' must be a numeric or matrix-like (not array).")
  }
  
  if (is.matrix(x)) {
    type <- "matrix"
  }
  else {
    x <- matrix(x, ncol = 1)
    type <- "vec"
  }
  
  if (is_null(at)) {
    at <- colMeans(x, na.rm = na.rm)
  }
  else if (length(at) %nin% c(1L, ncol(x))) {
    stop("'at' is not the right length.")
  }
  
  out <- x - matrix(at, byrow = TRUE, ncol = ncol(x), nrow = nrow(x))
  
  switch(type,
         "df" = as.data.frame.matrix(out),
         "vec" = drop(out),
         out)
}
w.m <- function(x, w = NULL, na.rm = TRUE) {
  if (is_null(w)) w <- rep.int(1, length(x))
  if (anyNA(x)) is.na(w)[is.na(x)] <- TRUE
  
  sum(x * w, na.rm = na.rm) / sum(w, na.rm = na.rm)
}
col.w.m <- function(mat, w = NULL, na.rm = TRUE) {
  if (is_null(w)) w <- rep.int(1, nrow(mat))
  
  w.sum <- {
    if (na.rm && anyNA(mat)) colSums(w * !is.na(mat))
    else sum(w)
  }
  
  colSums(mat * w, na.rm = na.rm) / w.sum
}
col.w.v <- function(mat, w = NULL, bin.vars = NULL, na.rm = TRUE) {
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
        stop("'mat' must be a numeric matrix.")
      }
      mat <- data.matrix(mat)
    }
    else if (is.numeric(mat)) {
      mat <- matrix(mat, ncol = 1)
    }
    else {
      stop("'mat' must be a numeric matrix.")
    }
  }
  
  if (is_null(bin.vars)) {
    bin.vars <- rep.int(FALSE, ncol(mat))
  }
  else if (length(bin.vars) != ncol(mat) || anyNA(as.logical(bin.vars))) {
    stop("'bin.vars' must be a logical vector with length equal to the number of columns of 'mat'.", call. = FALSE)
  }
  
  bin.var.present <- any(bin.vars)
  non.bin.vars.present <- !all(bin.vars)
  
  var <- setNames(numeric(ncol(mat)), colnames(mat))
  if (is_null(w)) {
    if (non.bin.vars.present) {
      den <- colSums(!is.na(mat[, !bin.vars, drop = FALSE])) - 1
      var[!bin.vars] <- colSums(center(mat[, !bin.vars, drop = FALSE])^2, na.rm = na.rm) / den
    }
    if (bin.var.present) {
      means <- colMeans(mat[, bin.vars, drop = FALSE], na.rm = na.rm)
      var[bin.vars] <- means * (1 - means)
    }
  }
  else if (na.rm && anyNA(mat)) {
    # n <- nrow(mat)
    w <- array(w, dim = dim(mat))
    is.na(w)[is.na(mat)] <- TRUE
    s <- colSums(w, na.rm = na.rm)
    w <- mat_div(w, s)
    if (non.bin.vars.present) {
      x <- center(mat[, !bin.vars, drop = FALSE],
                  at = colSums(w[, !bin.vars, drop = FALSE] * mat[, !bin.vars, drop = FALSE], na.rm = na.rm))
      var[!bin.vars] <- colSums(w[, !bin.vars, drop = FALSE] * x * x, na.rm = na.rm) / (1 - colSums(w[, !bin.vars, drop = FALSE]^2, na.rm = na.rm))
    }
    if (bin.var.present) {
      means <- colSums(w[, bin.vars, drop = FALSE] * mat[, bin.vars, drop = FALSE], na.rm = na.rm)
      var[bin.vars] <- means * (1 - means)
    }
  }
  else {
    if (is_null(w)) w <- rep.int(1, nrow(mat))
    w <- w / sum(w)
    if (non.bin.vars.present) {
      x <- center(mat[, !bin.vars, drop = FALSE],
                  at = colSums(w * mat[, !bin.vars, drop = FALSE], na.rm = na.rm))
      var[!bin.vars] <- colSums(w * x * x, na.rm = na.rm) / (1 - sum(w^2))
    }
    if (bin.var.present) {
      means <- colSums(w * mat[, bin.vars, drop = FALSE], na.rm = na.rm)
      var[bin.vars] <- means * (1 - means)
    }
  }
  
  var
}
col.w.cov <- function(mat, y, w = NULL, na.rm = TRUE) {
  if (!is.matrix(mat)) {
    if (is_null(w)) {
      return(cov(mat, y, use = if (na.rm) "pair" else "everything"))
    }
    mat <- matrix(mat, ncol = 1)
  }
  
  if (is_null(w)) {
    y <- array(y, dim = dim(mat))
    if (anyNA(mat)) is.na(y)[is.na(mat)] <- TRUE
    if (anyNA(y)) is.na(mat)[is.na(y)] <- TRUE
    den <- colSums(!is.na(mat * y)) - 1
    cov <- colSums(center(mat, na.rm = na.rm) * center(y, na.rm = na.rm), na.rm = na.rm) / den
  }
  else if (na.rm && anyNA(mat)) {
    w <- array(w, dim = dim(mat))
    is.na(w)[is.na(mat)] <- TRUE
    s <- colSums(w, na.rm = na.rm)
    w <- mat_div(w, s)
    x <- w * center(mat, at = colSums(w * mat, na.rm = na.rm))
    cov <- colSums(x * y, na.rm = na.rm) / (1 - colSums(w^2, na.rm = na.rm))
  }
  else {
    w <- w / sum(w)
    x <- w * center(mat, at = colSums(w * mat, na.rm = na.rm))
    cov <- colSums(x * y, na.rm = na.rm) / (1 - sum(w^2))
  }
  
  cov
}
col.w.r <- function(mat, y, w = NULL, s.weights = NULL, bin.vars = NULL, na.rm = TRUE) {
  if (is_null(w) && is_null(s.weights)) {
    return(cor(mat, y, w, use = if (na.rm) "pair" else "everything"))
  }
  
  cov <- col.w.cov(mat, y = y, w = w, na.rm = na.rm)
  den <- sqrt(col.w.v(mat, w = s.weights, bin.vars = bin.vars, na.rm = na.rm)) *
    sqrt(col.w.v(y, w = s.weights, na.rm = na.rm))
  
  cov / den
}
.mean_abs_dev <- function(x) {
  mean_fast(abs(x - mean_fast(x, TRUE)), TRUE)
}
rms <- function(x) {
  sqrt(mean_fast(x^2))
}
.geom_mean <- function(y) {
  exp(mean_fast(log(y[is.finite(log(y))]), TRUE))
}
mat_div <- function(mat, vec) {
  mat / vec[col(mat)]
}
abs_ <- function(x, ratio = FALSE) {
  if (ratio) pmax(x, 1 / x)
  else abs(x)
}
mean_fast <- function(x, nas.possible = FALSE) {
  #Equal to mean(x, na.rm = TRUE) but faster
  #Set no.nas = FALSE if it's possible there are NAs
  if (nas.possible && anyNA(x)) {
    s <- sum(x, na.rm = TRUE)
    n <- sum(!is.na(x))
  }
  else {
    s <- sum(x)
    n <- length(x)
  }
  
  s / n
}
bw.nrd <- function(x) {
  #R's bw.nrd doesn't always work, but bw.nrd0 does
  stats::bw.nrd0(x) * 1.06 / .9
}
ave_w.m <- function(x, ..., w = NULL) {
  #ave() version of w.m() since ave() doesn't do well with multiple variables to split
  if (missing(...)) {
    x[] <- w.m(x, w)
  }
  else {
    g <- interaction(...)
    split(x, g) <- lapply(levels(g), function(i) w.m(x[g == i], w[g == i]))
  }
  
  x
}

#Formulas
subbars <- function(term) {
  if (is.name(term) || !is.language(term)) {
    return(term)
  }
  
  if (length(term) == 2L) {
    term[[2L]] <- subbars(term[[2L]])
    return(term)
  }
  
  if (is.call(term) && (term[[1L]] == as.name("|") || term[[1L]] == as.name("||"))) {
    term[[1L]] <- as.name("+")
  }
  
  for (j in seq_len(term)[-1L]) {
    term[[j]] <- subbars(term[[j]])
  }
  
  term
}

#treat/covs
assign.treat.type <- function(treat, use.multi = FALSE) {
  #Returns treat with treat.type attribute
  nunique.treat <- nunique(treat)
  
  if (nunique.treat < 2L) {
    .err("the treatment must have at least two unique values")
  }
  
  if (!use.multi && nunique.treat == 2L) {
    treat.type <- "binary"
  }
  else if (use.multi || chk::vld_character_or_factor(treat)) {
    treat.type <- "multinomial"
    if (!inherits(treat, "processed.treat")) {
      treat <- factor(treat)
    }
  }
  else {
    treat.type <- "continuous"
  }
  
  attr(treat, "treat.type") <- treat.type
  treat
}
get.treat.type <- function(treat) {
  out <- attr(treat, "treat.type")
  
  if (identical(out, "multi-category")) {
    return("multinomial")
  }
  
  out
}
has.treat.type <- function(treat) {
  is_not_null(get.treat.type(treat))
}
get_treated_level <- function(treat, estimand = NULL, focal = NULL) {
  if (is_not_null(attr(treat, "control", TRUE)) &&
      is_not_null(attr(treat, "treated", TRUE))) {
    return(attr(treat, "treated"))
  }
  
  if (is_not_null(focal)) {
    if (length(focal) > 1L || focal %nin% treat) {
      .err("the argument supplied to `focal` must be the name of a level of treatment")
    }
    
    if (is_null(estimand) || !isTRUE(estimand == "ATC")) {
      return(focal)
    }
    
    unique.vals <- {
      if (chk::vld_character_or_factor(treat))
        levels(factor(treat, nmax = 2L))
      else
        sort(unique(treat, nmax = 2L))
    }
    
    return(setdiff(unique.vals, focal))
  }
  
  treated <- attr(treat, "treated", TRUE)
  
  if (is_not_null(treated)) {
    return(treated)
  }
  
  control <- attr(treat, "control", TRUE)
  
  if (is_not_null(control)) {
    unique.vals <- {
      if (chk::vld_character_or_factor(treat))
        levels(factor(treat, nmax = 2L))
      else
        sort(unique(treat, nmax = 2L))
    }
    
    return(setdiff(unique.vals, control))
  }
  
  if (is.logical(treat)) {
    return(TRUE)
  }
  
  unique.vals <- {
    if (chk::vld_character_or_factor(treat))
      levels(factor(treat, nmax = 2L))
    else
      sort(unique(treat, nmax = 2L))
  }
  
  if (is.numeric(unique.vals) && any(unique.vals == 0)) {
    return(unique.vals[unique.vals != 0])
  }
  
  if (can_str2num(unique.vals)) {
    unique.vals.numeric <- str2num(unique.vals)
    
    if (any(unique.vals.numeric == 0)) {
      return(unique.vals[unique.vals.numeric != 0])
    }
    
    return(setdiff(unique.vals, unique.vals[which.min(unique.vals.numeric)]))
  }
  
  treated_options <- c("t", "tr", "treat", "treated", "exposed")
  t_match <- which(unique.vals %in% treated_options)
  
  if (length(t_match) == 1L) {
    return(unique.vals[t_match])
  }
  
  control_options <- c("c", "co", "ctrl", "control", "unexposed")
  c_match <- which(unique.vals %in% control_options)
  
  if (length(c_match) == 1L) {
    return(setdiff(unique.vals, unique.vals[c_match]))
  }
  
  unique.vals[2L]
}

#Uniqueness
nunique <- function(x, nmax = NA, na.rm = TRUE) {
  if (is_null(x)) {
    return(0L)
  }
  
  if (is.factor(x) && all(seq_len(nlevels(x)) %in% unclass(x))) {
    return(nlevels(x))
  }
  
  if (na.rm && anyNA(x)) {
    x <- na.rem(x)
  }
  
  length(unique(x, nmax = nmax))
}
nunique.gt <- function(x, n, na.rm = TRUE) {
  if (missing(n)) {
    stop("'n' must be supplied.")
  }
  
  if (n < 0) {
    stop("'n' must be non-negative.")
  }
  
  if (is_null(x)) {
    return(FALSE)
  }
  
  if (n == 1) {
    return(!all_the_same(x, na.rm))
  }
  
  if (length(x) < 2000L) {
    return(nunique(x, na.rm = na.rm) > n)
  }
  
  tryCatch(nunique(x, nmax = n, na.rm = na.rm) > n,
           error = function(e) TRUE)
}
all_the_same <- function(x, na.rm = TRUE) {
  if (anyNA(x)) {
    x <- na.rem(x)
    if (!na.rm) {
      return(is_null(x))
    }
  }
  
  if (is.numeric(x)) {
    return(check_if_zero(max(x) - min(x)))
  }
  
  all(x == x[1L])
}
is_binary <- function(x, na.rm = TRUE) {
  if (na.rm && anyNA(x)) x <- na.rem(x)
  !all_the_same(x) && all_the_same(x[x != x[1L]])
}
is_binary_col <- function(dat, na.rm = TRUE) {
  if (length(dim(dat)) != 2L) {
    stop("is_binary_col() cannot be used with objects that don't have 2 dimensions.")
  }
  
  apply(dat, 2L, is_binary)
}
is_0_1 <- function(x) {
  is_not_null(x) && 
    all(x == 1 | x == 0)
}

#R Processing
make_list <- function(n) {
  if (length(n) == 1L && is.numeric(n)) {
    vector("list", as.integer(n))
  }
  else if (length(n) > 0L && is.atomic(n)) {
    setNames(vector("list", length(n)), as.character(n))
  }
  else {
    stop("'n' must be an integer(ish) scalar or an atomic variable.")
  }
}
make_df <- function(ncol, nrow = 0L, types = "numeric") {
  if (missing(ncol) || is_null(ncol)) {
    ncol <- 0L
  }
  
  if (length(ncol) == 1L && is.numeric(ncol)) {
    col_names <- NULL
    ncol <- as.integer(ncol)
  }
  else if (is.atomic(ncol)) {
    col_names <- as.character(ncol)
    ncol <- length(ncol)
  }
  
  if (is_null(nrow)) {
    nrow <- 0L
  }
  
  if (length(nrow) == 1L && is.numeric(nrow)) {
    row_names <- NULL
    nrow <- as.integer(nrow)
  }
  else if (is.atomic(nrow)) {
    row_names <- as.character(nrow)
    nrow <- length(nrow)
  }
  
  df <- as.data.frame.matrix(matrix(NA_real_, nrow = nrow, ncol = ncol))
  
  names(df) <- col_names
  rownames(df) <- row_names
  
  if (is_null(types)) {
    return(df)
  }
  
  if (length(types) %nin% c(1L, ncol)) {
    stop("'types' must be equal to the number of columns.")
  }
  
  if (!is.character(types) ||
      !all(types %in% c("numeric", "integer", "logical", "character", NA))) {
    stop("'types' must be an acceptable type. For factors, use NA.")
  }
  
  if (length(types) == 1L) {
    types <- rep.int(types, ncol)
  }
  
  for (i in which(!is.na(types))) {
    if (types[i] != "numeric") {
      df[[i]] <- get(types[i])(nrow)
    }
  }
  
  df
}
rep_with <- function(x, y) {
  #Helper function to fill named vectors with x and given names of y
  setNames(rep.int(x, length(y)), names(y))
}
ifelse_ <- function(...) {
  dotlen <- ...length()
  if (dotlen %% 2 == 0) {
    stop("`ifelse_()` must have an odd number of arguments: pairs of test/yes, and one no.")
  }
  
  out <- ...elt(dotlen)
  
  if (dotlen <= 1 && !is.atomic(out)) {
    stop("The first entry to `ifelse_()` must be atomic.")
  }
  
  if (!is.atomic(out)) {
    stop("The last entry to `ifelse_()` must be atomic.")
  }
  
  if (length(out) == 1L) {
    out <- rep.int(out, length(..1))
  }
  
  n <- length(out)
  
  for (i in seq_len((dotlen - 1) / 2)) {
    test <- ...elt(2 * i - 1)
    yes <- ...elt(2 * i)
    if (length(yes) == 1) yes <- rep.int(yes, n)
    if (length(yes) != n || length(test) != n) stop("All entries must have the same length.")
    if (!is.logical(test)) stop(sprintf("The %s entry to `ifelse_()` must be logical.", ordinal(2 * i - 1)))
    if (!is.atomic(yes)) stop(sprintf("The %s entry to `ifelse_()` must be atomic.", ordinal(2 * i)))
    pos <- which(test)
    out[pos] <- yes[pos]
  }
  
  out
}
is_ <- function(x, types, stop = FALSE, arg.to = FALSE) {
  s1 <- deparse1(substitute(x))
  
  if (is_not_null(x)) {
    for (i in types) {
      if (i == "list") {
        it.is <- is.list(x) && !is.data.frame(x)
      }
      else if (is_not_null(get0(paste0("is_", i)))) {
        it.is <- get0(paste0("is_", i))(x)
      }
      else if (is_not_null(get0(paste.("is", i)))) {
        it.is <- get0(paste.("is", i))(x)
      }
      else {
        it.is <- inherits(x, i)
      }
      
      if (it.is) break
    }
  }
  else {
    it.is <- FALSE
  }
  
  if (stop && !it.is) {
    s0 <- if (arg.to) "the argument to " else ""
    s2 <- if (any(types %in% c("factor", "character", "numeric", "logical"))) "vector" else ""
    
    .err(sprintf("%s'%s' must be a %s %s",
                 s0, s1, word_list(types, and.or = "or"), s2))
  }
  
  it.is
}
is_mat_like <- function(x) {
  length(dim(x)) == 2L
}
is_null <- function(x) length(x) == 0L
is_not_null <- function(x) !is_null(x)
if_null_then <- function(x1 = NULL, x2 = NULL, ...) {
  if (is_not_null(x1)) {
    return(x1)
  }
  
  if (is_not_null(x2)) {
    return(x2)
  }
  
  if (...length() > 0L) {
    for (k in seq_len(...length())) {
      if (is_not_null(...elt(k))) {
        return(...elt(k))
      }
    }
  }
  
  x1
}
clear_null <- function(x) {
  x[vapply(x, function(i) {
    is_null(i) &&
      (is_null(dim(i)) || all(dim(i) == 0L))
  }, logical(1L))] <- NULL
  x
}
clear_attr <- function(x, all = FALSE) {
  if (all) {
    attributes(x) <- NULL
  }
  else {
    dont_clear <- c("names", "class", "dim", "dimnames", "row.names")
    attributes(x)[names(attributes(x)) %nin% dont_clear] <- NULL
  }
  x
}
probably.a.bug <- function() {
  fun <- paste(deparse(sys.call(-1L)), collapse = "\n")
  .err(paste0("An error was produced and is likely a bug. Please let the maintainer know a bug was produced by the function\n",
              fun), tidy = FALSE)
}
`%nin%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))
`%pin%` <- function(x, table) {
  #Partial in. TRUE if x uniquely identifies values in table.
  !is.na(pmatch(x, table))
}
`%cin%` <- function(x, table) {
  #Partial in w/ charmatch. TRUE if x at all in table.
  !is.na(charmatch(x, table))
}
is_error <- function(x) {inherits(x, "try-error")}
null_or_error <- function(x) {is_null(x) || is_error(x)}

#Extract variables from ..., similar to ...elt(), by name without evaluating list(...)
...get <- function(x, ifnotfound = NULL) {
  expr <- quote({
    .m1 <- match(.x, ...names())
    if (anyNA(.m1)) {
      .ifnotfound
    }
    else {
      .m2 <- ...elt(.m1[1L])
      if (is_not_null(.m2)) .m2
      else .ifnotfound
    }
  })
  
  eval(expr,
       pairlist(.x = x[1L], .ifnotfound = ifnotfound),
       parent.frame(1L))
}
...mget <- function(x) {
  found <- match(x, eval(quote(...names()), parent.frame(1L)))
  
  not_found <- is.na(found)
  
  if (all(not_found)) {
    return(list())
  }
  
  setNames(lapply(found[!not_found], function(z) {
    eval(quote(...elt(.z)),
         pairlist(.z = z),
         parent.frame(3L))
  }), x[!not_found])
}

#More informative and cleaner version of base::match.arg(). Uses chk.
match_arg <- function(arg, choices, several.ok = FALSE) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg)) {
    stop("No argument was supplied to match_arg.")
  }
  
  arg.name <- deparse1(substitute(arg), width.cutoff = 500L)
  
  if (missing(choices)) {
    sysP <- sys.parent()
    formal.args <- formals(sys.function(sysP))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }
  
  if (is_null(arg)) {
    return(choices[1L])
  }
  
  if (several.ok) {
    chk::chk_character(arg, x_name = add_quotes(arg.name, "`"))
  }
  else {
    chk::chk_string(arg, x_name = add_quotes(arg.name, "`"))
    if (identical(arg, choices)) {
      return(arg[1L])
    }
  }
  
  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    .err(sprintf("the argument to `%s` should be %s%s",
                 arg.name,
                 ngettext(length(choices), "", if (several.ok) "at least one of " else "one of "),
                 word_list(choices, and.or = "or", quotes = 2L)))
  i <- i[i > 0L]
  
  choices[i]
}

grab <- function(x, what) {
  lapply(x, function(z) z[[what]])
}
last <- function(x) {
  x[[length(x)]]
}
`last<-` <- function(x, value) {
  `[[<-`(x, length(x), value)
}
len <- function(x, recursive = TRUE) {
  if (is_null(x)) 0L
  else if (length(dim(x)) > 1) NROW(x)
  else if (is.list(x) && recursive) vapply(x, len, numeric(1L), recursive = FALSE)
  else length(x)
}
seq_row <- function(x) {
  if (is_null(x)) {
    return(integer(0L))
  }
  
  if (length(dim(x)) != 2L) {
    stop("dim(x) must have length 2")
  }
  
  seq_len(nrow(x))
}
seq_col <- function(x) {
  if (is_null(x)) {
    return(integer(0L))
  }
  
  if (length(dim(x)) != 2L) {
    stop("dim(x) must have length 2")
  }
  
  seq_len(ncol(x))
}
na.rem <- function(x) {
  #A faster na.omit for vectors
  x[!is.na(x)]
}
anyNA_col <- function(x) {
  colSums(is.na(x)) > 0L
}
check_if_call_from_fun <- function(fun) {
  # Check if called from within function f
  if (missing(fun) || !exists(deparse1(substitute(fun)), mode = "function")) {
    return(FALSE)
  }
  
  for (x in sys.parents()) {
    if (identical(fun, sys.function(x), ignore.environment = TRUE)) {
      return(TRUE)
    }
  }
  
  FALSE
}
has_method <- function(class, fun) {
  if (!is.character(fun) || length(fun) != 1) {
    stop("'fun' must be a string of length 1.")
  }
  
  if (!is.character(class)) {
    stop("'class' must be a character vector.")
  }
  
  vapply(class, function(cl) is_not_null(utils::getS3method(fun, cl, optional = TRUE,
                                                            envir = asNamespace(utils::packageName()))),
         logical(1L))
}
set_class <- function(x, .class, .replace = TRUE, .last = TRUE) {
  if (missing(.class)) {
    return(x)
  }
  
  if (.replace) {
    class(x) <- .class
  }
  else if (.last) {
    class(x) <- c(class(x), .class)
  }
  else {
    class(x) <- c(.class, class(x))
  }
  
  x
}

#Efficient versions of any(vapply(...)) and all(vapply(...))
any_apply <- function(X, FUN, ...) {
  FUN <- match.fun(FUN)
  if (!is.vector(X) || is.object(X)) {
    X <- as.list(X)
  }
  
  for (x in X) {
    if (isTRUE(FUN(x, ...))) {
      return(TRUE)
    }
  }
  
  FALSE
}
all_apply <- function(X, FUN, ...) {
  FUN <- match.fun(FUN)
  if (!is.vector(X) || is.object(X)) {
    X <- as.list(X)
  }
  
  for (x in X) {
    if (isFALSE(FUN(x, ...))) {
      return(FALSE)
    }
  }
  
  TRUE
}
