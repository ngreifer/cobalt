# #' @exportS3Method as.data.frame bal.tab
as.data.frame.bal.tab <- function(x, row.names = NULL, optional = FALSE, agg.fun = NULL, var.names = NULL, ...) {
    
    #Regular bal.tab
    
    B <- x[["Balance"]]
    B[["variable"]] <- rownames(B)
    
    p <- attr(x, "print.options")
    
    B_long <- reshape(B, direction = "long", 
                      idvar = c("variable", "Type"),
                      varying = list(Diff = c("Diff.Un", "Diff.Adj"),
                                     KS = c("KS.Un", "KS.Adj")),
                      v.names = c("Diff", "KS"),
                      timevar = "sample", times = c("Un", "Adj")) |>
        reshape(direction = "long",
                varying = c("Diff", "KS"),
                v.names = "value",
                timevar = "statistic",
                times = c("Diff", "KS"))
}