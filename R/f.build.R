f.build <- function(y, rhs) {
  f<-as.formula(paste(deparse(substitute(y))," ~ ",paste(names(rhs),collapse="+")))
  return(f)
}
