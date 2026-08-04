#Shared body of every `bal.tab()` method: gather the method's own arguments and
#`...`, hand them to the named `x2base` method, and dispatch on the resulting `X`.
#
#`.x2base` is a method name rather than the `x2base` generic. Several `bal.tab`
#methods are aliases across classes -- `bal.tab.matrix` is `bal.tab.data.frame`,
#and there is no `x2base.matrix` -- so generic dispatch would land on
#`x2base.default` for those.
#
#Both formals are dot-prefixed: `...` carries the user's own arguments, and a named
#element of `...` matches a formal of the same name even when the caller supplied
#that formal positionally.
.bal.tab_dispatch <- function(.x, .x2base, .env = parent.frame()) {
  #The caller's *matched* call, so that an abbreviated argument name is normalized
  #before being compared with `names(args)` below.
  .call <- match.call(definition = sys.function(-1L), call = sys.call(-1L))

  args <- try_arg(c(as.list(.env), eval(quote(list(...)), .env)))
  args[["x"]] <- NULL

  #A formal the caller left missing has no value to pass on. One whose value is
  #zero-length is dropped as well, unless the caller named it explicitly.
  args[vapply(args, rlang::is_missing, logical(1L))] <- NULL
  args[lengths(args) == 0L & names(args) %nin% names(.call)[-1L]] <- NULL

  X <- do.call(.x2base, c(list(.x), args), quote = TRUE)

  #Everything `x2base()` consumed is removed: `base.bal.tab()` receives those
  #values inside `X`, and leaving them in `args` too would make `do.call()` reject
  #the call for matching a formal twice. `names(X)` is the full slot list including
  #the empty slots, which is what makes this complete.
  args[names(X)] <- NULL

  do.call("base.bal.tab", c(list(.assign_X_class(X)), args), quote = TRUE)
}
