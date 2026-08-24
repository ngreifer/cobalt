## Test environments

* local macOS (aarch64), R 4.6.1
* win-builder, R-devel
* R-hub

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

One reverse dependency changes to worse. `mvGPS` gains a warning when it is
installed:

```
Warning: replacing previous import 'WeightIt::.cens' by 'cobalt::.cens'
  when loading 'mvGPS'
```

`mvGPS` imports both *WeightIt* and *cobalt* in full, and this version of
*cobalt* exports `.cens()` for the first time. It is the only name the two
packages export in common, so this is the first time the two wholesale imports
have collided.

`mvGPS` does not call `.cens()` anywhere, so which of the two functions it ends
up bound to makes no difference to it, and it installs and checks otherwise as
before. I have notified its maintainer.

The warning will resolve on its own shortly. *WeightIt*, which I also maintain,
is being updated to re-export `cobalt::.cens()` rather than define its own copy,
after which the two bindings are the same object and R does not warn. That
update cannot be submitted until this version of *cobalt* is on CRAN, since it
imports from it.

## Method references

The methods implemented are described in the references given on the help pages
of the functions that implement them.
