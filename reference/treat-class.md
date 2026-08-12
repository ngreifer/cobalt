# The `treat` Class

cobalt represents a processed treatment variable as an object of class
`treat`, carrying attributes that record what kind of treatment it is
and what its groups are called. Every
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.md)
method builds one, and every balance computation reads it rather than
re-deriving the same facts from the values.

The attributes are

- `treat.type`:

  one of `"binary"`, `"multinomial"`, `"continuous"`, or `"censoring"`.
  This is what decides which balance statistics apply and how the
  treatment is compared.

- `treat.name`:

  the name the treatment was given, if it had one.

- `treat_names`:

  for a binary treatment, the display names of the two groups, named
  `control` and `treated`; for a multi-category treatment, one name per
  level; for a censoring indicator, `Uncensored` and `Censored`. `NULL`
  for a continuous treatment.

- `treat_vals`:

  the values in the data that those names correspond to.

- `group_labels`:

  how the treatment's groups are labeled in the column names of a
  balance table. `c("0", "1")` for a binary treatment, which is cobalt's
  positional convention (`M.0` is the control group's mean); otherwise
  whatever names the groups have.

`[` preserves all of them, so a subset of a `treat` is still a `treat`.

[`.cens()`](https://ngreifer.github.io/cobalt/reference/cens.md) is the
only way to build one directly; it returns a censoring indicator.
[`WeightIt::.cens()`](https://ngreifer.github.io/WeightIt/reference/dot-cens.html)
returns an object with the same class and a compatible `treat.type`, so
an indicator tagged by either package is accepted here.

The class is shared, so carrying it means only "this is a treatment",
not "this carries the attributes above": WeightIt tags its own
treatments with it and supplies `treat.type` alone, and even cobalt's
[`.cens()`](https://ngreifer.github.io/cobalt/reference/cens.md) has no
group names until the treatment is processed. Anything that reads
`treat_names` or `treat_vals` therefore has to check for them rather
than for the class, and process the treatment when they are absent.

## See also

- [`.cens()`](https://ngreifer.github.io/cobalt/reference/cens.md)

- [`class-bal.tab.cens`](https://ngreifer.github.io/cobalt/reference/class-bal.tab.cens.md)
