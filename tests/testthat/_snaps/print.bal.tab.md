# error messages from print are stable

    Code
      print(b, disp = "bogus")
    Condition
      Error:
      ! "bogus" is not allowed in `disp`.

---

    Code
      print(b, stats = "ks.statistics")
    Condition
      Error:
      ! `stats` cannot contain "ks.statistics" when it was not requested in the original call to `bal.tab()`.

---

    Code
      print(b, un = "yes")
    Condition
      Error:
      ! `un` must be a logical value (TRUE or FALSE).

