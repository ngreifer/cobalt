# No adjustment

    Code
      bal.tab(covs, treat = lalonde$treat, m.threshold = 0.1, v.threshold = 2,
      imbalanced.only = T)
    Message
      Note: 's.d.denom' not specified; assuming pooled.
    Output
      Balance Measures
                    Type Diff.Un     M.Threshold.Un V.Ratio.Un   V.Threshold.Un
      age        Contin. -0.2419 Not Balanced, >0.1     0.4400 Not Balanced, >2
      educ       Contin.  0.0448     Balanced, <0.1     0.4959 Not Balanced, >2
      race_black  Binary  0.6404 Not Balanced, >0.1          .                 
      race_white  Binary -0.5577 Not Balanced, >0.1          .                 
      married     Binary -0.3236 Not Balanced, >0.1          .                 
      nodegree    Binary  0.1114 Not Balanced, >0.1          .                 
      re74       Contin. -0.5958 Not Balanced, >0.1     0.5181     Balanced, <2
      re75       Contin. -0.2870 Not Balanced, >0.1     0.9563     Balanced, <2
      
      Balance tally for mean differences
                         count
      Balanced, <0.1         2
      Not Balanced, >0.1     7
      
      Variable with the greatest mean difference
         Variable Diff.Un     M.Threshold.Un
       race_black  0.6404 Not Balanced, >0.1
      
      Balance tally for variance ratios
                       count
      Balanced, <2         2
      Not Balanced, >2     2
      
      Variable with the greatest variance ratio
       Variable V.Ratio.Un   V.Threshold.Un
            age       0.44 Not Balanced, >2
      
      Sample sizes
          Control Treated
      All     429     185

---

    Code
      bal.tab(f.build("treat", covs), data = lalonde, ks.threshold = 0.1)
    Message
      Note: 's.d.denom' not specified; assuming pooled.
    Output
      Balance Measures
                     Type Diff.Un  KS.Un    KS.Threshold.Un
      age         Contin. -0.2419 0.1577 Not Balanced, >0.1
      educ        Contin.  0.0448 0.1114 Not Balanced, >0.1
      race_black   Binary  0.6404 0.6404 Not Balanced, >0.1
      race_hispan  Binary -0.0827 0.0827     Balanced, <0.1
      race_white   Binary -0.5577 0.5577 Not Balanced, >0.1
      married      Binary -0.3236 0.3236 Not Balanced, >0.1
      nodegree     Binary  0.1114 0.1114 Not Balanced, >0.1
      re74        Contin. -0.5958 0.4470 Not Balanced, >0.1
      re75        Contin. -0.2870 0.2876 Not Balanced, >0.1
      
      Balance tally for KS statistics
                         count
      Balanced, <0.1         1
      Not Balanced, >0.1     8
      
      Variable with the greatest KS statistic
         Variable  KS.Un    KS.Threshold.Un
       race_black 0.6404 Not Balanced, >0.1
      
      Sample sizes
          Control Treated
      All     429     185

# MatchIt with PS

    Code
      bal.tab(m1, int = T, v.threshold = 2, imbalanced.only = F)
    Output
      Call
       MatchIt::matchit(formula = treat ~ log(age) * married + educ + 
          race + nodegree + re74 + re75, data = lalonde, discard = "both", 
          reestimate = TRUE, replace = T, ratio = 2)
      
      Balance Measures
                                   Type Diff.Adj V.Ratio.Adj      V.Threshold
      distance                 Distance  -0.0004      0.9889     Balanced, <2
      log(age)                  Contin.   0.0849      0.4519 Not Balanced, >2
      married                    Binary   0.0000           .                 
      educ                      Contin.  -0.0221      0.5851     Balanced, <2
      race_black                 Binary  -0.0083           .                 
      race_hispan                Binary  -0.0250           .                 
      race_white                 Binary   0.0333           .                 
      nodegree                   Binary   0.0472           .                 
      re74                      Contin.  -0.0548      1.1790     Balanced, <2
      re75                      Contin.  -0.0511      1.0739     Balanced, <2
      log(age) * married_0      Contin.   0.0079      0.9589     Balanced, <2
      log(age) * married_1      Contin.   0.0089      1.0216     Balanced, <2
      log(age) * educ           Contin.   0.0456      0.6614     Balanced, <2
      log(age) * race_black     Contin.  -0.0202      0.9854     Balanced, <2
      log(age) * race_hispan    Contin.  -0.1001      0.7395     Balanced, <2
      log(age) * race_white     Contin.   0.1244      1.6089     Balanced, <2
      log(age) * nodegree_0     Contin.  -0.0970      0.9298     Balanced, <2
      log(age) * nodegree_1     Contin.   0.1131      0.9089     Balanced, <2
      log(age) * re74           Contin.  -0.0679      1.0360     Balanced, <2
      log(age) * re75           Contin.  -0.0356      1.1955     Balanced, <2
      married_0 * educ          Contin.   0.0050      0.8862     Balanced, <2
      married_0 * race_black     Binary  -0.0194           .                 
      married_0 * race_hispan    Binary  -0.0194           .                 
      married_0 * race_white     Binary   0.0389           .                 
      married_0 * nodegree_0     Binary   0.0000           .                 
      married_0 * nodegree_1     Binary   0.0000           .                 
      married_0 * re74          Contin.  -0.0291      0.8570     Balanced, <2
      married_0 * re75          Contin.  -0.0675      0.6186     Balanced, <2
      married_1 * educ          Contin.  -0.0160      0.9128     Balanced, <2
      married_1 * race_black     Binary   0.0111           .                 
      married_1 * race_hispan    Binary  -0.0056           .                 
      married_1 * race_white     Binary  -0.0056           .                 
      married_1 * nodegree_0     Binary  -0.0472           .                 
      married_1 * nodegree_1     Binary   0.0472           .                 
      married_1 * re74          Contin.  -0.0451      1.6016     Balanced, <2
      married_1 * re75          Contin.  -0.0145      1.3939     Balanced, <2
      educ * race_black         Contin.  -0.0317      0.8896     Balanced, <2
      educ * race_hispan        Contin.  -0.1085      0.6776     Balanced, <2
      educ * race_white         Contin.   0.1031      1.3348     Balanced, <2
      educ * nodegree_0         Contin.  -0.1320      0.8414     Balanced, <2
      educ * nodegree_1         Contin.   0.1560      0.9629     Balanced, <2
      educ * re74               Contin.  -0.0835      1.0430     Balanced, <2
      educ * re75               Contin.  -0.1052      0.8354     Balanced, <2
      race_black * nodegree_0    Binary  -0.0528           .                 
      race_black * nodegree_1    Binary   0.0444           .                 
      race_black * re74         Contin.  -0.0222      1.2171     Balanced, <2
      race_black * re75         Contin.  -0.0675      1.0063     Balanced, <2
      race_hispan * nodegree_0   Binary  -0.0139           .                 
      race_hispan * nodegree_1   Binary  -0.0111           .                 
      race_hispan * re74        Contin.  -0.0951      0.6470     Balanced, <2
      race_hispan * re75        Contin.  -0.0097      1.3203     Balanced, <2
      race_white * nodegree_0    Binary   0.0194           .                 
      race_white * nodegree_1    Binary   0.0139           .                 
      race_white * re74         Contin.  -0.0525      0.6864     Balanced, <2
      race_white * re75         Contin.   0.0673      2.8908 Not Balanced, >2
      nodegree_0 * re74         Contin.  -0.0674      1.2444     Balanced, <2
      nodegree_0 * re75         Contin.  -0.3959      0.3324 Not Balanced, >2
      nodegree_1 * re74         Contin.  -0.0030      1.0111     Balanced, <2
      nodegree_1 * re75         Contin.   0.1559      2.6239 Not Balanced, >2
      re74 * re75               Contin.   0.0520      2.9528 Not Balanced, >2
      
      Balance tally for variance ratios
                       count
      Balanced, <2        34
      Not Balanced, >2     5
      
      Variable with the greatest variance ratio
                Variable V.Ratio.Adj      V.Threshold
       nodegree_0 * re75      0.3324 Not Balanced, >2
      
      Sample sizes
                           Control Treated
      All                   429.       185
      Matched (ESS)          67.71     180
      Matched (Unweighted)  118.       180
      Unmatched             239.         0
      Discarded              72.         5

# MatchIt with distance

    Code
      bal.tab(m1.1, data = lalonde)
    Output
      Call
       MatchIt::matchit(formula = f.build("treat", covs), data = lalonde, 
          distance = runif(nrow(lalonde)))
      
      Balance Measures
                      Type Diff.Adj
      distance    Distance   0.0043
      age          Contin.  -0.2553
      educ         Contin.   0.0511
      race_black    Binary   0.5946
      race_hispan   Binary  -0.0486
      race_white    Binary  -0.5459
      married       Binary  -0.2973
      nodegree      Binary   0.1243
      re74         Contin.  -0.7496
      re75         Contin.  -0.2864
      
      Sample sizes
                Control Treated
      All           429     185
      Matched       185     185
      Unmatched     244       0

# MatchIt with complicated formula

    Code
      bal.tab(m1.2, addl = ~ rms::rcs(educ, 3) + poly(age, 2))
    Output
      Call
       MatchIt::matchit(formula = treat ~ log(age) * married + poly(educ, 
          2) + rms::rcs(age, 3) + race + factor(nodegree) + re74 + 
          re75, data = lalonde, replace = T, ratio = 2)
      
      Balance Measures
                           Type Diff.Adj
      distance         Distance   0.0014
      log(age)          Contin.   0.0496
      married            Binary   0.0757
      poly(educ, 2)_1   Contin.   0.0094
      poly(educ, 2)_2   Contin.   0.0613
      age               Contin.  -0.0011
      age'              Contin.  -0.1261
      race_black         Binary  -0.0135
      race_hispan        Binary  -0.0081
      race_white         Binary   0.0216
      factor(nodegree)   Binary   0.0108
      re74              Contin.   0.1383
      re75              Contin.   0.1371
      educ'             Contin.   0.0102
      poly(age, 2)_2    Contin.  -0.2459
      
      Sample sizes
                           Control Treated
      All                   429.       185
      Matched (ESS)          36.88     185
      Matched (Unweighted)  108.       185
      Unmatched             321.         0

# MatchIt with Malahanobis

    Code
      bal.tab(m2, data = lalonde, int = T, v.threshold = 2, addl = "race")
    Output
      Call
       MatchIt::matchit(formula = f.build("treat", covs), data = lalonde, 
          distance = "mahalanobis")
      
      Balance Measures
                                  Type Diff.Adj V.Ratio.Adj      V.Threshold
      age                      Contin.   0.1254      0.5333     Balanced, <2
      educ                     Contin.  -0.0457      0.8306     Balanced, <2
      race_black                Binary   0.3784           .                 
      race_hispan               Binary   0.0000           .                 
      race_white                Binary  -0.3784           .                 
      married                   Binary  -0.0595           .                 
      nodegree                  Binary   0.0486           .                 
      re74                     Contin.  -0.2482      0.8710     Balanced, <2
      re75                     Contin.  -0.1331      1.0081     Balanced, <2
      age * educ               Contin.   0.0776      0.6363     Balanced, <2
      age * race_black         Contin.   0.8498      0.6054     Balanced, <2
      age * race_hispan        Contin.   0.0118      1.1057     Balanced, <2
      age * race_white         Contin.  -1.1167      0.3474 Not Balanced, >2
      age * married_0          Contin.   0.2640      0.8431     Balanced, <2
      age * married_1          Contin.  -0.1814      0.7026     Balanced, <2
      age * nodegree_0         Contin.  -0.0995      0.8893     Balanced, <2
      age * nodegree_1         Contin.   0.1654      0.8485     Balanced, <2
      age * re74               Contin.  -0.3121      0.5767     Balanced, <2
      age * re75               Contin.  -0.1297      1.0007     Balanced, <2
      educ * race_black        Contin.   0.9414      0.6058     Balanced, <2
      educ * race_hispan       Contin.  -0.0023      0.9786     Balanced, <2
      educ * race_white        Contin.  -1.2219      0.3613 Not Balanced, >2
      educ * married_0         Contin.   0.1279      0.8115     Balanced, <2
      educ * married_1         Contin.  -0.1580      0.7993     Balanced, <2
      educ * nodegree_0        Contin.  -0.1147      0.8939     Balanced, <2
      educ * nodegree_1        Contin.   0.1238      0.9564     Balanced, <2
      educ * re74              Contin.  -0.2464      0.8729     Balanced, <2
      educ * re75              Contin.  -0.1636      0.9095     Balanced, <2
      race_black * married_0    Binary   0.3514           .                 
      race_black * married_1    Binary   0.0270           .                 
      race_black * nodegree_0   Binary   0.0649           .                 
      race_black * nodegree_1   Binary   0.3135           .                 
      race_black * re74        Contin.   0.0738      1.4717     Balanced, <2
      race_black * re75        Contin.   0.1283      1.6861     Balanced, <2
      race_hispan * married_0   Binary   0.0000           .                 
      race_hispan * married_1   Binary   0.0000           .                 
      race_hispan * nodegree_0  Binary   0.0000           .                 
      race_hispan * nodegree_1  Binary   0.0000           .                 
      race_hispan * re74       Contin.   0.0125      2.3313 Not Balanced, >2
      race_hispan * re75       Contin.   0.0429      1.3644     Balanced, <2
      race_white * married_0    Binary  -0.2919           .                 
      race_white * married_1    Binary  -0.0865           .                 
      race_white * nodegree_0   Binary  -0.1135           .                 
      race_white * nodegree_1   Binary  -0.2649           .                 
      race_white * re74        Contin.  -1.7053      0.0495 Not Balanced, >2
      race_white * re75        Contin.  -1.1031      0.1021 Not Balanced, >2
      married_0 * nodegree_0    Binary   0.0054           .                 
      married_0 * nodegree_1    Binary   0.0541           .                 
      married_0 * re74         Contin.  -0.1213      0.8948     Balanced, <2
      married_0 * re75         Contin.  -0.0188      1.0963     Balanced, <2
      married_1 * nodegree_0    Binary  -0.0541           .                 
      married_1 * nodegree_1    Binary  -0.0054           .                 
      married_1 * re74         Contin.  -0.2140      0.7090     Balanced, <2
      married_1 * re75         Contin.  -0.1377      0.8936     Balanced, <2
      nodegree_0 * re74        Contin.  -0.2048      0.7711     Balanced, <2
      nodegree_0 * re75        Contin.  -0.3461      0.3803 Not Balanced, >2
      nodegree_1 * re74        Contin.  -0.1269      0.8231     Balanced, <2
      nodegree_1 * re75        Contin.   0.0401      1.5570     Balanced, <2
      re74 * re75              Contin.  -0.0904      1.1162     Balanced, <2
      
      Balance tally for variance ratios
                       count
      Balanced, <2        32
      Not Balanced, >2     6
      
      Variable with the greatest variance ratio
                Variable V.Ratio.Adj      V.Threshold
       race_white * re74      0.0495 Not Balanced, >2
      
      Sample sizes
                Control Treated
      All           429     185
      Matched       185     185
      Unmatched     244       0

# MatchIt with subclassification

    Code
      bal.tab(m3, int = F, v.threshold = 2, disp.subclass = T, ks.threshold = 0.1)
    Output
      Call
       MatchIt::matchit(formula = f.build("treat", covs), data = lalonde, 
          method = "subclass")
      
      Balance by subclass
       - - - Subclass 1 - - - 
                      Type Diff.Adj V.Ratio.Adj  V.Threshold KS.Adj
      distance    Distance   0.2785      1.5087              0.4866
      age          Contin.  -0.4024      0.4163 Balanced, <2 0.1600
      educ         Contin.   0.1142      0.5625 Balanced, <2 0.0959
      race_black    Binary   0.0823           .              0.0823
      race_hispan   Binary   0.1492           .              0.1492
      race_white    Binary  -0.2315           .              0.2315
      married       Binary  -0.2877           .              0.2877
      nodegree      Binary  -0.0003           .              0.0003
      re74         Contin.  -0.5864      1.1631 Balanced, <2 0.4293
      re75         Contin.  -0.1729      1.0676 Balanced, <2 0.2623
                        KS.Threshold
      distance                      
      age         Not Balanced, >0.1
      educ            Balanced, <0.1
      race_black      Balanced, <0.1
      race_hispan Not Balanced, >0.1
      race_white  Not Balanced, >0.1
      married     Not Balanced, >0.1
      nodegree        Balanced, <0.1
      re74        Not Balanced, >0.1
      re75        Not Balanced, >0.1
      
       - - - Subclass 2 - - - 
                      Type Diff.Adj V.Ratio.Adj  V.Threshold KS.Adj
      distance    Distance   0.1873      1.3761              0.3858
      age          Contin.  -0.7473      0.3404 Balanced, <2 0.4220
      educ         Contin.   0.1183      0.6596 Balanced, <2 0.1035
      race_black    Binary   0.0094           .              0.0094
      race_hispan   Binary  -0.0094           .              0.0094
      race_white    Binary   0.0000           .              0.0000
      married       Binary  -0.2473           .              0.2473
      nodegree      Binary  -0.0121           .              0.0121
      re74         Contin.  -0.0352      0.9023 Balanced, <2 0.1882
      re75         Contin.  -0.0970      0.9876 Balanced, <2 0.2325
                        KS.Threshold
      distance                      
      age         Not Balanced, >0.1
      educ        Not Balanced, >0.1
      race_black      Balanced, <0.1
      race_hispan     Balanced, <0.1
      race_white      Balanced, <0.1
      married     Not Balanced, >0.1
      nodegree        Balanced, <0.1
      re74        Not Balanced, >0.1
      re75        Not Balanced, >0.1
      
       - - - Subclass 3 - - - 
                      Type Diff.Adj V.Ratio.Adj  V.Threshold KS.Adj
      distance    Distance  -0.0140      1.0943              0.1197
      age          Contin.   0.0524      0.4654 Balanced, <2 0.2677
      educ         Contin.   0.1372      0.3845 Balanced, <2 0.2191
      race_black    Binary   0.0000           .              0.0000
      race_hispan   Binary   0.0000           .              0.0000
      race_white    Binary   0.0000           .              0.0000
      married       Binary   0.3550           .              0.3550
      nodegree      Binary   0.2191           .              0.2191
      re74         Contin.  -0.2669      0.2635 Balanced, <2 0.3225
      re75         Contin.  -0.0970      0.6585 Balanced, <2 0.1460
                        KS.Threshold
      distance                      
      age         Not Balanced, >0.1
      educ        Not Balanced, >0.1
      race_black      Balanced, <0.1
      race_hispan     Balanced, <0.1
      race_white      Balanced, <0.1
      married     Not Balanced, >0.1
      nodegree    Not Balanced, >0.1
      re74        Not Balanced, >0.1
      re75        Not Balanced, >0.1
      
       - - - Subclass 4 - - - 
                      Type Diff.Adj V.Ratio.Adj      V.Threshold KS.Adj
      distance    Distance  -0.0003      1.4231                  0.1860
      age          Contin.  -0.0499      0.2048     Balanced, <2 0.2902
      educ         Contin.  -0.1436      1.0063     Balanced, <2 0.1711
      race_black    Binary   0.0000           .                  0.0000
      race_hispan   Binary   0.0000           .                  0.0000
      race_white    Binary   0.0000           .                  0.0000
      married       Binary  -0.1116           .                  0.1116
      nodegree      Binary  -0.0417           .                  0.0417
      re74         Contin.  -0.0073      2.2051 Not Balanced, >2 0.3051
      re75         Contin.  -0.0801      1.7074     Balanced, <2 0.2902
                        KS.Threshold
      distance                      
      age         Not Balanced, >0.1
      educ        Not Balanced, >0.1
      race_black      Balanced, <0.1
      race_hispan     Balanced, <0.1
      race_white      Balanced, <0.1
      married     Not Balanced, >0.1
      nodegree        Balanced, <0.1
      re74        Not Balanced, >0.1
      re75        Not Balanced, >0.1
      
       - - - Subclass 5 - - - 
                      Type Diff.Adj V.Ratio.Adj      V.Threshold KS.Adj
      distance    Distance  -0.0224      0.7370                  0.2832
      age          Contin.   0.2640      0.7413     Balanced, <2 0.2885
      educ         Contin.  -0.2977      0.5077     Balanced, <2 0.1774
      race_black    Binary   0.0000           .                  0.0000
      race_hispan   Binary   0.0000           .                  0.0000
      race_white    Binary   0.0000           .                  0.0000
      married       Binary   0.0000           .                  0.0000
      nodegree      Binary   0.0376           .                  0.0376
      re74         Contin.   0.0190      2.3025 Not Balanced, >2 0.1398
      re75         Contin.   0.1233      4.2576 Not Balanced, >2 0.1613
                        KS.Threshold
      distance                      
      age         Not Balanced, >0.1
      educ        Not Balanced, >0.1
      race_black      Balanced, <0.1
      race_hispan     Balanced, <0.1
      race_white      Balanced, <0.1
      married         Balanced, <0.1
      nodegree        Balanced, <0.1
      re74        Not Balanced, >0.1
      re75        Not Balanced, >0.1
      
       - - - Subclass 6 - - - 
                      Type Diff.Adj V.Ratio.Adj      V.Threshold KS.Adj
      distance    Distance   0.0143      2.8080                  0.3441
      age          Contin.   0.5245      0.4039     Balanced, <2 0.6022
      educ         Contin.   0.2781      6.1419 Not Balanced, >2 0.1398
      race_black    Binary   0.0000           .                  0.0000
      race_hispan   Binary   0.0000           .                  0.0000
      race_white    Binary   0.0000           .                  0.0000
      married       Binary   0.0000           .                  0.0000
      nodegree      Binary  -0.1290           .                  0.1290
      re74         Contin.  -0.0152      1.3924     Balanced, <2 0.2043
      re75         Contin.  -0.2407      2.6542 Not Balanced, >2 0.7419
                        KS.Threshold
      distance                      
      age         Not Balanced, >0.1
      educ        Not Balanced, >0.1
      race_black      Balanced, <0.1
      race_hispan     Balanced, <0.1
      race_white      Balanced, <0.1
      married         Balanced, <0.1
      nodegree    Not Balanced, >0.1
      re74        Not Balanced, >0.1
      re75        Not Balanced, >0.1
      

# MatchIt: full matching

    Code
      bal.tab(m4, int = T, ks.threshold = 0.05)
    Output
      Call
       MatchIt::matchit(formula = f.build("treat", covs), data = lalonde, 
          method = "full", distance = "glm", link = "probit", estimand = "ATE")
      
      Balance Measures
                                   Type Diff.Adj KS.Adj        KS.Threshold
      distance                 Distance   0.0043 0.0987                    
      age                       Contin.  -0.1864 0.2046 Not Balanced, >0.05
      educ                      Contin.   0.1588 0.1016 Not Balanced, >0.05
      race_black                 Binary   0.0029 0.0029     Balanced, <0.05
      race_hispan                Binary  -0.0015 0.0015     Balanced, <0.05
      race_white                 Binary  -0.0014 0.0014     Balanced, <0.05
      married                    Binary  -0.0045 0.0045     Balanced, <0.05
      nodegree                   Binary  -0.1016 0.1016 Not Balanced, >0.05
      re74                      Contin.  -0.2312 0.2455 Not Balanced, >0.05
      re75                      Contin.  -0.2711 0.2308 Not Balanced, >0.05
      age * educ                Contin.  -0.0290 0.1442 Not Balanced, >0.05
      age * race_black          Contin.   0.0288 0.1064 Not Balanced, >0.05
      age * race_hispan         Contin.  -0.0658 0.0439     Balanced, <0.05
      age * race_white          Contin.  -0.1170 0.1314 Not Balanced, >0.05
      age * married_0           Contin.   0.0260 0.1961 Not Balanced, >0.05
      age * married_1           Contin.  -0.1354 0.1393 Not Balanced, >0.05
      age * nodegree_0          Contin.   0.1473 0.1671 Not Balanced, >0.05
      age * nodegree_1          Contin.  -0.2570 0.1238 Not Balanced, >0.05
      age * re74                Contin.  -0.3229 0.2455 Not Balanced, >0.05
      age * re75                Contin.  -0.2871 0.2236 Not Balanced, >0.05
      educ * race_black         Contin.   0.0014 0.0216     Balanced, <0.05
      educ * race_hispan        Contin.   0.0155 0.0247     Balanced, <0.05
      educ * race_white         Contin.   0.0755 0.1218 Not Balanced, >0.05
      educ * married_0          Contin.  -0.0054 0.0517 Not Balanced, >0.05
      educ * married_1          Contin.   0.0867 0.1245 Not Balanced, >0.05
      educ * nodegree_0         Contin.   0.1901 0.1016 Not Balanced, >0.05
      educ * nodegree_1         Contin.  -0.1651 0.0974 Not Balanced, >0.05
      educ * re74               Contin.  -0.1654 0.2405 Not Balanced, >0.05
      educ * re75               Contin.  -0.2567 0.2482 Not Balanced, >0.05
      race_black * married_0     Binary  -0.0066 0.0066     Balanced, <0.05
      race_black * married_1     Binary   0.0096 0.0096     Balanced, <0.05
      race_black * nodegree_0    Binary  -0.0146 0.0146     Balanced, <0.05
      race_black * nodegree_1    Binary   0.0176 0.0176     Balanced, <0.05
      race_black * re74         Contin.   0.0328 0.0770 Not Balanced, >0.05
      race_black * re75         Contin.   0.0007 0.0683 Not Balanced, >0.05
      race_hispan * married_0    Binary   0.0150 0.0150     Balanced, <0.05
      race_hispan * married_1    Binary  -0.0165 0.0165     Balanced, <0.05
      race_hispan * nodegree_0   Binary  -0.0056 0.0056     Balanced, <0.05
      race_hispan * nodegree_1   Binary   0.0040 0.0040     Balanced, <0.05
      race_hispan * re74        Contin.  -0.1263 0.0485     Balanced, <0.05
      race_hispan * re75        Contin.  -0.1024 0.0434     Balanced, <0.05
      race_white * married_0     Binary  -0.0038 0.0038     Balanced, <0.05
      race_white * married_1     Binary   0.0024 0.0024     Balanced, <0.05
      race_white * nodegree_0    Binary   0.1218 0.1218 Not Balanced, >0.05
      race_white * nodegree_1    Binary  -0.1231 0.1231 Not Balanced, >0.05
      race_white * re74         Contin.  -0.2602 0.1387 Not Balanced, >0.05
      race_white * re75         Contin.  -0.3474 0.1716 Not Balanced, >0.05
      married_0 * nodegree_0     Binary   0.0190 0.0190     Balanced, <0.05
      married_0 * nodegree_1     Binary  -0.0144 0.0144     Balanced, <0.05
      married_0 * re74          Contin.  -0.0453 0.1208 Not Balanced, >0.05
      married_0 * re75          Contin.  -0.0836 0.1320 Not Balanced, >0.05
      married_1 * nodegree_0     Binary   0.0826 0.0826 Not Balanced, >0.05
      married_1 * nodegree_1     Binary  -0.0871 0.0871 Not Balanced, >0.05
      married_1 * re74          Contin.  -0.2223 0.1322 Not Balanced, >0.05
      married_1 * re75          Contin.  -0.2402 0.1683 Not Balanced, >0.05
      nodegree_0 * re74         Contin.   0.0055 0.0963 Not Balanced, >0.05
      nodegree_0 * re75         Contin.  -0.2082 0.1340 Not Balanced, >0.05
      nodegree_1 * re74         Contin.  -0.3151 0.2519 Not Balanced, >0.05
      nodegree_1 * re75         Contin.  -0.1537 0.1878 Not Balanced, >0.05
      re74 * re75               Contin.  -0.1808 0.1970 Not Balanced, >0.05
      
      Balance tally for KS statistics
                          count
      Balanced, <0.05        21
      Not Balanced, >0.05    38
      
      Variable with the greatest KS statistic
                Variable KS.Adj        KS.Threshold
       nodegree_1 * re74 0.2519 Not Balanced, >0.05
      
      Sample sizes
                           Control Treated
      All                    429.   185.  
      Matched (ESS)          257.3   24.68
      Matched (Unweighted)   429.   185.  

---

    Code
      bal.tab(m4, pairwise = FALSE, un = T)
    Output
      Call
       MatchIt::matchit(formula = f.build("treat", covs), data = lalonde, 
          method = "full", distance = "glm", link = "probit", estimand = "ATE")
      
      Balance summary across all treatment pairs
                      Type Max.Diff.Un Max.Diff.Adj
      distance    Distance      1.2347       0.0033
      age          Contin.      0.1690       0.2023
      educ         Contin.      0.0313       0.1712
      race_black    Binary      0.4475       0.0016
      race_hispan   Binary      0.0578       0.0049
      race_white    Binary      0.3897       0.0046
      married       Binary      0.2261       0.0126
      nodegree      Binary      0.0778       0.1110
      re74         Contin.      0.4162       0.2329
      re75         Contin.      0.2005       0.2630
      
      Sample sizes
                           Control Treated
      All                    429.   185.  
      Matched (ESS)          257.3   24.68
      Matched (Unweighted)   429.   185.  

# MatchIt: genetic matching, using matching weights

    Code
      bal.tab(m5, method = "m", estimand = "ATT")
    Output
      Call
       MatchIt::matchit(formula = f.build("treat", covs), data = lalonde, 
          method = "genetic", replace = F, ratio = 1, pop.size = 50)
      
      Balance Measures
                      Type Diff.Adj
      distance    Distance   1.0083
      age          Contin.  -0.0272
      educ         Contin.   0.0941
      race_black    Binary   0.3730
      race_hispan   Binary  -0.2216
      race_white    Binary  -0.1514
      married       Binary  -0.1081
      nodegree      Binary   0.0649
      re74         Contin.  -0.1967
      re75         Contin.  -0.0866
      
      Sample sizes
                Control Treated
      All           429     185
      Matched       185     185
      Unmatched     244       0

