#sourcing
for (i in dir("R/")) source(paste0("R/", i))
library(ggplot2)
library(crayon)
stop("Done sourcing.", call. = FALSE)

#Tests things quickly
#library("cobalt")
data("lalonde", package = "cobalt")
covs <- subset(lalonde, select = -c(re78, treat))
lalonde_ <- splitfactor(lalonde, "race", drop.first = F)
covs_ <- subset(lalonde_, select = -c(re78, treat))
#No Adjustment
bal.tab(covs, treat = lalonde$treat, m.threshold = .1, v.threshold = 2, imbalanced.only = T)
bal.tab(f.build("treat",covs), data = lalonde, ks.threshold = .1)

#MatchIt: matching w/ PS
library("MatchIt")
m1 <- matchit(f.build("treat", covs), data = lalonde, replace = T, ratio = 2, 
              discard = "both", reestimate = TRUE)
bal.tab(m1, int = T, v.threshold = 2, imbalanced.only = T)
#MatchIt: matching w/ distance
m1.1 <- matchit(f.build("treat", covs), data = lalonde, distance = runif(nrow(lalonde)))
bal.tab(m1.1, data = lalonde)
#MatchIt: matching w/Mahalanobis
m2 <- matchit(f.build("treat", covs[,-3]), data = lalonde, distance = "mahalanobis")
bal.tab(m2, data = lalonde, int = T, v.threshold = 2, addl = "race")
#MatchIt: subclassification
m3 <- matchit(f.build("treat", covs), data = lalonde, method = "subclass")
bal.tab(m3, int = F, v.threshold = 2, disp.subclass = T, ks.threshold = .1)
#Matchit: full matching
m4 <- matchit(f.build("treat", covs), data = lalonde, method = "full", distance = "probit")
bal.tab(m4, int = T, quick = T, ks.threshold = .05)
#Matchit: genetic matching, using matching weights
m5 <- matchit(f.build("treat", covs), data = lalonde, method = "genetic", replace = F,
              ratio = 1, print.level = 0, pop.size = 50)
bal.tab(m5, method = "m", estimand = "ATT")
#twang
library("twang")
ps.out <- ps(f.build("treat", covs), data = lalonde, 
   stop.method = c("ks.max", "es.max"), 
   estimand = "ATT", verbose = FALSE, n.trees = 2000, bag.fraction = .6)
bal.tab(ps.out, thresholds = c(k=.05, m = .05), stats = "k")
sampw <- sample(c(1.25, 0.75), nrow(covs), TRUE, c(.5, .5))
ps.out.s <- ps(f.build("treat", covs), data = lalonde, 
             stop.method = c("ks.max"), 
             estimand = "ATT", verbose = FALSE,
             sampw = sampw, n.trees = 1000)
bal.tab(ps.out.s, un = T)
bal.tab(covs, lalonde$treat, weights = get.w(ps.out.s), s.weights = sampw,
        un = T, distance = ps.out.s$ps)
#CBPS: binary
library("CBPS")
cbps.out <- CBPS(f.build("treat", covs), data = lalonde, ATT = F)
bal.tab(cbps.out, disp.bal.tab = T)
cbps.out.e <- CBPS(f.build("treat", covs), data = lalonde, method = "exact", ATT = F)
bal.tab(cbps.out, weights = list("exact" = (cbps.out.e)),
        stats = c("m", "v", "k"))
#Matching
library("Matching")
p.score <- glm(f.build("treat", covs), 
               data = lalonde, family = "binomial")$fitted.values
Match.out <- Match(Tr = lalonde$treat, X = p.score, estimand = "ATE",
                   M = 1, replace = T, CommonSupport = F)
bal.tab(Match.out, formula = f.build("treat", covs), data = lalonde)
gen <- GenMatch(Tr = lalonde$treat, X = covs_,
                   M = 1, replace = F, pop.size = 10, 
                print.level = 0, ties = F)
Gen.Match.out <- Match(Tr=lalonde$treat, X=covs_, Weight.matrix=gen,
                       M = 2, replace = F, ties = F)
bal.tab(Gen.Match.out, formula = f.build("treat", covs), data = lalonde,
        addl = data.frame(age2 = covs$age^2))
gen2 <- GenMatch(Tr = lalonde$treat, X = covs_,
                M = 1, replace = F, pop.size = 1000, 
                print.level = 0, ties = F, BalanceMatrix = cbind(covs_, age2 = covs_$age^2))
Gen.Match.out2 <- Match(Tr=lalonde$treat, X=covs_, Weight.matrix=gen2,
                       M = 1, replace = F, ties = F)
bal.tab(Gen.Match.out2, formula = f.build("treat", covs), data = lalonde,
        addl = data.frame(age2 = covs$age^2))

Matchby.out <- Matchby(Tr = lalonde$treat, X = p.score, by = lalonde$married)
bal.tab(Matchby.out, formula = f.build("treat", covs), data = lalonde,
        addl = ~I(age^2) + married*race)

#optmatch
library("optmatch")
ls <- lalonde[sample(1:614, 614, F),]
p.score <- glm(treat ~ age + educ + race + 
                              married + nodegree + re74 + re75, 
                          data = ls, family = binomial)$fitted.values
pm <- pairmatch(treat ~ qlogis(p.score), data = ls)

## Using formula and data (treatment not needed)
bal.tab(pm,  ~ age + educ + race + 
            married + nodegree + re74 + re75, data = ls,
        distance = p.score)

#WeightIt
library(WeightIt)
W <- weightit(f.build("treat", covs), data = lalonde,
              method = "ps", estimand = "ATT")
bal.tab(W)
W.cont <- weightit(f.build("re75", covs[-7]), data = lalonde,
              method = "ebal")
bal.tab(W.cont, stats = c("c", "sp"), disp = "means", un = T)
W.mult <- weightit(f.build("race", covs[-3]), data = lalonde,
                  method = "ps", estimand = "ATE",
                  focal = "black")
bal.tab(W.mult, disp = c("m", "sds"), which.treat = .all)

#Data frame/formula: weighting
glm1 <- glm(treat ~ age + educ + race, data = lalonde, 
            family = "binomial")
lalonde$distance <- glm1$fitted.values
lalonde$iptw.weights <- ifelse(lalonde$treat==1, 
                               1/lalonde$distance, 
                               1/(1-lalonde$distance))
bal.tab(f.build("treat", covs), data = lalonde, 
        weights = "iptw.weights", method = "weighting", 
        addl = data.frame(age2 = covs$age^2),
        distance = "distance")
#Data frame/formula: subclassification
lalonde$subclass <- findInterval(lalonde$distance, 
                                 quantile(lalonde$distance, (0:6)/6), all.inside = T)

bal.tab(covs, treat = lalonde$treat, subclass = lalonde$subclass, 
        disp.subclass = TRUE, addl = ~I(age^2),
        m.threshold = .1, v.threshold = 2)

bal.tab(covs, treat = lalonde$re75, subclass = lalonde$subclass, 
        disp.subclass = TRUE, addl = ~I(age^2))

#Entropy balancing
library("ebal")
e.out <- ebalance(lalonde$treat, cbind(covs_[,-3], age_2 = covs_$age^2, re74_2 = covs_$re74^2/1000, 
                                       re75_2 = covs_$re75^2/1000, educ_2 = covs_$educ^2))
bal.tab(e.out, treat = lalonde$treat, covs = covs, disp.ks = T, disp.v.ratio = T)
e.out.trim <- ebalance.trim(e.out)
bal.tab(e.out.trim, treat = lalonde$treat, covs = covs, disp.ks = T, disp.v.ratio = T)
bal.tab(covs, lalonde$treat, weights = list(e = e.out,
                                            e.trim = e.out.trim),
        disp.ks = T, disp.v.ratio = T)
#Continuous treatment (CBPS)
cbps.out2 <- CBPS(f.build("re78", covs), data = lalonde, method = "exact")
bal.tab(cbps.out2, stats = c("c", "sp"))
cbps.out2.e <- CBPS(f.build("re78", covs), data = do.call(data.frame, lapply(lalonde, rank)), method = "exact")
bal.tab(cbps.out2, weights = list(CBPS.sp = cbps.out2.e),
        stats = c("c", "sp"))
#Clustering with MatchIt
lalonde$school <- sample(LETTERS[1:4], nrow(lalonde), replace = T)
m5 <- matchit(f.build("treat", covs), data = lalonde, exact = "school")
bal.tab(m5, cluster = lalonde$school, disp.v.ratio = T, cluster.summary = T)
bal.tab(m5, cluster = lalonde$school, disp.ks = T, abs = T, cluster.fun = c("max"))
#Clustering w/ continuous treatment
bal.tab(cbps.out2, cluster = lalonde$school, un = T, cluster.fun = "mean")
#Multiple imputation
data("lalonde_mis", package = "cobalt")
covs_mis <- subset(lalonde_mis, select = -c(re78, treat))
lalonde_mis_ <- splitfactor(lalonde_mis, "race")
covs_mis_ <- subset(lalonde_mis_, select = -c(re78, treat))

library("mice")
imp <- mice(lalonde_mis, m = 3)
imp.data <- complete(imp, "long", include = FALSE)
imp.data <- imp.data[with(imp.data, order(.imp, .id)),]

ps <- match.weight <- rep(0, nrow(imp.data))
for (i in unique(imp.data$.imp)) {
    in.imp <- imp.data$.imp == i
    ps[in.imp] <- glm(f.build("treat", covs_mis), data = imp.data[in.imp,], 
                      family = "binomial")$fitted.values
    m.out <- matchit(treat ~ age, data = imp.data[in.imp,], distance = ps[in.imp])
    match.weight[in.imp] <- m.out$weights
}
imp.data <- cbind(imp.data, ps = ps, match.weight = match.weight)

bal.tab(f.build("treat", covs_mis), data = imp.data, weights = "match.weight", 
        method = "matching", imp = ".imp", distance = "ps", which.imp = .all,
        disp.ks = T)

bal.tab(f.build("treat", covs_mis), data = imp, weights = match.weight, 
        method = "matching", distance = ps, which.imp = NULL,
        cluster = "race", disp.v = T, abs = T, which.cluster = "black")

#With WeightIt
library(WeightIt)
W.imp <- weightit(f.build("treat", covs_mis), data = imp.data, estimand = "ATT",
                  by = ".imp", method = "optweight", moments = 2)
bal.tab(W.imp, imp = ".imp", disp.v.ratio = TRUE, abs = F)

#With MatchThem
library(MatchThem)
mt.out <- matchthem(f.build("treat", covs_mis), datasets = imp,
                    approach = "within")
bal.tab(mt.out)

wt.out <- weightthem(f.build("race", covs_mis[-3]), datasets = imp,
                    approach = "within", method = "optweight")
bal.tab(wt.out)
    
#With continuous treatment
imp.data <- complete(imp, "long", include = FALSE)
imp.data <- imp.data[with(imp.data, order(.imp, .id)),]
w <- rep(0, nrow(imp.data))
for (i in unique(imp.data$.imp)) {
    in.imp <- imp.data$.imp == i
    w[in.imp] <- get.w(CBPS(f.build("re78", covs_mis), data = imp.data[in.imp,],
                            method = "exact"))
}
imp.data <- cbind(imp.data, cbps.w = w)
bal.tab(f.build("re78", covs_mis), data = imp, weights = w, 
        which.imp = .all, un = T, thresholds = .1)

bal.tab(f.build("re78", covs_mis), data = imp.data, weights = "cbps.w", 
        method = "w", imp = ".imp", which.imp = NULL, cluster = "race")

#Missingness indicators
data("lalonde_mis", package = "cobalt")
covs_mis <- subset(lalonde_mis, select = -c(re78, treat))
lalonde_mis_ <- splitfactor(lalonde_mis, "race")
covs_mis_ <- subset(lalonde_mis_, select = -c(re78, treat))
library(twang)

ps.out <- ps(f.build("treat", covs_mis), data = lalonde_mis, 
             stop.method = c("es.max"), 
             estimand = "ATE", verbose = FALSE, n.trees = 1000)
bal.tab(ps.out, int = T)



#love.plot
v <- data.frame(old = c("age", "educ", "race_black", "race_hispan", 
                        "race_white", "married", "nodegree", "re74", "re75", "distance"),
                new = c("Age", "Years of Education", "Black", 
                        "Hispanic", "White", "Married", "No Degree Earned", 
                        "Earnings 1974", "Earnings 1975", "Propensity Score"))
plot(bal.tab(m1), threshold = .1, var.names = v)
love.plot(m1, stat = c("m", "v"), stars = "std", drop.distance = T)

love.plot(bal.tab(m1, cluster = lalonde$school), stat = "ks", agg.fun = "range", which.cluster = NA)
love.plot(cbps.out2.e, drop.distance = F, line = T, abs = T,
          var.order = "u")
love.plot(bal.tab(f.build("treat", covs_mis), data = imp.data, weights = "match.weight", 
                  method = "matching", imp = ".imp", distance = "ps"),
          agg.fun = "range")
love.plot(f.build("treat", covs_mis), data = imp.data, weights = "match.weight", 
                  method = "matching", imp = ".imp", distance = "ps",
                  cluster = "race", which.cluster = 1:3,
          agg.fun = "range", stat = c("m", "ks"))
#bal.plot
bal.plot(m1, "age", which = "both")
bal.plot(m1, "race")
bal.plot(cbps.out, "age", mirror = TRUE, type = "h", bins = 20, colors = c("white", "black"))
bal.plot(cbps.out, "race", which = "both")
bal.plot(cbps.out2, "age", which = "both")
bal.plot(cbps.out2, "race", which = "u")
bal.plot(m3, "age", which.sub = 2)
bal.plot(m3, "race", which.sub = 1:3, which = "both")
bal.plot(m5, "age", cluster = lalonde$school, which.cluster = c("A", "B"), which = "both")
bal.plot(m5, "race", cluster = lalonde$school, which.cluster = .none)
bal.plot(cbps.out2, "age", cluster = lalonde$school, which.cluster = c("A", "B"), which = "both")
bal.plot(cbps.out2, "race", cluster = lalonde$school)
bal.plot(f.build("treat", covs_mis), data = imp.data, weights = "match.weight", 
         method = "matching", imp = ".imp", var.name = "age", which = "b")
bal.plot(f.build("treat", covs_mis), data = imp.data, weights = "match.weight", 
         method = "matching", imp = ".imp", cluster = "race", which.imp = 1,
         var.name = "age", which = "b")

#Other packages

#ATE
library("ATE")
ate.att <- ATE(Y = lalonde_$re78, lalonde_$treat, covs_, ATT = T)
ate.att$weights.q[lalonde$treat == 1] <- 1
bal.tab(covs, lalonde$treat, weights = ate.att$weights.q, method = "w", estimand = "att", disp.v.ratio = T)

ate.ate <- ATE(Y = rep(0, nrow(lalonde)), lalonde_$treat, covs_, ATT = F, theta = 1)
ate.ate$weights <- ate.ate$weights.q + ate.ate$weights.p
bal.tab(covs, lalonde$treat, weights = ate.ate$weights, method = "w", estimand = "ate", disp.v.ratio = T)

#CMatching
library("CMatching")
lalonde$school <- sample(1:2, nrow(lalonde), replace = T)
ls <- lalonde[order(lalonde$school),]
p.score <- glm(treat ~ age + educ + race + married + nodegree + re74 + re75 + factor(school) - 1, 
               data = ls, family = "binomial")$fitted.values
MW <- MatchPW(Tr = ls$treat, X = p.score,
                   M = 1, replace = T, Group = ls$school,
             caliper = 2, estimand = "ATT")
bal.tab(MW, formula = f.build("treat", covs), data = ls, cluster = "school")

#designmatch
library(designmatch)
ls <- lalonde[order(lalonde$treat, decreasing = T),]
cs <- covs[order(lalonde$treat, decreasing = T),]
dmout <- bmatch(ls$treat,
                dist_mat = NULL,
                subset_weight = 1,
                total_groups = 100,
                mom = list(covs = as.matrix(cs[-(3:5)]),
                           tols = absstddif(as.matrix(cs[-(3:5)]), ls$treat, .0001)),
                # ks = list(covs = as.matrix(covs_[-c(3:7)]), 
                #           n_grid = 7,
                #            tols = rep(.05, 4)),
                n_controls = 1,
                solver = list(name = "glpk", approximate = 1)
                # , total_groups = 185
                , fine = list(covs = cs[c(3,4,5)])
                )

bal.tab(dmout, covs = cs, treat = ls$treat)
bal.tab(dmout, formula = treat ~ covs, data = ls)

#optweight (default method)
library(optweight)
W <- optweight(f.build("treat", covs), data = lalonde,
              estimand = "ATT")
bal.tab(W)
W.cont <- optweight(f.build("re75", covs[-7]), data = lalonde)
bal.tab(W.cont)
W.mult <- optweight(f.build("race", covs[-3]), data = lalonde,
                   estimand = "ATE",
                   focal = "black")
bal.tab(W.mult)

#Multinomial
lalonde$treat3 <- factor(ifelse(lalonde$treat == 1, "A", sample(c("B", "C"), nrow(lalonde), T)))
bal.tab(f.build("treat3", covs), data = lalonde, focal = 1, which.treat = 1:3, m.threshold = .1)
mnps3.out <- mnps(f.build("treat3", covs), data = lalonde, 
                  stop.method = c("es.mean", "es.max"), 
                  estimand = "ATE", verbose = FALSE,
                  n.trees = 200)
bal.tab(mnps3.out, which.treat = 1:3)
bal.plot(mnps3.out, var.name = "age")
mnps3.att <- mnps(f.build("treat3", covs), data = lalonde, 
                  stop.method = c("ks.max"), 
                  estimand = "ATT", verbose = FALSE,
                  treatATT = "B")
bal.tab(mnps3.att, which.treat = .all)
bal.plot(mnps3.att, var.name = "age")
cbps.out3 <- CBPS(f.build("treat3", covs), data = lalonde)
bal.tab(cbps.out3)
bal.plot(cbps.out3, var.name = "age")
bal.plot(f.build("treat3", covs), data = lalonde, var.name = "age",
         weights = list((cbps.out3),
                              (mnps3.out)),
         method = "w", which = "both")
bal.tab(f.build("treat3", covs), data = lalonde, 
         weights = data.frame(cbps = get.w(cbps.out3),
                              gbm = get.w(mnps3.out)),
         method = "w")
ate3.out <- ATE(rep(0, nrow(lalonde)), as.numeric(lalonde$treat3)-1,
                covs_, ATT = TRUE)
ate3.out$weights <- apply(ate3.out$weights.mat, 2, sum)
bal.plot(f.build("treat3", covs), data = lalonde, var.name = "age",
         weights = data.frame(ate = ate3.out$weights),
         method = "w", which = "both")
bal.tab(f.build("treat3", covs), data = lalonde,
         weights = data.frame(ate = ate3.out$weights),
         method = "w")
#MSMs
library("twang")
data(iptwExWide, package = "twang")
bal.tab(list(iptwExWide[c("use0", "gender", "age")], iptwExWide[c("use0", "gender", "age", "use1", "tx1")]), treat.list = iptwExWide[c("tx1", "tx2")])
bal.tab(list(tx1 ~ use0 + gender + age,
             tx2 ~ use1 + use0 + tx1 + gender + age,
             tx3 ~ use2 + use1 + use0 + tx2 + tx1 + gender + age),
        data = iptwExWide)

iptw.Ex <- iptw(list(tx1 ~ use0 + gender + age,
                     tx2 ~ use1 + use0 + tx1 + gender + age,
                     tx3 ~ use2 + use1 + use0 + tx2 + tx1 + gender + age),
                timeInvariant ~ gender + age,
                data = iptwExWide,
                cumulative = FALSE,
                priorTreatment = FALSE,
                verbose = FALSE,
                stop.method = "es.max",
                n.trees = 2000)
bal.tab(iptw.Ex)
data("iptwExLong")
iptw.l <- iptw(tx ~ gender + age + use, data = iptwExLong$covariates, 
               timeIndicators = iptwExLong$covariates$time, ID = iptwExLong$covariates$ID,
               n.trees = 200, stop.method = "es.max",
               verbose = FALSE)
bal.tab(iptw.l)

library("WeightIt")
Wmsm <- weightitMSM(list(tx1 ~ use0 + gender + age,
                     tx2 ~ use1 + use0 + tx1 + gender + age,
                     tx3 ~ use2 + use1 + use0 + tx2 + tx1 + gender + age),
                data = iptwExWide,
                method = "ps")
bal.tab(Wmsm)

library("CBPS")
data(iptwExLong, package = "twang")
cbps.msm <- CBMSM(tx ~ age + use,
                data = iptwExLong$covariates,
                id = iptwExLong$covariates$ID,
                time = iptwExLong$covariates$time,
                time.vary = T, msm.variance = "full")
bal.tab(cbps.msm)
W.msm <- weightitMSM(list(tx1 ~ use0 + gender + age,
                         tx2 ~ use1 + use0 + tx1 + gender + age,
                         tx3 ~ use2 + use1 + use0 + tx2 + tx1 + gender + age),
                    data = iptwExWide, stabilize = F,
                    verbose = FALSE,
                    method = "ps")
bal.tab(W.msm)

#Target checking
library(WeightIt)
w1 <- weightit(lalonde$treat ~ covs, estimand = "ATE")
target.bal.tab(w1, which.treat = NULL)

w2 <- weightit(lalonde$race ~ covs[-3], estimand = "ATE", method = "cbps")
target.bal.tab(w2, which.treat = NULL)

w3 <- weightit(lalonde$re78 ~ covs, estimand = "ATE", method = "cbps", over = FALSE)
target.bal.tab(w3, which.treat = NULL)
target.bal.tab(w3, disp.means = T, int = 2, imbalanced.only = T, v.threshold = .5)

w4 <- weightitMSM(list(tx1 ~ use0 + gender + age,
                          tx2 ~ use1 + use0 + tx1 + gender + age,
                          tx3 ~ use2 + use1 + use0 + tx2 + tx1 + gender + age),
                     data = iptwExWide,
                     verbose = FALSE,
                     method = "ps")
target.bal.tab(w4, which.treat = NULL)


a <- function() {
    d <- lalonde
lapply(1:3, function(i) {
    lapply(1, function(x) {
        b <- do.call(bal.tab, list(f.build("treat",covs), data = d, cluster = "race", which.cluster = NULL), quote = T)
    })
})
}
