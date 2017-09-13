#Tests things quickly
library("cobalt")
data("lalonde", package = "cobalt")
covs <- subset(lalonde, select = -c(re78, treat))
lalonde_ <- splitfactor(lalonde, "race")
covs_ <- subset(lalonde_, select = -c(re78, treat))
#No Adjustment
bal.tab(covs, treat = lalonde$treat, m.threshold = .1, v.threshold = 2)
bal.tab(f.build("treat",covs), lalonde, ks.threshold = .1)

#MatchIt: matching w/ PS
library("MatchIt")
m1 <- matchit(f.build("treat", covs), data = lalonde, replace = T, ratio = 2, discard = "control")
bal.tab(m1, int = T, quick = T, v.threshold = 2)
#MatchIt: matching w/Mahalanobis
m2 <- matchit(f.build("treat", covs_), data = lalonde_, distance = "mahalanobis", 
              discard = "hull.control")
bal.tab(m2, int = T, quick = T, v.threshold = 2)
#MatchIt: subclassification
m3 <- matchit(f.build("treat", covs), data = lalonde, method = "subclass")
bal.tab(m3, int = T, quick = T, v.threshold = 2, disp.subclass = T, ks.threshold = .1)
#Matchit: full matching
m4 <- matchit(f.build("treat", covs), data = lalonde, method = "full", distance = "probit")
bal.tab(m4, int = T, quick = T, ks.threshold = .05)
#Matchit: genetic matching, using matching weights
m5 <- matchit(f.build("treat", covs), data = lalonde, method = "genetic", replace = F,
              ratio = 2, print.level = 0, pop.size = 1000)
bal.tab(f.build("treat", covs), data = lalonde, weights = m5$weights, method = "m",
        estimand = "ATT")
#twang
library("twang")
ps.out <- ps(f.build("treat", covs), data = lalonde, 
   stop.method = c("ks.max", "es.max"), 
   estimand = "ATT", verbose = FALSE)
bal.tab(ps.out, disp.ks = T)
sampw <- sample(c(1.25, 0.75), nrow(covs), TRUE, c(.5, .5))
ps.out.s <- ps(f.build("treat", covs), data = lalonde, 
             stop.method = c("ks.max"), 
             estimand = "ATT", verbose = FALSE,
             sampw = sampw)
bal.tab(ps.out.s, un = T)
bal.tab(covs, lalonde$treat, weights = get.w(ps.out.s), s.weights = sampw,
        un = T, distance = ps.out.s$ps)
#CBPS: binary
library("CBPS")
cbps.out <- CBPS(f.build("treat", covs), data = lalonde)
bal.tab(cbps.out)
cbps.out.e <- CBPS(f.build("treat", covs), data = lalonde, method = "exact")
bal.tab(covs, lalonde$treat, weights = list("ps" = get.w(cbps.out),
                                            "exact" = get.w(cbps.out.e)),
        method = "w", estimand = "att", disp.ks = T, disp.v.ratio = T)
#Matching
library("Matching")
p.score <- glm(f.build("treat", covs), 
               data = lalonde, family = "binomial")$fitted.values
Match.out <- Match(Tr = lalonde$treat, X = p.score,
                   M = 1, replace = F, CommonSupport = T)
bal.tab(Match.out, formula = f.build("treat", covs), data = lalonde)
gen <- GenMatch(Tr = lalonde$treat, X = covs_,
                   M = 2, replace = F, pop.size = 1000, 
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
                                 quantile(lalonde$distance[lalonde$treat==1], (0:8)/8), all.inside = T)

bal.tab(covs, treat = lalonde$treat, subclass = lalonde$subclass, 
        method = "subclassification", disp.subclass = TRUE, addl = covs$age^2,
        m.threshold = .1, v.threshold = 2)
#Entropy balancing
library("ebal")
e.out <- ebalance(lalonde$treat, cbind(covs_, age_2 = covs_$age^2, re74_2 = covs_$re74^2/1000, 
                                       re75_2 = covs_$re75^2/1000, educ_2 = covs_$educ^2,
                                       covs_$re75^3/1000000, covs_$age^3.5/1000))
bal.tab(e.out, treat = lalonde$treat, covs = covs, disp.ks = T, disp.v.ratio = T)
e.out.trim <- ebalance.trim(e.out)
class(e.out.trim) <- "ebalance"
bal.tab(e.out.trim, treat = lalonde$treat, covs = covs, disp.ks = T, disp.v.ratio = T)
bal.tab(covs, lalonde$treat, weights = list(e = get.w(e.out, treat = lalonde$treat),
                                            e.trim = get.w(e.out.trim, treat = lalonde$treat)),
        disp.ks = T, disp.v.ratio = T)
#Continuous treatment (CBPS)
cbps.out2 <- CBPS(f.build("re78", covs), data = lalonde)
bal.tab(cbps.out2)
cbps.out2.e <- CBPS(f.build("re78", covs), data = lalonde,method = "exact")
bal.tab(covs, lalonde$re78, weights = list(c = get.w(cbps.out2),
                                            c.e = get.w(cbps.out2.e)),
        disp.ks = T)
#Clustering with MatchIt
lalonde$school <- sample(LETTERS[1:4], nrow(lalonde), replace = T)
m5 <- matchit(f.build("treat", covs), data = lalonde, exact = "school")
bal.tab(m5, quick = T, cluster = lalonde$school, disp.v.ratio = T)
print(bal.tab(m5, quick = T, cluster = lalonde$school, disp.ks = T), cluster.fun = c("mean", "max"))
#Clustering w/ continuous treatment
bal.tab(cbps.out2.e, quick = T, cluster = lalonde$school, un = T)
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
for (i in levels(imp.data$.imp)) {
    in.imp <- imp.data$.imp == i
    ps[in.imp] <- glm(f.build("treat", covs_mis), data = imp.data[in.imp,], 
                      family = "binomial")$fitted.values
    m.out <- matchit(treat ~ age, data = imp.data[in.imp,], distance = ps[in.imp])
    match.weight[in.imp] <- m.out$weights
}
imp.data <- cbind(imp.data, ps = ps, match.weight = match.weight)

bal.tab(f.build("treat", covs_mis), data = imp.data, weights = "match.weight", 
        method = "matching", imp = ".imp", distance = "ps", which.imp = NULL,
        disp.ks = T)

bal.tab(f.build("treat", covs_mis), data = imp.data, weights = "match.weight", 
        method = "matching", imp = ".imp", distance = "ps", which.imp = NULL,
        cluster = "race", disp.ks = T)
#With continuous treatment
imp.data <- complete(imp, "long", include = FALSE)
imp.data <- imp.data[with(imp.data, order(.imp, .id)),]
w <- rep(0, nrow(imp.data))
for (i in levels(imp.data$.imp)) {
    in.imp <- imp.data$.imp == i
    w[in.imp] <- get.w(CBPS(f.build("re78", covs_mis), data = imp.data[in.imp,],
                            method = "exact"))
}
imp.data <- cbind(imp.data, cbps.w = w)
bal.tab(f.build("re78", covs_mis), data = imp.data, weights = "cbps.w", 
        method = "w", imp = ".imp", which.imp = NULL, un = T, r.threshold = .1)

bal.tab(f.build("re78", covs_mis), data = imp.data, weights = "cbps.w", 
        method = "w", imp = ".imp", which.imp = NULL, cluster = "race")

#love.plot
plot(bal.tab(m1), threshold = .1, limits = c(-1, 1.5))
love.plot(bal.tab(m5, cluster = lalonde$school), stat = "ks", agg.fun = "range")
love.plot(bal.tab(cbps.out2.e), drop.distance = F, line = T, abs = T,
          var.order = "u")
love.plot(bal.tab(f.build("treat", covs_mis), data = imp.data, weights = "match.weight", 
                  method = "matching", imp = ".imp", distance = "ps"),
          agg.fun = "range")
love.plot(bal.tab(f.build("treat", covs_mis), data = imp.data, weights = "match.weight", 
                  method = "matching", imp = ".imp", distance = "ps",
                  cluster = "race", which.cluster = 1:3),
          agg.fun = "range", stat = "ks")
#bal.plot
bal.plot(m1, "age", which = "both")
bal.plot(m1, "race")
bal.plot(cbps.out, "age")
bal.plot(cbps.out, "race", which = "both")
bal.plot(cbps.out2, "age", which = "both")
bal.plot(cbps.out2, "race", which = "u")
bal.plot(m3, "age", which.sub = 2)
bal.plot(m3, "race", which.sub = 1:3, which = "both")
bal.plot(m5, "age", cluster = lalonde$school, which.cluster = c("A", "B"), which = "both")
bal.plot(m5, "race", cluster = lalonde$school, which.cluster = NULL)
bal.plot(cbps.out2, "age", cluster = lalonde$school, which.cluster = c("A", "B"), which = "both")
bal.plot(cbps.out2, "race", cluster = lalonde$school, which.cluster = NA)
bal.plot(f.build("treat", covs_mis), data = imp.data, weights = "match.weight", 
         method = "matching", imp = ".imp", var.name = "age", which = "b")
bal.plot(f.build("treat", covs_mis), data = imp.data, weights = "match.weight", 
         method = "matching", imp = ".imp", cluster = "race", which.imp = 1,
         var.name = "age", which = "b")

#Add MI tests
#Other packages
#sbw
library("sbw")
s <- sbw(splitfactor(lalonde, drop.first = F), "treat", 
         names(splitfactor(covs, drop.first = F)), .001, target = "treated", 
         solver = "quadprog")
s$w <- s$data_frame_weights$weights; s$w[lalonde$treat==1] <- 1; s$w[s$w < 0] <- 0
bal.tab(covs, lalonde$treat, weights = s$w, method = "w", estimand = "att", disp.v.ratio = T)
#ATE
library("ATE")
ate.att <- ATE(Y = lalonde_$re78, lalonde_$treat, covs_, ATT = T)
ate.att$weights.q[lalonde$treat == 1] <- 1
bal.tab(covs, lalonde$treat, weights = ate.att$weights.q, method = "w", estimand = "att", disp.v.ratio = T)

ate.ate <- ATE(Y = rep(0, nrow(lalonde)), lalonde_$treat, covs_, ATT = F, theta = 1)
ate.ate$weights <- ate.ate$weights.q + ate.ate$weights.p
bal.tab(covs, lalonde$treat, weights = ate.ate$weights, method = "w", estimand = "ate", disp.v.ratio = T)

#Multinomial
lalonde$treat4 <- factor(ifelse(lalonde$treat == 1, 1, sample(2:4, nrow(lalonde), T)))
mnps4.out <- mnps(f.build("treat4", covs), data = lalonde, 
                  stop.method = c("ks.max"), 
                  estimand = "ATE", verbose = FALSE)
bal.tab(mnps4.out)
bal.plot(mnps4.out, var.name = "age")
mnps4.att <- mnps(f.build("treat4", covs), data = lalonde, 
                  stop.method = c("ks.max"), 
                  estimand = "ATT", verbose = FALSE,
                  treatATT = 4)
bal.tab(mnps4.att)
bal.plot(mnps4.att, var.name = "age")
cbps.out4 <- CBPS(f.build("treat4", covs), data = lalonde)
bal.tab(cbps.out4)
bal.plot(cbps.out4, var.name = "age")
bal.plot(f.build("treat4", covs), data = lalonde, var.name = "age",
         weights = data.frame(cbps = get.w(cbps.out4),
                              gbm = get.w(mnps4.out)),
         method = "w", which = "both")

#sourcing
source('R/x2base.R')
source('R/bal.tab.R')
source('R/functions_for_processing.R')
source('R/print.bal.tab.R')
source('R/utilities.R')
source('R/love.plot.R')
source('R/bal.plot.R')
library(ggplot2)
