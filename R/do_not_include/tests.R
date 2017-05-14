#Tests things quickly
library("cobalt")
data("lalonde", package = "cobalt")
covs <- subset(lalonde, select = -c(re78, treat))
lalonde_ <- split.factor("race", lalonde)
covs_ <- subset(lalonde_, select = -c(re78, treat))
#MatchIt: matching w/ PS
library("MatchIt")
m1 <- matchit(f.build("treat", covs), data = lalonde, ratio = 2)
bal.tab(m1, int = T, quick = T, v.threshold = 2)
#MatchIt: matching w/Mahalanobis
m2 <- matchit(f.build("treat", covs), data = lalonde, distance = "mahalanobis")
bal.tab(m2, int = T, quick = T, v.threshold = 2)
#MatchIt: subclassification
m3 <- matchit(f.build("treat", covs), data = lalonde, method = "subclass")
bal.tab(m3, int = T, quick = T, v.threshold = 2, disp.subclass = T)
#Matchit: full matching
m4 <- matchit(f.build("treat", covs), data = lalonde, method = "full")
bal.tab(m4, int = T, quick = T, v.threshold = 2)
#Matchit: genetic matching, using matching weights
m5 <- matchit(f.build("treat", covs), data = lalonde, method = "genetic", replace = T,
              ratio = 2, print.level = 0, pop.size = 1000)
bal.tab(f.build("treat", covs_), data = lalonde_, weights = m5$weights, method = "m",
        estimand = "ATT")
#twang
library("twang")
ps.out <- ps(f.build("treat", covs), data = lalonde, 
   stop.method = c("ks.max"), 
   estimand = "ATT", verbose = FALSE)
bal.tab(ps.out)
#CBPS: binary
library("CBPS")
cbps.out <- CBPS(f.build("treat", covs), data = lalonde)
bal.tab(cbps.out)
#Matching
library("Matching")
p.score <- glm(f.build("treat", covs), 
               data = lalonde, family = "binomial")$fitted.values
Match.out <- Match(Tr = lalonde$treat, X = p.score,
                   M = 2, replace = F)
bal.tab(Match.out, formula = f.build("treat", covs), data = lalonde)
gen <- GenMatch(Tr = lalonde$treat, X = covs_,
                   M = 2, replace = F, pop.size = 1000, 
                print.level = 0, ties = F)
Gen.Match.out <- Match(Tr=lalonde$treat, X=covs_, Weight.matrix=gen,
                       M = 2, replace = F, ties = F)
bal.tab(Gen.Match.out, formula = f.build("treat", covs), data = lalonde,
        addl = data.frame(age2 = covs$age^2))
#Data frame/formula: weighting
glm1 <- glm(treat ~ age + educ + race, data = lalonde, 
            family = "binomial")
lalonde$distance <- glm1$fitted.values
lalonde$iptw.weights <- ifelse(lalonde$treat==1, 
                               1/lalonde$distance, 
                               1/(1-lalonde$distance))
bal.tab(covs, treat = "treat", data = lalonde, 
        weights = "iptw.weights", method = "weighting")
bal.tab(f.build("treat", covs), data = lalonde, 
        weights = "iptw.weights", method = "weighting", 
        addl = data.frame(age2 = covs$age^2))
#Data frame/formula: subclassification
lalonde$subclass <- findInterval(lalonde$distance, 
                                 quantile(lalonde$distance[lalonde$treat==1], (0:6)/6), all.inside = T)

bal.tab(covs, treat = lalonde$treat, subclass = lalonde$subclass, 
        method = "subclassification", disp.subclass = TRUE, addl = covs$age^2,
        m.threshold = .1, v.threshold = 2)
#Entropy balancing
library("ebal")
e.out <- ebalance(lalonde$treat, covs_[,-3])
bal.tab(e.out, treat = lalonde$treat, covs = covs)
#Continuous treatment (CBPS)
cbps.out2 <- CBPS(f.build("re78", covs), data = lalonde)
bal.tab(cbps.out2)
#Clustering with MatchIt
lalonde$school <- sample(LETTERS[1:4], nrow(lalonde), replace = T)
m5 <- matchit(f.build("treat", covs), data = lalonde, exact = "school")
bal.tab(m5, quick = T, cluster = lalonde$school)
#Clustering w/ continuous treatment
bal.tab(cbps.out2, quick = T, cluster = lalonde$school)
#love.plot
love.plot(bal.tab(m1), threshold = .1)
love.plot(bal.tab(m5, cluster = lalonde$school), cluster.fun = "range", limits = c(0, 2))
love.plot(bal.tab(cbps.out2))
#bal.plot
bal.plot(m1, "age")
bal.plot(m1, "race")
bal.plot(cbps.out, "age")
bal.plot(cbps.out, "race")
bal.plot(cbps.out2, "age")
bal.plot(cbps.out2, "race")
bal.plot(m3, "age", which.sub = 2)
bal.plot(m3, "race", which.sub = 2)
bal.plot(m5, "age", cluster = lalonde$school, which.cluster = "A")
bal.plot(m5, "race", cluster = lalonde$school, which.cluster = 2)
bal.plot(cbps.out2, "age", cluster = lalonde$school, which.cluster = "A")
bal.plot(cbps.out2, "race", cluster = lalonde$school, which.cluster = 2)
#dd MI tests
#Other packages
#sbw
library("sbw")
s <- sbw(lalonde, "treat", names(cov), .001, "treated", solver = "quadprog")
s$w <- s$data_frame_weights$weights; s$w[lalonde$treat==1] <- 1
bal.tab(cov, lalonde$treat, weights = s$w, method = "w", estimand = "att", disp.v.ratio = T)
#ATE
library("ATE")
ate <- ATE(Y = lalonde$re78, lalonde$treat, cov, ATT = TRUE)
ate$weights.q[lalonde$treat == 1] <- 1
bal.tab(cov, lalonde$treat, weights = ate$weights.q, method = "w", estimand = "att", disp.v.ratio = T)

#sourcing
source('~/Dropbox (Personal)/Research/R/cobalt/R/x2base.R')
source('~/Dropbox (Personal)/Research/R/cobalt/R/bal.tab.R')
source('~/Dropbox (Personal)/Research/R/cobalt/R/functions_for_processing.R')
source('~/Dropbox (Personal)/Research/R/cobalt/R/print.bal.tab.R')
source('~/Dropbox (Personal)/Research/R/cobalt/R/utilities.R')