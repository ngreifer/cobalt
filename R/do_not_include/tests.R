#Tests things quickly
library("cobalt")
data("lalonde", package = "cobalt")
covs <- subset(lalonde, select = -c(re78, treat))
#MatchIt: matching w/ PS
library("MatchIt")
m1 <- matchit(f.build("treat", covs), data = lalonde, ratio = 2)
bal.tab(m1, int = T, quick = T, disp.v.ratio = T)
#MatchIt: matching w/Mahalanobis
m2 <- matchit(f.build("treat", covs), data = lalonde, distance = "mahalanobis")
bal.tab(m2, int = T, quick = T, disp.v.ratio = T)
#MatchIt: subclassification
m3 <- matchit(f.build("treat", covs), data = lalonde, method = "subclass")
bal.tab(m3, int = T, quick = T, disp.v.ratio = T, disp.subclass = T)
#Matchit: full matching
m4 <- matchit(f.build("treat", covs), data = lalonde, method = "full")
bal.tab(m4, int = T, quick = T, disp.v.ratio = T)
#twang
library("twang")
ps.out <- ps(f.build("treat", covs), data = lalonde, 
   stop.method = c("ks.max"), 
   estimand = "ATT", verbose = FALSE)
bal.tab(ps.out)
#CBPS: binary
library("CBPS")
cbps.out <- CBPS(treat ~ age + educ + black + hispan + married + nodegree + re74 + 
                     re75, data = lalonde, 
             standardize = FALSE)
bal.tab(cbps.out)
#Matching
library("Matching")
Match.out <- Match(Tr = lalonde$treat, X = glm(treat ~ age + educ + black + hispan + 
                                                   married + nodegree + re74 + re75, 
                                               data = lalonde, family = binomial)$fitted.values)
bal.tab(Match.out, treat ~ age + educ + black + hispan + 
            married + nodegree + re74 + re75, data = lalonde)
#Data frame/formula: weighting
glm1 <- glm(treat ~ age + educ + black + hispan, data = lalonde, 
            family = "binomial")
lalonde$distance <- glm1$fitted.values
lalonde$iptw.weights <- ifelse(lalonde$treat==1, 
                               1/lalonde$distance, 
                               1/(1-lalonde$distance))
bal.tab(covs, treat = "treat", data = lalonde, 
        weights = "iptw.weights", method = "weighting", 
        s.d.denom = "pooled")
bal.tab(f.build("treat", covs), data = lalonde, 
        weights = "iptw.weights", method = "weighting", 
        s.d.denom = "pooled")
#Data frame/formula: subclassification
lalonde$subclass <- findInterval(lalonde$distance, 
                                 quantile(lalonde$distance[lalonde$treat==1], (0:6)/6), all.inside = T)

bal.tab(covs, treat = lalonde$treat, subclass = lalonde$subclass, 
        method = "subclassification", disp.subclass = TRUE)
#Entropy balancing
library("ebal")
e.out <- ebalance(lalonde$treat, covs)
bal.tab(e.out, treat = lalonde$treat, covs = covs)
#Continuous treatment (CBPS)
cbps.out2 <- CBPS(re78 ~ age + educ + black + hispan + married + nodegree + re74 + 
                     re75, data = lalonde, 
                 standardize = FALSE)
bal.tab(cbps.out2)
#Clustering with MatchIt
lalonde$school <- sample(LETTERS[1:4], nrow(lalonde), replace = T)
m5 <- matchit(f.build("treat", covs), data = lalonde, exact = "school")
bal.tab(m5, quick = T, cluster = lalonde$school)
#Clustering w/ continuous treatment
bal.tab(cbps.out2, quick = T, cluster = lalonde$school)
#love.plot
love.plot(bal.tab(m1))
love.plot(bal.tab(m5, cluster = lalonde$school), cluster.fun = "range")
love.plot(bal.tab(cbps.out2))
#bal.plot
bal.plot(m1, "age")
bal.plot(m1, "black")
bal.plot(cbps.out, "age")
bal.plot(cbps.out, "black")
bal.plot(cbps.out2, "age")
bal.plot(cbps.out2, "black")
bal.plot(m3, "age", which.sub = 2)
bal.plot(m3, "black", which.sub = 2)
bal.plot(m5, "age", cluster = lalonde$school, which.cluster = "A")
bal.plot(m5, "black", cluster = lalonde$school, which.cluster = 2)
bal.plot(cbps.out2, "age", cluster = lalonde$school, which.cluster = "A")
bal.plot(cbps.out2, "black", cluster = lalonde$school, which.cluster = 2)