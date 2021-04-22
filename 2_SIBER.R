#Using SIBER to look at ellipses among "groups"
# I'll do an arbitrary M-F one first I think

library(tidyverse)
library(SIBER, quietly = TRUE,
        verbose = FALSE,
        logical.return = FALSE)

library(viridis)

palette(viridis(4))

ttt <- read.csv("data/TTT_fish_SI.csv", header = T, stringsAsFactors = F)
names(ttt)

fishmeta <- read.csv("data/TTT_fishdat.csv", header = T, stringsAsFactors = F)
names(fishmeta)

ttt <- ttt %>%
  left_join(fishmeta, by = "SampleID")

#### SIBER FOR M vs F ####
# prep data for SIBER
mydata <- ttt %>%
  select(iso1 = d13C, iso2 = d15N, group = Sex) %>%
  mutate(community = 1)

# create the siber object
siber.example <- createSiberObject(mydata)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, 
                             lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")


# ellipses and group.hulls are set to TRUE or T for short to force
# their plotting. 
par(mfrow=c(1,1))
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = FALSE, community.hulls.args, 
                ellipses = TRUE, group.ellipses.args,
                group.hulls = TRUE, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(siber.example, n = 100, p.interval = 0.50,
                  lty = 1, lwd = 2)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(siber.example, n = 100, p.interval = 0.95,
                  ci.mean = TRUE, lty = 1, lwd = 2)


# Calculate sumamry statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.example)
print(group.ML)

# add a legend
legend("topright", colnames(group.ML), 
       pch = c(21, 24), col = c(1:2), lty = 1)

## Run SIBER
# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for (adds n.iter to the list)
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors - default vague problems that you should probably never change
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method. with SIBER multivariate
ellipses.posterior <- siberMVN(siber.example, parms, priors)

# ----------------------------------------------------------------
# Plot out some of the data and results
# ----------------------------------------------------------------

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(SEA.B.credibles)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

print(SEA.B.modes)

## Compare posterior distributions
#this is to compare ellipse 1 to ellipse 2 of your dataset
Pg1.lt.g2 <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
print(Pg1.lt.g2)

# Overlap between ellipses
# the choice of ellipse size does matter here
# bigger ellipses increases chance for overlap
# area is in units of per mil squared
overlap.1M.1F <- maxLikOverlap("1.M", "1.F", siber.example, p = 0.95, n =100)
overlap.1M.1F

# Bayesian overlap
bayes.overlap.1M.1F <- bayesianOverlap("1.M", "1.F", ellipses.posterior, 
                                       draws = 10, p.interval = 0.95,
                                       n = 360)
print(bayes.overlap.1M.1F)

# get credible intervals
overlap.credibles <- lapply(
  as.data.frame(bayes.overlap.1M.1F), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(overlap.credibles)

#### SIBER FOR YEAR ####
# prep data for SIBER
mydata <- ttt %>%
  select(iso1 = d13C, iso2 = d15N, group = Year) %>%
  mutate(community = 1)

# create the siber object
siber.example <- createSiberObject(mydata)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, 
                             lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")


# ellipses and group.hulls are set to TRUE or T for short to force
# their plotting. 
par(mfrow=c(1,1))
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = FALSE, community.hulls.args, 
                ellipses = TRUE, group.ellipses.args,
                group.hulls = TRUE, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(siber.example, n = 100, p.interval = 0.50,
                  lty = 1, lwd = 2)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(siber.example, n = 100, p.interval = 0.95,
                  ci.mean = TRUE, lty = 1, lwd = 2)


# Calculate sumamry statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.example)
print(group.ML)

# add a legend
legend("topright", colnames(group.ML), 
       pch = c(21, 24, 25), col = c(1:3), lty = 1)

## Run SIBER
# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for (adds n.iter to the list)
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors - default vague problems that you should probably never change
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method. with SIBER multivariate
ellipses.posterior <- siberMVN(siber.example, parms, priors)

# ----------------------------------------------------------------
# Plot out some of the data and results
# ----------------------------------------------------------------

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(SEA.B.credibles)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

print(SEA.B.modes)

## Compare posterior distributions
#this is to compare ellipse 1 to ellipse 2 of your dataset
Pg1.lt.g3 <- sum( SEA.B[,1] < SEA.B[,3] ) / nrow(SEA.B)
print(Pg1.lt.g3)

# Overlap between ellipses
# the choice of ellipse size does matter here
# bigger ellipses increases chance for overlap
# area is in units of per mil squared
overlap.2014.2016 <- maxLikOverlap("1.2014", "1.2016", siber.example, p = 0.95, n =100)
overlap.2014.2016

# Bayesian overlap
bayes.overlap.2014.2016 <- bayesianOverlap("1.2014", "1.2016", ellipses.posterior, 
                                       draws = 10, p.interval = 0.95,
                                       n = 360)
print(bayes.overlap.2014.2016)

# get credible intervals
overlap.credibles <- lapply(
  as.data.frame(bayes.overlap.1M.1F), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(overlap.credibles)


#### SIBER FOR MONTH ####
table(ttt$Mo)
# 5  6  7  8  9 
# 24  8 51 12  1

#### SIBER FOR AREA ####
table(ttt$Area)# let's remove unknown and GOM for now and just see what happens

# prep data for SIBER
mydata <- ttt %>%
  filter(!Area %in% c("unknown", "GOM")) %>%
  mutate(Estuary = case_when(Area == "E_MS_Sound" ~ "MS_Sound",
                             Area == "W_MS_Sound" ~ "MS_Sound",
                             Area == "N_Mob_Bay" ~ "Mobile_Bay",
                             Area == "S_Mob_Bay" ~ "Mobile_Bay",
                             TRUE ~ Area)) %>%
  select(iso1 = d13C, iso2 = d15N, group = Estuary) %>%
  mutate(community = 1)

# create the siber object
siber.example <- createSiberObject(mydata)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, 
                             lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")


# ellipses and group.hulls are set to TRUE or T for short to force
# their plotting. 
par(mfrow=c(1,1))
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = FALSE, community.hulls.args, 
                ellipses = TRUE, group.ellipses.args,
                group.hulls = TRUE, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(siber.example, n = 100, p.interval = 0.50,
                  lty = 1, lwd = 2)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(siber.example, n = 100, p.interval = 0.95,
                  ci.mean = TRUE, lty = 1, lwd = 2)


# Calculate sumamry statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.example)
print(group.ML)

# add a legend
legend("topright", colnames(group.ML), 
       pch = c(1, 1), col = c(1:2), lty = 1)

## Run SIBER
# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for (adds n.iter to the list)
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors - default vague problems that you should probably never change
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method. with SIBER multivariate
ellipses.posterior <- siberMVN(siber.example, parms, priors)

# ----------------------------------------------------------------
# Plot out some of the data and results
# ----------------------------------------------------------------

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)
#THIS SEEMS REALLY STRANGE OUTPUT TO ME WITH TINY DENSITY PLOTS

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(SEA.B.credibles)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

print(SEA.B.modes)

## Compare posterior distributions
#this is to compare ellipse 1 to ellipse 2 of your dataset
Pg1.lt.g2 <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
print(Pg1.lt.g2)

# Overlap between ellipses
# the choice of ellipse size does matter here
# bigger ellipses increases chance for overlap
# area is in units of per mil squared
overlap.1M.1F <- maxLikOverlap("1.M", "1.F", siber.example, p = 0.95, n =100)
overlap.1M.1F

# Bayesian overlap
bayes.overlap.1M.1F <- bayesianOverlap("1.M", "1.F", ellipses.posterior, 
                                       draws = 10, p.interval = 0.95,
                                       n = 360)
print(bayes.overlap.1M.1F)

# get credible intervals
overlap.credibles <- lapply(
  as.data.frame(bayes.overlap.1M.1F), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(overlap.credibles)

#### SIBER FOR AREA C and S ####
table(ttt$Area)# let's remove unknown and GOM for now and just see what happens

# prep data for SIBER
mydata <- ttt %>%
  filter(!Area %in% c("unknown", "GOM")) %>%
  mutate(Estuary = case_when(Area == "E_MS_Sound" ~ "MS_Sound",
                             Area == "W_MS_Sound" ~ "MS_Sound",
                             Area == "N_Mob_Bay" ~ "Mobile_Bay",
                             Area == "S_Mob_Bay" ~ "Mobile_Bay",
                             TRUE ~ Area)) %>%
  select(iso1 = d13C, iso2 = d34S, group = Estuary) %>%
  mutate(community = 1)

# create the siber object
siber.example <- createSiberObject(mydata)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, 
                             lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")


# ellipses and group.hulls are set to TRUE or T for short to force
# their plotting. 
par(mfrow=c(1,1))
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = FALSE, community.hulls.args, 
                ellipses = TRUE, group.ellipses.args,
                group.hulls = TRUE, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^34*S~'\u2030')
)

# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(siber.example, n = 100, p.interval = 0.50,
                  lty = 1, lwd = 2)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(siber.example, n = 100, p.interval = 0.95,
                  ci.mean = TRUE, lty = 1, lwd = 2)


# Calculate sumamry statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.example)
print(group.ML)

# add a legend
legend("topright", colnames(group.ML), 
       pch = c(1, 1), col = c(1:2), lty = 1)

## Run SIBER
# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for (adds n.iter to the list)
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors - default vague problems that you should probably never change
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method. with SIBER multivariate
ellipses.posterior <- siberMVN(siber.example, parms, priors)

# ----------------------------------------------------------------
# Plot out some of the data and results
# ----------------------------------------------------------------

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)
#THIS SEEMS REALLY STRANGE OUTPUT TO ME WITH TINY DENSITY PLOTS

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(SEA.B.credibles)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

print(SEA.B.modes)

## Compare posterior distributions
#this is to compare ellipse 1 to ellipse 2 of your dataset
Pg1.lt.g2 <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
print(Pg1.lt.g2)

# Overlap between ellipses
# the choice of ellipse size does matter here
# bigger ellipses increases chance for overlap
# area is in units of per mil squared
overlap.1M.1F <- maxLikOverlap("1.M", "1.F", siber.example, p = 0.95, n =100)
overlap.1M.1F

# Bayesian overlap
bayes.overlap.1M.1F <- bayesianOverlap("1.M", "1.F", ellipses.posterior, 
                                       draws = 10, p.interval = 0.95,
                                       n = 360)
print(bayes.overlap.1M.1F)

# get credible intervals
overlap.credibles <- lapply(
  as.data.frame(bayes.overlap.1M.1F), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(overlap.credibles)