# April 20, 2021

# Start Looking at SIA data, after taking A. Parnell and A. Jackson's online
#  course for using simmr and MixSIAR

library(tidyverse)
library(simmr)

ttt <- read.csv("data/TTT_fish_SI.csv", header = T, stringsAsFactors = F)
names(ttt)

fishmeta <- read.csv("data/TTT_fishdat.csv", header = T, stringsAsFactors = F)
names(fishmeta)

ttt <- ttt %>%
  left_join(fishmeta, by = "SampleID")

sources <- read.csv("data/TTT_sources_SI.csv", header = T, stringsAsFactors = F)
names(sources)

all_spp <- ttt %>%
  full_join(sources)

p1 <- ggplot(data = all_spp,
             aes(x = d13C,
                 y = d15N)) +
  geom_point(aes(color = Species), size = 5)
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(text = element_text(size=20)) +
  theme_classic() +
  coord_equal()
p1

fish_summary <- all_spp %>%
  group_by(Species) %>%
  summarise(count = n(),
            mC = mean(d13C),
            sdC = sd(d13C),
            mN = mean(d15N),
            sdN = sd(d15N),
            mS = mean(d34S),
            sdS = sd(d34S))

p2 <- p1 +
  geom_errorbar(data = fish_summary,
                mapping = aes(x = mC, y = mN,
                              ymin = mN - 1.96*sdN,
                              ymax = mN + 1.96*sdN),
                        width = 0) +
  geom_errorbarh(data = fish_summary, 
                 mapping = aes(x = mC, y = mN,
                               xmin = mC - 1.96*sdC,
                               xmax = mC + 1.96*sdC),
                 height = 0) + 
  geom_point(data = fish_summary, aes(x = mC, 
                             y = mN,
                             fill = Species), 
             color = "black", shape = 22, size = 5,
             alpha = 0.7, show.legend = FALSE) 
p2                  


# Ellipses plot instead of error bars

p.ell <- 0.95

p_ellipse <- p1 +
  stat_ellipse(aes(group = Species, 
                   fill = Species, 
                   color = Species), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon")
p_ellipse

#### Let's Run simmr ####

# Need to prep the data:
Tripletail <- all_spp %>%
  filter(Species == "Atlantic Tripletail") %>%
  select(d13C, d15N, d34S) %>%
  as.matrix()
s_names <- all_spp %>%
  select(Species) %>%
  filter(Species != "Atlantic Tripletail") %>%
  distinct(Species) %>%
  pull()
s_means <- all_spp %>%
  filter(Species != "Atlantic Tripletail") %>%
  group_by(Species) %>%
  summarize(d13C = mean(d13C), d15N=mean(d15N), d34S = mean(d34S)) %>%
  #need to drop column of species names
  select(-Species) %>%
  as.matrix()
s_sds <- all_spp %>%
  filter(Species != "Atlantic Tripletail") %>%
  group_by(Species) %>%
  summarize(d13C = sd(d13C), d15N=sd(d15N), d34S = sd(d34S)) %>%
  #need to drop column of species names
  select(-Species) %>%
  as.matrix()

# let's try to estimate some TDFs or TEFs from our observed data
# library(tRophicPosition)
# TDFs_prep <- all_spp %>%
#   group_by(Species) %>%
#   summarise(nN = n(),
#             meanN = mean(d15N),
#             sdN = sd(d15N),
#             nC = n(),
#             meanC = mean(d13C),
#             sdC = sd(d13C),
#             seed = 3)
# 
# AB <- simulateTDF(nN = TDFs_prep[1,2],
#                   meanN = TDFs_prep[1,3],
#                   sdN = TDFs_prep[1,4])
# 
# AB <- simulateTDF(nN = 12,
#                   meanN = 14.8,
#                   sdN = 0.851)

# Load it in 
simmr_in = simmr_load(mixtures=Tripletail,
                      source_names=s_names,
                      source_means=s_means,
                      source_sds=s_sds)
simmr_in #gives you a quick summary; should give an error if you've mis-specified something

# Want help?
#help(simmr_load)

# Default iso-space
plot(simmr_in, tracers = c(1,2)) #plot C & N
plot(simmr_in, tracers = c(1,3)) #plot C & S
plot(simmr_in, tracers = c(2,3)) #plot N & S

# Better iso-space plots, for all pairs of isotopes
plot(simmr_in, tracers = c(1,2),
     xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")), 
     title='Isospace plot of Atlantic Tripletail and diet sources')

plot(simmr_in, tracers = c(1,3),
     xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     ylab=expression(paste(delta^34, "S (\u2030)",sep="")), 
     title='Isospace plot of Atlantic Tripletail and diet sources')

plot(simmr_in, tracers = c(2,3),
     xlab=expression(paste(delta^15, "N (\u2030)",sep="")), 
     ylab=expression(paste(delta^34, "S (\u2030)",sep="")), 
     title='Isospace plot of Atlantic Tripletail and diet sources')

# Getting help on this is a bit fiddlier
help(plot.simmr_input)

# Run simmr!
simmr_out = simmr_mcmc(simmr_in)
simmr_out
# Look at output ----------------------------------------------------------

# Always start with convergence
summary(simmr_out, type='diagnostics')
# If these are close to 1 (e.g. less than 1.1) you can continue; if not, try longer run of MCMC

# Now you can explore ...
help(summary.simmr_output)
summary(simmr_out, type = 'statistics')
summary(simmr_out, type = 'quantiles')
summary(simmr_out, type = 'correlations') # good if you have overlapping sources (those should be correlated)
#  typically negative because as one goes up, the other goes down; if positive then both props go up to keep consumers in middle of sources in plot
# Plot the output
help(plot.simmr_output)
plot(simmr_out, type = 'histogram')
plot(simmr_out, type = 'density')
plot(simmr_out, type = 'boxplot')
plot(simmr_out, type = 'matrix')

# To look at prior vs. posterior distributions
#  How to interpret:
#  You don't necessarily want the prior and posterior to line up or converge
#  If you’re using the default vague priors you want them to be different. 
#  If they look too similar then the data isn’t overwhelming the vague prior.
#  If you’ve got informative priors then you want them to be reasonably similar
#  because you want the prior to play a role in the estimation of the parameters.
#  A. Parnell says my plot with the vague default priors looks ok to him
prior_viz(simmr_out)

# plot the posterior predictive power 
posterior_predictive(simmr_out) # this fails saying differing number of rows 192 & 288 in matrix y_post_pred_ci???

# Do more with the output -------------------------------------------------

# Do more with the output -------------------------------------------------

# Compare two sources
compare_sources(simmr_out,source_names=c('Brown Shrimp','Northern White Shrimp'))

# Compare multiple sources
compare_sources(simmr_out,source_names=c('Brown Shrimp','Northern White Shrimp', 'Atlantic Croaker'))

# IF you want to combine sources a posteriori
# combine sources C and D which are in positions 3 and 4
simmr_out_a_posteriori <- combine_sources(
  simmr_out, 
  to_combine = simmr_out$input$source_names[c(3,4)], 
  new_source_name = "CD")

# Plot the a posteriori aggregated diet estimatess
plot(simmr_out_a_posteriori, type = "density")

# -----------------------------------------------------------------------
# Multiple groups ---------------------------------------------------------
# -----------------------------------------------------------------------
#### Let's try area groupings:

all_spp <- all_spp %>%
  mutate(Estuary = case_when(Area == "E_MS_Sound" ~ "MS_Sound",
                             Area == "W_MS_Sound" ~ "MS_Sound",
                             Area == "N_Mob_Bay" ~ "Mobile_Bay",
                             Area == "S_Mob_Bay" ~ "Mobile_Bay",
                             TRUE ~ Area))

# Need to prep the data:
Tripletail <- all_spp %>%
  filter(Species == "Atlantic Tripletail") %>%
  select(d13C, d15N, d34S) %>%
  as.matrix()
s_names <- all_spp %>%
  select(Species) %>%
  filter(Species != "Atlantic Tripletail") %>%
  distinct(Species) %>%
  pull()
s_means <- all_spp %>%
  filter(Species != "Atlantic Tripletail") %>%
  group_by(Species) %>%
  summarize(d13C = mean(d13C), d15N=mean(d15N), d34S = mean(d34S)) %>%
  #need to drop column of species names
  select(-Species) %>%
  as.matrix()
s_sds <- all_spp %>%
  filter(Species != "Atlantic Tripletail") %>%
  group_by(Species) %>%
  summarize(d13C = sd(d13C), d15N=sd(d15N), d34S = sd(d34S)) %>%
  #need to drop column of species names
  select(-Species) %>%
  as.matrix()
groups <- all_spp %>%
  filter(Species == "Atlantic Tripletail") %>%
  select(Estuary) %>%
  as.matrix()

# Load it in 
simmr_in_area = simmr_load(mixtures=Tripletail,
                      source_names=s_names,
                      source_means=s_means,
                      source_sds=s_sds,
                      group=groups)
simmr_in_area #gives you a quick summary; should give an error if you've mis-specified something

# Want help?
#help(simmr_load)

# Default iso-space
plot(simmr_in_area, tracers = c(1,2)) #plot C & N
plot(simmr_in_area, tracers = c(1,3)) #plot C & S
plot(simmr_in_area, tracers = c(2,3)) #plot N & S

# Better iso-space plots, for all pairs of isotopes
plot(simmr_in_area, tracers = c(1,2),
     xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")), 
     title='Isospace plot of Atlantic Tripletail and diet sources')

plot(simmr_in_area, tracers = c(1,3),
     xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     ylab=expression(paste(delta^34, "S (\u2030)",sep="")), 
     title='Isospace plot of Atlantic Tripletail and diet sources')

plot(simmr_in_area, tracers = c(2,3),
     xlab=expression(paste(delta^15, "N (\u2030)",sep="")), 
     ylab=expression(paste(delta^34, "S (\u2030)",sep="")), 
     title='Isospace plot of Atlantic Tripletail and diet sources')

# Getting help on this is a bit fiddlier
help(plot.simmr_input)

# Run simmr!
simmr_out_area = simmr_mcmc(simmr_in_area)
simmr_out_area
# Look at output ----------------------------------------------------------

# Always start with convergence
summary(simmr_out_area, type='diagnostics', group = 1)
summary(simmr_out_area, type='diagnostics', group = 2) #etc. here
summary(simmr_out_area, type='diagnostics', group = 3) 
# If these are close to 1 (e.g. less than 1.1) you can continue; if not, try longer run of MCMC

# Now you can explore ...
help(summary.simmr_output)
summary(simmr_out_area, type = 'statistics')
summary(simmr_out_area, type = 'quantiles')
summary(simmr_out_area, type = 'correlations') # good if you have overlapping sources (those should be correlated)
#  typically negative because as one goes up, the other goes down; if positive then both props go up to keep consumers in middle of sources in plot
# Plot the output
help(plot.simmr_output)
plot(simmr_out_area, type = 'histogram', group = 1)
plot(simmr_out_area, type = 'histogram', group = 2)
plot(simmr_out_area, type = 'histogram', group = 3)
plot(simmr_out_area, type = 'density')
plot(simmr_out_area, type = 'boxplot')
plot(simmr_out_area, type = 'matrix')

# To look at prior vs. posterior distributions
#  How to interpret:
#  You don't necessarily want the prior and posterior to line up or converge
#  If you’re using the default vague priors you want them to be different. 
#  If they look too similar then the data isn’t overwhelming the vague prior.
#  If you’ve got informative priors then you want them to be reasonably similar
#  because you want the prior to play a role in the estimation of the parameters.
#  A. Parnell says my plot with the vague default priors looks ok to him
prior_viz(simmr_out_area)

# plot the posterior predictive power 
posterior_predictive(simmr_out_area) # this fails saying differing number of rows 192 & 288 in matrix y_post_pred_ci???

# Do more with the output -------------------------------------------------

# Do more with the output -------------------------------------------------

# Compare source proportions between groups
compare_groups(simmr_out_area ,source = "Atlantic Bumper", groups = c(1,3))
compare_groups(simmr_out_area ,source = "Gulf Menhaden", groups = c(1,3))
compare_groups(simmr_out_area ,source = "Brown Shrimp", groups = c(1,3))
compare_groups(simmr_out_area ,source = "Northern White Shrimp", groups = c(1,3))

