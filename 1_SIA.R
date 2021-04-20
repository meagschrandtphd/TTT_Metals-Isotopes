# April 20, 2021

# Start Looking at SIA data, after taking A. Parnell and A. Jackson's online
#  course for using simmr and MixSIAR

library(tidyverse)

ttt <- read.csv("data/TTT_fish_SI.csv", header = T, stringsAsFactors = F)
names(ttt)

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
library(simmr)

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
plot(simmr_in)

# Better iso-space plot
plot(simmr_in,
     xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")), 
     title='Isospace plot of Atlantic Tripletail')

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
# Do more with the output -------------------------------------------------

# Do more with the output -------------------------------------------------

# Compare two sources
compare_sources(simmr_out,source_names=c('Brown Shrimp','Northern White Shrimp'))

# Compare multiple sources
compare_sources(simmr_out,source_names=c('Brown Shrimp','Northern White Shrimp', 'Atlantic Croaker'))

# -----------------------------------------------------------------------
# Multiple groups ---------------------------------------------------------
# -----------------------------------------------------------------------