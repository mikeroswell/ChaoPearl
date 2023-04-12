# script to generate "piano plots" and "slug plots" to visually analyze
# performance of Chao et al. uncertainty estimates for Hill diversity.

# Using highly idealized case
# Random sampling of individuals from community with known, finite diversity and
# infinite abundance

# This recapitulates some work from 2019 and 2020 by Dushoff and Roswell but
# here uses a cleaner/more stable codebase and sits in an easy to navigate repo.

# First, load libraries
library(dplyr)
library(furrr)
library(purrr)
library(MeanRarity)
library(checkPlotR)

# and load functions
source("scripts/diversity_checkplot_helpers.R")

# Generate communities

SADs_list<-map(c("lnorm", "gamma"), function(distr){
    map(c(100, 200), function(rich){
        # map(c(0.05, .15,.25,.5,.75,.85), function(simp_Prop){
        
        map(c(0.01,0.05, .15,.25,.5,.75,.85), function(simp_Prop){ #this is the special sads loop
            MeanRarity::fit_SAD(rich = rich, simpson = simp_Prop*(rich-1)+1, dstr = distr)
        })
    })
})

nc <- 9 # set number of cores
future::plan(strategy = "multisession", workers = nc)

cpDat_full <- map_dfr(flatten(flatten(SADs_list)), function(SAD){
    map_dfr(floor(10^seq(2, 5, 0.25)), function(ss){
        cpDat <- checkplot_inf(SAD, l = 0, inds = ss, reps = 10000)
        data.frame(cpDat
                   , dist = SAD$distribution_info[1]
                   , param = SAD$distribution_info[2]
                   , rich = SAD$community_info[1]
                   , HillShannon = SAD$community_info[2]
                   , HillSimpson = SAD$community_info[3]
        )
        })
})

# save the data
data.table::fwrite(cpDat_full, "data/Hill_Shannon_checkplot_data.csv", row.names = FALSE)

# read data back in
cpDat_full <- read.csv("data/Hill_Shannon_checkplot_data.csv")

# example of piano plots
cpDat_full %>% 
    filter(param == unique(param)[[3]]) %>%     checkPlot(facets = 13) + 
    facet_wrap(~inds) + 
    theme_classic() +
    ylim(c(0, 5000))

# and slugplot
cpDat_full %>% 
    filter(param == unique(param)[[3]]
           , inds == 1000) %>% 
    mutate(est = chaoest) %>% 
    rangePlot(title = "slugplot for Chao & Jost Hill-Shannon CI estimator.\nN= 1000; lognormal SAD with 100 species and Hill-Simpson of 15.85")
