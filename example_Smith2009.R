# set working directory, load data, packages and source files #
#_____________________________________________________________#

rm(list = ls())

setwd("~/Dropbox/PhD Nottingham/Academic/ppc paper additional work/Files for github")
source(file = "functions.R")
load(file="example_Smith2009.RData")




# distance method #
#_________________#

distance.method(R.obs = R.obs, R.STAR = R.STAR)




# position-time method #
#______________________#

position.time.method(R.obs = R.obs, R.STAR = R.STAR)
