###

#Rhine snow - General settings
#Erwin Rottler, University of Potsdam
#Spring, 2019

###

#packages----

# devtools::install_github('ERottler/alptempr')

remove.packages("meltimr")
devtools::install_github('ERottler/meltimr')
library("meltimr")

pacman::p_load(ncdf4, ncdf4.helpers, PCICt, dplyr, readr, tidyr, rgeos, ggplot2, 
               sp, viridis, rgdal, leaflet, ggmap, zoo, zyp, alptempr, lmomco, 
               raster, foreach, rfs, dismo, XML, parallel, doParallel, Lmoments,
               shape, devtools, pbapply, profvis, RColorBrewer, viridis, Rcpp, rEchseSnow,
               Rlibeemd, xts, emdbook, rfs, meltimr)

#directories----

#set base direcoty
base_dir <- "u:/RhineFlow/rhine_snow/"

#gridded climate data
# file_dir <- "e:/mhm_data/04_Daten/lobith_6435060/input/"
file_dir <- "d:/nrc_user/rottler/toErwin1/6435060/"

#GRDC discharge data
grdc_dir <- "d:/nrc_user/rottler/GRDC_DAY/"

#functions----

#load functions
source(paste0(base_dir, "R/melTim/2_functions.R"))


#cluster----

#Cluster for parallel computing

#stop cluster
stopCluster(my_clust)

n_cores <- 45 #number of cores used for parallel computing

#Make cluster for parallel computing
my_clust <- makeCluster(n_cores)
clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr, lmomco, ncdf4, rEchseSnow, sp, raster, betareg, rfs, meltimr))
registerDoParallel(my_clust)


