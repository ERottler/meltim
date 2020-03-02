###

#Rhine snow - General settings
#Erwin Rottler, University of Potsdam
#Spring, 2019

###

#packages----

# remove.packages("meltimr")
# devtools::install_github('ERottler/meltimr')
# devtools::install_github('ERottler/rfs')
# devtools::install_github("laubblatt/phaselag")

# #ECHSE snow functions in package
# Rcpp::sourceCpp(paste0(base_dir, "R/meltim/echse_snow.cpp"))
# Rcpp::Rcpp.package.skeleton(name = "rEchseSnow", cpp_files = paste0(base_dir, "R/meltim/echse_snow.cpp"))
# install.packages(paste0(base_dir, "R/meltim/rEchseSnow"), repos=NULL, type="source")
# library(rEchseSnow)

pacman::p_load(ncdf4, ncdf4.helpers, PCICt, dplyr, readr, tidyr, rgeos, ggplot2, 
               sp, viridis, rgdal, leaflet, ggmap, zoo, zyp, alptempr, lmomco, 
               raster, foreach, rfs, dismo, XML, parallel, doParallel, Lmoments,
               shape, devtools, pbapply, profvis, RColorBrewer, viridis, Rcpp, rEchseSnow,
               Rlibeemd, xts, emdbook, rfs, meltimr, readr, tmap, sf, hydroGOF)
# library(phaselag)

#directories----

#set base direcoty
base_dir <- "U:/rhine_snow/"

#gridded climate data
# file_dir <- "e:/mhm_data/04_Daten/lobith_6435060/input/"
file_dir <- "d:/nrc_user/rottler/toErwin1/6435060/"

#snow cover data from EURAC
scf_eurac_dir <- "D:/nrc_user/rottler/SCF_data/snow_eurac/zenodo_02_cloudremoval.tar/zenodo_02_cloudremoval/05_temporal_complete_max10d/" 

#snow cover data from DLR
scf_dlr_dir <- "D:/nrc_user/rottler/SCF_data/snow_dlr/SnowPack_DLR.tar/SnowPack_DLR/" 

#catchments
ezg_dir <- "D:/nrc_user/rottler/basin_data/EZG_Schweiz_BAFU/"


#cluster----

#Cluster for parallel computing

# stopCluster(my_clust)

n_cores <- 30 #number of cores used for parallel computing

#Make cluster for parallel computing
my_clust <- makeCluster(n_cores)
clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr, lmomco, ncdf4, rEchseSnow, sp, raster, betareg, viridis, rfs, meltimr))
registerDoParallel(my_clust)

#save_files----

save(date_snow, grid_points_d_in, snows_d, temps_d, precs_d, elevs_d,
     basin, basin_buf,
     
     scf_eurac, 
     

     file = "U:/rhine_snow/R/draft_snow_17_12.RData")


# 
# load(file = "U:/rhine_snow/R/draft_snow.RData")
# 
# save(snows_d_band, my_elev_bands, date_snow, svolu_d_band, 
#      smea_band, vmea_band, sslo_band, vslo_band, vdif_band, vdis_band,
#      tmea_band, tmea_band_mea, tslo_band, tslo_band_mea,
#      pmea_band, pmea_band_mea, pslo_band, pslo_band_mea,
#      meta_grid_bands, sno_vol_basin, prec_basin,
#      
#      file = "U:/rhine_snow/R/draft_snow_17_12.RData")
#      # file = "U:/rhine_snow/R/draft_snow_100.RData")


# smea_band_diep <- smea_band
# sslo_band_diep <- sslo_band
# vmea_band_diep <- vmea_band
# vslo_band_diep <- vslo_band
# vdif_band_diep <- vdif_band
# vdis_band_diep <- vdis_band
# my_elev_bands_diep <- my_elev_bands
# basin_84_diep <- basin_84
# snows_d_band_diep <- snows_d_band
# 
# smea_band_reus <- smea_band
# sslo_band_reus <- sslo_band
# vmea_band_reus <- vmea_band
# vslo_band_reus <- vslo_band
# vdif_band_reus <- vdif_band
# vdis_band_reus <- vdis_band
# my_elev_bands_reus <- my_elev_bands
# basin_84_reus <- basin_84
# snows_d_band_reus <- snows_d_band
# 
# smea_band_aare <- smea_band
# sslo_band_aare <- sslo_band
# vmea_band_aare <- vmea_band
# vslo_band_aare <- vslo_band
# vdif_band_aare <- vdif_band
# vdis_band_aare <- vdis_band
# my_elev_bands_aare <- my_elev_bands
# basin_84_aare <- basin_84
# snows_d_band_aare <- snows_d_band
# 
# smea_band_base <- smea_band
# sslo_band_base <- sslo_band
# vmea_band_base <- vmea_band
# vslo_band_base <- vslo_band
# vdif_band_base <- vdif_band
# vdis_band_base <- vdis_band
# my_elev_bands_base <- my_elev_bands
# basin_84_base <- basin_84
# snows_d_band_base <- snows_d_band




# load("U:/rhine_snow/R/meltim_snow.Rdata")
# 
# save(smea_band_diep, sslo_band_diep, vmea_band_diep, vslo_band_diep, vdif_band_diep, vdis_band_diep,
#      snows_d_band_diep, my_elev_bands_diep, basin_84_diep,
#      smea_band_reus, sslo_band_reus, vmea_band_reus, vslo_band_reus, vdif_band_reus, vdis_band_reus,
#      snows_d_band_reus, my_elev_bands_reus, basin_84_reus,
#      smea_band_aare, sslo_band_aare, vmea_band_aare, vslo_band_aare, vdif_band_aare, vdis_band_aare,
#      snows_d_band_aare, my_elev_bands_aare, basin_84_aare,
#      smea_band_base, sslo_band_base, vmea_band_base, vslo_band_base, vdif_band_base, vdis_band_base,
#      snows_d_band_base, my_elev_bands_base, basin_84_base,
#      date_snow,
#      
#      file = "U:/rhine_snow/R/meltim_snow.RData", version = 2)
