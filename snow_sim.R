###

#Snow simulations + Snow/meteo analysis on sub-basin scale
#Erwin Rottler, University of Potsdam

###

sta_yea_bas <- 1950
end_yea_bas <- 2014
# basin_sel <- "basel"        # basel, alp_rhine, reuss, aare, 
high_stat_thresh <- 2000 #1900
middle_stat_thresh <- 1000 #900

snow_params <- read.table(paste0(base_dir, "R/melTim/snow_param.txt"), header = T, sep = ";")

#Update snow parameter with calibrated values
dds_log <-  read.table(file = paste0(base_dir, "R/meltim/melt_calib/dds.log"), sep = "\t", header = T)

best_run_ind <- min_na(which(dds_log$objective_function == min_na(dds_log$objective_function)))

# snow_params$a0                <- dds_log[best_run_ind, which(colnames(dds_log) == "a0")]
# snow_params$a1                <- dds_log[best_run_ind, which(colnames(dds_log) == "a1")]
snow_params$kSatSnow          <- dds_log[best_run_ind, which(colnames(dds_log) == "kSatSnow")]
snow_params$densDrySnow       <- dds_log[best_run_ind, which(colnames(dds_log) == "densDrySnow")]
snow_params$specCapRet        <- dds_log[best_run_ind, which(colnames(dds_log) == "specCapRet")]
snow_params$emissivitySnowMin <- dds_log[best_run_ind, which(colnames(dds_log) == "emissivitySnowMin")]
snow_params$emissivitySnowMax <- dds_log[best_run_ind, which(colnames(dds_log) == "emissivitySnowMax")]
snow_params$tempAir_crit      <- dds_log[best_run_ind, which(colnames(dds_log) == "tempAir_crit")]
snow_params$albedoMin         <- dds_log[best_run_ind, which(colnames(dds_log) == "albedoMin")]
snow_params$albedoMax         <- dds_log[best_run_ind, which(colnames(dds_log) == "albedoMax")]
snow_params$agingRate_tAirPos <- dds_log[best_run_ind, which(colnames(dds_log) == "agingRate_tAirPos")]
snow_params$agingRate_tAirNeg <- dds_log[best_run_ind, which(colnames(dds_log) == "agingRate_tAirNeg")]
snow_params$weightAirTemp     <- dds_log[best_run_ind, which(colnames(dds_log) == "weightAirTemp")]

#data prep----

#Load ncdf E-OBS gridded datasets
nc_temp_file <- paste0(file_dir, "meteo/tavg.nc")
nc_prec_file <- paste0(file_dir, "meteo/pre.nc")
nc_petr_file <- paste0(file_dir, "meteo_HS/pet.nc")

nc_temp <- ncdf4::nc_open(nc_temp_file)
nc_prec <- ncdf4::nc_open(nc_prec_file)
nc_petr <- ncdf4::nc_open(nc_petr_file)

#Projections used later
crswgs84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
epsg3035 <- sp::CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 
                    +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

#Load DEM
# dem = raster(paste0(file_dir, "morph/dem.asc"), crs = epsg3035)
dem = raster(paste0(base_dir, "data/eu_dem/processed/eu_dem_1000.tif"))

#Load basin boundaries (shapefile delineated beforehand using Q-GIS)
# my_basin <- paste0(base_dir, "data/basins_shp/", basin_sel, ".shp")
# basin <- rgdal::readOGR(dsn = my_basin)
basins <-  rgdal::readOGR(dsn = "D:/nrc_user/rottler/basin_data/EZG_Schweiz_BAFU/ezg_kombiniert.shp", encoding = "UTF8")
basin <- spTransform(basins[basins@data$Ort == "Basel, Rheinhalle",], CRS = crs(dem, asText = T))
basin_buf <- buffer(basin, width = 8000)

#corp DEM sub-basin area
dem_cro <- raster::crop(dem, extent(basin))
dem_sub <- mask(dem_cro, basin)

#get elevations of cells cropped dem
dem_ele_NA <- dem_sub@data@values
dem_ele <- dem_ele_NA[!is.na(dem_ele_NA)]

#Area of basin in m2
area_m2 <- area(basin)

# plot(dem_sub)
# plot(basin, add =T)
# hist(dem_ele, nclass = 100)

#Transform projection WGS84
basin_buf_84 <- sp::spTransform(basin_buf, CRS = crswgs84)
# dem_84    <- raster::projectRaster(dem, crs = crswgs84)

#get lat/lon/time of .nc meteo data
lon <- ncdf4::ncvar_get(nc_temp, varid = "lon2D")
lat <- ncdf4::ncvar_get(nc_temp, varid = "lat2D")
date <- as.Date(as.character(ncdf4.helpers::nc.get.time.series(nc_temp, time.dim.name = "time")))

#define date start index and count for extraction from .nc file
date_min_index <- which(date == paste0(sta_yea_bas, "-01-01"))
date_max_index <- which(date == paste0(end_yea_bas, "-12-31"))
date_count <- date_max_index - date_min_index +1
meteo_date <- date[date_min_index:date_max_index]

#extract coordinates from basin shapefile
basin_coords <- extract_coords(basin_buf_84)

#get min/max of lat/lon for extraction from .nc file
lon_min <- min(basin_coords[ ,1]) - 0.5
lon_max <- max(basin_coords[ ,1]) + 0.5
lat_min <- min(basin_coords[ ,2]) - 0.5
lat_max <- max(basin_coords[ ,2]) + 0.5

lon_min_index <- which(abs(lon - lon_min) == min(abs(lon - lon_min)), arr.ind = T)[1,1]
lat_min_index <- which(abs(lat - lat_max) == min(abs(lat - lat_max)), arr.ind = T)[1,2]

lon_max_index <- which(abs(lon - lon_max) == min(abs(lon - lon_max)), arr.ind = T)[1,1]
lat_max_index <- which(abs(lat - lat_min) == min(abs(lat - lat_min)), arr.ind = T)[1,2]

lon_count <- abs(lon_max_index - lon_min_index) 
lat_count <- abs(lat_max_index - lat_min_index) 

#extract meteo data from .nc file using previously defined course cube extensions

print(paste0(Sys.time(), " Extract temperature from ncdf-file."))

temps_cube <- ncdf4::ncvar_get(nc_temp, start =c(lon_min_index, lat_min_index, date_min_index), 
                               count = c(lon_count, lat_count, date_count), varid = "tavg")

print(paste0(Sys.time(), " Extract precipitation from ncdf-file."))

precs_cube <- ncdf4::ncvar_get(nc_prec, start =c(lon_min_index, lat_min_index, date_min_index), 
                               count = c(lon_count, lat_count, date_count), varid = "pre")

#get lat/lon values of extracted cube
lon2D <- lon[lon_min_index : (lon_min_index+lon_count-1), lat_min_index : (lat_min_index+lat_count-1)]
lat2D <- lat[lon_min_index : (lon_min_index+lon_count-1), lat_min_index : (lat_min_index+lat_count-1)]

#spatial grid points from lat/lon info
grid_points_cube_84 <-  sp::SpatialPoints(data.frame(lon = c(lon2D), lat = c(lat2D)), proj4string =  crswgs84)
grid_points_cube     <- sp::spTransform(grid_points_cube_84, CRS = crs(basin, asText = T))

#get grid points inside watershed
inside <- !is.na(sp::over(grid_points_cube, as(basin_buf, "SpatialPolygons")))
grid_points <- grid_points_cube[which(inside == T)]
grid_points_84 <- sp::spTransform(grid_points, CRS = crs(grid_points_cube_84, asText = T))

#get index in cube from points inside sub-basin
cube_index_col <- sapply(grid_points_84@coords[,1], get_cube_index_col)
cube_index_row <- sapply(grid_points_84@coords[,1], get_cube_index_row)

#get time series from grid points inside basin (with buffer)
#temperature
for (i in 1:length(cube_index_col)) {
  
  temp_sing <- temps_cube[cube_index_col[i], cube_index_row[i], ]
  
  if(i == 1){
    temps <- temp_sing
  }else{
    temps <- cbind(temps, temp_sing)
  }
}
#precipitation
for (i in 1:length(cube_index_col)) {
  
  precs_sing <- precs_cube[cube_index_col[i], cube_index_row[i], ]
  
  if(i == 1){
    precs <- precs_sing
  }else{
    precs <- cbind(precs, precs_sing)
  }
}

#sub-basin average values for rangking of weather types
temp_basin <- apply(temps, 1, med_na)
prec_basin <- apply(precs, 1, mea_na) 

#Elevation grid points: 5km squared buffer
my_elev_buff <- function(point_index){
  
  elev_buff(point_in = grid_points[point_index], dem_in = dem)

  }

elevs <- foreach(i = 1:length(grid_points), .combine = 'c') %dopar% {
  
  my_elev_buff(i) #[m]
  
}
elevs_ord <- elevs[order(elevs)]

# #elevation of grip points from DEM
# elevs <- raster::extract(dem, grid_points)

#Syntetic meta data info for grid points to use alptempr functions
HS_amount <- length(which(elevs_ord > high_stat_thresh))
MS_amount <- length(which(elevs_ord > middle_stat_thresh)) - HS_amount
LS_amount <- length(elevs_ord) - MS_amount - HS_amount

meta_grid <- data.frame(stn = paste0("point",1:length(grid_points)),
                        alt = elevs_ord,
                        category = c(rep("low", LS_amount), rep("middle", MS_amount), rep("high", HS_amount)),
                        data_qual = rep("quality-checked", length(grid_points)),
                        clim_reg = rep("Jura", length(grid_points)))

#Radiation data for snow simulations
#Measured daily solar radiation from station Napf
#Mean annual cycle as simplified input for snow simulations

radi_sno <- read.table(paste0(base_dir, "data/idaweb/order62894/order_62894_data.txt"), sep = ";", skip = 2, header = T)
radi_sno$date <- as.Date(strptime(radi_sno$time, "%Y%m%d", tz="UTC"))
radi_sno$radi <- radi_sno$gre000d0

radi_snow_day <- ord_day(data_in = radi_sno$gre000d0,
                         date= radi_sno$date,
                         start_y = 1981,
                         end_y = 2017)

radi_mea <- apply(radi_snow_day, 2, mea_na)

radi_mea_smo <- smoothFFT(radi_mea, sd = 15)
# plot(radi_mea_smo, type = "l")

#Syntetic radiation time series based on mean values
radi_mea_seri <- rep( c(radi_mea_smo, radi_mea_smo[365], rep(radi_mea_smo, 3)), 50)[1:nrow(temps)]


#downscale----

#Downscale temperature grid using simple lapse-rate based approach
#lapse rate in basin on daily bases used

#downscale points
down_points <- function(lon_lat_point_in){
  
  res_new <- 1000 # new desirted resolution in [m]
  
  d_00 <- lon_lat_point_in
  d_01 <- d_00 + res_new * 1
  d_02 <- d_00 + res_new * 2
  d_03 <- d_00 - res_new * 1
  d_04 <- d_00 - res_new * 2
  d_05 <- d_00
  d_05[1] <- d_05[1] + res_new * 1
  d_05[2] <- d_05[2] - res_new * 1
  d_06 <- d_00
  d_06[1] <- d_06[1] + res_new * 2
  d_06[2] <- d_06[2] - res_new * 2
  d_07 <- d_00
  d_07[1] <- d_07[1] - res_new * 1
  d_07[2] <- d_07[2] + res_new * 1
  d_08 <- d_00
  d_08[1] <- d_08[1] - res_new * 2
  d_08[2] <- d_08[2] + res_new * 2
  d_09 <- d_00
  d_09[1] <- d_09[1] + res_new * 1
  d_09[2] <- d_09[2] + res_new * 0
  d_10 <- d_00
  d_10[1] <- d_10[1] + res_new * 2
  d_10[2] <- d_10[2] + res_new * 0
  d_11 <- d_00
  d_11[1] <- d_11[1] - res_new * 1
  d_11[2] <- d_11[2] + res_new * 0
  d_12 <- d_00
  d_12[1] <- d_12[1] - res_new * 2
  d_12[2] <- d_12[2] + res_new * 0
  d_13 <- d_00
  d_13[1] <- d_13[1] + res_new * 0
  d_13[2] <- d_13[2] - res_new * 1
  d_14 <- d_00
  d_14[1] <- d_14[1] + res_new * 0
  d_14[2] <- d_14[2] - res_new * 2
  d_15 <- d_00
  d_15[1] <- d_15[1] + res_new * 0
  d_15[2] <- d_15[2] + res_new * 1
  d_16 <- d_00
  d_16[1] <- d_16[1] + res_new * 0
  d_16[2] <- d_16[2] + res_new * 2
  d_17 <- d_00
  d_17[1] <- d_17[1] + res_new * 2
  d_17[2] <- d_17[2] + res_new * 1
  d_18 <- d_00
  d_18[1] <- d_18[1] + res_new * 2
  d_18[2] <- d_18[2] - res_new * 1
  d_19 <- d_00
  d_19[1] <- d_19[1] + res_new * 1
  d_19[2] <- d_19[2] - res_new * 2
  d_20 <- d_00
  d_20[1] <- d_20[1] - res_new * 1
  d_20[2] <- d_20[2] - res_new * 2
  d_21 <- d_00
  d_21[1] <- d_21[1] - res_new * 2
  d_21[2] <- d_21[2] - res_new * 1
  d_22 <- d_00
  d_22[1] <- d_22[1] - res_new * 2
  d_22[2] <- d_22[2] + res_new * 1
  d_23 <- d_00
  d_23[1] <- d_23[1] - res_new * 1
  d_23[2] <- d_23[2] + res_new * 2
  d_24 <- d_00
  d_24[1] <- d_24[1] + res_new * 1
  d_24[2] <- d_24[2] + res_new * 2
  
  d_points_raw <- rbind(d_00, d_01, d_02, d_03, d_04, d_05, d_06, d_07, d_08, d_09, d_10, d_11, d_12,
                        d_13, d_14, d_15, d_16, d_17, d_18, d_19, d_20, d_21, d_22, d_23, d_24)
  
  
  #spatial grip points from lat/lon info
  d_points <- sp::SpatialPoints(data.frame(lon = d_points_raw[, 1], lat = d_points_raw[, 2]), proj4string =  sp::CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000
                    +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  
  return(d_points)
  
}

#Lapse rate on daily basis from all data points selected
my_lin_trend <- function(data_in){
  
  lin_trend(data_in, index_in = elevs_ord)
  
}

temps_ord <- temps[, order(elevs)]

lapse_temp <- foreach(i = 1:nrow(temps_ord), .combine = 'c') %dopar% {
  
  my_lin_trend(temps_ord[i, ]) #[°C/m]
  
}
lapse_prec <- rep(0, nrow(precs)) #no lapse rate approach for precipitation; rain stays the same

#downscale meteo on 1 km grid
f_temps_laps <- function(data_in){
  
  f_laps_mod <- function(index_in){
    
    data_laps <- data_in + lapse_temp* d_points_elevs_dif[index_in]
    return(data_laps)
    
  }
  
  temps_down <- sapply(1:length(d_points_elevs_dif), f_laps_mod)
  
  return(temps_down)
  
}

f_precs_laps <- function(data_in){
  
  f_laps_mod <- function(index_in){
    
    data_laps <- data_in + lapse_prec* d_points_elevs_dif[index_in]
    return(data_laps)
    
  }
  
  precs_down <- sapply(1:length(d_points_elevs_dif), f_laps_mod)
  return(precs_down)
  
}

grid_points_d <- foreach(i = 1:length(grid_points), .combine = 'rbind') %dopar% {
  
  d_points <- down_points(grid_points@coords[i, ])
  
  d_points@coords
  
}

#spatial grid points from lat/lon info
grid_points_d     <- sp::SpatialPoints(data.frame(lon = grid_points_d[, 1], lat = grid_points_d[, 2]), proj4string =  epsg3035)
grid_points_d_84  <- sp::spTransform(grid_points_d, CRS = crswgs84)

temps_d <- foreach(i = 1:ncol(temps), .combine = 'cbind') %dopar% {
  
  d_points <- down_points(grid_points@coords[i, ])
  
  d_points_elevs <- raster::extract(dem, d_points)
  
  d_points_elevs_dif <- d_points_elevs - elevs[i]
  
  f_temps_laps(temps[, i]) #[°C]
  
}

precs_d <- foreach(i = 1:ncol(precs), .combine = 'cbind') %dopar% {
  
  d_points <- down_points(grid_points@coords[i, ])
  
  d_points_elevs <- raster::extract(dem, d_points)
  
  d_points_elevs_dif <- d_points_elevs - elevs[i]
  
  f_precs_laps(precs[, i]) #[°C]
  
}

elevs_d <-  raster::extract(dem, grid_points_d)

#snow_simu----

# #if no downscaling > rename
# temps_d <- temps
# precs_d <- precs
# elevs_d <- elevs

#Sonw simulations using physically-based ECHSE snow routine

#Simulations using stand alone ECHSE snow model
if(F){
  
  # for(i in 1:ncol(precs)){
  for(i in 1:ncol(precs)){  
    print(paste(Sys.time(), "Processing grid point", i, "out of", ncol(precs)))
    
    #Rainfall [mm]
    rain <- precs[, i]
    
    #Temperature [?C]
    temp = temps[, i]
    
    #Raidation [W/m2]
    radi=rep(100, nrow(precs))
    
    #Humidity [%]
    humi=rep(70, nrow(precs))
    
    #Wind [m/s]
    wind=rep(1, nrow(precs))
    
    #Could Coverage [-]
    cloud=rep(0.5, nrow(precs))
    
    #Air Pressure [hPa]
    pressAir = rep(1000, nrow(precs))
    
    ### Export table
    output = cbind(temp, rain, radi, humi, pressAir, wind, cloud)
    write.table(file = paste0(base_dir, "snow_sim/snowAlone/input/input_meteo.dat"), x=output, 
                na = "-9999", sep="\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
    
    system(snow_exe)
    
    swe <- read.table(paste0(base_dir, "snow_sim/snowAlone/output/swe.out"), header = T)
    
    if(i == 1){
      snows <- swe$swe
    }else{
      snows <- cbind(snows, swe$swe)
    }
    
  }
  
  #Save data as .Rdata
  save(file = paste0(base_dir, "snow_sim/data_sim/", "snows_", basin_sel, ".RData"), list="snows")
  
}

#Simulations using cpp functions in R
block_size <- 1000
block_stas <- c(1, seq(block_size+1, ncol(temps_d), by = block_size))
block_ends <- c(seq(block_size, ncol(temps_d), by = block_size), ncol(temps_d))

#disturb input signal
if(F){
  # temps_d <- temps_d - 2
  precs_d <- precs_d * 0.8
}

for(b in 1:length(block_stas)){
  
  print(paste(Sys.time(),"Snow simulations", "Block:", b, "out of", length(block_stas)))
  
  temps_simu <- temps_d[, block_stas[b]:block_ends[b]]
  precs_simu <- precs_d[, block_stas[b]:block_ends[b]]
  
  snows <- foreach(k = 1:ncol(temps_simu), .combine = 'cbind') %dopar% {
    
    swe_sim   <- rep(NA, nrow(temps_simu))
    sec_sim   <- rep(NA, nrow(temps_simu))
    alb_sim   <- rep(NA, nrow(temps_simu))
    
    swe_init <- .0
    sec_init <- .0
    alb_init <- snow_params$albedoMax
    
    swe_sim[1] <- swe_init
    sec_sim[1] <- sec_init
    alb_sim[1] <- alb_init
    
    for(i in 2:nrow(temps_simu)){
      
      sim_out <- snowModel_inter(
        #Forcings
        precipSumMM = precs_simu[i, k],
        shortRad = radi_mea_seri[i],
        tempAir = temps_simu[i, k],
        pressAir = 1000,
        relHumid = 70,
        windSpeed = 1,
        cloudCoverage = 0.5,
        #Parameters
        precipSeconds = snow_params$precipSeconds,
        a0 = snow_params$a0,
        a1 = snow_params$a1,
        kSatSnow = snow_params$kSatSnow,
        densDrySnow = snow_params$densDrySnow,
        specCapRet = snow_params$specCapRet,
        emissivitySnowMin = snow_params$emissivitySnowMin,
        emissivitySnowMax = snow_params$emissivitySnowMax,
        tempAir_crit = snow_params$tempAir_crit,
        albedoMin = snow_params$albedoMin,
        albedoMax = snow_params$albedoMax,
        agingRate_tAirPos = snow_params$agingRate_tAirPos,
        agingRate_tAirNeg = snow_params$agingRate_tAirNeg,
        soilDepth = snow_params$soilDepth,
        soilDens = snow_params$soilDens,
        soilSpecHeat = snow_params$soilSpecHeat,
        weightAirTemp = snow_params$weightAirTemp,
        tempMaxOff = snow_params$tempMaxOff,
        tempAmpli = snow_params$tempAmpli,
        #States
        snowEnergyCont = sec_sim[i-1],
        snowWaterEquiv = swe_sim[i-1],
        albedo = alb_sim[i-1],
        #Outputs
        TEMP_MEAN = NA,
        TEMP_SURF = NA,
        LIQU_FRAC = NA,
        flux_R_netS = NA,
        flux_R_netL = NA,
        flux_R_soil = NA,
        flux_R_sens = NA,
        stoi_f_prec = NA,
        stoi_f_subl = NA,
        stoi_f_flow = NA,
        flux_M_prec = NA,
        flux_M_subl = NA,
        flux_M_flow = NA,
        rate_G_alb = NA
      )
      
      sec_sim[i] <- sim_out[1]
      swe_sim[i] <- sim_out[2]
      alb_sim[i] <- sim_out[3]
      
    }
    
    swe_sim
    
  }
  
  if(b == 1){
    snows_d <- snows
  }else{
    snows_d <- cbind(snows_d, snows)
  }
  
}  

#first year as warm up to get snow dynamics going
#remove first year
index_next_1Jan <- which(meteo_date == paste0(sta_yea_bas+1,"-01-01"))-1
snows_d <- snows_d[-(1:index_next_1Jan), ]
temps_d <- temps_d[-(1:index_next_1Jan), ]
precs_d <- precs_d[-(1:index_next_1Jan), ]
date_snow <- meteo_date[-(1:index_next_1Jan)]

#Get grid points inside watershed (calculations with buffer around)
sim_inside <- !is.na(sp::over(grid_points_d, as(basin, "SpatialPolygons")))
grid_points_d_in <- grid_points_d[which(sim_inside == T)]
grid_points_d_in_84 <- sp::spTransform(grid_points_d_in, CRS = crs(grid_points_cube_84, asText = T))

points_outside <- unique(which(sim_inside == F))
snows_d <- snows_d[, -points_outside]
temps_d <- temps_d[, -points_outside]
precs_d <- precs_d[, -points_outside]
elevs_d <- elevs_d[-points_outside]


plot(grid_points_d_in[1:30])
raster::extract(dem, grid_points_d_in[1:30])
elevs_d[1:10]

#analysis----

snow_max <- apply(snows_d, 2, max_na)
#length(which(snow_max > 3))

grid_m2 <- 1000*1000 #grid 1 km resolution

block_size <- 1000
block_stas <- c(1, seq(block_size+1, ncol(snows_d), by = block_size))
block_ends <- c(seq(block_size, ncol(snows_d), by = block_size), ncol(snows_d))

#Snow height: mean average
for(b in 1:length(block_stas)){
  
  snows_calc <- snows_d[, block_stas[b]:block_ends[b]]
  
  print(paste(Sys.time(),"Average mean swe", "Block:", b, "out of", length(block_stas)))
  
  smea_block <- foreach(i = 1:ncol(snows_calc), .combine = 'cbind') %dopar% {
    
    day_ana(snows_calc[, i], 
            date = date_snow, 
            start_year = sta_yea_bas+1, 
            end_year = end_yea_bas,
            break_day = 274,
            do_ma = F,
            window_width = 30, 
            method_ana = "mean") #[m]
    
  }
  
  if(b == 1){
    smea <- smea_block
  }else{
    smea <- cbind(smea, smea_block)
  }
  
}

# smea[, which(snows_d[nrow(snows_d), ] > 5)] <- NA #remove 'glacier' points

smea[, which(snow_max > 3)] <- NA #remove snow towers

#Snow volume: mean average
vmea <- smea * grid_m2

#Snow height: trends 30DMA
for(b in 1:length(block_stas)){

  snows_calc <- snows_d[, block_stas[b]:block_ends[b]]

  print(paste(Sys.time(),"Trends 30DMA swe", "Block:", b, "out of", length(block_stas)))
  sslo_block <- foreach(i = 1:ncol(snows_calc), .combine = 'cbind') %dopar% {

    day_ana(snows_calc[, i], 
            date = date_snow, 
            start_year = sta_yea_bas, 
            end_year = end_yea_bas,
            break_day = 274,
            do_ma = T,
            window_width = 30, 
            method_ana = "sens_slope") * 10 #[m/dec]

  }

  if(b == 1){
    sslo <- sslo_block
  }else{
    sslo <- cbind(sslo, sslo_block)
  }

}

sslo[, which(snow_max > 3)] <- NA #remove snow towers

#Snow volume: trends 30 DMA
vslo <- sslo * grid_m2

#Snow accumulation and melt water outflow
sv_diff <- function(snow_volume_in){
  
  sv_diff <- c(NA, diff(snow_volume_in))
  
  return(sv_diff)
  
}

snows_d_dif <- apply(snows_d, 2, sv_diff)

#Snow height diff: mean average
for(b in 1:length(block_stas)){
  
  snows_calc <- snows_d_dif[, block_stas[b]:block_ends[b]]
  
  print(paste(Sys.time(),"Average mean diff swe", "Block:", b, "out of", length(block_stas)))
  
  sdif_block <- foreach(i = 1:ncol(snows_calc), .combine = 'cbind') %dopar% {
    
    day_ana(snows_calc[, i], 
            date = date_snow, 
            start_year = sta_yea_bas, 
            end_year = end_yea_bas,
            break_day = 274,
            do_ma = F,
            window_width = 30, 
            method_ana = "mean") #[m]
    
  }
  
  if(b == 1){
    sdif <- sdif_block
  }else{
    sdif <- cbind(sdif, sdif_block)
  }
  
}

sdif[, which(snow_max > 3)] <- NA #remove snow towers

#Snow volume diff: mean average
vdif <- sdif * grid_m2

#Snow height diff: trends 30DMA
for(b in 1:length(block_stas)){

  snows_calc <- snows_d_dif[, block_stas[b]:block_ends[b]]

  print(paste(Sys.time(),"Trends 30 DMA diff swe", "Block:", b, "out of", length(block_stas)))
  sdis_block <- foreach(i = 1:ncol(snows_calc), .combine = 'cbind') %dopar% {

    day_ana(snows_calc[, i], 
            date = date_snow, 
            start_year = sta_yea_bas, 
            end_year = end_yea_bas,
            break_day = 274,
            do_ma = T,
            window_width = 30, 
            method_ana = "sens_slope") * 10#[m/dec]

  }

  if(b == 1){
    sdis <- sdis_block
  }else{
    sdis <- cbind(sdis, sdis_block)
  }

}

sdis[, which(snow_max > 3)] <- NA #remove snow towers

#Snow volume diff: trends 30 DMA
vdis <- sdis * grid_m2

#Temperature: mean average
for(b in 1:length(block_stas)){

  temps_calc <- temps_d[, block_stas[b]:block_ends[b]]

  print(paste(Sys.time(),"Average mean temperature", "Block:", b, "out of", length(block_stas)))

  tmea_block <- foreach(i = 1:ncol(temps_calc), .combine = 'cbind') %dopar% {

    day_ana(temps_calc[, i],
            date = date_snow,
            start_year = sta_yea_bas,
            end_year = end_yea_bas,
            break_day = 274,
            do_ma = F,
            window_width = 30,
            method_ana = "mean") #[°C]

  }

  if(b == 1){
    tmea <- tmea_block
  }else{
    tmea <- cbind(tmea, tmea_block)
  }

}

#Temperature: mean average 5km x 5km
block_size_5km <- 1000
block_stas_5km <- c(1, seq(block_size+1, ncol(temps), by = block_size))
block_ends_5km <- c(seq(block_size, ncol(temps), by = block_size), ncol(temps))

for(b in 1:length(block_stas_5km)){
  
  temps_calc <- temps[, block_stas_5km[b]:block_ends_5km[b]]
  
  print(paste(Sys.time(),"Average mean temperature", "Block:", b, "out of", length(block_stas_5km)))
  
  tmea_block <- foreach(i = 1:ncol(temps_calc), .combine = 'cbind') %dopar% {
    
    day_ana(temps_calc[, i],
            date = meteo_date,
            start_year = sta_yea_bas,
            end_year = end_yea_bas,
            break_day = 274,
            do_ma = F,
            window_width = 30,
            method_ana = "mean") #[°C]
    
  }
  
  if(b == 1){
    tmea_5km <- tmea_block
  }else{
    tmea_5km <- cbind(tmea_5km, tmea_block)
  }
  
}

#Temperature: 30 DMA trends
for(b in 1:length(block_stas)){

  temps_calc <- temps_d[, block_stas[b]:block_ends[b]]

  print(paste(Sys.time(),"30DMA trends temperature", "Block:", b, "out of", length(block_stas)))

  tslo_block <- foreach(i = 1:ncol(temps_calc), .combine = 'cbind') %dopar% {

    day_ana(temps_calc[, i],
            date = date_snow,
            start_year = sta_yea_bas,
            end_year = end_yea_bas,
            break_day = 274,
            do_ma = T,
            window_width = 30,
            method_ana = "sens_slope") * 10 #[°C/dec]

  }

  if(b == 1){
    tslo <- tslo_block
  }else{
    tslo <- cbind(tslo, tslo_block)
  }

}

#Precipitation: mean average
for(b in 1:length(block_stas)){

  precs_calc <- precs_d[, block_stas[b]:block_ends[b]]

  print(paste(Sys.time(),"Average mean precipitation", "Block:", b, "out of", length(block_stas)))

  pmea_block <- foreach(i = 1:ncol(precs_calc), .combine = 'cbind') %dopar% {

    day_ana(precs_calc[, i],
            date = date_snow,
            start_year = sta_yea_bas,
            end_year = end_yea_bas,
            break_day = 274,
            do_ma = F,
            window_width = 30,
            method_ana = "mean") #[mm]

  }

  if(b == 1){
    pmea <- pmea_block
  }else{
    pmea <- cbind(pmea, pmea_block)
  }

}

#Precipitation: 30 DMA trends
for(b in 1:length(block_stas)){

  precs_calc <- precs_d[, block_stas[b]:block_ends[b]]

  print(paste(Sys.time(),"30DMA trends precipitation", "Block:", b, "out of", length(block_stas)))

  pslo_block <- foreach(i = 1:ncol(precs_calc), .combine = 'cbind') %dopar% {

    day_ana(precs_calc[, i],
            date = date_snow,
            start_year = sta_yea_bas,
            end_year = end_yea_bas,
            break_day = 274,
            do_ma = T,
            window_width = 30,
            method_ana = "sens_slope")  * 10 #[mm/dec]

  }

  if(b == 1){
    pslo <- pslo_block
  }else{
    pslo <- cbind(pslo, pslo_block)
  }

}
# 
# #Evaportranspiration: median average
# for(b in 1:length(block_stas)){
# 
#   petrs_calc <- petrs_d[, block_stas[b]:block_ends[b]]
# 
#   print(paste(Sys.time(),"Average mean evapotranspiration", "Block:", b, "out of", length(block_stas)))
# 
#   emea_block <- foreach(i = 1:ncol(petrs_calc), .combine = 'cbind') %dopar% {
# 
#     day_ana(petrs_calc[, i],
#             date = meteo_date,
#             start_year = sta_yea_bas,
#             end_year = end_yea_bas,
#             do_ma = F,
#             window_width = 30,
#             method_ana = "mean") #[mm]
# 
#   }
# 
#   if(b == 1){
#     emea <- emea_block
#   }else{
#     emea <- cbind(emea, emea_block)
#   }
# 
# }
# 
# #Evapotranspiration: 30 DMA trends
# for(b in 1:length(block_stas)){
# 
#   petrs_calc <- petrs_d[, block_stas[b]:block_ends[b]]
# 
#   print(paste(Sys.time(),"30DMA trends evapotranspiration", "Block:", b, "out of", length(block_stas)))
# 
#   eslo_block <- foreach(i = 1:ncol(petrs_calc), .combine = 'cbind') %dopar% {
# 
#     day_ana(petrs_calc[, i],
#             date = meteo_date,
#             start_year = sta_yea_bas,
#             end_year = end_yea_bas,
#             do_ma = T,
#             window_width = 30,
#             method_ana = "sens_slope") * 10 #[mm/dec]
# 
#   }
# 
#   if(b == 1){
#     eslo <- eslo_block
#   }else{
#     eslo <- cbind(eslo, eslo_block)
#   }
# 
# }


#elev_ranges----

f_elev_bands <- function(data_in, elev_bands = my_elev_bands, func_aggr = "mean", meta_dat = elevs_d){
  
  for(i in 1:(length(elev_bands) - 1)){
    # print(i)
    data_points_range <- which(meta_dat > elev_bands[i] & meta_dat < elev_bands[i+1])
    
    if(length(data_points_range) == 1){
      data_range_sing <- data_in[, data_points_range]
    }else{
      if(func_aggr == "mean"){
        data_range_sing <- apply(data_in[, data_points_range], 1, mea_na)
      }
      if(func_aggr == "medi"){
        data_range_sing <- apply(data_in[, data_points_range], 1, med_na)
      }
      if(func_aggr == "sum"){
        data_range_sing <- apply(data_in[, data_points_range], 1, sum_na)
      }
    }
    
    
    if(i == 1){
      
      data_range <- data_range_sing
      
    }else{
      
      data_range <- cbind(data_range, data_range_sing)
      
    }
    
  }
  
  return(data_range)
  
}

range(elevs_d)
# hist(elevs_d)
# range(elevs)
# my_elev_bands <- c(seq(400, 3000, 50), 4000) #Diepoldsau
my_elev_bands <- seq(250, 3200, 50) #Basel 1km
# my_elev_bands <- c(seq(250, 3000, 50), 4000) #Basel 5km
# my_elev_bands <- c(seq(350, 3150, 50), 4000) #Reuss
# my_elev_bands <- c(seq(350, 3550, 50), 4000) #Aare

smea_band <- f_elev_bands(data_in = smea, func_aggr = "mean")
sslo_band <- f_elev_bands(data_in = sslo, func_aggr = "mean")
vmea_band <- f_elev_bands(data_in = vmea, func_aggr = "sum")
vslo_band <- f_elev_bands(data_in = vslo, func_aggr = "sum")
vdif_band <- f_elev_bands(data_in = vdif, func_aggr = "sum")
vdis_band <- f_elev_bands(data_in = vdis, func_aggr = "sum")

snows_d[, which(snow_max > 3)] <- NA #remove snow towers

snows_d_band <- f_elev_bands(data_in = snows_d, func_aggr = "mean")
temps_d_band <- f_elev_bands(data_in = temps_d, func_aggr = "mean")
precs_d_band <- f_elev_bands(data_in = precs_d, func_aggr = "mean")

svolu_d_band <- f_elev_bands(data_in = snows_d, func_aggr = "sum")*grid_m2

#Syntetic meta data info for grid points to use alptempr functions
HS_amount <- length(which(my_elev_bands[-length(my_elev_bands)] > high_stat_thresh))
MS_amount <- length(which(my_elev_bands[-length(my_elev_bands)] > middle_stat_thresh)) - HS_amount
LS_amount <- length(my_elev_bands[-length(my_elev_bands)]) - MS_amount - HS_amount

meta_grid_bands <- data.frame(stn = paste0("band", 1:(length(my_elev_bands)-1)),
                              alt = my_elev_bands[-length(my_elev_bands)],
                              category = c(rep("low", LS_amount), rep("middle", MS_amount), rep("high", HS_amount)),
                              data_qual = rep("quality-checked", (length(my_elev_bands)-1)),
                              clim_reg = rep("Jura", (length(my_elev_bands)-1)))

tmea_band <- f_elev_bands(data_in = tmea, func_aggr = "mean")
colnames(tmea_band) <- paste0("band", 1:ncol(tmea_band))
tmea_band_mea <- apply(tmea_band, 2, mea_na)#annual average values

tslo_band <- f_elev_bands(data_in = tslo, func_aggr = "mean")
colnames(tslo_band) <- paste0("band", 1:ncol(tslo_band))
tslo_band_mea <- apply(tslo_band, 2, mea_na)#annual average values

pmea_band <- f_elev_bands(data_in = pmea, func_aggr = "mean")
colnames(pmea_band) <- paste0("band", 1:ncol(pmea_band))
pmea_band_mea <- apply(pmea_band, 2, mea_na)#annual average values

pslo_band <- f_elev_bands(data_in = pslo, func_aggr = "mean")
colnames(pslo_band) <- paste0("band", 1:ncol(pslo_band))
pslo_band_mea <- apply(pslo_band, 2, mea_na)#annual average values

tmea_band_5km <- f_elev_bands(data_in = tmea_5km, func_aggr = "mean", meta_dat = elevs)
colnames(tmea_band_5km) <- paste0("band", 1:ncol(tmea_band_5km))
tmea_band_mea_5km <- apply(tmea_band_5km, 2, mea_na)#annual average values

# emea_band <- f_elev_bands(data_in = emea, func_aggr = "mean")
# colnames(emea_band) <- paste0("band", 1:ncol(emea_band))
# emea_band_mea <- apply(emea_band, 2, mea_na)#annual average values
# 
# eslo_band <- f_elev_bands(data_in = eslo, func_aggr = "mean")
# colnames(eslo_band) <- paste0("band", 1:ncol(eslo_band))
# eslo_band_mea <- apply(eslo_band, 2, mea_na)#annual average values


#visu_snow_bands----

#Plot: Snow variables

plot_snow <- vdis_band #smea_band, sslo_band, vmea_band, vslo_band, vdif_band, vdis_band
col_zero <- T #color range center at zero

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(   46,74,105,135,166,196,227,258,288,319,349)-15

my_col <- colorRampPalette(c("white", viridis(9, direction = 1)[c(3,4)], "cadetblue3", "grey80",
                             "yellow2","gold", "orange2", "orangered2"))(200)

if(col_zero){
  
  n_max <- round(abs(max_na(plot_snow[, ])) / (max_na(plot_snow[, ]) + abs(min_na(plot_snow[, ]))), digits = 2) * 200
  n_min <- 200 - n_max
  
  cols_min <- colorRampPalette(c(viridis(9, direction = 1)[1:4], "cadetblue3", "grey90"))(n_min)
  cols_max <- colorRampPalette(c("grey90", "yellow2", "gold", "orange2", "orangered2"))(n_max)
  my_col <- c(cols_min, cols_max)
  
}

my_bre <- seq(min_na(plot_snow), max_na(plot_snow), length.out = length(my_col)+1)

par(mar = c(1.6, 3, 0.6, 0))

layout(matrix(c(1,1,1,1,1,1,1,2),
              1, 8), widths=c(), heights=c())

image(x = 1:365,
      y = my_elev_bands[-length(my_elev_bands)],
      z = plot_snow, col =my_col, breaks = my_bre,
      ylab = "", xlab = "", axes = F)
axis(1, at = x_axis_tic, c("","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S","O","N","D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.15, 0))#plot labels
mtext("Elevation [m]", side = 2, line = 1.5, cex = 0.8)
axis(2, mgp=c(3, 0.15, 0), tck = -0.001)
box()

par(mar = c(1.6, 0.5, 0.6, 2.7))

image_scale(as.matrix(plot_snow), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.15, 0), tck = -0.08)
mtext("Snow", side = 4, line = 1.5, cex = 0.8)

box()




#visu_meteo_comp----

#Plot melt compensation: Trend snow volume diff over elevation
melt_comp <- apply(vdis_band, 1, sum_na)

par(mfrow = c(1, 1))
par(mar = c(1.6, 3.0, 0.6, 0.6))

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15

plot(smoothFFT(melt_comp, sd = 5), type = "l", col = "black", axes = F,
     ylab = "", xlab = "", lwd = 2, ylim = rev(range(smoothFFT(melt_comp, sd = 5))))
abline(h = 0, lty = "dashed", lwd = 0.7)
abline(v = x_axis_tic, lty = "dashed", lwd = 0.5)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S","O","N","D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.15, 0))#plot labels
mtext("Trend total snow melt volume [m³/dec]", side = 2, line = 1.5, cex = 0.8)
axis(2, mgp=c(3, 0.15, 0), tck = -0.001)

box(lwd = 0.7)

#visu_meteo----

par(mfrow = c(1, 2))

plot_cycl_elev(data_in = tmea_band, data_mk = tmea_band, data_in_me = tmea_band_mea,
               data_meta = meta_grid_bands, main_text = paste0("Temperature [°C]"),
               margins_1 = c(1.4,1.8,1.8,0.2), margins_2 = c(1.4,0.2,1.8,3.5),
               no_col = F, show_mk = F, aggr_cat_mean = T, with_hom_dat = F,
               smooth_val = 0.2, mk_sig_level = 0.05, add_st_num = T)

# plot_cycl_elev(data_in = tmea_band_5km, data_mk = tmea_band_5km, data_in_me = tmea_band_mea_5km,
#                data_meta = meta_grid_bands, main_text = paste0("Temperature [°C]"),
#                margins_1 = c(1.4,1.8,1.8,0.2), margins_2 = c(1.4,0.2,1.8,3.5),
#                no_col = F, show_mk = F, aggr_cat_mean = T, with_hom_dat = F,
#                smooth_val = 0.2, mk_sig_level = 0.05, add_st_num = T)

plot_cycl_elev(data_in = tslo_band, data_mk = tslo_band, data_in_me = tslo_band_mea,
               data_meta = meta_grid_bands, main_text = paste0("Temperature [°C/dec]"),
               margins_1 = c(1.4,1.8,1.8,0.2), margins_2 = c(1.4,0.2,1.8,3.5),
               no_col = F, show_mk = F, aggr_cat_mean = T, with_hom_dat = F,
               smooth_val = 0.1, mk_sig_level = 0.05, add_st_num = T)

plot_cycl_elev(data_in = pmea_band, data_mk = pmea_band, data_in_me = pmea_band_mea,
               data_meta = meta_grid_bands, main_text= paste0("Precipitation [mm]"),
               margins_1 = c(1.4,1.8,1.8,0.2), margins_2 = c(1.4,0.2,1.8,3.5),
               no_col = F, show_mk = F, aggr_cat_mean = T, with_hom_dat = F,
               smooth_val = 0.2, mk_sig_level = 0.05, add_st_num = T)

plot_cycl_elev(data_in = pslo_band, data_mk = pslo_band, data_in_me = pslo_band_mea,
               data_meta = meta_grid_bands, main_text = paste0("Precipitation [mm/dec]"),
               margins_1 = c(1.4,1.8,1.8,0.2), margins_2 = c(1.4,0.2,1.8,3.5),
               no_col = F, show_mk = F, aggr_cat_mean = T, with_hom_dat = F,
               smooth_val = 0.2, mk_sig_level = 0.05, add_st_num = T)

plot_cycl_elev(data_in = emea_band, data_mk = emea_band, data_in_me = emea_band_mea,
               data_meta = meta_grid_bands, main_text= paste0("Evapotranspiration [mm]"),
               margins_1 = c(1.4,1.8,1.8,0.2), margins_2 = c(1.4,0.2,1.8,3.5),
               no_col = F, show_mk = F, aggr_cat_mean = T, with_hom_dat = F,
               smooth_val = 0.2, mk_sig_level = 0.05, add_st_num = T)

plot_cycl_elev(data_in = eslo_band, data_mk = eslo_band, data_in_me = eslo_band_mea,
               data_meta = meta_grid_bands, main_text = paste0("Evapotranspiration [mm/dec]"),
               margins_1 = c(1.4,1.8,1.8,0.2), margins_2 = c(1.4,0.2,1.8,3.5),
               no_col = F, show_mk = F, aggr_cat_mean = T, with_hom_dat = F,
               smooth_val = 0.2, mk_sig_level = 0.05, add_st_num = F)


#visu_spatial----

snows_d_mea <- apply(snows_d, 2, mean)
temps_d_mea <- apply(temps_d, 2, mean)
precs_d_mea <- apply(precs_d, 2, mean)

sno_vol_basin <- apply(snows_d, 1, sum_na) * 1000 * 1000
prec_basin <- apply(precs_d, 1, mea_na) * area_m2 / 1000 #[m³]

#corp DEM sub-basin area
my_ext <- extent(basin)

my_box <- as(my_ext, 'SpatialPolygons')
dem_cro_swiss <- raster::crop(dem, extent(my_box))
dem_sub_swiss <- mask(dem_cro, my_box)

val2col <- function(val_in){
  
  # dat_ref <- snows_d_mea
  # dat_ref <- temps_d_mea
  # dat_ref <- precs_d_mea
  dat_ref <- elevs_d
  
  # val_in <- log(val_in)
  dat_ref_log <- dat_ref
  
  # my_col <- colorRampPalette(c("white", viridis(9, direction = 1)[c(3,4)], "cadetblue3", "grey80",
  #                              "yellow2","gold", "orange2", "orangered2"))(200)
  my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
  # my_col <- colorRampPalette(c("grey95", viridis::viridis(9, direction = 1)[4:1]))(200)
      
  col_ind <- round((val_in-min_na(dat_ref_log)) / (max_na(dat_ref_log)-min_na(dat_ref_log)) * 200)  
    
  if(col_ind == 0){#for minimum and very small values
    
    col_ind <- 1
    
  }
  
  
  col_out <- my_col[col_ind]
  
  if(length(col_out) < 1){
    
     col_out <- "red"
    
  }
  
  return(col_out)
}

cols_spat <- foreach(i = 1:length(elevs_d), .combine = 'cbind') %dopar% {
  
  # val2col(temps_d_mea[i])
  # val2col(precs_d_mea[i])
  # val2col(snows_d_mea[i])
  val2col(elevs_d[i])
  
}

pdf(paste0(base_dir,"R/figs_exp/simu_spat.pdf"), width = 12, height = 6)

# plot(dem_sub_swiss, axes = F, legend = F,  col = colorRampPalette(c("white", "black"))(200), box = F)
plot(dem_sub, axes = F, legend = F,  col = colorRampPalette(c("white", "black"))(200), box = F)

plot(basin, add =T)
# points(grid_points_d_in@coords[, 1], grid_points_d_in@coords[, 2], pch = 19, col = cols_spat, cex = 0.30)
points(grid_points_d@coords[, 1], grid_points_d@coords[, 2], pch = 19, col = cols_spat, cex = 0.30)
points(grid_points_d@coords[, 1], grid_points_d@coords[, 2], pch = 19, cex = 0.20, col = "red3")

dev.off()

elevs_d[which(cols_spat == "red")]

range(dem_ele)
range(elevs_d)

hist(elevs_d, breaks = my_elev_bands)
hist(dem_ele, breaks = my_elev_bands)



plot(grid_points_d_in[1:10505], pch = 19)




swi_region <- st_bbox(c(xmin = 3950000, xmax = 4400000,
                        ymin = 2500000, ymax = 2850000))

swi_ele_map = tm_shape(dem, bbox = swi_region) +
  tm_raster(style = "cont", palette = "Greys", n = 20, legend.show = F) +
  tm_shape(basin) + tm_borders(col = "red3", alpha = 1, lwd = 2)

swi_ele_map


snows_d_mea <- apply(snows_d, 2, mean)

plot(grid_points_d_84[1:100])
plot(test[1:5])

hist(dem_ele, nclass = 50)






#melt_ext----

#Test spatial coherence results
points_sel <- which(grid_points_d@coords[, 1] > med_na(grid_points_d@coords[, 1]))
points_sel <- which(grid_points@coords[, 1] > min_na(grid_points@coords[, 1]))

range(elevs_d[points_sel])
my_elev_bands_sel <- seq(400, 3000, 50)

svolu_d_band_sel <- f_elev_bands(data_in = snows_d[, points_sel], 
                                 elev_bands = my_elev_bands_sel,
                                 func_aggr = "sum", 
                                 meta_dat = elevs_d[points_sel])*grid_m2 

#Bands calculate
yea_min_mag_bands_14 <- NULL
yea_min_doy_bands_14 <- NULL

for(i in 1:ncol(svolu_d_band_sel)){
  
  swe_band <- svolu_d_band_sel[, i]
  
  swe_band_dif <- c(NA, diff(swe_band))
  
  swe_band_dif[which(swe_band_dif > 0)] <- 0
  
  #Moving average filter
  swe_band_dif_ma_14 <- rollapply(data = swe_band_dif, width = 14,
                                  FUN = sum_na, align = "center", fill = NA)
  
  #Order data by day
  data_day_14 <- ord_day(data_in = swe_band_dif_ma_14,
                         date = date_snow,
                         start_y = 1954,
                         end_y = 2014,
                         break_day = 274,
                         do_ma = F,
                         window_width = 30)
  
  min_doy <- function(data_in){
    
    doy_min <- as.numeric(which(data_in == min_na(data_in)))
    
    if(length(doy_min) > 1){
      doy_min <- mea_na(doy_min)
    }
    
    return(doy_min)
  }
  
  yea_min_mag_14 <- apply(data_day_14, 1, min_na)
  yea_min_doy_14 <- apply(data_day_14, 1, min_doy)
  
  yea_min_mag_bands_14 <- cbind(yea_min_mag_bands_14, yea_min_mag_14)
  yea_min_doy_bands_14 <- cbind(yea_min_doy_bands_14, yea_min_doy_14)
  
}

freq_bands_1 <- NULL
freq_bands_2 <- NULL
freq_diff_all <- NULL
length_bands_all <- NULL
width_over <- 14

for(b in 1:ncol(svolu_d_band_sel)){
  
  band_sel <- b  
  
  bands_frequ_1 <- NULL
  bands_lengt_1 <- NULL
  for(y in 1:30){
    
    mid_day <- yea_min_doy_bands_14[y, band_sel]
    min_day <- mid_day - width_over
    max_day <- mid_day + width_over
    melt_window <- min_day:max_day
    
    bands_asso_1 <- which(yea_min_doy_bands_14[y, ] > min_day & yea_min_doy_bands_14[y, ] < max_day)
    
    bands_frequ_1 <- c(bands_frequ_1, bands_asso_1)
    bands_lengt_1 <- c(bands_lengt_1, length(bands_asso_1))
    
  }
  hist_1 <- hist(bands_frequ_1, breaks = seq(from = 0.5, to = ncol(svolu_d_band_sel)+0.5, by = 1), plot = F)
  
  bands_frequ_2 <- NULL
  bands_lengt_2 <- NULL
  for(x in 31:60){
    
    # print(y)
    mid_day <- yea_min_doy_bands_14[x, band_sel]
    min_day <- mid_day - width_over
    max_day <- mid_day + width_over
    melt_window <- min_day:max_day
    
    bands_asso_2 <- which(yea_min_doy_bands_14[x, ] > min_day & yea_min_doy_bands_14[x, ] < max_day)
    
    bands_frequ_2 <- c(bands_frequ_2, bands_asso_2)
    bands_lengt_2 <- c(bands_lengt_2, length(bands_asso_2))
    
  }
  hist_2 <- hist(bands_frequ_2, breaks = seq(from = 0.5, to = ncol(svolu_d_band_sel)+0.5, by = 1), plot = F)
  
  freq_diff <- hist_2$counts - hist_1$counts
  
  freq_bands_1 <- cbind(freq_bands_1, hist_1$counts)
  freq_bands_2 <- cbind(freq_bands_2, hist_2$counts)
  freq_diff_all <- cbind(freq_diff_all, freq_diff)
  length_bands_all <- cbind(length_bands_all, bands_lengt_1, bands_lengt_2)
  
}


pdf("/home/rottler/ownCloud/RhineFlow/rhine_snow/manus/meltim_v1/figures/bands_asso.pdf", width = 8, height = 5)

band_test <- 28 #my_elev_bands[band_test]
col_1 <- viridis(9, direction = 1)[4] 
col_2 <- "darkred"

par(mar = c(2.5, 2.5, 2, 0.2))

plot(freq_bands_1[, band_test], type = "h", col = alpha(col_1, alpha = 0.5), lwd = 7, ylab = "", xlab = "", 
     axes = F, lend = 2, xaxs = "i", yaxs = "i", ylim = c(0, 32), xlim = c(0, 60))
abline(v = band_test, col = "black", lwd = 1.5)
legend("topleft", c("1985-2014", "1954-1984"), col = c(col_2, col_1), pch = 19)
par(new = T)
plot(freq_bands_2[, band_test], type = "h", col = alpha(col_2, alpha = 0.5), lwd = 7, ylab = "", xlab = "", 
     axes = F, lend = 2, xaxs = "i", yaxs = "i", ylim = c(0, 32), xlim = c(0, 60))
elevs_sel <- c(6, 16, 26, 36, 46, 56)
axis(1, at = elevs_sel, labels = my_elev_bands[elevs_sel], mgp=c(3, 0.15, 0), tck = -0.01, cex.axis = 0.9)
axis(2, mgp=c(3, 0.15, 0), tck = -0.01, cex.axis = 0.9)
mtext("Elevation band [m]", side = 1, adj = 0.5, line = 1.3, cex = 1.1)
mtext("Frequency concurrent melt [-]", side = 2, adj = 0.5, line = 1.3, cex = 1.1)
mtext("a) Frequency concurrent melt (elevation band 1550-1600 m)", side = 3, adj = 0.0, line = 0.2, cex = 1.4)
box()

dev.off()


pdf("/home/rottler/ownCloud/RhineFlow/rhine_snow/manus/meltim_v1/figures/bands_chang.pdf", width = 8, height = 5)

par(mar = c(2.5, 2.5, 2.0, 0.2))

plot(freq_diff_all[band_test, ], type = "n", axes = F, ylab = "", xlab = "", xaxs = "i", yaxs = "i", 
     xlim = c(0, 60), ylim = c(-9, 9))
abline(v = band_test, col = "black", lwd = 1.5)
lines(freq_diff_all[band_test, ], lwd = 2, col = "grey55")
polygon(x = c(0.999, 1:59) , y = c(0, freq_diff_all[band_test, ]), col = alpha("grey55", alpha = 0.7), border = F)
elevs_sel <- c(6, 16, 26, 36, 46, 56)
axis(1, at = elevs_sel, labels = my_elev_bands[elevs_sel], mgp=c(3, 0.15, 0), tck = -0.01, cex.axis = 0.9)
axis(2, mgp=c(3, 0.15, 0), tck = -0.01, cex.axis = 0.9)
abline(h = 0, lty = "dotted")
mtext("Elevation band [m]", side = 1, adj = 0.5, line = 1.3, cex = 1.1)
mtext("Change frequency concurrent melt [-]", side = 2, adj = 0.5, line = 1.3, cex = 1.1)
mtext("b) Changes in concurrent melt (elevation band 1550-1600 m)", side = 3, adj = 0.0, line = 0.2, cex = 1.4)
box()

dev.off()


pdf("/home/rottler/ownCloud/RhineFlow/rhine_snow/manus/meltim_v1/figures/forgot_dim.pdf", width = 8, height = 5)

cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "cadetblue3", "grey90"))(50)
cols_max <- colorRampPalette(c("grey90", "gold3",  "orange3", "orangered4", "orangered4", "darkred"))(50)
cols_mel <- c(cols_min, cols_max)

layout(matrix(c(rep(1, 7), 2),
              1, 8), widths=c(), heights=c())

par(mar = c(3.5, 3.5, 2.5, 0.2))

image(x = 1:nrow(freq_diff_all),
      y = 1:ncol(freq_diff_all),
      z = freq_diff_all, 
      col = cols_mel, breaks = seq(from = -12, to = 12, length.out = 101), axes = F, ylab = "", xlab = "")
elevs_sel <- c(6, 16, 26, 36, 46, 56)-3
axis(1, at = elevs_sel, labels = my_elev_bands_sel[elevs_sel], mgp=c(3, 0.55, 0), tck = -0.005, cex.axis = 1.4)
axis(2, at = elevs_sel, labels = my_elev_bands_sel[elevs_sel], mgp=c(3, 0.25, 0), tck = -0.005, cex.axis = 1.4)
mtext("Elevation band [m]", side = 1, adj = 0.5, line = 2.1, cex = 1.1)
mtext("Elevation band [m]", side = 2, adj = 0.5, line = 2.1, cex = 1.1)
mtext("c) Changes in concurrent melt (all elevation bands)", side = 3, adj = 0.0, line = 0.2, cex = 1.4)
mtext("[-]", side = 3, adj = 1.0, line = 0.2, cex = 1.1)
box()

# points(c(18, 27, 35, 43, 51), c(18, 27, 35, 43, 51), pch = 21, cex = 2, col = "black")

par(mar = c(3.5, 0.2, 2.5, 4.0))

alptempr::image_scale(as.matrix(freq_diff_all), col = cols_mel, breaks = seq(from = -12, to = 12, length.out = 101), horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.35, 0), tck = -0.08, cex.axis = 1.4)
# mtext(lab_unit, side = 3, line = 0.3, cex = 1)
box()

dev.off()


#sc_frac----

sc <- snows_d
for(i in 1:nrow(sc)){
 
  print(i)
  row_sel <- i
  
  sc[row_sel, which(sc[row_sel, ] > 0.01)] <- 1
  sc[row_sel, which(sc[row_sel, ] < 0.01)] <- 0
  
}

sc_frac <- apply(sc, 1, sum_na)/ncol(sc)

plot(sc_frac[100:2000], type = "l")

save(sc_frac, file = "U:/rhine_snow/R/sc_frac.RData")

