###

#Rhine snow - Snow simulations + Meteo/snow analysis on sub-basin scale
#Erwin Rottler, University of Potsdam
#Spring 2019

###

sta_yea_bas <- 1961
end_yea_bas <- 2014
basin_sel <- "basel"        # alp_rhine,  reuss,     aare,  moselle, nahe,      neckar,   main,      lahn, basel
high_stat_thresh <- 1900 #1900
middle_stat_thresh <- 900 #900

snow_params <- read.table(paste0(base_dir, "R/melTim/snow_param.txt"), header = T, sep = ";")

# #ECHSE snow functions in package
# Rcpp::sourceCpp(paste0(base_dir, "R/meltim/echse_snow.cpp"))
# Rcpp::Rcpp.package.skeleton(name = "rEchseSnow", cpp_files = paste0(base_dir, "R/meltim/echse_snow.cpp"))
# install.packages(paste0(base_dir, "R/meltim/rEchseSnow"), repos=NULL, type="source")
# library(rEchseSnow)

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
dem = raster(paste0(base_dir, "data/basins_shp/eu_dem_500_fil.tif"))

#Load basin boundaries (shapefile delineated beforehand using Q-GIS)
my_basin <- paste0(base_dir, "data/basins_shp/", basin_sel, ".shp")
basin <- rgdal::readOGR(dsn = my_basin)

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
basin_84 <- sp::spTransform(basin, CRS = crswgs84)
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
basin_coords <- extract_coords(basin_84)

#get min/max of lat/lon for extraction from .nc file
lon_min <- min(basin_coords[ ,1]) - 0.15
lon_max <- max(basin_coords[ ,1]) + 0.15
lat_min <- min(basin_coords[ ,2]) - 0.15
lat_max <- max(basin_coords[ ,2]) + 0.15

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

print(paste0(Sys.time(), " Extract evapotranspiration from ncdf-file."))

petr_cube <- ncdf4::ncvar_get(nc_petr, start =c(lon_min_index, lat_min_index, date_min_index),
                              count = c(lon_count, lat_count, date_count), varid = "pet")

#get lat/lon values of extracted cube
lon2D <- lon[lon_min_index : (lon_min_index+lon_count-1), lat_min_index : (lat_min_index+lat_count-1)]
lat2D <- lat[lon_min_index : (lon_min_index+lon_count-1), lat_min_index : (lat_min_index+lat_count-1)]

#spatial grid points from lat/lon info
grid_points_cube_84 <-  sp::SpatialPoints(data.frame(lon = c(lon2D), lat = c(lat2D)), proj4string =  crswgs84)
grid_points_cube     <- sp::spTransform(grid_points_cube_84, CRS = crs(basin, asText = T))

#get grid points inside watershed
inside <- !is.na(sp::over(grid_points_cube, as(basin, "SpatialPolygons")))
grid_points <- grid_points_cube[which(inside == T)]
grid_points_84 <- sp::spTransform(grid_points, CRS = crs(grid_points_cube_84, asText = T))

#get index in cube from points inside sub-basin
cube_index_col <- sapply(grid_points_84@coords[,1], get_cube_index_col)
cube_index_row <- sapply(grid_points_84@coords[,1], get_cube_index_row)

#get time series from grid points inside sub-basin
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
#evapotranspiration
for (i in 1:length(cube_index_col)) {
  
  petr_sing <- petr_cube[cube_index_col[i], cube_index_row[i], ]
  
  if(i == 1){
    petrs <- petr_sing
  }else{
    petrs <- cbind(petrs, petr_sing)
  }
}

#sub-basin average values for rangking of weather types
temp_basin <- apply(temps, 1, med_na)
prec_basin <- apply(precs, 1, mea_na) 
petr_basin <- apply(petrs, 1, med_na) 

#Average mean snow water equivalent
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

radi_mea_smo <- smoothFFT(radi_mea, sd = 7)
# plot(radi_mea_smo, type = "l")

#Syntetic radiation time series based on mean values
radi_mea_seri <- rep( c(radi_mea_smo, radi_mea_smo[365], rep(radi_mea_smo, 3)), 50)[1:nrow(temps)]


#downscale----

#Downscale temperature/precipitation grid using simple lapse-rate based approach
#lapse rate in basin on daily bases used

#Lapse rate on daily basis from all data points selected
my_lin_trend <- function(data_in){
  
  lin_trend(data_in, index_in = elevs_ord)
  
}

temps_ord <- temps[, order(elevs)]
precs_ord <- precs[, order(elevs)]
petrs_ord <- petrs[, order(elevs)]

lapse_temp <- foreach(i = 1:nrow(temps_ord), .combine = 'c') %dopar% {
  
  my_lin_trend(temps_ord[i, ]) #[°C/m]
  
}
lapse_prec <- foreach(i = 1:nrow(precs_ord), .combine = 'c') %dopar% {
  
  my_lin_trend(precs_ord[i, ]) #[°C/m]
  
}
lapse_petr <- foreach(i = 1:nrow(petrs_ord), .combine = 'c') %dopar% {
  
  my_lin_trend(petrs_ord[i, ]) #[°C/m]
  
}
lapse_prec[which(is.na(lapse_prec))] <- 0 #when no rain at all recoreded NA; put to zero

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
  
  temps_down <- sapply(1:length(d_points_elevs_dif), f_laps_mod)
  return(temps_down)
}
f_petrs_laps <- function(data_in){
  
  f_laps_mod <- function(index_in){
    
    data_laps <- data_in + lapse_petr* d_points_elevs_dif[index_in]
    return(data_laps)
    
  }
  
  temps_down <- sapply(1:length(d_points_elevs_dif), f_laps_mod)
  return(temps_down)
}

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

petrs_d <- foreach(i = 1:ncol(petrs), .combine = 'cbind') %dopar% {
  
  d_points <- down_points(grid_points@coords[i, ])
  
  d_points_elevs <- raster::extract(dem, d_points)
  
  d_points_elevs_dif <- d_points_elevs - elevs[i]
  
  f_petrs_laps(petrs[, i]) #[mm]
  
}

elevs_d <- foreach(i = 1:length(grid_points), .combine = 'c') %dopar% {
  
  d_points <- down_points(grid_points@coords[i, ])
  
  raster::extract(dem, d_points)
  
}


#snow_simu----

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
        shortRad = radi_mea_seri[k],
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


#analysis----

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
            date = meteo_date, 
            start_year = sta_yea_bas, 
            end_year = end_yea_bas,
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

smea[, which(snows_d[nrow(snows_d), ] > 5)] <- NA #remove 'glacier' points

#Snow volume: mean average
vmea <- smea * grid_m2

#Snow height: trends 30DMA
for(b in 1:length(block_stas)){

  snows_calc <- snows_d[, block_stas[b]:block_ends[b]]

  print(paste(Sys.time(),"Trends 30DMA swe", "Block:", b, "out of", length(block_stas)))
  sslo_block <- foreach(i = 1:ncol(snows_calc), .combine = 'cbind') %dopar% {

    day_ana(snows_calc[, i], 
            date = meteo_date, 
            start_year = sta_yea_bas, 
            end_year = end_yea_bas,
            do_ma = T,
            window_width = 30, 
            method_ana = "sens_slope") #[m]

  }

  if(b == 1){
    sslo <- sslo_block
  }else{
    sslo <- cbind(sslo, sslo_block)
  }

}

sslo[, which(snows_d[nrow(snows_d), ] > 5)] <- NA #remove 'glacier' points

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
            date = meteo_date, 
            start_year = sta_yea_bas, 
            end_year = end_yea_bas,
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

sdif[, which(snows_d[nrow(snows_d), ] > 5)] <- NA #remove 'glacier' points

#Snow volume diff: mean average
vdif <- sdif * grid_m2

#Snow height diff: trends 30DMA
for(b in 1:length(block_stas)){

  snows_calc <- snows_d_dif[, block_stas[b]:block_ends[b]]

  print(paste(Sys.time(),"Trends 30 DMA diff swe", "Block:", b, "out of", length(block_stas)))
  sdis_block <- foreach(i = 1:ncol(snows_calc), .combine = 'cbind') %dopar% {

    day_ana(snows_calc[, i], 
            date = meteo_date, 
            start_year = sta_yea_bas, 
            end_year = end_yea_bas,
            do_ma = T,
            window_width = 30, 
            method_ana = "sens_slope") #[m]

  }

  if(b == 1){
    sdis <- sdis_block
  }else{
    sdis <- cbind(sdis, sdis_block)
  }

}

sdis[, which(snows_d[nrow(snows_d), ] > 5)] <- NA #remove 'glacier' points

#Snow volume diff: trends 30 DMA
vdis <- sdis * grid_m2

#Temperature: mean average
for(b in 1:length(block_stas)){
  
  temps_calc <- temps_d[, block_stas[b]:block_ends[b]]
  
  print(paste(Sys.time(),"Average median temperature", "Block:", b, "out of", length(block_stas)))
  
  tmea_block <- foreach(i = 1:ncol(temps_calc), .combine = 'cbind') %dopar% {
    
    day_ana(temps_calc[, i], 
            date = meteo_date, 
            start_year = sta_yea_bas, 
            end_year = end_yea_bas,
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

#Temperature: 30 DMA trends
for(b in 1:length(block_stas)){

  temps_calc <- temps_d[, block_stas[b]:block_ends[b]]

  print(paste(Sys.time(),"30DMA trends temperature", "Block:", b, "out of", length(block_stas)))

  tslo_block <- foreach(i = 1:ncol(temps_calc), .combine = 'cbind') %dopar% {

    day_ana(temps_calc[, i], 
            date = meteo_date, 
            start_year = sta_yea_bas, 
            end_year = end_yea_bas,
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
            date = meteo_date, 
            start_year = sta_yea_bas, 
            end_year = end_yea_bas,
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
            date = meteo_date, 
            start_year = sta_yea_bas, 
            end_year = end_yea_bas,
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

#Evaportranspiration: median average
for(b in 1:length(block_stas)){

  petrs_calc <- petrs_d[, block_stas[b]:block_ends[b]]

  print(paste(Sys.time(),"Average mean evapotranspiration", "Block:", b, "out of", length(block_stas)))

  emea_block <- foreach(i = 1:ncol(petrs_calc), .combine = 'cbind') %dopar% {

    day_ana(petrs_calc[, i], 
            date = meteo_date, 
            start_year = sta_yea_bas, 
            end_year = end_yea_bas,
            do_ma = F,
            window_width = 30, 
            method_ana = "mean") #[mm]

  }

  if(b == 1){
    emea <- emea_block
  }else{
    emea <- cbind(emea, emea_block)
  }

}

#Evapotranspiration: 30 DMA trends
for(b in 1:length(block_stas)){

  petrs_calc <- petrs_d[, block_stas[b]:block_ends[b]]

  print(paste(Sys.time(),"30DMA trends evapotranspiration", "Block:", b, "out of", length(block_stas)))

  eslo_block <- foreach(i = 1:ncol(petrs_calc), .combine = 'cbind') %dopar% {

    day_ana(petrs_calc[, i], 
            date = meteo_date, 
            start_year = sta_yea_bas, 
            end_year = end_yea_bas,
            do_ma = T,
            window_width = 30, 
            method_ana = "sens_slope") * 10 #[mm/dec]

  }

  if(b == 1){
    eslo <- eslo_block
  }else{
    eslo <- cbind(eslo, eslo_block)
  }

}


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
my_elev_bands <- c(seq(400, 3000, 50), 4000) #Diepoldsau

smea_band <- f_elev_bands(data_in = smea, func_aggr = "mean")
vmea_band <- f_elev_bands(data_in = vmea, func_aggr = "sum")
vslo_band <- f_elev_bands(data_in = vslo, func_aggr = "sum")
vdif_band <- f_elev_bands(data_in = vdif, func_aggr = "sum")
vdis_band <- f_elev_bands(data_in = vdis, func_aggr = "sum")

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

emea_band <- f_elev_bands(data_in = emea, func_aggr = "mean")
colnames(emea_band) <- paste0("band", 1:ncol(emea_band))
emea_band_mea <- apply(emea_band, 2, mea_na)#annual average values

eslo_band <- f_elev_bands(data_in = eslo, func_aggr = "mean")
colnames(eslo_band) <- paste0("band", 1:ncol(eslo_band))
eslo_band_mea <- apply(eslo_band, 2, mea_na)#annual average values



#visu----

#Plot: Snow variables

plot_snow <- vdis_band #smea_band, vmea_band, vslo_band, vdif_band, vdis_band
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

#Plot meteo

par(mfrow = c(1, 2))

plot_cycl_elev(data_in = tmea_band, data_mk = tmea_band, data_in_me = tmea_band_mea,
               data_meta = meta_grid_bands, main_text = paste0("Temperature [°C]"),
               margins_1 = c(1.4,1.8,1.8,0.2), margins_2 = c(1.4,0.2,1.8,3.5),
               no_col = F, show_mk = F, aggr_cat_mean = T, with_hom_dat = F,
               smooth_val = 0.2, mk_sig_level = 0.05, add_st_num = T)

plot_cycl_elev(data_in = tslo_band, data_mk = tslo_band, data_in_me = tslo_band_mea,
               data_meta = meta_grid_bands, main_text = paste0("Temperature [°C/dec]"),
               margins_1 = c(1.4,1.8,1.8,0.2), margins_2 = c(1.4,0.2,1.8,3.5),
               no_col = F, show_mk = F, aggr_cat_mean = T, with_hom_dat = F,
               smooth_val = 0.2, mk_sig_level = 0.05, add_st_num = T)

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










