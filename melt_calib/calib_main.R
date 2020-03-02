###

#Snow model calibration using using PSO / DDS
#Erwin Rottler, Summer 2019

###

#set up----

# devtools::install_github("TillF/ppso")
pacman::p_load(ppso, rEchseSnow, alptempr, rgeos, raster, meltimr, rfs, 
               parallel, doParallel, lhs, hydroGOF)

#set base direcoty
base_dir <- "U:/rhine_snow/"
# base_dir <- "/home/erwin/Documents/"

#gridded climate data
file_dir <- "d:/nrc_user/rottler/toErwin1/6435060/"

#Wrapper for snow simulation and calculation of bjective function
source("melt_calib/optim_wrapper.R")

# #Read snow station data
# snow_data_all_1 <- read.table(paste0(base_dir,"data/snow/order_73106_data.txt"), sep = ";", skip = 2, stringsAsFactors = F, na.strings = "-")
# snow_data_all_2 <- read.table(paste0(base_dir,"data/snow/order_73164_data.txt"), sep = ";", skip = 2, stringsAsFactors = F, na.strings = "-")
# snow_data_all_3 <- read.table(paste0(base_dir,"data/snow/order_74297_data.txt"), sep = ";", skip = 2, stringsAsFactors = F, na.strings = "-")
# snow_data_all <- rbind(snow_data_all_1, snow_data_all_2, snow_data_all_3)
# 
# #Read station to calibrate
# stat_calib <- read.table(paste0(base_dir, "R/meltim/melt_calib/station_calib.txt"), sep = ";", header = T)

#Read initial snow parameters
snow_params <- read.table(paste0(base_dir, "R/meltim/snow_param.txt"), header = T, sep = ";")
snow_params <- read.table(paste0(base_dir, "R/meltim/snow_param.txt"), header = T, sep = ";")

sta_yea_sno <- 1984 #start year snow simulation
end_yea_sno <- 2015 #end year snow simulation

sta_yea_cal <- 1985 #start year calibration period (compare simulated end measured values)
end_yea_cal <- 2014 #end year calibration period

#Cluster for parallel computing

# stopCluster(my_clust)

n_cores <- 35 #number of cores used for parallel computing

#Make cluster for parallel computing
# my_clust <- makeCluster(n_cores)
# clusterEvalQ(my_clust, pacman::p_load(zoo, zyp, alptempr, lmomco, ncdf4, rEchseSnow, sp, raster, betareg, rfs, meltimr))
# registerDoParallel(my_clust)
my_clust <- makeCluster(n_cores)
clusterEvalQ(my_clust, pacman::p_load(alptempr, rEchseSnow, meltimr))
registerDoParallel(my_clust)

#data_prep_grid----

#Load ncdf E-OBS gridded datasets
nc_temp_file <- paste0(file_dir, "meteo/tavg.nc")
nc_prec_file <- paste0(file_dir, "meteo/pre.nc")

nc_temp <- ncdf4::nc_open(nc_temp_file)
nc_prec <- ncdf4::nc_open(nc_prec_file)

#Projections used later
crswgs84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
epsg3035 <- sp::CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000
                    +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

#get lat/lon/time of .nc meteo data
lon <- ncdf4::ncvar_get(nc_temp, varid = "lon2D")
lat <- ncdf4::ncvar_get(nc_temp, varid = "lat2D")
date <- as.Date(as.character(ncdf4.helpers::nc.get.time.series(nc_temp, time.dim.name = "time")))

#define date start index and count for extraction from .nc file
date_min_index <- which(date == paste0(sta_yea_sno, "-01-01"))
date_max_index <- which(date == paste0(end_yea_sno, "-12-31"))
date_count <- date_max_index - date_min_index +1
meteo_date <- date[date_min_index:date_max_index]

#spatial grid points from lat/lon info
grid_points_84 <-  sp::SpatialPoints(data.frame(lon = c(lon), lat = c(lat)), proj4string =  crswgs84)
grid_points    <- sp::spTransform(grid_points_84, CRS = epsg3035)

stat_points_84 <- sp::SpatialPoints(data.frame(lon = stat_calib$longitude, lat = stat_calib$latitude), proj4string =  crswgs84)
stat_points    <- sp::spTransform(stat_points_84, CRS = epsg3035)

#Extract meteo time series from grid
for(k in 1:nrow(stat_calib)){
  
  print(paste(Sys.time(), "Extract meteo data of station", k, "out of", nrow(stat_calib)))
  stat_point <- stat_points[k]
  
  dist_point <- rep(NA, length(grid_points))
  
  for(i in 1:length(grid_points)){
    # print(i)
    dist_point[i] <- gDistance(spgeom1 = grid_points[i], stat_point) 
  }
  
  point_index <- which(dist_point == min_na(dist_point))
  
  point_sel <- grid_points_84@coords[which(dist_point == min_na(dist_point)), ]
  
  get_row <- function(row_num){
    
    point_sel[2] %in% lat[row_num,]
    
  }
  
  get_col <- function(col_num){
    
    point_sel[2] %in% lat[, col_num]
    
  }
  
  row_sel <- which(sapply(1:nrow(lat), get_row) == T)
  col_sel <- which(sapply(1:ncol(lat), get_col) == T)
  
  
  #Extract meteo data for determined grid cell(s)
  
  temp <- ncdf4::ncvar_get(nc_temp, start =c(row_sel, col_sel, date_min_index), 
                           count = c(1, 1, date_count), varid = "tavg")
  prec <- ncdf4::ncvar_get(nc_prec, start =c(row_sel, col_sel, date_min_index), 
                           count = c(1, 1, date_count), varid = "pre")
  
  if(k ==1){
    
    temps_all <- temp
    precs_all <- prec
    
  }else{
    
    temps_all <- cbind(temps_all, temp)
    precs_all <- cbind(precs_all, prec)
    
  }
  
}

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
radi_mea_seri <- rep( c(radi_mea_smo, radi_mea_smo[365], rep(radi_mea_smo, 3)), 50)[1:nrow(temps_all)]


#Get snow measurements data

for(i in 1:nrow(stat_calib)){
  
  station_id <- stat_calib$ID[i]
  sel_ind <- which(snow_data_all$V1 == as.character(station_id))
  snow_valu <- as.numeric(snow_data_all$V3[sel_ind])
  snow_date <- as.Date(strptime(snow_data_all$V2[sel_ind], "%Y%m%d", tz = "UTC"))
  
  data_full <- data.frame(date = snow_date, value = snow_valu)
  
  input_data <- data_full[as.numeric(format(data_full$date, "%Y")) >= sta_yea_sno, ]
  
  input_data <- input_data[as.numeric(format(input_data$date,  "%Y")) <= end_yea_sno, ]
  
  start_date <- as.POSIXct(strptime(paste0(sta_yea_sno, "-01-01"), "%Y-%m-%d", tz = "UTC"))
  end_date <- as.POSIXct(strptime(paste0(end_yea_sno, "-12-31"), "%Y-%m-%d", tz = "UTC"))
  full_date <- seq(start_date, end_date, by = "day")
  
  input_data <- data.frame(dates = full_date, values = with(input_data, 
                                                            value[match(as.Date(full_date), as.Date(date))]))
  
  if(i == 1){
    
    snows_stat_all <- input_data$value
    
  }else{
    
    snows_stat_all <- cbind(snows_stat_all, input_data$values)
    
  }
  
}



#data_prep_stat----

#Read station data
rain_data_all_1 <- read.table(paste0(base_dir,"data/idaweb/order74624/order_74624_data.txt"), sep = ";", 
                              skip = 2, stringsAsFactors = F, na.strings = "-")
rain_data_all_2 <- read.table(paste0(base_dir,"data/idaweb/order74631/order_74631_data.txt"), sep = ";", 
                              skip = 2, stringsAsFactors = F, na.strings = "-")
rain_data_all <- rbind(rain_data_all_1, rain_data_all_2)

temp_data_all_1 <- read.table(paste0(base_dir,"data/idaweb/order74625/order_74625_data.txt"), sep = ";", 
                              skip = 2, stringsAsFactors = F, na.strings = "-")
temp_data_all_2 <- read.table(paste0(base_dir,"data/idaweb/order74632/order_74632_data.txt"), sep = ";", 
                              skip = 2, stringsAsFactors = F, na.strings = "-")
temp_data_all <- rbind(temp_data_all_1, temp_data_all_2)

#Get measurements data

for(i in 1:nrow(stat_calib)){
  
  print(i)
  
  station_id <- stat_calib$ID[i]
  sel_ind_sno <- which(snow_data_all$V1 == as.character(station_id))
  sel_ind_rai <- which(rain_data_all$V1 == as.character(station_id))
  sel_ind_tem <- which(temp_data_all$V1 == as.character(station_id))
  
  snow_valu <- as.numeric(snow_data_all$V3[sel_ind_sno])
  rain_valu <- as.numeric(rain_data_all$V3[sel_ind_rai])
  temp_valu <- as.numeric(temp_data_all$V3[sel_ind_tem])
  
  snow_date <- as.Date(strptime(snow_data_all$V2[sel_ind_sno], "%Y%m%d", tz = "UTC"))
  rain_date <- as.Date(strptime(rain_data_all$V2[sel_ind_rai], "%Y%m%d", tz = "UTC"))
  temp_date <- as.Date(strptime(temp_data_all$V2[sel_ind_tem], "%Y%m%d", tz = "UTC"))
  
  data_full_sno <- data.frame(date = snow_date, value = snow_valu)
  data_full_rai <- data.frame(date = rain_date, value = rain_valu)
  data_full_tem <- data.frame(date = temp_date, value = temp_valu)
  
  input_data_sno <- data_full_sno[as.numeric(format(data_full_sno$date, "%Y")) >= sta_yea_sno, ]
  input_data_rai <- data_full_rai[as.numeric(format(data_full_rai$date, "%Y")) >= sta_yea_sno, ]
  input_data_tem <- data_full_tem[as.numeric(format(data_full_tem$date, "%Y")) >= sta_yea_sno, ]
  
  input_data_sno <- input_data_sno[as.numeric(format(input_data_sno$date,  "%Y")) <= end_yea_sno, ]
  input_data_rai <- input_data_rai[as.numeric(format(input_data_rai$date,  "%Y")) <= end_yea_sno, ]
  input_data_tem <- input_data_tem[as.numeric(format(input_data_tem$date,  "%Y")) <= end_yea_sno, ]
  
  start_date <- as.POSIXct(strptime(paste0(sta_yea_sno, "-01-01"), "%Y-%m-%d", tz = "UTC"))
  end_date <- as.POSIXct(strptime(paste0(end_yea_sno, "-12-31"), "%Y-%m-%d", tz = "UTC"))
  full_date <- seq(start_date, end_date, by = "day")
  
  input_data_sno <- data.frame(dates = full_date, values = with(input_data_sno, 
                                                            value[match(as.Date(full_date), as.Date(date))]))
  input_data_rai <- data.frame(dates = full_date, values = with(input_data_rai, 
                                                                value[match(as.Date(full_date), as.Date(date))]))
  input_data_tem <- data.frame(dates = full_date, values = with(input_data_tem, 
                                                                value[match(as.Date(full_date), as.Date(date))]))
  
  if(i == 1){
    
    snows_stat_all <- input_data_sno$value
    precs_all <- input_data_rai$value
    temps_all <- input_data_tem$value
    
  }else{
    
    snows_stat_all <- cbind(snows_stat_all, input_data_sno$values)
    precs_all <- cbind(precs_all, input_data_rai$values)
    temps_all <- cbind(temps_all, input_data_tem$values)
    
  }
  
}

cal_sta_ind <- min_na(which(format(full_date, "%Y") == sta_yea_cal))
cal_end_ind <- min_na(which(format(full_date, "%Y") == end_yea_cal))

for(i in 1:nrow(stat_calib)){
  print(i)
  print(as.character(stat_calib$name[i]))
  print(summary(precs_all[cal_sta_ind:cal_end_ind, i]))
  print(summary(temps_all[cal_sta_ind:cal_end_ind, i]))
  print(summary(snows_stat_all[cal_sta_ind:cal_end_ind, i]))
  
}

# #Remove stations with NAs in temperature or precipitation
# stat_rem <- c(5, 13, 17)
# precs_all <- precs_all[, -stat_rem]
# temps_all <- temps_all[, -stat_rem]
# snows_stat_all <- snows_stat_all[, -stat_rem]

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

radi_mea_smo <- smoothFFT(radi_mea, sd = 19)
# plot(radi_mea_smo, type = "l")

#Syntetic radiation time series based on mean values
radi_mea_seri <- rep( c(radi_mea_smo, radi_mea_smo[365], rep(radi_mea_smo, 3)), 50)[1:nrow(temps_all)]

meteo_date <- as.Date(full_date) #when preparing grid input, date called meteo_date


#data_prep_swe----

f_read_swe <- function(slf_swe_file, sta_year = sta_yea_sno, end_year = end_yea_sno){
  
  swe_data <- read.table(slf_swe_file, header = T, sep = ",")
  swe_data$DATUM <- as.Date(swe_data$DATUM, "%d.%m.%Y")
  
  start_date <- as.POSIXct(strptime(paste0(1900, "-01-01"), "%Y-%m-%d", tz = "UTC"))
  end_date <- as.POSIXct(strptime(paste0(2020, "-12-31"), "%Y-%m-%d", tz = "UTC"))
  full_date <- seq(start_date, end_date, by = "day")
  
  swe_data_full <- data.frame(dates = full_date, 
                              values = with(swe_data, SWE..mm.[match(as.Date(full_date), as.Date(DATUM))]))
  
  swe_data_full <- swe_data_full[format(swe_data_full$dates, '%Y') >= sta_year, ]
  swe_data_full <- swe_data_full[format(swe_data_full$dates, '%Y') <= end_year, ]
  
  return(swe_data_full)
  
}

swe_and <- f_read_swe(slf_swe_file = paste0(base_dir, "data/slf_swe/sp_2AN.txt"))
swe_dav <- f_read_swe(slf_swe_file = paste0(base_dir, "data/slf_swe/sp_5DF.txt"))      
swe_wfj <- f_read_swe(slf_swe_file = paste0(base_dir, "data/slf_swe/sp_5WJ.txt")) 
# swe_sta <- f_read_swe(slf_swe_file = paste0(base_dir, "data/slf_swe/sp_7ST.txt"))
# swe_zer <- f_read_swe(slf_swe_file = paste0(base_dir, "data/slf_swe/sp_4ZE.txt"))
# swe_ulr <- f_read_swe(slf_swe_file = paste0(base_dir, "data/slf_swe/sp_4UL.txt"))  

swe_all <- cbind(swe_and$values, swe_dav$values, swe_wfj$values)

swe_all <- swe_all / 1000 # [m]

#Read temperature and precipitation data

meteo_data_1 <- read.table(paste0(base_dir,"data/idaweb/order78285/order_78285_data.txt"), sep = ";", nrows = 112947, 
                        skip = 2, stringsAsFactors = F, na.strings = "-", blank.lines.skip = T, header = T)

meteo_data_2 <- read.table(paste0(base_dir,"data/idaweb/order78285/order_78285_data.txt"), sep = ";", nrows = 49060,
                        skip = 112952, stringsAsFactors = F, na.strings = "-", blank.lines.skip = T, header = T)

meteo_data_3 <- read.table(paste0(base_dir,"data/idaweb/order78285/order_78285_data.txt"), sep = ";", nrows = 148505,
                        skip = (112952+49062), stringsAsFactors = F, na.strings = "-", blank.lines.skip = T, header = T)

meteo_data_4 <- read.table(paste0(base_dir,"data/idaweb/order78285/order_78285_data.txt"), sep = ";", nrows = 21190,
                          skip = (112952+49062+148513), stringsAsFactors = F, na.strings = "-", blank.lines.skip = T, header = T)

meteo_data_5 <- read.table(paste0(base_dir,"data/idaweb/order78285/order_78285_data.txt"), sep = ";",
                          skip = (112952+49062+148513+21200), stringsAsFactors = F, na.strings = "-", blank.lines.skip = T, header = T)

colnames(meteo_data_5)
unique(meteo_data_5$stn)[which(unique(meteo_data_5$stn) != "stn")]

f_meteo_ext <- function(data_in, ID_in, var_in, sta_year = sta_yea_sno, end_year = end_yea_sno){

  row_sel <- which(data_in$stn == ID_in)
  data_sel <- data_in[row_sel, ]
  
  data_sel$time <- as.Date(data_sel$time, "%Y%m%d")
  
  start_date <- as.POSIXct(strptime(paste0(1900, "-01-01"), "%Y-%m-%d"))
  end_date <- as.POSIXct(strptime(paste0(2020, "-12-31"), "%Y-%m-%d"))
  full_date <- seq(start_date, end_date, by = "day")
  
  if(var_in == "tre200d0"){
    
    data_sel$tre200d0 <- as.numeric(data_sel$tre200d0)
    
    data_sel_full <- data.frame(dates = full_date, 
                                values = with(data_sel, tre200d0[match(as.Date(full_date), as.Date(time))]))
    
  }
  
  if(var_in == "rre150d0"){
    
    data_sel$rre150d0 <- as.numeric(data_sel$rre150d0)
    
    data_sel_full <- data.frame(dates = full_date, 
                                values = with(data_sel, rre150d0[match(as.Date(full_date), as.Date(time))]))
    
  }
  
  data_sel_full$values <- as.numeric(data_sel_full$values)
  data_sel_full <- data_sel_full[format(data_sel_full$dates, '%Y') >= sta_year, ]
  data_sel_full <- data_sel_full[format(data_sel_full$dates, '%Y') <= end_year, ]
  
  return(data_sel_full)
  
}

temp_ant <- f_meteo_ext(data_in = meteo_data_1,
                        ID_in = "ANT",
                        var_in = "tre200d0")

rain_ant <- f_meteo_ext(data_in = meteo_data_1,
                        ID_in = "ANT",
                        var_in = "rre150d0")

temp_dav <- f_meteo_ext(data_in = meteo_data_1,
                        ID_in = "DAV",
                        var_in = "tre200d0")

rain_dav <- f_meteo_ext(data_in = meteo_data_1,
                        ID_in = "DAV",
                        var_in = "rre150d0")

temp_wfj <- f_meteo_ext(data_in = meteo_data_3,
                        ID_in = "WFJ",
                        var_in = "tre200d0")

rain_wfj <- f_meteo_ext(data_in = meteo_data_3,
                        ID_in = "WFJ",
                        var_in = "rre150d0")

# temp_zer <- f_meteo_ext(data_in = meteo_data_3,
#                         ID_in = "ZER",
#                         var_in = "tre200d0")
# 
# rain_zer <- f_meteo_ext(data_in = meteo_data_3,
#                         ID_in = "ZER",
#                         var_in = "rre150d0")

# temp_ulr <- f_meteo_ext(data_in = meteo_data_3,
#                         ID_in = "ULR",
#                         var_in = "tre200d0")
# 
# rain_ulr <- f_meteo_ext(data_in = meteo_data_3,
#                         ID_in = "ULR",
#                         var_in = "rre150d0")
# 
# temp_smm <- f_meteo_ext(data_in = meteo_data_3,
#                         ID_in = "SMM",
#                         var_in = "tre200d0")
# 
# rain_smm <- f_meteo_ext(data_in = meteo_data_3,
#                         ID_in = "SMM",
#                         var_in = "rre150d0")

temps_all <- cbind(temp_ant$values, temp_dav$values, temp_wfj$values)
precs_all <- cbind(rain_ant$values, rain_dav$values, rain_wfj$values)


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

radi_mea_smo <- smoothFFT(radi_mea, sd = 19)
# plot(radi_mea_smo, type = "l")

#Syntetic radiation time series based on mean values
radi_mea_seri <- rep( c(radi_mea_smo, radi_mea_smo[365], rep(radi_mea_smo, 3)), 50)[1:nrow(temps_all)]

start_date <- as.POSIXct(strptime(paste0(sta_yea_sno, "-01-01"), "%Y-%m-%d"))
end_date <- as.POSIXct(strptime(paste0(end_yea_sno, "-12-31"), "%Y-%m-%d"))
meteo_date <- seq(start_date, end_date, by = "day")

#calib_sing----

for(i in 1:ncol(temps_all)){
  
  # temps <- temps_all[, i]
  # precs <- precs_all[, i]
  # snows_stat <- snows_stat_all[, i]
  temps <- temps_all
  precs <- precs_all
  snows_stat <- swe_all
  
  save(temps, precs, radi_mea_seri, snows_stat, snow_params, meteo_date, sta_yea_cal, end_yea_cal,
       file = paste0(base_dir, "R/meltim/melt_calib/calib_thread/calib_snow_data.Rdata"), version = 2)
  
  param_ranges=rbind(#define parameter ranges
    
    # a0                  = c(0, 0.01), #empirical coefficient (m/s); linear dependence of turbulent transfer coefficient (D) in sensible heat flux: D = a0 + a1 * WindSpeed
    # a1                  = c(0, 0.01),	#empirical coefficient; linear dependence of turbulent transfer coefficient (D) in sensible heat flux: D = a0 + a1 * WindSpeed
    kSatSnow            = c(0.000001, 0.001), #Saturated hydraulic conductivity of snow [m/s]
    densDrySnow         = c(200,600), #Snow dry density [kg/m³]
    specCapRet          = c(0.01, 0.8), #Capill. retention volume as fraction of solid SWE [-]
    emissivitySnowMin   = c(0.30, 0.90), #Minimum emissivity old snow [-]
    emissivitySnowMax   = c(0.90, 0.99), #Maximum emissivity of fresh snow [-]
    tempAir_crit        = c(-3, 3), #Threshold temperature for rain-/snowfall [°C]
    albedoMin           = c(0.3, 0.7), #Minimum snow albedo after a long time without snowfall
    albedoMax           = c(0.8, 0.99), #Maximum snow albedo right after snowfall
    agingRate_tAirPos   = c(0.0000001,  0.0001), #Rate to describe the intensity of snow aging process at positive temperatures
    agingRate_tAirNeg   = c(0.00000001, 0.0001), #Rate to describe the intensity of snow aging process at negative temperatures
    # soilDepth           = c(0.05, 0.2), #Depth upper soil (thermally interacts with snow cover) [m]
    # soilDens            = c(1000, 1500), #Density upper soil [kg/m³]
    # soilSpecHeat        = c(1.8, 2.9 #Specific heat capacity upper soil [kJ/kg/K]
    weightAirTemp       = c(0, 1) #Weighting paramameter for air temp.  and snow mean temp. to get surface temp snow (-) 0...1
    # precipSeconds       = c(86400, 86400) #Seconds per time step
    # tempMaxOff          = c(0, 2), #Offset maximum hourly temperatures to 12:00
    # tempAmpli           = c(1, 15) #Amplitude daily temperature cycle
    
  )
  
  
  starting.values=rbind(#initial estimates
    
    a0                  = 0.002,
    a1                  = 0.0008,
    kSatSnow            = 0.00004,
    densDrySnow         = 450,
    specCapRet          = 0.05, 
    emissivitySnowMin   = 0.84,
    emissivitySnowMax   = 0.99,
    tempAir_crit        = 0.2,
    albedoMin           = 0.55,
    albedoMax           = 0.88,
    agingRate_tAirPos   = 0.00000111,
    agingRate_tAirNeg   = 0.000000462,
    # soilDepth           = 0.1, 
    # soilDens            = 1300,
    # soilSpecHeat        = 2.18,
    weightAirTemp       = 0.5
    # precipSeconds       = 86400, #Seconds per time step
    # tempMaxOff          = 2, 
    # tempAmpli           = 8
    
  )
  
  max_number_function_calls = 1000
  
  # source("2_optim_wrapper.R")
  
  res <- optim_dds(objective_function=optim_wrapper, number_of_particles=3, number_of_parameters=NROW(param_ranges), parameter_bounds=param_ranges, initial_estimates=starting.values, lhc_init=T,
                   logfile="melt_calib/dds.log", projectfile="melt_calib/dds.pro", load_projectfile="no", break_file=NULL, max_wait_iterations=1000, max_number_function_calls=max_number_function_calls, tryCall=F)
  
  # plot_optimization_progress(logfile = "dds.log", projectfile = "dds.pro")
  # 
  # savePlot(filename = paste0("dds_progress1_", i, ".png"), type = "png", device = dev.prev())
  # savePlot(filename = paste0("dds_progress2_", i, ".png"), type = "png", device = dev.cur())
  # graphics.off()
  
  file.rename("melt_calib/dds.log", paste0("melt_calib/dds", i, ".log"))
  
  }

#calib_all----

# base_dir <- "/home/erwin/Documents/"

temps <- temps_all
precs <- precs_all
# snows_stat <- snows_stat_all
snows_stat <- swe_all

# #snow parameters start
# snow_params=cbind(
#   
#   a0                  = 0.002,
#   a1                  = 0.0008,
#   kSatSnow            = 0.00004,
#   densDrySnow         = 450,
#   specCapRet          = 0.05, 
#   emissivitySnowMin   = 0.84,
#   emissivitySnowMax   = 0.99,
#   tempAir_crit        = 0.2,
#   albedoMin           = 0.55,
#   albedoMax           = 0.88,
#   agingRate_tAirPos   = 0.00000111,
#   agingRate_tAirNeg   = 0.000000462,
#   soilDepth           = 0.1,
#   soilDens            = 1300,
#   soilSpecHeat        = 2.18,
#   weightAirTemp       = 0.5,
#   precipSeconds       = 86400, #Seconds per time step
#   tempMaxOff          = 2,
#   tempAmpli           = 8
#   
# )
# 
# snow_params <- as.data.frame(snow_params)
  
save(temps, precs, radi_mea_seri, snows_stat, snow_params, meteo_date, sta_yea_cal, end_yea_cal,
     file = paste0(base_dir, "R/meltim/melt_calib/calib_thread/calib_snow_data.Rdata"), version = 2)

param_ranges=rbind(#define parameter ranges
  
  # a0                  = c(0, 0.01), #empirical coefficient (m/s); linear dependence of turbulent transfer coefficient (D) in sensible heat flux: D = a0 + a1 * WindSpeed
  # a1                  = c(0, 0.01),	#empirical coefficient; linear dependence of turbulent transfer coefficient (D) in sensible heat flux: D = a0 + a1 * WindSpeed
  kSatSnow            = c(0.000001, 0.001), #Saturated hydraulic conductivity of snow [m/s]
  densDrySnow         = c(200,600), #Snow dry density [kg/m³]
  specCapRet          = c(0.01, 0.8), #Capill. retention volume as fraction of solid SWE [-]
  emissivitySnowMin   = c(0.30, 0.90), #Minimum emissivity old snow [-]
  emissivitySnowMax   = c(0.90, 0.99), #Maximum emissivity of fresh snow [-]
  tempAir_crit        = c(-3, 3), #Threshold temperature for rain-/snowfall [°C]
  albedoMin           = c(0.3, 0.7), #Minimum snow albedo after a long time without snowfall
  albedoMax           = c(0.8, 0.99), #Maximum snow albedo right after snowfall
  agingRate_tAirPos   = c(0.0000001,  0.0001), #Rate to describe the intensity of snow aging process at positive temperatures
  agingRate_tAirNeg   = c(0.00000001, 0.0001), #Rate to describe the intensity of snow aging process at negative temperatures
  # soilDepth           = c(0.05, 0.2), #Depth upper soil (thermally interacts with snow cover) [m]
  # soilDens            = c(1000, 1500), #Density upper soil [kg/m³]
  # soilSpecHeat        = c(1.8, 2.9 #Specific heat capacity upper soil [kJ/kg/K]
  weightAirTemp       = c(0, 1) #Weighting paramameter for air temp.  and snow mean temp. to get surface temp snow (-) 0...1
  # precipSeconds       = c(86400, 86400) #Seconds per time step
  # tempMaxOff          = c(0, 2), #Offset maximum hourly temperatures to 12:00
  # tempAmpli           = c(1, 15) #Amplitude daily temperature cycle
  
)


starting.values=rbind(#initial estimates
  
  # a0                  = 0.002,
  # a1                  = 0.0008,
  kSatSnow            = 0.00004,
  densDrySnow         = 450,
  specCapRet          = 0.05, 
  emissivitySnowMin   = 0.84,
  emissivitySnowMax   = 0.99,
  tempAir_crit        = 0.2,
  albedoMin           = 0.55,
  albedoMax           = 0.88,
  agingRate_tAirPos   = 0.00000111,
  agingRate_tAirNeg   = 0.000000462,
  # soilDepth           = 0.1, 
  # soilDens            = 1300,
  # soilSpecHeat        = 2.18,
  weightAirTemp       = 0.5
  # precipSeconds       = 86400, #Seconds per time step
  # tempMaxOff          = 2, 
  # tempAmpli           = 8
  
)

max_number_function_calls = 5000

#Save results indivudial stations
results_individual <- NULL
save(results_individual, file =  paste0(base_dir, "R/meltim/melt_calib/results_stations.Rdata"), version = 2)

res <- optim_dds(objective_function=optim_wrapper, number_of_particles=3, number_of_parameters=NROW(param_ranges), parameter_bounds=param_ranges, initial_estimates=starting.values, lhc_init=T,
                 logfile="melt_calib/dds.log", projectfile="melt_calib/dds.pro", load_projectfile="no", break_file=NULL, max_wait_iterations=1000, max_number_function_calls=max_number_function_calls, tryCall=F)

plot_optimization_progress(logfile = "melt_calib/dds.log", projectfile = "melt_calib/dds.pro")

savePlot(filename = paste0("melt_calib/dds_progress1", ".png"), type = "png", device = dev.prev())
savePlot(filename = paste0("melt_calib/dds_progress2", ".png"), type = "png", device = dev.cur())
graphics.off()


#visu_sing----

for(i in 1:ncol(temps_all)){
  
  dds_log <- read.table(paste0(base_dir, "R/meltim/melt_calib/dds", i, ".log"), header = T, sep = "\t")
  
  best_row <- min_na(which(dds_log$objective_function == min_na(dds_log$objective_function)))
  
  best_run <- dds_log[best_row, ]
  
  if(i == 1){
    
    best_all <- best_run
    
  }else{
    best_all <- rbind(best_all, best_run)
  }
  
  
}

best_all

for(i in 2:(ncol(best_all)-1)){
  plot(best_all[, i], pch = 19, main = colnames(best_all)[i])
}


#visu_all----

dds_log <- read.table(paste0(base_dir, "meltim/melt_calib/dds", ".log"), header = T, sep = "\t")

best_row <- min_na(which(dds_log$objective_function == min_na(dds_log$objective_function)))

best_run <- dds_log[best_row, ]


load(file =  paste0(base_dir, "R/meltim/melt_calib/results_stations.Rdata"))

results_individual

plot(results_individual[best_row, ], stat_calib$altitude)
plot(stat_calib$altitude, results_individual[best_row, ])

#simu_best----

dds_log <- read.table(paste0(base_dir, "R/meltim/melt_calib/dds", ".log"), header = T, sep = "\t")

best_row <- min_na(which(dds_log$objective_function == min_na(dds_log$objective_function)))

best_run <- dds_log[best_row, ]; best_run

#Modify model parameters
# snow_params$a0                <- params[which(names(params) == "a0")]
# snow_params$a1                <- params[which(names(params) == "a1")]
snow_params$kSatSnow          <- best_run$kSatSnow
snow_params$densDrySnow       <- best_run$densDrySnow
snow_params$specCapRet        <- best_run$specCapRet
snow_params$emissivitySnowMin <- best_run$emissivitySnowMin
snow_params$emissivitySnowMax <- best_run$emissivitySnowMax
snow_params$tempAir_crit      <- best_run$tempAir_crit
snow_params$albedoMin         <- best_run$albedoMin
snow_params$albedoMax         <- best_run$albedoMax
snow_params$agingRate_tAirPos <- best_run$agingRate_tAirPos
snow_params$agingRate_tAirNeg <- best_run$agingRate_tAirNeg
snow_params$weightAirTemp     <- best_run$weightAirTemp
# # snow_params$tempAmpli         <- params[which(names(params) == "tempAmpli")]

if(is.vector(temps)){
  numb_stat <- 1
  meteo_length <- length(temps)
}else{
  numb_stat <- ncol(temps)
  meteo_length <- nrow(temps) 
}

snows <- foreach(k = 1:ncol(temps), .combine = 'cbind') %dopar% {
  
  swe_sim   <- rep(NA, meteo_length)
  sec_sim   <- rep(NA, meteo_length)
  alb_sim   <- rep(NA, meteo_length)
  sde_sim   <- rep(NA, meteo_length) #snow depth [m]
  lfr_sim   <- rep(NA, meteo_length) #snow depth [m]
  
  swe_init <- .0
  sec_init <- .0
  alb_init <- snow_params$albedoMax
  sde_init <- .0
  lfr_init <- .0
  
  swe_sim[1] <- swe_init
  sec_sim[1] <- sec_init
  alb_sim[1] <- alb_init
  sde_sim[1] <- sde_init
  lfr_sim[1] <- lfr_init
  
  if(numb_stat > 1){
    temps_sel <- temps[, k]
    precs_sel <- precs[, k]
  }else{
    temps_sel <- temps
    precs_sel <- precs
  }
  
  
  for(i in 2:meteo_length){
    
    sim_out <- snowModel_inter(
      #Forcings
      precipSumMM = precs_sel[i],
      shortRad = radi_mea_seri[i],
      tempAir = temps_sel[i],
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
    lfr_sim[i] <- sim_out[6]
    #Densitiy of dry snow = Density of snow * (1 - Liquid/Solid fraction)
    
    # densSnow <- snow_params$densDrySnow / (1 - sim_out[6])  
    # 
    # #maximum density = density of water
    # if(densSnow > 1000){
    #   densSnow <- 1000
    # }
    # 
    # sde_sim[i] <- sim_out[2] * (1/(densSnow/1000))
    
  }
  
  swe_sim
  
}

min_ind <- which(meteo_date == paste0(sta_yea_cal, "-01-01"))
max_ind <- which(meteo_date == paste0(end_yea_cal, "-12-31"))

snows_cal <- snows[min_ind:max_ind, ]
snows_stat_cal <- snows_stat[min_ind:max_ind, ]
meteo_date_cal <- meteo_date[min_ind:max_ind]

load(file =  paste0(base_dir, "R/meltim/melt_calib/results_stations.Rdata"))
results_individual[best_row+1, ]
mean(results_individual[best_row+1, ])

val_all <- NULL
for(i in 1:ncol(snows_cal)){
  
  val_cal <- obj_func(snows_cal[, i], snows_stat_cal[, i])
  val_all <- c(val_all, val_cal)
  print(val_cal)
  
}
mean(val_all)
for(i in 1:3){
  plot(snows_cal[, i], snows_stat_cal[, i])
  print(obj_func(snows_cal[, i], snows_stat_cal[, i]))
  abline(0, 1)
}

pdf(paste0(base_dir, "R/figs_exp/calib_res.pdf"), width = 16, height = 7.5)

par(mfrow = c(3, 1))
par(mar = c(2, 5.0, 2.5, 0.5))

main_plots <- c("a) Andermatt (1440 m)", "b) Davos (1560 m)", "c) Weissfluhjoch (2540 m)")

for(i in 1:ncol(snows_cal)){

  y_max <- max_na(c(snows_cal[, i], snows_stat_cal[, i]))
  y_min <- min_na(c(snows_cal[, i], snows_stat_cal[, i]))
  x_tics <- which(format(meteo_date_cal, '%m-%d') == '01-01')
  x_labs <- which(format(meteo_date_cal, '%m-%d') == '07-01')
  x_labels <- unique(format(meteo_date_cal, '%Y'))
  
  plot(snows_cal[, i], type = "n", ylim = c(y_min, y_max), axes = F, ylab = "", 
       xlab = "", xaxs = "i")
  abline(v = c(x_tics), lty = "dashed", col = "grey50", lwd = 0.8)
  lines(snows_cal[, i], col = scales::alpha("black", alpha = 1.0), lwd = 1.5)
  points(snows_stat_cal[, i], pch = 1, col = scales::alpha("darkblue", alpha = 1.0), cex = 1.5)
  axis(2, cex = 1.5, mgp=c(3, 0.50, 0), cex.axis = 1.6, tck = -0.02)
  axis(1, at = x_tics, labels = rep("", length(x_tics)), tck = -0.07)
  axis(1, at = x_labs, labels = x_labels, tick = F, mgp=c(3, 0.50, 0), cex.axis = 1.4)
  mtext(main_plots[i], side= 3, line = 0.3, cex = 1.4, adj = 0.0)
  mtext("SWE [m]", side= 2, line = 2.5, cex = 1.2, adj = 0.5)
  box()
  
}

dev.off()




plot(meteo_date, swe_all[, 1])
summary(swe_all)
length(swe_all)
length
length(meteo_date)
