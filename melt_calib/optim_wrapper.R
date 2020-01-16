###

#Calibration Snow model - Optim wrapper
#Erwin Rottler, Summer 2019

###

optim_wrapper <- function(params = NULL){
  

  #Snow simulations----
  
  load(paste0(base_dir, "R/meltim/melt_calib/calib_thread/calib_snow_data.Rdata"))
  
  #Modify model parameters
  
    # snow_params$a0                <- params[which(names(params) == "a0")]
    # snow_params$a1                <- params[which(names(params) == "a1")]
    snow_params$kSatSnow          <- params[which(names(params) == "kSatSnow")]
    snow_params$densDrySnow       <- params[which(names(params) == "densDrySnow")]
    snow_params$specCapRet        <- params[which(names(params) == "specCapRet")]
    snow_params$emissivitySnowMin <- params[which(names(params) == "emissivitySnowMin")]
    snow_params$emissivitySnowMax <- params[which(names(params) == "emissivitySnowMax")]
    snow_params$tempAir_crit      <- params[which(names(params) == "tempAir_crit")]
    snow_params$albedoMin         <- params[which(names(params) == "albedoMin")]
    snow_params$albedoMax         <- params[which(names(params) == "albedoMax")]
    snow_params$agingRate_tAirPos <- params[which(names(params) == "agingRate_tAirPos")]
    snow_params$agingRate_tAirNeg <- params[which(names(params) == "agingRate_tAirNeg")]
    snow_params$weightAirTemp     <- params[which(names(params) == "weightAirTemp")]
    # # snow_params$tempAmpli         <- params[which(names(params) == "tempAmpli")]
  
    if(is.vector(temps)){
      numb_stat <- 1
      meteo_length <- length(temps)
    }else{
      numb_stat <- ncol(temps)
      meteo_length <- nrow(temps) 
    }

  snows <- foreach(k = 1:numb_stat, .combine = 'cbind') %dopar% {
    
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
  
  #Goodness_fit----
  
  source("melt_calib/obj_function.R")

  #selecte calibration period
  min_ind <- which(meteo_date == paste0(sta_yea_cal, "-01-01"))
  max_ind <- which(meteo_date == paste0(end_yea_cal, "-12-31"))
  
  obj_measure_all <- rep(NA, numb_stat)
  
  for(i in 1:numb_stat){
    
    if(numb_stat > 1){
      obj_measure_all[i] <- obj_func(snows[min_ind:max_ind, i], snows_stat[min_ind:max_ind, i])
    }else{
      obj_measure_all <- obj_func(snows[min_ind:max_ind], snows_stat[min_ind:max_ind])
    }
    
  }
    
  obj_measure <- mean(obj_measure_all)
  
  #Save individual results
  load(file =  paste0(base_dir, "R/meltim/melt_calib/results_stations.Rdata"))
  results_individual <- rbind(results_individual, obj_measure_all)
  save(results_individual, file =  paste0(base_dir, "R/meltim/melt_calib/results_stations.Rdata"), version = 2)

  return(obj_measure)
 
}

