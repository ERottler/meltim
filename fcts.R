###

#Additional Functions for melTim

###

# Extract coordinates from Polygon
# Extract coordinates from SpatialPolyognsDataFrame
extract_coords <- function(sp_df){
  
  results <- list()
  
  for(i in 1:length(sp_df@polygons[[1]]@Polygons)){
    results[[i]] <- sp_df@polygons[[1]]@Polygons[[i]]@coords
  }
  
  results <- Reduce(rbind, results)
  
  return(results)
  
}


# Buffer mean value
# Get mean elevation of squared buffer around grid point.
elev_buff <- function(point_in, radius = 2500, dem_in){
  
  y_plu <-point_in@bbox[2,1] + radius
  x_plu <-point_in@bbox[1,1] + radius
  y_min <-point_in@bbox[2,1] - radius
  x_min <-point_in@bbox[1,1] - radius
  
  square <- cbind(x_min, y_plu, #NW corner
                  x_plu, y_plu, #NE corner
                  x_plu, y_min, #SW corner
                  x_min, y_min, #SW corner
                  x_min, y_plu) #NW corner again to close polygon
  
  pol_sq <- sp::Polygon(matrix(square, ncol = 2, byrow = T))
  sq_buf <- sp::SpatialPolygons(list(sp::Polygons(list(pol_sq), ID = "sq_buf")),
                                proj4string = sp::CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  
  
  #get mean elevation in squared buffer...
  
  elev_square <- raster::extract(dem_in, sq_buf, fun = mean, na.rm = T)
  
  return(elev_square)
  
}


# Grid around spatial point.
# Add points as 1 km grid around spatial point. One point turn into a spatial grind with 1 km  resolution. In total 25 spatial points returned.
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
  d_09[2] <- d_07[2] + res_new * 0
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

# Cube index from points inside polygon.
# Get column number of matching point.
get_cube_index_col <- function(val_in, lons_in = lon2D, col_or_row = "col"){
  
  if(col_or_row == "col"){
    index_out <- which(round(lons_in, digits =6) == round(val_in, digits =6), arr.ind = T)[1,1]
  }
  
  if(col_or_row == "row"){
    index_out <- which(round(lons_in, digits =6) == round(val_in, digits =6), arr.ind = T)[1,2]
  }
  
  return(index_out)
  
}


# Cube index from points inside polygon.
# Get row number of matching point.
get_cube_index_row <- function(val_in, lons_in = lon2D, col_or_row = "row"){
  if(col_or_row == "col"){
    index_out <- which(round(lons_in, digits =6) == round(val_in, digits =6), arr.ind = T)[1,1]
  }
  
  if(col_or_row == "row"){
    index_out <- which(round(lons_in, digits =6) == round(val_in, digits =6), arr.ind = T)[1,2]
  }
  
  return(index_out)
}

#Phase lags
#Calculations of phase lags between mean hydrograph and individual years using package: devtools::install_github("laubblatt/phaselag").
f_phase_lag <- function(grdc_data, break_day, end_yea_cla, sta_yea_cla){
  sta_yea_cla <- as.numeric(format(grdc_data$date[1], "%Y"))
  end_yea_cla <- as.numeric(format(grdc_data$date[nrow(grdc_data)], "%Y"))
  
  #Order data by day (including break day to set start hydrologica year)
  data_day <- ord_day(data_in = grdc_data$value,
                      date = grdc_data$date,
                      start_y = sta_yea_cla,
                      end_y = end_yea_cla,
                      break_day = break_day_pha,
                      do_ma = do_ma,
                      window_width = ma_window)
  
  #Mean seasonal cycle discharge
  dis_mea <- apply(data_day, 2, mea_na)
  
  #calculate lag between mean cycle and all other yearly cycles
  dis_lag <- function(dat, ref, year_in = 2010, smo_ref = -1, smo_dat = -1){
    
    ind_sel <- which(format(dat$date, "%Y") == year_in)
    
    if(length(ind_sel) < 365 ){
      
      my_phaselag <- NA
      
    }else{
      
      if(length(which(is.na(dat$value[ind_sel])) > 0)){
        
        my_phaselag <- NA
        
      }else{
        
        dis_y <- smoothFFT(dat$value[ind_sel], sd = smo_dat)
        ref <- smoothFFT(ref, sd = smo_ref)
        
        if(length(dis_y) > 365){
          ref <- c(ref, ref[length(ref)])
        }
        
        dref <- c(NA, diff(ref))
        
        my_slope_1 <- coef(lm(dis_y ~ ref))[2]
        my_slope_2 <- coef(lm(dis_y ~ dref))[2]
        my_day <- 365 #number of measurements per period
        my_unit <- 365 #number of time steps per period
        
        my_phaselag <- phaselag_time(slope1 = my_slope_1, slope2 = my_slope_2, nday = my_day, timeunitperday = my_unit)
        # plot(dis_y, type = "l")
      }
      
    }
    
    return(my_phaselag * -1) #multiply with -1 so that positive values = earlier
    
  }
  
  lag_years <-sta_yea_cla:end_yea_cla
  my_lags <- rep(NA, length(lag_years))
  
  for(i in 1:length(lag_years)){
    print(lag_years[i])
    my_lags[i] <- dis_lag(dat = grdc_data, ref = dis_mea, year_in = lag_years[i],
                          smo_ref = -1, smo_dat = -1)
    
  }
  
  return(my_lags)
  
}