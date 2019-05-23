###

#Rhine snow - Runoff timing: phaselag
#Erwin Rottler, University of Potsdam
#Spring 2019

###

break_day_pha <- 0
do_ma <- F
ma_window <- 30

#calc_visu----

sta_yea_cla <- as.numeric(format(grdc_data$date[1], "%Y"))
end_yea_cla <- as.numeric(format(grdc_data$date[nrow(grdc_data)], "%Y"))

data_meta  <- paste0(gauge_sel, " ", sta_yea_cla, "-", end_yea_cla)

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
# plot(dis_mea, type = "l")
# lines(smoothFFT(dis_mea, sd = 5), col = "red3", lwd = 2)

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

slo_dec <- round(sens_slope(my_lags)*10, 2)
slo_all <- round(sens_slope(my_lags)*length(lag_years), 2)

par(mar = c(2.2, 2.2, 2, 0.2))

cex_dat <- 0.6

plot(my_lags, type = "l", axes = F, ylab = "", xlab = "")
lines(smoothFFT(my_lags, sd = 5), col = "blue3", lwd = 2)
axis(2, mgp=c(3, 0.15, 0), cex.axis = 0.7)
axis(1, at = 1:length(my_lags),labels = 1:length(my_lags), mgp=c(3, 0.15, 0), cex.axis = 0.7)
mtext("Percentage Discharge", side = 1, line = 1.3)
mtext("Days earlier [day]", side = 2, line = 1.3)
mtext(data_meta, line = 0.2, side = 3, cex = 0.8)
abline(h = 0, lty = "dashed")
mtext(paste(slo_dec, "day/dec"), side = 3, line = 0.2, adj = 0, cex = cex_dat)
mtext(paste(slo_all, "days"), side = 3, line = 0.2, adj = 1, cex = cex_dat)
mtext(data_meta, line = 0.2, side = 3, cex = 0.8)
box(lwd = 0.7)

# ind_sel <- which(format(grdc_data$date, "%Y") == lag_years[which(my_lags == min_na(my_lags))]) 
# plot(grdc_data$value[ind_sel], type = "l")
# length(my_lags)
