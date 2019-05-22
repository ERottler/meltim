###

#Rhine snow - Runoff timing: phaselag
#Erwin Rottler, University of Potsdam
#Spring 2019

###

break_day <- 0
do_ma <- F
ma_window <- 30


#calc----

sta_yea_cla <- as.numeric(format(grdc_data$date[1], "%Y"))
end_yea_cla <- as.numeric(format(grdc_data$date[nrow(grdc_data)], "%Y"))


#Order data by day (including break day to set start hydrologica year)
data_day <- ord_day(data_in = grdc_data$value,
                    date = grdc_data$date,
                    start_y = sta_yea_cla,
                    end_y = end_yea_cla,
                    break_day = break_day,
                    do_ma = do_ma,
                    window_width = ma_window)

#Mean seasonal cycle discharge
dis_mea <- apply(data_day, 2, mea_na)
plot(dis_mea, type = "l")
lines(smoothFFT(dis_mea, sd = 15), col = "red3", lwd = 2)




#package_test----


# Set number of time steps per period (day)
nday = 24
#' set time unit minutes
timeunitperday = 24 * 60
dt = data.table(time = seq(0.5,timeunitperday, by = timeunitperday/nday))
dt[ , Date := as.IDate("2007-07-07")]
dt[ , ptime := time * 2* pi / timeunitperday]
dt[ , ptimecos := cos(ptime  - pi) ]
dt[ , Rsd := ifelse(ptimecos<0,0,800 * ptimecos)]
dt
# create a second time series without a time lag
dt[ , LEnolag :=  Rsd/2 + 30 * rnorm(1) , by = time ]
# create another series which has a lag of one hour
(lagLE = 60 *  2 * pi / timeunitperday)
dt[ , ptimecos_lag1h := cos(ptime  - pi - lagLE)  ]
dt[ , LElag1h :=  ifelse(ptimecos_lag1h<0,0, 400 * ptimecos_lag1h ) + 30 * rnorm(1)  , by = time ]

#Plot the generated series
plot(Rsd ~ time, data = dt)
lines(LEnolag ~ time, data = dt, col = 3)
lines(LElag1h ~ time, data = dt, col = 4)

## Hysteresis loops appear when plotted against Reference
plot(LEnolag ~ Rsd, data = dt, col = 3, type = "b")
lines(LElag1h ~ Rsd, data = dt, col = 4, type = "b", pch = 0)


# first estimate the time derivative of the reference series
dt[ ,dRsd := c(NA, diff(Rsd))]
# then call the estimation of the phase lags
dt[ , camuffo_phaselag_time(Y = LElag1h, X = Rsd, dX = dRsd, nday = nday, timeunitperday = timeunitperday)]
dt[ , camuffo_phaselag_time(Y = LEnolag, X = Rsd, dX = dRsd, nday = nday, timeunitperday = timeunitperday)]
## compare  phaselagtime with defined phase lag
## the significance of the lag can be assessed by slope2_pvalue

(camuffo_nolag =  dt[ , as.list(coef(lm(LEnolag ~ Rsd + dRsd))), by = Date])
(camuffo_lag1h =  dt[ , as.list(coef(lm(LElag1h ~ Rsd + dRsd))), by = Date])
### compute the phase lag in time units (here time)
phaselag_time(slope1 = camuffo_nolag[ ,Rsd], slope2 = camuffo_nolag[ , dRsd], nday = nday, timeunitperday = timeunitperday)
phaselag_time(slope1 = camuffo_lag1h[ ,Rsd], slope2 = camuffo_lag1h[ , dRsd], nday = nday, timeunitperday = timeunitperday)



#my_data----

sta_yea <- 1869
end_yea <- 2017
stat_sel <- "Basel_Rheinhalle_2" # Diepoldsau_2, Basel_Rheinhalle_2

#load discharge data
load(paste0(base_dir, "data/discharge/dis_new.RData"))
dat_sel_all <- dis_new[, which(names(dis_new) == stat_sel)]
dat_dis_all <- data.frame(date = dis_new$date, values = dat_sel_all)
dat_dis <- dat_dis_all[which((format(dat_dis_all$date, '%Y') >= sta_yea) & 
                               (format(dat_dis_all$date, '%Y') <= end_yea)), ]

#load temperature data
temp_sel <- read.table(paste0(base_dir, "/data/order64387/order_64387_data.txt"), sep = ";", skip = 2, header = T)
temp_sel$date <- as.Date(strptime(temp_sel$time, "%Y%m%d", tz="UTC"))
dat_temp_all <- data.frame(date = temp_sel$date, values  = temp_sel$ths200d0)
dat_temp <- dat_temp_all[which((format(dat_temp_all$date, '%Y') >= sta_yea) & (format(dat_temp_all$date, '%Y') <= end_yea)), ]

break_day <- 0
data_day <- order_day(date_in = dat_temp$date,
                      valu_in = dat_temp$values,
                      sta_yea = sta_yea,
                      end_yea = end_yea,
                      break_day = break_day)

#Mean seasonal cycle discharge
tem_mea <- apply(data_day, 2, mea_na)
#plot(tem_mea, type = "l")

f_dis_lag <- function(dat, ref, year_in = 2010, smooth_val = 10){
  
  row_sel <- which(format(dat$date, "%Y") == year_in)
  dis_y <- smoothFFT(dat$values[row_sel], sd = smooth_val)
  tem_y <- smoothFFT(ref, sd = smooth_val)
  
  if(length(dis_y) > 365){
    tem_y <- c(tem_y, tem_y[length(tem_y)])
  }
  
  
  dtem_y <- c(NA, diff(tem_y))
  
  my_slope_1 <- coef(lm(dis_y ~ tem_y))[2]
  my_slope_2 <- coef(lm(dis_y ~ dtem_y))[2]
  my_day <- 365 #number of measurements (per day)
  my_unit <- 365 #number of timesteps per day
  
  my_phaselag <- phaselag_time(slope1 = my_slope_1, slope2 = my_slope_2, nday = my_day, timeunitperday = my_unit)
  
  return(my_phaselag)
  
}

lag_years <-sta_yea:end_yea
my_lags <- rep(NA, length(lag_years))

for(i in 1:length(lag_years)){
  print(i)
  my_lags[i] <- f_dis_lag(dat = dat_dis, ref = tem_mea, 
                          year_in = lag_years[i], smooth_val = 5)
  
}

plot(my_lags, type = "l")
lines(smoothFFT(my_lags, sd = 5), col = "blue3")
abline(h = 0, lty = "dashed")

abline(v = 135)

sens_slope(my_lags)*10
sens_slope(my_lags)*length(lag_years)

