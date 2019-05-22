###

#Rhine snow - Runoff timing: classical approach
#Erwin Rottler, University of Potsdam
#Spring 2019

###

break_day <- 304 #1.Oct = 274, 1.Nov = 304
mid_yea_cal <- 1960
gauge_sel <- "Thoerishaus" # Domat, Martinsbruck, Porte_du_Scex, Wasserburg, Burghausen, 
                         # Diepoldsau, Neuhausen, Andelfingen, Appenzell, Rekingen, Mellingen, Bern_Schoenau, 
                         # Emmenmatt, Montier, Murgenthal, Brugg, Untersiggenthal, Basel, Maxau, Koeln, Wuerzburg, Schwaibach, Thoerishaus
do_ma <-  F
ma_window <- 30

#calc----
if(gauge_sel == "Domat"){          file_sel <- "6935145_Q_Day.Cmd.txt"}
if(gauge_sel == "Martinsbruck"){   file_sel <- "6943100_Q_Day.Cmd.txt"}
if(gauge_sel == "Porte_du_Scex"){  file_sel <- "6939200_Q_Day.Cmd.txt"}
if(gauge_sel == "Wasserburg"){     file_sel <- "6343100_Q_Day.Cmd.txt"}
if(gauge_sel == "Burghausen"){     file_sel <- "6343500_Q_Day.Cmd.txt"}
if(gauge_sel == "Diepoldsau"){     file_sel <- "6935500_Q_Day.Cmd.txt"}
if(gauge_sel == "Mellingen"){      file_sel <- "6935310_Q_Day.Cmd.txt"}
if(gauge_sel == "Untersiggenthal"){file_sel <- "6935300_Q_Day.Cmd.txt"}
if(gauge_sel == "Neuhausen"){      file_sel <- "6935055_Q_Day.Cmd.txt"}
if(gauge_sel == "Andelfingen"){    file_sel <- "6935400_Q_Day.Cmd.txt"}
if(gauge_sel == "Rekingen"){       file_sel <- "6935054_Q_Day.Cmd.txt"}
if(gauge_sel == "Bern_Schoenau"){  file_sel <- "6935020_Q_Day.Cmd.txt"}
if(gauge_sel == "Montier"){        file_sel <- "6935060_Q_Day.Cmd.txt"}
if(gauge_sel == "Appenzell"){      file_sel <- "6935412_Q_Day.Cmd.txt"}
if(gauge_sel == "Murgenthal"){     file_sel <- "6935302_Q_Day.Cmd.txt"}
if(gauge_sel == "Brugg"){          file_sel <- "6935301_Q_Day.Cmd.txt"}
if(gauge_sel == "Basel"){          file_sel <- "6935051_Q_Day.Cmd.txt"}
if(gauge_sel == "Schwaibach"){     file_sel <- "6335125_Q_Day.Cmd.txt"}
if(gauge_sel == "Maxau"){          file_sel <- "6335200_Q_Day.Cmd.txt"}
if(gauge_sel == "Koeln"){          file_sel <- "6335060_Q_Day.Cmd.txt"}
if(gauge_sel == "Wuerzburg"){      file_sel <- "6335500_Q_Day.Cmd.txt"}

grdc_data <- read_grdc(paste0(grdc_dir, file_sel)) ; summary(grdc_data$value)

#define start year
sta_yea_cla <- as.numeric(format(grdc_data$date[1], "%Y"))

# if(length(which(is.na(grdc_data$value))) > 0){
#   sta_yea_cla <- as.numeric(format(grdc_data$date[max_na(which(is.na(grdc_data$value)))], "%Y")) + 1
# }else{
#   sta_yea_cla <- as.numeric(format(grdc_data$date[1], "%Y"))
# }

# if(format(grdc_data$date[min_na(which(format(grdc_data$date, "%Y") == sta_yea_cla))], "%m-%d") == "01-01"){
#   sta_yea_cla <- sta_yea_cla
# }else{
#   sta_yea_cla <- sta_yea_cla + 1
# }

end_yea_cla <- as.numeric(format(grdc_data$date[nrow(grdc_data)], "%Y"))

# if(format(grdc_data$date[nrow(grdc_data)], "%m-%d") == "12-31"){
#   end_yea_cla <- as.numeric(format(grdc_data$date[nrow(grdc_data)], "%Y"))
# }else{
#   end_yea_cla <- as.numeric(format(grdc_data$date[nrow(grdc_data)], "%Y")) - 1
# }


#Order data by day (including break day to set start hydrologica year)
data_day <- ord_day(data_in = grdc_data$value,
                    date = grdc_data$date,
                    start_y = sta_yea_cla,
                    end_y = end_yea_cla,
                    break_day = break_day,
                    do_ma = do_ma,
                    window_width = ma_window)

#Mean seasonal cycle discharge
data_mea <- apply(data_day, 2, mea_na)

# break_year <- floor(nrow(data_day) / 2) #Mean cycle fist/second half of time frame
break_year <- nrow(data_day) - (end_yea_cla - mid_yea_cal)

data_mea_1 <- apply(data_day[c(1:break_year), ], 2, mea_na)
data_mea_2 <- apply(data_day[c((break_year+1):nrow(data_day)), ], 2, mea_na)
data_meta  <- paste0(gauge_sel, " ", sta_yea_cla, "-", end_yea_cla)


#only year with complete data recordings
na_check <- function(data_in){
  
  if(length(which(is.na(data_in))) > 0){
    data_in <- rep(NA, length(data_in))
  }
  
  return(data_in)

}

for(i in 1:nrow(data_day)){
  
  data_day[i, ] <- na_check(data_day[i, ])
  
}


#Cumulative sum discharges per year
data_cumsum <- apply(data_day, 1, cumsum)

# plot(data_cumsum[, 1], type = "l", ylim = c(0, max_na(data_cumsum)))
# for(i in 1:ncol(data_cumsum)){
# 
#   lines(data_cumsum[, i], col = "black")
# 
# }

#Scale cummulative sums

#Devide array by last element (Scale cumsums)
cumsum_scale <- function(data_in){
  
  if(length(which(is.na(data_in))) > 0){
    data_out <- rep(NA, length(data_in))
  }else{
    data_out <- data_in / data_in[length(data_in)]
  }
    
  return(data_out)
  
}

data_cumsum_scale <- apply(data_cumsum, 2, cumsum_scale)

# plot(data_cumsum_scale[, 1], type = "l", ylim = c(0, max_na(data_cumsum_scale)))
# for(i in 1:ncol(data_cumsum_scale)){
# 
#   lines(data_cumsum_scale[, i], col = "black")
# 
# }

#DOY percentage discharge through

percents <- seq(0.01, 0.99, 0.01)
day_cross <- matrix(nrow = length(percents), ncol = ncol(data_cumsum_scale))

for(p in 1:length(percents)){
  
  for(i in 1:ncol(data_cumsum_scale)){
    
    if(length(which(is.na(data_cumsum_scale[, i]))) > 0){
      day_cross[p, i] <- NA
    }else{
      day_cross[p, i] <- day_cross[p, i] <- min(which(data_cumsum_scale[, i] > percents[p]))
    }
  }
}

#Slope and mean of crossing days
#Mann-Kendall significance test (with prewithening)
mk_sig <- function(data_in, cover_thresh = 0.9){
  
  if(length(which(is.na(data_in))) / length(data_in) > (1-cover_thresh)){
    mk_sig <-  NA
  }else{
    time_step <- 1:length(data_in)
    # sens_slo <- as.numeric(zyp.sen(data_in~time_step)$coefficients[2])
    mk_sig <- as.numeric(zyp.trend.vector(data_in, method = "zhang", conf.intervals = F)[6])
  }
  
  return(mk_sig)
  
}

day_cross_slo <- apply(day_cross, 1, sens_slo) * 10 * -1 # [day/dec]
day_cross_sig <- apply(day_cross, 1, mk_sig)
day_cross_mea <- apply(day_cross, 1, mea_na)

#visu_1----

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15

y_max <- max(c(smoothFFT(data_mea, sd =5), smoothFFT(data_mea_1, sd =5), smoothFFT(data_mea_2, sd =5)))
y_min <- min(c(smoothFFT(data_mea, sd =5), smoothFFT(data_mea_1, sd =5), smoothFFT(data_mea_2, sd =5)))

plot(smoothFFT(data_mea_1, sd = 5), type = "n", col ="blue3", 
     axes = F, ylab = "", xlab = "", ylim = c(y_min, y_max))
lines(smoothFFT(data_mea_1, sd =5), col = "blue3", lwd = 2)
lines(smoothFFT(data_mea, sd =5), col = "black", lwd = 2)
lines(smoothFFT(data_mea_2, sd =5), col = "red3", lwd = 2)
axis(2, mgp=c(3, 0.15, 0))
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
if(break_day == 304){
  axis(1, at = x_axis_lab, c("N","D","J","F","M","A","M","J","J","A","S","O"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.15, 0))#plot labels
}
if(break_day == 374){
  axis(1, at = x_axis_lab, c("O","N","D","J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.15, 0))#plot labels
}

box(lwd = 0.7)
legend("topleft", c("start - 1960", "all", "1961 - end"), pch = 19, col = c("blue3", "black", "red3"))
mtext(data_meta, line = 0.2, side = 3, cex = 1.5)

#visu_2----

#Plot: Selected percentage values
smo_val <- 0.25
cex_dat <- 0.6
decs <- length(sta_yea_cla:end_yea_cla)/10

par(mar = c(1.5, 1.5, 4, 0.1))
plot(day_cross[75, ], type = "n", ylim = c(min_na(day_cross[25, ]), max_na(day_cross[75, ])),
     axes = F, ylab = "", xlab = "")
                                           
lines(day_cross[75, ], col = "blue3")
lines(day_cross[50, ], col = "black")
lines(day_cross[25, ], col = "red3")
lines(loess_NA_restore(day_cross[75, ], smoo_val = smo_val), col = "blue3")
lines(loess_NA_restore(day_cross[50, ], smoo_val = smo_val), col = "black")
lines(loess_NA_restore(day_cross[25, ], smoo_val = smo_val), col = "red3")
axis(2, mgp=c(3, 0.15, 0), cex.axis = 0.7)
axis(1, at = 1:length(sta_yea_cla:end_yea_cla),labels = sta_yea_cla:end_yea_cla, mgp=c(3, 0.15, 0), cex.axis = 0.7)
box(lwd = 0.6)
mtext(data_meta, side = 3, line = 3, cex = 0.8)

mtext(paste0(round(day_cross_slo[75], 2), " day/dec"), side = 3, line = 1.0, adj = 0.0, col = "blue3", cex = cex_dat)
mtext(paste0(round(day_cross_slo[50], 2), " day/dec"), side = 3, line = 1.0, adj = 0.5, col = "black", cex = cex_dat)
mtext(paste0(round(day_cross_slo[25], 2), " day/dec"), side = 3, line = 1.0, adj = 1.0, col = "red3", cex = cex_dat)
mtext(paste0(round(day_cross_mea[75], 2), " DOY"), side = 3, line = 2.0, adj = 0.0, col = "blue3", cex = cex_dat)
mtext(paste0(round(day_cross_mea[50], 2), " DOY"), side = 3, line = 2.0, adj = 0.5, col = "black", cex = cex_dat)
mtext(paste0(round(day_cross_mea[25], 2), " DOY"), side = 3, line = 2.0, adj = 1.0, col = "red3", cex = cex_dat)
mtext(paste0(round(day_cross_slo[75], 2) * decs, " days"), side = 3, line = 0.1, adj = 0.0, col = "blue3", cex = cex_dat)
mtext(paste0(round(day_cross_slo[50], 2) * decs, " days"), side = 3, line = 0.1, adj = 0.5, col = "black", cex = cex_dat)
mtext(paste0(round(day_cross_slo[25], 2) * decs, " days"), side = 3, line = 0.1, adj = 1.0, col = "red3", cex = cex_dat)


#visu_3----

#Plot: All precentages with trends
par(mar = c(2.5, 2.5, 1.5, 1.7))

plot(day_cross_slo * decs, type = "l", axes = F, ylab = "", xlab = "")
points(day_cross_slo * decs, pch = 19, cex = 0.7)
abline(h = 0, lty = "dashed")
axis(2, mgp=c(3, 0.15, 0), cex.axis = 0.7)
axis(1, mgp=c(3, 0.15, 0), cex.axis = 0.7)
par(new = T)
plot(day_cross_mea, type = "l", col = "red3", axes = F, ylab = "", xlab = "")
points(day_cross_mea, pch = 19, cex = 0.7, col = "red3")
box()
axis(4,  mgp=c(3, 0.15, 0), cex.axis = 0.7, col = "red3", col.axis = "red3")
mtext("Percentage Discharge", side = 1, line = 1.3)
mtext("Days earlier [day]", side = 2, line = 1.3)
mtext("Average DOY", side = 4, line = 1.3, col = "red3")
mtext(data_meta, side = 3, line = 0.1, cex = 1.0)

#visu_4----

day_sel <- 30 + 31 + 31 + 28 #1 Feb is 120 with break day of 1 Nov
year_sel <- which(sta_yea_cla:end_yea_cla == 2006)
smo_val <- 0.25
visu_day <- data_day[, day_sel]
visu_yea <- data_day[year_sel, ]

par(mar = c(2.5, 2.5, 1.5, 0.1))

plot(visu_day, type = "n", ylab = "", xlab = "", axes = F)
lines(visu_day)
lines(loess_NA_restore(visu_day, smoo_val = smo_val), lwd = 2)
axis(2, mgp=c(3, 0.15, 0), cex.axis = 0.7)
axis(1, mgp=c(3, 0.15, 0), cex.axis = 0.7)
box()
mtext("Year", side = 1, line = 1.3)
mtext("Discharge [m³/s]", side = 2, line = 1.3)
mtext(data_meta, side = 3, line = 0.1, cex = 1.0)
mtext("1 Feb", side = 3, line = 0.1, cex = 0.7, adj = 0.0)


plot(visu_yea, type = "n", ylab = "", xlab = "", axes = F)
lines(visu_yea)
lines(loess_NA_restore(visu_yea, smoo_val = smo_val), lwd = 2)
axis(2, mgp=c(3, 0.15, 0), cex.axis = 0.7)
axis(1, mgp=c(3, 0.15, 0), cex.axis = 0.7)
box()
mtext("DOY", side = 1, line = 1.3)
mtext("Discharge [m³/s]", side = 2, line = 1.3)
mtext(data_meta, side = 3, line = 0.1, cex = 1.0)
mtext("01.11.2006- 31.10.2007", side = 3, line = 0.1, cex = 0.7, adj = 0.0)


#export----

dat19   <- data_mea
dat19_1 <- data_mea_1
dat19_2 <- data_mea_2
meta_19 <- data_meta

day_cross_2 <- day_cross
day_cross_mea_2 <- day_cross_mea
day_cross_slo_2 <- day_cross_slo
meta_cross_2 <- data_meta

load( paste0(base_dir, "R/meltim/riv_flow_new.Rdata"))

save(dat1, dat1_1, dat1_2, meta_1, dat2, dat2_1, dat2_2, meta_2,
     dat3, dat3_1, dat3_2, meta_3, dat4, dat4_1, dat4_2, meta_4,
     dat5, dat5_1, dat5_2, meta_5, dat6, dat6_1, dat6_2, meta_6,
     dat7, dat7_1, dat7_2, meta_7, dat8, dat8_1, dat8_2, meta_8,
     dat9, dat9_1, dat9_2, meta_9, 
     dat10, dat10_1, dat10_2, meta_10, dat11, dat11_1, dat11_2, meta_11,
     dat12, dat12_1, dat12_2, meta_12, dat13, dat13_1, dat13_2, meta_13,
     dat14, dat14_1, dat14_2, meta_14, dat15, dat15_1, dat15_2, meta_15,
     dat16, dat16_1, dat16_2, meta_16, dat17, dat17_1, dat17_2, meta_17,
     dat18, dat18_1, dat18_2, meta_18, dat19, dat19_1, dat19_2, meta_19,
     
     day_cross_1, day_cross_mea_1, day_cross_slo_1, meta_cross_1,
     day_cross_2, day_cross_mea_2, day_cross_slo_2, meta_cross_2,
     
     file = paste0(base_dir, "R/meltim/riv_flow_new.Rdata"))



plot(data_day[100, ], type = "l")
abline(v = c(54,62))
