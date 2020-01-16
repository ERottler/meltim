###

#Rhine snow - Runoff timing: classical approach
#Erwin Rottler, University of Potsdam
#Spring 2019

###

break_day <- 274 #1.Oct = 274, 1.Nov = 304
mid_yea_cal <- 1960
do_ma <-  F
ma_window <- 30

#calc----

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
day_cross_day <- apply(day_cross, 1, sens_slo) * 10 * -1 *length(sta_yea_cla:end_yea_cla)/10 # [days]
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
if(break_day == 274){
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

plot(day_cross_day, type = "l", axes = F, ylab = "", xlab = "")
points(day_cross_day, pch = 19, cex = 0.7)
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


#rasthy----

cols_min <- grDevices::colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "cadetblue3", "white"))(100)
cols_max <- grDevices::colorRampPalette(c("white", "yellow2","gold2", "orange2", "orangered3", "orangered4", "red4"))(100)
cols_hydro <- c(cols_min, cols_max)

max_break <- max_na(data_day)
min_break <- min_na(data_day)
qua_break <- quantile(data_day, probs = 0.6, type = 8, na.rm = T)

breaks_1 <- lseq(min_break, qua_break, length.out = length(cols_hydro)/2)
breaks_2 <- lseq(qua_break+0.01, max_break, length.out = length(cols_hydro)/2 + 1)
breaks_2[length(breaks_2)] <- breaks_2[length(breaks_2)] + 0.1

breaks_hydro <- c(breaks_1, breaks_2)

par(mar = c(1.6, 3, 1.6, 0.1))

layout(matrix(c(1,1,1,1,1,1,1,2),
              1, 8), widths=c(), heights=c())

image(x = 1:ncol(data_day),
      y = 1:nrow(data_day),
      z = t(data_day), 
      col = cols_hydro, 
      breaks = breaks_hydro,
      ylab = "", xlab = "", axes = F)

axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.02)#plot ticks
if(break_day == 304){
  axis(1, at = x_axis_lab, c("N","D","J","F","M","A","M","J","J","A","S","O"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.15, 0))#plot labels
}
if(break_day == 274){
  axis(1, at = x_axis_lab, c("O","N","D","J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.15, 0))#plot labels
}
mtext("Year", side = 2, line = 1.5, cex = 0.8)
axis(2, mgp=c(3, 0.15, 0), tck = -0.001)
mtext(data_meta, side = 3, line = 0.1, cex = 1.0)
box()

par(mar = c(1.6, 0.1, 1.6, 1.6))

alptempr::image_scale(as.matrix(data_day), col = cols_hydro, breaks = breaks_hydro, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.15, 0), tck = -0.08)
mtext("[m³/s]", side = 3, line = 0.1, cex = 0.8)

box()


#export----

data_bru <- data_mea
dat1_bru <- data_mea_1
dat2_bru <- data_mea_2
meta_bru <- data_meta

day_cross_bru <- day_cross
day_cross_mea_bru <- day_cross_mea
day_cross_slo_bru <- day_cross_slo
day_cross_day_bru <- day_cross_day

my_lags_bru <- my_lags
slo_dec_bru <- slo_dec
slo_all_bru <- slo_all
  
smea_band_die <- smea_band
sslo_band_die <- sslo_band
vmea_band_die <- vmea_band
vslo_band_die <- vslo_band
vdif_band_die <- vdif_band
vdis_band_die <- vdis_band
my_elev_bands_die <- my_elev_bands
  
tmea_band_die     <- tmea_band
tmea_band_mea_die <- tmea_band_mea
tslo_band_die     <- tslo_band
tslo_band_mea_die <- tslo_band_mea
pmea_band_die     <- pmea_band 
pmea_band_mea_die <- pmea_band_mea
pslo_band_die     <- pslo_band 
pslo_band_mea_die <- pslo_band_mea
emea_band_die     <- emea_band
emea_band_mea_die <- emea_band_mea
eslo_band_die     <- eslo_band
eslo_band_mea_die <- eslo_band_mea

visu_day_mar <- visu_day
visu_yea_mar <- visu_yea

load(paste0(base_dir, "R/meltim/meltim.Rdata"))

save(data_dom, dat1_dom, dat2_dom, meta_dom,
     day_cross_dom, day_cross_mea_dom, day_cross_slo_dom, day_cross_day_dom,
     my_lags_dom, slo_dec_dom, slo_all_dom,
     
     data_mar, dat1_mar, dat2_mar, meta_mar,
     day_cross_mar, day_cross_mea_mar, day_cross_slo_mar, day_cross_day_mar,
     my_lags_mar, slo_dec_mar, slo_all_mar,
     
     data_por, dat1_por, dat2_por, meta_por,
     day_cross_por, day_cross_mea_por, day_cross_slo_por, day_cross_day_por,
     my_lags_por, slo_dec_por, slo_all_por,
     
     data_die, dat1_die, dat2_die, meta_die,
     day_cross_die, day_cross_mea_die, day_cross_slo_die, day_cross_day_die,
     my_lags_die, slo_dec_die, slo_all_die,
     
     data_was, dat1_was, dat2_was, meta_was,
     day_cross_was, day_cross_mea_was, day_cross_slo_was, day_cross_day_was,
     my_lags_was, slo_dec_was, slo_all_was,
     
     data_bur, dat1_bur, dat2_bur, meta_bur,
     day_cross_bur, day_cross_mea_bur, day_cross_slo_bur, day_cross_day_bur,
     my_lags_bur, slo_dec_bur, slo_all_bur,
     
     data_neu, dat1_neu, dat2_neu, meta_neu,
     day_cross_neu, day_cross_mea_neu, day_cross_slo_neu, day_cross_day_neu,
     my_lags_neu, slo_dec_neu, slo_all_neu,
     
     data_rek, dat1_rek, dat2_rek, meta_rek,
     day_cross_rek, day_cross_mea_rek, day_cross_slo_rek, day_cross_day_rek,
     my_lags_rek, slo_dec_rek, slo_all_rek,
     
     data_bas, dat1_bas, dat2_bas, meta_bas,
     day_cross_bas, day_cross_mea_bas, day_cross_slo_bas, day_cross_day_bas,
     my_lags_bas, slo_dec_bas, slo_all_bas,
     
     data_bru, dat1_bru, dat2_bru, meta_bru,
     day_cross_bru, day_cross_mea_bru, day_cross_slo_bru, day_cross_day_bru,
     my_lags_bru, slo_dec_bru, slo_all_bru,
     
     data_ber, dat1_ber, dat2_ber, meta_ber,
     day_cross_ber, day_cross_mea_ber, day_cross_slo_ber, day_cross_day_ber,
     my_lags_ber, slo_dec_ber, slo_all_ber,
     
     data_mur, dat1_mur, dat2_mur, meta_mur,
     day_cross_mur, day_cross_mea_mur, day_cross_slo_mur, day_cross_day_mur,
     my_lags_mur, slo_dec_mur, slo_all_mur,
     
     data_unt, dat1_unt, dat2_unt, meta_unt,
     day_cross_unt, day_cross_mea_unt, day_cross_slo_unt, day_cross_day_unt,
     my_lags_unt, slo_dec_unt, slo_all_unt,
     
     data_mel, dat1_mel, dat2_mel, meta_mel,
     day_cross_mel, day_cross_mea_mel, day_cross_slo_mel, day_cross_day_mel,
     my_lags_mel, slo_dec_mel, slo_all_mel,
     
     data_max, dat1_max, dat2_max, meta_max,
     day_cross_max, day_cross_mea_max, day_cross_slo_max, day_cross_day_max,
     my_lags_max, slo_dec_max, slo_all_max,
     
     smea_band_bas, sslo_band_bas, vmea_band_bas, vslo_band_bas, 
     vdif_band_bas, vdis_band_bas, my_elev_bands_bas,
     
     tmea_band_bas, tmea_band_mea_bas, tslo_band_bas, tslo_band_mea_bas,
     pmea_band_bas, pmea_band_mea_bas, pslo_band_bas, pslo_band_mea_bas,
     emea_band_bas, emea_band_mea_bas, eslo_band_bas, eslo_band_mea_bas,
     
     smea_band_die, sslo_band_die, vmea_band_die, vslo_band_die, 
     vdif_band_die, vdis_band_die, my_elev_bands_die,
     
     tmea_band_die, tmea_band_mea_die, tslo_band_die, tslo_band_mea_die,
     pmea_band_die, pmea_band_mea_die, pslo_band_die, pslo_band_mea_die,
     emea_band_die, emea_band_mea_die, eslo_band_die, eslo_band_mea_die,
     
     meta_sel,
     
     visu_day_dom, visu_yea_dom, 
     visu_day_mar, visu_yea_mar,
  
     file = paste0(base_dir, "R/meltim/meltim.Rdata"))



plot(data_day[100, ], type = "l")
abline(v = c(54,62))


#exp_app----

save(data_day,
     
     file = paste0("U:/RhineFlow/rhine_snow/R/meltimr/inst/shiny_app/meltApp/data/meltApp.Rdata"))




