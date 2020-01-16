###

#Export figures for manusript melTim
#Erwin Rottler, University of Potsam

###

pacman::p_load(devtools, leaflet, raster, tmap, sf, prettymapr, meltimr, alptempr, rfs, viridis, 
               shape, scales, emdbook, zoo, zyp)

#set base direcoty
base_dir <- "U:/rhine_snow/"

#Load results snow simulations
load(paste0(base_dir, "R/draft_snow_17_12.RData"))

#snow_stats----

sta_yea <- 1959
end_yea <- 2018

#Read snow station data
snow_data_1 <- read.table(paste0(base_dir,"data/snow/order_73106_data.txt"), sep = ";", skip = 2, stringsAsFactors = F, header = T, na.strings = "-")
snow_data_2 <- read.table(paste0(base_dir,"data/snow/order_73164_data.txt"), sep = ";", skip = 2, stringsAsFactors = F, header = T, na.strings = "-")
snow_data_3 <- read.table(paste0(base_dir,"data/snow/order_74297_data.txt"), sep = ";", skip = 2, stringsAsFactors = F, header = T, na.strings = "-")
snow_data_all <- rbind(snow_data_1, snow_data_2, snow_data_3)

id_all <- unique(snow_data_all$stn); id_all <- id_all[-which(id_all == "stn")]
snow_all <- NULL
stat_col <- NULL

for(i in 1: length(id_all)){
  
  print(i)
  
  sel_ind <- which(snow_data_all$stn == id_all[i])
  snow_valu <- as.numeric(snow_data_all$hto000d0[sel_ind])
  snow_date <- as.Date(strptime(snow_data_all$time[sel_ind], "%Y%m%d", tz = "UTC"))
  
  snow_data <- data.frame(date = snow_date,
                          valu = snow_valu)
  
  #start/end year selected snow recordings
  sta_yea_sno <- as.numeric(format(snow_date[1], "%Y"))
  end_yea_sno <- as.numeric(format(snow_date[length(snow_date)], "%Y"))
  
  #1.Criteria: Recordings at least from 1969
  if(sta_yea_sno <= 1969){
    
    #Fill possible gaps
    sta_day <- paste0(sta_yea, "-01-01")
    end_day <- paste0(end_yea, "-12-31")
    
    sta_date <- as.POSIXct(strptime(sta_day, "%Y-%m-%d", tz="UTC"))
    end_date <- as.POSIXct(strptime(end_day, "%Y-%m-%d", tz="UTC"))
    snow_stat_date  <- as.Date(seq(sta_date, end_date, by="day"))
    
    snow_sel <- data.frame(date  = snow_stat_date,
                           values = with(snow_data, valu[match(as.Date(snow_stat_date), as.Date(date))]))
    
    if(id_all[i] == "SAE"){#Data Saentis only from 1965 on
      
      snow_sel$values[which(format(snow_sel$date, "%Y") < 1965)] <- NA
      
    }
    
    #2.Criteria: not more than 15 years in total missing
    if(length(which(is.na(snow_sel$values))) < length(snow_sel$values) *(15/60)){
      
      snow_all <- cbind(snow_all, snow_sel$values)
      stat_col <- c(stat_col, id_all[i])
      
    }
    
  }
  
}

colnames(snow_all) <- stat_col


#Plot observational stations

pdf("U:/rhine_snow/R/figs_exp/snow_stats.pdf", width = 12, height = 14.5)

par(family = "serif")

sing_snow <- function(stat_sel, do_head = F, sta_yea = 1958, bre_yea = 1988, end_yea = 2018, bre_day = 274){
  
  snow_stat <- snow_all[, which(colnames(snow_all) == stat_sel)]
  
  #Calculate mean yearly cycles
  data_day_1 <- ord_day(data_in = snow_stat,
                        date = snow_stat_date,
                        start_y = sta_yea,
                        end_y = bre_yea,
                        break_day = bre_day,
                        do_ma = T,
                        window_width = 30)
  
  data_day_2 <- ord_day(data_in = snow_stat,
                        date = snow_stat_date,
                        start_y = bre_yea,
                        end_y = end_yea,
                        break_day = bre_day,
                        do_ma = T,
                        window_width = 30)
  
  sno_mea_1 <- apply(data_day_1, 2, mea_na)
  sno_mea_2 <- apply(data_day_2, 2, mea_na)
  
  #Snow meanaccumulation/melt rates
  snow_stat_diff <- c(diff(snow_stat), NA)
  
  data_day_diff_1 <- ord_day(data_in = snow_stat_diff,
                             date = snow_stat_date,
                             start_y = sta_yea,
                             end_y = bre_yea,
                             break_day = bre_day,
                             do_ma = T,
                             window_width = 30)
  
  data_day_diff_2 <- ord_day(data_in = snow_stat_diff,
                             date = snow_stat_date,
                             start_y = bre_yea,
                             end_y = end_yea,
                             break_day = bre_day,
                             do_ma = T,
                             window_width = 30)
  
  diff_mea_1 <- apply(data_day_diff_1, 2, mea_na)
  diff_mea_2 <- apply(data_day_diff_2, 2, mea_na)
  
  #vertical ablines at maximum/zero
  doy_max_1 <- which(sno_mea_1 == max_na(sno_mea_1))
  doy_max_2 <- which(sno_mea_2 == max_na(sno_mea_2))
  
  
  #Plot 1: Raster graph
  
  data_day <- ord_day(data_in = snow_stat,
                      date = snow_stat_date,
                      start_y = sta_yea,
                      end_y = end_yea,
                      break_day = bre_day,
                      do_ma = F,
                      window_width = 30)
  
  x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
  x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
  
  # cols_hydro <- c(colorRampPalette(c("white", viridis::viridis(20, direction = -1)))(200))
  # cols_hydro <- colorRampPalette(c("grey85", "cadetblue3", viridis::viridis(9, direction = 1)[4:1]))(200)
  cols_hydro <- colorRampPalette(c("grey90", viridis::viridis(9, direction = 1)[4:1]))(200)
  breaks_hydro <- seq(alptempr::min_na(data_day), alptempr::max_na(data_day), length.out = length(cols_hydro)+1)
  
  par(mar = c(1.7, 1.7, 0.5, 0.2))
  
  image(x = 1:ncol(data_day),
        y = 1958:2017,
        z = t(data_day),
        col = cols_hydro,
        breaks = breaks_hydro,
        ylab = "", xlab = "", axes = F)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.08)#plot ticks
  axis(1, at = x_axis_lab, c("O","N","D","J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.2)#plot labels
  axis(2, mgp=c(3, 0.10, 0), tck = -0.015, cex.axis = 1.2)
  if(do_head){
    mtext("a) Snow depth raster", side= 3, line = 0.5, cex = 1.2, adj = 0.0)
    mtext("[cm]", side= 3, line = 0.2, cex = 0.9, adj = 1.0)
    # legend("topleft", c("1959-1988", "1989-2018"), pch = 19, cex = 1.0, col = c(col_1, col_2), bg = "white")
  }
  box()
  
  par(mar = c(1.7, 0.2, 0.5, 1.5))
  
  alptempr::image_scale(as.matrix(data_day), col = cols_hydro, breaks = breaks_hydro, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
  axis(4, mgp=c(3, 0.20, 0), tck = -0.10, cex.axis = 1.2)
  # mtext("cm", side = 3, line = 0.2, cex = 1.0)
  box()
  
  
  #Plot 2: Mean annual cycle 
  
  par(mar =c(1.7, 1.7, 0.5, 1.0))
  
  x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
  x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
  
  col_1 <- "steelblue4"
  col_2 <- "darkred"
  col_3 <- "black"
  my_lwd <- 2.5
  lwd_ab <- 2.5
  
  y_max <- max(c(sno_mea_1, sno_mea_2))
  y_min <- min(c(sno_mea_1, sno_mea_2))
  
  plot(sno_mea_1, type = "n", col =col_1, axes = F, ylab = "", xlab = "", ylim = c(y_min, y_max))
  abline(v = x_axis_tic, col = "grey55",  lty = "dashed", lwd = 0.7)
  lines(sno_mea_1, col = col_1, lwd = my_lwd)
  lines(sno_mea_2, col = col_2, lwd = my_lwd)
  # abline(v = doy_max_1, col = col_1, lty = "dashed", lwd = lwd_ab)
  # abline(v = doy_max_2, col = col_2,  lty = "dashed", lwd = lwd_ab)
  axis(2, mgp=c(3, 0.05, 0), tck = -0.015, cex.axis = 1.2)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.08)#plot ticks
  axis(1, at = x_axis_lab, c("O","N","D","J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.2)#plot labels
  if(do_head){
    mtext("b) Snow depth mean", side= 3, line = 0.5, cex = 1.2, adj = 0.0)
    mtext("[cm]", side= 3, line = 0.2, cex = 0.9, adj = 1.0)
    legend("topleft", c("1958-1987", "1988-2017"), pch = 19, cex = 1.1, col = c(col_1, col_2), bg = "white")
  }
  box(lwd = 0.7)
  
  
  #Plot 3: Annaul cycle meand diff.
  
  y_max <- max(c(diff_mea_1, diff_mea_2))
  y_min <- min(c(diff_mea_1, diff_mea_2))
  
  plot(diff_mea_1, type = "n", col = col_1, axes = F, ylab = "", xlab = "", ylim = c(y_min, y_max))
  abline(h = 0, col = "grey55", lty = "dashed", lwd = 0.7)
  # abline(v = doy_max_1, col = col_1, lty = "dashed", lwd = lwd_ab)
  # abline(v = doy_max_2, col = col_2, lty = "dashed", lwd = lwd_ab)
  abline(v = x_axis_tic, col = "grey55",  lty = "dashed", lwd = 0.7)
  lines(diff_mea_1, col = col_1, lwd = my_lwd)
  lines(diff_mea_2, col = col_2, lwd = my_lwd)
  axis(2, mgp=c(3, 0.05, 0), tck = -0.015, cex.axis = 1.2)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.08)#plot ticks
  axis(1, at = x_axis_lab, c("O","N","D","J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.2)#plot labels
  if(do_head){
    mtext("c) Acc./Melt rate", side= 3, line = 0.5, cex = 1.2, adj = 0.0)
    mtext("[cm/1d]", side= 3, line = 0.2, cex = 0.9, adj = 1.0)
  }
  # abline(v = x_axis_tic, col = "grey55", lty = "dashed", lwd = 0.8)
  # grid(nx = 0, ny = 5, lty = "dashed", col = "grey55", lwd = 0.8)
  box(lwd = 0.7)
  
}

layout(matrix(c(46, rep(45, 11),  
                rep(c(46, 1, 5,  9, 13, 17, 21, 25, 29, 33, 37, 41), 8),
                46, 2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42,
                rep(c(46, 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43), 9), 
                rep(c(46, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44), 9)
),
12, 28), widths=c(), heights=c(0.2, rep(1, 10)))

x_axis_lab <- c(15,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(   46,74,105,135,166,196,227,258,288,319,349)-15

sing_snow("WFJ", do_head = T)
sing_snow("ARO")
sing_snow("SIA")
sing_snow("GRC")
sing_snow("DAV")
sing_snow("ANT")
sing_snow("SMM")
sing_snow("ABO")
sing_snow("DIS")
sing_snow("ELM")
sing_snow("SMA")

#Station names
cex_header <- 1.0
par(mar = c(0,0,0,0))

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("WFJ 2691 m",  side = 2, line = -2.6, cex = cex_header, adj = 0.968, outer = T)
mtext("ARO 1878 m",  side = 2, line = -2.6, cex = cex_header, adj = 0.872, outer = T)
mtext("SIA 1804 m",  side = 2, line = -2.6, cex = cex_header, adj = 0.777, outer = T)
mtext("GRC 1606 m",  side = 2, line = -2.6, cex = cex_header, adj = 0.682, outer = T)
mtext("DAV 1594 m",  side = 2, line = -2.6, cex = cex_header, adj = 0.586, outer = T)
mtext("ANT 1438 m",  side = 2, line = -2.6, cex = cex_header, adj = 0.492, outer = T)
mtext("SMM 1386 m",  side = 2, line = -2.6, cex = cex_header, adj = 0.401, outer = T)
mtext("ABO 1322 m",  side = 2, line = -2.6, cex = cex_header, adj = 0.302, outer = T)
mtext("DIS 1197 m",   side = 2, line = -2.6, cex = cex_header, adj = 0.211, outer = T)
mtext("ELM 957 m",   side = 2, line = -2.6, cex = cex_header, adj = 0.117, outer = T)
mtext("SMA 555 m",   side = 2, line = -2.6, cex = cex_header, adj = 0.025, outer = T)

#Analytical method

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
# mtext("a) SD [cm]", 
#       side = 3, line = -2.35, cex = cex_header+0.2, adj = 0.17)
# mtext("b) AMR [cm/day]", 
#       side = 3, line = -2.35, cex = cex_header+0.2, adj = 0.53)
# mtext("c) Sub-AMR [cm/day]",   
#       side = 3, line = -2.35, cex = cex_header+0.2, adj = 0.90)
# mtext("1959-1988",   
#       side = 3, line = -1.60, cex = cex_header-0.1, adj = 0.07)
# mtext("1989-2018",   
#       side = 3, line = -2.95, cex = cex_header-0.1, adj = 0.07)
# col_1 <- "steelblue4"
# col_2 <- "darkgoldenrod4"
# points(3.2, 10, pch = 19, col = col_2)
# points(3.2, 65, pch = 19, col = col_1)

dev.off()

#sfdfsdf----

