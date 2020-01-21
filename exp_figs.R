###

#Export figures for manusript melTim
#Erwin Rottler, University of Potsam

###

pacman::p_load(devtools, leaflet, raster, tmap, sf, prettymapr, meltimr, alptempr, rfs, viridis, 
               shape, scales, emdbook, zoo, zyp, maptools, rworldmap)

#set base direcoty
base_dir <- "U:/rhine_snow/"
data_dir <- "D:/nrc_user/rottler/"
grdc_dir <- "D:/nrc_user/rottler/GRDC_DAY/"

#Load results snow simulations
load(paste0(base_dir, "R/draft_snow_17_12.RData"))
load(paste0(base_dir, "R/draft_meltim.RData"))

#map_over----

dem = raster(paste0(data_dir, "basin_data/eu_dem/processed/eu_dem_500.tif"))

rivers <- rgdal::readOGR(dsn = paste0(data_dir, "basin_data/eu_hydro/wise_rivers_lakes/Large_rivers.shp"))
tribus <- rgdal::readOGR(dsn = paste0(data_dir, "basin_data/eu_hydro/wise_rivers_lakes/Other_large_rivers_and_tributaries.shp"))
lakes  <- rgdal::readOGR(dsn = paste0(data_dir, "basin_data/eu_hydro/wise_rivers_lakes/Large_lakes.shp"))

rhin_riv <- rivers[rivers@data$NAME == "Rhine",]
aare_riv <- tribus[tribus@data$NAME == "Aare",]
brienzersee <- lakes[which(lakes@data$WPLKNM == "BRIENZER SEE" ),]
thunersee <- lakes[which(lakes@data$WPLKNM == "THUNER SEE" ),]
bielersee <- lakes[which(lakes@data$WPLKNM == "BIELER SEE" ),]
bodensee <- lakes[which(lakes@data$WPLKNM == "BODENSEE" ),]

basins <-  rgdal::readOGR(dsn = paste0(data_dir, "basin_data/EZG_Schweiz_BAFU/ezg_kombiniert.shp"), encoding = "UTF8")

basin_base <- spTransform(basins[basins@data$Ort == "Basel, Rheinhalle",], CRS = crs(dem, asText = T))
basin_unte <- spTransform(basins[basins@data$Ort == "Untersiggenthal, Stilli",], CRS = crs(dem, asText = T))
basin_brug <- spTransform(basins[basins@data$Ort == "Brugg",], CRS = crs(dem, asText = T))
basin_bern <- spTransform(basins[basins@data$Ort == "Bern, SchÃ¶nau",], CRS = crs(dem, asText = T))
basin_neuh <- spTransform(basins[basins@data$Ort == "Neuhausen, FlurlingerbrÃ¼cke",], CRS = crs(dem, asText = T))
basin_reki <- spTransform(basins[basins@data$Ort == "Rekingen",], CRS = crs(dem, asText = T))
basin_diep <- spTransform(basins[basins@data$Ort == "Diepoldsau, RietbrÃ¼cke",], CRS = crs(dem, asText = T))

#7.6167, 47.5594 Basel
#8.33, 47.5704 Rekingen
#8.6259, 47.6823 Neuhausen
#8.2348, 47.5166 Untersiggenthal
#8.1949, 47.4825 Brugg
#7.448, 46.9331 Bern
#9.6409, 47.3831 Diepoldsau
crswgs84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
gauges_84 <-  sp::SpatialPoints(data.frame(lon = c(7.6167, 8.33, 8.6259, 8.2348, 8.1949, 7.448, 9.6409),
                                           lat = c(47.5594, 47.5704, 47.6823, 47.5166, 47.4825, 46.9331, 47.3831)),
                                proj4string =  crswgs84)
gauges    <- sp::spTransform(gauges_84, CRS = crs(basin_base, asText = T))
# basin_base_84 <- sp::spTransform(basin_base, CRS = crswgs84)

#corp DEM sub-basin area
my_ext <- extent(basin_base)
my_ext_buf <- my_ext + c(-20000, +30000, -30000, +12000) #xmin, xmax, ymin, ymax

my_box <- as(my_ext_buf, 'SpatialPolygons')
dem_cro_swiss <- raster::crop(dem, extent(my_box))
dem_sub_swiss <- mask(dem_cro_swiss, my_box)

#Snow stations

stat_id <- c("WFJ", "SAE", "ARO", "BIV", "SIA", "GRC", "DAV", "ANT", "MVE", "NAP", "SMM", "ABO", "SCU", "DIS", "CHM",
             "ROB", "GTT", "CDF", "CHD", "ELM", "EIN", "SMA", "BER", "BAS")
stat_na <- c("WeiÃŸfluhjoch", "Saentis", "Arosa", "Bivio", "Segl-Maria", "Graechen", "Davos", "Andermatt", "Montana", "Napf",
             "Sta. Maria, Val Muestair", "Adelboden", "Scuol", "Disentis", "Chaumant", "Poschiavo/Robbia",
             "Guttannen", "La Chaux-de-Fonds", "Chateaux-d'Oez", "Elm", "Einsiedeln", "Zuerich/Fluntern",
             "Bern, Zollikofen", "Basel Binningen")
stat_al <- c(2691, 2502, 1878, 1856, 1804, 1606, 1594, 1438, 1427, 1404, 1386, 1322, 1303, 1197, 1136, 1078, 1055, 1017,
             1028,  957, 910, 555, 552, 316)
stat_lo <- c(9.8000, 9.3500, 9.6833, 9.6666, 9.7666, 7.8333, 9.8500, 8.5833, 7.4666, 7.94, 10.4333, 7.5666, 10.2833, 8.8500, 6.9833,
             10.0666, 8.3000, 6.8000, 7.1333, 9.1833, 8.7500, 8.5666, 7.4666, 7.5836)
stat_la <- c(46.8333, 47.2500, 46.8000, 46.4666, 46.4333, 46.2000, 46.8166, 46.6333,  46.3000, 47.0047222, 46.6000, 46.5000, 46.8000, 46.7000,
             47.0500, 46.3500, 46.6500, 47.0833, 46.4833, 46.9166, 47.1333, 47.3833, 46.9833, 47.5411)
stat_nu <- 1:length(stat_id)

stat_meta <- data.frame(V0 = stat_nu,
                        V1 = stat_id,
                        V2 = stat_na,
                        V3 = stat_al,
                        V4 = stat_lo,
                        V5 = stat_la)
colnames(stat_meta) <- c("Number", "ID", "Name", "Alt. [m]", "Longitude", "Latitude")

#Select stations
stats_used <- c("WFJ", "ARO", "SIA", "GRC", "DAV", "ANT", "SMM", "ABO", "DIS", "ELM", "SMA", "NAP", "EIN")

stat_meta <- stat_meta[stat_meta$ID %in% stats_used, ]

snow_84 <-  sp::SpatialPoints(data.frame(lon = stat_meta$Longitude,
                                         lat = stat_meta$Latitude),
                              proj4string =  crswgs84)
snow    <- sp::spTransform(snow_84, CRS = crs(basin_base, asText = T))


pdf(paste0(base_dir,"R/figs_exp/map_over_raw.pdf"), width = 10, height = 6)

par(mar = c(0,0,0,0))
par(family = "serif")

plot(dem_sub_swiss, axes = F, legend = F,  col = colorRampPalette(c("white", "black"))(200), box = F)
plot(basin_bern, add =T, lwd = 0.025, col = alpha("steelblue4", alpha = 0.35))
plot(basin_brug, add =T, lwd = 0.025, col = alpha("steelblue4", alpha = 0.25))
plot(basin_unte, add =T, lwd = 0.025, col = alpha("steelblue4", alpha = 0.45))
plot(basin_diep, add =T, lwd = 0.025, col = alpha("darkolivegreen", alpha = 0.35))
plot(basin_neuh, add =T, lwd = 0.025, col = alpha("darkolivegreen", alpha = 0.25))
plot(basin_reki, add =T, lwd = 0.025, col = alpha("darkolivegreen", alpha = 0.45))
plot(basin_base, add =T, lwd = 1.8, border = "darkgoldenrod4")
plot(rhin_riv, add = T, col = "blue4")
plot(aare_riv, add = T, col = "blue4")
plot(brienzersee, add = T, col = "blue4")
plot(thunersee, add = T, col = "blue4")
plot(bielersee, add = T, col = "blue4")
plot(bodensee, add = T, col = "blue4")
plot(gauges, add = T, pch = 23, cex = 1.7,
     bg = alpha(c("darkgoldenrod4", "darkolivegreen", "darkolivegreen", "steelblue4",
                   "steelblue4", "steelblue4", "darkolivegreen"), alpha = 1.0))
plot(gauges, add = T, pch = 19, cex = 0.5)
plot(snow, add = T, pch = "*", cex = 1.3)
lab_mov <- 5200
lab_pos_1 <- c(lab_mov, -lab_mov, rep(lab_mov, 2), +lab_mov, rep(lab_mov, 7), +lab_mov)
lab_pos_2 <- c(lab_mov, -lab_mov, rep(lab_mov, 2), -lab_mov, rep(lab_mov, 7), -lab_mov)
text(snow@coords[, 1]+lab_pos_1, snow@coords[, 2]+lab_pos_2, labels = stat_meta$ID, col = "black", cex = 0.9)
addscalebar(plotunit = "m", widthhint = 0.2, htin = 0.15, pos = "topleft",
            padin = c(0.3, 0.5))

dev.off()


#map_ins----

dem = raster(paste0(data_dir, "basin_data/eu_dem/processed/eu_dem_1000.tif"))

pdf(paste0(base_dir,"R/figs_exp/map_ins_raw.pdf"), width = 10, height = 6)

#Create polygons for extends overview maps
boxes <- rbind(c(5.0, 40.0,  5.0, 55.0,  10.0, 55.0,  10.0, 40.0,  5.0, 40.0),
               c(my_ext@xmin, my_ext@ymin,  my_ext@xmin, my_ext@ymax,  my_ext@xmax, my_ext@ymax,  
                 my_ext@xmax, my_ext@ymin,  my_ext@xmin, my_ext@ymin))

matrix_1 <- matrix(boxes[1, ], ncol=2, byrow=TRUE)
matrix_2 <- matrix(boxes[2, ], ncol=2, byrow=TRUE)

ID <- c("box_1", "box_2")

polys <- SpatialPolygons(list(
  Polygons(list(Polygon(matrix_1)), ID[1]),
  Polygons(list(Polygon(matrix(boxes[2, ], ncol=2, byrow=TRUE))), ID[2])
))


plot(dem, col = colorRampPalette(c("grey85", "black"))(200), axes = F, legend = F, box = F)
plot(polys[2], col = alpha("white", alpha = 0.1), border = "red3", add = T, lwd = 3)

dev.off()


#elev_dist----

pdf(paste0(base_dir,"R/figs_exp/elev_dist.pdf"), width = 6, height = 2)

#Load DEM
dem = raster(paste0(data_dir, "basin_data/eu_dem/processed/eu_dem_500.tif"))

#Load basin boundaries
basins <-  rgdal::readOGR(dsn = paste0(data_dir, "basin_data/EZG_Schweiz_BAFU/ezg_kombiniert.shp"), encoding = "UTF8")
basin_base <- spTransform(basins[basins@data$Ort == "Basel, Rheinhalle",], CRS = crs(dem, asText = T))

#corp DEM sub-basin area
dem_cro <- raster::crop(dem, extent(basin_base))
dem_sub <- mask(dem_cro, basin_base)

#get elevations of cells cropped dem
dem_ele_NA <- dem_sub@data@values
dem_ele <- dem_ele_NA[!is.na(dem_ele_NA)]

hist_breaks <- seq(4000, 200, -50)

par(mar = c(1, 2.2, 2.0, 1))

par(family = "serif")

hist_res <- hist(dem_ele, breaks = hist_breaks, plot = F)

ylabs <- (hist_res$counts / length(dem_ele))  * 100

areal_perc <- c(round(length(which(dem_ele < 1000)) / length(dem_ele) * 100, digits = 2),
                round(length(which(dem_ele > 1000 & dem_ele < 2000)) / length(dem_ele) * 100, digits = 2),
                round(length(which(dem_ele > 2000 & dem_ele < 3000)) / length(dem_ele) * 100, digits = 2),
                round(length(which(dem_ele > 3000)) / length(dem_ele) * 100, digits = 2))

hist(dem_ele, breaks = hist_breaks, axes = F, ylab = "", xlab = "", main = "",
     col = alpha("darkgoldenrod4", alpha = 0.8), yaxs = "i", ylim = c(0,16000), xlim = c(4000, 0))
# rect(xleft =  6000, xright = -300, ytop = 150000, ybottom = 0 ,col = alpha("grey", alpha = 0.30), border = NA)
box()
axis(3, mgp=c(3, 0.08, 0), tck = -0.01, cex.axis = 0.8)
axis(2, mgp=c(3, 0.08, 0), tck = -0.01, cex.axis = 0.8,
     at = c(length(dem_ele) * c(0.01, 0.025, 0.05, 0.075, 0.1)), labels = c("1.0", "2.5", "5.0", "7.5", "10.0"))
abline(v = c(0, 1000, 2000, 3000, 4000), col = "black", lty = "dashed", lwd = 1.0)
mtext("Elevation [m]", side = 3, line = 1.1, adj = 0.5, cex = 1.2)
mtext("Areal fraction [%]", side = 2, line = 1.2, adj = 0.5, cex = 1.2)
mtext(paste(areal_perc[1], "%"), side = 2, line = -21.0, adj = 0.65, col = "black")
mtext(paste(areal_perc[2], "%"), side = 2, line = -17.0, adj = 0.65, col = "black")
mtext(paste(areal_perc[3], "%"), side = 2, line = -11.0, adj = 0.65, col = "black")
mtext(paste(areal_perc[4], "%"), side = 2, line = -4.5, adj = 0.65, col = "black")

dev.off()


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

pdf(paste0(base_dir,"R/figs_exp/snow_stats.pdf"), width = 12, height = 14.5)

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

#disc_perc----

grdc_data_base <- read_grdc(paste0(grdc_dir, "6935051_Q_Day.Cmd.txt"))
grdc_data_unte <- read_grdc(paste0(grdc_dir, "6935300_Q_Day.Cmd.txt"))
grdc_data_reki <- read_grdc(paste0(grdc_dir, "6935054_Q_Day.Cmd.txt"))
grdc_data_brug <- read_grdc(paste0(grdc_dir, "6935301_Q_Day.Cmd.txt"))
grdc_data_neuh <- read_grdc(paste0(grdc_dir, "6935055_Q_Day.Cmd.txt"))
grdc_data_bern <- read_grdc(paste0(grdc_dir, "6935020_Q_Day.Cmd.txt"))
grdc_data_diep <- read_grdc(paste0(grdc_dir, "6935500_Q_Day.Cmd.txt"))

pdf(paste0(base_dir,"R/figs_exp/disc_perc.pdf"), width = 12, height = 8)
# tiff("/home/rottler/ownCloud/RhineFlow/rhine_snow/manus/meltim_v1/figures/disc_perc.tiff", width = 12, height = 8,
#      units = "in", res = 800)

layout(matrix(c(rep(c(15, 3, 7, 11), 4),
                rep(c(1, 3, 7, 11), 3),
                rep(c(1, 4, 8, 12), 1),
                rep(c(1, 5, 9, 13), 3),
                rep(c(2, 5, 9, 13), 1),
                rep(c(16, 5, 9, 13), 3),
                rep(c(16, 6, 10, 14), 1)
),
4, 16), widths=c(), heights=c())

# layout.show(n = 16)

sta_yea_cla <- 1927
yea_cla_1 <- sta_yea_cla
yea_cla_2 <- 1966
yea_cla_3 <- 1977
yea_cla_4 <- 2016

perce_plot <- function(data_in, date_in, main_in = "", year_1 = yea_cla_1, year_2 = yea_cla_2,
                       year_3 = yea_cla_3, year_4 = yea_cla_4){
  
  data_day <- ord_day(data_in = data_in,
                      date = date_in,
                      start_y = year_1,
                      end_y = year_4)
  
  #Mean seasonal cycles
  ind_break_1 <- (year_1 - year_1) + 1
  ind_break_2 <- (year_2 - year_1) + 1
  ind_break_3 <- (year_3 - year_1) + 1
  ind_break_4 <- (year_4 - year_1) + 1
  
  #Calculation percentile graph
  jan_cols <- 1:31
  feb_cols <- 32:59
  mar_cols <- 60:90
  apr_cols <- 91:120
  may_cols <- 121:151
  jun_cols <- 152:181
  jul_cols <- 182:212
  aug_cols <- 213:243
  sep_cols <- 244:273
  oct_cols <- 274:304
  nov_cols <- 305:334
  dec_cols <- 335:365
  
  month_cols <- list(jan_cols, feb_cols, mar_cols, apr_cols, may_cols, jun_cols, jul_cols, aug_cols, sep_cols, oct_cols, nov_cols, dec_cols)
  
  f_quants_all <- function(data_in){
    
    probs <- seq(0.01, 0.99, 0.01)
    quants_all <- quantile(data_in,  probs = probs, type = 8, na.rm = T)
    
    return(quants_all)
    
  }
  
  qmon_1 <- matrix(NA, ncol = 12, nrow = length(seq(0.01, 0.99, 0.01)))
  qmon_2 <-  matrix(NA, ncol = 12, nrow = length(seq(0.01, 0.99, 0.01)))
  
  for(i in 1:12){
    qmon_1[ , i] <- f_quants_all(data_day[ind_break_1:ind_break_2, month_cols[[i]]])
    qmon_2[ , i] <- f_quants_all(data_day[ind_break_3:ind_break_4, month_cols[[i]]])
  }
  
  qdif <- qmon_2 - qmon_1
  
  
  #Percentile graph
  x_axis_lab <- 1:12
  x_axis_tic <- (1:13)-0.5
  
  n_max <- round(abs(max_na(qdif[, ])) / (max_na(qdif[, ]) + abs(min_na(qdif[, ]))), digits = 2) * 200
  n_min <- 200 - n_max
  # n_max <- 100
  # n_min <- 100

  
  # cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[c(1,1,2,3,4)], "grey98"))(n_min)
  # cols_max <- colorRampPalette(c("grey98", "gold3",  "orangered4", "firebrick4", "firebrick4", "darkred"))(n_max)
  
  cols_min <- colorRampPalette(c("darkred", "firebrick4", "orangered4", "gold3", "grey98"))(n_min)
  cols_max <- colorRampPalette(c("white", "azure3", viridis::viridis(9, direction = 1)[c(4,3,2,1)]))(n_max)
  my_col <- colorRampPalette(c(cols_min, cols_max))(200)
  # my_bre <- seq(min_na(qdif[, ]), max_na(qdif[, ]), length.out = 201)
  my_bre <- c(seq(-max_na(abs(qdif)), 0, length.out = n_min),
              seq(0, max_na(abs(qdif)), length.out = n_max+1))

  # plot(1:n_max, 1:n_max, pch = 19, cex = 2, col = cols_max)
  # plot(1:n_min, 1:n_min, pch = 19, cex = 2, col = cols_min)
  
  par(mar = c(1.6, 3.5, 2.0, 0.2))
  par(family = "serif")
  
  image(x = 1:12,
        y = seq(0.01, 0.99, 0.01),
        z = t(qdif), col = my_col, breaks = my_bre,
        ylab = "", xlab = "", axes = F)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.05)#plot ticks
  axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S","O", "N", "D"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.5)#plot labels
  axis(2, mgp=c(3, 0.25, 0), tck = -0.025, cex.axis = 1.5)
  mtext("Prob. level", side = 2, line = 1.7, cex = 1.25)
  mtext(main_in, side = 3, line = 0.25, cex = 1.5, adj = 0.0)
  mtext("[m³/s]", side = 3, line = 0.2, cex = 1.2, adj = 1.0)
  
  box()
  
  par(mar = c(2.0, 0.2, 2.5, 3.5))
  
  alptempr::image_scale(as.matrix(qdif), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
  axis(4, mgp=c(3, 0.45, 0), tck = -0.1, cex.axis = 1.5)
  box()
  
}

perce_plot(data_in = grdc_data_base$value,
           date_in = grdc_data_base$date,
           main_in = "a) Basel")

perce_plot(data_in = grdc_data_unte$value,
           date_in = grdc_data_unte$date,
           main_in = "b) Untersiggenthal")

perce_plot(data_in = grdc_data_reki$value,
           date_in = grdc_data_reki$date,
           main_in = "c) Rekingen")

perce_plot(data_in = grdc_data_brug$value,
           date_in = grdc_data_brug$date,
           main_in = "d) Brugg")

perce_plot(data_in = grdc_data_neuh$value,
           date_in = grdc_data_neuh$date,
           main_in = "e) Neuhausen")

perce_plot(data_in = grdc_data_bern$value,
           date_in = grdc_data_bern$date,
           main_in = "f) Bern")

perce_plot(data_in = grdc_data_diep$value,
           date_in = grdc_data_diep$date,
           main_in = "g) Diepoldsau")

dev.off()


#disc_rast----

pdf(paste0(base_dir,"R/figs_exp/disc_rast.pdf"), width = 12, height = 8)

layout(matrix(c(rep(c(15, 3, 7, 11), 4),
                rep(c(1, 3, 7, 11), 3),
                rep(c(1, 4, 8, 12), 1),
                rep(c(1, 5, 9, 13), 3),
                rep(c(2, 5, 9, 13), 1),
                rep(c(16, 5, 9, 13), 3),
                rep(c(16, 6, 10, 14), 1)
),
4, 16), widths=c(), heights=c())

sta_yea_cla <- 1927
yea_cla_1 <- sta_yea_cla
yea_cla_2 <- 1966
yea_cla_3 <- 1977
yea_cla_4 <- 2016


raster_plot <- function(data_in, date_in, main_in = "", year_1  = yea_cla_1, year_2 = yea_cla_4,
                        quant_break = 0.50){
  
  data_day <- ord_day(data_in = data_in,
                      date = date_in,
                      start_y = year_1,
                      end_y = year_2)
  
  #Preparation Raster graph
  x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
  x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
  
  # cols_min <- grDevices::colorRampPalette(c(viridis::viridis(9, direction = 1)[c(1,1, 2:4)], "grey98"))(100)
  # cols_max <- grDevices::colorRampPalette(c("grey98", "gold3", "orangered4", "darkred"))(100)
  cols_min <- colorRampPalette(c("darkred", "orangered4", "goldenrod3", "gold3", "white"))(100)
  cols_max <- colorRampPalette(c("white", "azure3", viridis::viridis(9, direction = 1)[c(4,3,2,1,1)]))(100)
  cols_hydro <- c(cols_min, rep("white", 16), cols_max)
  # cols_hydro <- grDevices::colorRampPalette(c(viridis::viridis(9, direction = 1)))(100)
  # plot(1:100, 1:100, pch = 19, col = cols_min)
  
  max_break <- max_na(data_day)
  min_break <- min_na(data_day)
  qua_break <- quantile(data_day, probs = quant_break, type = 8, na.rm = T)
  
  breaks_1 <- lseq(min_break, qua_break, length.out = length(cols_hydro)/2)
  breaks_2 <- lseq(qua_break+0.01, max_break, length.out = length(cols_hydro)/2 + 1)
  breaks_2[length(breaks_2)] <- breaks_2[length(breaks_2)] + 0.1
  breaks_hydro <- c(breaks_1, breaks_2)
  
  #Raster graph
  par(mar = c(1.6, 3.5, 2.0, 0.2))
  par(family = "serif")
  
  image(x = 1:ncol(data_day),
        y = 1927:2016,
        z = t(data_day),
        col = cols_hydro,
        breaks = breaks_hydro,
        ylab = "", xlab = "", axes = F)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.05)#plot ticks
  axis(1, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S","O", "N", "D"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.5)#plot labels
  axis(2, mgp=c(3, 0.15, 0), tck = -0.025, cex.axis = 1.5)
  mtext("Year", side = 2, line = 1.7, cex = 1.25)
  mtext(main_in, side = 3, line = 0.25, cex = 1.5, adj = 0.0)
  mtext("[m³/s]", side = 3, line = 0.2, cex = 1.2, adj = 1.0)
  
  box()
  
  par(mar = c(2.0, 0.2, 2.0, 3.5))
  
  alptempr::image_scale(as.matrix(data_day), col = cols_hydro, breaks = breaks_hydro, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
  axis(4, mgp=c(3, 0.45, 0), tck = -0.1, cex.axis = 1.5)
  box()
  
  
  
  
  
}

raster_plot(data_in = grdc_data_base$value,
            date_in = grdc_data_base$date,
            main_in = "a) Basel")

raster_plot(data_in = grdc_data_unte$value,
            date_in = grdc_data_unte$date,
            main_in = "b) Untersiggenthal")

raster_plot(data_in = grdc_data_reki$value,
            date_in = grdc_data_reki$date,
            main_in = "c) Rekingen")

raster_plot(data_in = grdc_data_brug$value,
            date_in = grdc_data_brug$date,
            main_in = "d) Brugg")

raster_plot(data_in = grdc_data_neuh$value,
            date_in = grdc_data_neuh$date,
            main_in = "e) Neuhausen")

raster_plot(data_in = grdc_data_bern$value,
            date_in = grdc_data_bern$date,
            main_in = "f) Bern")

raster_plot(data_in = grdc_data_diep$value,
            date_in = grdc_data_diep$date,
            main_in = "g) Diepoldsau")

dev.off()














#comp_illu----

pdf(paste0(base_dir,"R/figs_exp/comp_illu_raw.pdf"), width = 7, height = 3)

par(mar = c(1.3, 2.5, 1.2, 0.2))
par(family = "serif")

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
col_1 <- alpha("steelblue4", alpha = 0.8)
col_2 <- alpha("darkred", alpha = 0.8)

plot(1:10, 1:10, xlim = c(1, 365), ylim = c(500, 3200), type = "n", axes = F, ylab = "", xlab = "")
#Points Present
points(c(142, 142, 143, 144, 145, 146), c(500, 600, 700, 800, 900, 1000),
       col = col_1, pch = 19, cex = 1.2)
points(c(170, 172, 174, 176, 178, 179, 180, 181), c(1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800),
       col = col_1, pch = 19, cex = 1.2)
points(c(214, 215, 216, 217, 219, 221, 224, 230, 235, 240), c(1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800),
       col = col_1, pch = 19, cex = 1.2)
points(c(256, 260, 265, 272), c(2900, 3000, 3100, 3200),
       col = col_1, pch = 19, cex = 1.2)
#Points Future
points(c(142, 142, 143, 144, 145, 146), c(900, 1000, 1100, 1200, 1300, 1400),
       col = col_2, pch = 19, cex = 1.2)
points(c(170, 172, 174, 176, 178, 179, 180, 181), c(1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200),
       col = col_2, pch = 19, cex = 1.2)
points(c(214, 215, 216, 217, 219, 221, 224, 230, 235, 240), c(2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200),
       col = col_2, pch = 19, cex = 1.2)
abline(h = c(500, 1000, 1500, 2000, 2500, 3000), lty = "dotted", col = "grey55", lwd = 0.8)
abline(v = x_axis_tic, lty = "dotted", col = "grey55", lwd = 0.8)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.04)#plot ticks
axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.10, 0), cex.axis = 0.9)#plot labels
axis(2, mgp=c(3, 0.15, 0), tck = -0.01, cex.axis = 1.0)
mtext("Elevation [m]", side = 2, line = 1.3, cex = 1.0, adj = 0.5)
mtext("b) Snow melt elevation compensation", side = 3, line = 0.15, cex = 1.2, adj = 0.0)
legend("topleft", c("future", "present"), col = c(col_2, col_1), pch = 19, bg = "white")
Arrows(x0 = 166, x1 =156, y0 = 1250, y1 = 1250, arr.type = "triangle")
Arrows(x0 = 205, x1 =192, y0 = 2000, y1 = 2000, arr.type = "triangle")
Arrows(x0 = 254, x1 =243, y0 = 3000, y1 = 3000, arr.type = "triangle")
Arrows(x0 = 130, x1 =130, y0 = 750,  y1 = 1100, arr.type = "triangle", col = "black")
Arrows(x0 = 162, x1 =162, y0 = 1650, y1 = 2000, arr.type = "triangle", col = "black")
Arrows(x0 = 208, x1 =208, y0 = 2650, y1 = 3000, arr.type = "triangle", col = "black")
box()



dev.off()


#snow_simu----

pdf(paste0(base_dir,"R/figs_exp/snow_simu.pdf"), width = 12, height = 6)

layout(matrix(c(rep(c(1, 5, 9), 8), 2, 6, 10,
                rep(c(3, 7, 11), 8), 4, 8, 12),
              3, 18), widths=c(), heights=c())
# layout.show(n=12)

snow_sim_plot <- function(data_plot, cols, breaks, header, lab_unit){
  
  par(mar = c(2.2, 3.5, 2.5, 0.2))
  par(family = "serif")
  
  x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
  x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
  
  image(x = 1:365,
        y = elev_bands[-length(elev_bands)],
        z = data_plot, col = cols, breaks = breaks,
        ylab = "", xlab = "", axes = F)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.04)#plot ticks
  axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.6)#plot labels
  axis(2, mgp=c(3, 0.25, 0), tck = -0.005, cex.axis = 1.6)
  mtext("Elevation", side = 2, line = 1.8, cex = 1.3)
  mtext(header, side = 3, line = 0.3, cex = 1.5, adj = 0.0)
  mtext(lab_unit, side = 3, line = 0.2, cex = 1.2, adj = 1.0)
  box()
  
  par(mar = c(2.2, 0.2, 2.5, 3.0))
  
  alptempr::image_scale(as.matrix(data_plot), col = cols, breaks = breaks, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
  axis(4, mgp=c(3, 0.35, 0), tck = -0.08, cex.axis = 1.4)
  box()
  
  
}

#SWE depth mean
my_col <- colorRampPalette(c("grey95", viridis::viridis(9, direction = 1)[4:1]))(200)
my_bre <- seq(alptempr::min_na(smea_band), alptempr::max_na(smea_band), length.out = length(my_col)+1)

snow_sim_plot(smea_band, cols = my_col, breaks = my_bre,
              header = "a) SWE depth mean", lab_unit = "[m]")

#SWE depth trend
cols_min <- colorRampPalette(c("darkred", "firebrick4", "firebrick4", "orange3", "darkgoldenrod3", "grey98"))(100)
cols_max <- colorRampPalette(c("grey98", viridis::viridis(9, direction = 1)[4:1]))(100)

my_col <- alpha(c(cols_min, cols_max), alpha = 1.0)
my_bre <- seq(-max_na(abs(sslo_band)), max_na(abs(sslo_band)), length.out = length(my_col)+1)

snow_sim_plot(sslo_band, cols = my_col, breaks = my_bre,
              header = "a) SWE depth mean", lab_unit = "[m/dec]")

#SWE volume mean
my_col <- colorRampPalette(c("grey95", viridis::viridis(9, direction = 1)[4:1]))(200)
my_bre <- seq(alptempr::min_na(vmea_band), alptempr::max_na(vmea_band), length.out = length(my_col)+1)

snow_sim_plot(vmea_band, cols = my_col, breaks = my_bre,
              header = "a) SWE volume mean", lab_unit = "[hm³]")

#SWE volume trend
cols_min <- colorRampPalette(c("darkred", "firebrick4", "orange3", "darkgoldenrod3", "grey98"))(100)
cols_max <- colorRampPalette(c("grey98", viridis::viridis(9, direction = 1)[4:1]))(100)

my_col <- c(cols_min, cols_max)
my_bre <- seq(-max_na(abs(vslo_band)), max_na(abs(vslo_band)), length.out = length(my_col)+1)

snow_sim_plot(vslo_band, cols = my_col, breaks = my_bre,
              header = "a) SWE depth mean", lab_unit = "[hm³/dec]")

#SWE volume diff mean
cols_max <- colorRampPalette(c("grey98", "darkgoldenrod3", "orange3", "firebrick4", "darkred"))(100)
cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "grey98"))(100)
my_col <- c(cols_min, cols_max)
my_bre <- seq(alptempr::min_na(vdif_band), alptempr::max_na(vdif_band), length.out = length(my_col)+1)

snow_sim_plot(vdif_band, cols = my_col, breaks = my_bre,
              header = "a) Accum./Melt mean", lab_unit = "[hm³/dec]")



n_max <- 100
n_min <- 100

cols_max <- colorRampPalette(c("grey98", "darkgoldenrod3", "orange3", "firebrick4", "darkred"))(n_max)
cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "grey98"))(n_min)
my_col <- alpha(c(cols_min, cols_max), alpha = 1.0)



# #SWE mean
# data_plot <- smea_band; col_zero <- F; lab_unit <- "[m]"
# 
# par(mar = c(2.2, 3.5, 2.5, 0.2))
# par(family = "serif")
# 
# x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
# x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
# 
# plot_test <- data_plot
# 
# my_col <- colorRampPalette(c("grey95", viridis::viridis(9, direction = 1)[4:1]))(200)
# my_col <- alpha(my_col, alpha = 1.0)
# 
# my_bre <- seq(alptempr::min_na(plot_test), alptempr::max_na(plot_test), length.out = length(my_col)+1)
# 
# image(x = 1:365,
#       y = elev_bands[-length(elev_bands)],
#       z = plot_test, col =my_col, breaks = my_bre,
#       ylab = "", xlab = "", axes = F)
# axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
#      col = "black", col.axis = "black", tck = -0.04)#plot ticks
# axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
#      col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.6)#plot labels
# axis(2, mgp=c(3, 0.25, 0), tck = -0.005, cex.axis = 1.6)
# mtext("Elevation", side = 2, line = 1.8, cex = 1.3)
# mtext("a) SWE depth mean", side = 3, line = 0.3, cex = 1.5, adj = 0.0)
# mtext("[m]", side = 3, line = 0.2, cex = 1.2, adj = 1.0)
# box()
# 
# par(mar = c(2.2, 0.2, 2.5, 3.0))
# 
# alptempr::image_scale(as.matrix(plot_test), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
# axis(4, mgp=c(3, 0.35, 0), tck = -0.08, cex.axis = 1.4)
# # mtext(lab_unit, side = 3, line = 0.3, cex = 1)
# box()


# #SWE slo
# data_plot <- sslo_band; col_zero <- T; lab_unit <- "[m/dec]"
# 
# par(mar = c(2.2, 3.5, 2.5, 0.2))
# 
# x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
# x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
# 
# plot_test <- data_plot
# 
# my_col <- c(colorRampPalette(c("white", viridis::viridis(20, direction = -1)))(200))
# 
# my_bre <- seq(alptempr::min_na(plot_test), alptempr::max_na(plot_test), length.out = length(my_col)+1)
# 
# if(col_zero){
#   
#   # n_max <- round(abs(alptempr::max_na(plot_test[, ])) / (alptempr::max_na(plot_test[, ]) + abs(alptempr::min_na(plot_test[, ]))), digits = 2) * 200
#   # n_min <- 200 - n_max
#   n_max <- 100
#   n_min <- 100
#   
#   # cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "cadetblue3", "white"))(n_min)
#   # cols_min <- colorRampPalette(c("darkred","orangered4", "orangered4", "orange3", "gold3", "grey90"))(n_min)
#   # cols_min <- colorRampPalette(c("darkred","orangered4", "orangered4", "orange3", "grey92"))(n_min)
#   # cols_max <- colorRampPalette(c("grey80", "cadetblue3"))(n_max)
#   # cols_max <- colorRampPalette(c("grey92", viridis::viridis(9, direction = 1)[4:1]))(n_max)
#   # cols_max <- colorRampPalette(c("grey90", viridis::viridis(9, direction = 1)[4:1]))(n_max)
#   cols_min <- colorRampPalette(c("darkred", "firebrick4", "firebrick4", "orange3", "darkgoldenrod3", "grey98"))(n_min)
#   # cols_max <- colorRampPalette(c("grey80", "cadetblue3"))(n_max)
#   cols_max <- colorRampPalette(c("grey98", viridis::viridis(9, direction = 1)[4:1]))(n_max)
#   
#   my_col <- alpha(c(cols_min, cols_max), alpha = 1.0)
#   
#   my_bre <- seq(-max_na(abs(plot_test)), max_na(abs(plot_test)), length.out = length(my_col)+1)
#   
# }
# 
# image(x = 1:365,
#       y = elev_bands[-length(elev_bands)],
#       z = plot_test, col =my_col, breaks = my_bre,
#       ylab = "", xlab = "", axes = F)
# axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
#      col = "black", col.axis = "black", tck = -0.04)#plot ticks
# axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
#      col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.6)#plot labels
# axis(2, mgp=c(3, 0.25, 0), tck = -0.005, cex.axis = 1.6)
# mtext("Elevation", side = 2, line = 1.8, cex = 1.3)
# mtext("b) SWE depth trend", side = 3, line = 0.3, cex = 1.5, adj = 0.0)
# mtext("[m/dec]", side = 3, line = 0.2, cex = 1.2, adj = 1.0)
# box()
# 
# par(mar = c(2.2, 0.2, 2.5, 3.0))
# 
# alptempr::image_scale(as.matrix(plot_test), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
# axis(4, mgp=c(3, 0.35, 0), tck = -0.08, cex.axis = 1.4)
# # mtext(lab_unit, side = 3, line = 0.3, cex = 1)
# box()


# #SWE volume mean
# data_plot <- vmea_band/1000000; col_zero <- F; lab_unit <- "[m³]"
# 
# par(mar = c(2.2, 3.5, 2.5, 0.2))
# 
# x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
# x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
# 
# plot_test <- data_plot
# 
# my_col <- colorRampPalette(c("grey95", viridis::viridis(9, direction = 1)[4:1]))(200)
# my_col <- alpha(my_col, alpha = 1.0)
# 
# my_bre <- seq(alptempr::min_na(plot_test), alptempr::max_na(plot_test), length.out = length(my_col)+1)
# 
# image(x = 1:365,
#       y = elev_bands[-length(elev_bands)],
#       z = plot_test, col =my_col, breaks = my_bre,
#       ylab = "", xlab = "", axes = F)
# axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
#      col = "black", col.axis = "black", tck = -0.04)#plot ticks
# axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
#      col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.6)#plot labels
# axis(2, mgp=c(3, 0.25, 0), tck = -0.005, cex.axis = 1.6)
# mtext("Elevation", side = 2, line = 1.8, cex = 1.3)
# mtext("c) SWE volume", side = 3, line = 0.3, cex = 1.5, adj = 0.0)
# mtext("[hm³]", side = 3, line = 0.2, cex = 1.2, adj = 1.0)
# box()
# 
# par(mar = c(2.2, 0.2, 2.5, 3.0))
# 
# alptempr::image_scale(as.matrix(plot_test), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
# axis(4, mgp=c(3, 0.35, 0), tck = -0.08, cex.axis = 1.4)
# box()


# #SWE volume slo
# data_plot <- vslo_band/1000000; col_zero <- T; lab_unit <- "[hm³/dec]"
# 
# par(mar = c(2.2, 3.5, 2.5, 0.2))
# 
# x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
# x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
# 
# plot_test <- data_plot
# 
# my_col <- c(colorRampPalette(c("white", viridis::viridis(20, direction = -1)))(200))
# 
# my_bre <- seq(alptempr::min_na(plot_test), alptempr::max_na(plot_test), length.out = length(my_col)+1)
# 
# if(col_zero){
#   
#   # n_max <- round(abs(alptempr::max_na(plot_test[, ])) / (alptempr::max_na(plot_test[, ]) + abs(alptempr::min_na(plot_test[, ]))), digits = 2) * 200
#   # n_min <- 200 - n_max
#   n_max <- 100
#   n_min <- 100
#   cols_min <- colorRampPalette(c("darkred", "firebrick4", "orange3", "darkgoldenrod3", "grey98"))(n_min)
#   # cols_max <- colorRampPalette(c("grey80", "cadetblue3"))(n_max)
#   cols_max <- colorRampPalette(c("grey98", viridis::viridis(9, direction = 1)[4:1]))(n_max)
#   
#   my_col <- alpha(c(cols_min, cols_max), alpha = 1.0)
#   
#   my_bre <- seq(-max_na(abs(plot_test)), max_na(abs(plot_test)), length.out = length(my_col)+1)
#   
# }
# 
# image(x = 1:365,
#       y = elev_bands[-length(elev_bands)],
#       z = plot_test, col =my_col, breaks = my_bre,
#       ylab = "", xlab = "", axes = F)
# axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
#      col = "black", col.axis = "black", tck = -0.04)#plot ticks
# axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
#      col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.6)#plot labels
# axis(2, mgp=c(3, 0.25, 0), tck = -0.005, cex.axis = 1.6)
# mtext("Elevation", side = 2, line = 1.8, cex = 1.3)
# mtext("d) SWE volume trend", side = 3, line = 0.3, cex = 1.5, adj = 0.0)
# mtext("[hm³/dec]", side = 3, line = 0.2, cex = 1.2, adj = 1.0)
# box()
# 
# par(mar = c(2.2, 0.2, 2.5, 3.0))
# 
# alptempr::image_scale(as.matrix(plot_test), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
# axis(4, mgp=c(3, 0.35, 0), tck = -0.08, cex.axis = 1.4)
# box()


#Snow diff mean
data_plot <- vdif_band/1000000; col_zero <- T; lab_unit <- "[hm³]"

par(mar = c(2.2, 3.5, 2.5, 0.2))

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15

plot_test <- data_plot

my_col <- c(colorRampPalette(c("white", viridis::viridis(20, direction = -1)))(200))

my_bre <- seq(alptempr::min_na(plot_test), alptempr::max_na(plot_test), length.out = length(my_col)+1)

if(col_zero){
  
  # n_max <- round(abs(alptempr::max_na(plot_test[, ])) / (alptempr::max_na(plot_test[, ]) + abs(alptempr::min_na(plot_test[, ]))), digits = 2) * 200
  # n_min <- 200 - n_max
  n_max <- 100
  n_min <- 100

  cols_max <- colorRampPalette(c("grey98", "darkgoldenrod3", "orange3", "firebrick4", "darkred"))(n_max)
  cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "grey98"))(n_min)
  my_col <- alpha(c(cols_min, cols_max), alpha = 1.0)
  
  my_bre <- seq(-max_na(abs(plot_test)), max_na(abs(plot_test)), length.out = length(my_col)+1)
  
}

image(x = 1:365,
      y = elev_bands[-length(elev_bands)],
      z = plot_test, col =my_col, breaks = my_bre,
      ylab = "", xlab = "", axes = F)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.04)#plot ticks
axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.6)#plot labels
axis(2, mgp=c(3, 0.25, 0), tck = -0.005, cex.axis = 1.6)
mtext("Elevation", side = 2, line = 1.8, cex = 1.3)
mtext("e) Accum./Melt mean", side = 3, line = 0.3, cex = 1.5, adj = 0.0)
mtext("[hm³]", side = 3, line = 0.2, cex = 1.2, adj = 1.0)
box()

par(mar = c(2.2, 0.2, 2.5, 3.0))

alptempr::image_scale(as.matrix(plot_test), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.35, 0), tck = -0.08, cex.axis = 1.4)
box()


#Snow diff slo
data_plot <- vdis_band/1000000; col_zero <- T; lab_unit <- "[hm³/dec]"

par(mar = c(2.2, 3.5, 2.5, 0.2))

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15

plot_test <- data_plot

my_col <- c(colorRampPalette(c("white", viridis::viridis(20, direction = -1)))(200))

my_bre <- seq(alptempr::min_na(plot_test), alptempr::max_na(plot_test), length.out = length(my_col)+1)

if(col_zero){
  
  # n_max <- round(abs(alptempr::max_na(plot_test[, ])) / (alptempr::max_na(plot_test[, ]) + abs(alptempr::min_na(plot_test[, ]))), digits = 2) * 200
  # n_min <- 200 - n_max
  n_max <- 100
  n_min <- 100

  cols_max <- colorRampPalette(c("grey98", "darkgoldenrod3", "orange3", "firebrick4", "firebrick4", "darkred", "darkred"))(n_max)
  cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "grey98"))(n_min)
  my_col <- alpha(c(cols_min, cols_max), alpha = 1.0)
  
  my_bre <- seq(-max_na(abs(plot_test)), max_na(abs(plot_test)), length.out = length(my_col)+1)
  
}

image(x = 1:365,
      y = elev_bands[-length(elev_bands)],
      z = plot_test, col =my_col, breaks = my_bre,
      ylab = "", xlab = "", axes = F)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.04)#plot ticks
axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.6)#plot labels
axis(2, mgp=c(3, 0.25, 0), tck = -0.005, cex.axis = 1.6)
mtext("Elevation", side = 2, line = 1.8, cex = 1.3)
mtext("f) Accum./Melt trend", side = 3, line = 0.3, cex = 1.5, adj = 0.0)
mtext("[hm³/dec]", side = 3, line = 0.2, cex = 1.2, adj = 1.0)
box()

par(mar = c(2.2, 0.2, 2.5, 3.0))

alptempr::image_scale(as.matrix(plot_test), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.35, 0), tck = -0.08, cex.axis = 1.4)
box()


dev.off()

