###

#Export figures for manusript melTim
#Erwin Rottler, University of Potsdam

###

pacman::p_load(devtools, leaflet, raster, tmap, sf, prettymapr, meltimr, alptempr, rfs, viridis, 
               shape, scales, emdbook, zoo, zyp, maptools, rworldmap, parallel, doParallel,
               DescTools, POT)

#set base direcoty
base_dir <- "U:/rhine_snow/"
data_dir <- "D:/nrc_user/rottler/"
grdc_dir <- "D:/nrc_user/rottler/GRDC_DAY/"

#Load results snow simulations
load("U:/rhine_snow/R/figs_exp/sim_scf_exp.RData")
load("U:/rhine_snow/R/figs_exp/temp_prec_exp.RData")

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
basin_bern <- spTransform(basins[basins@data$Ort == "Bern, Schönau",], CRS = crs(dem, asText = T))
basin_neuh <- spTransform(basins[basins@data$Ort == "Neuhausen, Flurlingerbrücke",], CRS = crs(dem, asText = T))
basin_reki <- spTransform(basins[basins@data$Ort == "Rekingen",], CRS = crs(dem, asText = T))
basin_diep <- spTransform(basins[basins@data$Ort == "Diepoldsau, Rietbrücke",], CRS = crs(dem, asText = T))

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
my_ext_buf <- my_ext + c(-20000, +30000, -40000, +13000) #xmin, xmax, ymin, ymax

my_box <- as(my_ext_buf, 'SpatialPolygons')
dem_cro_swiss <- raster::crop(dem, extent(my_box))
dem_sub_swiss <- mask(dem_cro_swiss, my_box)

#Meteo stations

stat_id <- c("WFJ", "SAE", "ARO", "BIV", "SIA", "GRC", "DAV", "ANT", "MVE", "NAP", "SMM", "ABO", "SCU", "DIS", "CHM",
             "ROB", "GTT", "CDF", "CHD", "ELM", "EIN", "SMA", "BER", "BAS", "ZER", "5WJ", "4ZE", "5DF", "2AN")
stat_na <- c("Weißfluhjoch", "Saentis", "Arosa", "Bivio", "Segl-Maria", "Graechen", "Davos", "Andermatt", "Montana", "Napf",
             "Sta. Maria, Val Muestair", "Adelboden", "Scuol", "Disentis", "Chaumant", "Poschiavo/Robbia",
             "Guttannen", "La Chaux-de-Fonds", "Chateaux-d'Oez", "Elm", "Einsiedeln", "Zuerich/Fluntern",
             "Bern, Zollikofen", "Basel Binningen", "Zermatt", "Weissfluhjoch", "Zermatt", "Davos Flueestr.", "Andermatt")
stat_al <- c(2691, 2502, 1878, 1856, 1804, 1606, 1594, 1438, 1427, 1404, 1386, 1322, 1303, 1197, 1136, 1078, 1055, 1017,
             1028,  957, 910, 555, 552, 316, 1638, 2540, 1600, 1560, 1440)
stat_lo <- c(9.8000, 9.3500, 9.6833, 9.6666, 9.7666, 7.8333, 9.8500, 8.5833, 7.4666, 7.94, 10.4333, 7.5666, 10.2833, 8.8500, 6.9833,
             10.0666, 8.3000, 6.8000, 7.1333, 9.1833, 8.7500, 8.5666, 7.4666, 7.5836, 7.7531, 9.8093, 7.7512, 9.8482, 8.5919)
stat_la <- c(46.8333, 47.2500, 46.8000, 46.4666, 46.4333, 46.2000, 46.8166, 46.6333,  46.3000, 47.0047222, 46.6000, 46.5000, 46.8000, 46.7000,
             47.0500, 46.3500, 46.6500, 47.0833, 46.4833, 46.9166, 47.1333, 47.3833, 46.9833, 47.5411, 46.0291667, 46.8294, 46.0234,
             46.8125, 46.6329)
stat_nu <- 1:length(stat_id)

stat_meta <- data.frame(V0 = stat_nu,
                        V1 = stat_id,
                        V2 = stat_na,
                        V3 = stat_al,
                        V4 = stat_lo,
                        V5 = stat_la)
colnames(stat_meta) <- c("Number", "ID", "Name", "Alt. [m]", "Longitude", "Latitude")

#Select stations
stats_used <- c("WFJ", "ARO", "SIA", "GRC", "DAV", "ANT", "SMM", "ABO", "DIS", "ELM", "SMA", "NAP", "EIN", "ZER",
                "BER", "BAS", "5WJ", "4ZE", "5DF", "2AN")

stat_meta <- stat_meta[stat_meta$ID %in% stats_used, ]

snow_84 <-  sp::SpatialPoints(data.frame(lon = stat_meta$Longitude,
                                         lat = stat_meta$Latitude),
                              proj4string =  crswgs84)
snow    <- sp::spTransform(snow_84, CRS = crs(basin_base, asText = T))


pdf(paste0(base_dir,"R/figs_exp/map_over_raw.pdf"), width = 10, height = 6)

par(bg=NA, mar=c(0,0,0,0), oma=c(0,0,0,0))
par(family = "serif")

col_elev <- colorRampPalette(c("white", "black"))(200)
# col_elev <- viridis::viridis(n = 200, direction = -1)
# col_elev <- terrain.colors(200)
# cols_min <- colorRampPalette(c("darkred", "firebrick4", "orangered4", "gold3", "grey98"))(n_min)
# cols_max <- colorRampPalette(c("grey98", "lightcyan3", viridis::viridis(9, direction = 1)[c(4,3,2,1)]))(n_max)

plot(dem_sub_swiss, axes = F, legend = F,  col = col_elev, box = F)
plot(basin_bern, add =T, lwd = 0.025, col = alpha("steelblue4", alpha = 0.45))
plot(basin_brug, add =T, lwd = 0.025, col = alpha("steelblue4", alpha = 0.25))
plot(basin_unte, add =T, lwd = 0.025, col = alpha("steelblue4", alpha = 0.35))
plot(basin_diep, add =T, lwd = 0.025, col = alpha("darkolivegreen", alpha = 0.45))
plot(basin_neuh, add =T, lwd = 0.025, col = alpha("darkolivegreen", alpha = 0.25))
plot(basin_reki, add =T, lwd = 0.025, col = alpha("darkolivegreen", alpha = 0.35))
plot(basin_bern, add =T, lwd = 0.8, border = "black")
plot(basin_brug, add =T, lwd = 0.8, border = "black")
plot(basin_unte, add =T, lwd = 0.8, border = "black")
plot(basin_diep, add =T, lwd = 0.8, border = "black")
plot(basin_neuh, add =T, lwd = 0.8, border = "black")
plot(basin_reki, add =T, lwd = 0.8, border = "black")
plot(basin_base, add =T, lwd = 2.5, border = "gold3")
plot(rhin_riv, add = T, col = "blue4", lwd = 1.5)
plot(aare_riv, add = T, col = "blue4", lwd = 1.5)
plot(brienzersee, add = T, col = "blue4")
plot(thunersee, add = T, col = "blue4")
plot(bielersee, add = T, col = "blue4")
plot(bodensee, add = T, col = "blue4")
plot(gauges, add = T, pch = 23, cex = 1.7,
     bg = alpha(c("gold3", "darkolivegreen", "darkolivegreen", "steelblue4",
                   "steelblue4", "steelblue4", "darkolivegreen"), alpha = 1.0))
plot(gauges, add = T, pch = 19, cex = 0.5)
plot(snow, add = T, pch = "*", cex = 1.3)
lab_mov <- 5200
lab_pos_1 <- c(lab_mov, -lab_mov, rep(lab_mov, 2), +lab_mov, rep(lab_mov, 7), +lab_mov, lab_mov,        0, -lab_mov, -lab_mov-4000,  lab_mov+3000, lab_mov+4000, -lab_mov)
lab_pos_2 <- c(lab_mov, -lab_mov, rep(lab_mov, 2), -lab_mov, rep(lab_mov, 7), -lab_mov, lab_mov, -lab_mov, lab_mov,  lab_mov,   lab_mov,       0,  -lab_mov)
text(snow@coords[, 1]+lab_pos_1, snow@coords[, 2]+lab_pos_2, labels = stat_meta$ID, col = "black", cex = 0.9)
prettymapr::addscalebar(plotunit = "m", widthhint = 0.25, htin = 0.15, pos = "topleft",
                        padin = c(0.15, 0.15))

dev.off()


#map_ins----

dem = raster(paste0(data_dir, "basin_data/eu_dem/processed/eu_dem_1000.tif"))

# pdf(paste0(base_dir,"R/figs_exp/map_ins_raw.pdf"), width = 10, height = 6)
# tiff(paste0(base_dir,"R/figs_exp/map_ins_raw.tiff"), width = 10, height = 6,
#      units = "in", res = 800)
png(paste0(base_dir,"R/figs_exp/map_ins_raw.png"), width = 1000, height = 600)

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


plot(dem, col = alpha("grey35", alpha = 1.0), axes = F, legend = F, box = F)
plot(polys[2], col = alpha("grey35", alpha = 1.0), border = "red3", add = T, lwd = 4)

dev.off()


#elev_dist----

pdf(paste0(base_dir,"R/figs_exp/elev_dist.pdf"), width = 6, height = 2)

#Load DEM
dem = raster(paste0(data_dir, "basin_data/eu_dem/processed/eu_dem_500.tif"))

#Load basin boundaries
basins <-  rgdal::readOGR(dsn = paste0(data_dir, "basin_data/EZG_Schweiz_BAFU/ezg_kombiniert.shp"), encoding = "UTF8")
basin_base <- spTransform(basins[basins@data$Ort == "Basel, Rheinhalle",], CRS = crs(dem, asText = T))
# basin_unte <- spTransform(basins[basins@data$Ort == "Untersiggenthal, Stilli",], CRS = crs(dem, asText = T))
# basin_brug <- spTransform(basins[basins@data$Ort == "Brugg",], CRS = crs(dem, asText = T))
# basin_bern <- spTransform(basins[basins@data$Ort == "Bern, Schönau",], CRS = crs(dem, asText = T))
# basin_neuh <- spTransform(basins[basins@data$Ort == "Neuhausen, Flurlingerbrücke",], CRS = crs(dem, asText = T))
# basin_reki <- spTransform(basins[basins@data$Ort == "Rekingen",], CRS = crs(dem, asText = T))
# basin_diep <- spTransform(basins[basins@data$Ort == "Diepoldsau, Rietbrücke",], CRS = crs(dem, asText = T))

#corp DEM sub-basin area
dem_cro_base <- raster::crop(dem, extent(basin_base))
dem_sub_base <- mask(dem_cro_base, basin_base)
# dem_cro_unte <- raster::crop(dem, extent(basin_unte))
# dem_sub_unte <- mask(dem_cro_unte, basin_unte)
# dem_cro_brug <- raster::crop(dem, extent(basin_brug))
# dem_sub_brug <- mask(dem_cro_brug, basin_brug)
# dem_cro_bern <- raster::crop(dem, extent(basin_bern))
# dem_sub_bern <- mask(dem_cro_bern, basin_bern)
# dem_cro_reki <- raster::crop(dem, extent(basin_reki))
# dem_sub_reki <- mask(dem_cro_reki, basin_reki)
# dem_cro_neuh <- raster::crop(dem, extent(basin_neuh))
# dem_sub_neuh <- mask(dem_cro_neuh, basin_neuh)
# dem_cro_diep <- raster::crop(dem, extent(basin_diep))
# dem_sub_diep <- mask(dem_cro_diep, basin_diep)


#get elevations of cells cropped dem
dem_ele_NA_base <- dem_sub_base@data@values
dem_ele_base <- dem_ele_NA_base[!is.na(dem_ele_NA_base)]
# dem_ele_NA_unte <- dem_sub_unte@data@values
# dem_ele_unte <- dem_ele_NA_unte[!is.na(dem_ele_NA_unte)]
# dem_ele_NA_brug <- dem_sub_brug@data@values
# dem_ele_brug <- dem_ele_NA_brug[!is.na(dem_ele_NA_brug)]
# dem_ele_NA_bern <- dem_sub_bern@data@values
# dem_ele_bern <- dem_ele_NA_bern[!is.na(dem_ele_NA_bern)]
# dem_ele_NA_reki <- dem_sub_reki@data@values
# dem_ele_reki <- dem_ele_NA_reki[!is.na(dem_ele_NA_reki)]
# dem_ele_NA_neuh <- dem_sub_neuh@data@values
# dem_ele_neuh <- dem_ele_NA_neuh[!is.na(dem_ele_NA_neuh)]
# dem_ele_NA_diep <- dem_sub_diep@data@values
# dem_ele_diep <- dem_ele_NA_diep[!is.na(dem_ele_NA_diep)]

hist_breaks <- seq(4000, 200, -50)

par(mar = c(1, 2.2, 2.0, 1))

par(family = "serif")

hist_res <- hist(dem_ele_base, breaks = hist_breaks, plot = F)

ylabs <- (hist_res$counts / length(dem_ele_base))  * 100

areal_perc <- c(round(length(which(dem_ele_base < 1000)) / length(dem_ele_base) * 100, digits = 2),
                round(length(which(dem_ele_base > 1000 & dem_ele_base < 2000)) / length(dem_ele_base) * 100, digits = 2),
                round(length(which(dem_ele_base > 2000 & dem_ele_base < 3000)) / length(dem_ele_base) * 100, digits = 2),
                round(length(which(dem_ele_base > 3000)) / length(dem_ele_base) * 100, digits = 2))

hist(dem_ele_base, breaks = hist_breaks, axes = F, ylab = "", xlab = "", main = "",
     col = alpha("gold3", alpha = 0.8), yaxs = "i", ylim = c(0,16000), xlim = c(4000, 0))
# rect(xleft =  6000, xright = -300, ytop = 150000, ybottom = 0 ,col = alpha("grey", alpha = 0.30), border = NA)
box()
axis(3, mgp=c(3, 0.08, 0), tck = -0.01, cex.axis = 0.8)
axis(2, mgp=c(3, 0.08, 0), tck = -0.01, cex.axis = 0.8,
     at = c(length(dem_ele_base) * c(0.01, 0.025, 0.05, 0.075, 0.1)), labels = c("1.0", "2.5", "5.0", "7.5", "10.0"))
abline(v = c(0, 1000, 2000, 3000, 4000), col = "black", lty = "dashed", lwd = 1.0)
mtext("Elevation [m]", side = 3, line = 1.1, adj = 0.5, cex = 1.2)
mtext("Areal fraction [%]", side = 2, line = 1.2, adj = 0.5, cex = 1.2)
mtext(paste(areal_perc[1], "%"), side = 2, line = -21.0, adj = 0.65, col = "black")
mtext(paste(areal_perc[2], "%"), side = 2, line = -17.0, adj = 0.65, col = "black")
mtext(paste(areal_perc[3], "%"), side = 2, line = -11.0, adj = 0.65, col = "black")
mtext(paste(areal_perc[4], "%"), side = 2, line = -4.5, adj = 0.65, col = "black")

dev.off()

#meteo_plots----



layout(matrix(c(rep(1, 8), 2, rep(3, 8), 4,
                rep(5, 8), 6, rep(7, 8), 8),
              2, 18, byrow = T), widths=c(), heights=c())

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15

#Plot: Temperature mean
data_day <- tmea_band

cols_max <- colorRampPalette(c("grey98", "gold3", "darkgoldenrod3", "darkorange4", "darkred"))(100)
cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "grey98"))(100)
cols_hydro <- c(cols_min, cols_max)
breaks_hydro <- c(seq(min_na(data_day), 0, length.out = 100),
                  seq(0, max_na(data_day), length.out = 100+1))

par(mar = c(1.7, 3.7, 2.5, 0.2))

image(x = 1:nrow(data_day),
      y = seq(250, 3150, 50),
      z = data_day,
      col = cols_hydro,
      breaks = breaks_hydro,
      ylab = "", xlab = "", axes = F)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("O","N","D","J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.4)#plot labels
axis(2, mgp=c(3, 0.20, 0), tck = -0.03, cex.axis = 1.5)
mtext("a) Temperature mean", side= 3, line = 0.3, cex = 1.6, adj = 0.0)
mtext("[°C]", side= 3, line = 0.3, cex = 1.6, adj = 1.0)
mtext("Elevation", side= 2, line = 1.9, cex = 1.6)

box()

par(mar = c(1.7, 0.2, 2.5, 1.7))

alptempr::image_scale(as.matrix(data_day), col = cols_hydro, breaks = breaks_hydro, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.20, 0), tck = -0.10, cex.axis = 1.4)
box()


#Plot: Temperature trend
data_day <- tslo_band

cols_max <- colorRampPalette(c("grey98", "gold3", "darkgoldenrod3", "darkorange4", "darkred"))(100)
cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "grey98"))(100)
cols_hydro <- c(cols_min, cols_max)
breaks_hydro <- seq(-max_na(abs(data_day)), max_na(abs(data_day)), length.out = length(cols_hydro)+1)

par(mar = c(1.7, 3.7, 2.5, 0.2))

image(x = 1:nrow(data_day),
      y = seq(250, 3150, 50),
      z = data_day,
      col = cols_hydro,
      breaks = breaks_hydro,
      ylab = "", xlab = "", axes = F)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("O","N","D","J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.4)#plot labels
axis(2, mgp=c(3, 0.20, 0), tck = -0.03, cex.axis = 1.5)
mtext("b) Temperature trend", side= 3, line = 0.3, cex = 1.6, adj = 0.0)
mtext(expression(paste("[°C ", "dec"^"-1", "]")), side= 3, line = 0.3, cex = 1.6, adj = 1.0)
mtext("Elevation", side= 2, line = 1.9, cex = 1.6)

box()

par(mar = c(1.7, 0.2, 2.5, 1.7))

alptempr::image_scale(as.matrix(data_day), col = cols_hydro, breaks = breaks_hydro, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.20, 0), tck = -0.10, cex.axis = 1.4)
box()


#Plot: Precipitation mean
data_day <- pmea_band

cols_hydro <- grDevices::colorRampPalette(c("grey98", viridis::viridis(9, direction = 1)[c(4,3,2,1,1)]))(100)

breaks_hydro <- c(seq(min_na(data_day), max_na(data_day),  length.out = 100+1))

cols_hydro <- colorRampPalette(c(viridis::viridis(9, direction = 1)[9:1]))(200)
breaks_hydro <- seq(0, alptempr::max_na(data_day), length.out = length(cols_hydro)+1)

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15

par(mar = c(1.7, 3.7, 2.5, 0.2))

image(x = 1:nrow(data_day),
      y = seq(250, 3150, 50),
      z = data_day,
      col = cols_hydro,
      breaks = breaks_hydro,
      ylab = "", xlab = "", axes = F)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("O","N","D","J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.4)#plot labels
axis(2, mgp=c(3, 0.20, 0), tck = -0.03, cex.axis = 1.5)
mtext("c) Precipitation mean", side= 3, line = 0.3, cex = 1.6, adj = 0.0)
mtext("[mm]", side= 3, line = 0.3, cex = 1.6, adj = 1.0)
mtext("Elevation", side= 2, line = 1.9, cex = 1.6)

box()

par(mar = c(1.7, 0.2, 2.5, 1.7))

alptempr::image_scale(as.matrix(data_day), col = cols_hydro, breaks = breaks_hydro, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.20, 0), tck = -0.10, cex.axis = 1.4)
# mtext("cm", side = 3, line = 0.2, cex = 1.0)
box()


#Plot: Precipitation trend
data_day <- pslo_band

cols_max <- colorRampPalette(c("grey98", "gold3", "darkgoldenrod3", "darkorange4", "darkred"))(100)
cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "grey98"))(100)
cols_hydro <- c(cols_min, cols_max)
breaks_hydro <- seq(-max_na(abs(data_day)), max_na(abs(data_day)), length.out = length(cols_hydro)+1)

par(mar = c(1.7, 3.7, 2.5, 0.2))

image(x = 1:nrow(data_day),
      y = seq(250, 3150, 50),
      z = data_day,
      col = cols_hydro,
      breaks = breaks_hydro,
      ylab = "", xlab = "", axes = F)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.08)#plot ticks
axis(1, at = x_axis_lab, c("O","N","D","J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.4)#plot labels
axis(2, mgp=c(3, 0.20, 0), tck = -0.03, cex.axis = 1.5)
mtext("b) Precipitation trend", side= 3, line = 0.3, cex = 1.6, adj = 0.0)
mtext(expression(paste("[mm ", "dec"^"-1", "]")), side= 3, line = 0.3, cex = 1.6, adj = 1.0)
mtext("Elevation", side= 2, line = 1.9, cex = 1.6)

box()

par(mar = c(1.7, 0.2, 2.5, 1.7))

alptempr::image_scale(as.matrix(data_day), col = cols_hydro, breaks = breaks_hydro, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.20, 0), tck = -0.10, cex.axis = 1.4)
box()





#snow_stats----

sta_yea <- 1954
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

# pdf(paste0(base_dir,"R/figs_exp/snow_stats.pdf"), width = 12, height = 14.5)
png(paste0(base_dir,"R/figs_exp/snow_stats.png"), width = 12, height = 14.5,
    units = "in", res = 300)
# tiff(paste0(base_dir,"R/figs_exp/snow_stats.tiff"), width = 12, height = 14.5, 
#     units = "in", res = 300)

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
       col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.4)#plot labels
  axis(2, mgp=c(3, 0.20, 0), tck = -0.03, cex.axis = 1.5)
  if(do_head){
    mtext("a) Snow depth raster", side= 3, line = 0.5, cex = 1.6, adj = 0.0)
    mtext("[cm]", side= 3, line = 0.3, cex = 1.4, adj = 1.0)
    # legend("topleft", c("1959-1988", "1989-2018"), pch = 19, cex = 1.0, col = c(col_1, col_2), bg = "white")
  }
  box()
  
  par(mar = c(1.7, 0.2, 0.5, 1.7))
  
  alptempr::image_scale(as.matrix(data_day), col = cols_hydro, breaks = breaks_hydro, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
  axis(4, mgp=c(3, 0.20, 0), tck = -0.10, cex.axis = 1.4)
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
  axis(2, mgp=c(3, 0.2, 0), tck = -0.03, cex.axis = 1.5)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.08)#plot ticks
  axis(1, at = x_axis_lab, c("O","N","D","J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.4)#plot labels
  if(do_head){
    mtext("b) Snow depth mean", side= 3, line = 0.5, cex = 1.6, adj = 0.0)
    mtext("[cm]", side= 3, line = 0.3, cex = 1.4, adj = 1.0)
    legend("topright", c("1958-1988", "1988-2018"), pch = 19, cex = 1.25, col = c(col_1, col_2), bg = "white")
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
  axis(2, mgp=c(3, 0.2, 0), tck = -0.03, cex.axis = 1.5)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.08)#plot ticks
  axis(1, at = x_axis_lab, c("O","N","D","J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.4)#plot labels
  if(do_head){
    mtext("c) Acc./Melt rate", side= 3, line = 0.5, cex = 1.6, adj = 0.0)
    mtext(expression(paste("[cm ","day"^"-1", "]")), side= 3, line = 0.15, cex = 1.4, adj = 1.0)
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
12, 28), widths=c(1.2, rep(1, 27)), heights=c(0.2, rep(1, 10)))

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
cex_header <- 1.4
par(mar = c(0,0,0,0))

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("WFJ",  side = 2, line = -1.8, cex = cex_header, adj = 0.958, outer = T)
mtext("2691 m",  side = 2, line = -3.5, cex = cex_header, adj = 0.967, outer = T)
mtext("ARO",  side = 2, line = -1.8, cex = cex_header, adj = 0.865, outer = T)
mtext("1878 m",  side = 2, line = -3.5, cex = cex_header, adj = 0.871, outer = T)
mtext("SIA",  side = 2, line = -1.8, cex = cex_header, adj = 0.771, outer = T)
mtext("1804 m",  side = 2, line = -3.5, cex = cex_header, adj = 0.776, outer = T)
mtext("GRC",  side = 2, line = -1.8, cex = cex_header, adj = 0.679, outer = T)
mtext("1606 m",  side = 2, line = -3.5, cex = cex_header, adj = 0.682, outer = T)
mtext("DAV",  side = 2, line = -1.8, cex = cex_header, adj = 0.588, outer = T)
mtext("1594 m",  side = 2, line = -3.5, cex = cex_header, adj = 0.588, outer = T)
mtext("ANT",  side = 2, line = -1.8, cex = cex_header, adj = 0.494, outer = T)
mtext("1438 m",  side = 2, line = -3.5, cex = cex_header, adj = 0.493, outer = T)
mtext("SMM",  side = 2, line = -1.8, cex = cex_header, adj = 0.400, outer = T)
mtext("1386 m",  side = 2, line = -3.5, cex = cex_header, adj = 0.397, outer = T)
mtext("ABO",  side = 2, line = -1.8, cex = cex_header, adj = 0.305, outer = T)
mtext("1322 m",  side = 2, line = -3.5, cex = cex_header, adj = 0.303, outer = T)
mtext("DIS",   side = 2, line = -1.8, cex = cex_header, adj = 0.214, outer = T)
mtext("1197 m",   side = 2, line = -3.5, cex = cex_header, adj = 0.209, outer = T)
mtext("ELM",   side = 2, line = -1.8, cex = cex_header, adj = 0.120, outer = T)
mtext("957 m",   side = 2, line = -3.5, cex = cex_header, adj = 0.117, outer = T)
mtext("SMA",   side = 2, line = -1.8, cex = cex_header, adj = 0.031, outer = T)
mtext("555 m",   side = 2, line = -3.5, cex = cex_header, adj = 0.027, outer = T)

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

# pdf(paste0(base_dir,"R/figs_exp/disc_perc.pdf"), width = 12, height = 15)
# tiff(paste0(base_dir,"R/figs_exp/disc_perc.tiff"), width = 12, height = 15,
#      units = "in", res = 300)
png(paste0(base_dir,"R/figs_exp/disc_perc.png"), width = 12, height = 15,
    units = "in", res = 300)

layout(matrix(c(32, rep(29, 18),
                32, rep(1, 8), 2,   rep(3, 8), 4,
                32, rep(30, 18),
                32, rep(5, 8), 6,   rep(7, 8), 8,
                32, rep(9, 8), 10,  rep(11, 8), 12,
                32, rep(13, 8), 14, rep(15, 8), 16,
                32, rep(31, 18),
                32, rep(17, 8), 18, rep(19, 8), 20,
                32, rep(21, 8), 22, rep(23, 8), 24,
                32, rep(25, 8), 26, rep(27, 8), 28), 
              10, 19, byrow = T), widths=c(0.8, rep(1, 18)), heights=c(0.20, 1, 0.15, 1, 1, 1, 0.15, 1, 1, 1))

# layout.show(n = 32)

perce_plot <- function(data_in, date_in, main_in = "", year_1, year_2, year_3, year_4, 
                       do_ylab = T, do_unit = T, pos_unit = 1.0){
  
  data_day <- ord_day(data_in = data_in,
                      date = date_in,
                      start_y = year_1,
                      end_y = year_4)
  
  #Mean seasonal cycles
  ind_break_1 <- (year_1 - year_1) + 1
  ind_break_2 <- (year_2 - year_1) + 1
  ind_break_3 <- (year_3 - year_1) + 1
  ind_break_4 <- (year_4 - year_1) + 0
  
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
  
  month_cols <- list(oct_cols, nov_cols, dec_cols, jan_cols, feb_cols, mar_cols, apr_cols, may_cols, jun_cols, jul_cols, aug_cols, sep_cols)
  
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
  cols_max <- colorRampPalette(c("grey98", "lightcyan3", viridis::viridis(9, direction = 1)[c(4,3,2,1)]))(n_max)
  # my_col <- colorRampPalette(c(cols_min, cols_max))(200)
  # my_bre <- seq(min_na(qdif[, ]), max_na(qdif[, ]), length.out = 201)
  my_bre <- c(seq(-max_na(abs(qdif)), 0, length.out = n_min),
              seq(0, max_na(abs(qdif)), length.out = n_max+1))
  my_col <- colorRampPalette(c(cols_min, cols_max))(length(my_bre)-1)

  # plot(1:n_max, 1:n_max, pch = 19, cex = 2, col = cols_max)
  # plot(1:n_min, 1:n_min, pch = 19, cex = 2, col = cols_min)
  
  par(mar = c(1.6, 3.5, 1.5, 0.2))
  par(family = "serif")
  
  image(x = 1:12,
        y = seq(0.01, 0.99, 0.01),
        z = t(qdif), col = my_col, breaks = my_bre,
        ylab = "", xlab = "", axes = F)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.05)#plot ticks
  axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.5)#plot labels
  axis(2, mgp=c(3, 0.25, 0), tck = -0.025, cex.axis = 1.5)
  if(do_ylab){
    mtext("Prob. level", side = 2, line = 1.7, cex = 1.25) 
  }
  mtext(main_in, side = 3, line = 0.25, cex = 1.5, adj = 0.0)
  if(do_unit){
    mtext(expression(paste("[m"^"3", " s"^"-1", "]")), side = 3, line = 0.2, cex = 1.2, adj = pos_unit)
  }
  
  box()
  
  par(mar = c(1.6, 0.2, 1.5, 3.5))
  
  alptempr::image_scale(as.matrix(qdif), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
  axis(4, mgp=c(3, 0.45, 0), tck = -0.1, cex.axis = 1.5)
  box()
  
}

perce_plot(data_in = grdc_data_base$value,
           date_in = grdc_data_base$date,
           main_in = "",
           year_1 = 1951, year_2 = 1982, year_3 = 1983, year_4 = 2014)

perce_plot(data_in = grdc_data_base$value,
           date_in = grdc_data_base$date,
           main_in = "",
           year_1 = 1919, year_2 = 1967, year_3 = 1968, year_4 = 2016)


perce_plot(data_in = grdc_data_unte$value,
           date_in = grdc_data_unte$date,
           main_in = "", do_unit = F,
           year_1 = 1951, year_2 = 1982, year_3 = 1983, year_4 = 2014)

perce_plot(data_in = grdc_data_unte$value,
           date_in = grdc_data_unte$date,
           main_in = "", do_ylab = F, do_unit = F,
           year_1 = 1919, year_2 = 1967, year_3 = 1968, year_4 = 2016)

perce_plot(data_in = grdc_data_brug$value,
           date_in = grdc_data_brug$date,
           main_in = "", do_unit = F,
           year_1 = 1951, year_2 = 1982, year_3 = 1983, year_4 = 2014)

perce_plot(data_in = grdc_data_brug$value,
           date_in = grdc_data_brug$date,
           main_in = "", do_ylab = F, do_unit = F,
           year_1 = 1919, year_2 = 1967, year_3 = 1968, year_4 = 2016)

perce_plot(data_in = grdc_data_bern$value,
           date_in = grdc_data_bern$date,
           main_in = "", do_unit = F,
           year_1 = 1951, year_2 = 1982, year_3 = 1983, year_4 = 2014)

perce_plot(data_in = grdc_data_bern$value,
           date_in = grdc_data_bern$date,
           main_in = "", do_ylab = F,  do_unit = F,
           year_1 = 1919, year_2 = 1967, year_3 = 1968, year_4 = 2016)



perce_plot(data_in = grdc_data_reki$value,
           date_in = grdc_data_reki$date,
           main_in = "", do_unit = F,
           year_1 = 1951, year_2 = 1982, year_3 = 1983, year_4 = 2014)

perce_plot(data_in = grdc_data_reki$value,
           date_in = grdc_data_reki$date,
           main_in = "", do_ylab = F, do_unit = F,
           year_1 = 1919, year_2 = 1967, year_3 = 1968, year_4 = 2016)

perce_plot(data_in = grdc_data_neuh$value,
           date_in = grdc_data_neuh$date,
           main_in = "", do_unit = F,
           year_1 = 1951, year_2 = 1982, year_3 = 1983, year_4 = 2014)

perce_plot(data_in = grdc_data_neuh$value,
           date_in = grdc_data_neuh$date,
           main_in = "", do_ylab = F,  do_unit = F,
           year_1 = 1919, year_2 = 1967, year_3 = 1968, year_4 = 2016)

perce_plot(data_in = grdc_data_diep$value,
           date_in = grdc_data_diep$date,
           main_in = "", do_unit = F,
           year_1 = 1951, year_2 = 1982, year_3 = 1983, year_4 = 2014)

perce_plot(data_in = grdc_data_diep$value,
           date_in = grdc_data_diep$date,
           main_in = "", do_ylab = F, do_unit = F,
           year_1 = 1919, year_2 = 1967, year_3 = 1968, year_4 = 2016)

#Time frames
par(mar = c(0,0,0,0))
plot(1:10, 1:10, axes = F, type = "n", ylab = "", xlab = "")
mtext("a) 1951-2014", side = 3, line = -2.8, adj = 0.20, cex = 1.6)
mtext("b) 1919-2016", side = 3, line = -2.8, adj = 0.78, cex = 1.6)

#Aare branch
par(mar = c(0,0,0,0))
plot(1:10, 1:10, axes = F, type = "n", ylab = "", xlab = "")
mtext("Aare branch", side = 3, line = -3.0, adj = 0.50, cex = 1.6)

#Rhine branch
par(mar = c(0,0,0,0))
plot(1:10, 1:10, axes = F, type = "n", ylab = "", xlab = "")
mtext("Rhine branch", side = 3, line = -3.0, adj = 0.50, cex = 1.6)

#Gauges
cex_header <- 1.6
par(mar = c(0,0,0,0))

plot(1:100, 1:100, axes = F, type = "n", xlab = "", ylab = "")
mtext("1. Basel",            side = 2, line = -2.8, cex = cex_header, adj = 0.929, outer = T)
mtext("2. Untersiggenthal",  side = 2, line = -2.8, cex = cex_header, adj = 0.789, outer = T)
mtext("3. Brugg",            side = 2, line = -2.8, cex = cex_header, adj = 0.627, outer = T)
mtext("4. Bern",             side = 2, line = -2.8, cex = cex_header, adj = 0.490, outer = T)
mtext("5. Rekingen",         side = 2, line = -2.8, cex = cex_header, adj = 0.315, outer = T)
mtext("6. Neuhausen",        side = 2, line = -2.8, cex = cex_header, adj = 0.170, outer = T)
mtext("7. Diepoldsau",       side = 2, line = -2.8, cex = cex_header, adj = 0.015, outer = T)

dev.off()


#disc_rast----

# pdf(paste0(base_dir,"R/figs_exp/disc_rast.pdf"), width = 12, height = 8)
# tiff(paste0(base_dir,"R/figs_exp/disc_rast.tiff"), width = 12, height = 8,
#      units = "in", res = 300)
png(paste0(base_dir,"R/figs_exp/disc_rast.png"), width = 12, height = 8,
    units = "in", res = 300)

layout(matrix(c(rep(c(15, 3, 7, 11), 4),
                rep(c(1, 3, 7, 11), 3),
                rep(c(1, 4, 8, 12), 1),
                rep(c(1, 5, 9, 13), 3),
                rep(c(2, 5, 9, 13), 1),
                rep(c(16, 5, 9, 13), 3),
                rep(c(16, 6, 10, 14), 1)
),
4, 16), widths=c(), heights=c())

sta_yea_cla <- 1919
yea_cla_1 <- sta_yea_cla
yea_cla_2 <- 1967
yea_cla_3 <- 1968
yea_cla_4 <- 2016


raster_plot <- function(data_in, date_in, main_in = "", year_1  = yea_cla_1, year_2 = yea_cla_4,
                        quant_break = 0.50){
  
  data_day <- ord_day(data_in = data_in,
                      date = date_in,
                      start_y = year_1,
                      end_y = year_2,
                      break_day = 274)
  
  #Preparation Raster graph
  x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
  x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
  
  # cols_min <- grDevices::colorRampPalette(c(viridis::viridis(9, direction = 1)[c(1,1, 2:4)], "grey98"))(100)
  # cols_max <- grDevices::colorRampPalette(c("grey98", "gold3", "orangered4", "darkred"))(100)
  cols_min <- colorRampPalette(c("darkred", "darkorange4", "goldenrod3", "gold3", "lightgoldenrod2", "grey98"))(50)
  cols_max <- colorRampPalette(c("grey98", "lightcyan3", viridis::viridis(9, direction = 1)[c(4,3,2,1)]))(50)
  cols_hydro <- colorRampPalette(c(cols_min, cols_max))(200)
  # cols_hydro <- grDevices::colorRampPalette(c(viridis::viridis(9, direction = -1)))(200)
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
        y = year_1:year_2,
        z = t(data_day),
        col = cols_hydro,
        breaks = breaks_hydro,
        ylab = "", xlab = "", axes = F)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.05)#plot ticks
  axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.30, 0), cex.axis = 1.5)#plot labels
  axis(2, mgp=c(3, 0.15, 0), tck = -0.025, cex.axis = 1.5)
  mtext("Year", side = 2, line = 1.7, cex = 1.25)
  mtext(main_in, side = 3, line = 0.25, cex = 1.5, adj = 0.0)
  mtext(expression(paste("[m"^"3", " s"^"-1", "]")), side = 3, line = 0.2, cex = 1.2, adj = 1.0)
  
  
  box()
  
  par(mar = c(1.6, 0.2, 2.0, 3.5))
  
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

pdf(paste0(base_dir,"R/figs_exp/comp_illu_raw.pdf"), width = 8, height = 3)

par(mar = c(1.3, 3.5, 2.0, 0.2))
par(family = "serif")

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
col_1 <- alpha("steelblue4", alpha = 0.8)
col_2 <- alpha("grey55", alpha = 0.8)

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
     col = "black", col.axis = "black", tck = -0.05)#plot ticks
axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.27, 0), cex.axis = 1.5)#plot labels
axis(2, mgp=c(3, 0.25, 0), tck = -0.01, cex.axis = 1.5)
mtext("Elevation [m]", side = 2, line = 1.8, cex = 1.7, adj = 0.5)
mtext("b) Snowmelt elevation compensation", side = 3, line = 0.25, cex = 1.7, adj = 0.0)
legend("topleft", c("future", "present"), col = c(col_2, col_1), pch = 19,  box.lwd = 1, box.col = "black", bg = "white", cex = 1.4)
Arrows(x0 = 166, x1 =156, y0 = 1250, y1 = 1250, arr.type = "triangle", arr.width = 0.15, arr.length = 0.2, lwd = 1.2)
Arrows(x0 = 205, x1 =192, y0 = 2000, y1 = 2000, arr.type = "triangle", arr.width = 0.15, arr.length = 0.2, lwd = 1.2)
Arrows(x0 = 254, x1 =243, y0 = 3000, y1 = 3000, arr.type = "triangle", arr.width = 0.15, arr.length = 0.2, lwd = 1.2)
Arrows(x0 = 130, x1 =130, y0 = 750,  y1 = 1100, arr.type = "triangle", col = "black", arr.width = 0.15, arr.length = 0.2, lwd = 1.2)
Arrows(x0 = 162, x1 =162, y0 = 1650, y1 = 2000, arr.type = "triangle", col = "black", arr.width = 0.15, arr.length = 0.2, lwd = 1.2)
Arrows(x0 = 208, x1 =208, y0 = 2650, y1 = 3000, arr.type = "triangle", col = "black", arr.width = 0.15, arr.length = 0.2, lwd = 1.2)
box()

dev.off()


#snow_simu----

# pdf(paste0(base_dir,"R/figs_exp/snow_simu.pdf"), width = 12, height = 10)
# tiff(paste0(base_dir,"R/figs_exp/snow_simu.tiff"), width = 12, height = 10,
#      units = "in", res = 300)
png(paste0(base_dir,"R/figs_exp/snow_simu.png"), width = 12, height = 10,
    units = "in", res = 300)


layout(matrix(c(rep(c(1, 5, 9, 13, 17), 8), 2, 6, 10, 14, 18,
                rep(c(3, 7, 11, 15, 19), 8), 4, 8, 12, 16, 20),
              5, 18), widths=c(), heights=c())
# layout.show(n=12)

snow_sim_plot <- function(data_plot, cols, breaks, header, lab_unit, elev_bands = my_elev_bands){
  
  par(mar = c(2.2, 3.5, 2.5, 0.2))
  par(family = "serif")
  
  x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
  x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
  
  image(x = 1:365,
        y = elev_bands[-length(elev_bands)],
        z = data_plot, col = cols, breaks = breaks,
        ylab = "", xlab = "", axes = F)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.06)#plot ticks
  axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.40, 0), cex.axis = 1.7)#plot labels
  axis(2, mgp=c(3, 0.25, 0), tck = -0.02, cex.axis = 1.6)
  mtext("Elevation", side = 2, line = 1.8, cex = 1.3)
  mtext(header, side = 3, line = 0.3, cex = 1.5, adj = 0.0)
  mtext(lab_unit, side = 3, line = 0.3, cex = 1.3, adj = 1.0)
  box()
  
  par(mar = c(2.2, 0.2, 2.5, 3.0))
  
  alptempr::image_scale(as.matrix(data_plot), col = cols, breaks = breaks, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
  axis(4, mgp=c(3, 0.35, 0), tck = -0.08, cex.axis = 1.5)
  box()
  
  
}

#SWE depth mean
my_col <- colorRampPalette(c("grey98", viridis::viridis(9, direction = 1)[4:1]))(200)
my_bre <- seq(alptempr::min_na(smea_band), alptempr::max_na(smea_band), length.out = length(my_col)+1)

snow_sim_plot(smea_band, cols = my_col, breaks = my_bre,
              header = "a) SWE depth mean", lab_unit = "[m]")

#SWE depth trend
cols_min <- colorRampPalette(c("darkred", "darkorange4", "darkgoldenrod3", "gold3","grey98"))(100)
cols_max <- colorRampPalette(c("grey98", viridis::viridis(9, direction = 1)[4:1]))(100)
# cols_max <- colorRampPalette(c("grey98", "gold3", "darkgoldenrod3", "darkorange4", "darkred"))(100)
# cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "grey98"))(100)

my_col <- alpha(c(cols_min, cols_max), alpha = 1.0)
my_bre <- seq(-max_na(abs(sslo_band)), max_na(abs(sslo_band)), length.out = length(my_col)+1)

snow_sim_plot(sslo_band, cols = my_col, breaks = my_bre,
              header = "b) SWE depth trend", lab_unit = expression(paste("[", "m", " decade"^"-1", "]"))) 

#SWE volume mean
my_col <- colorRampPalette(c("grey98", viridis::viridis(9, direction = 1)[4:1]))(200)
my_bre <- seq(alptempr::min_na(vmea_band), alptempr::max_na(vmea_band), length.out = length(my_col)+1)

snow_sim_plot(vmea_band, cols = my_col, breaks = my_bre,
              header = "c) SWE volume mean", lab_unit = expression(paste("[", "10"^"6", "m"^"3", "]"))) #[hm?]

#SWE volume trend
cols_min <- colorRampPalette(c("darkred", "darkorange4", "darkgoldenrod3", "gold3","grey98"))(100)
cols_max <- colorRampPalette(c("grey98", viridis::viridis(9, direction = 1)[4:1]))(100)
# cols_max <- colorRampPalette(c("grey98", "gold3", "darkgoldenrod3", "darkorange4", "darkred"))(100)
# cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "grey98"))(100)

my_col <- c(cols_min, cols_max)
my_bre <- seq(-max_na(abs(vslo_band)), max_na(abs(vslo_band)), length.out = length(my_col)+1)

snow_sim_plot(vslo_band, cols = my_col, breaks = my_bre,
              header = "d) SWE volume trend", lab_unit = expression(paste("[", "10"^"6", "m"^"3", " decade"^"-1", "]")))

#SWE volume diff mean
# cols_max <- colorRampPalette(c("grey98", "gold3", "darkgoldenrod3", "darkorange4", "darkred"))(100)
# cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "grey98"))(100)
cols_min <- colorRampPalette(c("darkred", "darkorange4", "darkgoldenrod3", "gold3","grey98"))(100)
cols_max <- colorRampPalette(c("grey98", viridis::viridis(9, direction = 1)[4:1]))(100)


my_col <- c(cols_min, cols_max)
my_bre <- seq(-max_na(abs(vdif_band)), max_na(abs(vdif_band)), length.out = length(my_col)+1)

snow_sim_plot(vdif_band, cols = my_col, breaks = my_bre,
              header = "e) Accum./Melt mean", lab_unit = expression(paste("[", "10"^"6", "m"^"3", " day"^"-1", "]")))

#SWE volume diff trend
# cols_max <- colorRampPalette(c("grey98", "gold3", "darkgoldenrod3", "darkorange4", "darkred"))(100)
# cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "grey98"))(100)
cols_min <- colorRampPalette(c("darkred", "darkorange4", "darkgoldenrod3", "gold3","grey98"))(100)
cols_max <- colorRampPalette(c("grey98", viridis::viridis(9, direction = 1)[4:1]))(100)


my_col <- c(cols_min, cols_max)
my_bre <- seq(-max_na(abs(vdis_band)), max_na(abs(vdis_band)), length.out = length(my_col)+1)

snow_sim_plot(vdis_band, cols = my_col, breaks = my_bre,
              header = "f) Accum./Melt trend", lab_unit = expression(paste("[", "10"^"6", "m"^"3", " day"^"-1", "decade"^"-1", "]"))) #

#Temperature mean
cols_max <- colorRampPalette(c("grey98", "gold3", "darkgoldenrod3", "darkorange4", "darkred"))(100)
cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "grey98"))(100)
my_col <- c(cols_min, cols_max)
my_bre <- c(seq(min_na(tmea_band), 0, length.out = 100), seq(0, max_na(tmea_band), length.out = 100+1))

snow_sim_plot(tmea_band, cols = my_col, breaks = my_bre,
              header = "g) Temperature mean", lab_unit = expression(paste("[", "°C", "]"))) 

#Temperature trend
cols_max <- colorRampPalette(c("grey98", "gold3", "darkgoldenrod3", "darkorange4", "darkred"))(100)
cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "grey98"))(100)
my_col <- c(cols_min, cols_max)
my_bre <- seq(-max_na(abs(tslo_band)), max_na(abs(tslo_band)), length.out = length(my_col)+1)

snow_sim_plot(tslo_band, cols = my_col, breaks = my_bre,
              header = "h) Temperature trend", lab_unit = expression(paste("[", "°C ", "decade"^"-1", "]"))) 

#Precipitation mean
my_col <- colorRampPalette(c("grey98", viridis::viridis(9, direction = 1)[4:1]))(200)
my_col <- viridis::viridis(200, direction = -1)
# alptempr::min_na(pmea_band)
my_bre <- seq(0, alptempr::max_na(pmea_band), length.out = length(my_col)+1)


snow_sim_plot(pmea_band, cols = my_col, breaks = my_bre,
              header = "i) Precipitation mean", lab_unit = expression(paste("[", "mm", "]"))) 

#Precipitation trend
# cols_max <- colorRampPalette(c("grey98", "gold3", "darkgoldenrod3", "darkorange4", "darkred"))(100)
# cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "grey98"))(100)
cols_min <- colorRampPalette(c("darkred", "darkorange4", "darkgoldenrod3", "gold3","grey98"))(100)
cols_max <- colorRampPalette(c("grey98", viridis::viridis(9, direction = 1)[4:1]))(100)

my_col <- c(cols_min, cols_max)
my_bre <- seq(-max_na(abs(pslo_band)), max_na(abs(pslo_band)), length.out = length(my_col)+1)

snow_sim_plot(pslo_band, cols = my_col, breaks = my_bre,
              header = "j) Precipitation trend", lab_unit = expression(paste("[", "mm ", "decade"^"-1", "]"))) 


dev.off()







#varia_asso----


plot(doy_mea_2-doy_mea_1, type = "l")

cinp <- NULL
for(i in 1:length(doy_mea)){

  band_doy_mea <- doy_mea[i]
  cinp <- c(cinp, mea_na(cdis_band[((band_doy_mea-14):(band_doy_mea)), i]))
  
}

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
lab_months <- c("O", "N", "D","J","F","M","A","M","J","J","A","S")

par(mar = c(1.5, 1.5, 0.2, 0.2))
par(family ="serif")

plot(doy_mea, cinp, type = "l", axes = F, xlab = "", ylab = "", 
     ylim = c(range(cinp, na.rm = T)), xlim = c(0, 365))
abline(v = x_axis_tic, lty = "dashed", lwd = 0.8, col = "grey55")
abline(h = 0, col = "grey65", lty = "dashed", lwd = 0.8)
axis(2, mgp = c(3, 0.3, 0), tck = -0.02, cex.axis = 1.7)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
for(i in 1:length(x_axis_lab)){
  axis(1, at = x_axis_lab[i], lab_months[i], tick = FALSE, col="black", col.axis="black", 
       mgp=c(4, 0.45, 0), cex.axis = 1.7)
}
# mtext(expression(paste("Temp. [?C", "decades"^"-1", "]")), side = 2, line = 1.6, cex = 1.5)
box()


#Plot: Temperature trends

pdf(paste0(base_dir,"R/figs_exp/varia_asso_1_14.pdf"), width = 5, height = 2)

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
lab_months <- c("O", "N", "D","J","F","M","A","M","J","J","A","S")

par(mar = c(1.5, 1.5, 0.2, 0.2))
par(family ="serif")

data_band_in <- tslo_band

plot(data_band_in[, 1], type = "n", axes = F, xlab = "", ylab = "", ylim = c(range(data_band_in, na.rm = T)))
abline(v = x_axis_tic, lty = "dashed", lwd = 0.8, col = "grey55")
abline(h = 0, col = "grey65", lty = "dashed", lwd = 0.8)
for(i in 1: ncol(data_band_in)){
  
  lines(data_band_in[, i], col = viridis(ncol(data_band_in))[i])
  
}
axis(2, mgp = c(3, 0.3, 0), tck = -0.02, cex.axis = 1.7)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
for(i in 1:length(x_axis_lab)){
  axis(1, at = x_axis_lab[i], lab_months[i], tick = FALSE, col="black", col.axis="black", 
       mgp=c(4, 0.45, 0), cex.axis = 1.7)
}
# mtext(expression(paste("Temp. [?C", "decades"^"-1", "]")), side = 2, line = 1.6, cex = 1.5)
box()

dev.off()



#Plot: Timing change

pdf(paste0(base_dir,"R/figs_exp/varia_asso_2.pdf"), width = 5, height = 2.0)

# my_elev_bands[12]
band_ind_min <- 1
cols_elev <- viridis(length(my_elev_bands[band_ind_min:length(my_elev_bands)]), direction = -1)

par(mar = c(1.5, 1.5, 0.2, 0.2))
par(family ="serif")

plot(1:10, 1:10, type = "n", xlim = c(0, 365), axes = F, ylab = "", xlab = "", ylim = range(yea_min_mag_1[band_ind_min:length(doy_mea_2)]-yea_min_mag_2[band_ind_min:length(doy_mea_2)]),
     yaxs = "i")
abline(v = x_axis_tic, lty = "dashed", lwd = 0.8, col = "grey55")
abline(h = 0, col = "grey65", lty = "dashed", lwd = 0.8)
lines(doy_mea[band_ind_min:length(doy_mea)],
       (yea_min_mag_1[band_ind_min:length(doy_mea_2)]-yea_min_mag_2[band_ind_min:length(doy_mea_2)]),
       col = "black", lwd = 1.2)
points(doy_mea[band_ind_min:length(doy_mea)],
      (doy_mea_2[band_ind_min:length(doy_mea_2)]-doy_mea_1[band_ind_min:length(doy_mea_2)]),
      col = cols_elev, pch = 19, cex = 0.9)
axis(2, at = c(0, -5, -10), labels = c("0", "-5", "-10"), mgp = c(3, 0.3, 0), 
     tck = -0.02, cex.axis = 1.7)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
for(i in 1:length(x_axis_lab)){
  axis(1, at = x_axis_lab[i], lab_months[i], tick = FALSE, col="black", col.axis="black", 
       mgp=c(4, 0.45, 0), cex.axis = 1.7)
}
box()

dev.off()


#Plot: Color scale elevations
pdf(paste0(base_dir,"R/figs_exp/varia_asso_3.pdf"), width = 0.5, height = 1.5)

par(mar = c(0.2, 0.2, 1.2, 1.1))

alptempr::image_scale(as.matrix(freq_diff_all), col = cols_elev, 
                      breaks = seq(from = my_elev_bands[band_ind_min], to = my_elev_bands[length(my_elev_bands)]+50, by = 50), 
                      horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.20, 0), tck = -0.00, cex.axis = 1.0,
     at = c(1500, 3000), labels = c("1500","3000"))
axis(4, mgp=c(3, 0.35, 0), tck = -0.08, cex.axis = 1.0,
     at = c(1000, 1500, 2000, 2500, 3000), labels = c("", "","", "", ""))

mtext("[m]", side = 3, line = 0.3, cex = 1)
box()

dev.off()


#Plot: Changes in melt integrated over elevation range

pdf(paste0(base_dir,"R/figs_exp/varia_asso_4.pdf"), width = 5, height = 2.0)

par(mar = c(1.5, 1.5, 0.2, 0.2))
par(family ="serif")

# my_elev_bands[12]
band_ind_min <- 1

vdims_inte <- apply(vdims_band[, band_ind_min:ncol(vdims_band)], 1, sum_na) / 1000000 #[hm?]
#re-arange so it starts 1st Jan
vdims_inte <- c(vdims_inte[93:length(vdims_inte)], vdims_inte[1:92])

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15

ylims = range(c(vdims_inte), na.rm = T)

plot(vdims_inte, type = "n", ylim = c(ylims[2], ylims[1]), axes = F, lwd = 2, ylab = "", xlab = "")
abline(v = x_axis_tic, lty = "dashed", lwd = 0.8, col = "grey55")
abline(h = 0, col = "grey65", lty = "dashed", lwd = 0.8)
lines(vdims_inte, lwd = 2)
axis(2, at = c(-0.5, 0, 0.5, 1), labels = c(-0.5, 0, 0.5, 1), mgp = c(3, 0.3, 0), 
     tck = -0.02, cex.axis = 1.7)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
for(i in 1:length(x_axis_lab)){
  axis(1, at = x_axis_lab[i], lab_months[i], tick = FALSE, col="black", col.axis="black", 
       mgp=c(4, 0.45, 0), cex.axis = 1.7)
}
box()

dev.off()



#Bands calculate
yea_min_mag_bands_14 <- NULL
yea_min_doy_bands_14 <- NULL

for(i in band_ind_min:ncol(svolu_d_band)){
  
  swe_band <- svolu_d_band[, i]
  
  swe_band_dif <- c(NA, diff(swe_band))
  
  swe_band_dif[which(swe_band_dif > 0)] <- 0
  
  #Moving average filter
  swe_band_dif_ma_14 <- rollapply(data = swe_band_dif, width = 14,
                                  FUN = sum_na, align = "center", fill = NA)
  
  #Order data by day
  data_day_14 <- ord_day(data_in = swe_band_dif_ma_14,
                         date = date_snow,
                         start_y = 1954,
                         end_y = 2014,
                         break_day = 274,
                         do_ma = F,
                         window_width = 30)
  
  min_doy <- function(data_in){
    
    doy_min <- as.numeric(which(data_in == min_na(data_in)))
    
    if(length(doy_min) > 1){
      doy_min <- mea_na(doy_min)
    }
    
    return(doy_min)
  }
  
  yea_min_mag_14 <- apply(data_day_14, 1, min_na)
  yea_min_doy_14 <- apply(data_day_14, 1, min_doy)
  
  yea_min_mag_bands_14 <- cbind(yea_min_mag_bands_14, yea_min_mag_14)
  yea_min_doy_bands_14 <- cbind(yea_min_doy_bands_14, yea_min_doy_14)
  
}

bands_doy_mea <- apply(yea_min_doy_bands_14, 2, mea_na)
bands_doy_slo <- apply(yea_min_doy_bands_14, 2, sens_slo)*-1
bands_doy_mea_1 <- apply(yea_min_doy_bands_14[1:30,], 2, med_na)
bands_doy_mea_2 <- apply(yea_min_doy_bands_14[31:60,], 2, med_na)

freq_bands_1 <- NULL
freq_bands_2 <- NULL
freq_diff_all <- NULL
length_bands_1 <- NULL
length_bands_2 <- NULL
width_over <- 14

for(b in 1:length(bands_doy_mea)){
  
  band_sel <- b
  
  bands_frequ_1 <- NULL
  band_lengt_1 <- NULL
  for(y in 1:30){
    
    mid_day <- yea_min_doy_bands_14[y, band_sel]
    min_day <- mid_day - width_over
    max_day <- mid_day + width_over
    melt_window <- min_day:max_day
    
    bands_asso_1 <- which(yea_min_doy_bands_14[y, ] > min_day & yea_min_doy_bands_14[y, ] < max_day)
    
    bands_frequ_1 <- c(bands_frequ_1, bands_asso_1)
    band_lengt_1 <- c(band_lengt_1, length(bands_asso_1))
    
  }
  hist_1 <- hist(bands_frequ_1, breaks = seq(from = 0.5, to = length(bands_doy_mea)+0.5, by = 1), plot = F)
  
  bands_frequ_2 <- NULL
  band_lengt_2 <- NULL
  for(x in 31:60){
    
    # print(y)
    mid_day <- yea_min_doy_bands_14[x, band_sel]
    min_day <- mid_day - width_over
    max_day <- mid_day + width_over
    melt_window <- min_day:max_day
    
    bands_asso_2 <- which(yea_min_doy_bands_14[x, ] > min_day & yea_min_doy_bands_14[x, ] < max_day)
    
    bands_frequ_2 <- c(bands_frequ_2, bands_asso_2)
    band_lengt_2 <- c(band_lengt_2, length(bands_asso_2))
    
  }
  hist_2 <- hist(bands_frequ_2, breaks = seq(from = 0.5, to = length(bands_doy_mea)+0.5, by = 1), plot = F)
  
  freq_diff <- hist_2$counts - hist_1$counts
  
  freq_bands_1 <- cbind(freq_bands_1, hist_1$counts)
  freq_bands_2 <- cbind(freq_bands_2, hist_2$counts)
  freq_diff_all <- cbind(freq_diff_all, freq_diff)
  length_bands_1 <- cbind(length_bands_1, band_lengt_1)
  length_bands_2 <- cbind(length_bands_2, band_lengt_2)
}

length_bands_1_mea <- apply(length_bands_1, 2, mea_na)
length_bands_2_mea <- apply(length_bands_2, 2, mea_na)

elevs_asso <- my_elev_bands[band_ind_min:length(my_elev_bands)]

elevs_sel <- c(6, 16, 26, 36, 46, 56)-1


#Plot Frequencies individual band

pdf(paste0(base_dir,"R/figs_exp/varia_asso_5.pdf"), width = 5, height = 2)

band_test <- 27 #elevs_asso[27]
col_1 <- viridis(9, direction = 1)[4]
col_2 <- "darkred"

par(mar = c(1.5, 1.5, 0.2, 0.2))
par(family ="serif")

plot(freq_bands_1[, band_test], type = "h", col = alpha(col_1, alpha = 0.6), lwd = 7, ylab = "", xlab = "",
     axes = F, lend = 2, xaxs = "i", yaxs = "i", ylim = c(0, 32), xlim = c(0, length(elevs_asso)))
legend("topleft", c("1985-2014", "1954-1984"), col = c(col_2, col_1), pch = 19, cex = 1.2)
par(new = T)
plot(freq_bands_2[, band_test], type = "h", col = alpha(col_2, alpha = 0.6), lwd = 7, ylab = "", xlab = "",
     axes = F, lend = 2, xaxs = "i", yaxs = "i", ylim = c(0, 32), xlim = c(0, length(elevs_asso)))
par(new = T)
plot(band_test, freq_bands_2[band_test, band_test], type = "h", col = alpha("black", alpha = 1.0), lwd = 7, ylab = "", xlab = "",
     axes = F, lend = 2, xaxs = "i", yaxs = "i", ylim = c(0, 32), xlim = c(0, length(elevs_asso)))
axis(1, at = elevs_sel, labels = elevs_asso[elevs_sel], tck = -0.02, cex.axis = 1.7,mgp = c(3, 0.4, 0))
axis(2, mgp = c(3, 0.3, 0), tck = -0.02, cex.axis = 1.7)
# mtext("Elevation band [m]", side = 1, adj = 0.5, line = 1.3, cex = 1.1)
# mtext("Frequency concurrent melt [-]", side = 2, adj = 0.5, line = 1.3, cex = 1.1)
# mtext("a) Frequency concurrent melt (elevation band 2050-2100 m)", side = 3, adj = 0.0, line = 0.2, cex = 1.4)
box()

dev.off()



pdf(paste0(base_dir,"R/figs_exp/varia_asso_6.pdf"), width = 5, height = 2)

par(mar = c(1.5, 1.5, 0.2, 0.2))
par(family ="serif")

plot(freq_diff_all[band_test, ], type = "n", axes = F, ylab = "", xlab = "", xaxs = "i", yaxs = "i",
     xlim = c(0,length(elevs_asso)), ylim = c(-12.9, 9))
abline(v = band_test, col = "black", lwd = 1.5)
lines(freq_diff_all[band_test, ], lwd = 2, col = "grey35")
polygon(x = c(0.999, 1:47, 58.0001) , y = c(0, freq_diff_all[band_test, ], 0), col = alpha("grey55", alpha = 0.7), border = F)
axis(1, at = elevs_sel, labels = elevs_asso[elevs_sel], tck = -0.02, cex.axis = 1.7,mgp = c(3, 0.4, 0))
axis(2, mgp = c(3, 0.3, 0), tck = -0.02, cex.axis = 1.7)
# abline(h = 0, lty = "dashed")
# mtext("Elevation band [m]", side = 1, adj = 0.5, line = 1.3, cex = 1.1)
# mtext("Change frequency concurrent melt [-]", side = 2, adj = 0.5, line = 1.3, cex = 1.1)
# mtext("b) Changes in concurrent melt (elevation band 2050-2100 m)", side = 3, adj = 0.0, line = 0.2, cex = 1.4)
box()

dev.off()



pdf(paste0(base_dir,"R/figs_exp/varia_asso_7.pdf"), width = 4.5, height = 3.5)

cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "cadetblue3", "grey90"))(50)
cols_max <- colorRampPalette(c("grey90", "gold3",  "orange3", "orangered4", "orangered4", "darkred"))(50)
cols_mel <- c(cols_min, cols_max)

layout(matrix(c(rep(1, 9), 2),
              1, 10), widths=c(), heights=c())

par(mar = c(2.3, 2.5, 0.7, 0.2))
par(family ="serif")

image(x = 1:nrow(freq_diff_all),
      y = 1:ncol(freq_diff_all),
      z = freq_diff_all,
      col = cols_mel, breaks = seq(from = -12, to = 12, length.out = 101), axes = F, ylab = "", xlab = "")
axis(1, at = elevs_sel, labels = elevs_asso[elevs_sel], mgp=c(3, 0.85, 0), tck = -0.01, cex.axis = 2.2)
axis(2, at = elevs_sel, labels = elevs_asso[elevs_sel], mgp=c(3, 0.55, 0), tck = -0.01, cex.axis = 2.2)
# mtext("Elevation band [m]", side = 1, adj = 0.5, line = 2.1, cex = 1.1)
# mtext("Elevation band [m]", side = 2, adj = 0.5, line = 2.1, cex = 1.1)
box()

par(mar = c(2.3, 0.2, 0.7, 2.2))

alptempr::image_scale(as.matrix(freq_diff_all), col = cols_mel, breaks = seq(from = -12, to = 12, length.out = 101), horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.85, 0), tck = -0.08, cex.axis = 2.2)
# mtext(lab_unit, side = 3, line = 0.3, cex = 1)
box()

dev.off()








#basin_flood----

#calculations

sno_vol_dif <- c(NA, diff(sno_vol_basin))

snow_vol_dif_mel <- sno_vol_dif
snow_vol_dif_mel[which(snow_vol_dif_mel > 0)] <- NA

sno_vol_dif_acc <- sno_vol_dif
sno_vol_dif_acc[which(sno_vol_dif < 0)] <- NA

#Moving average sum

sno_vol_dif_mel_ma_14 <- rollapply(data = snow_vol_dif_mel, width = 14,
                                   FUN = sum_na, align = "center", fill = NA)

sno_vol_dif_acc_ma_14 <- rollapply(data = sno_vol_dif_acc, width = 14,
                                   FUN = sum_na, align = "center", fill = NA)

#Order data by day

data_day_mel_14 <- ord_day(data_in = sno_vol_dif_mel_ma_14,
                           date = date_snow,
                           start_y = 1954,
                           end_y = 2014,
                           break_day = 274,
                           do_ma = F,
                           window_width = 30)

data_day_acc_14 <- ord_day(data_in = sno_vol_dif_acc_ma_14,
                           date = date_snow,
                           start_y = 1954,
                           end_y = 2014,
                           break_day = 274,
                           do_ma = F,
                           window_width = 30)

#Discharge data gauge Basel
grdc_base <- read_grdc(paste0(grdc_dir, "6935051_Q_Day.Cmd.txt"))

grdc_day_base <- ord_day(data_in = grdc_base$value,
                         date = grdc_base$date,
                         start_y = 1954,
                         end_y = 2014,
                         break_day = 274)

#Basin precipitation

prec_basin_ms_3 <- rollapply(data = prec_basin, width = 3,
                             FUN = sum_na, align = "center", fill = NA)

prec_basin_ms_3_sol <- rollapply(data = prec_basin_sol, width = 3,
                                 FUN = sum_na, align = "center", fill = NA)

prec_day_ms_3 <- ord_day(data_in = prec_basin_ms_3,
                         date = date_snow,
                         start_y = 1954,
                         end_y = 2014,
                         break_day = 274)

prec_day_ms_3_sol <- ord_day(data_in = prec_basin_ms_3_sol,
                             date = date_snow,
                             start_y = 1954,
                             end_y = 2014,
                             break_day = 274)


#Timing + magnitude 14-day melt max all evelation bands

yea_min_mag_14_all <- NULL
yea_min_doy_14_all <- NULL

for (i in 1:ncol(svolu_d_band)) {
  
  print(paste0(i, " out of ", ncol(svolu_d_band)))
  
  svol_stat <- svolu_d_band[, i]
  
  svol_stat_dif <- c(NA,diff(svol_stat))
  
  svol_stat_dif[which(svol_stat_dif > 0)] <- NA
  
  #Moving average filter
  svol_stat_dif_ma_14 <- rollapply(data = svol_stat_dif, width = 14,
                                   FUN = sum_na, align = "center", fill = NA)
  
  #Order data by day
  data_day_14 <- ord_day(data_in = svol_stat_dif_ma_14,
                         date = date_snow,
                         start_y = 1954,
                         end_y = 2014,
                         break_day = 274,
                         do_ma = F,
                         window_width = 30)
  
  yea_min_mag_14 <- apply(data_day_14, 1, min_na)
  
  min_doy <- function(data_in){
    
    doy_min <- which(data_in == min_na(data_in))[1]
    
    return(doy_min)
  }
  
  yea_min_doy_14 <- apply(data_day_14, 1, min_doy)
  
  
  #Collect results
  
  yea_min_mag_14_all <- cbind(yea_min_mag_14_all, yea_min_mag_14)
  yea_min_doy_14_all <- cbind(yea_min_doy_14_all, yea_min_doy_14)
  
}

#Export plot

# pdf(paste0(base_dir,"R/figs_exp/basin_flood.pdf"), width = 8, height = 4*2.6)
# tiff(paste0(base_dir,"R/figs_exp/basin_flood.tiff"), width = 8, height = 4*2.6,
#      units = "in", res = 300)
png(paste0(base_dir,"R/figs_exp/basin_flood.png"), width = 8, height = 4*2.6,
    units = "in", res = 300)

melt_plot <- function(i, head_ind = "", do_labs = F){
  
  # layout(matrix(c(1,2, 3,3
  # ),
  # 2, 2), widths=c(), heights=c())
  col_rain <- alpha("steelblue4", alpha = 0.7)
  col_rain_sol <- "black"
  col_accu_basin <- "grey35" #viridis(100)[40]
  col_melt_basin <- "grey35" #"steelblue4" #"darkgoldenrod4"
  col_dis <- "steelblue4" #"grey25"
  alpha_melt <- 0.8
  alpha_accu <- 0.6
  alpha_disc <- 0.6
  lwd_lines <- 1.5
  
  snow_acc_plot <-  data_day_acc_14/1000000
  snow_mel_plot <-  data_day_mel_14/1000000
  disc_data_plot <- grdc_day_base *24*3600/1000000
  rain_data_plot <- prec_day_ms_3/1000000
  rain_data_plot_sol <- prec_day_ms_3_sol/1000000
  # rain_data_plot <- prec_day_ms_5/1000000
  
  x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
  x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
  my_years <- 1954:2014

  #Snow accumulation/melt
  x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
  x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
  my_years <- 1954:2014
  
  par(mar = c(1.3, 4.2, 1.6, 2.8))
  par(family = "serif")
  
  snow_mel_plot[i, is.na(snow_mel_plot[i, ])] <- 0
  snow_acc_plot[i, is.na(snow_acc_plot[i, ])] <- 0
  #Snow melt/acc axis labels
  snow_labs_all <- c(-3000, -2500, -2000, -1500, -1000, -500, 500, 1000, 1500, 2000, 2500, 3000)
  snow_max <- max_na(abs(c(snow_mel_plot[i, ], snow_acc_plot[i, ])))
  snow_labs <- snow_labs_all[which(abs(snow_labs_all) < snow_max)]
  #Rain labels
  rain_labs_all <- c(500, 1000, 1500, 2000, 2500, 3000, 3500)
  rain_max <- max_na(rain_data_plot[i, ])
  rain_labs <- rain_labs_all[which(abs(rain_labs_all) < rain_max)]
  #Discharge labels
  disc_labs_all <- c(0, 50, 100, 150, 200, 300, 400, 600, 800)
  disc_max <- max_na(disc_data_plot[i, ])
  disc_labs <- disc_labs_all[which(abs(disc_labs_all) < disc_max)]
  
  #Snow volume basin
  plot(snow_acc_plot[i, ], type = "n", axes = F, lwd = 2, ylab = "", xlab = "",
       ylim = c(max_na(abs(c(snow_mel_plot[i, ], snow_acc_plot[i, ])))*-2,
                max_na(abs(c(snow_mel_plot[i, ], snow_acc_plot[i, ])))*2))
  # abline(v = x_axis_tic, col = "black", lwd = 0.8, lty = "dotted")
  lines(snow_mel_plot[i, ], col = alpha(col_melt_basin, alpha = 0.9), lwd = 1.8)
  polygon(x = c(0.99, 1:length(snow_mel_plot[i, ]), length(snow_mel_plot[i, ])+0.01), y = c(0, snow_mel_plot[i, ], 0),
          col = alpha(col_melt_basin, alpha = alpha_melt), border = F)
  lines(snow_acc_plot[i, ], col = alpha(col_accu_basin, alpha = 1.0), lwd = 1.8)
  polygon(x = c(0.99, 1:length(snow_acc_plot[i, ]), length(snow_acc_plot[i, ])+0.01),
          y = c(0, snow_acc_plot[i, ], 0), col = alpha(col_accu_basin, alpha = alpha_accu),
          border = F)
  # abline(v = c(yea_min_doy_14[i]-7, yea_min_doy_14[i]+7) , col= "black")
  abline(h = 0)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.04)#plot ticks
  axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.06, 0), cex.axis = 1.2)#plot labels
  axis(4, at = snow_labs, labels = snow_labs, mgp=c(3, 0.08, 0), tck = -0.005, cex.axis = 1.2)
  mtext(expression(paste("SWE [10"^"6", "m"^"3", " 14day"^"-1","]")), side = 4, line = 1.5, adj = 0.5, cex = 0.75)
  mtext(paste0(head_ind, "1) Precip, snowmelt and discharge"), side = 3, line = 0.11, adj = 0.0, cex = 1.1)
  mtext(paste0(my_years[i], " - ", my_years[i+1]), side = 2, line = 2.9, adj = 0.5, cex = 1.1)
  if(do_labs){
    text(140, -1290, "melt", cex = 1.2)
    text(60, 2650, "accumulation", cex = 1.2)
    text(245, 5400, "solid", cex = 1.2, col = col_rain_sol)
    text(245, 4200, "liquid", cex = 1.2, col = col_rain)
  }
  
  box()
  
  par(new = T)
  
  #Discharge Basel
  grdc_max <- which(disc_data_plot[i, ] == max_na(disc_data_plot[i, ]))
  plot(disc_data_plot[i, ], type = "n", axes = F, lwd = 2, yaxs = "i",
       ylim = c(25, max_na(disc_data_plot[i, ])*4), ylab = "", xlab = "")
  lines(disc_data_plot[i, ], col = col_dis, lwd = 1.8)
  polygon(x = c(0.99, 1:length(disc_data_plot[i, ]), length(disc_data_plot[i, ])+0.01),
          y = c(0, disc_data_plot[i, ], 0), col = alpha(col_dis, alpha = alpha_disc),
          border = F)
  points(grdc_max, disc_data_plot[i, grdc_max], cex = 2, pch = 21, col = "black", bg = NULL)
  axis(2, at = disc_labs, labels = disc_labs, mgp=c(3, 0.15, 0), tck = -0.005, cex.axis = 1.1)
  mtext(expression(paste("Q [10"^"6", "m"^"3"," 1day"^"-1", "]")), side = 2, line = 1.1, adj = 0.0, cex = 0.75)
  par(new = T)
  
  #Precipitation High Rhine
  plot(rain_data_plot[i, ], type = "h", ylim = c(max_na(rain_data_plot[i, ])*4, 0), col = col_rain,
       axes = F, yaxs = "i", lend =2, lwd = 0.7, ylab = "", xlab = "")
  lines(rain_data_plot_sol[i, ], type = "h", col = col_rain_sol, lend = 2, lwd = 0.7)
  axis(2, at = rain_labs, labels = rain_labs, mgp=c(3, 0.15, 0), tck = -0.005, cex.axis = 1.1)
  mtext(expression(paste("P [10"^"6", "m"^"3"," 3day"^"-1", "]")), side = 2, line = 1.1, adj = 1.0, cex = 0.75)
  box()
  
  rain_data_plot_sol[i, 145] / rain_data_plot[i, 145]
  
  #Snow volume basin
  par(mar = c(1.3, 2.8, 1.6, 2.8))
  
  doy_min <- which(data_day_mel_14[i, ] == min_na(data_day_mel_14[i, ])) - 7
  doy_max <- which(data_day_mel_14[i, ] == min_na(data_day_mel_14[i, ])) + 7
  
  snow_mel_plot[i, is.na(snow_mel_plot[i, ])] <- 0
  
  plot(snow_mel_plot[i, ], type = "n", axes = F, lwd = 2, yaxs = "i",
       ylim = c(min_na(snow_mel_plot[i, ])*1.2, 0))
  lines(snow_mel_plot[i, ], col = alpha(col_melt_basin, alpha = 0.9), lwd = 1.8)
  polygon(x = c(0.99, 1:length(snow_mel_plot[i, ]), length(snow_mel_plot[i, ])+0.01), y = c(0, snow_mel_plot[i, ], 0),
          col = alpha(col_melt_basin, alpha = alpha_melt), border = F)
  # abline(v = c(doy_min, doy_max), col= "orange3", lty = "solid")
  # abline(v = x_axis_tic, col = "black", lty = "dotted", lwd = 0.8)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.06)#plot ticks
  axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.06, 0), cex.axis = 1.2)#plot labels
  axis(2, mgp=c(3, 0.05, 0), tck = -0.02, cex.axis = 1.1)
  mtext(expression(paste("SWE [10"^"6", "m"^"3", " 14day"^"-1","]")), side = 2, line = 1.2, cex = 0.75, adj = 0.5)
  # mtext(paste0(my_years[i], "/",my_years[i+1]), side = 3, line = 0.02, adj = 1.0, cex = 1.0)
  mtext(paste0(head_ind, "2) Basin melt  and elevation band timing"), side = 3, line = 0.11, adj = 0.0, cex = 1.1)
  box()
  
  
  #Plot:DOY melt min bands
  par(new = T)
  par(mar = c(1.3, 2.8, 1.6, 2.8))
  
  cex_mags <- (abs(yea_min_mag_14_all[i, ]) / max_na(abs(yea_min_mag_14_all[i, ]))) * 2.5
  # elevs_ins <- my_elev_bands[which(yea_min_doy_14_all[i, ] > doy_min & yea_min_doy_14_all[i, ] < doy_max)]
  
  plot(yea_min_doy_14_all[i, ], my_elev_bands[-length(my_elev_bands)], type = "p", pch = 19, cex = cex_mags,
       col = alpha("black", alpha = 0.6), xlim = c(0,365), axes = F)
  # abline(v = c(doy_min, doy_max), col= "red3", lty = "solid")
  # abline(v = x_axis_tic, col = "black", lty = "dotted", lwd = 0.8)
  # abline(h = c(1000, 2000, 3000), col = "grey55", lty = "dotted")
  # abline(h = range(elevs_ins), col= "red3")
  # axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
  #      col = "black", col.axis = "black", tck = -0.04)#plot ticks
  # axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
  #      col="black", col.axis="black", mgp=c(3, 0.20, 0), cex.axis = 0.9)#plot labels
  axis(4, mgp=c(3, 0.03, 0), tck = -0.02, cex.axis = 1.1)
  # mtext(paste0(my_years[i], "/",my_years[i+1]), side = 3, line = 0.02, adj = 1.0, cex = 1.0)
  # mtext("2. Timing max. melt", side = 3, line = 0.02, adj = 0.5, cex = 1.0)
  mtext("Elevation [m]", side = 4, line = 1.3, cex = 0.75, adj = 0.5)
  box()
  
  #Snow stations
  
  sel_sno_dat <- function(statID, index_year){
    
    snow_stat <- snow_all[, which(colnames(snow_all) == statID)]
    
    data_day_sno <- ord_day(data_in = snow_stat,
                            date = snow_stat_date,
                            start_y = 1954,
                            end_y = 2014,
                            break_day = 274,
                            do_ma = F,
                            window_width = 30)
    
    sno_sel <- data_day_sno[index_year, ]
    
    return(sno_sel)
    
  }
  
  snow_wfj <- sel_sno_dat("WFJ", i)
  snow_aro <- sel_sno_dat("ARO", i)
  snow_sia <- sel_sno_dat("SIA", i)
  snow_dav <- sel_sno_dat("DAV", i)
  snow_gtt <- sel_sno_dat("GTT", i)
  snow_elm <- sel_sno_dat("ELM", i)
  snow_ein <- sel_sno_dat("EIN", i)
  snow_sma <- sel_sno_dat("SMA", i)
  
  #One month NA summer 2011 in WFJ
  #put to NA for polygon drawing
  snow_wfj[is.na(snow_wfj)] <- 0
  
  ssmo <- -1 #negative value is no smoothing
  
  x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
  x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
  
  par(mar = c(1.3, 2.8, 1.3, 2.8))
  
  # col_a <- viridis(100)[1]; col_a <- viridis(100)[18]
  # col_b <- viridis(100)[22]; col_b <- viridis(100)[35]
  # col_c <- viridis(100)[59]; col_c <- "grey50"
  # col_d <- viridis(100)[82]; col_d <- "darkgoldenrod4"
  # col_e <- viridis(100)[100]; col_e <- "darkred"
  
  col_a <- viridis(100)[1]; col_a <-  viridis(100)[27]
  col_b <- viridis(100)[22]; col_b <- viridis(100)[35]
  col_c <- viridis(100)[59]; col_c <- viridis(100)[45]
  col_d <- viridis(100)[82]; col_d <- "grey40"
  col_e <- viridis(100)[100]; col_e <- "gold3"
  
  plot(snow_wfj, type = "n", axes = F, ylab = "", xlab = "")
  lines(smoothFFT(snow_wfj, sd = ssmo), col = col_a, lwd = lwd_lines)
  lines(smoothFFT(snow_aro, sd = ssmo), col = alpha(col_b, alpha = 1.0), lwd = lwd_lines)
  # lines(smoothFFT(snow_sia, sd = ssmo), col = alpha("blue3", alpha = 0.7), lwd = lwd_lines)
  lines(smoothFFT(snow_dav, sd = ssmo), col = alpha(col_c, alpha = 1.0), lwd = lwd_lines)
  # lines(smoothFFT(snow_gtt, sd = ssmo), col = alpha("black", alpha = 0.7), lwd = lwd_lines)
  # lines(smoothFFT(snow_elm, sd = ssmo), col = alpha("black", alpha = 0.7), lwd = lwd_lines)
  lines(smoothFFT(snow_ein, sd = ssmo), col = alpha(col_d, alpha = 1.0), lwd = lwd_lines)
  lines(smoothFFT(snow_sma, sd = ssmo), col = alpha(col_e, alpha = 1.0), lwd = lwd_lines)
  polygon(x = c(1:365, 365:1) , y = c(smoothFFT(snow_wfj, sd = ssmo), rev(smoothFFT(snow_aro, sd = ssmo))),
          col = alpha(col_a, alpha = 0.7), border = F)
  polygon(x = c(1:365, 365:1) , y = c(smoothFFT(snow_aro, sd = ssmo), rev(smoothFFT(snow_dav, sd = ssmo))),
          col = alpha(col_b, alpha = 0.8), border = F)
  polygon(x = c(1:365, 365:1) , y = c(smoothFFT(snow_dav, sd = ssmo), rev(smoothFFT(snow_ein, sd = ssmo))),
          col = alpha(col_c, alpha = 0.8), border = F)
  polygon(x = c(1:365, 365:1) , y = c(smoothFFT(snow_ein, sd = ssmo), rev(smoothFFT(snow_sma, sd = ssmo))),
          col = alpha(col_d, alpha = 0.7), border = F)
  polygon(x = c(1:365) , y = c(smoothFFT(snow_sma, sd = ssmo)),
          col = alpha(col_e, alpha = 1.0), border = F)
  # abline(v = x_axis_tic, col = "black", lty = "dotted", lwd = 0.8)
  axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
       col = "black", col.axis = "black", tck = -0.06)#plot ticks
  axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
       col="black", col.axis="black", mgp=c(3, 0.06, 0), cex.axis = 1.1)#plot labels
  axis(2, mgp=c(3, 0.05, 0), tck = -0.02, cex.axis = 1.1)
  mtext("Depth [cm]", side = 2, line = 1.2, cex = 0.75, adj = 0.5)
  mtext(paste0(head_ind, "3) Snow observations"), side = 3, line = 0.04, adj = 0.0, cex = 1.1)
  legend("topright", c("WFJ", "ARO", "DAV", "EIN", "SMA"), col = c(col_a, col_b, col_c, col_d, col_e),
         pch = 19, cex = 0.85, bg = "white")
  box()
  
  
}

layout(matrix(c(1, 1, 13, 4, 4, 14, 7, 7, 15, 10, 10,
                2, 3, 13, 5, 6, 14, 8, 9, 15, 11, 12),
              11, 2), widths=c(), heights=c(1, 1, 0.35, 1, 1, 0.35, 1, 1, 0.35, 1, 1, 0.35, 1, 1))

# melt_plot(i = 13, head_ind = "a", do_labs = T)

melt_plot(i = 16, head_ind = "a", do_labs = T)

melt_plot(i = 25, head_ind = "b")

melt_plot(i = 27, head_ind = "c")

melt_plot(i = 34, head_ind = "d")

dev.off()




#calib_res----

# pdf(paste0(base_dir, "R/figs_exp/calib_res.pdf"), width = 16, height = 4*2.5)
png(paste0(base_dir,"R/figs_exp/calib_res.png"), width = 16, height = 4*2.5,
    units = "in", res = 300)
# tiff(paste0(base_dir,"R/figs_exp/calib_res.tiff"), width = 16, height = 4*2.5,
#      units = "in", res = 300)

par(mfrow = c(4, 1))
par(mar = c(2, 5.5, 2.7, 0.5))
par(family = "serif")

main_plots <- c("a) Andermatt (1440 m)", "b) Davos (1560 m)", "c) Zermatt (1600 m)","d) Weissfluhjoch (2540 m)", "Zermatt")

for(i in 1:ncol(snows_cal)){
  
  y_max <- max_na(c(snows_cal[, i], snows_stat_cal[, i]))
  y_min <- min_na(c(snows_cal[, i], snows_stat_cal[, i]))
  x_tics <- which(format(meteo_date_cal, '%m-%d') == '01-01')
  x_labs <- which(format(meteo_date_cal, '%m-%d') == '07-01')
  x_labels <- unique(format(meteo_date_cal, '%Y'))
  
  plot(snows_cal[, i], type = "n", ylim = c(y_min, y_max), axes = F, ylab = "", 
       xlab = "", xaxs = "i")
  abline(v = c(x_tics), lty = "dashed", col = "grey50", lwd = 0.8)
  lines(snows_cal[, i], col = scales::alpha("black", alpha = 1.0), lwd = 1.7)
  points(snows_stat_cal[, i], pch = 21, bg = scales::alpha("steelblue4", alpha = 0.7), 
         cex = 1.5, col = "steelblue")
  axis(2, cex = 1.5, mgp=c(3, 0.50, 0), cex.axis = 2.0, tck = -0.02)
  axis(1, at = x_tics, labels = rep("", length(x_tics)), tck = -0.02)
  axis(1, at = x_labs, labels = x_labels, tick = F, mgp=c(3, 0.65, 0), cex.axis = 2.0)
  mtext(main_plots[i], side= 3, line = 0.3, cex = 1.7, adj = 0.0)
  mtext("SWE [m]", side= 2, line = 2.7, cex = 1.7, adj = 0.5)
  box()
  if(i == 1){
    legend("topright", c("obs", "sim."), pch = 19, cex = 1.8, col = c("steelblue4", "black"), bg = "white")
    
  }
  
}

dev.off()

#scd_maps----

# stopCluster(clust_exp)

n_cores <- 10 #number of cores used for parallel computing

#Make cluster for parallel computing
clust_exp <- makeCluster(n_cores)
clusterEvalQ(clust_exp, pacman::p_load(alptempr, meltimr))
registerDoParallel(clust_exp)

val2col <- function(val_in, dat_ref, do_log = F, do_bicol = T, col_na = "white"){
  
  if(do_log){
    
    val_in <- log(val_in)
    dat_ref <- log(dat_ref)
    
  }
  
  if(is.na(val_in)){#set NAs to mean to keep script running; later back to NA
    val_in <- mea_na(dat_ref)
    set2NA_1 <- T
  }else{
    set2NA_1 <- F
  }
  
  if(do_bicol){
    
    col_ind <- round((abs(val_in) / max_na(abs(dat_ref))) * 100)
    
    if(val_in < 0){
      my_col  <- colorRampPalette(c("grey80", "lemonchiffon2", "lightgoldenrod2", "gold3", "goldenrod3", "orangered4", "darkred"))(100)
    }else{
      my_col  <- colorRampPalette(c("grey80", "lightcyan3", viridis::viridis(9, direction = 1)[c(4,3,2,1,1)]))(100)
    }
    
  }else{
    col_ind <- round((val_in-min_na(dat_ref)) / (max_na(dat_ref)-min_na(dat_ref)) * 200)  
    my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
  }
  
  
  if(is.na(col_ind)){
    set2NA_2 <- T
    col_ind <- 1 #set to one to keep script running; later set to NA color
  }else{
    set2NA_2 = F
  }
  
  if(col_ind == 0){#for minimum and very small values
    
    col_ind <- 1
    
  }
  
  col_out <- my_col[col_ind]
  
  if(length(col_out) < 1){
    
    col_out <- col_na
    
  }
  
  if(set2NA_1 | set2NA_2){
    
    col_out <- col_na
    
  }
  
  return(col_out)
  
}

#Values to colors simulation
sc_doy_simu_ann <- sc_doy_simu / round(length(date_vali) / 365)
cols_spat_sim <- foreach(i = 1:length(sc_doy_simu_ann), .combine = 'cbind') %dopar% {
  
  val2col(val_in = sc_doy_simu_ann[i],
          dat_ref = sc_doy_simu_ann,
          do_bicol = F)
  
}

#Values to colors difference
scd_dif <- (sc_doy_simu - scd_eurac) / round(length(date_vali) / 365) #Calculate difference Obs. and Sim.
cols_spat_dif <- foreach(i = 1:length(scd_dif), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_dif[i], 
          dat_ref = scd_dif,
          do_log = F,
          do_bicol = T)
  
}

#Values to colors observations
scd_eurac_ann <- scd_eurac / round(length(date_vali) / 365)
cols_spat_eur <- foreach(i = 1:length(scd_eurac_ann), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_eurac_ann[i],
          dat_ref = scd_eurac_ann,
          do_bicol = F)
  
}


# pdf(paste0(base_dir,"R/figs_exp/scd_maps.pdf"), width = 16, height = 3.2)
# tiff(paste0(base_dir,"R/figs_exp/scd_maps.tiff"), width = 16, height = 3.2,
#      units = "in", res = 300)
png(paste0(base_dir,"R/figs_exp/scd_maps.png"), width = 16, height = 3.2,
    units = "in", res = 300)

#Plot maps
layout(matrix(c(rep(1, 7), 2, rep(3, 7), 4, rep(5, 7), 6),
              1, 24, byrow = T), widths=c(), heights=c())
# layout.show(n = 7)

par(family = "serif")
cex_pch <- 0.27

#Map Simulations
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_base)
points(grid_points_d_in@coords[, 1], grid_points_d_in@coords[, 2], pch = 15, col = cols_spat_sim, cex = cex_pch)
plot(basin_base, add =T, lwd = 1.5)
mtext("a) Snow simulations", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(0, max_na(abs(sc_doy_simu_ann)), length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(sc_doy_simu), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

#Map Difference
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_base)
points(grid_points_d_in@coords[, 1], grid_points_d_in@coords[, 2], pch = 15, col = cols_spat_dif, cex = cex_pch)
plot(basin_base, add = T, lwd = 1.5)
mtext("b) Difference (a - c)", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
cols_min <- colorRampPalette(c("darkred", "darkorange4", "goldenrod3", "gold3", "lightgoldenrod2", "lemonchiffon2", "grey80"))(100)
cols_max <- colorRampPalette(c("grey80", "lightcyan3", viridis::viridis(9, direction = 1)[c(4,3,2,1,1)]))(100)
my_col <- colorRampPalette(c(cols_min, cols_max))(200)
my_bre <- seq(-max_na(abs(scd_dif)), max_na(abs(scd_dif)), length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_dif), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

#Map EURAC
par(mar = c(0.5, 0.5, 1.0, 0.5))
plot(basin_base)
points(grid_points_d_in@coords[, 1], grid_points_d_in@coords[, 2], pch = 15, col = cols_spat_eur, cex = cex_pch)
plot(basin_base, add =T, lwd = 1.5)
mtext("c) MODIS snow cover", side = 3, line = -1.0, cex = 1.7)

par(mar = c(2.0, 0.2, 5.0, 3.0))
my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
my_bre <- seq(0, max_na(abs(scd_eurac_ann)), length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_eurac_ann), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.50, 0), tck = -0.1, cex.axis = 1.6)
mtext(expression(paste("[", "day ", "year"^"-1", "]")), side = 3, line = 0.8, cex = 1.3)
box()

# #Header above maps
# par(mar = c(0.2, 0.2, 0.2, 0.2))
# plot(1:100, 1:100, type = "n", axes = F, ylab = "", xlab = "")
# mtext("Average annual total days with snow cover  for the time frame 08/2002 to 07/2012", 
#       side = 3, line = -1.8, cex = 1.4)

dev.off()



#melt_boxes----

#Calculate Timing and Magnitude
yea_min_mag_bands_14 <- NULL
yea_min_doy_bands_14 <- NULL

for(i in 1:ncol(svolu_d_band)){
  
  swe_band <- svolu_d_band[, i]
  
  swe_band_dif <- c(NA, diff(swe_band))
  
  swe_band_dif[which(swe_band_dif > 0)] <- 0
  
  #Moving average filter
  swe_band_dif_ma_14 <- rollapply(data = swe_band_dif, width = 30,
                                  FUN = sum_na, align = "center", fill = NA)
  
  #Order data by day
  data_day_14 <- ord_day(data_in = swe_band_dif_ma_14,
                         date = date_snow,
                         start_y = 1954,
                         end_y = 2014,
                         break_day = 274,
                         do_ma = F)
  
  min_doy <- function(data_in){
    
    doy_min <- as.numeric(which(data_in == min_na(data_in)))
    
    if(length(doy_min) > 1){
      doy_min <- mea_na(doy_min)
    }
    
    return(doy_min)
  }
  
  yea_min_mag_14 <- apply(data_day_14, 1, min_na)
  yea_min_doy_14 <- apply(data_day_14, 1, min_doy)
  
  yea_min_mag_bands_14 <- cbind(yea_min_mag_bands_14, yea_min_mag_14)
  yea_min_doy_bands_14 <- cbind(yea_min_doy_bands_14, yea_min_doy_14)
  
}

yea_min_doy_1 <- apply(yea_min_doy_bands_14[1:30, ], 2, mea_na)
yea_min_doy_2 <- apply(yea_min_doy_bands_14[31:60, ], 2, mea_na)
yea_min_mag_1 <- apply(yea_min_mag_bands_14[1:30, ], 2, mea_na) * -1
yea_min_mag_2 <- apply(yea_min_mag_bands_14[31:60, ], 2, mea_na) * -1

doy_mea <- round(apply(yea_min_doy_bands_14, 2, mea_na))
doy_mea_1 <- yea_min_doy_1
doy_mea_2 <- yea_min_doy_2
range_doy_mea <- range(doy_mea)

#Calculate melt elevation
melt_elevs_ts <- NULL
for(d in range_doy_mea[1]:range_doy_mea[2]){
  
  time_melt <- (d-14):(d+15)
  elev_collector <- NULL
  
  for(i in 1:nrow(yea_min_doy_bands_14)){
    
    elev_collector <- c(elev_collector, mea_na(my_elev_bands[which(yea_min_doy_bands_14[i, ] %in% time_melt)]))
    
  }
  
  melt_elevs_ts <- cbind(melt_elevs_ts, elev_collector)
  
}

melt_elevs_1 <- apply(melt_elevs_ts[1:30,], 2, mea_na)
melt_elevs_2 <- apply(melt_elevs_ts[31:60,], 2, mea_na)


#Plot Boxplots

pdf(paste0(base_dir,"R/figs_exp/melt_boxes_raw.pdf"), width = 6, height = 3)

par(mfrow = c(1, 3))
par(mar = c(0.3, 6, 1.0, 1))
par(family = "serif")
box_col <- "grey85"

doy_ranges <- range(doy_mea_2-doy_mea_1)
boxplot(doy_mea_2-doy_mea_1, col = box_col, axes = F, lwd = 1.8, outpch = 19,
        ylim = c(doy_ranges[2:1]), whisklty = "solid", medlwd = 2.5, medcol = "black")
axis(2, mgp=c(3, 0.25, 0), tck = -0.011, cex.axis = 1.8)
mtext("Timing change [day]", side = 2, line = 2.4, cex = 1.6)
box()

boxplot((yea_min_mag_2 - yea_min_mag_1)/1000000, col = box_col, axes = F, lwd = 1.8,
        ylim = range((yea_min_mag_2 - yea_min_mag_1)/1000000)[c(2,1)], 
        whisklty = "solid", medlwd = 2.5, medcol = "black", outpch = 19)
axis(2, mgp=c(3, 0.25, 0), tck = -0.011, cex.axis = 1.8)
mtext(expression(paste("Magnitude change [10"^"6","m"^"3","]")), side = 2, line = 2.4, cex = 1.6)
box()

boxplot(melt_elevs_2-melt_elevs_1, col = box_col, axes = F, lwd = 1.8,
        whisklty = "solid", medlwd = 2.5, medcol = "black", outpch = 19)
axis(2, mgp=c(3, 0.25, 0), tck = -0.011, cex.axis = 1.8)
mtext("Elevation change [m]", side = 2, line = 2.4, cex = 1.6)
box()

dev.off()



#Raster graph: Elevation band

# pdf("/home/rottler/ownCloud/RhineFlow/rhine_snow/manus/meltim_v1/figures/raster_elev_band.pdf", width = 6, height = 3)
tiff(paste0(base_dir,"R/figs_exp/raster_elev_band.tiff"), width = 6, height = 3,
     units = "in", res = 800)

band_sel <- 31 #my_elev_bands[31]

snow_band <- svolu_d_band[, band_sel]/ 1000000 # [hm³]

data_day <- ord_day(data_in = snow_band,
                    date = date_snow,
                    start_y = 1954,
                    end_y = 2014,
                    break_day = 274)

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15

# cols_hydro <- c(colorRampPalette(c("white", viridis::viridis(20, direction = -1)))(200))
cols_hydro <- colorRampPalette(c("grey95", viridis::viridis(9, direction = 1)[4:1]))(200)

breaks_hydro <- seq(alptempr::min_na(data_day), alptempr::max_na(data_day), length.out = length(cols_hydro)+1)

layout(matrix(c(rep(1, 8), 2),
              1, 9), widths=c(), heights=c())

par(mar = c(1.5, 3.5, 2.5, 0.2))
par(family = "serif")

image(x = 1:ncol(data_day),
      y = 1954:2013,
      z = t(data_day),
      col = cols_hydro,
      breaks = breaks_hydro,
      ylab = "", xlab = "", axes = F)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.05)#plot ticks
axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.50, 0), cex.axis = 1.7)#plot labels
axis(2, mgp=c(3, 0.12, 0), tck = -0.01, cex.axis = 1.5)
abline(h = 1983.5, col = "black", lwd = 2)
# mtext("Year", side = 2, line = 2.1, cex = 1.2)
mtext("a) SWE volume elevation band", side = 3, line = 0.3, cex = 1.4, adj = 0.0)
mtext("[hm³]", side = 3, line = 0.3, cex = 1.2, adj = 1.0)
box()

par(mar = c(1.5, 0.2, 2.5, 2.2))

alptempr::image_scale(as.matrix(data_day), col = cols_hydro, breaks = breaks_hydro, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.40, 0), tck = -0.10, cex.axis = 1.3)
box()

dev.off()







#band_asso----

#Bands calculate
yea_min_mag_bands_14 <- NULL
yea_min_doy_bands_14 <- NULL

for(i in 1:ncol(svolu_d_band)){
  
  swe_band <- svolu_d_band[, i]
  
  swe_band_dif <- c(NA, diff(swe_band))
  
  swe_band_dif[which(swe_band_dif > 0)] <- 0
  
  #Moving average filter
  swe_band_dif_ma_14 <- rollapply(data = swe_band_dif, width = 14,
                                  FUN = sum_na, align = "center", fill = NA)
  
  #Order data by day
  data_day_14 <- ord_day(data_in = swe_band_dif_ma_14,
                         date = date_snow,
                         start_y = 1954,
                         end_y = 2014,
                         break_day = 274,
                         do_ma = F,
                         window_width = 30)
  
  min_doy <- function(data_in){
    
    doy_min <- as.numeric(which(data_in == min_na(data_in)))
    
    if(length(doy_min) > 1){
      doy_min <- mea_na(doy_min)
    }
    
    return(doy_min)
  }
  
  yea_min_mag_14 <- apply(data_day_14, 1, min_na)
  yea_min_doy_14 <- apply(data_day_14, 1, min_doy)
  
  yea_min_mag_bands_14 <- cbind(yea_min_mag_bands_14, yea_min_mag_14)
  yea_min_doy_bands_14 <- cbind(yea_min_doy_bands_14, yea_min_doy_14)
  
}

bands_doy_mea <- apply(yea_min_doy_bands_14, 2, mea_na)
bands_doy_slo <- apply(yea_min_doy_bands_14, 2, sens_slo)*-1
bands_doy_mea_1 <- apply(yea_min_doy_bands_14[1:30,], 2, med_na)
bands_doy_mea_2 <- apply(yea_min_doy_bands_14[31:60,], 2, med_na)

freq_bands_1 <- NULL
freq_bands_2 <- NULL
freq_diff_all <- NULL
length_bands_1 <- NULL
length_bands_2 <- NULL
width_over <- 14

for(b in 1:ncol(svolu_d_band)){
  
  band_sel <- b
  
  bands_frequ_1 <- NULL
  band_lengt_1 <- NULL
  for(y in 1:30){
    
    mid_day <- yea_min_doy_bands_14[y, band_sel]
    min_day <- mid_day - width_over
    max_day <- mid_day + width_over
    melt_window <- min_day:max_day
    
    bands_asso_1 <- which(yea_min_doy_bands_14[y, ] > min_day & yea_min_doy_bands_14[y, ] < max_day)
    
    bands_frequ_1 <- c(bands_frequ_1, bands_asso_1)
    band_lengt_1 <- c(band_lengt_1, length(bands_asso_1))
    
  }
  hist_1 <- hist(bands_frequ_1, breaks = seq(from = 0.5, to = ncol(svolu_d_band)+0.5, by = 1), plot = F)
  
  bands_frequ_2 <- NULL
  band_lengt_2 <- NULL
  for(x in 31:60){
    
    # print(y)
    mid_day <- yea_min_doy_bands_14[x, band_sel]
    min_day <- mid_day - width_over
    max_day <- mid_day + width_over
    melt_window <- min_day:max_day
    
    bands_asso_2 <- which(yea_min_doy_bands_14[x, ] > min_day & yea_min_doy_bands_14[x, ] < max_day)
    
    bands_frequ_2 <- c(bands_frequ_2, bands_asso_2)
    band_lengt_2 <- c(band_lengt_2, length(bands_asso_2))
    
  }
  hist_2 <- hist(bands_frequ_2, breaks = seq(from = 0.5, to = ncol(svolu_d_band)+0.5, by = 1), plot = F)
  
  freq_diff <- hist_2$counts - hist_1$counts
  
  freq_bands_1 <- cbind(freq_bands_1, hist_1$counts)
  freq_bands_2 <- cbind(freq_bands_2, hist_2$counts)
  freq_diff_all <- cbind(freq_diff_all, freq_diff)
  length_bands_1 <- cbind(length_bands_1, band_lengt_1)
  length_bands_2 <- cbind(length_bands_2, band_lengt_2)
}

length_bands_1_mea <- apply(length_bands_1, 2, mea_na)
length_bands_2_mea <- apply(length_bands_2, 2, mea_na)

#Plot Frequencies individual band
pdf(paste0(base_dir,"R/figs_exp/bands_asso_1.pdf"), width = 8, height = 5)

band_test <- 38 #my_elev_bands[38]
col_1 <- viridis(9, direction = 1)[4]
col_2 <- "darkred"

par(mar = c(2.5, 2.5, 2, 0.2))

plot(freq_bands_1[, band_test], type = "h", col = alpha(col_1, alpha = 0.6), lwd = 7, ylab = "", xlab = "",
     axes = F, lend = 2, xaxs = "i", yaxs = "i", ylim = c(0, 32), xlim = c(0, 60))
legend("topleft", c("1985-2014", "1954-1984"), col = c(col_2, col_1), pch = 19)
par(new = T)
plot(freq_bands_2[, band_test], type = "h", col = alpha(col_2, alpha = 0.6), lwd = 7, ylab = "", xlab = "",
     axes = F, lend = 2, xaxs = "i", yaxs = "i", ylim = c(0, 32), xlim = c(0, 60))
par(new = T)
plot(band_test, freq_bands_2[band_test, band_test], type = "h", col = alpha("black", alpha = 1.0), lwd = 7, ylab = "", xlab = "",
     axes = F, lend = 2, xaxs = "i", yaxs = "i", ylim = c(0, 32), xlim = c(0, 60))
elevs_sel <- c(6, 16, 26, 36, 46, 56)
axis(1, at = elevs_sel, labels = my_elev_bands[elevs_sel], mgp=c(3, 0.15, 0), tck = -0.01, cex.axis = 0.9)
axis(2, mgp=c(3, 0.15, 0), tck = -0.01, cex.axis = 0.9)
mtext("Elevation band [m]", side = 1, adj = 0.5, line = 1.3, cex = 1.1)
mtext("Frequency concurrent melt [-]", side = 2, adj = 0.5, line = 1.3, cex = 1.1)
mtext("a) Frequency concurrent melt (elevation band 2050-2100 m)", side = 3, adj = 0.0, line = 0.2, cex = 1.4)
box()

dev.off()



pdf(paste0(base_dir,"R/figs_exp/bands_asso_2.pdf"), width = 8, height = 5)

par(mar = c(2.5, 2.5, 2.0, 0.2))

plot(freq_diff_all[band_test, ], type = "n", axes = F, ylab = "", xlab = "", xaxs = "i", yaxs = "i",
     xlim = c(0, 60), ylim = c(-12.9, 9))
abline(v = band_test, col = "black", lwd = 1.5)
lines(freq_diff_all[band_test, ], lwd = 2, col = "grey35")
polygon(x = c(0.999, 1:58, 58.0001) , y = c(0, freq_diff_all[band_test, ], 0), col = alpha("grey55", alpha = 0.7), border = F)
elevs_sel <- c(6, 16, 26, 36, 46, 56)
axis(1, at = elevs_sel, labels = my_elev_bands[elevs_sel], mgp=c(3, 0.15, 0), tck = -0.01, cex.axis = 0.9)
axis(2, mgp=c(3, 0.15, 0), tck = -0.01, cex.axis = 0.9)
# abline(h = 0, lty = "dashed")
mtext("Elevation band [m]", side = 1, adj = 0.5, line = 1.3, cex = 1.1)
mtext("Change frequency concurrent melt [-]", side = 2, adj = 0.5, line = 1.3, cex = 1.1)
mtext("b) Changes in concurrent melt (elevation band 2050-2100 m)", side = 3, adj = 0.0, line = 0.2, cex = 1.4)
box()

dev.off()


pdf(paste0(base_dir,"R/figs_exp/bands_asso_3.pdf"), width = 8, height = 5)

cols_min <- colorRampPalette(c(viridis::viridis(9, direction = 1)[1:4], "cadetblue3", "grey90"))(50)
cols_max <- colorRampPalette(c("grey90", "gold3",  "orange3", "orangered4", "orangered4", "darkred"))(50)
cols_mel <- c(cols_min, cols_max)

layout(matrix(c(rep(1, 9), 2),
              1, 10), widths=c(), heights=c())

par(mar = c(3.5, 3.5, 2.5, 0.2))

image(x = 1:nrow(freq_diff_all),
      y = 1:ncol(freq_diff_all),
      z = freq_diff_all,
      col = cols_mel, breaks = seq(from = -12, to = 12, length.out = 101), axes = F, ylab = "", xlab = "")
elevs_sel <- c(6, 16, 26, 36, 46, 56)
axis(1, at = elevs_sel, labels = my_elev_bands[elevs_sel], mgp=c(3, 0.55, 0), tck = -0.005, cex.axis = 1.4)
axis(2, at = elevs_sel, labels = my_elev_bands[elevs_sel], mgp=c(3, 0.25, 0), tck = -0.005, cex.axis = 1.4)
mtext("Elevation band [m]", side = 1, adj = 0.5, line = 2.1, cex = 1.1)
mtext("Elevation band [m]", side = 2, adj = 0.5, line = 2.1, cex = 1.1)
mtext("c) Changes in concurrent melt (all elevation bands)", side = 3, adj = 0.0, line = 0.2, cex = 1.4)
mtext("[-]", side = 3, adj = 1.0, line = 0.2, cex = 1.1)
box()

par(mar = c(3.5, 0.2, 2.5, 2.0))

alptempr::image_scale(as.matrix(freq_diff_all), col = cols_mel, breaks = seq(from = -12, to = 12, length.out = 101), horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.35, 0), tck = -0.08, cex.axis = 1.4)
# mtext(lab_unit, side = 3, line = 0.3, cex = 1)
box()

dev.off()





#wtc_change----


#Calculate changes in frequencies of weather types

gwt26_data <- read.table("U:/rhine_snow/data/idaweb/order75953/order_75953_data.txt",
                         sep = ";", skip = 1, header = TRUE, na.strings = "-")

gwt26_data$time <- as.POSIXct(strptime(gwt26_data$time, "%Y%m%d", tz="UTC"))

start_year <- 1959; start_day <- "1959-01-01"
end_year <- 2018; end_day   <- "2018-12-31"

start_date <- as.POSIXct(strptime(start_day, "%Y-%m-%d", tz="UTC"))
end_date   <- as.POSIXct(strptime(end_day,   "%Y-%m-%d", tz="UTC"))
full_date  <- seq(start_date, end_date, by="day")

data_gwt26 <- data.frame(date = full_date,
                         value = with(gwt26_data, gwt26_data$wkwtg3d0[match(full_date, time)]))

gwt_low_tem0 <- c(1:8, 25)
gwt_high_tem0 <- c(9:16, 26)
window_width <- 30
cover_thres <- 0.50

sens_slope <- function(data_in, cover_thresh = 0.9){
  
  if(length(which(is.na(data_in))) / length(data_in) > (1-cover_thresh)){
    sens_slo <-  NA
  }else{
    time_step <- 1:length(data_in)
    sens_slo <- as.numeric(zyp.sen(data_in~time_step)$coefficients[2])
    #sens_slo <- as.numeric(zyp.trend.vector(data_in, method = "zhang", conf.intervals = F)[2])
  }
  return(sens_slo)
}

moving_analys_wt <- function (dates, values, start_year, end_year, window_width,
                              cover_thresh, method_analys, weather_type = 1) {
  
  f_weatherYN <- function(data_in, gwts = weather_type) {
    if (data_in %in% gwts) {
      data_in <- 1
    }
    else {
      data_in <- 0
    }
    return(data_in)
  }
  
  valuesYN <- sapply(values, f_weatherYN)
  valuesYN_mw <- rollapply(data = valuesYN,
                           width = window_width, FUN = mea_na, align = "center",
                           fill = NA)
  
  data_day <- ord_day(data_in = valuesYN_mw,
                      date = dates,
                      start_y = start_year,
                      end_y = end_year,
                      break_day = 274,
                      do_ma = F,
                      window_width = 30)
  
  f_sens_slope <- function(data_in) {
    sens_slope(data_in = data_in, cover_thresh = cover_thresh)
  }
  
  mov_res <- apply(data_day[, -1], 2, f_sens_slope)
  
  for (i in 1:365) {
    if ((length(which(is.na(data_day[, i])))/nrow(data_day)) < (1 - cover_thresh)) {
      if (length(which(data_day[, i] == 1)) >= ((nrow(data_day) - 2) - length(which(is.na(data_day[, i]))))) {
        mov_res[i] <- 0
      }
      if (length(which(data_day[, i] == 0)) >=  ((nrow(data_day) - 2) - length(which(is.na(data_day[, i]))))) {
        mov_res[i] <- 0
      }
    }
  }
  
  return(mov_res)
}

gwt_tem0_high <- moving_analys_wt(dates = data_gwt26$date, values = data_gwt26$value, start_year = start_year,
                                  end_year = end_year, window_width = window_width,
                                  cover_thresh= cover_thres, method_analys = "weather_type_window_likeli_sens_slope",
                                  weather_type = gwt_high_tem0)*100*10# [%/dec]

gwt_tem0_low  <- moving_analys_wt(dates = data_gwt26$date, values = data_gwt26$value, start_year = start_year,
                                  end_year = end_year, window_width = window_width,
                                  cover_thresh= cover_thres, method_analys = "weather_type_window_likeli_sens_slope",
                                  weather_type = gwt_low_tem0)*100*10 # [%/dec]


#Calculate melt elevation changes

#Calculate melt elevation
melt_elevs_ts <- NULL
for(d in range_doy_mea[1]:range_doy_mea[2]){
  
  time_melt <- (d-14):(d+15)
  elev_collector <- NULL

  for(i in 1:nrow(yea_min_doy_bands_14)){
    
    elev_collector <- c(elev_collector, mea_na(my_elev_bands[which(yea_min_doy_bands_14[i, ] %in% time_melt)]))

  }
  
  melt_elevs_ts <- cbind(melt_elevs_ts, elev_collector)
  
}

melt_elevs_1 <- apply(melt_elevs_ts[1:30,], 2, mea_na)
melt_elevs_2 <- apply(melt_elevs_ts[31:60,], 2, mea_na)


#Plot: Weather types, temperature and elevation shift

pdf(paste0(base_dir, "R/figs_exp/temp_wt_elev.pdf"), width = 8, height = 4.5)

par(mfrow = c(2, 1))
par(family = "serif")

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15

days_sel <- 1:364

par(mar = c(2.0, 3.0, 2.0, 3.0))

plot(x = 1:364, y = tslo_band[days_sel, 28], type = "n", axes = F, ylim = range(tslo_band), ylab = "", xlab = "")
for(i in 1:ncol(tslo_band)){
  lines(x = days_sel, y = tslo_band[days_sel, i], col = alpha("black", alpha = 0.2))
}
abline(h = 0, lty = "dashed", col = "grey55", lwd = 1)
axis(2, mgp=c(3, 0.15, 0), tck = -0.01, cex.axis = 1.0)
abline(v = x_axis_tic, lty = "dashed", col = "grey55", lwd = 1)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
axis(1, at = x_axis_lab, c("O","N","D","J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.05, 0), cex.axis = 1.0)#plot labels
mtext("a) Temperature trends and changes in melt elevation", side = 3, cex = 1.2, line = 0.3, adj = 0)
mtext("Temp. [°C/dec]", side = 2, cex = 1.0, line = 1.5)

par(new = T)

col_elev <- "darkblue"
ele_range <- c(-max_na(abs(melt_elevs_2 - melt_elevs_1)), max_na(abs(melt_elevs_2 - melt_elevs_1)))
plot(range_doy_mea[1]:range_doy_mea[2], melt_elevs_2 - melt_elevs_1, type = "l", axes = F, lwd = 2,
     xlim = c(1,365), ylim = ele_range, col = alpha(col_elev, alpha = 0.7), ylab = "", xlab = "")
points(range_doy_mea[1]:range_doy_mea[2], melt_elevs_2 - melt_elevs_1, pch = 19, cex = 0.3, 
       col = alpha(col_elev, alpha = 0.7))
axis(4, mgp=c(3, 0.15, 0), tck = -0.03, cex.axis = 1.0, col = col_elev, col.axis = col_elev)
mtext("Elev. change [m]", side = 4, cex = 1.0, line = 1.4, col = col_elev)
box()


#Plot: WTC frequencies
col_1 <- "steelblue4"
col_2 <- "darkred"

days_sel <- 1:365
plot(gwt_tem0_high[days_sel], type = "l", col = col_2, ylim = c(range(c(gwt_tem0_high-gwt_tem0_low))), 
     axes = F, lwd = 2, ylab = "", xlab = "")
lines(gwt_tem0_low[days_sel], col = col_1, lwd = 2)
lines(gwt_tem0_high[days_sel] - gwt_tem0_low[days_sel], lwd = 2)
abline(h = 0, lty = "dashed", col = "grey55", lwd = 1)
abline(v = x_axis_tic, lty = "dashed", col = "grey55", lwd = 1)
axis(2, mgp=c(3, 0.15, 0), tck = -0.01, cex.axis = 1.0)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
axis(1, at = x_axis_lab, c("O","N","D","J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.05, 0), cex.axis = 1.0)#plot labels
mtext("b) Changes in warm/cold weather types", side = 3, cex = 1.2, line = 0.3, adj = 0)
mtext("Window prob. [%/dec]", side = 2, cex = 1.0, line = 1.5)
legend("topleft", c("              ","     "), box.lwd = 0, box.col = "white",bg = "white", cex = 0.8)
mtext("Warm GWTs", side = 3, cex = 0.75, line = -1.0, adj = 0.01, col = col_2)
mtext("Cold GWTs", side = 3, cex = 0.75, line = -1.7, adj = 0.01, col = col_1)
mtext("WTE index", side = 3, cex = 0.75, line = -2.4, adj = 0.01)
box()

dev.off()




#melt_basin_old----

pdf(paste0(base_dir,"R/figs_exp/melt_basin.pdf"), width = 16, height = 3.5)
# tiff("/home/rottler/ownCloud/RhineFlow/rhine_snow/manus/meltim_v1/figures/melt_ext.tiff", width = 16, height = 10,
#      units = "in", res = 800)


col_melt_14 <- viridis(9, direction = 1)[4]
col_melt_7 <- "black"
par(family = "serif")

layout(matrix(c(rep(1, 8), 2, rep(3, 9), rep(4, 9)),
              1, 27), widths=c(), heights=c())

#Raster graph SWE volume basin
sno_vol_basin_h <- sno_vol_basin / 1000000 # [hm³]

data_day <- ord_day(data_in = sno_vol_basin_h,
                    date = date_snow,
                    start_y = 1954,
                    end_y = 2014,
                    break_day = 274)

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15

# cols_hydro <- c(colorRampPalette(c("white", viridis::viridis(20, direction = -1)))(200))
cols_hydro <- colorRampPalette(c("grey95", viridis::viridis(9, direction = 1)[4:1]))(200)

breaks_hydro <- seq(alptempr::min_na(data_day), alptempr::max_na(data_day), length.out = length(cols_hydro)+1)

par(mar = c(1.8, 4.5, 2.5, 0.2))

image(x = 1:ncol(data_day),
      y = 1954:2013,
      z = t(data_day),
      col = cols_hydro,
      breaks = breaks_hydro,
      ylab = "", xlab = "", axes = F)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.04)#plot ticks
axis(1, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.52, 0), cex.axis = 1.6)#plot labels
axis(2, mgp=c(3, 0.25, 0), tck = -0.01, cex.axis = 1.6)
mtext("Year", side = 2, line = 2.0, cex = 1.4)
mtext("a) SWE volume basin", side = 3, line = 0.3, cex = 1.6, adj = 0.0)
mtext("[hm³]", side = 3, line = 0.3, cex = 1.2, adj = 1.0)

box()

par(mar = c(1.5, 0.2, 2.5, 2.2))

alptempr::image_scale(as.matrix(data_day), col = cols_hydro, breaks = breaks_hydro, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.40, 0), tck = -0.08, cex.axis = 1.2)
box()


#Basin calculate
sno_vol_basin_h <- sno_vol_basin / 1000000 # [hm³]

sno_vol_basin_h_dif <- c(NA, diff(sno_vol_basin_h))

sno_vol_basin_h_dif[which(sno_vol_basin_h_dif > 0)] <- 0

#Moving average filter
snow_basin_dif_ma_7 <- rollapply(data = sno_vol_basin_h_dif, width = 7,
                                 FUN = sum_na, align = "center", fill = NA)
snow_basin_dif_ma_14 <- rollapply(data = sno_vol_basin_h_dif, width = 14,
                                  FUN = sum_na, align = "center", fill = NA)

#Order data by day
data_day_7 <- ord_day(data_in = snow_basin_dif_ma_7,
                      date = date_snow,
                      start_y = 1954,
                      end_y = 2014,
                      break_day = 274,
                      do_ma = F,
                      window_width = 30)

data_day_14 <- ord_day(data_in = snow_basin_dif_ma_14,
                       date = date_snow,
                       start_y = 1954,
                       end_y = 2014,
                       break_day = 274,
                       do_ma = F,
                       window_width = 30)

yea_min_mag_7 <- apply(data_day_7, 1, min_na)* -1 #melt positive value
yea_min_mag_14 <- apply(data_day_14, 1, min_na)* -1 #melt positive value

min_doy <- function(data_in){
  
  doy_min <- which(data_in == min_na(data_in))[1]
  
  return(doy_min)
}

yea_min_doy_7 <- apply(data_day_7, 1, min_doy)
yea_min_doy_14 <- apply(data_day_14, 1, min_doy)

slo_mag_7  <- sens_slo(yea_min_mag_7)
int_mag_7  <- as.numeric(zyp.trend.vector(yea_min_mag_7, x = 1:length(yea_min_mag_7),
                                          method = "zhang", conf.intervals = F)[11])
sig_mag_7 <- as.numeric(zyp.trend.vector(yea_min_mag_7, x = 1:length(yea_min_mag_7),
                                         method = "zhang", conf.intervals = F)[6])

slo_mag_14 <- sens_slo(yea_min_mag_14)
int_mag_14 <- as.numeric(zyp.trend.vector(yea_min_mag_14, x = 1:length(yea_min_mag_14),
                                          method = "zhang", conf.intervals = F)[11])
sig_mag_14 <- as.numeric(zyp.trend.vector(yea_min_mag_14, x = 1:length(yea_min_mag_14),
                                          method = "zhang", conf.intervals = F)[6])

slo_doy_7 <- sens_slo(yea_min_doy_7)
int_doy_7 <- as.numeric(zyp.trend.vector(yea_min_doy_7, x = 1:length(yea_min_doy_7),
                                         method = "zhang", conf.intervals = F)[11])
sig_doy_7 <- as.numeric(zyp.trend.vector(yea_min_doy_7, x = 1:length(yea_min_doy_7),
                                         method = "zhang", conf.intervals = F)[6])

slo_doy_14 <- sens_slo(yea_min_doy_14)
int_doy_14 <- as.numeric(zyp.trend.vector(yea_min_doy_14, x = 1:length(yea_min_doy_14),
                                          method = "zhang", conf.intervals = F)[11])
sig_doy_14 <- as.numeric(zyp.trend.vector(yea_min_doy_14, x = 1:length(yea_min_doy_14),
                                          method = "zhang", conf.intervals = F)[6])


#DOY basin

par(mar = c(1.8, 5.5, 2.5, 0.2))

plot(1:length(yea_min_doy_7), yea_min_doy_7, type = "n", axes = F,
     ylab = "", xlab = "", ylim = c(90, 250))
lines(1:length(yea_min_doy_7), yea_min_doy_7, type = "l", col = col_melt_7, lwd = 1.8)
lines(1:length(yea_min_doy_14), yea_min_doy_14, type = "l", col = col_melt_14, lwd = 1.8)
segments(x0 = 1, y0 = (1*slo_doy_7+int_doy_7), x1 = 59, y1 =(59*slo_doy_7+int_doy_7),
         col = col_melt_7, lwd = 1.2, lty = "dashed")
segments(x0 = 1, y0 = (1*slo_doy_14+int_doy_14), x1 = 59, y1 =(59*slo_doy_14+int_doy_14),
         col = col_melt_14, lwd = 1.2, lty = "dashed")
abline(h = x_axis_tic, col = "grey65", lwd = 0.7, lty = "dotted")
abline(v =c(7, 17, 27, 37, 47, 57), col = "grey65", lwd = 0.7, lty = "dotted")
axis(1, at = c(7, 17, 27, 37, 47, 57) , labels = c(1960, 1970, 1980, 1990, 2000, 2010),
     mgp=c(3, 0.60, 0), tck = -0.01, cex.axis = 1.6)
x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
axis(2, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.04)#plot ticks
axis(2, at = x_axis_lab, c("O", "N", "D", "J","F","M","A","M","J","J","A","S"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.25, 0), cex.axis = 1.6)#plot labels
mtext("b) Timing basin", side = 3, line = 0.3, adj = 0.0, cex = 1.6)
mtext("Month", side = 2, line = 2.0, adj = 0.5, cex = 1.4)
mtext(paste0(round(slo_doy_7, digits = 2)*10,  " day/dec"), col = col_melt_7, side = 3, line = 1.1, cex = 0.9, adj = 0.75)
mtext(paste0(round(slo_doy_14, digits = 2)*10, " day/dec"), col = col_melt_14, side = 3, line = 0.01, cex = 0.9, adj = 0.75)
mtext(paste0("p = ", round(sig_doy_7, digits = 2)), col = col_melt_7, side = 3, line = 1.1, cex = 0.9, adj = 1.0)
mtext(paste0("p = ", round(sig_doy_14, digits = 2)), col = col_melt_14, side = 3, line = 0.01, cex = 0.9, adj = 1.0)
legend("bottomleft", c("7-day sum", "14-day sum"), pch = 19, col = c(col_melt_7, col_melt_14), cex = 1.5, bg = "white")
box()


#Magnitude Basin
ylims <- range(yea_min_mag_14, yea_min_mag_7)

plot(1:length(yea_min_mag_14), yea_min_mag_14, type = "n", ylim = ylims, col = viridis(9, direction = 1)[4], axes = F,
     ylab = "", xlab = "")
lines(1:length(yea_min_mag_7), yea_min_mag_7, type = "l", ylim = ylims, col = col_melt_7, lwd = 1.8)
lines(1:length(yea_min_mag_14), yea_min_mag_14, type = "l", ylim = ylims, col = col_melt_14, lwd = 1.8)
segments(x0 = 1, y0 = (1*slo_mag_7+int_mag_7), x1 = 59, y1 =(59*slo_mag_7+int_mag_7),
         col = col_melt_7, lwd = 1.2, lty = "dashed")
segments(x0 = 1, y0 = (1*slo_mag_14+int_mag_14), x1 = 59, y1 =(59*slo_mag_14+int_mag_14),
         col = col_melt_14, lwd = 1.2, lty = "dashed")
axis(1, at = c(7, 17, 27, 37, 47, 57) ,labels = c(1960, 1970, 1980, 1990, 2000, 2010),
     mgp=c(3, 0.60, 0), tck = -0.01, cex.axis = 1.6)
axis(2, mgp=c(3, 0.25, 0), tck = -0.01, cex.axis = 1.6)
abline(v = c(7, 17, 27, 37, 47, 57), col = "grey65", lwd = 0.7, lty = "dotted")
abline(h = c(1000, 2000, 3000), col = "grey65", lwd = 0.7, lty = "dotted")
mtext("c) Magnitude basin", side = 3, line = 0.3, adj = 0.0, cex = 1.6)
mtext("Melt volume [hm³]", side = 2, line = 2.0, adj = 0.5, cex = 1.4)
mtext(paste0(round(slo_mag_7, digits = 2)*10,  " hm³/dec"), col = col_melt_7, side = 3, line = 1.1, cex = 0.9, adj = 0.75)
mtext(paste0(round(slo_mag_14, digits = 2)*10, " hm³/dec"), col = col_melt_14, side = 3, line = 0.01, cex = 0.9, adj = 0.75)
mtext(paste0("p = ", round(sig_mag_7, digits = 2)), col = col_melt_7, side = 3, line = 1.1, cex = 0.9, adj = 1.0)
mtext(paste0("p = ", round(sig_mag_14, digits = 2)), col = col_melt_14, side = 3, line = 0.01, cex = 0.9, adj = 1.0)
legend("topright", c("7-day sum", "14-day sum"), pch = 19, col = c(col_melt_7, col_melt_14), cex = 1.5, bg = "white")
box()

dev.off()


#melt_basin_new----

pdf(paste0(base_dir,"R/figs_exp/melt_basin.pdf"), width = 8, height = 3.0)
# tiff("/home/rottler/ownCloud/RhineFlow/rhine_snow/manus/meltim_v1/figures/melt_ext.tiff", width = 16, height = 10,
#      units = "in", res = 800)

col_melt_14 <- "black"
par(family = "serif")
par(mfrow = c(1, 2))

#Basin calculate
sno_vol_basin_h <- sno_vol_basin / 1000000 # [hm³]

sno_vol_basin_h_dif <- c(NA, diff(sno_vol_basin_h))

sno_vol_basin_h_dif[which(sno_vol_basin_h_dif > 0)] <- 0

#Moving average filter
snow_basin_dif_ma_14 <- rollapply(data = sno_vol_basin_h_dif, width = 14,
                                  FUN = sum_na, align = "center", fill = NA)

data_day_14 <- ord_day(data_in = snow_basin_dif_ma_14,
                       date = date_snow,
                       start_y = 1954,
                       end_y = 2014,
                       break_day = 274,
                       do_ma = F,
                       window_width = 30)

yea_min_mag_14 <- apply(data_day_14, 1, min_na)* -1 #melt positive value

min_doy <- function(data_in){
  
  doy_min <- which(data_in == min_na(data_in))[1]
  
  return(doy_min)
}

yea_min_doy_14 <- apply(data_day_14, 1, min_doy)

slo_mag_14 <- sens_slo(yea_min_mag_14)
int_mag_14 <- as.numeric(zyp.trend.vector(yea_min_mag_14, x = 1:length(yea_min_mag_14),
                                          method = "zhang", conf.intervals = F)[11])
sig_mag_14 <- as.numeric(zyp.trend.vector(yea_min_mag_14, x = 1:length(yea_min_mag_14),
                                          method = "zhang", conf.intervals = F)[6])


slo_doy_14 <- sens_slo(yea_min_doy_14)
int_doy_14 <- as.numeric(zyp.trend.vector(yea_min_doy_14, x = 1:length(yea_min_doy_14),
                                          method = "zhang", conf.intervals = F)[11])
sig_doy_14 <- as.numeric(zyp.trend.vector(yea_min_doy_14, x = 1:length(yea_min_doy_14),
                                          method = "zhang", conf.intervals = F)[6])


#DOY basin

par(mar = c(1.8, 3.5, 2.2, 0.2))

plot(1:length(yea_min_doy_14), yea_min_doy_14, type = "n", axes = F,
     ylab = "", xlab = "", ylim = c(93, 273), yaxs = "i")
abline(h = x_axis_tic, col = "grey65", lwd = 0.9, lty = "dotted")
abline(v =c(7, 17, 27, 37, 47, 57), col = "grey65", lwd = 0.9, lty = "dotted")
lines(1:length(yea_min_doy_14), yea_min_doy_14, type = "l", col = col_melt_14, lwd = 1.8)
segments(x0 = 1, y0 = (1*slo_doy_14+int_doy_14), x1 = 59, y1 =(59*slo_doy_14+int_doy_14),
         col = col_melt_14, lwd = 1.2, lty = "dashed")
axis(1, at = c(7, 17, 27, 37, 47, 57) , labels = c(1960, 1970, 1980, 1990, 2000, 2010),
     mgp=c(3, 0.40, 0), tck = -0.02, cex.axis = 1.3)
x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)+92
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15+92
axis(2, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.04)#plot ticks
axis(2, at = x_axis_lab, c("J","F","M","A","M","J","J","A","S", "O", "N", "D"), tick = FALSE,
     col="black", col.axis="black", mgp=c(3, 0.25, 0), cex.axis = 1.2)#plot labels
mtext("Timing basin", side = 3, line = 0.3, adj = 0.0, cex = 1.6)
mtext("Month", side = 2, line = 1.5, adj = 0.5, cex = 1.6)
mtext(paste(round(slo_doy_14, digits = 2)*10, "day/decade"), col = col_melt_14, side = 3, line = 1.00, cex = 1.0, adj = 1.0)
mtext(paste0("p-value = ", round(sig_doy_14, digits = 2)),           col = col_melt_14, side = 3, line = 0.01, cex = 1.0, adj = 1.0)
box()


#Magnitude Basin
ylims <- range(yea_min_mag_14)

plot(1:length(yea_min_mag_14), yea_min_mag_14, type = "n", ylim = ylims, col = viridis(9, direction = 1)[4], axes = F,
     ylab = "", xlab = "")
abline(v = c(7, 17, 27, 37, 47, 57), col = "grey65", lwd = 0.9, lty = "dotted")
abline(h = c(1000, 2000, 3000), col = "grey65", lwd = 0.9, lty = "dotted")
lines(1:length(yea_min_mag_14), yea_min_mag_14, type = "l", ylim = ylims, col = col_melt_14, lwd = 1.8)
segments(x0 = 1, y0 = (1*slo_mag_14+int_mag_14), x1 = 59, y1 =(59*slo_mag_14+int_mag_14),
         col = col_melt_14, lwd = 1.2, lty = "dashed")
axis(1, at = c(7, 17, 27, 37, 47, 57) ,labels = c(1960, 1970, 1980, 1990, 2000, 2010),
     mgp=c(3, 0.40, 0), tck = -0.02, cex.axis = 1.3)
axis(2, mgp=c(3, 0.25, 0), tck = -0.01, cex.axis = 1.2)
mtext("Magnitude basin", side = 3, line = 0.3, adj = 0.0, cex = 1.6)
mtext(expression(paste("Melt volume [10"^"6","m"^"3","]")), side = 2, line = 1.5, adj = 0.5, cex = 1.6)
mtext(paste(round(slo_mag_14, digits = 2)*10, " 106", "m3", " decade-1"), col = col_melt_14, side = 3, line = 1.00, cex = 1.0, adj = 1.0)
mtext(paste0("p-value = ", round(sig_mag_14, digits = 2)),              col = col_melt_14, side = 3, line = 0.01, cex = 1.0, adj = 1.0)
box()

dev.off()



#flood_gene----

grdc_base_full <- read_grdc(paste0(grdc_dir, "6935051_Q_Day.Cmd.txt"))

#clip 1951-2014
sta_ind <- which(grdc_base_full$date == "1951-01-01")
end_ind <- which(grdc_base_full$date == "2014-12-31")

grdc_base <- grdc_base_full[sta_ind:end_ind, ]

#Peak over threshold discharge
pot_thre_base <- quantile(grdc_base$value, 0.95, na.rm = T)

pot_data_base <- data.frame(obs = grdc_base$value,
                              time = grdc_base$date)

pot_peaks_base <- clust(data = pot_data_base, u = pot_thre_base, tim.cond = 14, clust.max = T, plot = F)


#14d snowmelt
sno_vol_dif <- c(NA, diff(sno_vol_basin))

snow_vol_dif_mel <- sno_vol_dif
snow_vol_dif_mel[which(snow_vol_dif_mel > 0)] <- NA

sno_vol_dif_mel_ma_14 <- rollapply(data = snow_vol_dif_mel, width = 14,
                                   FUN = sum_na, align = "right", fill = NA)

#Basin precipitation
prec_basin_ms_3 <- rollapply(data = prec_basin, width = 3,
                             FUN = sum_na, align = "right", fill = NA)

prec_basin_liq_ms_3 <- rollapply(data = prec_basin_liq, width = 3,
                             FUN = sum_na, align = "right", fill = NA)

melt_peak <- sno_vol_dif_mel_ma_14[pot_peaks_base[, 3]] * -1 #melt is positive
prec_peak_liq <- prec_basin_liq_ms_3[pot_peaks_base[, 3]]
prec_peak <- prec_basin_ms_3[pot_peaks_base[, 3]]

peak_frac_mel <- melt_peak / (melt_peak + prec_peak_liq)

mel_thres_1 <- 0.25
mel_thres_2 <- 0.50
mel_thres_3 <- 0.75

firs_ind <- which(grdc_base$date[pot_peaks_base[, 3]] <= "1982-12-31")
seco_ind <- which(grdc_base$date[pot_peaks_base[, 3]] > "1982-12-31")
firs_len <- length(firs_ind)
seco_len <- length(seco_ind)

firs_len_mel_1 <- length(which(peak_frac_mel[firs_ind] > mel_thres_1))
seco_len_mel_1 <- length(which(peak_frac_mel[seco_ind] > mel_thres_1))
firs_len_mel_2 <- length(which(peak_frac_mel[firs_ind] > mel_thres_2))
seco_len_mel_2 <- length(which(peak_frac_mel[seco_ind] > mel_thres_2))
firs_len_mel_3 <- length(which(peak_frac_mel[firs_ind] > mel_thres_3))
seco_len_mel_3 <- length(which(peak_frac_mel[seco_ind] > mel_thres_3))

#between February and June
firs_ind_fj <- which(grdc_base$date[pot_peaks_base[, 3]] <= "1982-12-31" &
                     as.numeric(format(grdc_base$date[pot_peaks_base[, 3]], "%m")) %in% c(2:6))
seco_ind_fj <- which(grdc_base$date[pot_peaks_base[, 3]] > "1982-12-31"&
                     as.numeric(format(grdc_base$date[pot_peaks_base[, 3]], "%m")) %in% c(2:6))
firs_len_fj <- length(firs_ind_fj)
seco_len_fj <- length(seco_ind_fj)

firs_len_mel_1_fj <- length(which(peak_frac_mel[firs_ind_fj] > mel_thres_1))
seco_len_mel_1_fj <- length(which(peak_frac_mel[seco_ind_fj] > mel_thres_1))
firs_len_mel_2_fj <- length(which(peak_frac_mel[firs_ind_fj] > mel_thres_2))
seco_len_mel_2_fj <- length(which(peak_frac_mel[seco_ind_fj] > mel_thres_2))
firs_len_mel_3_fj <- length(which(peak_frac_mel[firs_ind_fj] > mel_thres_3))
seco_len_mel_3_fj <- length(which(peak_frac_mel[seco_ind_fj] > mel_thres_3))


pdf(paste0(base_dir, "R/figs_exp/flood_gene_raw.pdf"), width = 7.5, height = 3.0)

col_1 <- alpha("black", alpha = 0.8)
col_2 <- "grey65"
col_3 <- "lightskyblue3"
col_4 <- "steelblue4"

par(family = "serif")
par(mar = c(1, 3.5, 1, 1))
lwd_bar <- 55

plot(1:10, 1:10, xlim = c(0.0, 6.5), ylim = c(0, max(c(firs_len, seco_len)+10)), type = "n",
     axes = F, ylab = "", xlab = "", xaxs = "i", yaxs = "i")
lines(c(1,2.3), c(firs_len, seco_len),             type = "h", col = col_1, lwd = lwd_bar, lend = 1)
lines(c(1,2.3), c(firs_len_mel_1, seco_len_mel_1), type = "h", col = col_2, lwd = lwd_bar, lend = 1)
lines(c(1,2.3), c(firs_len_mel_2, seco_len_mel_2), type = "h", col = col_3, lwd = lwd_bar, lend = 1)
lines(c(1,2.3), c(firs_len_mel_3, seco_len_mel_3), type = "h", col = col_4, lwd = lwd_bar, lend = 1)
lines(c(4,5.3), c(firs_len_fj, seco_len_fj),       type = "h", col = col_1, lwd = lwd_bar, lend = 1)
lines(c(4,5.3), c(firs_len_mel_1_fj, seco_len_mel_1_fj), type = "h", col = col_2, lwd = lwd_bar, lend = 1)
lines(c(4,5.3), c(firs_len_mel_2_fj, seco_len_mel_2_fj), type = "h", col = col_3, lwd = lwd_bar, lend = 1)
lines(c(4,5.3), c(firs_len_mel_3_fj, seco_len_mel_3_fj), type = "h", col = col_4, lwd = lwd_bar, lend = 1)
axis(2,  mgp=c(3, 0.40, 0), tck = -0.02, cex.axis = 1.5)
mtext("Frequency [-]", side = 2, line = 2.2, cex = 1.5)
box()

dev.off()


firs_len; firs_len - firs_len_mel_1; firs_len_mel_1 - firs_len_mel_2; firs_len_mel_2 - firs_len_mel_3; firs_len_mel_3
seco_len; seco_len - seco_len_mel_1; seco_len_mel_1 - seco_len_mel_2; seco_len_mel_2 - seco_len_mel_3; seco_len_mel_3 

firs_len_fj; firs_len_fj - firs_len_mel_1_fj; firs_len_mel_1_fj - firs_len_mel_2_fj; firs_len_mel_2_fj - firs_len_mel_3_fj; firs_len_mel_3_fj
seco_len_fj; seco_len_fj - seco_len_mel_1_fj; seco_len_mel_1_fj - seco_len_mel_2_fj; seco_len_mel_2_fj - seco_len_mel_3_fj; seco_len_mel_3_fj 

#POT magnitudes

# max(melt_peak/1000000, na.rm = T)
# max(prec_peak_liq/1000000, na.rm = T)

pdf(paste0(base_dir, "R/figs_exp/flood_mags_raw.pdf"), width = 7.5, height = 3.0)

par(mfrow = c(1, 4))
par(family = "serif")
box_width <- 1.2

par(mar = c(1, 4.0, 2, 0))

boxplot(melt_peak[firs_ind]/1000000, ylim = c(0, 4000), width = box_width,
        notch = T, axes = F, col = "grey35", pch = 19)
axis(2,  mgp=c(3, 0.60, 0), tck = -0.02, cex.axis = 2.5)

par(mar = c(1, 0, 2, 4.0))

boxplot(melt_peak[seco_ind]/1000000, ylim = c(0, 4000), width = box_width, 
        notch = T, axes = F, col = "grey35", pch = 19)

par(mar = c(1, 4.0, 2, 0))

boxplot(prec_peak_liq[firs_ind]/1000000, ylim = c(0, 4000), width = box_width,
        notch = T, axes = F, col = alpha("steelblue4", alpha = 0.8), pch = 19)

par(mar = c(1, 0, 2, 4.0))

boxplot(prec_peak_liq[seco_ind]/1000000, ylim = c(0, 4000), width = box_width,
        notch = T, axes = F, col = alpha("steelblue4", alpha = 0.8), pch = 19)
axis(4,  mgp=c(3, 1.15, 0), tck = -0.02, cex.axis = 2.5)

dev.off()




#Simple calculation protective effect
prec_basin_day <- ord_day(data_in = prec_basin_ms_3,
                          date = date_snow,
                          start_y = 1951,
                          end_y = 2014,
                          break_day = 274)

prec_basin_liq_day <- ord_day(data_in = prec_basin_liq_ms_3,
                          date = date_snow,
                          start_y = 1951,
                          end_y = 2014,
                          break_day = 274)

prec_basin_max <- apply(prec_basin_day[, 1:182], 1, max_na)
prec_basin_liq_max <- apply(prec_basin_liq_day[, 1:182], 1, max_na)

prot_eff <- (1-(prec_basin_liq_max / prec_basin_max))*100

hist(prot_eff)
mea_na(prot_eff)


#reservoirs----

#High Rhine reservoirs
#Wildenhahn und Klaholz 1996: Gro?e Speicherseen im Einzugsgebiet des Rheins, Internationale Kommission fur die Hydrologie des Rheingebietes (KHR)
# https://www.chr-khr.org/de/veroffentlichung/grosse-speicherseen-im-einzugsgebiet-des-rheins?position=16&list=zjXE9vdtnkjt146rdN7VhvwPTUHgowCdFNWGTSyOXtU


#Vorderrhein
name_vr <- c("Cumera", "Nalps", "Santa Maria", "Runcahez", "Zervreila", "Egschi")
volu_vr <- c(40.80, 44.50, 67.00, 0.44, 100.00, 0.40)
year_vr <- c(1966, 1962, 1968, 1961, 1957, 1949)
rese_vr <- cbind(name_vr, volu_vr, year_vr)

#Hinterrehin
name_hr <- c("Surner", "Valle di Lei", "B?reburg", "Marmorera", "Davoser See", "Solis", "Isel")
volu_hr <- c(18.30, 197.00, 1.00, 60.00, 11.30, 1.46, 0.3) 
year_hr <- c(1962, 1961, 1960, 1954, 1925, 1985, 1969)
rese_hr <- cbind(name_hr, volu_hr, year_hr)

#Tamina
name_ta <- c("Gigerwald", "Mapprag")
volu_ta <- c(33.40, 5.10)
year_ta <- c(1976, 1976)
rese_ta <- cbind(name_ta, volu_ta, year_ta)

volu_tot <- sum_na(volu_vr) + sum_na(volu_hr) + sum_na(volu_ta)

grdc_data_diep <- read_grdc(paste0(grdc_dir, "6935500_Q_Day.Cmd.txt"))

dis_diep_day <- ord_day(data_in = grdc_data_diep$value,
                        date = grdc_data_diep$date,
                        start_y = 1927,
                        end_y = 2016)

dis_diep_ann <- apply(dis_diep_day, 1, sum_na)
dis_diep_tot <- mea_na(dis_diep_ann)*3600*24/1000000 #hm?

volu_tot / dis_diep_tot * 100


#Thur
name_tu <- c("Seealpsee")
volu_tu <- c(0.60)
year_tu <- c(1905)
rese_tu <- cbind(name_tu, volu_tu, year_tu)

#Aare
name_aa <- c("Sanetsch", "Amensee", "Rossiniere", "Lac d'Hogrin", "Lessoc", "Montsalvens", "Rossens",  "Lac de Perolles",
             "Schiffenen", "Oberaar", "Tr?btensee", "Totensee", "Grimsel", "Raetherichboden", "Gelmer", "Mattenalp", "Engstlensee",
             "Wohlensee", "Stausee Niederried")
volu_aa <- c(2.70, 10.30, 1.70, 52.10, 0.75, 11.00, 180.00, 0.30, 35.50, 56.00, 1.00, 2.50, 98.70, 25.00, 13.40, 2.00, 2.00, 1.60, 0.40)
year_aa <- c(1965, 1942, 1972, 1968, 1973, 1920, 1947, 1870, 1963, 1953, 1950, 1950, 1932, 1950, 1929, 1950, 1961, 1920, 1913)
rese_aa <- cbind(name_aa, volu_aa, year_aa)

#Reuss
name_re <- c("Lucendro", "Oberalpsee", "G?scheneralpsee", "Schl?ttli (Seglis)", "Bannalpsee", "Lungernsee", "Wichelsee")
volu_re <- c(25.00, 0.83, 75.00, 0.35, 1.63, 50.00, 0.38)
year_re <- c(1947, 1963, 1960, 1965, 1976, 1921, 1957)
rese_re <- cbind(name_re, volu_re, year_re)

#Limmat
name_li <- c("Limmern", "Garichte", "Kloental", "Chapfensee", "Murgtal", "W?gitall/Schr?h", "Rempen", "Sihlsee", "Wettingen")
volu_li <- c(92.00, 2.90, 39.80, 0.50, 1.20, 80.30, 0.36, 91.80, 6.00)
year_li <- c(1963, 1931, 1910, 1948, 1925, 1924, 1924, 1936, 1933)
rese_li <- cbind(name_li, volu_li, year_li)

volu_tot <- sum_na(volu_aa) + sum_na(volu_re) + sum_na(volu_li)

grdc_data_unte <- read_grdc(paste0(grdc_dir, "6935300_Q_Day.Cmd.txt"))

dis_unte_day <- ord_day(data_in = grdc_data_unte$value,
                        date = grdc_data_unte$date,
                        start_y = 1927,
                        end_y = 2016)

dis_unte_ann <- apply(dis_diep_day, 1, sum_na)
dis_unte_tot <- mea_na(dis_diep_ann)*3600*24/1000000 #hm?

volu_tot / dis_diep_tot * 100



#All resevoirs
rese_all_ <- rbind(rese_aa, rese_bo, rese_br, rese_hr, rese_il, 
                  rese_li, rese_or, rese_re, rese_ta, rese_tu, rese_vr)
reses <- data.frame(name = rese_all[, 1],
                    volu = as.numeric(rese_all[, 2]),
                    year = as.numeric(rese_all[, 3]))

reses_agg <- aggregate(reses$volu, by = list(years = reses$year), FUN = sum)
reses_agg$cum <- cumsum(reses_agg$x)

perc_vol <- c(
  round(reses_agg$cum[which(reses_agg$years == 1929)] / max(cumsum(reses_agg$x)), 2)*100,
  round(reses_agg$cum[which(reses_agg$years == 1936)] / max(cumsum(reses_agg$x)), 2)*100,
  round(reses_agg$cum[which(reses_agg$years == 1950)] / max(cumsum(reses_agg$x)), 2)*100,
  round(reses_agg$cum[which(reses_agg$years == 1960)] / max(cumsum(reses_agg$x)), 2)*100,
  round(reses_agg$cum[which(reses_agg$years == 1969)] / max(cumsum(reses_agg$x)), 2)*100
)

tota_vol <- c(
  reses_agg$cum[which(reses_agg$years == 1929)],
  reses_agg$cum[which(reses_agg$years == 1936)],
  reses_agg$cum[which(reses_agg$years == 1950)],
  reses_agg$cum[which(reses_agg$years == 1960)],
  reses_agg$cum[which(reses_agg$years == 1969)]
)





#liq_sol----

liq_frac_prec <- prec_basin_liq / prec_basin

liq_prec_day <- ord_day(data_in = prec_basin_liq,
                        date = date_snow,
                        start_y = 1951,
                        end_y = 2014,
                        break_day = 274,
                        do_ma = T,
                        window_width = 3)

liq_prec_mea <- apply(liq_prec_day, 2, mea_na) * 100 #[%]
liq_prec_slo <- apply(liq_prec_day, 2, sens_slo) * 100 *10 #[%/decade]

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
lab_months <- c("O", "N", "D", "J","F","M","A","M","J","J","A","S")

par(mfrow = c(2, 1))

plot(liq_prec_mea, type = "l", axes = F, ylab = "", xlab = "")
abline(h = 0)
axis(2, mgp=c(3, 0.12, 0), tck = -0.01, cex.axis = 1.5)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
for(i in 1:length(x_axis_lab)){
  axis(1, at = x_axis_lab[i], lab_months[i], tick = FALSE, col="black", col.axis="black", 
       mgp=c(4, 0.45, 0), cex.axis = 1.7)
}
# mtext(expression(paste("Temp. [?C", "decades"^"-1", "]")), side = 2, line = 1.6, cex = 1.5)
box()


plot(liq_prec_slo, type = "l", axes = F, ylab = "", xlab = "")

abline(h = 0)
axis(2, mgp=c(3, 0.12, 0), tck = -0.01, cex.axis = 1.5)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
for(i in 1:length(x_axis_lab)){
  axis(1, at = x_axis_lab[i], lab_months[i], tick = FALSE, col="black", col.axis="black", 
       mgp=c(4, 0.45, 0), cex.axis = 1.7)
}
# mtext(expression(paste("Temp. [?C", "decades"^"-1", "]")), side = 2, line = 1.6, cex = 1.5)
box()







liq_area_day <- ord_day(data_in = liq_area,
                        date = date_snow,
                        start_y = 1951,
                        end_y = 2014,
                        break_day = 274,
                        do_ma = T,
                        window_width = 30)

liq_area_mea <- apply(liq_area_day, 2, mea_na) * 100 #[%]
liq_area_slo <- apply(liq_area_day, 2, sens_slo) * 100 *10 #[%/decade]

x_axis_lab <- c(16,46,74,105,135,166,196,227,258,288,319,349)
x_axis_tic <- c(16,46,74,105,135,166,196,227,258,288,319,349,380)-15
lab_months <- c("O", "N", "D", "J","F","M","A","M","J","J","A","S")

par(mfrow = c(2, 1))

plot(liq_area_mea, type = "l", axes = F, ylab = "", xlab = "", ylim = c(0, 100))
abline(h = 100)
axis(2, mgp=c(3, 0.12, 0), tck = -0.01, cex.axis = 1.5)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
for(i in 1:length(x_axis_lab)){
  axis(1, at = x_axis_lab[i], lab_months[i], tick = FALSE, col="black", col.axis="black", 
       mgp=c(4, 0.45, 0), cex.axis = 1.7)
}
# mtext(expression(paste("Temp. [?C", "decades"^"-1", "]")), side = 2, line = 1.6, cex = 1.5)
box()


plot(liq_area_slo, type = "l", axes = F, ylab = "", xlab = "")

abline(h = 0)
axis(2, mgp=c(3, 0.12, 0), tck = -0.01, cex.axis = 1.5)
axis(1, at = x_axis_tic, c("","","","","","","","","","","","",""), tick = TRUE,
     col = "black", col.axis = "black", tck = -0.06)#plot ticks
for(i in 1:length(x_axis_lab)){
  axis(1, at = x_axis_lab[i], lab_months[i], tick = FALSE, col="black", col.axis="black", 
       mgp=c(4, 0.45, 0), cex.axis = 1.7)
}
# mtext(expression(paste("Temp. [?C", "decades"^"-1", "]")), side = 2, line = 1.6, cex = 1.5)
box()







#weatype_calc----

start_year <- 1959
end_year <- 2017

start_day <- "1959-01-01"
end_day   <- "2017-12-31"

wtc_sel <- "gwt26_msl" #gwt26_msl, gwt26_500, cap27_msl

my_wtc_sel <- 1:26

# base_dir <- "/home/rottler/ownCloud/RhineFlow/rhine_obs/"
# source(paste0(base_dir, "R/rhine_flow/2_rhine_functions.R"))

cli_stats <- c("BAS", "BER", "SMA")

for(i in 1:3){
  
  cli_sel <- cli_stats[i]
  #sel_varies
  
  
  if(wtc_sel == "gwt26_msl"){
    
    #Weather type classification data: GWT26 based on mean sea level pressure
    wtc_data <- read.table("U:/rhine_snow/data/idaweb/order64866/order_64866_data.txt",
                           sep = ";", skip = 1, header = TRUE, na.strings = "-")
    wtc_data$date <- as.POSIXct(strptime(wtc_data$time, "%Y%m%d", tz="UTC"))
    wtc_data$valu <- wtc_data$wkwtp3d0
    
  }
  
  if(wtc_sel == "gwt26_500"){
    
    #Weather type classification data: GWT26 based on mean sea level pressure
    wtc_data <- read.table("U:/rhine_snow/data/idaweb/order75953/order_75953_data.txt",
                           sep = ";", skip = 1, header = TRUE, na.strings = "-")
    wtc_data$date <- as.POSIXct(strptime(wtc_data$time, "%Y%m%d", tz="UTC"))
    wtc_data$valu <- wtc_data$wkwtg3d0
    
  }
  
  if(wtc_sel == "cap27_msl"){
    
    #Weather type classification data: GWT26 based on mean sea level pressure
    wtc_data <- read.table(paste0(base_dir, "data/idaweb/order64973/order_64973_data.txt"),
                           sep = ";", skip = 1, header = TRUE, na.strings = "-")
    wtc_data$date <- as.POSIXct(strptime(wtc_data$time, "%Y%m%d", tz="UTC"))
    wtc_data$valu <- wtc_data$wkcap3d0
    
  }
  
  if(cli_sel == "BAS"){
    
    # Basel / Binningen
    data_cli <- read.table("U:/rhine_snow/data/idaweb/order64388/order_64388_data.txt",
                           sep = ";", skip = 2, header = T, na.strings = c("-"))
    data_cli$date <- as.POSIXct(strptime(data_cli$time, "%Y%m%d", tz="UTC"))
    data_cli$valu <- data_cli$rhs150d0
    
  }
  
  if(cli_sel == "BER"){
    
    # Bern
    data_cli <- read.table("U:/rhine_snow/data/idaweb/order64387/order_64387_data.txt",
                           sep = ";", skip = 2, header = T, na.strings = c("-"))
    data_cli$date <- as.POSIXct(strptime(data_cli$time, "%Y%m%d", tz="UTC"))
    data_cli$valu <- data_cli$rhs150d0
    
  }
  
  if(cli_sel == "SMA"){
    
    # Zuerich
    data_cli <- read.table("U:/rhine_snow/data/idaweb/order64389/order_64389_data.txt",
                           sep = ";", skip = 2, header = T, na.strings = c("-"))
    data_cli$date <- as.POSIXct(strptime(data_cli$time, "%Y%m%d", tz="UTC"))
    data_cli$valu <- data_cli$rhs150d0
    
  }
  
  
  #calc_wtc_cli
  
  #Analysis within-type changes weather types using climate station data
  f_wtc_cli <- function(wtc_data_in, clim_data_in, annu_analy,
                        wtc_sel = my_wtc_sel, method_analy){
    
    start_date <- as.POSIXct(strptime(start_day, "%Y-%m-%d", tz="UTC"))
    end_date   <- as.POSIXct(strptime(end_day,   "%Y-%m-%d", tz="UTC"))
    full_date  <- seq(start_date, end_date, by="day")
    
    data_wtc <- data.frame(date = full_date,
                           value = with(wtc_data_in, wtc_data_in$valu[match(full_date, date)]))
    
    data_cli <- data.frame(date = full_date,
                           value = with(clim_data_in, clim_data_in$valu[match(full_date, date)]))
    
    n_wtcs <- length(wtc_sel)
    n_year <- length(start_year:end_year)
    
    
    f_wtc_annu <- function(input_wtc, input_cli, gwt_sel, my_annu_analy = annu_analy){
      
      #Remove 29th of February
      input_wtc <- input_wtc[-which(format(input_wtc$date, "%m%d") == "0229"),]
      input_cli <- input_cli[-which(format(input_cli$date, "%m%d") == "0229"),]
      
      #Vector with the 365 days of the year
      days <- seq(as.Date('2014-01-01'), to=as.Date('2014-12-31'), by='days')
      days <- format(days,"%m-%d")
      
      #Order data by day
      data_day_wtc <-  matrix(NA, nrow = length(start_year:end_year), ncol = 366)
      data_day_cli <-  matrix(NA, nrow = length(start_year:end_year), ncol = 366)
      colnames(data_day_wtc) <- c("year", days)
      colnames(data_day_cli) <- c("year", days)
      data_day_wtc[ ,1] <- start_year:end_year
      data_day_cli[ ,1] <- start_year:end_year
      
      for(i in 0:(length(start_year:end_year)-1)) {
        
        data_day_wtc[i+1, 2:366] <- input_wtc$valu[(i*365+1):((i+1)*365)]
        data_day_cli[i+1, 2:366] <- input_cli$valu[(i*365+1):((i+1)*365)]
        
      }
      
      wtc_out <- rep(NA, nrow(data_day_cli))
      
      for(i in 1:nrow(data_day_cli)){
        
        if(my_annu_analy == "mean"){
          
          wtc_out[i] <- mea_na(data_day_cli[i, which(data_day_wtc[i, ] == gwt_sel)])
          
        }
        
        if(my_annu_analy == "sum"){
          
          wtc_out[i] <- sum_na(data_day_cli[i, which(data_day_wtc[i, ] == gwt_sel)])
          
        }
        
      }
      
      return(wtc_out)
      
    }
    
    wtc_out <- matrix(data=rep(NA, n_wtcs*n_year), ncol = n_wtcs)
    
    for(i in 1:n_wtcs){
      
      wtc_out[, i] <- f_wtc_annu(input_wtc = data_wtc,
                                 input_cli = data_cli,
                                 gwt_sel = i)
    }
    
    if(method_analy == "mean"){
      
      wtc_return <- apply(wtc_out, 2, mea_na)
      
    }
    
    if(method_analy == "sens_slope"){
      wtc_sens_slope <- function(wtc_data_in, wtc_cover = 0.1){
        sens_slo(data_in = wtc_data_in, cover_thresh = wtc_cover)
      }
      
      wtc_return <- apply(wtc_out, 2, wtc_sens_slope) * 10 # per decade
    }
    
    return(wtc_return)
    
  }
  
  #Analysis frequency and  change in frequency weather types
  f_wtc_fre <- function(wtc_data_in, wtc_sel = my_wtc_sel, method_analy){
    
    start_date <- as.POSIXct(strptime(start_day, "%Y-%m-%d", tz="UTC"))
    end_date   <- as.POSIXct(strptime(end_day,   "%Y-%m-%d", tz="UTC"))
    full_date  <- seq(start_date, end_date, by="day")
    
    data_wtc <- data.frame(date = full_date,
                           value = with(wtc_data_in, wtc_data_in$valu[match(full_date, date)]))
    
    n_wtcs <- length(wtc_sel)
    n_year <- length(start_year:end_year)
    
    f_annu_fre <- function(input_wtc, gwt_sel){
      
      #Remove 29th of February
      input_wtc <- input_wtc[-which(format(input_wtc$date, "%m%d") == "0229"),]
      
      #Vector with the 365 days of the year
      days <- seq(as.Date('2014-01-01'), to=as.Date('2014-12-31'), by='days')
      days <- format(days,"%m-%d")
      
      #Order data by day
      data_day_wtc <-  matrix(NA, nrow = length(start_year:end_year), ncol = 366)
      colnames(data_day_wtc) <- c("year", days)
      data_day_wtc[ ,1] <- start_year:end_year
      
      for(i in 0:(length(start_year:end_year)-1)) {
        
        data_day_wtc[i+1, 2:366] <- input_wtc$valu[(i*365+1):((i+1)*365)]
        
      }
      
      fre_wtc <- rep(NA, nrow(data_day_wtc))
      
      for(i in 1:nrow(data_day_wtc)){
        
        fre_wtc[i] <-length(which(data_day_wtc[i, ] == gwt_sel))
        
      }
      
      return(fre_wtc)
      
    }
    
    wtc_fre <- matrix(data=rep(NA, n_wtcs*n_year), ncol = n_wtcs)
    
    k <- 1
    for(i in wtc_sel){
      
      wtc_fre[, k] <- f_annu_fre(input_wtc = data_wtc,
                                 gwt_sel = i)
      k <- k+1
    }
    
    #plot(wtc_fre[, 24], type = "l")
    
    if(method_analy == "mean"){
      
      wtc_fre_out <- apply(wtc_fre, 2, mea_na)
      
    }
    
    if(method_analy == "sens_slope"){
      
      f_sens_slope <- function(data_in, cover_thresh = 0.9){
        
        if(length(which(is.na(data_in))) / length(data_in) > (1-cover_thresh)){
          sens_slo <-  NA
        }else{
          time_step <- 1:length(data_in)
          sens_slo <- as.numeric(zyp.sen(data_in~time_step)$coefficients[2])
          # sens_slo <- as.numeric(zyp.trend.vector(data_in, method = "zhang", conf.intervals = F)[2])
        }
        return(sens_slo)
      }
      
      wtc_fre_out <- apply(wtc_fre, 2, f_sens_slope) * 10 # per decaade
      
    }
    
    return(wtc_fre_out)
    
  }
  
  
  wtc_sum_mea <- f_wtc_cli(wtc_data_in = wtc_data,
                           clim_data_in = data_cli,
                           annu_analy = "sum",
                           wtc_sel = my_wtc_sel,
                           method_analy = "mean")
  
  wtc_sum_slo <- f_wtc_cli(wtc_data_in = wtc_data,
                           clim_data_in = data_cli,
                           annu_analy = "sum",
                           wtc_sel = my_wtc_sel,
                           method_analy = "sens_slope")
  
  wtc_mea_mea <- f_wtc_cli(wtc_data_in = wtc_data,
                           clim_data_in = data_cli,
                           annu_analy = "mean",
                           wtc_sel = my_wtc_sel,
                           method_analy = "mean")
  
  wtc_mea_slo <- f_wtc_cli(wtc_data_in = wtc_data,
                           clim_data_in = data_cli,
                           annu_analy = "mean",
                           wtc_sel = my_wtc_sel,
                           method_analy = "sens_slope")
  
  if(cli_sel == "BAS"){
    
    wtc_sum_mea_bas <- wtc_sum_mea
    wtc_sum_slo_bas <- wtc_sum_slo
    wtc_mea_mea_bas <- wtc_mea_mea
    wtc_mea_slo_bas <- wtc_mea_slo
    
  }
  
  if(cli_sel == "BER"){
    
    wtc_sum_mea_ber <- wtc_sum_mea
    wtc_sum_slo_ber <- wtc_sum_slo
    wtc_mea_mea_ber <- wtc_mea_mea
    wtc_mea_slo_ber <- wtc_mea_slo
    
  }
  
  if(cli_sel == "SMA"){
    
    wtc_sum_mea_sma <- wtc_sum_mea
    wtc_sum_slo_sma <- wtc_sum_slo
    wtc_mea_mea_sma <- wtc_mea_mea
    wtc_mea_slo_sma <- wtc_mea_slo
    
  }
  
}

#calc_wtc

wtc_fre_mea <- f_wtc_fre(wtc_data_in = wtc_data,
                         method_analy = "mean")

wtc_fre_slo <- f_wtc_fre(wtc_data_in = wtc_data,
                         method_analy = "sens_slope")

#weatype_visu----

pdf(paste0(base_dir,"R/figs_exp/weatype.pdf"), width = 14, height = 8)

layout(matrix(c(1,3,5, 2,4,6),
              3, 2), widths=c(), heights=c(1,1,1,1))
par(oma = c(0,0,0,0))

mar_1 <- c(2.2, 2.0, 2.0, 0.5)
mar_2 <- c(3.2, 2.2, 2.2, 0.5)

gap_lenght <- 2
lwd_bar <- 2.5
gaps_wtc_plot <- 0:25 * gap_lenght
size_y_labs <- 1.5
size_x_labs <- 1.5
size_main <- 1.5
line_main <- 0.3
x_lab_posi <- c(1:8, 10, 12, 14, 16, 18, 20, 22, 24, 26)

par(family = "serif")

#Plot a: Mean annual frequency of weather types

par(mar = mar_1)

my_ylim <- c(0, max_na(wtc_fre_mea) + 5)
my_xlim <- c(0.5, 26.5)

plot(wtc_fre_mea, pch = 19, type = "h", lwd = 8, lend = 1, axes = F, xlab = "", ylab = "",
     xaxs = "i", yaxs = "i", ylim = my_ylim, xlim = my_xlim)
axis(1, at = (1:27)-0.5, labels = rep("", 27), tick = TRUE,
     col="black", col.axis="black", tck=-0.04)#plot ticks
axis(1, at = x_lab_posi, labels = x_lab_posi, tick = FALSE,
     col = "black", col.axis = "black", mgp = c(3, 0.4, 0), cex.axis = size_x_labs)
axis(2, mgp = c(3, 0.3, 0), tck = -0.015, cex.axis = size_y_labs)
abline(h = 0, lty = "dashed")
abline(v = c(8.5, 16.5, 24.5), lty = "dashed")
mtext("a) Annual GWT frequency", side = 3, line = line_main, cex = size_main, adj = 0)
mtext("[day]", side = 3, line = line_main, cex = 1.2, adj = 1)
box(lwd = 1.2)

directs <- rep(c("W", "SW", "NW", "N", "NE", "E", "SE", "S"), 3)
pos_labs <- (1:26)
for (i in 1:8){
  mtext(text = directs[i], at = pos_labs[i], cex = 0.8, side = 3, line = - 1.2)
}
for (i in 9:16){
  mtext(text = directs[i], at = pos_labs[i], cex = 0.8, side = 3, line = - 1.2)
}
for (i in 17:24){
  mtext(text = directs[i], at = pos_labs[i], cex = 0.8, side = 3, line = - 1.2)
}
mtext("cyclonic",                  side = 3, line = -1.8, adj = 0.12, padj = 1, cex = 0.9)
mtext("anticyclonic",              side = 3, line = -1.8, adj = 0.45, padj = 1, cex = 0.9)
mtext("indifferent",               side = 3, line = -1.8, adj = 0.79, padj = 1, cex = 0.9)
mtext("low pressure",              side = 4, line = -3.8, adj = 0.90, padj = 1, cex = 0.9)
mtext("high pressure",             side = 4, line = -2.1, adj = 0.925, padj = 1, cex = 0.9)

#Plot b: Trend annual frequency of weather types

my_ylim <- c(min_na(wtc_fre_slo) - 0.2, max_na(wtc_fre_slo) + 0.2)
my_xlim <- c(0.5, 26.5)

plot(wtc_fre_slo, pch = 19, type = "h", lwd = 8, lend = 1, axes = F, xlab = "", ylab = "",
     xaxs = "i", yaxs = "i", ylim = my_ylim, xlim = my_xlim)
axis(1, at = (1:27)-0.5, labels = rep("", 27), tick = TRUE,
     col="black", col.axis="black", tck=-0.04)#plot ticks
axis(1, at = x_lab_posi, labels = x_lab_posi, tick = FALSE,
     col = "black", col.axis = "black", mgp = c(3, 0.4, 0), cex.axis = size_x_labs)
axis(2, mgp = c(3, 0.3, 0), tck = -0.015, cex.axis = size_y_labs)
abline(h = 0, lty = "dashed")
abline(v = c(8.5, 16.5, 24.5), lty = "dashed")
mtext("b) Trend annual GWT frequency", side = 3, line = line_main, cex = size_main, adj = 0)
mtext("[day/dec]", side = 3, line = line_main, cex = 1.2, adj = 1)
box(lwd = 1.2)

directs <- rep(c("W", "SW", "NW", "N", "NE", "E", "SE", "S"), 3)
pos_labs <- (1:26)
for (i in 1:8){
  mtext(text = directs[i], at = pos_labs[i], cex = 0.8, side = 3, line = - 1.2)
}
for (i in 9:16){
  mtext(text = directs[i], at = pos_labs[i], cex = 0.8, side = 1, line = - 1.2)
}
for (i in 17:24){
  mtext(text = directs[i], at = pos_labs[i], cex = 0.8, side = 3, line = - 1.2)
}
mtext("cyclonic",                  side = 3, line = -1.8, adj = 0.12, padj = 1, cex = 0.9)
mtext("anticyclonic",              side = 1, line = -3.2, adj = 0.45, padj = 1, cex = 0.9)
mtext("indifferent",               side = 3, line = -1.8, adj = 0.79, padj = 1, cex = 0.9)
mtext("low pressure",              side = 4, line = -3.8, adj = 0.10, padj = 1, cex = 0.9)
mtext("high pressure",             side = 4, line = -2.1, adj = 0.10, padj = 1, cex = 0.9)


# Plot c: Mean precipitation per GWT

gap_lenght <- 2
lwd_bar <- 2.5
gaps_wtc_plot <- 0:25 * gap_lenght

my_ylim <- c(0, max_na(c(wtc_mea_mea_bas, wtc_mea_mea_ber, wtc_mea_mea_sma)) + 1.2)
my_xlim <- c(-0.5,(3 * 26 + gap_lenght*25) + gap_lenght - 0.5)

plot(((1:26) * 3 - 2) + gaps_wtc_plot, wtc_mea_mea_bas, type = "h", col = "black", lwd = 3, lend = 2,
     xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)
par(new = T)
plot(((1:26) * 3 - 1) + gaps_wtc_plot, wtc_mea_mea_ber, type = "h", col = "black", lwd = 3, lend = 2,
     xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)
par(new = T)
plot(((1:26) * 3 - 0) + gaps_wtc_plot, wtc_mea_mea_sma, type = "h", col = "black", lwd = 3, lend = 2,
     xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)
# par(new = T)
# plot(((1:26) * 4 - 0) + gaps_wtc_plot, wtc_mea_mea_hoh, type = "h", col = "black", lwd = 3, lend = 2,
#      xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)
axis(1, at = c(-0.5, ((1:26) * 3 + 1.5) + gaps_wtc_plot), labels = rep("", 27), tick = TRUE,
     col="black", col.axis="black", tck=-0.04)#plot ticks

axis(1, at = (x_lab_posi * 3) + gaps_wtc_plot[x_lab_posi] -1.0, labels = x_lab_posi, tick = FALSE,
     col = "black", col.axis = "black", mgp = c(3, 0.4, 0), cex.axis = size_x_labs)
axis(2, mgp = c(3, 0.3, 0), tck = -0.015, cex.axis = size_y_labs)
abline(h = 0, lty = "dashed", lwd = 0.7)
abline(v = c(8, 16, 24) * 3 + 1.5 + gaps_wtc_plot[c(8, 16, 24)], lty = "dashed", lwd = 0.7)
mtext("c) Precipitation per occurence", side = 3, line = line_main, cex = size_main, adj = 0)
mtext("[mm/1d]", side = 3, line = line_main, cex = 1.2, adj = 1)
box(lwd = 1.2)

# mtext("Basel",     side = 3, line = -3.0 + 2.5, adj = 0.005, padj = 1, cex = 0.6)
# mtext("Bern",      side = 3, line = -3.6 + 2.5, adj = 0.016, padj = 1, cex = 0.6)
# mtext("Zuerich",   side = 3, line = -4.2 + 2.5, adj = 0.022,  padj = 1, cex = 0.6)
# # mtext("Hohenp.",   side = 3, line = -4.8 + 2.5, adj = 0.029, padj = 1, cex = 0.6)
# 
# lines(c(1,1), c(9.3, 10.5), type = "l", lwd = 0.5)
# lines(c(2,2), c(9.3, 10.1), type = "l", lwd = 0.5)
# lines(c(3,3), c(9.3, 9.7), type = "l", lwd = 0.5)
# lines(c(4,4), c(9.9, 10.00), type = "l", lwd = 0.5)


# Plot: d: Trend mean precipiation per day

my_ylim <- c(min_na(c(wtc_mea_slo_bas, wtc_mea_slo_ber, wtc_mea_slo_sma)) - 0.03,
             max_na(c(wtc_mea_slo_bas, wtc_mea_slo_ber, wtc_mea_slo_sma)) + 0.03)
my_xlim <- c(-0.5,(3 * 26 + gap_lenght*25) + gap_lenght - 0.5)

plot(((1:26) * 3 - 2) + gaps_wtc_plot, wtc_mea_slo_bas, type = "h", col = "black", lwd = 3, lend = 2,
     xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)
par(new = T)
plot(((1:26) * 3 - 1) + gaps_wtc_plot, wtc_mea_slo_ber, type = "h", col = "black", lwd = 3, lend = 2,
     xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)
par(new = T)
plot(((1:26) * 3 - 0) + gaps_wtc_plot, wtc_mea_slo_sma, type = "h", col = "black", lwd = 3, lend = 2,
     xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)
# par(new = T)
# plot(((1:26) * 4 - 0) + gaps_wtc_plot, wtc_mea_slo_hoh, type = "h", col = "black", lwd = 3, lend = 2,
#      xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)
axis(1, at = c(-0.5, ((1:26) * 3 + 1.5) + gaps_wtc_plot), labels = rep("", 27), tick = TRUE,
     col="black", col.axis="black", tck=-0.04)#plot ticks

axis(1, at = (x_lab_posi * 3) + gaps_wtc_plot[x_lab_posi] -1.0, labels = x_lab_posi, tick = FALSE,
     col = "black", col.axis = "black", mgp = c(3, 0.4, 0), cex.axis = size_x_labs)
axis(2, mgp = c(3, 0.3, 0), tck = -0.015, cex.axis = size_y_labs)
abline(h = 0, lty = "dashed", lwd = 0.7)
abline(v = c(8, 16, 24) * 3 + 1.5 + gaps_wtc_plot[c(8, 16, 24)], lty = "dashed", lwd = 0.7)
mtext("d) Trend precipiation per occur.", side = 3, line = line_main, cex = size_main, adj = 0)
mtext("[(mm/1d)/dec]", side = 3, line = line_main, cex = 1.2, adj = 1)
box(lwd = 1.2)


# Plot e: Mean annual sum

par(mar = mar_2)

my_ylim <- c(0, max_na(c(wtc_sum_mea_bas, wtc_sum_mea_ber, wtc_sum_mea_sma)) + 5)
my_xlim <- c(-0.5,(3 * 26 + gap_lenght*25) + gap_lenght - 0.5)

plot(((1:26) * 3 - 2) + gaps_wtc_plot, wtc_sum_mea_bas, type = "h", col = "black", lwd = 3, lend = 2,
     xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)
par(new = T)
plot(((1:26) * 3 - 1) + gaps_wtc_plot, wtc_sum_mea_ber, type = "h", col = "black", lwd = 3, lend = 2,
     xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)
par(new = T)
plot(((1:26) * 3 - 0) + gaps_wtc_plot, wtc_sum_mea_sma, type = "h", col = "black", lwd = 3, lend = 2,
     xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)

axis(1, at = c(-0.5, ((1:26) * 3 + 1.5) + gaps_wtc_plot), labels = rep("", 27), tick = TRUE,
     col="black", col.axis="black", tck=-0.04)#plot ticks

axis(1, at = (x_lab_posi * 3) + gaps_wtc_plot[x_lab_posi] -1.0, labels = x_lab_posi, tick = FALSE,
     col = "black", col.axis = "black", mgp = c(3, 0.4, 0), cex.axis = size_x_labs)
axis(2, mgp = c(3, 0.3, 0), tck = -0.015, cex.axis = size_y_labs)
abline(h = 0, lty = "dashed", lwd = 0.7)
abline(v = c(8, 16, 24) * 3 + 1.5 + gaps_wtc_plot[c(8, 16, 24)], lty = "dashed", lwd = 0.7)
mtext("e) Annual total precipiation", side = 3, line = line_main, cex = size_main, adj = 0)
mtext("[mm/year]", side = 3, line = line_main, cex = 1.2, adj = 1)
mtext("GWT26 weather type", side = 1, line = 2, cex = 1.2, adj = 0.5)
box(lwd = 1.2)

mtext("Basel",     side = 3, line = -2.2 - 0.8, adj = 0.005, padj = 1, cex = 0.8)
mtext("Bern",      side = 3, line = -3.1 - 0.8, adj = 0.017, padj = 1, cex = 0.8)
mtext("Zuerich",   side = 3, line = -4.1 - 0.8, adj = 0.023, padj = 1, cex = 0.8)
# mtext("Hohenp.",   side = 3, line = -4.8 - 0.8, adj = 0.029, padj = 1, cex = 0.6)

lines(c(1,1), c(75+0, 132-20), type = "l", lwd = 0.5)
lines(c(2,2), c(75+0, 122-20), type = "l", lwd = 0.5)
lines(c(3,3), c(75+0, 116-20), type = "l", lwd = 0.5)
# lines(c(4,4), c(75+0, 110-20), type = "l", lwd = 0.5)

# Plot f: Trend annual sum

my_ylim <- c(min_na(c(wtc_sum_slo_bas, wtc_sum_slo_ber, wtc_sum_slo_sma)) - 0.5,
             max_na(c(wtc_sum_slo_bas, wtc_sum_slo_ber, wtc_sum_slo_sma)) + 0.5)
my_xlim <- c(-0.5,(3 * 26 + gap_lenght*25) + gap_lenght - 0.5)

plot(((1:26) * 3 - 2) + gaps_wtc_plot, wtc_sum_slo_bas, type = "h", col = "black", lwd = 3, lend = 2,
     xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)
par(new = T)
plot(((1:26) * 3 - 1) + gaps_wtc_plot, wtc_sum_slo_ber, type = "h", col = "black", lwd = 3, lend = 2,
     xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)
par(new = T)
plot(((1:26) * 3 - 0) + gaps_wtc_plot, wtc_sum_slo_sma, type = "h", col = "black", lwd = 3, lend = 2,
     xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)
# par(new = T)
# plot(((1:26) * 3 - 0) + gaps_wtc_plot, wtc_sum_slo_hoh, type = "h", col = "black", lwd = 3, lend = 2,
#      xaxs = "i", yaxs = "i", axes = F, ylab = "", xlab = "", ylim = my_ylim, xlim = my_xlim)
axis(1, at = c(-0.5, ((1:26) * 3 + 1.5) + gaps_wtc_plot), labels = rep("", 27), tick = TRUE,
     col="black", col.axis="black", tck=-0.04)#plot ticks

axis(1, at = (x_lab_posi * 3) + gaps_wtc_plot[x_lab_posi] -1.0, labels = x_lab_posi, tick = FALSE,
     col = "black", col.axis = "black", mgp = c(3, 0.4, 0), cex.axis = size_x_labs)
axis(2, mgp = c(3, 0.3, 0), tck = -0.015, cex.axis = size_y_labs)
abline(h = 0, lty = "dashed", lwd = 0.7)
abline(v = c(8, 16, 24) * 3 + 1.5 + gaps_wtc_plot[c(8, 16, 24)], lty = "dashed", lwd = 0.7)
mtext("f) Trend annual total precipiation", side = 3, line = line_main, cex = size_main, adj = 0)
mtext("[(mm/year)/dec]", side = 3, line = line_main, cex = 1.2, adj = 1)
mtext("GWT26 weather type", side = 1, line = 2, cex = 1.2, adj = 0.5)
box(lwd = 1.2)

dev.off()


#scf_series----

pdf(paste0(base_dir, "R/figs_exp/sc_frac.pdf"), width = 8, height = 2.5)

par(mar = c(1.5, 2.0, 1.5, 0.8))
par(family = "serif")

x_labs <- 2003:2012
x_tics <- which(as.character(scf_eurac$date) %in% c(paste0(x_labs, "-01-01")) )

plot(scf_eurac$scf, type = "n", axes = F, ylim = c(0, 1), xlab = "", ylab ="")
lines(scf_eurac$scf, col = "steelblue4")
lines(scf_simu_vali$scf, col = "darkred")
axis(1, at = x_tics, mgp=c(3, 0.15, 0), labels = rep("", length(x_tics)), tick = TRUE)
axis(1, at = x_tics+182, mgp=c(3, 0.15, 0), labels = x_labs, tick = FALSE, cex.axis = 0.8)
axis(2, mgp=c(3, 0.15, 0), tck = -0.011, cex.axis = 0.7)
mtext("Snow cover fraction [-]", side = 2, line = 1.0, cex = 0.9) 
legend("topleft", c("observed", "similated"), pch = 19, cex = 0.6, col = c("steelblue4", "darkred"), bty = "n")
box()

dev.off()

ind_sel <- 366:(9*365)
plot(date[ind_sel], lapse_temp[ind_sel], type = "l")

abline(h= median(lapse_temp[ind_sel]))

       