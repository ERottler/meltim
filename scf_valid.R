###

#Snow cover fractions from MODIS data
#Erwin Rottler, University of Potsdam

###

#set_up----

#Select validation period
sta_date <- as.POSIXct(strptime("2002-08-01", "%Y-%m-%d", tz = "UTC"))
end_date <- as.POSIXct(strptime("2012-07-31", "%Y-%m-%d", tz = "UTC"))
date_vali <- as.Date(seq(sta_date, end_date, by = "day"))

#select basin
basins <-  rgdal::readOGR(dsn = paste0(ezg_dir,"ezg_kombiniert.shp"))
basin_base <- spTransform(basins[basins@data$Ort == "Basel, Rheinhalle",], CRS = crs(scf_file, asText = T))
basin_base_buf <- buffer(basin_base, width = 8000)

#functs----

scf_extr <- function(file_path, basin_in, aggr_fac = 2, snow_val = 50, provider = "DLR"){
  
  #Read file
  scf <- raster(file_path)
  
  #corp file to basin area
  scf_cro <- raster::crop(scf, extent(basin_base))
  scf_sub <- mask(scf_cro, basin_base)
  
  scf_sub <- aggregate(scf_sub, fact = aggr_fac, fun = modal, na.rm = TRUE)
  
  #get values of cells cropped
  scf_sub_val_NA <- scf_sub@data@values
  scf_sub_val <- scf_sub_val_NA[!is.na(scf_sub_val_NA)]
  
  #Calculate snow cover fraction
  scf_valu <- length(which(scf_sub_val == snow_val)) / length(scf_sub_val)  
  
  #Extract date from file name
  if(provider == "DLR"){
    
    doy <- as.numeric(substr(file_path, nchar(file_path)-6, nchar(file_path)-4))
    yea <- substr(file_path, nchar(file_path)-11, nchar(file_path)-8)
    date <- as.character(as.Date(doy, origin = paste0(yea, "-01-01")))
    
  }
  
  if(provider == "EURAC"){
    
    day <- substr(file_path, nchar(file_path)-12, nchar(file_path)-11)
    mon <- substr(file_path, nchar(file_path)-14, nchar(file_path)-13)
    yea <- substr(file_path, nchar(file_path)-18, nchar(file_path)-15)
    date <- paste0(yea, "-", mon, "-", day)
    
  }
  
  return(c(date, scf_valu))
  
}

scf_ext_dlr <- function(file_path){scf_extr(file_path = file_path, 
                                            basin_in = basin_base,
                                            aggr_fac = 2,
                                            snow_val = 50,
                                            provider = "DLR")}

scf_ext_eurac <- function(file_path){scf_extr(file_path = file_path, 
                                            basin_in = basin_base,
                                            aggr_fac = 4,
                                            snow_val = 1,
                                            provider = "EURAC")}

#calc_dlr----

file_names <- dir(path = scf_dlr_dir, recursive = T)
# file_names <- file_names[-which(file_names == "README.txt")]
file_names <- paste0(scf_dlr_dir, "/", file_names[which(nchar(file_names) == 28)])
gsp_crs <- raster(file_names[1])

basins <-  rgdal::readOGR(dsn = paste0(ezg_dir,"ezg_kombiniert.shp"))
basin_base <- spTransform(basins[basins@data$Ort == "Basel, Rheinhalle",], CRS = crs(gsp_crs, asText = T))

scf_out <- foreach(i = 1:length(file_names), .combine = 'cbind') %dopar% {
  
   scf_ext_dlr(file_names[i])
  
}

date <- as.Date(as.character(scf_out[1, ]), "%Y-%m-%d")
scf <- as.numeric(scf_out[2, ])

scf_dlr <- data.frame(date = date[order(date)],
                      scf = scf[order(date)])

#calc_eurac----

#get file names
file_names <- dir(path = scf_eurac_dir, recursive = T)
scf_file <- raster(paste0(scf_eurac_dir , file_names[1]))

#get dates from file names
f_scf_date <- function(file_path, provider = "EURAC"){
  
  #Extract date from file name
  if(provider == "DLR"){
    
    doy <- as.numeric(substr(file_path, nchar(file_path)-6, nchar(file_path)-4))
    yea <- substr(file_path, nchar(file_path)-11, nchar(file_path)-8)
    date <- as.character(as.Date(doy, origin = paste0(yea, "-01-01")))
    
  }
  
  if(provider == "EURAC"){
    
    day <- substr(file_path, nchar(file_path)-12, nchar(file_path)-11)
    mon <- substr(file_path, nchar(file_path)-14, nchar(file_path)-13)
    yea <- substr(file_path, nchar(file_path)-18, nchar(file_path)-15)
    date <- paste0(yea, "-", mon, "-", day)
    
  }
  
  return(date)
  
}

scf_dates <- foreach(i = 1:length(file_names), .combine = 'cbind') %dopar% {
  
  f_scf_date(paste0(scf_eurac_dir, file_names[i]))
  
}

scf_date <- as.Date(as.character(scf_dates[1, ]), "%Y-%m-%d")

#select files covering validation period
file_names_sel <- which(scf_date %in% date_vali)

#SCF times series for selected basin
scf_out <- foreach(i = 1:length(file_names_sel), .combine = 'cbind') %dopar% {
  
  scf_ext_eurac(paste0(scf_eurac_dir, file_names[i]))
  
}

date <- as.Date(as.character(scf_out[1, ]), "%Y-%m-%d")
scf <- as.numeric(scf_out[2, ])

scf_eurac <- data.frame(date = date[order(date)],
                        scf = scf[order(date)])

#Spatial map with DOY with snow cover for selected basin and time frame

sc_doy_extr <- function(file_path, snow_val = 1, provider = "EURAC"){
  
  #Read file
  scf <- raster(file_path)
  
  #corp file to basin area (with buffer)
  scf_cro <- raster::crop(scf, extent(basin_base_buf))
  scf_sub <- mask(scf_cro, basin_base_buf)
  
  #get values and set to 0 if not snow
  scf_val_NA <- values(scf_sub)
  scf_val_NA[which(scf_val_NA != 1)] <- 0
  
  #Extract date from file name
  if(provider == "DLR"){
    
    doy <- as.numeric(substr(file_path, nchar(file_path)-6, nchar(file_path)-4))
    yea <- substr(file_path, nchar(file_path)-11, nchar(file_path)-8)
    date <- as.character(as.Date(doy, origin = paste0(yea, "-01-01")))
    
  }
  
  if(provider == "EURAC"){
    
    day <- substr(file_path, nchar(file_path)-12, nchar(file_path)-11)
    mon <- substr(file_path, nchar(file_path)-14, nchar(file_path)-13)
    yea <- substr(file_path, nchar(file_path)-18, nchar(file_path)-15)
    date <- paste0(yea, "-", mon, "-", day)
    
  }
  
  return(scf_val_NA)
  
}

for(i in 1:length(file_names_sel)){
  
  print(i)
  
  sc_doy_sing <- sc_doy_extr(paste0(scf_eurac_dir, file_names[i]))
  
  if(i == 1){
    sc_doy_valus <- sc_doy_sing
  }else{
    sc_doy_valus <- sc_doy_valus + sc_doy_sing 
  }
  
}

sc_doy_out <- foreach(i = 1:length(file_names_sel), .combine = 'rbind') %dopar% {
  
  sc_doy_extr(paste0(scf_eurac_dir, file_names[i]))
  
}

sc_doy_valus <- apply(sc_doy_out, 2, sum_na)

#fill dummy raster with calculated data values
scf_file_cro <- raster::crop(scf_file, extent(basin_base_buf))
scf_file_sub <- mask(scf_file_cro, basin_base_buf)
scf_file_sub@data@values <- sc_doy_valus

#visu_time----

scf_simu_vali  <- scf_simu [which(scf_simu$date  %in% date_vali), ]

#Calculate metrics
my_kge <- KGE(scf_eurac$scf, scf_simu_vali$scf)
my_cor <- cor(scf_eurac$scf, scf_simu_vali$scf, method = "pearson")

pdf("/home/rottler/ownCloud/RhineFlow/rhine_snow/manus/meltim_v1/figures/sc_frac.pdf", width = 6, height = 2.5)

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
mtext("Snow cover fraction Rhine basin until gauge Basel", side = 3, line = 0.1, cex = 1.0, adj = 0.00) 
legend("topleft", c("observed", "similated"), pch = 19, cex = 0.6, col = c("steelblue4", "darkred"), bty = "n")
box()

dev.off()

#visu_space----

plot(scf_file_sub, col = viridis(365, direction = -1))

#clip validation period from simulations
snow_simu_sel <- which(date_snow %in% date_vali)
snows_d_vali <- snows_d[snow_simu_sel, ]

#calculate sum days with snow for validation period
f_snow2sc <- function(snow_in, snow_thresh = 0.02){
  
  snow_in[which(snow_in >= snow_thresh)] <- 1
  snow_in[which(snow_in <  snow_thresh)] <- 0
  snow_in[which(is.na(snow_in))] <- NA
  
  snow_days_sum <- sum(snow_in, na.rm = T)
  
  return(snow_days_sum)
  
}

sc_doy_simu <- apply(snows_d_vali, 2, f_snow2sc)

val2col <- function(val_in, dat_ref, do_log = F){
  
  if(do_log){
    
    val_in <- log(val_in)
    dat_ref <- log(dat_ref)
    
  }

  my_col <- c(colorRampPalette(c(viridis::viridis(20, direction = -1)))(200))
  
  col_ind <- round((val_in-min_na(dat_ref)) / (max_na(dat_ref)-min_na(dat_ref)) * 200)  
  
  if(is.na(col_ind)){
    set2NA <- T
    col_ind <- 1 #set to one to keep script running; later set to NA color
  }else{
    set2NA = F
  }
  
  if(col_ind == 0){#for minimum and very small values
    
    col_ind <- 1
    
  }
  
  col_out <- my_col[col_ind]
  
  if(length(col_out) < 1){
    
    col_out <- "red"
    
  }
  
  if(set2NA){
    
    col_out <- "red"
    
  }
  
  return(col_out)
}
cols_spat <- foreach(i = 1:length(sc_doy_simu), .combine = 'cbind') %dopar% {
  
  val2col(val_in = sc_doy_simu[i],
          dat_ref = sc_doy_simu)
  
}

plot(dem_sub, axes = F, legend = F,  col = colorRampPalette(c("white", "black"))(200), box = F)
plot(basin_base, add =T)
points(grid_points_d_in@coords[, 1], grid_points_d_in@coords[, 2], pch = 19, col = cols_spat, cex = 0.30)

#Resample EURAC data fit simulations
scf_file_sub_aggr <- aggregate(scf_file_sub, fact = 4, fun = modal, na.rm = TRUE)
sc_doy_eurac <- raster::extract(scf_file_sub_aggr, grid_points_d_in)

cols_spat <- foreach(i = 1:length(sc_doy_eurac), .combine = 'cbind') %dopar% {
  
  val2col(val_in = sc_doy_eurac[i],
          dat_ref = sc_doy_eurac)
  
}

plot(dem_sub, axes = F, legend = F,  col = colorRampPalette(c("white", "black"))(200), box = F)
plot(basin_base, add =T)
points(grid_points_d_in@coords[, 1], grid_points_d_in@coords[, 2], pch = 19, col = cols_spat, cex = 0.30)

#Calculate difference Obs. and Sim.
scd_dif <- sc_doy_simu - sc_doy_eurac

scd_dif[which(abs(scd_dif) > 300)] <- NA


val2col <- function(val_in, dat_ref, do_log = F, do_bicol = T, col_na = "green"){
  
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
      my_col  <- colorRampPalette(c("grey98", "gold3", "orangered4", "firebrick4", "darkred"))(100)
    }else{
      my_col  <- colorRampPalette(c("grey98", "azure3", viridis::viridis(9, direction = 1)[c(4,3,2,1)]))(100)
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

cols_spat <- foreach(i = 1:length(scd_dif), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_dif[i], 
          dat_ref = scd_dif,
          do_log = F,
          do_bicol = T)
  
}

layout(matrix(c(rep(1, 7), 2),
              1, 8), widths=c(), heights=c())
# layout.show(n = 2)

par(mar = c(2.0, 0.5, 0.5, 0.5))
# plot(dem_sub, axes = F, legend = F,  col = colorRampPalette(c("white", "black"))(200), box = F)
plot(basin_base)
points(grid_points_d_in@coords[, 1], grid_points_d_in@coords[, 2], pch = 19, col = cols_spat, cex = 0.30)
plot(basin_base, add = T)

par(mar = c(2.0, 0.2, 2.5, 3.5))
cols_min <- colorRampPalette(c("darkred", "firebrick4", "orangered4", "gold3", "grey98"))(100)
cols_max <- colorRampPalette(c("grey98", "azure3", viridis::viridis(9, direction = 1)[c(4,3,2,1)]))(100)
my_col <- colorRampPalette(c(cols_min, cols_max))(200)
my_bre <- seq(-max_na(abs(scd_dif)), max_na(abs(scd_dif)), length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_dif), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.45, 0), tck = -0.1, cex.axis = 1.5)
box()

