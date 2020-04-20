###

#Snow cover fractions from MODIS data
#Erwin Rottler, University of Potsdam

###

#set_up----

#Select validation period
sta_date <- as.POSIXct(strptime("2002-08-01", "%Y-%m-%d", tz = "UTC"))
end_date <- as.POSIXct(strptime("2014-07-31", "%Y-%m-%d", tz = "UTC"))
date_vali <- as.Date(seq(sta_date, end_date, by = "day"))

#calc_eurac----

#get file names
file_names <- dir(path = scf_eurac_dir, recursive = T)
scf_file <- raster(paste0(scf_eurac_dir , file_names[1]))

#select basin
basins <-  rgdal::readOGR(dsn = paste0(ezg_dir,"ezg_kombiniert.shp"))
basin_base <- spTransform(basins[basins@data$Ort == "Basel, Rheinhalle",], CRS = crs(scf_file, asText = T))
basin_base_buf <- buffer(basin_base, width = 8000)

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
file_names_sel <- file_names[which(scf_date %in% date_vali)]

#sum up days with snow cover for selected basin (with buffer) and time frame

f_scd_extr <- function(file_path, snow_val, provider){
  
  #Read file
  scf <- raster(file_path)
  
  #corp file to basin area (with buffer)
  scf_cro <- raster::crop(scf, extent(basin_base_buf))
  scf_sub <- mask(scf_cro, basin_base_buf)
  
  #get values and set to 0 if not snow, to 1 if snow
  scf_val_NA <- values(scf_sub)
  scf_val_NA[which(scf_val_NA != snow_val)] <- 0
  scf_val_NA[which(scf_val_NA == snow_val)] <- 1
  
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

block_size <- 1000
block_stas <- c(1, seq(block_size+1, length(file_names_sel), by = block_size))
block_ends <- c(seq(block_size, length(file_names_sel), by = block_size), length(file_names_sel))

for(b in 1:length(block_stas)){
  
  file_names_calc <- file_names_sel[block_stas[b]:block_ends[b]]
  
  print(paste(Sys.time(),"Spatial analysis: Days with snow cover", "Block:", b, "out of", length(block_stas)))
  
  scd_out <- foreach(i = 1:length(file_names_calc), .combine = 'cbind') %dopar% {
    
    f_scd_extr(file_path = paste0(scf_eurac_dir, file_names_calc[i]),
               snow_val = 1,
               provider = "EURAC")
    
  }
  
  if(b == 1){
    scd_eurac_buf_all <- scd_out
  }else{
    scd_eurac_buf_all <- cbind(scd_eurac_buf_all, scd_out)
  }
}

scd_eurac_buf_sum <- apply(scd_eurac_buf_all, 1, sum_na)
rm(scd_eurac_buf_all) #remove file after calculation as very big
gc() #colltect some garbage

#fill dummy raster with calculated snow cover fraction values
scf_buf_crop <- raster::crop(scf_file, extent(basin_base_buf))
scf_buf <- mask(scf_buf_crop, basin_base_buf)
scf_buf@data@values <- scd_eurac_buf_sum
scf_buf_aggr <- aggregate(scf_buf, fact = 4, fun = mean, na.rm = TRUE)
plot(scf_buf_aggr, col = viridis(200, direction = -1))

#get lake sufaces
block_size <- 1000
block_stas <- c(1, seq(block_size+1, length(file_names_sel), by = block_size))
block_ends <- c(seq(block_size, length(file_names_sel), by = block_size), length(file_names_sel))

for(b in 1:length(block_stas)){
  
  file_names_calc <- file_names_sel[block_stas[b]:block_ends[b]]
  
  print(paste(Sys.time(),"Spatial analysis: Lakes", "Block:", b, "out of", length(block_stas)))
  
  lak_out <- foreach(i = 1:length(file_names_calc), .combine = 'cbind') %dopar% {
    
    f_scd_extr(file_path = paste0(scf_eurac_dir, file_names_calc[i]),
               snow_val = 5,
               provider = "EURAC")
    
  }
  
  if(b == 1){
    lak_eurac_buf_all <- lak_out
  }else{
    lak_eurac_buf_all <- cbind(lak_eurac_buf_all, lak_out)
  }
  
}

lak_eurac_buf_sum <- apply(lak_eurac_buf_all, 1, sum_na)
rm(lak_eurac_buf_all) #remove file after calculation as very big
gc() #colltect some garbage

#fill dummy raster with calculated lake surfaces
lak_buf_crop <- raster::crop(scf_file, extent(basin_base_buf))
lak_buf <- mask(lak_buf_crop, basin_base_buf)
lak_buf@data@values <- lak_eurac_buf_sum
lak_buf_aggr <- aggregate(lak_buf, fact = 4, fun = mean, na.rm = TRUE)
laks_ind <- which(lak_buf_aggr@data@values > (length(date_vali)*0.90))#define lake surfaces
lak_buf_aggr@data@values <- 0
lak_buf_aggr@data@values[laks_ind] <- 100
plot(lak_buf_aggr)
plot(basin_base_buf, add = T)

#calculate glacier surfaces
gla_buf_crop <- raster::crop(scf_file, extent(basin_base_buf))
gla_buf <- mask(gla_buf_crop, basin_base_buf)
gla_buf@data@values <- scd_eurac_buf_sum
gla_buf_aggr <- aggregate(gla_buf, fact = 4, fun = mean, na.rm = TRUE)
glac_ind <- which(gla_buf_aggr@data@values > (length(date_vali)*0.90))#define lake surfaces
gla_buf_aggr@data@values <- 0
gla_buf_aggr@data@values[glac_ind] <- 100
plot(gla_buf_aggr)
plot(basin_base_buf, add = T)

#Remove lake (and glacier surfaces)

scf_buf_aggr@data@values[laks_ind] <- NA
# scf_buf_aggr@data@values[glac_ind] <- NA
plot(scf_buf_aggr, col = viridis(n=200, direction = -1))

#Get values grid points simulated
scd_eurac <- raster::extract(scf_buf_aggr, grid_points_d_in)

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
    
    col_out <- "white"
    
  }
  
  if(set2NA){
    
    col_out <- "white"
    
  }
  
  return(col_out)
}

cols_spat <- foreach(i = 1:length(scd_eurac), .combine = 'cbind') %dopar% {
  
  val2col(val_in = scd_eurac[i],
          dat_ref = scd_eurac)
  
}

plot(basin_base)
points(grid_points_d_in@coords[, 1], grid_points_d_in@coords[, 2], pch = 19, col = cols_spat, cex = 0.30)
plot(basin_base, add =T)

#SCF times series for selected basin
f_scf <- function(file_path, basin_in, basin_in_buf, aggr_fac, snow_val, provider){
  
  #Read file
  scf <- raster(file_path)
  
  #clip basin area with buffer
  scf_cro <- raster::crop(scf, extent(basin_in_buf))
  scf_buf <- mask(scf_cro, basin_in_buf)
  
  #aggregate to simulation resolution
  scf_agg <- aggregate(scf_buf, fact = aggr_fac, fun = modal, na.rm = TRUE) #modal: most frequent value in a set of values
  
  #get values for grid points simulated
  scf_values <- raster::extract(scf_agg, grid_points_d_in)
  
  #remove lakes areas (previously determined)
  scf_values[which(is.na(scd_eurac))] <- NA
  
  #Calculate snow cover fraction
  scf_out <- length(which(scf_values == snow_val)) / length(scf_values)
  
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
  
  return(c(date, scf_out))
  
}

block_size <- 1000
block_stas <- c(1, seq(block_size+1, length(file_names_sel), by = block_size))
block_ends <- c(seq(block_size, length(file_names_sel), by = block_size), length(file_names_sel))

for(b in 1:length(block_stas)){
  
  file_names_calc <- file_names_sel[block_stas[b]:block_ends[b]]
  
  print(paste(Sys.time(),"Temporal analysis snow cover fraction", "Block:", b, "out of", length(block_stas)))
  
  scf_out <- foreach(i = 1:length(file_names_calc), .combine = 'cbind') %dopar% {
    
    f_scf(file_path = paste0(scf_eurac_dir, file_names_calc[i]),
          basin_in = basin_base,
          basin_in_buf = basin_base_buf,
          aggr_fac = 4,
          snow_val = 1,
          provider = "EURAC")
    
  }
  
  if(b == 1){
    scf_out_all <- scf_out
  }else{
    scf_out_all <- cbind(scf_out_all, scf_out)
  }
  
}

date <- as.Date(as.character(scf_out_all[1, ]), "%Y-%m-%d")
scf <- as.numeric(scf_out_all[2, ])

scf_eurac <- data.frame(date = date[order(date)],
                        scf = scf[order(date)])

plot(scf_eurac, type = "l")



#visu_time----

#Snow cover fraction simulations

scf_simu_val <- NULL

for(i in 1:nrow(snows_d)){
  
  print(i)
  
  scf_dummy <- snows_d[i, ]
  
  scf_dummy[which(scf_dummy >= 0.02)] <- 1
  scf_dummy[which(scf_dummy <  0.02)] <- 0
  
  #remove lake areas
  scf_dummy[which(is.na(scd_eurac))] <- NA
  
  scf_row <- sum(scf_dummy, na.rm = T)/  length(scf_dummy)
  
  scf_simu_val <- c(scf_simu_val, scf_row)
  
}

scf_simu <- data.frame(date = date_snow,
                       scf =  scf_simu_val)
scf_simu_vali  <- scf_simu[which(scf_simu$date  %in% date_vali), ]

#Calculate metrics
my_kge <- KGE(scf_eurac$scf, scf_simu_vali$scf)
my_cor <- cor(scf_eurac$scf, scf_simu_vali$scf, method = "pearson")

pdf(paste0(base_dir, "R/figs_exp/sc_frac.pdf"), width = 6, height = 2.5)

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

#put lake sufaces to NA
sc_doy_simu[which(is.na(scd_eurac))] <- NA

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
    
    col_out <- "white"
    
  }
  
  if(set2NA){
    
    col_out <- "white"
    
  }
  
  return(col_out)
}
cols_spat <- foreach(i = 1:length(sc_doy_simu), .combine = 'cbind') %dopar% {
  
  val2col(val_in = sc_doy_simu[i],
          dat_ref = sc_doy_simu)
  
}

plot(basin_base)
points(grid_points_d_in@coords[, 1], grid_points_d_in@coords[, 2], pch = 19, col = cols_spat, cex = 0.30)
plot(basin_base, add =T)


#Calculate difference Obs. and Sim.
scd_dif <- (sc_doy_simu - scd_eurac) / round(length(date_vali) / 365)

# scd_dif[which(abs(scd_dif) > 200)] <- NA

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
plot(basin_base)
points(grid_points_d_in@coords[, 1], grid_points_d_in@coords[, 2], pch = 19, col = cols_spat, cex = 0.30)
plot(basin_base, add = T)

par(mar = c(2.0, 0.2, 2.5, 3.5))
cols_min <- colorRampPalette(c("darkred", "darkorange4", "goldenrod3", "gold3", "lightgoldenrod2", "lemonchiffon2", "grey80"))(100)
cols_max <- colorRampPalette(c("grey80", "lightcyan3", viridis::viridis(9, direction = 1)[c(4,3,2,1,1)]))(100)
my_col <- colorRampPalette(c(cols_min, cols_max))(200)
my_bre <- seq(-max_na(abs(scd_dif)), max_na(abs(scd_dif)), length.out = length(my_col)+1)
alptempr::image_scale(as.matrix(scd_dif), col = my_col, breaks = my_bre, horiz=F, ylab="", xlab="", yaxt="n", axes=F)
axis(4, mgp=c(3, 0.45, 0), tck = -0.1, cex.axis = 1.5)
box()




  par(mfrow = c(1,1))
par(mar = c(2,2,2,2))
range(scd_dif, na.rm = T)
hist(scd_dif, breaks = seq(-160, 155, 1), col = "grey52", border = F)
abline(v = c(-30, 20))

quantile(scd_dif, probs = 0.9, na.rm = T)
boxplot(scd_dif)

length(which(abs(scd_dif) < 28)) / length(scd_dif)

plot(elevs_d, scd_dif, pch = 19, cex = 0.5, col = alpha("grey50", alpha = 0.5))


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
