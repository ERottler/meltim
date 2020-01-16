###

#Rhine snow - Runoff timing: classical approach
#Erwin Rottler, University of Potsdam
#Spring 2019

###

gauge_sel <- "Domat" 
# Domat, Martinsbruck, Porte_du_Scex, Wasserburg, Burghausen, 
# Diepoldsau, Neuhausen, Andelfingen, Appenzell, Rekingen, Mellingen, Bern_Schoenau, 
# Emmenmatt, Montier, Murgenthal, Brugg, Untersiggenthal, Basel, Maxau, Koeln, Wuerzburg, Schwaibach, Thoerishaus,
# Landsberg, Kempten, Kochel, Chancy_aux_Ripes

if(gauge_sel == "Domat"){           file_sel <- "6935145_Q_Day.Cmd.txt"}
if(gauge_sel == "Martinsbruck"){    file_sel <- "6943100_Q_Day.Cmd.txt"}
if(gauge_sel == "Porte_du_Scex"){   file_sel <- "6939200_Q_Day.Cmd.txt"}
if(gauge_sel == "Wasserburg"){      file_sel <- "6343100_Q_Day.Cmd.txt"}
if(gauge_sel == "Burghausen"){      file_sel <- "6343500_Q_Day.Cmd.txt"}
if(gauge_sel == "Diepoldsau"){      file_sel <- "6935500_Q_Day.Cmd.txt"}
if(gauge_sel == "Mellingen"){       file_sel <- "6935310_Q_Day.Cmd.txt"}
if(gauge_sel == "Untersiggenthal"){ file_sel <- "6935300_Q_Day.Cmd.txt"}
if(gauge_sel == "Neuhausen"){       file_sel <- "6935055_Q_Day.Cmd.txt"}
if(gauge_sel == "Andelfingen"){     file_sel <- "6935400_Q_Day.Cmd.txt"}
if(gauge_sel == "Rekingen"){        file_sel <- "6935054_Q_Day.Cmd.txt"}
if(gauge_sel == "Bern_Schoenau"){   file_sel <- "6935020_Q_Day.Cmd.txt"}
if(gauge_sel == "Montier"){         file_sel <- "6935060_Q_Day.Cmd.txt"}
if(gauge_sel == "Appenzell"){       file_sel <- "6935412_Q_Day.Cmd.txt"}
if(gauge_sel == "Murgenthal"){      file_sel <- "6935302_Q_Day.Cmd.txt"}
if(gauge_sel == "Brugg"){           file_sel <- "6935301_Q_Day.Cmd.txt"}
if(gauge_sel == "Basel"){           file_sel <- "6935051_Q_Day.Cmd.txt"}
if(gauge_sel == "Schwaibach"){      file_sel <- "6335125_Q_Day.Cmd.txt"}
if(gauge_sel == "Maxau"){           file_sel <- "6335200_Q_Day.Cmd.txt"}
if(gauge_sel == "Koeln"){           file_sel <- "6335060_Q_Day.Cmd.txt"}
if(gauge_sel == "Wuerzburg"){       file_sel <- "6335500_Q_Day.Cmd.txt"}
if(gauge_sel == "Landsberg"){       file_sel <- "6342513_Q_Day.Cmd.txt"}
if(gauge_sel == "Kempten"){         file_sel <- "6342200_Q_Day.Cmd.txt"}
if(gauge_sel == "Kochel"){          file_sel <- "6342930_Q_Day.Cmd.txt"}
if(gauge_sel == "Chancy_aux_Ripes"){file_sel <- "6939050_Q_Day.Cmd.txt"}

grdc_data <- read_grdc(paste0(grdc_dir, file_sel)) ; summary(grdc_data$value)


# #meta_map----
# 
# file_names <- list.files(path = grdc_dir, pattern = "*.Cmd", full.names = F)
# file_paths <- list.files(path = grdc_dir, pattern = "*.Cmd", full.names = T)
# 
# #read_meta----
# 
# for(i in 1:length(file_paths)){
#   
#   print(i)
#   
#   #get rows with meta information
#   meta_rows <- read_lines(file_paths[i], n_max = 32)
#   meta_rows <- iconv(meta_rows, "UTF-8", "ASCII", "")
#   #Name
#   row_name <- meta_rows[11]
#   sta_name <- substr(row_name, 26, nchar(row_name))
#   #Longitude
#   row_long <- meta_rows[14]
#   sta_long <- substr(row_long, 24, nchar(row_long))
#   #Latitude
#   row_lati <- meta_rows[13]
#   sta_lati <- substr(row_lati, 24, nchar(row_lati))
#   #Altitude
#   row_alti <- meta_rows[16]
#   sta_alti <- substr(row_alti, 28, nchar(row_alti))
#   #Start/End time series 
#   row_seri <- meta_rows[24]
#   sta_seri <- substr(row_seri, 26, nchar(row_seri)-13)
#   end_seri <- substr(row_seri, 36, nchar(row_seri)-3)
#   #Catchment area
#   row_catc <- meta_rows[15]
#   sta_catc <- substr(row_catc, 30, nchar(row_catc))
#   #Number of years
#   row_year <- meta_rows[25]
#   sta_year <- substr(row_year, 26, nchar(row_year))
#   
#   meta_sing <- c(sta_name, sta_lati, sta_long, sta_alti, sta_catc, sta_seri, end_seri, sta_year, file_names[i])
#   
#   
#   if(i == 1){
#     
#     grdc_meta <- meta_sing
#     
#   }else{
#     
#     grdc_meta <- rbind(grdc_meta, meta_sing)
#     
#   }
#   
#   
# }
# 
# colnames(grdc_meta) <- c("name", "latitude", "longitude", "altitude", "cath_area", "start_series", "end_series", "n_years", "file")
# rownames(grdc_meta) <- NULL
# grdc_meta <- as.data.frame(grdc_meta)
# grdc_meta$latitude   <- as.numeric(levels(grdc_meta$latitude))[grdc_meta$latitude]
# grdc_meta$longitude  <- as.numeric(levels(grdc_meta$longitude))[grdc_meta$longitude]
# grdc_meta$altitude   <- as.numeric(levels(grdc_meta$altitude))[grdc_meta$altitude]
# grdc_meta$cath_area  <- as.numeric(levels(grdc_meta$cath_area))[grdc_meta$cath_area]
# grdc_meta$start_series  <- as.numeric(levels(grdc_meta$start_series))[grdc_meta$start_series]
# grdc_meta$end_series  <- as.numeric(levels(grdc_meta$end_series))[grdc_meta$end_series]
# grdc_meta$n_years  <- as.numeric(levels(grdc_meta$n_years))[grdc_meta$n_years]
# 
# sel_ind <- which(grdc_meta[, 1] %in% c("DOMAT/EMS", "MARTINSBRUCK", "PORTE DU SCEX", "DIEPOLDSAU, RIETBRUECKE", "WASSERBURG",
#                                        "BURGHAUSEN", "NEUHAUSEN, FLURLINGERBRUECKE", "REKINGEN", "BASEL, RHEINHALLE", 
#                                        "BERN-SCHOENAU", "MURGENTHAL", "BRUGG", "MELLINGEN", "UNTERSIGGENTHAL, STILLI",
#                                        "MAXAU"))
# 
# meta_sel <- grdc_meta[sel_ind, ]
