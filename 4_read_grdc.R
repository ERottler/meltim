###

#Rhine snow - Runoff timing: classical approach
#Erwin Rottler, University of Potsdam
#Spring 2019

###

gauge_sel <- "Diepoldsau" 
# Domat, Martinsbruck, Porte_du_Scex, Wasserburg, Burghausen, 
# Diepoldsau, Neuhausen, Andelfingen, Appenzell, Rekingen, Mellingen, Bern_Schoenau, 
# Emmenmatt, Montier, Murgenthal, Brugg, Untersiggenthal, Basel, Maxau, Koeln, Wuerzburg, Schwaibach, Thoerishaus

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