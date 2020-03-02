###

#Calibration Snow model - Objective function
#Erwin Rottler, Summer 2019

###

obj_func <- function(modelled, measured, obj_method = "RMSE"){
  
  if(obj_method == "Correlation"){
    
    cor_calc <- cor(modelled, measured, method = "pearson", use = "pairwise.complete.obs")
    
    #ppso algorithm is minimizing result of objective runction
    #we want as as close to 1 as possible
    val_out <- 1-cor_calc
    
  }
  
  if(obj_method == "RMSE"){
    
    rmse_calc <- nrmse(modelled, measured, na.rm = T)
    
    #ppso algorithm is minimizing result of objective runction
    #smaller rmse value indicates better model performance
    val_out <- rmse_calc
    
  }
  # if(obj_method == "KGE"){
  #   
  #   kge_calc <- NSE(modelled, measured, na.rm = T, method = "2012")
  #   
  #   #ppso algorithm is minimizing result of objective runction
  #   #KGE takes vlues between -Inf and 1
  #   #we want as as close to 1 as possible
  #   val_out <- -(kge_calc-1)
  #   
  # }
  
  return(val_out)
  
}


