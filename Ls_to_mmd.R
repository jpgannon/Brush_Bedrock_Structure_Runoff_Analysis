Ls_to_mmd <- function(dischargeLS, WSsizeHa){
  #Ha to mm2
  #1 Ha = 10000000000 mm2 (1e10)
  WSsizemm2 <- WSsizeHa * 1e10
  
  #L/s to mm3/s
  #1L = 1000000 mm3 (1e6)
  dischargemm3s <- dischargeLS * 1e6
  
  #Discharge in mm3/s / WS Size in mm2
  # = Area normalized discharge in mm/s
  Discharge_mms <- dischargemm3s / WSsizemm2
  
  #sec to day
  #1 day = 86400 sec
  Discharge_mmd <- Discharge_mms * 86400
  
  #output
  Discharge_mmd
}

MDH <- 62.84
Beast <- 82.85

#single values
Ls_to_mmd(dischargeLS = 34,WSsizeHa = Beast)
Ls_to_mmd(dischargeLS = 30,WSsizeHa = MDH)

#vector of data
Qdat <- c(10,12,15,20) #totally made up
Ls_to_mmd(dischargeLS = Qdat, WSsizeHa = Beast)
