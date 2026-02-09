# HBV Hydrological Model
# Based on HBV-light structure as described in Seibert & Vis (2012)
# Adapted for Virginia Tech Hydroinformatics course

HBV <- function(pars, P, Temp, PET, routing = 0) {
  
  # Extract parameters from pars vector
  FC     <- pars[1]   # Max soil moisture storage (field capacity)
  beta   <- pars[2]   # Shape coefficient for soil moisture
  LP     <- pars[3]   # Threshold for reduction of evap
  SFCF   <- pars[4]   # Snowfall correction factor
  TT     <- pars[5]   # Threshold temperature
  CFMAX  <- pars[6]   # Degree-day factor
  k0     <- pars[7]   # Recession constant (upper storage, near surface)
  k1     <- pars[8]   # Recession constant (upper storage)
  k2     <- pars[9]   # Recession constant (lower storage)
  UZL    <- pars[10]  # Threshold for shallow storage
  PERC   <- pars[11]  # Percolation, max flow from upper to lower storage
  MAXBAS <- pars[12]  # Base of triangular routing function (days)
  
  n <- length(P)
  
  # Initialize state variables
  SWE <- 0      # Snow water equivalent
  SM  <- 0      # Soil moisture
  SUZ <- 0      # Upper zone storage
  SLZ <- 0      # Lower zone storage
  
  # Initialize output vectors
  q       <- numeric(n)  # Total discharge
  qs      <- numeric(n)  # Surface/fast runoff
  qi      <- numeric(n)  # Interflow
  qb      <- numeric(n)  # Baseflow
  Storage <- numeric(n)  # Total storage
  SWE_out <- numeric(n)  # Snow water equivalent
  AET     <- numeric(n)  # Actual evapotranspiration
  SF      <- numeric(n)  # Snowfall
  S1      <- numeric(n)  # Upper zone storage
  S2      <- numeric(n)  # Lower zone storage
  soil    <- numeric(n)  # Soil moisture
  w       <- numeric(n)  # Water input to soil
  
  for (i in 1:n) {
    
    # Snow routine
    if (Temp[i] < TT) {
      # It's snowing
      snowfall <- P[i] * SFCF
      rainfall <- 0
      SWE <- SWE + snowfall
      SF[i] <- snowfall
    } else {
      # It's raining
      snowfall <- 0
      rainfall <- P[i]
      # Snowmelt
      melt <- min(CFMAX * (Temp[i] - TT), SWE)
      SWE <- SWE - melt
      rainfall <- rainfall + melt
      SF[i] <- 0
    }
    
    # Total water input to soil
    water_input <- rainfall
    w[i] <- water_input
    
    # Soil moisture routine
    # Evapotranspiration
    if (SM < (LP * FC)) {
      AET[i] <- PET[i] * (SM / (LP * FC))
    } else {
      AET[i] <- PET[i]
    }
    AET[i] <- min(AET[i], SM)
    SM <- SM - AET[i]
    
    # Recharge to upper zone
    if (SM + water_input > FC) {
      recharge <- water_input - (FC - SM)
      SM <- FC
    } else {
      recharge <- water_input * ((SM + water_input) / FC)^beta
      SM <- SM + water_input - recharge
    }
    
    # Upper zone storage
    SUZ <- SUZ + recharge
    
    # Percolation to lower zone
    perc_actual <- min(PERC, SUZ)
    SUZ <- SUZ - perc_actual
    SLZ <- SLZ + perc_actual
    
    # Runoff from upper zone
    if (SUZ > UZL) {
      q0 <- k0 * (SUZ - UZL)
      SUZ <- SUZ - q0
    } else {
      q0 <- 0
    }
    
    q1 <- k1 * SUZ
    SUZ <- SUZ - q1
    
    # Runoff from lower zone
    q2 <- k2 * SLZ
    SLZ <- SLZ - q2
    
    # Total discharge
    q_total <- q0 + q1 + q2
    
    # Store values
    qs[i] <- q0
    qi[i] <- q1
    qb[i] <- q2
    q[i] <- q_total
    
    SWE_out[i] <- SWE
    S1[i] <- SUZ
    S2[i] <- SLZ
    soil[i] <- SM
    Storage[i] <- SM + SUZ + SLZ + SWE
  }
  
  # Apply routing if requested
  if (routing == 1 && MAXBAS > 1) {
    # Triangular unit hydrograph routing
    q_routed <- numeric(n)
    uh_length <- ceiling(MAXBAS)
    
    # Create triangular unit hydrograph
    uh <- numeric(uh_length)
    peak_time <- MAXBAS / 2
    for (j in 1:uh_length) {
      if (j <= peak_time) {
        uh[j] <- j / peak_time
      } else {
        uh[j] <- 2 - (j / peak_time)
      }
    }
    uh <- uh / sum(uh)  # Normalize
    
    # Convolve discharge with unit hydrograph
    for (i in 1:n) {
      for (j in 1:min(uh_length, i)) {
        q_routed[i] <- q_routed[i] + q[i - j + 1] * uh[j]
      }
    }
    q <- q_routed
  }
  
  # Return output as data frame
  output <- data.frame(
    q = q,
    qs = qs,
    qi = qi,
    qb = qb,
    Storage = Storage,
    SWE = SWE_out,
    AET = AET,
    SF = SF,
    S1 = S1,
    S2 = S2,
    soil = soil,
    w = w
  )
  
  return(output)
}
