library(dplyr)
library(readr)
library(lubridate)
library(zoo)
library(DEoptim)  # for optimization

# Source the HBV model function
# Make sure HBV.R is in your working directory
source('HBV.R')

start <- as.Date("2021-06-01")
end   <- as.Date("2022-01-31")

area_km2 <- 0.82

Lat_deg <- 37.2296
Lat_rad <- Lat_deg * pi / 180

q_file   <- "q_all.csv"
prcp_file <- "4161697.csv"

# Weather data download
asos_url <- function(station, start, end, tz = "America/New_York") {
  base <- "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?"
  sprintf(
    paste0(
      base,
      "station=%s",
      "&data=tmpf",
      "&year1=%d&month1=%d&day1=%d",
      "&year2=%d&month2=%d&day2=%d",
      "&tz=%s",
      "&format=onlycomma",
      "&latlon=no",
      "&missing=null"
    ),
    station,
    year(start), month(start), day(start),
    year(end),   month(end),   day(end),
    URLencode(tz)
  )
}

urls <- c(
  asos_url("BCB",  start, end),
  asos_url("KBCB", start, end)
)

# Temperature data
wx_raw <- NULL
for (u in urls) {
  message("Trying URL: ", u)
  tmp <- try(readr::read_csv(u, show_col_types = FALSE, na = c("M", "null")), silent = TRUE)
  if (!inherits(tmp, "try-error") && nrow(tmp) > 0) {
    wx_raw <- tmp
    break
  }
}
if (is.null(wx_raw) || nrow(wx_raw) == 0) stop("No data returned from IEM for BCB/KBCB.")

# Process temperature data
wx_daily <- wx_raw |>
  rename(timestamp = valid) |>
  mutate(
    timestamp = ymd_hms(timestamp, tz = "America/New_York"),
    date      = as.Date(timestamp)
  ) |>
  filter(!is.na(tmpf), date >= start, date <= end) |>
  group_by(date) |>
  summarise(
    tmean_f = mean(tmpf, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    Tmean = (tmean_f - 32) * 5/9
  ) |>
  select(date, Tmean) |>
  arrange(date)

temp_out <- "blacksburg_va_daily_Tmean_C.csv"
write_csv(wx_daily, temp_out)
message("Wrote: ", temp_out)

# Read discharge and precipitation
q_all <- read_csv(q_file, show_col_types = FALSE)
prcp  <- read_csv(prcp_file, show_col_types = FALSE)

q_daily <- q_all |>
  mutate(date = as.Date(dt)) |>
  group_by(date) |>
  summarise(Q_ls = mean(q.ls, na.rm = TRUE), .groups = "drop") |>
  mutate(Q_mm = Q_ls * 86.4 / area_km2) |>
  select(date, Q_mm)

prcp_daily <- prcp |>
  mutate(date = as.Date(DATE)) |>
  group_by(date) |>
  summarise(PRCP = sum(PRCP, na.rm = TRUE), .groups = "drop") |>
  mutate(P_mm = PRCP * 25.4) |>
  select(date, P_mm)

# Combine all data
df_all <- q_daily |>
  inner_join(prcp_daily, by = "date") |>
  inner_join(wx_daily,   by = "date") |>
  arrange(date)

dup_check <- df_all |>
  count(date) |>
  filter(n > 1)
if (nrow(dup_check) > 0) stop("Still have duplicate dates after daily aggregation. Check inputs.")

full_dates <- tibble(date = seq(min(df_all$date), max(df_all$date), by = "day"))

# Fill missing values
df_all_full <- full_dates |>
  left_join(df_all, by = "date") |>
  arrange(date) |>
  mutate(
    P_mm_filled = if_else(is.na(P_mm), 0, P_mm),
    Tmean_filled = zoo::na.approx(Tmean, x = date, na.rm = FALSE)
  )

df_all_full <- df_all_full |>
  mutate(
    Tmean_filled = zoo::na.locf(Tmean_filled, na.rm = FALSE),
    Tmean_filled = zoo::na.locf(Tmean_filled, fromLast = TRUE, na.rm = FALSE)
  )

if (any(is.na(df_all_full$Tmean_filled))) stop("Still have NA temperature after fill. Check temp series range.")

# Calculate PET using Hamon method
lat <- Lat_deg
latrad <- (lat/360) * 2 * pi

df_all_full <- df_all_full |>
  mutate(
    DOY = yday(date),
    tempvar = (2 * pi / 365) * DOY,
    delta_h = 0.4093 * sin(tempvar - 1.405),
    daylen = (2 * acos(-tan(delta_h) * tan(latrad)) / 0.2618),
    PET_mm = 29.8 * daylen * 0.611 * exp(17.3 * Tmean_filled / 
               (Tmean_filled + 237.3)) / (Tmean_filled + 273.2)
  ) |>
  select(date, P_mm_filled, Tmean_filled, PET_mm, Q_mm)

stopifnot(all(colSums(is.na(df_all_full[, c("P_mm_filled","Tmean_filled","PET_mm")])) == 0))

# Split data into warmup, calibration, and validation periods
warm_len <- 60
warmup_start <- min(df_all_full$date)
warmup_end   <- warmup_start + warm_len

calib_start  <- warmup_end + 1
calib_end    <- max(df_all_full$date)

valid_start  <- warmup_end + 1
valid_end    <- max(df_all_full$date)

# Define NSE function
calc_NSE <- function(obs, sim) {
  # Remove NAs
  valid_idx <- !is.na(obs) & !is.na(sim)
  obs <- obs[valid_idx]
  sim <- sim[valid_idx]
  
  if (length(obs) == 0) return(NA)
  
  numerator <- sum((sim - obs)^2)
  denominator <- sum((obs - mean(obs))^2)
  
  NSE <- 1 - (numerator / denominator)
  return(NSE)
}

# HBV model wrapper for calibration
run_hbv_model <- function(params, P, Temp, PET, Q_obs, 
                          calib_indices, routing = 0) {
  
  # Run HBV model
  model_output <- HBV(params, P, Temp, PET, routing)
  
  # Extract simulated discharge for calibration period
  Q_sim <- model_output$q[calib_indices]
  Q_obs_cal <- Q_obs[calib_indices]
  
  # Calculate NSE
  nse <- calc_NSE(Q_obs_cal, Q_sim)
  
  # Return negative NSE (for minimization)
  return(-nse)
}

# Prepare data vectors
P <- df_all_full$P_mm_filled
Temp <- df_all_full$Tmean_filled
PET <- df_all_full$PET_mm
Q <- df_all_full$Q_mm

# Define calibration and validation indices
calib_indices <- which(df_all_full$date >= calib_start & df_all_full$date <= calib_end)
valid_indices <- which(df_all_full$date >= valid_start & df_all_full$date <= valid_end)
all_indices <- which(df_all_full$date >= calib_start & df_all_full$date <= valid_end)

# Fill any NA values in Q for calibration period
Q_filled <- Q
if (any(is.na(Q_filled[calib_indices]))) {
  Q_filled <- zoo::na.approx(Q_filled, na.rm = FALSE)
  Q_filled <- zoo::na.locf(Q_filled, na.rm = FALSE)
  Q_filled <- zoo::na.locf(Q_filled, fromLast = TRUE, na.rm = FALSE)
}

# HBV parameter bounds
# Parameter order: FC, beta, LP, SFCF, TT, CFMAX, k0, k1, k2, UZL, PERC, MAXBAS
param_lower <- c(40,   1,   0.3,  0.4,  -1.5, 1,    0.05,  0.01,  0.001, 0,   0,    1)
param_upper <- c(400,  6,   1,    1.2,  1.2,  8,    0.5,   0.3,   0.15,  70,  4,    3)

message("\n--- Calibrating HBV Model ---")
message("This may take several minutes...")

# Run optimization using DEoptim
set.seed(123)  # for reproducibility
optim_result <- DEoptim(
  fn = run_hbv_model,
  lower = param_lower,
  upper = param_upper,
  control = DEoptim.control(
    NP = 50,          # population size
    itermax = 1000,    # max iterations
    trace = 10,       # print progress every 10 iterations
    parallelType = 0  # no parallelization
  ),
  P = P,
  Temp = Temp,
  PET = PET,
  Q_obs = Q_filled,
  calib_indices = calib_indices,
  routing = 0
)

# Extract best parameters
best_params <- optim_result$optim$bestmem
names(best_params) <- c("FC", "beta", "LP", "SFCF", "TT", "CFMAX", 
                        "k0", "k1", "k2", "UZL", "PERC", "MAXBAS")

message("\nOptimal HBV Parameters:")
print(round(best_params, 3))

# Run model with best parameters for full period
model_output <- HBV(best_params, P, Temp, PET, routing = 0)

# Add simulated discharge to dataframe
df_all_full$Qsim_HBV <- model_output$q

# Calculate performance metrics
nse_calib <- calc_NSE(Q_filled[calib_indices], df_all_full$Qsim_HBV[calib_indices])
nse_valid <- calc_NSE(Q_filled[valid_indices], df_all_full$Qsim_HBV[valid_indices])

# Create results summary
results <- data.frame(
  model = "HBV",
  nse_cal = nse_calib,
  nse_val = nse_valid,
  stringsAsFactors = FALSE
)

message("\n--- HBV Model Performance ---")
print(results)

# Save model output
write_csv(df_all_full, "hbv_model_output.csv")
message("\nWrote model output to: hbv_model_output.csv")

# Save plots to PDF files to avoid margin errors
message("\nCreating plots...")

# Plot 1: Full time series
pdf("hbv_timeseries.pdf", width = 10, height = 6)
par(mar = c(5, 4, 4, 2) + 0.1)
plot(df_all_full$date, df_all_full$Q_mm,
     type = "l", col = "black", lwd = 2,
     xlab = "Date", ylab = "Discharge (mm/day)",
     main = "Observed vs Simulated Discharge (HBV Model)")

lines(df_all_full$date, df_all_full$Qsim_HBV, col = "blue")

# Add vertical lines for period boundaries
abline(v = calib_start, col = "red", lty = 2)
abline(v = calib_end, col = "red", lty = 2)
abline(v = valid_start, col = "darkgreen", lty = 2)

legend("topright",
       legend = c("Observed", "HBV Simulated", "Calib. Period", "Valid. Period"),
       col = c("black", "blue", "red", "darkgreen"),
       lty = c(1, 1, 2, 2),
       lwd = c(2, 1, 1, 1),
       cex = 0.8)
dev.off()
message("Saved: hbv_timeseries.pdf")

# Plot 2: Scatter plot
pdf("hbv_scatter.pdf", width = 7, height = 7)
par(mar = c(5, 4, 4, 2) + 0.1)
plot(df_all_full$Q_mm, df_all_full$Qsim_HBV,
     xlab = "Observed Q (mm/day)", ylab = "Simulated Q (mm/day)",
     main = sprintf("HBV Model: NSE_cal = %.3f, NSE_val = %.3f", nse_calib, nse_valid),
     pch = 16, col = rgb(0, 0, 1, 0.4))
abline(0, 1, col = "black", lwd = 2)
dev.off()
message("Saved: hbv_scatter.pdf")

# Plot 3: Calibration and Validation periods
pdf("hbv_cal_val_periods.pdf", width = 10, height = 8)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

calib_data <- df_all_full[calib_indices, ]
plot(calib_data$date, calib_data$Q_mm,
     type = "l", col = "black", lwd = 2,
     xlab = "Date", ylab = "Discharge (mm/day)",
     main = sprintf("Calibration Period (NSE = %.3f)", nse_calib))
lines(calib_data$date, calib_data$Qsim_HBV, col = "blue")
legend("topright", legend = c("Observed", "Simulated"), 
       col = c("black", "blue"), lty = 1, lwd = c(2, 1), cex = 0.8)

valid_data <- df_all_full[valid_indices, ]
plot(valid_data$date, valid_data$Q_mm,
     type = "l", col = "black", lwd = 2,
     xlab = "Date", ylab = "Discharge (mm/day)",
     main = sprintf("Validation Period (NSE = %.3f)", nse_valid))
lines(valid_data$date, valid_data$Qsim_HBV, col = "blue")
legend("topright", legend = c("Observed", "Simulated"), 
       col = c("black", "blue"), lty = 1, lwd = c(2, 1), cex = 0.8)

dev.off()
message("Saved: hbv_cal_val_periods.pdf")

message("\n--- Analysis Complete ---")
