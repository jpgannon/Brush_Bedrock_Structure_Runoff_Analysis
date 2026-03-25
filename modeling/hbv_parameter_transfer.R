# ============================================================
# HBV Calibration + Parameter Transfer
# McDonald (HydroVu) -> Beast
# ============================================================

# -------------------------------
# Libraries
# -------------------------------
library(dplyr)
library(readr)
library(lubridate)
library(zoo)
library(DEoptim)
library(daymetr)
library(ggplot2)
library(stringr)

source("HBV.R")

# -------------------------------
# USER SETTINGS
# -------------------------------
hydrovu_dir <- "data/hydrovu"
q_all_path  <- "data/q_all.csv"

area_km2_mcd   <- 0.6284
area_km2_beast <- 0.82

mcd_lat   <- 37.265
mcd_lon   <- -80.428
beast_lat <- 37.2296
beast_lon <- -80.4101

tz_local <- "America/New_York"
warm_len <- 60

# -------------------------------
# Helper functions
# -------------------------------
Ls_to_mmd <- function(Q_ls, area_km2) Q_ls * 86.4 / area_km2

calc_NSE <- function(obs, sim) {
  idx <- !is.na(obs) & !is.na(sim)
  obs <- obs[idx]; sim <- sim[idx]
  if (length(obs) < 5) return(NA_real_)
  1 - sum((sim - obs)^2) / sum((obs - mean(obs))^2)
}

calc_RMSE <- function(obs, sim) {
  idx <- !is.na(obs) & !is.na(sim)
  obs <- obs[idx]; sim <- sim[idx]
  if (length(obs) < 5) return(NA_real_)
  sqrt(mean((sim - obs)^2))
}

fill_inputs <- function(df) {
  df %>%
    arrange(date) %>%
    mutate(
      P_mm = if_else(is.na(P_mm), 0, P_mm),
      Tmean = zoo::na.approx(Tmean, x = date, na.rm = FALSE),
      PET_mm = zoo::na.approx(PET_mm, x = date, na.rm = FALSE)
    ) %>%
    mutate(
      Tmean = zoo::na.locf(Tmean, na.rm = FALSE),
      Tmean = zoo::na.locf(Tmean, fromLast = TRUE, na.rm = FALSE),
      PET_mm = zoo::na.locf(PET_mm, na.rm = FALSE),
      PET_mm = zoo::na.locf(PET_mm, fromLast = TRUE, na.rm = FALSE)
    )
}

calc_hamon_pet <- function(date, Tmean_C, lat_deg) {
  latrad <- (lat_deg / 360) * 2 * pi
  DOY <- yday(date)
  delta_h <- 0.4093 * sin((2 * pi / 365) * DOY - 1.405)
  daylen <- (2 * acos(-tan(delta_h) * tan(latrad)) / 0.2618)
  29.8 * daylen * 0.611 * exp(17.3 * Tmean_C / (Tmean_C + 237.3)) /
    (Tmean_C + 273.2)
}

get_daymet_forcings <- function(lat, lon, start_date, end_date) {
  yrs <- seq(year(start_date), year(end_date))
  dm <- download_daymet(
    site = "HBV_site",
    lat = lat, lon = lon,
    start = min(yrs), end = max(yrs),
    internal = TRUE,
    silent = TRUE,
    simplify = FALSE
  )
  dm_tbl <- if (is.list(dm) && !is.null(dm$data)) dm$data else dm
  dm_tbl %>%
    mutate(
      date = as.Date(paste(year, yday, sep = "-"), format = "%Y-%j"),
      Tmean = (tmin..deg.c. + tmax..deg.c.) / 2
    ) %>%
    filter(date >= start_date, date <= end_date) %>%
    transmute(date, P_mm = prcp..mm.day., Tmean)
}

# -------------------------------
# Robust HydroVu reader (NO parse_date_time)
# -------------------------------
read_hydrovu_csv <- function(path) {
  
  raw <- readLines(path, warn = FALSE)
  
  header_i <- which(grepl("date time|datetime|timestamp", tolower(raw)))
  if (length(header_i) == 0) {
    stop("No datetime header found in: ", basename(path))
  }
  
  df <- read_csv(path, skip = header_i[1] - 1, show_col_types = FALSE)
  
  time_col <- names(df)[grepl("date|time", tolower(names(df)))][1]
  if (is.na(time_col)) {
    stop("No datetime column found in: ", basename(path))
  }
  
  df <- df %>%
    rename(datetime_raw = !!sym(time_col)) %>%
    mutate(
      datetime_chr = as.character(datetime_raw),
      datetime_chr = sub("Z$", "", datetime_chr),
      datetime_chr = sub("\\+\\d\\d:?\\d\\d$", "", datetime_chr),
      datetime_chr = sub("-\\d\\d:?\\d\\d$", "", datetime_chr),
      datetime_utc = suppressWarnings(ymd_hms(datetime_chr, tz = "UTC", quiet = TRUE)),
      datetime_loc = with_tz(datetime_utc, tz_local)
    ) %>%
    filter(!is.na(datetime_loc))
  
  if (nrow(df) == 0) {
    stop("Datetime parsing failed for: ", basename(path),
         "\nCheck the datetime format in the CSV.")
  }
  
  df
}

# -------------------------------
# Build McDonald stage (STRICT + unit-safe)
# -------------------------------
build_mcd_stage_daily <- function(dir) {
  
  files <- list.files(dir, pattern="\\.csv$", full.names=TRUE)
  if (length(files) == 0) stop("No CSVs found in: ", dir)
  
  dfs <- lapply(files, read_hydrovu_csv)
  
  # Combine all HydroVu data into one long table
  all <- bind_rows(lapply(seq_along(dfs), function(i) {
    df <- dfs[[i]]
    df$source_file <- basename(files[i])
    df
  }))
  
  nms <- names(all)
  nms_low <- tolower(nms)
  
  # ----- Candidate columns -----
  # Depth/Stage/Level columns (only accept if units suggest meters or cm)
  depth_m_col <- nms[
    str_detect(nms_low, "depth|water level|stage|\\blevel\\b") &
      (str_detect(nms_low, "\\(m\\)") | str_detect(nms_low, " m") | str_detect(nms_low, "_m\\b"))
  ][1]
  
  depth_cm_col <- nms[
    str_detect(nms_low, "depth|water level|stage|\\blevel\\b") &
      str_detect(nms_low, "cm")
  ][1]
  
  # Baro psi (must contain baro and psi)
  baro_col <- nms[str_detect(nms_low, "baro") & str_detect(nms_low, "psi")][1]
  
  # Pressure psi (must contain psi or pressure, but NOT baro)
  press_candidates <- nms[str_detect(nms_low, "press|psi")]
  press_candidates <- press_candidates[!str_detect(tolower(press_candidates), "baro")]
  press_psi_col <- press_candidates[str_detect(tolower(press_candidates), "psi")][1]
  
  # ----- Build stage time series -----
  if (!is.na(depth_m_col)) {
    stage_ts <- tibble(
      datetime_loc = all$datetime_loc,
      stage_m = suppressWarnings(as.numeric(all[[depth_m_col]]))
    )
    
  } else if (!is.na(depth_cm_col)) {
    stage_ts <- tibble(
      datetime_loc = all$datetime_loc,
      stage_m = suppressWarnings(as.numeric(all[[depth_cm_col]])) / 100
    )
    
  } else {
    # Pressure-based approach
    if (is.na(press_psi_col)) stop("No Depth/Stage column with units found AND no Pressure (psi) column found.")
    if (is.na(baro_col)) stop("Pressure (psi) found but no Baro (psi) column found. Export/upload the baro CSV.")
    
    press_ts <- tibble(
      datetime_loc = all$datetime_loc,
      press_psi = suppressWarnings(as.numeric(all[[press_psi_col]]))
    ) %>%
      filter(!is.na(datetime_loc), is.finite(press_psi)) %>%
      arrange(datetime_loc)
    
    baro_ts <- tibble(
      datetime_loc = all$datetime_loc,
      baro_psi = suppressWarnings(as.numeric(all[[baro_col]]))
    ) %>%
      filter(!is.na(datetime_loc), is.finite(baro_psi)) %>%
      arrange(datetime_loc)
    
    if (nrow(press_ts) < 2) stop("Not enough pressure points (<2).")
    if (nrow(baro_ts) < 2)  stop("Not enough baro points (<2).")
    
    baro_interp <- approx(
      x = as.numeric(baro_ts$datetime_loc),
      y = baro_ts$baro_psi,
      xout = as.numeric(press_ts$datetime_loc),
      rule = 2
    )$y
    
    stage_ts <- press_ts %>%
      mutate(
        gauge_psi = press_psi - baro_interp,
        stage_m = (gauge_psi * 6894.757) / (1000 * 9.80665)
      ) %>%
      select(datetime_loc, stage_m)
  }
  
  # Daily aggregation
  mcd_stage_daily <- stage_ts %>%
    filter(!is.na(datetime_loc), is.finite(stage_m)) %>%
    mutate(date = as.Date(datetime_loc, tz = tz_local)) %>%
    group_by(date) %>%
    summarise(Stage_m = mean(stage_m, na.rm = TRUE), .groups = "drop") %>%
    arrange(date)
  
  # HARD sanity check: stage should be in meters for a small stream
  s <- summary(mcd_stage_daily$Stage_m)
  if (is.finite(s[["Median"]]) && s[["Median"]] > 3) {
    stop(
      "Stage magnitude looks wrong (median > 3 m). ",
      "Your HydroVu extraction is still picking a non-meter column.\n",
      "Median Stage_m = ", round(s[["Median"]], 3),
      ".\nFix by hardcoding the correct HydroVu column name."
    )
  }
  
  mcd_stage_daily
}

# -------------------------------
# McDonald stage -> discharge
# -------------------------------
mcd_stage_daily <- build_mcd_stage_daily(hydrovu_dir)

mcd_q <- mcd_stage_daily %>%
  mutate(
    Q_ls = 1.737801e-21 * (Stage_m * 100)^14.4394,
    Q_mm = Ls_to_mmd(Q_ls, area_km2_mcd)
  ) %>%
  select(date, Q_mm)

# -------------------------------
# Beast discharge
# -------------------------------
q_all <- read_csv(q_all_path, show_col_types = FALSE)

beast_q <- q_all %>%
  mutate(date = as.Date(dt)) %>%
  group_by(date) %>%
  summarise(Q_mm = Ls_to_mmd(mean(q.ls, na.rm = TRUE), area_km2_beast), .groups = "drop")

# -------------------------------
# Date windows (separate per basin)
# -------------------------------
mcd_start <- min(mcd_q$date)
mcd_end   <- max(mcd_q$date)

beast_start <- min(beast_q$date)
beast_end   <- max(beast_q$date)

# -------------------------------
# Forcings
# -------------------------------
forc_mcd <- get_daymet_forcings(mcd_lat, mcd_lon, mcd_start, mcd_end) %>%
  mutate(PET_mm = calc_hamon_pet(date, Tmean, mcd_lat))

forc_beast <- get_daymet_forcings(beast_lat, beast_lon, beast_start, beast_end) %>%
  mutate(PET_mm = calc_hamon_pet(date, Tmean, beast_lat))

df_mcd   <- forc_mcd %>% left_join(mcd_q, by = "date") %>% fill_inputs()
df_beast <- forc_beast %>% left_join(beast_q, by = "date") %>% fill_inputs()

# -------------------------------
# Calibration indices
# -------------------------------
warm_end <- min(df_mcd$date) + warm_len
cal_idx <- which(df_mcd$date > warm_end & !is.na(df_mcd$Q_mm))

# -------------------------------
# HBV calibration
# -------------------------------
run_hbv <- function(params, P, T, PET, Q_obs, idx) {
  sim <- HBV(params, P, T, PET, routing = 0)$q
  -calc_NSE(Q_obs[idx], sim[idx])
}

lower <- c(40,1,0.3,0.4,-1.5,1,0.05,0.01,0.001,0,0,1)
upper <- c(400,6,1.0,1.2,1.2,8,0.5,0.3,0.15,70,4,3)

opt <- DEoptim(
  run_hbv, lower, upper,
  P = df_mcd$P_mm,
  T = df_mcd$Tmean,
  PET = df_mcd$PET_mm,
  Q_obs = df_mcd$Q_mm,
  idx = cal_idx,
  control = DEoptim.control(itermax = 600, NP = 120, trace = 10)
)

best_params <- opt$optim$bestmem

# -------------------------------
# Run model
# -------------------------------
sim_mcd   <- HBV(best_params, df_mcd$P_mm, df_mcd$Tmean, df_mcd$PET_mm)$q
sim_beast <- HBV(best_params, df_beast$P_mm, df_beast$Tmean, df_beast$PET_mm)$q

# -------------------------------
# Plots
# -------------------------------
df_mcd %>% mutate(Q_sim = sim_mcd) %>%
  ggplot(aes(date)) +
  geom_line(aes(y = Q_mm, color = "Observed")) +
  geom_line(aes(y = Q_sim, color = "Simulated")) +
  theme_bw() +
  labs(title = "McDonald HBV Calibration", y = "Q (mm/day)", color = "")

df_beast %>% mutate(Q_sim = sim_beast) %>%
  ggplot(aes(date)) +
  geom_line(aes(y = Q_mm, color = "Observed")) +
  geom_line(aes(y = Q_sim, color = "Simulated")) +
  theme_bw() +
  labs(title = "Beast Parameter Transfer", y = "Q (mm/day)", color = "")

# -------------------------------
# Diagnostics
# -------------------------------
cat("\n--- Stage diagnostics (McDonald) ---\n")
print(summary(mcd_stage_daily$Stage_m))
print(range(mcd_stage_daily$Stage_m, na.rm = TRUE))

mcd_q_debug <- mcd_stage_daily %>%
  mutate(
    Q_ls = 1.737801e-21 * (Stage_m * 100)^14.4394,
    Q_mm = Ls_to_mmd(Q_ls, area_km2_mcd)
  )

cat("\n--- Q diagnostics (McDonald) ---\n")
print(summary(mcd_q_debug$Q_mm))
print(range(mcd_q_debug$Q_mm, na.rm = TRUE))

