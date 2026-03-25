library(dplyr)
library(readr)
library(ggplot2)
library(lubridate)
library(stringr)

# -------------------------------
# Read discharge data
# -------------------------------
q_file <- "q_all.csv"
q_all  <- read_csv(q_file, show_col_types = FALSE)

# -------------------------------
# Process to daily discharge
# -------------------------------
q_daily <- q_all |>
  mutate(
    date = as.Date(dt),
    stream = str_trim(toupper(stream)),
    stream = case_when(
      stream == "BEAST" ~ "Scarp",
      stream == "MDH"   ~ "Dip",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(stream)) |>
  group_by(date, stream) |>
  summarise(Q_ls = mean(q.ls, na.rm = TRUE), .groups = "drop")

# -------------------------------
# Keep only dates shared by both streams
# -------------------------------
shared_dates <- q_daily |>
  count(date) |>
  filter(n == 2) |>
  pull(date)

q_plot <- q_daily |>
  filter(date %in% shared_dates) |>
  mutate(stream = factor(stream, levels = c("Scarp", "Dip")))

# -------------------------------
# Plot
# -------------------------------
ggplot(q_plot, aes(x = date, y = Q_ls, color = stream)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("Scarp" = "#6A3D9A", "Dip" = "#E66101"),
    drop = FALSE
  ) +
  labs(
    title = "Observed Discharge During Shared Monitoring Period",
    x = "Date",
    y = "Discharge (L/s)",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("shared_discharge_timeseries.png", width = 10, height = 5, dpi = 300)