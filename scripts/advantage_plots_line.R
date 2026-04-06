library(ggplot2)
library(dplyr)
library(readr)

library(rstudioapi)

data_path <- getActiveDocumentContext()$path
setwd(dirname(data_path))
path <- paste0(dirname(data_path),"/Data/Model_results/Crop_yield_final/csvs/")


scenario_map <- data.frame(
  file = c(
    "C4Rice_baseline_2099.csv",
    "C4Rice_SSP1-2.6_2050.csv", "C4Rice_SSP1-2.6_2099.csv",
    "C4Rice_SSP2-4.5_2050.csv", "C4Rice_SSP2-4.5_2099.csv",
    "C4Rice_SSP4-6.0_2050.csv", "C4Rice_SSP4-6.0_2099.csv",
    "C4Rice_SSP5-8.5_2050.csv", "C4Rice_SSP5-8.5_2099.csv"
  ),
  SSP = c(
    "Baseline",
    "SSP1-2.6", "SSP1-2.6",
    "SSP2-4.5", "SSP2-4.5",
    "SSP4-6.0", "SSP4-6.0",
    "SSP5-8.5", "SSP5-8.5"
  ),
  Year = c(
    "Baseline",
    "2050", "2099",
    "2050", "2099",
    "2050", "2099",
    "2050", "2099"
  )
)

get_net_advantage <- function(filename) {
  full_path <- paste0(path, filename)
  if (!file.exists(full_path)) {
    warning(paste("File missing:", full_path))
    return(NA_real_) # Ensure we return a double NA
  }

  df <- read_csv(full_path, show_col_types = FALSE)

  pos_perc <- mean(df$diff > 0, na.rm = TRUE) * 100
  neg_perc <- mean(df$diff < 0, na.rm = TRUE) * 100
  return(pos_perc - neg_perc)
}

baseline_val <- get_net_advantage("C4Rice_baseline_2099.csv")

scenarios <- c("SSP1-2.6", "SSP2-4.5", "SSP4-6.0", "SSP5-8.5")
years <- c("2050", "2099")

plot_data <- expand.grid(
  SSP = c("Baseline", scenarios), 
  Year = years, 
  stringsAsFactors = FALSE
)

plot_data <- plot_data %>%
  rowwise() %>%
  mutate(Net_Value = if_else(
    SSP == "Baseline",
    baseline_val,
    get_net_advantage(paste0("C4Rice_", SSP, "_", Year, ".csv"))
  )) %>%
  ungroup()

plot_data$SSP <- factor(plot_data$SSP, levels = c("Baseline", "SSP1-2.6", "SSP2-4.5", "SSP4-6.0", "SSP5-8.5"))
plot_data$Year <- factor(plot_data$Year, levels = c("2050", "2099"))

ggplot(plot_data, aes(x = SSP, y = Net_Value, color = Year, group = Year)) +
  geom_hline(yintercept = 0, color = "black", size = 0.8) +
  geom_line(size = 1.3) +
  geom_point(size = 4) +
  scale_y_continuous(
    limits = c(-105, 105),
    breaks = seq(-100, 100, by = 25),
    labels = function(x) paste0(ifelse(x > 0, "+", ""), x, "%")
  ) +
  scale_x_discrete(labels = c("Baseline" = "Historical baseline\n(1982-2013)")) +
  # Legend and Colors
  scale_color_manual(values = c("2050" = "#66c2a5", "2099" = "#fc8d62")) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Projection Scenario",
    y = "Net % Area Advantage (C4 vs. C3)",
    color = "Time Period"
  ) +
  theme_light(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    axis.text.x = element_text(face = "bold", size = 11),
    axis.title.y = element_text(face = "bold")
  ) +
  annotate("text", x = 1.5, y = 80, label = "C4 Advantage", color = "gray30", fontface = "italic") +
  annotate("text", x = 1.5, y = -80, label = "C3 Advantage", color = "gray30", fontface = "italic")

# Save the plot
ggsave("C4_Rice_percAreaAdvantage_line.jpeg", width = 9, height = 6, dpi = 300)
