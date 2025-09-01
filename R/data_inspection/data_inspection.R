library(ampvis2)
library(here)

#setwd("C:/Users/com/Documents/Msc project/msc_R")

source("functions/rclr_transform.R")

d <- readRDS(here("data", "d_initial_Simon_subset.rds"))

d_subset <-  amp_subset_samples(
  d,
  SampleSite %in% "Randers")

d_subset_rclr <- rclr_transform(d_subset)

p_gg <- amp_timeseries(
  d_subset_rclr,
  tax_show = c("ASV1360", "ASV461"),
  normalise = FALSE
  #plotly = TRUE
) 

# 2. Define your date limits using as.Date()
date_limits <- c(as.Date("2020-02-01"), as.Date("2020-04-30"))

# 3. Add the scale_x_date() layer with the limits
p_zoomed <- p_gg +
  scale_x_date(
    limits = date_limits,
    date_breaks = "2 months",      # Adjust breaks for the new range
    date_labels = "%b %Y"          # Format the labels
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 4. Print the plot with the new date range
print(p_zoomed)
