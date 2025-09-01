library(dplyr)
library(ampvis2)
library(ggplot2)
library(here)

d <- readRDS(here("data", "d_initial_Simon_subset.rds"))

d_subset <-  amp_subset_samples(
  d,
  SampleSite %in% "Randers")

time_vector <- data.frame(values = d_subset$metadata$SampleDate)

# 2. Calculate difference and create color flag
time_vector <- time_vector %>%
  arrange(values) %>% # Ensure data is sorted by date
  mutate(
    days_from_start = as.numeric(difftime(values, min(values), units = "days")),
    diff_from_prior = as.numeric(difftime(values, lag(values), units = "days")),
    is_seven_days = ifelse(!is.na(diff_from_prior) & diff_from_prior == 7, "Yes", "No")
  )

# 3. Create the 1D plot with conditional color
ggplot(time_vector, aes(x = days_from_start, y = 0, color = is_seven_days)) +
  geom_point(size = 2) +
  scale_color_manual(
    name = "7 Days After Prior Point?", # Legend title
    values = c("Yes" = "green", "No" = "black")
  ) +
  labs(
    title = "Data Points with 7-Day Gaps",
    x = "Days from Start Date",
    subtitle = "Points in green are exactly 7 days after the previous point."
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
