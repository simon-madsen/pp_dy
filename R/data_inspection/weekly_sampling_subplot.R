library(dplyr)
library(ampvis2)
library(ggplot2)
library(here)

d <- readRDS(here("data", "d_initial_Simon_subset.rds"))

# Get unique sample sites
sample_sites <- unique(d$metadata$SampleSite)

# Create an empty list to store data for each site
plot_data_list <- list()

# Loop over each sample site to prepare the data
for (site in sample_sites) {
  
  # Filter the data for the current site
  d_subset <- amp_filter_samples(d, SampleSite == site)
  
  # Continue if there is data for the site
  if (nrow(d_subset$metadata) > 1) {
    
    time_vector <- data.frame(values = d_subset$metadata$SampleDate,
                              SampleSite = site)
    
    # Calculate difference and create color flag
    time_vector <- time_vector %>%
      arrange(values) %>% # Ensure data is sorted by date
      mutate(
        days_from_start = as.numeric(difftime(values, min(values), units = "days")),
        diff_from_prior = as.numeric(difftime(values, lag(values), units = "days")),
        is_seven_days = ifelse(!is.na(diff_from_prior) & diff_from_prior <= 7, "Yes", "No")
      )
    
    # Calculate the ratio
    total_days <- max(time_vector$days_from_start)
    num_samples <- nrow(time_vector)
    ratio <- round(total_days / num_samples, 2)
    
    # Add annotation text
    time_vector$annotation <- paste("Avg. days between samples:", ratio)
    
    plot_data_list[[site]] <- time_vector
    
  } else {
    message(paste("Skipping site", site, "due to insufficient data points."))
  }
}

# Combine all data into a single data frame
plot_data <- do.call(rbind, plot_data_list)

# Create the plot with facets
ggplot(plot_data, aes(x = days_from_start, y = 0, color = is_seven_days)) +
  geom_point(size = 2, shape = 18) +
  scale_color_manual(
    name = "\u22647 Days After Prior Point?", # Legend title
    values = c("Yes" = "chartreuse2", "No" = "black")
  ) +
  labs(
    title = "Data Points with less than or equal to 7-Day Gaps by Site",
    x = "Days from Start Date"
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  facet_wrap(~ SampleSite, scales = "free_x") +
  geom_text(aes(x = Inf, y = Inf, label = annotation), 
            hjust = 1.1, vjust = 2, size = 3, color = "black",
            check_overlap = TRUE)

ggsave(
  here("output", "weekly_sampling_subplot.jpg"),
  dpi = 300
)
