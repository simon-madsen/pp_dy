# Load necessary libraries
library(readr)
library(dplyr)
library(here)
library(ampvis2)
library(ggplot2)

# --- Part 1: Identifying the Matching Pairs ---

# Read and filter the main dataframe 'd' from the .lsa file
d <- read_tsv(here("output", "elsa_results", "Aalborg_W_rclr_abund.lsa")) %>%
  filter(Q <= 0.05,
         # Delay != 0,
         # LS >= 0,
         # SPCC > PCC,
         # SSCC > SCC,
         # SSCC > 0.8,
  )
d$diffPCC <- d$SPCC - d$PCC
d$diffSCC <- d$SSCC - d$SCC

# Read the periods dataframe from the .csv file
periods <- read.csv(here("data", "all_sig_periods.csv")) %>%
  filter(Period <= max(Period)-1 )

# Join 'd' with periods for both ASV X and ASV Y to create all combinations
d_with_all_period_combinations <- d %>%
  left_join(periods %>% select(asv, Period), by = c("X" = "asv")) %>%
  rename(Period_X = Period) %>%
  left_join(periods %>% select(asv, Period), by = c("Y" = "asv")) %>%
  rename(Period_Y = Period)

# Filter these combinations based on the 1% difference criteria
d_filtered <- d_with_all_period_combinations %>%
  filter(!is.na(Period_X) & !is.na(Period_Y)) %>%
  filter(abs(Period_X - Period_Y) < (0.02 * Period_X))

# --- NEW: Create a dataframe of unique matching pairs ---
# This dataframe will have two columns, X and Y, for easy iteration.
matching_pairs_df <- d_filtered %>%
  select(X, Y) %>%
  distinct()

print(paste("Found", nrow(matching_pairs_df), "unique ASV pairs with similar periods."))

# --- Part 2: Time Series Plotting with ampvis2 ---

# 1. Load the initial abundance data
d_initial <- readRDS(here("data", "d_initial_Simon_subset.rds"))%>%
  amp_filter_samples(SampleSite == "Aalborg W")

d_tax <- d_initial$tax
# 2. Source the rclr_transform function
source(here("R", "functions", "rclr_transform.R"))

d_rclr <- rclr_transform(d_initial)

# Filter non-transformed data to get a list of ASVs to keep
#d_otu_filter <- filter_otus(d_initial, filter_otus = 0.3)
#taxa_filter_vector <- as.vector(d_otu_filter$tax$OTU)

# Subset the RCLR transformed data
#d_processed <- amp_subset_taxa(d_rclr, tax_vector = taxa_filter_vector)

# create timeseries plot with amp_timerseries()

for (i in 1:2) {
  row <- matching_pairs_df[i,]
  p<- amp_time_series(d_rclr,
                  time_variable = "SampleDate",
                  tax_show = c(row[1], row[2]),
                  normalise = FALSE#,
                  #tax_aggregate = "Species"
                  ) 
  show(p)
}

unique_taxa <- as.data.frame(unique(matching_pairs_df))

unique_taxa_vector <- c(as.vector(unique_taxa[,1]), as.vector(unique_taxa[,2]))

d_processed <- amp_subset_taxa(d_rclr, tax_vector = unique_taxa_vector)

species_names <- d_processed$tax$Species
