# 1. Load necessary libraries
library(here)
library(tidyverse)
library(tools) # For file_path_sans_ext()

# 2. Specify paths and column names
input_dir <- here("output", "elsa_results")    # 📁 Folder with .lsa files
output_dir <- here("output", "wavepal_input")  # 📂 Folder to save the result
cols_to_extract <- c("X", "Y")                 # ⚙️ Columns you want

# 3. Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 4. Get the list of input files
all_files <- list.files(path = input_dir,
                        pattern = "\\off.lsa$",
                        full.names = TRUE)

# 5. Iterate through each file, process it, and save it individually 💾
purrr::walk(all_files, ~{
  # Construct the new output filename
  base_name <- file_path_sans_ext(basename(.x))
  output_path <- file.path(output_dir, paste0(base_name, "_wavepal.txt"))
  
  # Read the input file, select columns, and write to the new file
  read_tsv(.x, show_col_types = FALSE) %>%
    filter(Delay != 0,
           Q < 0.001,
           X == "ASV2") %>%
    select(all_of(cols_to_extract)) %>%
    write_tsv(output_path)
})

# Print a confirmation message
print(paste("Processing complete. Files saved in:", output_dir))