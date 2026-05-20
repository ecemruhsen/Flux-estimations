### MS data converter for INCAWrapper ###

#All sections have to be modified based on the number of replicates, naming of the sample column and the time points 
library(dplyr)
library(tidyr)
library(stringr)

# LOAD METABOLITE ABBREVIATION LOOKUP TABLE
abbreviation_file <- "C:/Users/HP/Desktop/Master's Project/INCAWrapper/Metabolite_abbreviations.csv"
abbrev_lookup <- read.csv(abbreviation_file)
abbrev_dict <- setNames(abbrev_lookup$abbreviation, abbrev_lookup$full_name)

data <- read.csv("D:/Ecem export named/Cell/Vehicle_combined_relative_cell.csv")

data <- data %>%
  cbind(experiment_id = "exp1", .) %>%  
  rename(met_id = Analyte) %>% 
  mutate(unlabelled_atoms = "") %>%
  select(-FC, -Formula)

# APPLY ABBREVIATIONS AND FILTER OUT UNRELATED METABOLITES
data <- data %>%
  filter(met_id %in% names(abbrev_dict)) %>%
  mutate(met_id = abbrev_dict[met_id])

#Using Analysis column to create replicates section
#Remove rows with blank01, blank02, or QC in analysis column
data_filtered <- data %>%
  filter(!grepl("blank01|blank02|blank03", Analysis, ignore.case = TRUE))

#Creating the isotopologues and the intensity columns
df_transposed <- data_filtered %>%
  pivot_longer(
    cols = starts_with("A."), 
    names_to = "mass_isotope",                   
    values_to = "intensity"                  
  ) %>%
  mutate(mass_isotope = as.numeric(gsub("A\\.", "", mass_isotope))) %>%
  filter(!is.na(intensity))

# Calculate maximum mass_isotope for each metabolite and create labelled_atom_ids
metabolite_isotope_ranges <- df_transposed %>%
  group_by(met_id) %>%
  summarize(max_mass_isotope = max(mass_isotope, na.rm = TRUE)) %>%
  mutate(labelled_atom_ids = purrr::map_chr(max_mass_isotope, function(x) {
    if (x > 0) {
      paste0("[", paste0(1:x, collapse = ","), "]")
    } else {
      "[]"
    }
  })) %>%
  select(-max_mass_isotope)

# Join the labelled_atom_ids back to the main dataframe
df_with_labelled_ids <- df_transposed %>%
  left_join(metabolite_isotope_ranges, by = "met_id") %>%
  relocate(labelled_atom_ids, .before = unlabelled_atoms)

# Extract numeric portion from Analysis column to create sample_num
df_with_labelled_ids <- df_with_labelled_ids %>%
  mutate(sample_num = as.numeric(str_extract(Analysis, "(?<=hilic-|RPneg-)\\d+")))

# Read sequence file
seq_data <- read.csv("D:/Ecem export named/Cell/Vehicle_hilicRP_combined_seq.csv")

# Create a clean lookup table with unique sample_num to Time mapping
time_lookup <- seq_data %>%
  mutate(sample_num = as.numeric(str_extract(Sample, "(?<=hilic-|RPneg-)\\d+"))) %>%
  select(sample_num, Time) %>%
  distinct(sample_num, .keep_all = TRUE)

# Now join - sample_num exists in df_with_labelled_ids
df_with_time <- df_with_labelled_ids %>%
  left_join(time_lookup, by = "sample_num") %>%
  filter(!is.na(Time)) %>%
  mutate(measurement_replicate = as.character((sample_num - 1) %% 3 + 1)) %>%
  select(-Analysis, -sample_num)

# Assign experiment_id based on measurement_replicate
df_with_experiment <- df_with_time %>%
  mutate(experiment_id = case_when(
    measurement_replicate == "1" ~ "exp1",
    measurement_replicate == "2" ~ "exp2",
    measurement_replicate == "3" ~ "exp3",
    TRUE ~ "exp1"
  ))

# Set default standard error to 0.002
df_with_std_error <- df_with_experiment %>%
  mutate(intensity_std_error = 0.002)

# Set ms_id to be the same as met_id
df_with_ms_id <- df_with_std_error %>%
  mutate(ms_id = met_id)

#Create dataframe with correct column names
df_final <- df_with_ms_id %>%
  rename(time = Time) %>%
  select(experiment_id, met_id, ms_id, measurement_replicate,
         labelled_atom_ids, unlabelled_atoms, mass_isotope,
         intensity, intensity_std_error, time) %>%
  arrange(experiment_id, met_id, time, mass_isotope)

# Replace 0 values with NaN
df_final <- df_final %>%
  mutate(intensity = ifelse(intensity == 0, NaN, intensity))

# Verify uniqueness
duplicate_check <- df_final %>%
  group_by(experiment_id, met_id, time, mass_isotope) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count > 1)

if(nrow(duplicate_check) > 0) {
  print("ERROR: Duplicates found!")
  print(duplicate_check)
} else {
  print("No duplicates found, each (experiment, metabolite, time, mass_isotope) is unique.")
}

# Print summary
print(paste("Total rows:", nrow(df_final)))
print(paste("Experiments:", paste(unique(df_final$experiment_id), collapse=", ")))
print(paste("Time points:", paste(sort(unique(df_final$time)), collapse=", ")))

# Save the CSV files
write.csv(df_final, "C:/Users/HP/Desktop/Master's Project/INCAWrapper/Vehicle_Parallel_MS.csv", row.names = FALSE)
