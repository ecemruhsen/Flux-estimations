#Script for converting the MS measurements data format for Incawrapper

#Analysis column of the excel file has to be sorted before using the script!
#All sections have to be modified based on the number of replicates, naming of the sample column and the time points as well as time intervals
library(dplyr)
library(tidyr)
library(stringr)

# 1. LOAD METABOLITE ABBREVIATION LOOKUP TABLE
abbreviation_file <- "C:/Users/HP/Desktop/Master's Project/INCAWrapper/Glycolysis_TCA_abbreviations.csv"
abbrev_lookup <- read.csv(abbreviation_file)
abbrev_dict <- setNames(abbrev_lookup$abbreviation, abbrev_lookup$full_name)

data <- read.csv("C:/Users/HP/Desktop/Master's Project/INCAWrapper/Rikkedata_ctrl.csv")

data <- data %>%
  cbind(experiment_id = "exp1", .) %>%  
  rename(met_id = Analyte) %>% 
  mutate(unlabelled_atoms = "") %>%
  select(-FC, -Formula)

# APPLY ABBREVIATIONS AND FILTER OUT UNKNOWN METABOLITES
# Keep only metabolites that are in the abbreviation dictionary
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

#Analysis column should be all numeric
# Extract numbers after "glu"
df_with_labelled_ids$Analysis <- str_extract(df_with_labelled_ids$Analysis, "(?<=glu)\\d+")

# Read your sequence file
seq_data <- read.csv("C:/Users/HP/Desktop/Master's Project/INCAWrapper/Rikkeseq_ctrl.csv")

# Extract sample number from both dataframes and merge
df_with_time <- df_with_labelled_ids %>%
  mutate(
    sample_num = as.numeric(gsub("13Cglu", "", Analysis))
  ) %>%
  left_join(seq_data %>% 
              mutate(sample_num = as.numeric(gsub("13Cglu", "", sample))),
            by = "sample_num") %>%
  mutate(measurement_replicate = as.character((sample_num - 1) %% 3 + 1)) %>%
  select(-sample_num, -Analysis, -sample, -batch, -group, -order) 

# Calculate standard error using replicates before creating ms_id
df_with_std_error <- df_with_time %>%
  group_by(experiment_id, met_id, time, mass_isotope) %>%
  mutate(
    intensity_std_error = sd(intensity, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()

# Set ms_id to be exactly the same as met_id
df_with_time <- df_with_std_error %>%
  mutate(ms_id = met_id)  

# Condense the dataset by using average intensities and set all replicates to 1
df_condensed <- df_with_time %>%
  group_by(experiment_id, met_id, ms_id, time, mass_isotope, labelled_atom_ids, unlabelled_atoms) %>%
  summarize(
    intensity_avg = mean(intensity, na.rm = TRUE),
    intensity_std_error = first(intensity_std_error), # Keep the std error (should be same for all in group)
    .groups = "drop"
  ) %>%
  mutate(measurement_replicate = "1")  # Set all replicates to 1

df_final <- df_condensed %>%
  rename(intensity = intensity_avg) %>%
  relocate(experiment_id, met_id, ms_id, measurement_replicate,
           labelled_atom_ids, unlabelled_atoms, mass_isotope,
           intensity, intensity_std_error, time)

# Replace 0 values in intensity_std_error column with 0.02
df_final <- df_final %>%
  mutate(intensity_std_error = ifelse(intensity_std_error == 0, 0.02, intensity_std_error))

# Save the updated CSV file
write.csv(df_final, "C:/Users/HP/Desktop/Master's Project/INCAWrapper/Complete_TCA_glycolysis_ms file.csv", row.names = FALSE)
