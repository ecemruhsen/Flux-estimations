library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(tidyverse)

IS_CONC_PMOL       <- 2500      
IS_CONC_Lactate    <- 2500

# =============================================================================
# 1. LOAD AND CLEAN DATA
# =============================================================================
DATA_FILE     <- read_csv("D:/Ecem export named/Cell/Hilic/Fluxomics_Ecem data cell (hilic)_ecem db for Flux all iso_intensity_corrected_2026-04-17-14-42-34+0200_ink coa.csv")
METADATA_FILE <- read_csv("D:/Ecem export named/Cell/Hilic/hilic_cell_combined-seq.csv")

# Remove QC and Blank from metadata
metadata_clean <- METADATA_FILE %>%
  filter(!grepl("QC|Blank", Group, ignore.case = TRUE))

DATA_FILE <- DATA_FILE %>%
  filter(Analysis %in% metadata_clean$Sample)

# =============================================================================
# 2. DEFINE METABOLITE — IS PAIRS
# =============================================================================
metabolite_istd <- data.frame(
  metabolite              = c("Alanine",  "Glucose CL",  "L-Glutamic acid"),
  internal_standard       = c("Alanine 15N", "Glucose d7 CL", "L-Glutamate 15N"),
  istd_isotopologue       = c("A+3", "A+0", "A+5"),
  stringsAsFactors = FALSE
)

METABOLITES <- metabolite_istd$metabolite
ISTDS       <- metabolite_istd$internal_standard

istd_isotopologue_map <- setNames(
  metabolite_istd$istd_isotopologue,
  metabolite_istd$internal_standard
)

# =============================================================================
# 3. EXTRACT DATA FOR REGULAR METABOLITES
# =============================================================================
a0_data_metabolites <- DATA_FILE %>%
  filter(Analyte %in% METABOLITES) %>%
  select(Analyte, Analysis, A0 = `A+0`) %>%
  mutate(A0 = as.numeric(A0))

a0_data_istds <- DATA_FILE %>%
  filter(Analyte %in% ISTDS) %>%
  rowwise() %>%
  mutate(
    required_isotopologue = istd_isotopologue_map[Analyte],
    IS_value = case_when(
      required_isotopologue == "A+0" ~ as.numeric(`A+0`),
      required_isotopologue == "A+3" ~ as.numeric(`A+3`),
      required_isotopologue == "A+5" ~ as.numeric(`A+5`),
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup() %>%
  select(Analyte, Analysis, A0 = IS_value)

# =============================================================================
# 4. PIVOT AND MERGE FOR REGULAR METABOLITES
# =============================================================================
metabolites_wide <- a0_data_metabolites %>%
  pivot_wider(names_from = Analyte, values_from = A0, values_fill = NA)

istds_wide <- a0_data_istds %>%
  pivot_wider(names_from = Analyte, values_from = A0, values_fill = NA)

merged_data <- metabolites_wide %>%
  left_join(istds_wide, by = "Analysis") %>%
  left_join(metadata_clean, by = c("Analysis" = "Sample")) %>%
  filter(Time == 0)

# =============================================================================
# 5. EXTRACT LACTATE DATA SEPARATELY
# =============================================================================

# Extract lactate A+0 and A+3
lactate_data <- DATA_FILE %>%
  filter(Analyte == "Lactic acid") %>%
  select(Analysis, A0 = `A+0`, A3 = `A+3`) %>%
  mutate(
    A0 = as.numeric(A0),
    A3 = as.numeric(A3)
  ) %>%
  left_join(metadata_clean, by = c("Analysis" = "Sample")) %>%
  filter(Time == 0)

# Check if lactate data exists
print("Lactate data at Time 0:")
print(lactate_data)

# =============================================================================
# 6. CALCULATE POOL SIZES
# =============================================================================

# Create empty vectors to store results
metabolite_name <- c()
condition <- c()
pool_value <- c()

# Calculate for Alanine
for(i in 1:nrow(merged_data)) {
  met <- "Alanine"
  cond <- merged_data$Group[i]
  num <- as.numeric(merged_data$Alanine[i])
  denom <- as.numeric(merged_data$`Alanine 15N`[i])
  if(!is.na(num) & !is.na(denom) & denom != 0) {
    pool <- (num / denom) * IS_CONC_PMOL 
    metabolite_name <- c(metabolite_name, met)
    condition <- c(condition, cond)
    pool_value <- c(pool_value, pool)
  }
}

# Calculate for Glucose
for(i in 1:nrow(merged_data)) {
  met <- "Glucose"
  cond <- merged_data$Group[i]
  num <- as.numeric(merged_data$`Glucose CL`[i])
  denom <- as.numeric(merged_data$`Glucose d7 CL`[i])
  if(!is.na(num) & !is.na(denom) & denom != 0) {
    pool <- (num / denom) * IS_CONC_PMOL 
    metabolite_name <- c(metabolite_name, met)
    condition <- c(condition, cond)
    pool_value <- c(pool_value, pool)
  }
}

# Calculate for Glutamic acid
for(i in 1:nrow(merged_data)) {
  met <- "Glutamic_acid"
  cond <- merged_data$Group[i]
  num <- as.numeric(merged_data$`L-Glutamic acid`[i])
  denom <- as.numeric(merged_data$`L-Glutamate 15N`[i])
  if(!is.na(num) & !is.na(denom) & denom != 0) {
    pool <- (num / denom) * IS_CONC_PMOL 
    metabolite_name <- c(metabolite_name, met)
    condition <- c(condition, cond)
    pool_value <- c(pool_value, pool)
  }
}

# Calculate for LACTATE - Using A+0 / A+3 (self-normalization)
if(nrow(lactate_data) > 0) {
  for(i in 1:nrow(lactate_data)) {
    met <- "Lactic_acid"
    cond <- lactate_data$Group[i]
    num <- lactate_data$A0[i]   # A+0
    denom <- lactate_data$A3[i]  # A+3
    if(!is.na(num) & !is.na(denom) & denom != 0) {
      pool <- (num / denom) * IS_CONC_Lactate 
      metabolite_name <- c(metabolite_name, met)
      condition <- c(condition, cond)
      pool_value <- c(pool_value, pool)
    }
  }
} else {
  print("WARNING: No lactate data found at Time 0!")
}

# Create data frame from the vectors
results_df <- data.frame(
  metabolite = metabolite_name,
  Group = condition,
  pool_size = pool_value,
  stringsAsFactors = FALSE
)

print("Individual pool size values:")
print(results_df)

# =============================================================================
# 7. CALCULATE SUMMARY STATISTICS
# =============================================================================

# Get unique metabolites and groups
unique_metabolites <- unique(results_df$metabolite)
unique_groups <- unique(results_df$Group)

# Create summary dataframe
summary_list <- list()

for(met in unique_metabolites) {
  for(grp in unique_groups) {
    values <- results_df$pool_size[results_df$metabolite == met & results_df$Group == grp]
    
    if(length(values) > 0) {
      n_vals <- length(values)
      mean_val <- mean(values)
      sd_val <- sd(values)
      se_val <- sd_val / sqrt(n_vals)
      
      # Create met_id
      met_id <- case_when(
        met == "Alanine" ~ "Ala",
        met == "Glucose" ~ "Glc",
        met == "Glutamic_acid" ~ "Glu",
        met == "Lactic_acid" ~ "Lac"
      )
      
      summary_list[[length(summary_list) + 1]] <- data.frame(
        met_id = met_id,
        pool_size = round(mean_val, 4),
        Group = grp,
        pool_size_std_error = round(se_val, 4),
        sd = round(sd_val, 6),
        n = n_vals,
        stringsAsFactors = FALSE
      )
    }
  }
}

pool_size_summary_t0 <- do.call(rbind, summary_list) %>%
  arrange(met_id, Group)

print("=== POOL SIZE AT TIME 0 WITH STANDARD ERROR ===")
print(pool_size_summary_t0)

# =============================================================================
# 8. CREATE VISUALIZATION
# =============================================================================
pool_plot_t0 <- pool_size_summary_t0 %>%
  ggplot(aes(x = met_id, y = pool_size, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = pool_size - pool_size_std_error, 
                    ymax = pool_size + pool_size_std_error),
                position = position_dodge(0.7), width = 0.25) +
  labs(
    title = "Intracellular Metabolite Pool Sizes at Time 0",
    subtitle = "Mean ± SEM (n=3 per condition)",
    x = "Metabolite",
    y = "Pool Size (pmol/10⁶ cells)",
    fill = "Condition"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

print(pool_plot_t0)

