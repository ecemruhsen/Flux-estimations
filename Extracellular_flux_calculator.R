### Extracellular Flux Calculator ###

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(tidyverse)

DATA_FILE <- read_csv("D:/Ecem export named/Media/Fluxomics_Ecem media data (hilic)_ecem db for tasq1 media flux1_intensity_corrected_2026-04-16-16-58-55+0200.csv")
METADATA_FILE <- read_csv("D:/Ecem export named/Media/media_seq_combined.csv")

# Remove rows containing "QC" or "Blank" 
metadata_clean <- METADATA_FILE %>%
  filter(!grepl("QC|Blank", Group))
DATA_FILE <- DATA_FILE %>%
  filter(Analysis %in% metadata_clean$Sample)

# Define metabolites and their corresponding internal standards with specific isotopologues
metabolite_istd <- data.frame(
  metabolite = c("Alanine", "Glucose CL", "Lactic acid", "L-Glutamic acid"),
  internal_standard = c("Alanine 13C15N", "Glucose d7 CL", "Lactic acid C3", "L-Glutamic acid 13C15N work"),
  metabolite_isotopologue = c("A+0", "A+0", "A+0", "A+0"),
  istd_isotopologue = c("A+3", "A+0", "A+3", "A+5"),
  stringsAsFactors = FALSE
)

METABOLITES <- metabolite_istd$metabolite
ISTDS <- metabolite_istd$internal_standard

# Create a named vector for isotopologue mapping
istd_isotopologue_map <- setNames(
  metabolite_istd$istd_isotopologue,
  metabolite_istd$internal_standard
)

# Known concentration of the spiked standards
IS_CONC_PMOL <- 2500
IS_CONC_Lactate <- 2050

# Extract A+0 for metabolites
a0_data_metabolites <- DATA_FILE %>%
  filter(Analyte %in% METABOLITES,
         Analysis %in% metadata_clean$Sample) %>%
  select(Analyte, Analysis, A0 = `A+0`) %>%
  mutate(A0 = as.numeric(A0),
         A0 = ifelse(is.na(A0) | A0 == 0, NA_real_, A0))

# Extract internal standards with their specific isotopologues
a0_data_istds <- DATA_FILE %>%
  filter(Analyte %in% ISTDS,
         Analysis %in% metadata_clean$Sample) %>%
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
  select(Analyte, Analysis, A0 = IS_value) %>%
  mutate(A0 = ifelse(is.na(A0) | A0 == 0, NA_real_, A0))

# Pivot metabolites
metabolites_wide <- a0_data_metabolites %>%
  pivot_wider(
    names_from = Analyte,
    values_from = A0,
    values_fill = NA
  )

# Pivot internal standards
istds_wide <- a0_data_istds %>%
  pivot_wider(
    names_from = Analyte,
    values_from = A0,
    values_fill = NA
  )

# Merge metabolites and internal standards
merged_data <- metabolites_wide %>%
  left_join(istds_wide, by = "Analysis")

# Add metadata to merged_data
merged_data <- merged_data %>%
  left_join(metadata_clean, by = c("Analysis" = "Sample"))


# Add replicate numbers for each condition and time point 
merged_data <- merged_data %>%
  group_by(Group, Time) %>%
  mutate(replicate = row_number()) %>%
  ungroup()

# SPECIAL HANDLING FOR LACTATE: Replicate-specific T0 normalization
# Calculate T0 factors for EACH replicate 
lactate_t0_factors <- merged_data %>%
  filter(Time == 0) %>%
  mutate(
    Y = (`Lactic acid` / `Lactic acid C3`) * IS_CONC_PMOL,
    lactate_A0_t0 = `Lactic acid`
  ) %>%
  select(Group, replicate, Y, lactate_A0_t0)  

print(lactate_t0_factors)

# Calculate absolute concentrations with replicate-specific normalization
calculations <- merged_data %>%
  left_join(lactate_t0_factors, by = c("Group", "replicate")) %>%
  mutate(
    Alanine_conc_pmol = (Alanine / `Alanine 13C15N`) * IS_CONC_PMOL,
  
    Glucose_CL_conc_pmol = (`Glucose CL` / `Glucose d7 CL`) * IS_CONC_PMOL,
    
    Glutamic_acid_conc_pmol = (`L-Glutamic acid` / `L-Glutamic acid 13C15N work`) * IS_CONC_PMOL,
    
    Lactic_acid_conc_pmol = case_when(
      # At T0: use the standard IS calculation
      Time == 0 ~ (`Lactic acid` / `Lactic acid C3`) * IS_CONC_Lactate,
      # For other time points: Y × (A+0 at time X / A+0 at T0) for this specific replicate
      Time != 0 & !is.na(Y) ~ Y * (`Lactic acid` / lactate_A0_t0),
      TRUE ~ NA_real_
    )
  ) %>%
  select(Analysis, Alanine_conc_pmol, Glucose_CL_conc_pmol, 
         Glutamic_acid_conc_pmol, Lactic_acid_conc_pmol)

# Rename for consistency
names(calculations) <- c("Sample_ID", "Alanine_pmol", "Glucose_CL_pmol", 
                         "Glutamic_acid_pmol", "Lactic_acid_pmol")

# Merge with metadata
results <- calculations %>%
  left_join(metadata_clean, by = c("Sample_ID" = "Sample"))

# Add replicate numbers to results
results_with_replicates <- results %>%
  group_by(condition = Group, time = Time) %>%
  mutate(replicate = row_number()) %>%
  ungroup()


write.csv(results_with_replicates, "C:/Users/HP/Downloads/Concantrations2_MS.csv", row.names = FALSE)

# Reshape to long format
results_long <- results_with_replicates %>%
  pivot_longer(
    cols = ends_with("_pmol"),
    names_to = "metabolite",
    values_to = "conc",
    names_pattern = "(.*)_pmol"
  ) %>%
  filter(!is.na(conc)) %>%
  select(condition, metabolite, replicate, time, conc) %>%
  mutate(
    time = as.numeric(time),
    conc = as.numeric(conc)
  ) %>%
  arrange(condition, metabolite, replicate, time)


# Generate trend plot
trend_plot <- results_long %>%
  ggplot(aes(x = time, y = conc, color = condition, group = interaction(condition, replicate))) +
  geom_point(size = 2.5) +
  geom_line(size = 0.8) +
  facet_wrap(~metabolite, scales = "free_y") +
  scale_x_continuous(breaks = sort(unique(results_long$time))) +
  scale_color_manual(values = c("inhibitor" = "#FFB6C1", "vehicle" = "#90EE90")) +
  labs(
    title = "Metabolite Concentration Trends: Vehicle vs Inhibitor",
    x = "Time (minutes)",
    y = "Concentration (pmol/sample)",
    color = "Condition"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

print(trend_plot)

# FLUX CALCULATION
calculate_flux <- function(df) {
  flux_per_rep <- df %>%
    group_by(condition, metabolite, replicate) %>%
    arrange(time, .by_group = TRUE) %>%
    summarise(
      flux = coef(lm(conc ~ time))[["time"]],
      .groups = "drop"
    )
  
  flux_summary <- flux_per_rep %>%
    group_by(condition, metabolite) %>%
    summarise(
      mean_flux = mean(flux, na.rm = TRUE),
      se_flux   = sd(flux, na.rm = TRUE) / sqrt(n()),
      n         = n(),
      .groups   = "drop"
    )
  
  return(list(
    per_replicate = flux_per_rep,
    summary       = flux_summary
  ))
}

flux_results <- calculate_flux(results_long)
print(flux_results$summary)

# Generate flux barplots
condition_names <- unique(flux_results$summary$condition)

flux_plot <- flux_results$summary %>%
  ggplot(aes(x = metabolite, y = mean_flux, fill = condition)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.7), 
           width = 0.65, 
           color = "black") +
  geom_errorbar(aes(ymin = mean_flux - se_flux, 
                    ymax = mean_flux + se_flux),
                position = position_dodge(width = 0.7),
                width = 0.25,
                linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = setNames(c("#FFB6C1", "#90EE90"), condition_names)) +
  labs(
    title = "Extracellular Flux Rates",
    subtitle = "Lactate: Replicate-specific T0 normalization",
    x = "Metabolite",
    y = "Flux (pmol/10⁶ cells/time unit)",
    fill = "Condition"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    legend.position = "top"
  )

print(flux_plot)

# Create boxplot faceted
flux_boxplot_faceted <- flux_results$per_replicate %>%
  ggplot(aes(x = condition, y = flux, fill = condition)) +
  geom_boxplot(alpha = 0.7, width = 0.6, outlier.size = 2) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "gray40") +
  facet_wrap(~metabolite, scales = "free_y") +
  scale_fill_manual(values = setNames(c("#FFB6C1", "#90EE90"), condition_names)) +
  labs(
    title = "Flux Rates by Metabolite",
    x = "Condition",
    y = "Flux (pmol/10⁶ cells/time unit)",
    fill = "Condition"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "none"
  )

print(flux_boxplot_faceted)
