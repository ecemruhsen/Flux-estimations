# ============================================================================
# COMPLETE R SCRIPT FOR METABOLITE POOL SIZE ANALYSIS
# Metabolites: Pyr, Asp, Mal, Fum, Cit, Glu
# ============================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(stringr)

# ============================================================================
# 1. READ AND PREPARE DATA
# ============================================================================

# Read the main data file
data <- read.csv("D:/Ecem export named/Cell/Vehicle_combined_relative_cell.csv", stringsAsFactors = FALSE)

# Read the sequence file with time points
seq_data <- read.csv("D:/Ecem export named/Cell/Vehicle_hilicRP_combined_seq.csv", stringsAsFactors = FALSE)

# ============================================================================
# 2. DATA CLEANING AND RESHAPING
# ============================================================================

# Remove rows with blank01, blank02, or QC in Analysis column if present
data_filtered <- data %>%
  filter(!grepl("blank01|blank02|blank03|QC", Analysis, ignore.case = TRUE))

# Reshape from wide to long format (isotopologues as rows)
df_long <- data_filtered %>%
  pivot_longer(
    cols = starts_with("A."), 
    names_to = "mass_isotope",                   
    values_to = "intensity"                  
  ) %>%
  mutate(mass_isotope = as.numeric(gsub("A\\.", "", mass_isotope))) %>%
  filter(!is.na(intensity) & intensity > 0)  # Keep only positive intensities

# Extract numeric sample number from Analysis column
df_long <- df_long %>%
  mutate(
    sample_num = as.numeric(str_extract(Analysis, "\\d+"))
  )

# Merge with sequence data to get time points
df_with_time <- df_long %>%
  left_join(seq_data, by = c("Analysis" = "Sample")) %>%
  rename(Time = Time) %>%
  # Create technical replicate number (1,2,3 per time point)
  group_by(Time) %>%
  mutate(technical_replicate = as.character(dense_rank(sample_num))) %>%
  ungroup()

# ============================================================================
# 3. FILTER FOR METABOLITES OF INTEREST
# ============================================================================

# Define your metabolites of interest
metabolites_of_interest <- c("Pyruvate", "L-Aspartic Acid", "Malate", 
                             "Fumaric acid (RP)", "(iso)Citrate", "L-Glutamic acid")

# Filter for these metabolites
df_filtered <- df_with_time %>%
  filter(Analyte %in% metabolites_of_interest)

# ============================================================================
# 4. SUM ISOTOPOLOGUES PER METABOLITE, TIME, AND REPLICATE
# ============================================================================

summed_data <- df_filtered %>%
  group_by(Analyte, Time, technical_replicate) %>%
  summarise(
    total_intensity = sum(intensity, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================================
# 5. NORMALIZE TO MEAN AT TIME ZERO
# ============================================================================

normalized_data <- summed_data %>%
  group_by(Analyte) %>%
  mutate(
    # Calculate mean at time 0 across the 3 replicates
    mean_t0 = mean(total_intensity[Time == 0], na.rm = TRUE),
    # Normalize: divide by mean at time 0
    normalized_intensity = total_intensity / mean_t0
  ) %>%
  ungroup()

# ============================================================================
# 6. CREATE PLOTS
# ============================================================================

# Define cleaner names for plotting
metabolite_labels <- c(
  "Pyruvate" = "Pyruvate (Pyr)",
  "L-Aspartic Acid" = "Aspartate (Asp)",
  "Malate" = "Malate (Mal)",
  "Fumaric acid (RP)" = "Fumarate (Fum)",
  "(iso)Citrate" = "Citrate (Cit)",
  "L-Glutamic acid" = "Glutamate (Glu)"
)

# Add a clean name column for plotting
normalized_data <- normalized_data %>%
  mutate(metabolite_clean = recode(Analyte, !!!metabolite_labels))

# ============================================================================
# 7. DEFINE TIME POINTS AND FACTOR LEVELS
# ============================================================================
# Ensure Time is numeric
normalized_data <- normalized_data %>%
  mutate(Time = as.numeric(Time))

# Define the order of metabolites for faceting
metabolite_order <- c("Pyruvate (Pyr)", "Aspartate (Asp)", "Malate (Mal)",
                      "Fumarate (Fum)", "Citrate (Cit)", "Glutamate (Glu)")

normalized_data <- normalized_data %>%
  mutate(metabolite_clean = factor(metabolite_clean, levels = metabolite_order))

# ============================================================================
# 8. CREATE POOL SIZE PLOT
# ============================================================================

pool_size_plot <- ggplot(normalized_data, 
                         aes(x = Time, y = normalized_intensity)) +
  
  # Horizontal reference line at y = 1
  geom_hline(yintercept = 1, 
             linetype = "solid", 
             color = "black", 
             linewidth = 0.4) +
  
  # Individual replicate points with jitter to avoid overlap
  geom_point(size = 1.8, 
             color = "black",
             alpha = 0.8,
             position = position_jitter(width = 2, height = 0)) +  # Added jitter
  # width = 2 spreads points horizontally, height = 0 keeps y values accurate
  
  # Facet by metabolite
  facet_wrap(~ metabolite_clean, 
             ncol = 2,  # Changed from 3 to 2 columns for wider individual plots
             scales = "free_y") +
  
  # X axis breaks matching your time points
  scale_x_continuous(
    breaks = sort(unique(normalized_data$Time)),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  
  # Y axis
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  
  # Labels
  labs(
    title = "Metabolite Pool Sizes Over Time",
    subtitle = "Vehicle",
    x = "Time (minutes)",
    y = "Normalized Signal"
  ) +
  
  # Theme matching the reference figure style
  theme_bw() +
  theme(
    # Facet strip styling
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    
    # Axis styling
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9, color = "black"),
    axis.text.x = element_text(angle = 0),
    
    # Panel styling - wider panels
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    
    # Title
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    
    # No legend needed
    legend.position = "none"
  )

print(pool_size_plot)

# ============================================================================
# 9. SAVE PLOT
# ============================================================================
ggsave(
  filename = "pool_size_normalized_t0.pdf",
  plot = pool_size_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "pool_size_normalized_t0.png",
  plot = pool_size_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================================
# 10. OPTIONAL - PRINT SUMMARY TABLE FOR CHECKING
# ============================================================================
summary_check <- normalized_data %>%
  group_by(metabolite_clean, Time) %>%
  summarise(
    n_replicates = n(),
    mean_normalized = round(mean(normalized_intensity, na.rm = TRUE), 3),
    sd_normalized = round(sd(normalized_intensity, na.rm = TRUE), 3),
    .groups = "drop"
  )

print(summary_check)
