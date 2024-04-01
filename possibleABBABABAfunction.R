library(dplyr)
library(tidyr)

setwd("/Users/pedromadeira/Library/CloudStorage/OneDrive-UniversidadedoAlgarve/EckSA_paper/EckSA_data")

#matrix <- fread("EcklSA_AUS_final_matrix2021_orderedR.txt", drop = 1) 
matrix <- read.table("EcklSA_AUS_final_matrix2021_orderedABBAtestPopSamplesGrouped.txt", header = TRUE)

# Add a new column "Population" based on the "Sample" column
matrix$Population <- gsub("\\d+", "", matrix$Population)

# Stack the data so that each row represents an allele for a sample at a specific locus
data_long <- reshape2::melt(matrix, id.vars = c("Population", "Sample"), variable.name = "Locus")

# Split each allele into two digits
data_long$Allele1 <- substr(data_long$value, 1, 3)
data_long$Allele2 <- substr(data_long$value, 4, 6)

# Remove the original value column
data_long <- data_long[, c("Population", "Sample", "Locus", "Allele1", "Allele2")]

# Step 1: Calculate the total count of each allele across all populations
total_counts <- data_long %>%
  group_by(Locus, Allele1) %>%
  summarise(TotalCount = n()) %>%
  ungroup() %>%
  rename(Allele = Allele1) %>%
  bind_rows(
    data_long %>%
      group_by(Locus, Allele2) %>%
      summarise(TotalCount = n()) %>%
      ungroup() %>%
      rename(Allele = Allele2)
  ) %>%
  group_by(Locus, Allele) %>%
  summarise(TotalCount = sum(TotalCount)) %>%
  ungroup()


# Step 2: Calculate the frequency of each allele within each population
allele_frequencies <- data_long %>%
  group_by(Population, Locus, Allele1) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  rename(Allele = Allele1) %>%
  bind_rows(
    data_long %>%
      group_by(Population, Locus, Allele2) %>%
      summarise(Count = n()) %>%
      ungroup() %>%
      rename(Allele = Allele2)
  ) %>%
  left_join(total_counts, by = c("Locus", "Allele")) %>%
  mutate(Frequency = if_else(is.na(Count), 0, Count) / TotalCount)

# Assuming you have allele frequency data in a data frame called allele_frequencies
# with columns: Population, Locus, Allele, Frequency

# Create all possible trios of populations
trio_combinations <- combn(unique(allele_frequencies$Population), 3)

# Function to calculate ABBA/BABA statistics for a given trio of populations and locus
calculate_ABBA_BABA <- function(population1, population2, population3, locus) {
  # Subset data for the given trio of populations and locus
  subset_data <- subset(allele_frequencies, Population %in% c(population1, population2, population3) & Locus == locus)
  
  # Calculate ABBA and BABA counts
  ABBA <- sum(subset_data$Frequency[subset_data$Population == population1] * 
                subset_data$Frequency[subset_data$Population == population2])
  BABA <- sum(subset_data$Frequency[subset_data$Population == population1] * 
                subset_data$Frequency[subset_data$Population == population3])
  
  # Calculate D statistic
  D <- (ABBA - BABA) / (ABBA + BABA)
  
  # Calculate Z-score and p-value (assuming normal distribution)
  Z_score <- D / sqrt((1 / ABBA) + (1 / BABA))
  p_value <- 2 * pnorm(-abs(Z_score))
  
  return(c(D, Z_score, p_value))
}

results2 <- calculate_ABBA_BABA("E.maxima","E.radiata","DHHybrids", "EC03")
# Apply the function to each trio of populations and locus
results2 <- lapply(1:ncol(trio_combinations), function(i) {
  trio <- trio_combinations[, i]
  locus_results <- lapply(unique(allele_frequencies$Locus), function(locus) {
    calculate_ABBA_BABA(trio[1], trio[2], trio[3], locus)
  })
  names(locus_results) <- unique(allele_frequencies$Locus)
  return(locus_results)
})

# Convert results to a data frame
results_df <- as.data.frame(t(simplify2array(results)))

# Add trio information
trio_names <- apply(trio_combinations, 2, paste, collapse = "_")
results_df$Population_Trio <- rep(trio_names, each = nrow(results_df))

# Output results data frame
print(results_df)