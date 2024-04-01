# Load required packages
library(vcfR)

# Function to calculate ABBA/BABA statistics for a trio of variants
ABBA_BABA_function <- function(trio_variants, trio_names) {
  # Print trio names and column names for debugging
  print("Trio Names:")
  print(trio_names)
  print("Column Names:")
  print(colnames(trio_variants))
  # Check if trio_names contain duplicate populations
  if (length(unique(trio_names)) != 3) {
    stop("Trio names should consist of three distinct populations.")
  }
  
  # Initialize counts
  ABBA_count <- 0
  BABA_count <- 0
  
  # Iterate over each variant
  for (i in seq_along(trio_variants)) {
    # Extract allele counts for each population
    pop1_counts <- trio_variants[i, paste0(trio_names[1], "_counts")]
    print(pop1_counts)
    pop2_counts <- trio_variants[i, paste0(trio_names[2], "_counts")]
    pop3_counts <- trio_variants[i, paste0(trio_names[3], "_counts")]
    
    # Debugging: Print allele counts for each population
    cat("Pop1 Allele Counts:", pop1_counts, "\n")
    cat("Pop2 Allele Counts:", pop2_counts, "\n")
    cat("Pop3 Allele Counts:", pop3_counts, "\n")
    
    # Check for missing data
    if (anyNA(c(pop1_counts, pop2_counts, pop3_counts))) {
      next  # Skip variant if there is missing data
    }
    
    # Check for ABBA pattern
    if (pop1_counts == pop2_counts && pop1_counts != pop3_counts && pop3_counts == trio_variants[i, paste0(trio_names[1], "_ref")]) {
      ABBA_count <- ABBA_count + 1
    }
    
    # Check for BABA pattern
    if (pop1_counts != pop2_counts && pop1_counts == pop3_counts && pop2_counts == trio_variants[i, paste0(trio_names[1], "_ref")]) {
      BABA_count <- BABA_count + 1
    }
  }
  
  # Debugging: Print ABBA and BABA counts
  cat("ABBA Count:", ABBA_count, "\n")
  cat("BABA Count:", BABA_count, "\n")
  
  # Avoid division by zero
  if (ABBA_count == 0 || BABA_count == 0) {
    ABBA_BABA_statistic <- NaN
    p_value <- NaN
    Z_score <- NaN
  } else {
    # Calculate ABBA/BABA statistic (D statistic)
    ABBA_BABA_statistic <- (ABBA_count - BABA_count) / (ABBA_count + BABA_count)
    
    # Calculate p-value (assuming large sample size)
    p_value <- 2 * pnorm(-abs(ABBA_BABA_statistic))
    
    # Calculate Z-score
    Z_score <- ABBA_BABA_statistic / sqrt((ABBA_count + BABA_count) / 4)
  }
  
  # Create a data frame with results and trio names
  result <- data.frame(
    Trio = paste(trio_names[1], trio_names[2], trio_names[3], sep = "-"),
    ABBA_BABA_statistic = ABBA_BABA_statistic,
    p_value = p_value,
    Z_score = Z_score
  )
  
  return(result)
}

# Specify paths to VCF file and sample population mapping file
vcf_file <- "ABBA_BABA/EckSA_AUS_COX1_ForABBANoArborea.vcf"
sample_population_file <- "ABBA_BABA/sets3.txt"

# Read the file matching sequences to populations
population_file <- read.table(sample_population_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                              col.names = c("Sample", "Population"))
population_file <- population_file[-1,]

# Generate all unique trio combinations from the populations
unique_populations <- unique(population_file$Population)
trio_combinations <- combn(unique_populations, 3)

# Filter out combinations where Population1 < Population2 < Population3 to avoid duplicate combinations
valid_combinations <- apply(trio_combinations, 2, function(x) all(diff(match(x, unique_populations)) > 0))
trio_combinations <- trio_combinations[, valid_combinations]

# Convert trio combinations to a data frame with three columns
trio_combinations_df <- as.data.frame(t(trio_combinations))
colnames(trio_combinations_df) <- c("Population1", "Population2", "Population3")

# Modify trio names in the trio_combinations_df dataframe
# trio_combinations_df$Population1 <- paste0(trio_combinations_df$Population1, "_counts")
# trio_combinations_df$Population2 <- paste0(trio_combinations_df$Population2, "_counts")
# trio_combinations_df$Population3 <- paste0(trio_combinations_df$Population3, "_counts")

# Load VCF file
vcf <- read.vcfR(vcf_file)

# Initialize results list
results <- list()

# Extract genotype information for each trio and calculate ABBA/BABA statistics
results <- lapply(1:nrow(trio_combinations_df), function(i) {
  trio <- trio_combinations_df[i, ]
  # Extract genotype information for the trio
  print(paste("Processing Trio:", paste(trio, collapse = "-")))
  
  # Extract population names for the trio
  pop1 <- as.character(trio$Population1)
  pop2 <- as.character(trio$Population2)
  pop3 <- as.character(trio$Population3)
  
  # Get the sample names for each population from the population_data table
  pop1_samples <- population_file$Sample[population_file$Population == pop1]
  pop2_samples <- population_file$Sample[population_file$Population == pop2]
  pop3_samples <- population_file$Sample[population_file$Population == pop3]
  
  # Combine sample names into a vector for each population
  trio_samples <- c(pop1_samples, pop2_samples, pop3_samples)
  #print(trio_samples)
  
  # Extract column indices for the trio samples
  sample_indices <- match(trio_samples, colnames(vcf@gt))
  #print(sample_indices)
  
  # Subset the genotype matrix to include only the trio samples
  trio_genotypes <- vcf@gt[, sample_indices]
  print(trio_genotypes)
  
  # Calculate ABBA/BABA statistics for the trio
  result <- ABBA_BABA_function(trio_genotypes, c(pop1, pop2, pop3))
  
  return(result)
})

# Combine results into a data frame
results_df <- do.call(rbind, results)

# Print results data frame
print(results_df)
