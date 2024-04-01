library(dplyr)
library(tidyr)

setwd("/Users/pedromadeira/Library/CloudStorage/OneDrive-UniversidadedoAlgarve/EckSA_paper/EckSA_data")

matrix <- read.table("EcklSA_AUS_final_matrix2021_orderedABBAtestK4.txt", header = TRUE)

# Add a new column "Population" based on the "Sample" column
#matrix$Population <- gsub("\\d+", "", matrix$Population) #use only if the following function doesn't work properly

# Stack the data so that each row represents an allele for a sample at a specific locus
data_long <- reshape2::melt(matrix, id.vars = c("Population", "Sample"), variable.name = "Locus")

# Split each allele into two digits
data_long$Allele1 <- substr(data_long$value, 1, 3)
data_long$Allele2 <- substr(data_long$value, 4, 6)

# Remove the original value column
data_long <- data_long[, c("Population", "Sample", "Locus", "Allele1", "Allele2")]

#Create frequency tables from the allele data
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

#Perform ABBA/BABA for trios of populations/sets
# Create all combinations of trios of populations and loci
trio_combinations <- expand.grid(Population1 = unique(allele_frequencies$Population),
                                 Population2 = unique(allele_frequencies$Population),
                                 Population3 = unique(allele_frequencies$Population),
                                 Locus = unique(allele_frequencies$Locus))

# Convert factor variables to character vectors
trio_combinations$Population1 <- as.character(trio_combinations$Population1)
trio_combinations$Population2 <- as.character(trio_combinations$Population2)
trio_combinations$Population3 <- as.character(trio_combinations$Population3)

# Filter out combinations where Population1, Population2, and Population3 are not all distinct
trio_combinations <- trio_combinations[trio_combinations$Population1 != trio_combinations$Population2 &
                                         trio_combinations$Population1 != trio_combinations$Population3 &
                                         trio_combinations$Population2 != trio_combinations$Population3, ]

# Filter out combinations where Population1 < Population2 < Population3 to avoid duplicate combinations
trio_combinations <- trio_combinations[trio_combinations$Population1 < trio_combinations$Population2 &
                                         trio_combinations$Population2 < trio_combinations$Population3, ]

# Initialize an empty list to store results
results_trios <- list()

# Loop through each combination of trios of populations and loci
for (i in 1:nrow(trio_combinations)) {
  population1 <- trio_combinations$Population1[i]
  population2 <- trio_combinations$Population2[i]
  population3 <- trio_combinations$Population3[i]
  locus <- trio_combinations$Locus[i]
  
  # Subset the data for the specific populations and locus
  subset_data <- allele_frequencies[allele_frequencies$Population %in% c(population1, population2, population3) & 
                                      allele_frequencies$Locus == locus, ]
  
  # Create subsets for each population
  subset_data_p1 <- subset_data[subset_data$Population == population1, ]
  subset_data_p2 <- subset_data[subset_data$Population == population2, ]
  subset_data_p3 <- subset_data[subset_data$Population == population3, ]
  
  # Calculate D statistic for the trio of populations
  D_trio <- (sum(subset_data_p1$Frequency) - sum(subset_data_p2$Frequency)) / 
    (sum(subset_data_p1$Frequency) + sum(subset_data_p2$Frequency) + sum(subset_data_p3$Frequency))
  
  # Calculate variance
  variance_trio <- sum(subset_data_p1$Frequency^2) / nrow(subset_data_p1) +
    sum(subset_data_p2$Frequency^2) / nrow(subset_data_p2) +
    sum(subset_data_p3$Frequency^2) / nrow(subset_data_p3)
  
  # Calculate Z-score
  Z_score_trio <- D_trio / sqrt(variance_trio)
  
  # Calculate p-value
  p_value_trio <- 2 * pnorm(-abs(Z_score_trio))
  
  # Store results in a data frame
  results_trios[[i]] <- data.frame(Population1 = population1, Population2 = population2, Population3 = population3, Locus = locus,
                                   D_trio = D_trio, Z_score_trio = Z_score_trio, p_value_trio = p_value_trio)
}

# Combine results into a single data frame
results_trios <- do.call(rbind, results_trios)

#Change scientifc notation to regular notation, in case it helps reading the data
results_trios["p_value_trio"]<- format(results_trios$p_value_trio, scientific = FALSE)

#filter the results from result_trios to include only statistically significant samples
relevant_trios <- results_trios %>% filter(p_value_trio < 0.055) %>% filter(abs(Z_score_trio) > 3)

#Write the table with the relevant results
write.csv(relavant_trios, "ABBA_BABA/DstatsAllTriosAllLociSpecies2.csv")

#If there are many populations/sets the data can become hard to read in the full table
#As such, we can filter the tables to look only into our populations of interest
#To ease up the job, we can filter the tables in two ways, one including only positive D stats
#other including only negative D stats (or Z-scores). 
#This works because a positive D stat shows introgression between Pop2 and Pop3, while a negative
#D stat shows introgression between Pop1 and Pop3. So we combine this information to create tables
#where we keep positive D stats and rows where our population of interest is in Pop2 or Pop3 and
#where we keep negative D stats and rows where our population of interest is in Pop1 or Pop3
#Again, if there are many populations we can have multiple trio combinations, making it harder 
#to parse the information, so we can also make it so we only keep one row for each time the
#combination appears. Here we do so by filtering for the instance where the absolute Z-score is highest

#Filter table to include positive Dstats and our population of interest in column Pop2 or Pop3
filtered_tableDHhp <- relevant_trios %>%
  filter(Population2 == "DHh" | Population3 == "DHh") %>% #Keep only rows where ERSDH appears in either Pop2 or Pop3
  filter(Z_score_trio > 0) %>% #Keep only positive Dstats/Z-scores
  arrange(desc(Z_score_trio)) %>% #Arrange the Z-scores from highest to lowest
  distinct(Population2, Population3, .keep_all = TRUE) #Keep only the first instance of the combination, which will have the highest Z-score due to the line above

#Filter table to include negative Dstats and our population of interest in column Pop1 or Pop3
filtered_tableDHhm <- relevant_trios %>%
  filter(Population1 == "DHh" | Population3 == "DHh") %>% #Keep only rows where ERSDH is in Pop1 or Pop3
  filter(Z_score_trio < 0) %>% #Keep only negative Dstats/Z-scores
  arrange(Z_score_trio) %>% #Arrange Z-scores from lowest to highest (we want the highest absolute value)
  distinct(Population1, Population3, .keep_all = TRUE) #Keep only the first instance of the combination, which will have the lowest Z-score due to the line above

#Repeat the process for any population of interest by changing the name of the population in the
#first filter line of the function

filtered_tableBBhp <- relevant_trios %>%
  filter(Population2 == "BBh" | Population3 == "BBh") %>%
  filter(Z_score_trio > 0) %>%
  arrange(desc(Z_score_trio)) %>%
  distinct(Population2, Population3, .keep_all = TRUE)

filtered_tableBBhm <- relevant_trios %>%
  filter(Population1 == "BBh" | Population3 == "BBh") %>%
  filter(Z_score_trio < 0) %>%
  arrange(Z_score_trio) %>%
  distinct(Population1, Population3, .keep_all = TRUE)

filtered_tableVMp <- relevant_trios %>%
  filter(Population2 == "VM" | Population3 == "VM") %>%
  filter(Z_score_trio > 0) %>%
  arrange(desc(Z_score_trio)) %>%
  distinct(Population2, Population3, .keep_all = TRUE)

filtered_tableVMm <- relevant_trios %>%
  filter(Population1 == "VM" | Population3 == "VM") %>%
  filter(Z_score_trio < 0) %>%
  arrange(Z_score_trio) %>%
  distinct(Population1, Population3, .keep_all = TRUE)

filtered_tablRBp <- relevant_trios %>%
  filter(Population2 == "RB" | Population3 == "RB") %>%
  filter(Z_score_trio > 0) %>%
  arrange(desc(Z_score_trio)) %>%
  distinct(Population2, Population3, .keep_all = TRUE)

filtered_tableRBm <- relevant_trios %>%
  filter(Population1 == "RB" | Population3 == "RB") %>%
  filter(Z_score_trio < 0) %>%
  arrange(Z_score_trio) %>%
  distinct(Population1, Population3, .keep_all = TRUE)

filtered_tableEXBBp <- result_trios %>%
  filter(Population2 == "EXBB" | Population3 == "EXBB") %>%
  filter(Z_score_trio > 0) %>%
  arrange(desc(Z_score_trio)) %>%
  distinct(Population2, Population3, .keep_all = TRUE)

filtered_tableERhBOTm <- result_trios %>%
  filter(Population1 == "ERhBOT" | Population3 == "ERhBOT") %>%
  filter(Z_score_trio < 0) %>%
  arrange(Z_score_trio) %>%
  distinct(Population1, Population3, .keep_all = TRUE)

filtered_tableERhBOTp <- result_trios %>%
  filter(Population2 == "ERhBOT" | Population3 == "ERhBOT") %>%
  filter(Z_score_trio > 0) %>%
  arrange(desc(Z_score_trio)) %>%
  distinct(Population2, Population3, .keep_all = TRUE)

filtered_tableERRBm <- result_trios %>%
  filter(Population1 == "ERRB" | Population3 == "ERRB") %>%
  filter(Z_score_trio < 0) %>%
  arrange(Z_score_trio) %>%
  distinct(Population1, Population3, .keep_all = TRUE)

filtered_tableERRBDp <- result_trios %>%
  filter(Population2 == "ERRB" | Population3 == "ERRB") %>%
  filter(Z_score_trio > 0) %>%
  arrange(desc(Z_score_trio)) %>%
  distinct(Population2, Population3, .keep_all = TRUE)

filtered_tableEMVMm <- result_trios %>%
  filter(Population1 == "EMVM" | Population3 == "EMVM") %>%
  filter(Z_score_trio < 0) %>%
  arrange(Z_score_trio) %>%
  distinct(Population1, Population3, .keep_all = TRUE)

filtered_tableEMVMp <- result_trios %>%
  filter(Population2 == "EMVM" | Population3 == "EMVM") %>%
  filter(Z_score_trio > 0) %>%
  arrange(Z_score_trio) %>%
  distinct(Population2, Population3, .keep_all = TRUE)


#Compare to some completely E.maxima and some completely E.radiata samples
filtered_tableEMJBp <- result_trios %>%
  filter(Population2 == "EMJB" | Population3 == "EMJB") %>%
  filter(Z_score_trio > 0) %>%
  arrange(desc(Z_score_trio)) %>%
  distinct(Population2, Population3, .keep_all = TRUE)

filtered_tableEMJBm <- result_trios %>%
  filter(Population1 == "EMJB" | Population3 == "EMJB") %>%
  filter(Z_score_trio < 0) %>%
  arrange(Z_score_trio) %>%
  distinct(Population1, Population3, .keep_all = TRUE)





