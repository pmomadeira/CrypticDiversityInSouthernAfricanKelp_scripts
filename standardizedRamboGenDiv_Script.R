source("https://raw.githubusercontent.com/jorgeassis/rambo/master/Script.R")

#Set working directory
setwd("F:/EckSA_paper/EckSA_data")

#Load libraries
library(adegenet)
library(pegas)

#Load the allele matrix
EckSA_gen <- read.genetix("EckSA_genetix/BigPops/Eck_AUS_final_matrix2021_BigPops.gtx")

#Define a clustering vector for the populations to be tested
clustering.vector <- c(1:25) #in the case of BigPops, there are 25 different populations

#Define the variables
missing.data <- 000
replace <- FALSE
ncode <- 3
resample.number.auto <- TRUE
resample.number <- 11
discard.pops <- NULL
number.interactions <- 1000
alpha.number <- 0.05
savefile <- FALSE
save.filname <- "normalizedBigPops"
#Run the normalized genetic diversity algorithm
BigPopRamb <- Rambo("EckSA_genetix/BigPops/Eck_AUS_final_matrix2021_BigPops.gtx", missing.data, ncode, replace, resample.number.auto, resample.number, discard.pops,
                   number.interactions, alpha.number, clustering.vector, savefile)

EckRamb <- Rambo("EckSA_genetix/BigPops/Eck_AUS_final_matrix2021_BigPops.gtx", 000, 3, replace = FALSE, resample.number.auto = FALSE, 
      resample.number = 20, number.iteractions = 10000, alfa.test = 0.05, clustering.vector, 
      savefile = FALSE)

write.table(EckRamb,"EckSA_normalized_diversity.txt" , sep = ";", row.names = FALSE, col.names = TRUE, na = "NA", dec = ".")



