install.packages("diveRsity")
library(diveRsity)

divMeasures <- fastDivPart("CONVERT_ECKSA_MAR2020_noSmallPops_MOZasOnePop_genepop", outfile = "fastDivPart_calcs", fst = TRUE
                    ,bs_locus = TRUE, bs_pairwise = TRUE, boots = 10000, plot = FALSE, para = TRUE)


fstEck <- diffCalc("CONVERT_ECKSA_AUS_2021_BigPops", outfile = "diffCalcFinal_BigPops",
                       fst = TRUE, pairwise = TRUE, bs_locus = TRUE, bs_pairwise = TRUE,
                       boots = 10000, alpha = 0.05, para = TRUE)

install.packages("xlsx")
install.packages("sendplot")
install.packages("foreach")


??diffCalc

citation("diveRsity")
