#Using adegenet for genetic diversity analysis + DAPC
library(adegenet)
library(pegas)

#Set working directory
setwd("F:/EckSA_paper/EckSA_data")


#Import file containing the allelic information
EckSAg <- read.genetix("EckSA_AUS.gtx")

EckSA <- read.structure("CONVERT_ECKSA_AUS_2021.str", onerowperind = FALSE, 
                        n.ind = 773 , n.loc = 8, col.lab = 1, col.pop = 2, NA.char = "999",
                        row.marknames = 1)

EckSA_noNAs <- tab(EckSA, freq = FALSE, NA.method = "mean") #replaces NAs with the mean Allele freq value

#EckSA_noNAs_test <- tab(EckSA, freq = TRUE, NA.method = "mean")
#EckSA_scale <- scaleGen(EckSA, center = FALSE, scale = TRUE, NA.method = "mean")

#Statistics
nInd(EckSA) #nº of individuals = 773
nLoc(EckSA) #nº of loci = 8
nAll(EckSA) #nº of alleles per loci
#EC03   EC07   EC10   EC12 Erad01 Erad04 Erad05 Erad09 
#20     29     32     37     19     19      9     16 

#All Summary Statistics
EckSumm <- summary(EckSA)

par(mfrow=c(2,2))
plot(EckSumm$n.by.pop, EckSumm$pop.n.all, xlab="Sample Size", ylab="Number of Alleles",
     main="Allele number and Sample Size", type="n")
text(EckSumm$n.by.pop, EckSumm$pop.n.all, lab=names(EckSumm$n.by.pop), col=funky(30))

barplot(EckSumm$loc.n.all, ylab="Number of Alleles", main = "Number of alleles per locus")

barplot(EckSumm$Hexp-EckSumm$Hobs, main="Heterozigosity: expected-observed", ylab="Hexp - Hobs")

barplot(EckSumm$n.by.pop, main = "Sample Size per Population", ylab="Number of genotypes", las = 3)

Eckpop <- seppop(EckSA)

#Testing if mean Hexpected is significantly different from Hobserved

bartlett.test(list(EckSumm$Hexp, EckSumm$Hobs)) #bartlett test

t.test(EckSumm$Hexp,EckSumm$Hobs,pair=T,var.equal=TRUE,alter="greater") #t test

EckHWT <- hw.test(EckSA, B=0)

#FST
library("hierfstat")
fstat(EckSA)

EckFST <- pairwise.fst(EckSA)
#DAPC test
set.seed(13) #for reproducibility use it before each find.cluster

groups8 <- find.clusters(EckSA, n.clust = NULL, method = "kmeans")
groups11 <- find.clusters(EckSA, n.clust = NULL, method = "kmeans")
groups16 <- find.clusters(EckSA, n.clust = NULL, method = "kmeans")
groups16 <- groups15

#names(groups8) lists objects that are part of groups
#min(groups15$Kstat) shows mininmal k-mean value

#List to which group each sample belongs and print it to a table

#groups8$grp #lists the samples and the cluster to which they belong

#turns the lists into data frames so we can write them out later
SvG8 <- as.data.frame(groups8$grp)
SvG11 <- as.data.frame(groups11$grp)
SvG16 <- as.data.frame(groups16$grp)

#Writing group attribution tables
write.csv2(SvG8, "SampleGroupsMatch_k=8.csv")
write.csv2(SvG11, "SampleGroupsMatch_k=11.csv")
write.csv2(SvG16, "SampleGroupsMatch_k=16.csv")



CT8 <- table(pop(EckSA), groups8$grp)  #making a table with poplulations x groups
CT11 <- table(pop(EckSA), groups11$grp)
CT16 <- table(pop(EckSA), groups16$grp)

#Writing the tables for posterior checking
write.csv2(CT8, "SpeciesvClusters_k=8.csv")
write.csv2(CT11, "SpeciesvClusters_k=11.csv")
write.csv2(CT16, "SpeciesvClusters_k=16.csv")

dapc8<- dapc(EckSA, groups8$grp)
dapc11 <- dapc(EckSA, groups11$grp)
dapc16 <- dapc(EckSA, groups16$grp)

scatter.dapc(dapc8)
scatter.dapc(dapc11)
scatter.dapc(dapc16)


#dapcColors <- c("blue", "blueviolet", "aquamarine", "burlywood4", "chartreuse", "chocolate1",
#                "coral3", "cyan3", "darkgoldenrod", "darkgreen", "darkred", "mediumpurple",
#                "navy", "palevioletred", "slategray3")


#Create a series of color vectors to pass through each clusters, one for each dapc (8, 11 and 16 clusters)
dpCl8 <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#9C9C9C')
dpCl11 <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#9C9C9C','#1DF0D7',
            '#ED1769','#A618F2')
dpCl16 <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#9C9C9C','#1DF0D7',
            '#ED1769','#A618F2','#11FF00', '#1D60A3','#FFE100', '#26FF00', '#2EC76E')
dp3col <- c('#00FF36','#FF00FF','#FF00FF','#FF00FF','#FF00FF','#FF00FF','#FF00FF','#9400D3','#FF00FF',
            '#FF00FF','#FF00FF', '#9400D3','#00FF36','#9400D3','#9400D3','#FF00FF' )



#Attribute a specific shape to each cluster. In this case I chose to give a shape to the putative species
#clusters. 17 = De Hoop Hybrids; 18 = E. maxima clusters, 20 = E. radiata clusters
shapes8 <- c(18,20,17,18,18,18,20,20) #shape per putative cluster species
shapes11 <-c(18,17,18,18,20,20,18,20,18,18,18)
shapes16 <- c(17,18,18,18,18,18,18,20,18,18,18,20,17,20,20,18)

#Draw the plots
scatter(dapc8, scree.da = FALSE, bg = "white", pch = shapes8, cell = 0, cstar = 0, col = dpCl8,
        solid = 0.4, cex = 1.5, clab = 0, leg = TRUE, txt.leg = paste("Cluster", 1:8), 
        posi.leg = "bottomright")

scatter(dapc11, scree.da = FALSE, bg = "white", pch = shapes11, cell = 0, cstar = 0, col = dpCl11,
        solid = 0.4, cex = 1.5, clab = 0, leg = TRUE, txt.leg = paste("Cluster", 1:11), 
        posi.leg = "bottomright")

scatter(dapc16, scree.da = FALSE, bg = "white", pch = shapes16, cell = 0, cstar = 0, col = dpCl16,
        solid = 0.4, cex = 1.5, clab = 0, leg = TRUE, txt.leg = paste("Cluster", 1:16), 

                posi.leg = "bottomright")
scatter(dapc16, scree.da = FALSE, bg = "white", pch = 20, cell = 0, cstar = 0, col = dp3col,
        solid = 0.4, cex = 1.5, clab = 0)

#Draw structure-like barplots
compoplot(dapc8, posi = "bottomright", txt.leg=paste("Cluster", 1:8), lab = "",
          ncol=1, xlab = "individuals", col=dpCl8)

compoplot(dapc11, posi = "bottomright", txt.leg=paste("Cluster", 1:11), lab = "",
          ncol=1, xlab = "individuals", col=dpCl11)
compoplot(dapc16, posi = "bottomright", txt.leg=paste("Cluster", 1:16), lab = "",
          ncol=1, xlab = "individuals", col=funky(16))

#Try crossvalidation

EckMat <- as.matrix(tab(EckSA, freq = FALSE, NA.method = "mean"))
grp <- pop(EckSA)

xval <- xvalDapc(EckMat, grp, n.pca.max = 150, training.set = 0.9, result = "groupMean",
                 center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 30, xval.plot = TRUE)

system.time(xval2 <- xvalDapc(EckMat, grp, n.pca.max = 150, training.set = 0.9, result = "groupMean",
                     center = TRUE, scale = FALSE, n.pca = 40:100, n.rep = 1000, xval.plot = TRUE,
                     parallel = "snow", ncpus = 4L))

~#test optimal PC's to keep. Don't make a great difference overall in the graphs, slighty changes the distances
temp8 <- optim.a.score(dapc8)
temp11 <- optim.a.score(dapc11)
temp16 <- optim.a.score(dapc16)
