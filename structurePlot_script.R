#Ploting Structure results from _Q matrices
#Install package
devtools::install_github('ericarcher/strataG')
#Load package
library(strataG)

#Load the Q matrix you want to plot
K4_str <- read.table("results_job_T249_q", header = FALSE)
#structurePlot function will act up if it can't find the "orig.pop" column, so rename all columns accordingly
#In this example first column are samples, second column the population number/name, and the third-fifth are
#The cluster membership percentages
colnames(K4_str) <- c("ind", "orig.pop", "P1", "P2", "P3")


#Run the structurePlot function to plot the Q matrix resuls. horiz = FALSE makes it so the graph is ploted vertically
#Things like colors and legend.position can be ajusted
structurePlot(K4_str,pop.col = 2, prob.col = 3, legend.position = "bottom", 
              horiz = FALSE ,plot = TRUE)


