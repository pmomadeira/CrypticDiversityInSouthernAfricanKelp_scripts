#PLotting AFC with result matrices from Genetix
#After exploring some different packages and how they process data into and from PCA's
#I've started to understand the "problem" behind plotting an AFC from Genetix in R
#For my specific dataset, the information in the AFC is transformed in order to be
#reduced to 4 axis/components. This can be seen in Genetix AFC itself, where we have the option
#to plot the data acording to different combinations of PC's (1-X, 2-Y; 3-x, 2-Y, 2-X, 4Y, etc)
#Usually the two first components are the ones that better explain data variation, so K1 and K2
#in this specific case. This is the equivalent of have K1 serve as the x axis and K2 as the y axis
#for the plot. If needed be, I can just run the AFC to check the significance percentages of each axis
#to add to the axis label in R

#Below there are two ways to plot AFC results on R, one using ggplot, which is more flexible/easier to use and another one
#using plotly, that creates an interactive graph were we can check de individual coordinates of each scatter dot


#Set working directory (D:/ at home, F:/ at work)
setwd("D:/EckSA_paper/EckSA_data")

#Read the table with the data
#By default, this file is called flig.tab by genetix. We need to remove the first two lines so the file can be
#read into R.
#(First line contains the number of individuals and PCs; Second includes de file name of the original matrix used in genetix)
#AFCdataAll <- read.table("EckSA_genetix/SamplePops/EckSA_AUS_alleles2.txt") We do not need this data to plot the AFC in R
AFCdataInd <- read.table("EckSA_genetix/SamplePops/EckSA_AUS_ind2.txt")

#Add column names to the AFC table so we can select which variables will be used to plot the points and to color them.
colnames(AFCdataInd) <- c("Sample", "Population", "K1", "K2", "K3", "K4")

#To edit the plot we can use custom colors or symbols for the dots. In this case, I tested different approaches to have the best
#possible plot. The options were a) regular dots all around vs b) symbols matching the sampling map (square for E.radiata, circle for E. maxima, triangle for contact zones)
#Another possible option was to color the symbols according to putative taxa or cluster.
#The options were a) one colour per species, matching the k=2 for structure; b) colour each population matching their cluster from
#K=8 in structure; or c) colour populations per mtDNA haplotypes

#To do so, the following vectors were created
#Symbols vector
plotShapes <- rep(c(16,17,15,16,17,15), times = c(264,97,22,46,78,266)) #matches one shape (15,16 or 17) to each sample according to sample site
#The idea is to repeat one of the shapes a certain amount of times, to match the number of samples. So we get 16, x264, 17, x97 and so on.

#To get K=2 colours
maximaPops <- AFCdataInd$Population > 22 #creates a vector were samples were population id number is 22 or lower retrieve FALSE, 23 or higher retrieve TRUE
#pop 22 is the last maxima population

length(maximaPops[maximaPops == FALSE]) #Counts the nº of FALSES, meaning maxima populations
length(maximaPops[maximaPops == TRUE]) #Counts the number of TRUES, meaing radiata populations

k2colors <- rep(c("#0099E6", "#FF6600"), times = c(507, 266)) #Creates a vector of colours to pass to the samples.
#In this case, the first 507 samples are Maxima and the last 266 are radiata, so we create a vector of "#0099E6" x 507 ... "#FF6600" x 266
#For k=8
k8colors <- rep(c("#2242D3","#0099E6", "#2EBFD2", "#0099E6", "#2EBFD2", "#6F65FA", "#15B79D", "#0099E6","#3CFFF6", "#EBB35C", "#FF6600"),
                times = c(96, 96, 24, 72, 25, 48, 23, 46, 77, 112, 154))



#Having the data and all the formating vectors we can then plot the AFC
#Load ggplot2

library(ggplot2)

#Use the ggplot function to use the AFC data, determine that PC1/K1 is the x axis and PC2/K2 is the y axis and that the dots
#should be coloured according to our color selection. Then add the point with geom_point and use shape to give the intended
#shape to each dot/square/triange/etc. To have a clean plot, we can delete the background and grids, adding only a border and
#lines across the X = 0 and Y = 0 values.
ggplot (AFCdataInd, aes(x = K1, y = K2, color = k2colors)) + geom_point(shape = plotShapes) +
  theme_classic() +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 0.1) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.1) +
  theme (axis.line =  element_line(colour = "black"), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

#Plot for k=8 colours
ggplot (AFCdataInd, aes(x = K1, y = K2, color = k8colors)) + geom_point(shape = plotShapes) +
  theme_classic() +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 0.1) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.1) +
  theme (axis.line =  element_line(colour = "black"), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5))


#######################################################################################################################
#Interactive AFC ploting
#This is done using the plotly library. Check the 2D PCA Scatter Plot section on https://plotly.com/r/pca-visualization/
#for a tutorial of the code using data from the stats package


#Load AFC "coordinates" table (Coordinates as in, contains the values of each PC for each sample)
#By default, this file is called flig.tab by genetix. We need to remove the first two lines so the file can be
#read into R.
#(First line contains the number of individuals and PCs; Second includes de file name of the original matrix used in genetix)

AFCdataInd <- read.table("EckSA_genetix/SamplePops/EckSA_AUS_ind2.txt")
#Add column names to the AFC table so we can select which variables will be used to plot the points and to color them.
colnames(AFCdataInd) <- c("Sample", "Population", "K1", "K2", "K3", "K4")

#Load plotly library
library(plotly)

#General function to plot the AFC using plotly.
fig2 <- plotly::plot_ly(AFCdataInd, x = ~K1, y = ~K2, color = ~AFCdataInd$Population, colors = k2colors, type = 'scatter', mode = 'markers')%>%
  layout(
    legend=list(title=list(text='color')),
    plot_bgcolor='white',
    xaxis = list(
      title = "Axis 1",
      zerolinecolor = "black",
      zerolinewidth = 0.2,
      gridcolor='#ffff'),
    yaxis = list(
      title = "Axis 2",
      zerolinecolor = "black",
      zerolinewidth = 0.2,
      gridcolor='#ffff'))
fig2
