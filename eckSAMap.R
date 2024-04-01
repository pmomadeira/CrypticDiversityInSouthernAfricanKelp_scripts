#The objective of this script is to map smapling points to a geographical map
#It also has the information necessary to plot STRUCTURE pie charts, to observe the percentages of
#genetic clusters in each population

#Load coordinates for sampling sites
maximaCoord <- read.table("maximaFCoordinates.txt", header = TRUE)
radiataCoord <- read.table("radiataFcoordinates.txt", header = TRUE)
contactCoord <- read.table("mapCoordinates_contactPoints.txt", header = TRUE)

#Load map raster and packages necessary for plotting
library(ggplot2)
library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(raster)
library(maps)

#Example using previously prepared raster files
#Load the raster you want
coastlines <- shapefile("../../Work_cloud/GitHub/courseMarineEcologicalModelling/Sessions [P]/Data/vectorShapefiles/coastlines/ne_10m_coastline.shp")
#Define geographic limitis with extent
maximaExtent <- extent(-15, 40, -40, -15)
#Crop the raster to fit the geographic limits you want to look at
maximaArea <- crop(coastlines, maximaExtent)

#Plot a simple map to check if you croped the region you wanted. Ajust extent values accordingly to you needs
plot(maximaArea)

#The same thing can be done by using the maps package

#world <- ne_countries(scale = "medium")
#maximaExtent <- extent(-15, 40, -40, -15)
#maximaArea <- crop(world, maximaExtent)
#plot(maximaArea)


#Plot the map while using ggplot
#In this case I loaded 3 different coordinates tables so I could use different symbols/colors for each species
#There might be a better or more elegant way to do it, but this works
ggplot() + 
  geom_polygon(data = maximaArea, aes(x = long, y = lat, group = group), fill = "white", colour = "black") +   #this line is to plot the map itself
  geom_point(data = maximaCoord, aes(x = Longitude, y = Latitude), color = "black", shape = 19, size = 3) +    #this line is to plot the points for maxima/species 1
  geom_text(data = maximaCoord, aes(x = Longitude, y = Latitude, label = CODE), nudge_x = 1, nudge_y = 0.25) + #This line adds the labels to maxima points
  geom_point(data = contactCoord, aes(x = Longitude, y = Latitude), color = "black", shape = 17, size = 3) +   #this line plots contact zones points/species 2
  geom_text(data = contactCoord,aes(x = Longitude, y = Latitude, label = CODE), nudge_x = 1, nudge_y = 0.25) + #this line adds the lables to contact points
  geom_point(data = radiataCoord, aes(x = Longitude, y = Latitude), color = "black", shape = 15, size = 3) +   #this line is to plot points for radiata/species 3
  geom_text(data = radiataCoord,aes(x = Longitude, y = Latitude, label = CODE), nudge_x = 1, nudge_y = 0.25) + #this line adds lables to radiata points
  scale_y_continuous(breaks = seq(-90,90, by=20)) +                                                            #this line defines the scale for longitude (20 by 20, meaning "0-20-40" etc)
  scale_x_continuous(breaks = seq(-180,180,by=20)) +                                                           #this line defines the scale for latitude (see line before)
  coord_map("mercator", xlim = c(5,41), ylim = c(-35,-21))+                                                    #this line defines the latitude and longitude limits of the map to be drawn. Points outside of that limit will not show
  theme_bw()


#Adding Structure results to the maps in the form of pie charts
#Load the STRUCTURE/ADMIXTURE tables
#In this situation I've used the output from the CLUMPAK server for my selected best deltaK
#The ClumppPopFile is a table with the percentage of membership for each cluster/K of each population 
#Load the ClumppPopFiles
maximaSTR <- read.table("C:/Users/Pedro Madeira/Desktop/EckSA_data/EckSA_str/maximaHyb/maximaHyb_CLUMPAK/1630945777/K=4/MajorCluster/CLUMPP.files/ClumppPopFile")
radiataSTR <- read.table("C:/Users/Pedro Madeira/Desktop/EckSA_data/EckSA_str/radiataHyb/radiataHyb_CLUMPAK/1631029663/K=3/MajorCluster/CLUMPP.files/ClumppPopFile")

#Rename the columns of each table to use in the graphics afterwards. 
#Having column names is important to select what will be plotted
#In this case, the first columns refer to each population, the next to the K membership percentages and the last one is population/sampling size.
#Having the number of samples allows us to draw pie charts with different sizes depending on sampling effort
colnames(maximaSTR) <- c("Pop", "K1", "K2", "K3", "K4", "Size")
colnames(radiataSTR) <- c("Pop", "K1", "K2", "K3", "Size")

#Initialize packages nedded to plot the maps
#Notice that the packages were already loaded before, plotrix being the only "new" package loaded
#library(ggplot2)
#library(RColorBrewer)
#library(rnaturalearth)
#library(rnaturalearthdata)
#library(rnaturalearthhires)
#library(raster)
#(maps)
library(plotrix)

#Plot the map of the study region (check earlier lines in this code to see how to do that)
#In this case I'm using the object maximaArea which is a raster previously cropped to fit my study area
plot(maximaArea)
#This loop will use the coordinates of you sampling sites, match it to the sites in your STRUCTURE/ADMIXTURE/Clumpp file,
#and the plot the pie charts on those coordinates. It's important that the Population order in the coordinates file
#matches the order of the populations in the Clumpp file

for (x in 1:nrow(maximaCoord)) #For each row of coordinates 
  { floating.pie(maximaCoord$Longitude[x],maximaCoord$Latitude[x], #Plot the pies in X/Y coordinates
                                      c(maximaSTR$K1[x],maximaSTR$K2[x],maximaSTR$K3[x],maximaSTR$K4[x]),radius=maximaSTR$Size[x]/8, #Using K's from 1 to 4 (1 to X, X being whatever your max K is)
                                      col=c("red","blue","orange","green")) } #And give the K's this colors (Order of colors matches K order, so K1 = red and so on)


#Repeat the process for each data set you want to plot
plot(maximaArea)

for (x in 1:nrow(radiataCoord)) { floating.pie(radiataCoord$Longitude[x],radiataCoord$Latitude[x], #You can modify your loop to reflect this
                                              c(radiataSTR$K1[x],radiataSTR$K2[x],radiataSTR$K3[x]),radius=radiataSTR$Size[x]/8,
                                              col=c("red","blue","orange")) }

#PS - note that the graphs coming out of this script are very rudimental. I've proceeded to edit them heavily
#On Affinity Designer (Inkscape and Illustrator will do the same) in order to get the final result
#There might be some ways to improve this output but I didn't have the time to find it yet
#To help anyone who uses this, the pie charts are printed in order, so the last pie chart (the one on top) is
#for the last population. So you can then resize/move them on you software of choice and not lose track of which pie chart
#Comes from which population




