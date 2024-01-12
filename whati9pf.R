# whati9pf.R
#
# whati9pf_github folder
# Author: Andrew J. Wiebe, 8 Jan 2024
#
# Plot the flow system for Whati peninsula with a permafrost boundary
# Read in the results from the https://github.com/ajwiebe77/ICOP24 Octave code for permafrost (pf)
# Use a slightly modified version of the holzbecher_v10.R function (https://github.com/ajwiebe77/theYWVC; Wiebe and McKenzie, 2022) to line up the streamlines on either side of branch cuts, following the sequential contouring approach by Holzbecher (2018).
#
# Instructions:
#    * User needs to clean up the csv files from Octave prior to loading them here (remove the header lines and remove the first column, which is blank) *
#    User must specify the working directory employed by the setwd(...) command
#
# Notes:
#    Modifies code from app.R in the following repository: https://github.com/ajwiebe77/theYWVC (Wiebe and McKenzie, 2022)
#    Relies on a slightly modified version of the holzbecher_v10.R function from https://github.com/ajwiebe77/theYWVC (Wiebe and McKenzie, 2022)
#    Relies on mergecheck.R from https://github.com/ajwiebe77/theYWVC (Wiebe and McKenzie, 2022)
#    Relies on mergecheck2.R from https://github.com/ajwiebe77/theYWVC (Wiebe and McKenzie, 2022)
#    Relies on the long2UTM.R function (O'Brien, 2012; employed by https://github.com/ajwiebe77/theYWVC, Wiebe and McKenzie [2022])
#    The holzbecher_v10pf.R function relies on the mergecheck.R and mergecheck2.R functions from https://github.com/ajwiebe77/theYWVC (Wiebe and McKenzie, 2022)
#    The holzbecher_v10pf.R function does not succeed at completely removing the branch cut but does line up the streamlines on either side of it.
#    The trimAtBoundariesIMPERM2.R function is a modified version of the trimAtBoundaries.R function from https://github.com/ajwiebe77/theYWVC (Wiebe and McKenzie, 2022)
#    Convert the points for the head_pf.shp, flowlines_pf.shp, and capzones_pf.shp shapefiles to polylines by using an operation such as "Points to path" in GIS software, and using the group and order fields
#    (The points along the branch cut in the flowlines_pf.shp shapefile can be removed manually in GIS software prior to converting to polylines.)
#    The contour interval of the flowlines can be changed in GIS software by selecting series of points with levels with a desired spacing
#    The longitude of the site must be passed to the long2UTM(...) function
#    
# Known problems:
#    The holzbecher_v10pf.R function does not succeed at identifying and returning the capture zone boundary for this case (though the original version from https://github.com/ajwiebe77/theYWVC works better for other cases).
#    
# References:
#    
#    Holzbecher, E. 2018. Streamline visualization of potential flow with branch cuts, with applications to groundwater. J. Flow Vis. Image Process. 25(2), 119-144. doi:10.1615/JFlowVisImageProc.2018025918.
#    O'Brien, 2012. Response to "Determining UTM zone (to convert) from longitude/latitude." https://stackoverflow.com/questions/9186496/determining-utm-zone-to-convert-from-longitude-latitude. Cited 19 Dec 2022.
#    Wiebe, A.J., and McKenzie, J.M., 2022. An Open-Source Web Tool for Visualizing Estimates of Well Capture Zones Near Surface Water Features. Poster presentation at: AGU Fall Meeting 2022, Chicago, IL, USA, 12-16 Dec 2022. http://dx.doi.org/10.22541/essoar.167267811.10671930/v1.
#    
#    Code ideas from the websites listed in the code.
# 

library(sp)
library(raster) # for exporting shapefiles (points)
library(foreign) # https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile
library(sf) ## transforming coordinates, loading shapefiles

# FUNCTIONS -------------------------------------------

### Set the working directory here:

setwd("C:/RStudio_working/whati9pf_github")

# List function files
# https://www.r-bloggers.com/2013/01/storing-a-function-in-a-separate-file-in-r/
## https://stackoverflow.com/questions/40716589/shiny-app-with-own-function
source("holzbecher_v10pf.R", local = TRUE)
source("mergecheck.R", local = TRUE)
source("mergecheck2.R", local = TRUE)
source("trimAtBoundariesIMPERM2.R", local = TRUE)
source("long2UTM.R", local = TRUE)


### -------------------------------------------------------------------------
# Common parameters

x_origin <- 485891
y_origin <- 7001519
alpha <- 62.62429 * pi/180
beta <- pi + alpha/2
r_wedge <- 500
rotate1 <- -6.867509 # degrees
beta <- beta + rotate1 * pi/180 # **** 2 Aug 2023 try switching the sign ****
K = 0.015
b <- 9
dh_dl <- 1E-5
q0 <- abs(- K * b * dh_dl)


### -------------------------------------------------------------------------
# xy coordinates
numpts <- 1000

z_x <- seq(0,r_wedge,r_wedge/(numpts - 1))
z_y <- seq(0,r_wedge,r_wedge/(numpts - 1))



### -------------------------------------------------------------------------
### Potential
phi <- read.csv(file="phi_z1000x1000real_pf.csv",head=FALSE,sep=",")

phi <- as.matrix(phi)
phi <- t(phi) # transpose

Cu <- -527

for(i in 1:numpts){ # plot the streamlines
  for(ii in 1:numpts){
    if(! is.nan(phi[ii,i])){
      phi[ii,i] <- sqrt(2*(phi[ii,i] - Cu) / K) # NaN produced if phi[ii,i] < 0
    }
  }
}

phi_list <- contourLines(z_x, z_y, phi, levels=seq(min(phi[!is.nan(phi)]), max(phi[!is.nan(phi)]),by=0.00001))# find the capture zone boundary


### -------------------------------------------------------------------------
### Streamlines
psi <- read.csv(file="phi_z1000x1000imag_pf.csv",head=FALSE,sep=",")

psi <- as.matrix(psi)
psi <- t(psi) # transpose

# PF case
wells_temp <- rbind(c(486009.97, 7001611.41, 0.00059)) # holzbecher_v10 function cannot handle permafrost well because it is outside the allotted domain

numwells <- nrow(wells_temp)

wells_temp[,1] <- (wells_temp[,1] - x_origin)# / r_wedge # note: also need to account for rotation!
wells_temp[,2] <- (wells_temp[,2] - y_origin)# / r_wedge

rot_matrix <- rbind(c(cos(- rotate1 * pi/180), - sin(- rotate1 * pi / 180)), c(sin(- rotate1 * pi / 180), cos(- rotate1 * pi / 180)))
wells_temp2 <- matrix(0,numwells,3)
for(i in 1:numwells){
  wells_temp2[i,1] <- rot_matrix[1,1] * wells_temp[i,1] + rot_matrix[1,2] * wells_temp[i,2]
  wells_temp2[i,2] <- rot_matrix[2,1] * wells_temp[i,1] + rot_matrix[2,2] * wells_temp[i,2]
}

wells_temp2[,3] <- wells_temp[,3]

consts <- cbind(c(r_wedge, alpha, beta, K, b, q0))

print(wells_temp2)

xy_stg <- matrix(0.0,1,2)
xy_stg[1,1] <-  360.36 # manually enter the value from the xy_stg*.txt file
xy_stg[1,2] <- 223.72 # manually enter the value from the xy_stg*.txt file

print(xy_stg)

str_holz <- holzbecher_v10pf(psi, cbind(z_x,z_y), wells_temp2, consts, xy_stg) # works but does not find capture zone boundary
# spacing is output to the terminal as the third parameter listed after this function is called

psi_list_cz <- contourLines(z_x, z_y, psi, levels=psi[(floor(numpts * xy_stg[1,1] / r_wedge)),(floor(numpts * xy_stg[1,2] / r_wedge))]) # find the capture zone boundary

psi_check <- contourLines(z_x,z_y,psi,30)

### -------------------------------------------------------------------------
# equipotential contours

phi_list <- trimAtBoundariesIMPERM2(phi_list, 2, 1, alpha, r_wedge, TRUE, TRUE)



### -------------------------------------------------------------------------
# streamline contours

psi_list <- str_holz[[1]]

# if(length(str_holz[[2]]) > 0){ # add this 3 Jan 2023
  # capzones <- str_holz[[2]] # empty list - try a different approach - psi_list_cz
# }

psi_list <- trimAtBoundariesIMPERM2(psi_list, 2, 1, alpha, r_wedge, TRUE, TRUE)
psi_list_cz <- trimAtBoundariesIMPERM2(psi_list_cz, 2,1, alpha, r_wedge, TRUE, TRUE)

ptA <- c(r_wedge * cos(alpha), r_wedge * sin(alpha))
ptB <- c(0, 0)
ptC <- c(r_wedge, 0)
boundary <- rbind(ptA, ptB, ptC)

xlimits <- c(0,r_wedge)
ylimits <- c(0,r_wedge)

plot(boundary[,1], boundary[,2], lwd = 3, type = "l", col = "black", xlim=xlimits, ylim=ylimits, xlab="x", ylab="y")

# Streamlines

for(i in 1:length(psi_list)){ # plot the streamlines
lines(psi_list[[i]]$x, psi_list[[i]]$y, col="red") # the branch cut has been successfully corrected for all streamlines within the wedge aquifer, except that the cut line is still present
}

# Capture zones

for(i in 1:length(psi_list_cz)){ # plot the streamlines
lines(psi_list_cz[[i]]$x, psi_list_cz[[i]]$y, col="black") # the branch cut has been successfully corrected for all streamlines within the wedge aquifer, except that the cut line is still present
}

# head contours

for(i in 1:length(phi_list)){ # plot the streamlines
  lines(phi_list[[i]]$x, phi_list[[i]]$y, col="blue") # the branch cut has been successfully corrected for all streamlines within the wedge aquifer, except that the cut line is still present
}

points(xy_stg[,1], xy_stg[,2])


### -------------------------------------------------------------------------


# ### rotate all points ---------------------

rot_matrix <- rbind(c(cos(rotate1 * pi/180), - sin(rotate1 * pi / 180)), c(sin(rotate1 * pi / 180), cos(rotate1 * pi / 180)))

xy_stg2 <- xy_stg

for(i in 1:nrow(xy_stg2)){
  xy_stg2[i,1] <- rot_matrix[1,1] * xy_stg[i,1] + rot_matrix[1,2] * xy_stg[i,2]
  xy_stg2[i,2] <- rot_matrix[2,1] * xy_stg[i,1] + rot_matrix[2,2] * xy_stg[i,2]
}

psi_list2 <- psi_list

for(i in 1:length(psi_list2)){
  psi_list2[[i]]$x <- rot_matrix[1,1] * psi_list[[i]]$x + rot_matrix[1,2] * psi_list[[i]]$y
  psi_list2[[i]]$y <- rot_matrix[2,1] * psi_list[[i]]$x + rot_matrix[2,2] * psi_list[[i]]$y
}


capzones2 <- psi_list_cz

for(i in 1:length(capzones2)){
  capzones2[[i]]$x <- rot_matrix[1,1] * psi_list_cz[[i]]$x + rot_matrix[1,2] * psi_list_cz[[i]]$y # 1 Aug 2023
  capzones2[[i]]$y <- rot_matrix[2,1] * psi_list_cz[[i]]$x + rot_matrix[2,2] * psi_list_cz[[i]]$y # 1 Aug 2023
}


phi_list2 <- phi_list

for(i in 1:length(phi_list2)){
  phi_list2[[i]]$x <- rot_matrix[1,1] * phi_list[[i]]$x + rot_matrix[1,2] * phi_list[[i]]$y
  phi_list2[[i]]$y <- rot_matrix[2,1] * phi_list[[i]]$x + rot_matrix[2,2] * phi_list[[i]]$y
}

### Translate all points

xy_stg2[,1] <- xy_stg2[,1] + x_origin
xy_stg2[,2] <- xy_stg2[,2] + y_origin

for(i in 1:length(psi_list2)){
  psi_list2[[i]]$x <- psi_list2[[i]]$x + x_origin
  psi_list2[[i]]$y <- psi_list2[[i]]$y + y_origin
}


for(i in 1:length(capzones2)){
  capzones2[[i]]$x <- capzones2[[i]]$x + x_origin # 1 Aug 2023
  capzones2[[i]]$y <- capzones2[[i]]$y + y_origin # 1 Aug 2023
}


for(i in 1:length(phi_list2)){
  phi_list2[[i]]$x <- phi_list2[[i]]$x + x_origin
  phi_list2[[i]]$y <- phi_list2[[i]]$y + y_origin
}


# -------------------------------------------------------------------------------
# Export shapfiles

z1 <- long2UTM(-117.27)
proj <- CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84"))

### export head contour lines

counter <- 1
xyphi <- matrix(data=0.0,numpts*numpts,5)

for(i in 1:length(phi_list2)){
  for(ii in 1:length(phi_list2[[i]]$x)){
    if(is.na(phi_list2[[i]]$x[ii]) == FALSE){ # points outside aquifer domain will have NA coordinates due to trimAtBoundaries() function
      xyphi[counter,1] <- phi_list2[[i]]$x[ii]
      xyphi[counter,2] <- phi_list2[[i]]$y[ii]
      xyphi[counter,3] <- i # group number
      xyphi[counter,4] <- ii # order number
      xyphi[counter,5] <- phi_list2[[i]]$level # 1 Jan 2024
      counter <- counter + 1
    }
  }
}

xyphi <- xyphi[c(1:counter-1), c(1,2,3,4,5),drop=FALSE] # remove zero rows

equipotpts <- SpatialPoints(xyphi, proj4string = proj)
shapefile(equipotpts, filename = "head_pf.shp", overwrite = TRUE) # use the code below in the download handler
### add attribute (https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile)
dbfdata <- read.dbf("head_pf.dbf", as.is = TRUE)
# ## add new attribute data (just the numbers 1 to the number of objects)
dbfdata$group <- xyphi[,3] ### creates a new field called "group"
dbfdata$order <- xyphi[,4] ### creates a new field called "order"
dbfdata$level <- xyphi[,5] ### creates a new field called "level" # 1 Jan 2024
write.dbf(dbfdata, "head_pf.dbf") ## overwrite the file with this new copy

### export flowlines

counter <- 1
xypsi <- matrix(data=0.0,numpts*numpts,5)

for(i in 1:length(psi_list2)){
  for(ii in 1:length(psi_list2[[i]]$x)){
    if(is.na(psi_list2[[i]]$x[ii]) == FALSE){ # points outside aquifer domain will have NA coordinates due to trimAtBoundaries() function
      xypsi[counter,1] <- psi_list2[[i]]$x[ii]
      xypsi[counter,2] <- psi_list2[[i]]$y[ii]
      xypsi[counter,3] <- i # group number
      xypsi[counter,4] <- ii # order number
      xypsi[counter,5] <- psi_list2[[i]]$level # level # 1 Jan 2024
      counter <- counter + 1
    }
  }
}

xypsi <- xypsi[c(1:counter-1), c(1,2,3,4,5),drop=FALSE] # remove zero rows

flowlinepts <- SpatialPoints(xypsi, proj4string = proj)
shapefile(flowlinepts, filename = "flowlines_pf.shp", overwrite = TRUE) # use the code below in the download handler
## add attribute (https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile)
dbfdata <- read.dbf("flowlines_pf.dbf", as.is = TRUE)
# ## add new attribute data (just the numbers 1 to the number of objects)
dbfdata$group <- xypsi[,3] ### creates a new field called "group"
dbfdata$order <- xypsi[,4] ### creates a new field called "order"
dbfdata$level <- xypsi[,5] ### creates a new field called "level" # 1 Jan 2024
write.dbf(dbfdata, "flowlines_pf.dbf") ## overwrite the file with this new copy

### export capture zone boundaries

counter <- 1
xycap <- matrix(data=0.0,numpts*numpts,4)

for(i in 1:length(capzones2)){
  for(ii in 1:length(capzones2[[i]]$x)){ # 1 Aug 2023
    if(is.na(capzones2[[i]]$x[ii]) == FALSE){ # points outside aquifer domain will have NA coordinates due to trimAtBoundaries() function
      xycap[counter,1] <- capzones2[[i]]$x[ii] # 1 Aug 2023
      xycap[counter,2] <- capzones2[[i]]$y[ii] # 1 Aug 2023
      xycap[counter,3] <- i # group number
      xycap[counter,4] <- ii # order number
      counter <- counter + 1
    }
  }
}

xycap <- xycap[c(1:counter-1), c(1,2,3,4),drop=FALSE] # remove zero rows

cappts <- SpatialPoints(xycap, proj4string = proj)
shapefile(cappts, filename = "capzones_pf.shp", overwrite = TRUE) # use the code below in the download handler
## add attribute (https://gis.stackexchange.com/questions/6839/adding-attribute-data-to-shapefile)
dbfdata <- read.dbf("capzones_pf.dbf", as.is = TRUE)
# ## add new attribute data (just the numbers 1 to the number of objects)
dbfdata$group <- xycap[,3] ### creates a new field called "group"
dbfdata$order <- xycap[,4] ### creates a new field called "order"
write.dbf(dbfdata, "capzones_pf.dbf") ## overwrite the file with this new copy


### export stagnation points

stgpts <- SpatialPoints(xy_stg2, proj4string = proj)
str(stgpts)
projLatLong <- CRS(paste(sep="","+proj=utm ", "+zone=", z1, " +datum=WGS84")) # EPSG:3857
stgptsLatLong <- SpatialPoints(xy_stg2, proj4string=CRS("+proj=longlat +ellips=WGS84"))
shapefile(stgpts, filename = "stagnationpts_pf.shp", overwrite = TRUE)

