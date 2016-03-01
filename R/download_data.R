### environment ----------------------------------------------------------------

## packages
library(Rsenal)
library(MODIS)

## working directory
Orcs::setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/", 
              path_ext = "kilimanjaro/evapotranspiration/")

## modis options
MODISoptions(localArcPath = "MODIS_ARC", 
             outDirPath = "MODIS_ARC/PROCESSED", 
             MODISserverOrder = c("LPDAAC", "LAADS"), 
             outProj = "+init=epsg:21037")

product <- "MYD09Q1"
collection <- getCollection(product, forceCheck = TRUE)
begin <- "2013001"; end <- "2015365"
tileH <- 21; tileV <- 9
job <- product


### myd09q1 --------------------------------------------------------------------

runGdal(product = product, collection = collection, begin = begin, end = end, 
        tileH = tileH, tileV = tileV, job = job)

# "110011000000" for MYD11A2
# "010000000000000001000" for MYD09GA