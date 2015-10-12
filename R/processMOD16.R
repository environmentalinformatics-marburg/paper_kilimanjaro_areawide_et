### environmental stuff

## workspace clearance
rm(list = ls(all = TRUE))

## working directory
# library(devtools)
# install_github("fdetsch/Orcs")
library(Orcs)
setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/", 
        path_ext = "kilimanjaro/evapotranspiration")

# packages
# install.packages("MODIS", repos="http://R-Forge.R-project.org")
lib <- c("MODIS", "doParallel", "Rsenal")
sapply(lib, function(...) require(..., character.only = TRUE))

library(doParallel)
library(MODIS)
library(Rsenal)

## parallelization
cl <- makeCluster(3)
registerDoParallel(cl)


### Data import

## 'MODIS' global settings
MODISoptions(localArcPath = paste0(getwd(), "/MODIS_ARC/"), 
             outDirPath = paste0(getwd(), "/MODIS_ARC/PROCESSED/"), 
             MODISserverOrder = c("LAADS","LPDAAC"), 
             outProj = "+init=epsg:21037", quiet = TRUE)

## geographic extent
rst_kili <- kiliAerial(rasterize = TRUE, minNumTiles = 20)

## evapotranspiration (et) data
sds <- c("ET_1km.tif$", "ET_QC_1km.tif$")
  
## crop
ls_crp <- foreach(pattern = sds, overwrite = rep(FALSE, 2)) %do% {                                      
    
  # available files
  fls <- list.files(options()$MODIS_outDirPath, pattern = pattern, 
                    full.names = TRUE, recursive = TRUE)  
  
  # stack, crop and store
  rst <- stack(fls)
  
  rst_crp <- do.call("stack", 
                     foreach(i = 1:nlayers(rst), 
                             .export = ls(envir = globalenv()), 
                             .packages = c("raster", "rgdal")) %dopar% {
                               
        # output filename                                    
        file_out <- paste0("mod16a2/crp/CRP_", names(rst[[i]]))
        
        # crop
        if (!file.exists(file_out) | overwrite) {
          crop(rst[[i]], extent(rst_kili), filename = file_out, 
               format = "GTiff", overwrite = TRUE)
          
        } else {
          raster(file_out)
        }
        
  })

  return(rst_crp)
}
  
## isotope plots
shp_plots <- readOGR(dsn = "/media/permanent/kilimanjaro/coordinates/coords", 
                     layer = "PlotPoles_ARC1960_mod_20140807_final", 
                     p4s = "+init=epsg:21037")

shp_plots_amp <- subset(shp_plots, PoleType == "AMP")
shp_plots_amp <- subset(shp_plots_amp, 
                        PlotID %in% c("sav5", "hom4", "fpd0", "flm1", 
                                      "foc6", "foc0", "fpo0", "fer0"))
