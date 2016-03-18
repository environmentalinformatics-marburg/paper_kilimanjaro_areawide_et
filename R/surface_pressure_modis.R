### environment -----

## clear workspace
rm(list = ls(all = TRUE))

## source functions
source("R/barometricFormula.R")

## set working directory
Orcs::setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/",
              path_ext = "kilimanjaro/evapotranspiration/")

## load packages
library(raster)
library(doParallel)

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


### pre-processing -----

## reference raster and digital elevation model
rst_ref <- raster("data/MYD09Q1.006/ndvi/NDVI_MYD09Q1.A2013001.sur_refl.tif")

rst_dem <- raster("data/dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- resample(rst_dem, rst_ref)

rst_ref <- trim(projectRaster(rst_ref, crs = "+init=epsg:4326"))
rst_dem <- trim(projectRaster(rst_dem, crs = "+init=epsg:4326"))

## list and import temperature files
fls_ta <- list.files("data/MOD07_L2.006/Retrieved_Temperature_Profile", 
                     full.names = TRUE, pattern = "Profile.tif$")

lst_ta <- foreach(i = fls_ta, .packages = "raster") %dopar% raster::stack(i)

## resample geopotential files
fls_gp <- list.files("data/MOD07_L2.006/Retrieved_Height_Profile", 
                     full.names = TRUE, pattern = "Profile.tif$")

lst_gp <- foreach(i = fls_gp, .packages = "raster") %dopar% raster::stack(i)


### processing: calculate surface pressure -----

## pressure levels
p <- c(5, 10, 20, 30, 50, 70, 100, 150, 200, 250, 
       300, 400, 500, 620, 700, 780, 850, 920, 950, 1000)

## raster template
rst_tmp <- setValues(rst_ref, values = NA)

## target folder files
fls_sp <- gsub("Retrieved_Height_Profile", "Surface_Pressure", fls_gp)
if (!dir.exists(unique(dirname(fls_sp)))) dir.create(unique(dirname(fls_sp)))

## loop over scenes
lst_sp <- foreach(i = lst_gp, j = lst_ta, k = as.list(fls_sp), 
                  .packages = "raster") %dopar% {
  
  if (!(all(is.na(raster::maxValue(i))) | all(is.na(raster::maxValue(j))))) {
    if (file.exists(k)) {
      raster::raster(k)
    } else {  
      
      mat_gp <- raster::as.matrix(i)
      mat_ta <- raster::as.matrix(j)
      
      num_sp <- sapply(1:nrow(mat_gp), function(l) {
        z <- rst_dem[l]
        gp <- mat_gp[l, ]
        ta <- mat_ta[l, ]
        
        if (is.na(z) | all(is.na(gp)))
          return(NA)
        
        barometricFormula(z, gp, ta, p)
      })
      
      rst_sp <- raster::setValues(rst_tmp, num_sp)
      raster::writeRaster(rst_sp, k, format = "GTiff", overwrite = TRUE)
    }
  }
}
                    
rst_sp <- stack(lst_sp)
                    
## aggregate 8-day intervals
dts_sp <- sapply(fls_sp, function(i) substr(basename(i), 1, 7))
dts_sp <- as.Date(dts_sp, "%Y%j")

lst_agg <- foreach(i = 2013:2015, .combine = "c", .packages = "raster") %dopar% {
  
  # output folder and files
  id <- grep(paste0("^", i), dts_sp)
  
  timestamps <- strftime(as.Date(levels(cut(dts_sp[id], "8 days"))), "%Y%j")
  fls_agg <- fls_sp[sapply(timestamps, function(j) grep(j, unlist(fls_sp)))]
  dir_agg <- paste0(unique(dirname(unlist(fls_agg))), "/agg")
  if (!dir.exists(dir_agg)) dir.create(dir_agg)
  fls_agg <- basename(unlist(fls_agg))
  
  # aggregate
  indices <- as.numeric(cut(dts_sp[id], "8 days"))
  
  raster::stackApply(rst_sp[[id]], indices = indices, fun = "mean", 
                     filename = paste0(dir_agg, "/AGG"), bylayer = TRUE, 
                     suffix = fls_agg, format = "GTiff", overwrite = TRUE)
}

rst_agg <- stack(lst_agg)

## deregister parallel backend
stopCluster(cl)
