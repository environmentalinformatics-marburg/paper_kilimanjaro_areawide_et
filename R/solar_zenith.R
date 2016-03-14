### environment -----

## set working directory
Orcs::setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/",
              path_ext = "kilimanjaro/evapotranspiration/")

## load packages
library(raster)
library(doParallel)

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


### pre-processing: auxiliary data -----

## reference resolution and extent (250-m, southern kilimanjaro region)
rst_ref <- raster("data/MYD09Q1.006/ndvi/NDVI_MYD09Q1.A2013001.sur_refl.tif")

## digital elevation model
rst_dem <- raster("data/dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- resample(rst_dem, rst_ref)


### processing: solar zenith -----

## loop over products
lst_theta <- foreach(product = c("MOD05_L2.006", "MYD05_L2.006")) %do% {
  
  ## status message
  cat("Commencing with", product, "...\n")
  
  ## list available files
  fls_theta <- list.files(paste0("data/", product, "/Solar_Zenith"), 
                          full.names = TRUE, pattern = "Solar_Zenith.tif$")
  
  ## loop over files
  lst_theta <- foreach(i = fls_theta, .packages = "raster") %dopar% {
    
    # target file
    dir_res <- paste0(dirname(i), "/res")
    if (!dir.exists(dir_res)) dir.create(dir_res)
    
    fls_res <- paste0(dir_res, Orcs::pureBasename(i, TRUE), ".tif")
    
    if (!file.exists(fls_res)) {
      # import file
      rst <- raster(i)
      
      # project to utm-37s
      rst_prj <- raster::projectRaster(rst, crs = raster::projection(rst_ref))
      
      # resample to reference grid (extent, resolution)
      raster::resample(rst_prj, rst_ref, 
                       filename = fls_res, format = "GTiff", overwrite = TRUE)
      
    } else {
      raster::raster(fls_res)
    }
  }
  
  ## stack and return scenes
  raster::stack(lst_theta)
}

lst_ta <- unlist(lst_ta)
