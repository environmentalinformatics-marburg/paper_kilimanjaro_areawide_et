### environment -----

## source functions
source("R/downloadECMWF.R")
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


### download data -----

## first and last start date, last end date
st <- "2013-01-01"; nd <- "2015-12-01"; last <- "2015-12-31"

## temperature 
fls_ta <- downloadECMWF(st, nd, last, prm = "temperature", 
                        dsn = "data/ecmwf/temperature")

## geopotential
fls_gp <- downloadECMWF(st, nd, last, prm = "geopotential", 
                        dsn = "data/ecmwf/geopotential")


### pre-processing -----

## reference resolution and extent (250-m, southern kilimanjaro region)
rst_ref <- raster("data/MYD09Q1.006/ndvi/NDVI_MYD09Q1.A2013001.sur_refl.tif")

## digital elevation model
rst_dem <- raster("data/dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- resample(rst_dem, rst_ref)

## resample temperature files
lst_ta <- foreach(i = fls_ta, .packages = c("raster", "rgdal")) %dopar% {
  
  # target file
  fls_res <- paste0(dirname(i), Orcs::pureBasename(i, TRUE), ".tif")
  
  if (!file.exists(fls_res)) {
    # import data including metadata
    sgr <- rgdal::readGDAL(i, as.is = TRUE, p4s = "+init=epsg:4326")
    rst <- raster::stack(sgr)
    
    info <- suppressWarnings(rgdal::GDALinfo(i))
    
    ## apply offset and scale factor (value * scale + offset)
    scoff <- unique(attr(info, "ScaleOffset"))
    if (any(scoff != matrix(c(1, 0), nc = 2)))
      rst <- rst * scoff[, "scale"] + scoff[, "offset"]
    
    ## project to utm-37s
    rst_prj <- raster::projectRaster(rst, crs = raster::projection(rst_ref))
    
    ## resample to reference grid (extent, resolution)
    rst_res <- raster::resample(rst_prj, rst_ref, 
                                filename = fls_res, 
                                format = "GTiff", overwrite = TRUE)
    
  } else {
    rst_res <- raster::stack(fls_res)
  }

  ## stack daily data
  lapply(seq(1, raster::nlayers(rst_res), 18), function(j) {
    rst_res[[j:(j+17)]]
  })
}
lst_ta <- unlist(lst_ta)

## resample geopotential files
lst_gp <- foreach(i = fls_gp, .packages = c("raster", "rgdal")) %dopar% {
  
  # target file
  fls_res <- paste0(dirname(i), Orcs::pureBasename(i, TRUE), ".tif")
  
  if (!file.exists(fls_res)) {
    # import data including metadata
    sgr <- rgdal::readGDAL(i, as.is = TRUE, p4s = "+init=epsg:4326")
    rst <- raster::stack(sgr)
    
    info <- suppressWarnings(rgdal::GDALinfo(i))
    
    ## divide by gravity constant
    rst <- rst / 9.80665
    
    ## project to utm-37s
    rst_prj <- raster::projectRaster(rst, crs = raster::projection(rst_ref))
    
    ## resample to reference grid (extent, resolution)
    rst_res <- raster::resample(rst_prj, rst_ref, 
                                filename = fls_res, 
                                format = "GTiff", overwrite = TRUE)
    
  } else {
    rst_res <- raster::stack(fls_res)
  }
  
  ## stack daily data
  lapply(seq(1, raster::nlayers(rst_res), 18), function(j) {
    rst_res[[j:(j+17)]]
  })
}
lst_gp <- unlist(lst_gp)


### processing: calculate surface pressure -----

## pressure levels
p_levels <- c(seq(400, 750, 50), seq(775, 1000, 25))

## raster template
rst_tmp <- rst_ref[[1]]
rst_tmp[] <- NA

## target folder files
dates <- strftime(seq(as.Date(st), as.Date(last), "day"), "%Y%j")
fls_sp <- gsub("geopotential", "surface_pressure", fls_gp[1])
fls_sp <- gsub(".nc$", ".tif", fls_sp)
fls_sp <- lapply(dates, function(i) {
  gsub(substr(basename(fls_sp), 1, 15), i, fls_sp)
})

dir_sp <- dirname(fls_sp[[1]])
if (!dir.exists(dir_sp)) dir.create(dir_sp)

## loop over days
lst_sp <- foreach(i = lst_gp, j = lst_ta, k = fls_sp, 
                  .packages = "raster") %dopar% {
                    
  mat_gp <- raster::as.matrix(i)
  mat_ta <- raster::as.matrix(j)
  
  num_sp <- sapply(1:nrow(mat_gp), function(l) {
    z <- rst_dem[l]
    gp <- mat_gp[l, ]
    ta <- mat_ta[l, ]
    
    if (is.na(z) | any(is.na(gp)))
      return(NA)
    
    barometricFormula(z, gp, ta, p_levels)
  })
  
  rst_sp <- raster::setValues(rst_tmp, num_sp)
  raster::writeRaster(rst_sp, k, format = "GTiff", overwrite = TRUE)
}
