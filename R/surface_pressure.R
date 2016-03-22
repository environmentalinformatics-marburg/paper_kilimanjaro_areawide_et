### environment -----

## source functions
source("R/downloadECMWF.R")

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
rst_ref_utm <- raster("data/MYD09Q1.006/ndvi/NDVI_MYD09Q1.A2013001.sur_refl.tif")
rst_ref <- trim(projectRaster(rst_ref_utm, crs = "+init=epsg:4326"))

## digital elevation model
rst_dem <- raster("data/dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- trim(projectRaster(resample(rst_dem, rst_ref_utm), crs = "+init=epsg:4326"))

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
    
    ## resample to reference grid (extent, resolution)
    rst_res <- raster::resample(rst, rst_ref, filename = fls_res, 
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
    
    ## resample to reference grid (extent, resolution)
    rst_res <- raster::resample(rst, rst_ref, filename = fls_res, 
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
             
  if (file.exists(k)) {
    raster::raster(k)
  } else {  
                    
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
