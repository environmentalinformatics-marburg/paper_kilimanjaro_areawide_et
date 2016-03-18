### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("reset", "doParallel")
Orcs::loadPkgs(lib)

## set working directory
Orcs::setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/",
              path_ext = "kilimanjaro/evapotranspiration/")

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


### processing: extract relevant sds -------------------------------------------

## reference extent
rst_ref <- raster("data/MYD09Q1.006/ndvi/NDVI_MYD09Q1.A2013001.sur_refl.tif")
rst_ref <- trim(projectRaster(rst_ref, crs = "+init=epsg:4326"))

rst_kili <- stack("data/kili_aerial_ll.tif")
spt_kili <- data.frame(t(sp::bbox(rst_kili)))
coordinates(spt_kili) <- ~ s1 + s2
projection(spt_kili) <- "+init=epsg:4326"
ext_kili <- extent(rst_kili)
spy_kili <- as(ext_kili, "SpatialPolygons")

## loop over products
prm <- c("Solar_Zenith", "Scan_Start_Time", "Water_Vapor_Near_Infrared", 
         "Cloud_Mask_QA", "Quality_Assurance_Near_Infrared")

lst_prm <- foreach(product = c("MOD05_L2.006", "MYD05_L2.006")) %do% {
  
  ## status message
  cat("Commencing with", product, "...\n")
  
  ## discard nighttime scenes (infrared-based, i.e. 5 km spatial resolution)
  fls <- list.files(paste0("data/", product, "/Swath"), pattern = ".hdf$", 
                    full.names = TRUE)
  
  tms <- as.numeric(getSwathTime(fls)) 
  fls <- fls[tms >= 700 & tms <= 1800]

  ## discard scenes not covering reference extent
  # lst_ext <- list()
  # for (i in 1:length(fls)) {
  #   cat("Processing file", i, "...\n")
  #   lst_ext[[i]] <- try(reset::getSwathExtent(fls[i]), silent = TRUE)
  # }
  # saveRDS(lst_ext, paste0("data/", product, "/Swath/extents.rds"))
  
  lst_ext <- readRDS(paste0("data/", product, "/Swath/extents.rds"))
  
  inside <- sapply(lst_ext, function(i) {
    if (class(i) != "try-error") {
      rgeos::gCovers(as(i, "SpatialPolygons"), spy_kili)
    } else {
      FALSE
    }
  })
  
  if (any(!inside)) fls <- fls[inside]

  # ## extract precipitable water incl. quality layers
  # lst_prw <- foreach(i = fls, .packages = "reset") %dopar%
  #   reset::getSwathSDS(i, prm = prm[3:5], ext = rst_kili,
  #                      dsn = paste0("data/", product))
  
  ## extract and resample sun zenith angle and scan start time
  lst_szn <- foreach(i = fls, .packages = "reset") %dopar%
    reset::getSwathSDS(i, prm = prm[1:2], ext = rst_kili, template = rst_ref,
                       dsn = paste0("data/", product))
  
  return(lst_prw)
}


### processing: quality control ------------------------------------------------

## loop over products
lst_prw <- foreach(product = c("MOD05_L2.006", "MYD05_L2.006")) %do% {
          
  ## status message
  cat("Commencing with", product, "...\n")
  
  # import required parameters
  lst_prm <- foreach(param = prm[2:4], .packages = "raster") %dopar% {
    fls_prm <- list.files(paste("data", product, param, sep = "/"), 
                          full.names = TRUE, pattern = paste0(param, ".tif$"))
    lapply(fls_prm, raster)
  }
  
  # apply cloud mask
  dir_qc1 <- paste("data", product, prm[2], "qc1", sep = "/")
  if (!dir.exists(dir_qc1)) dir.create(dir_qc1)
  
  lst_qc1 <- foreach(i = lst_prm[[1]], j = lst_prm[[2]], 
                     .packages = "Rsenal") %dopar% {
                       
    fls_qc1 <- paste0(dir_qc1, "/QC1_", names(i), ".tif")
    if (file.exists(fls_qc1)) {
      raster::raster(fls_qc1)      
    } else {
                         
    raster::overlay(i, j, fun = function(x, y) {
      id <- sapply(y[], function(l) {
        bin <- Rsenal::number2binary(l, 8)
        status <- bin[8] == 1
        cloudiness <- paste(bin[6:7], collapse = "") %in% c("10", "11")
        day_night <- bin[5] == 1
        sunglint <- bin[4] == 1
        
        all(status, cloudiness, day_night, sunglint)
      })
      
      x[!id] <- NA
      return(x)
    }, filename = fls_qc1, overwrite = TRUE, format = "GTiff")
    }
  }
  
  # consider near-infrared qa
  dir_qc2 <- paste("data", product, prm[2], "qc2", sep = "/")
  if (!dir.exists(dir_qc2)) dir.create(dir_qc2)
  
  lst_qc2 <- foreach(i = lst_qc1, j = lst_prm[[3]], 
                     .packages = "Rsenal") %dopar% {
                       
    fls_qc2 <- paste0(dir_qc2, "/QC2_", names(i), ".tif")
    if (file.exists(fls_qc2)) {
      raster::raster(fls_qc2)      
    } else {
                       
    raster::overlay(i, j, fun = function(x, y) {
      
      id <- sapply(y[], function(l) {
        bin <- Rsenal::number2binary(l, 8)
        usefulness <- bin[8] == 1
        confidence <- paste(bin[5:7], collapse = "") %in% c("010", "011")

        all(usefulness, confidence)
      })
      
      x[!id] <- NA
      return(x)
    }, filename = fls_qc2, overwrite = TRUE, format = "GTiff")
    }
  }
  
  return(lst_qc2)
}


### processing: resample -------------------------------------------------------

## loop over products
lst_res <- foreach(product = c("MOD05_L2.006", "MYD05_L2.006")) %do% {
  
  ## status message
  cat("Commencing with", product, "...\n")
  
  # import required parameters
  fls_ptw <- list.files(paste("data", product, prm[2], "qc2", sep = "/"), 
                        full.names = TRUE, pattern = paste0(prm[2], ".tif$"))

  # target folder and files  
  dir_res <- paste("data", product, prm[2], "res", sep = "/")
  if (!dir.exists(dir_res)) dir.create(dir_res)
  
  fls_res <- paste0(dir_res, "/RES_", basename(fls_ptw))
  
  lst_res <- foreach(i = fls_ptw, j = fls_res, .packages = "raster") %dopar% {
    if (file.exists(j)) {
      raster::raster(j)
    } else {
      rst <- raster::raster(i)
      rst_prj <- raster::projectRaster(rst, crs = raster::projection(rst_ref))
      raster::resample(rst_prj, rst_ref, filename = j, 
                       format = "GTiff", overwrite = TRUE)
    }
  }
  
  raster::stack(lst_res)
}
  