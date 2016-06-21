### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("reset", "satellite", "doParallel", "Rsenal", "MODIS")
Orcs::loadPkgs(lib)

## set working directory
Orcs::setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "H:/",
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
         "Cloud_Mask_QA", "Quality_Assurance_Near_Infrared", "Sensor_Zenith")

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

  ## extract precipitable water incl. quality layers
  foreach(i = fls, .packages = "reset") %dopar%
    getSwathSDS(i, prm = prm, ext = rst_kili, template = rst_ref,
                dsn = paste0("data/", product))
}


### processing: quality control ------------------------------------------------

## loop over products
lst_prw <- foreach(product = c("MOD05_L2.006", "MYD05_L2.006")) %do% {
          
  ## status message
  cat("Commencing with", product, "...\n")
  
  # import required parameters
  lst_prm <- foreach(param = prm[3:5], .packages = lib) %dopar% {
    fls_prm <- list.files(paste("data", product, param, sep = "/"), 
                          full.names = TRUE, pattern = paste0(param, ".tif$"))
    lapply(fls_prm, raster)
  }
  
  ## apply cloud mask
  dir_qc1 <- paste("data", product, prm[3], "qc1", sep = "/")
  if (!dir.exists(dir_qc1)) dir.create(dir_qc1)
  
  lst_qc1 <- foreach(i = lst_prm[[1]], j = lst_prm[[2]], 
                     .packages = lib) %dopar% {
                       
    fls_qc1 <- paste0(dir_qc1, "/QC1_", names(i), ".tif")
    if (file.exists(fls_qc1)) {
      raster(fls_qc1)      
    } else {
                         
    overlay(i, j, fun = function(x, y) {
      id <- sapply(y[], function(l) {
        bin <- Rsenal::number2binary(l, 8)
        status <- bin[8] == 1
        cloudiness <- paste(bin[6:7], collapse = "") != "00"
        day_night <- bin[5] == 1
        sunglint <- bin[4] == 1
        
        all(status, cloudiness, day_night, sunglint)
      })
      
      x[!id] <- NA
      return(x)
    }, filename = fls_qc1, overwrite = TRUE, format = "GTiff")
    }
  }
  
  ## consider near-infrared qa
  dir_qc2 <- paste("data", product, prm[3], "qc2", sep = "/")
  if (!dir.exists(dir_qc2)) dir.create(dir_qc2)
  
  lst_qc2 <- foreach(i = lst_qc1, j = lst_prm[[3]], 
                     .packages = "Rsenal") %dopar% {
                       
    fls_qc2 <- paste0(dir_qc2, "/QC2_", names(i), ".tif")
    if (file.exists(fls_qc2)) {
      raster(fls_qc2)      
    } else {
                       
    overlay(i, j, fun = function(x, y) {
      
      id <- sapply(y[], function(l) {
        bin <- Rsenal::number2binary(l, 8)
        usefulness <- bin[8] == 1
        confidence <- paste(bin[5:7], collapse = "") != "000"

        all(usefulness, confidence)
      })
      
      x[!id] <- NA
      return(x)
    }, filename = fls_qc2, overwrite = TRUE, format = "GTiff")
    }
  }
  
  ## discard neighboring pixels
  dir_adj <- paste("data", product, prm[3], "adj", sep = "/")
  if (!dir.exists(dir_adj)) dir.create(dir_adj)
  
  lst_adj <- foreach(i = lst_qc2, .packages = lib) %dopar% {
    fls_adj <- paste0(dir_adj, "/", names(i), ".tif")
    if (file.exists(fls_adj)) {
      raster(fls_adj)
    } else {
      msk <- focal(i, w = matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3),
                   fun = function(...) any(is.na(...)))
      
      overlay(i, msk, fun = function(x, y) {
        x[y[] == 1] <- NA
        return(x)
      }, filename = fls_adj, format = "GTiff", overwrite = TRUE)
    }
  }

  return(lst_adj)
}


### processing: resample -------------------------------------------------------

## loop over products
lst_res <- foreach(product = c("MOD05_L2.006", "MYD05_L2.006")) %do% {
  
  ## status message
  cat("Commencing with", product, "...\n")
  
  # import required parameters
  fls_ptw <- list.files(paste("data", product, prm[3], "adj", sep = "/"), 
                        full.names = TRUE, pattern = paste0(prm[3], ".tif$"))

  # target folder and files  
  dir_res <- paste("data", product, prm[3], "res", sep = "/")
  if (!dir.exists(dir_res)) dir.create(dir_res)
  
  fls_res <- paste0(dir_res, "/RES_", basename(fls_ptw))
  
  stack(foreach(i = fls_ptw, j = fls_res, .packages = lib) %dopar% {
    if (file.exists(j)) {
      raster(j)
    } else {
      rst <- raster(i)
      
      # reproject
      if (!compareCRS(rst_ref, rst))
        rst <- projectRaster(rst, crs = projection(rst_ref))
      
      # resample
      resample(rst, rst_ref, filename = j, format = "GTiff", overwrite = TRUE)
    }
  })
}
  