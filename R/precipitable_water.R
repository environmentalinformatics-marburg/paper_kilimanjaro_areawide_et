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
rst_kili <- stack("data/kili_aerial_ll.tif")
spt_kili <- data.frame(t(sp::bbox(rst_kili)))
coordinates(spt_kili) <- ~ s1 + s2
projection(spt_kili) <- "+init=epsg:4326"
ext_kili <- extent(rst_kili)
spy_kili <- as(ext_kili, "SpatialPolygons")

## loop over products
prm <- c("Solar_Zenith", "Water_Vapor_Near_Infrared", 
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
  #   lst_ext[[i]] <- try(reset::getAtmosProfBbox(fls[i]), silent = TRUE)
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
  
  if (any(!inside)) {
    fls <- fls[inside]
    ext <- lst_ext[inside]
  }
  
  ## extract sds
  foreach(i = fls, .packages = "reset") %dopar%
    reset::getAtmosProfParam(i, prm = prm, ext = rst_kili, 
                             dsn = paste0("data/", product))
}


### processing: quality control ------------------------------------------------

## loop over products
lst_prw <- foreach(product = c("MOD05_L2.006", "MYD05_L2.006"), 
                   .packages = "raster") %dopar% {
          
  # import required parameters                              
  lst_prm <- foreach(param = prm[2:4]) %do% {
    fls_prm <- list.files(paste("data", product, param, sep = "/"), 
                          full.names = TRUE, pattern = paste0(param, ".tif$"))
    raster::stack(fls_prm)
  }
  
  # apply cloud mask
  rst_qc1 <- raster::overlay(lst_prm[[1]], lst_prm[[2]], fun = function(x, y) {
    ...
  })
  
  # consider near-infrared qa
  raster::overlay(rst_qc1, lst_prm[[3]], fun = function(x, y) {
    ...
  })
}