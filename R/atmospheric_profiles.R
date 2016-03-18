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


### rearrange data -------------------------------------------------------------

## loop over products
for (product in c("MOD07_L2", "MYD07_L2")) {
  
  ## status message
  cat("Commencing with processing of", product, "...\n")
  
  ## list available files and remove duplicates
  fls <- list.files("md07/raw006", full.names = TRUE, recursive = TRUE, 
                    pattern = paste(product, "006", ".hdf$", sep = ".*"))
  
  if (anyDuplicated(basename(fls)) != 0)
    fls <- fls[-which(duplicated(basename(fls)))]
  
  ## sort by datetime and stop if any date is missing
  meta <- Orcs::list2df(strsplit(basename(fls), "\\."), stringsAsFactors = FALSE)
  names(meta) <- c("product", "date", "time", "collection", "processing", "filetype")
  
  meta$datetime <- strptime(paste(meta$date, meta$time), format = "A%Y%j %H%M")
  fls <- fls[order(meta$datetime)]; # meta <- meta[order(meta$datetime), ]
  
  if (length(unique(meta$date)) != 
      length(seq(as.Date("2013-01-01"), as.Date("2015-12-31"), "day")))
    stop("There are dates missing. Please check available files.\n")
  
  ## copy files to destination folder
  file.copy(fls, paste0("data/", product, ".006/", basename(fls)))
}

## remove non-required files
for (product in c("MOD07_L2", "MYD07_L2")) {
  
  ## list available files and remove duplicates
  fls <- list.files("md07/raw006", full.names = TRUE, recursive = TRUE, 
                    pattern = paste(product, "006", ".hdf$", sep = ".*"))
  
  file.remove(fls)
}


### process coordinates --------------------------------------------------------

## reference extent
rst_ref <- raster("data/MYD09Q1.006/ndvi/NDVI_MYD09Q1.A2013001.sur_refl.tif")
rst_ref <- trim(projectRaster(rst_ref, crs = "+init=epsg:4326"))

rst_kili <- stack("data/kili_aerial_ll.tif")
spt_kili <- data.frame(t(sp::bbox(rst_kili)))
coordinates(spt_kili) <- ~ s1 + s2
projection(spt_kili) <- "+init=epsg:4326"
ext_kili <- extent(rst_kili)
spy_kili <- as(ext_kili, "SpatialPolygons")

## required parameters
prm <- c("Skin_Temperature", "Retrieved_Height_Profile",
         "Retrieved_Temperature_Profile", "Retrieved_Moisture_Profile")

## loop over products
for (product in c("MOD07_L2.006", "MYD07_L2.006")) {
  
  ## status message
  cat("Commencing with", product, "...\n")
  
  ## discard nighttime scenes (infrared-based, i.e. 5 km spatial resolution)
  fls <- list.files(paste0("data/", product), pattern = ".hdf$", 
                    full.names = TRUE)

  tms <- as.numeric(getSwathTime(fls)) 
  fls <- fls[tms >= 600 & tms <= 1800]
  
  ## discard scenes not covering reference extent
  # lst_ext <- lapply(1:length(fls), function(i) {
  #   try(reset::getSwathExtent(fls[i]), silent = TRUE)
  # })
  # saveRDS(lst_ext, paste0("data/", product, "/extents.rds"))
  lst_ext <- readRDS(paste0("data/", product, "/extents.rds"))
  
  inside <- sapply(lst_ext, function(i) {
    if (class(i) != "try-error") {
      rgeos::gCovers(as(i, "SpatialPolygons"), spy_kili)
    } else {
      FALSE
    }
  })
  
  if (any(!inside)) fls <- fls[inside]
  
  # ## relevant atmospheric paramters
  # lst_out <- getSwathSDS(fls, prm = prm, ext = rst_kili, 
  #                        dsn = "data/MYD07_L2.006", verbose = TRUE, 
  #                        cores = 3L)
  
  ## geopotential heights
  foreach(i = fls, .packages = "reset") %dopar%
    getSwathSDS(i, prm = "Retrieved_Temperature_Profile", ext = rst_kili, 
                template = rst_ref, dsn = paste0("data/", product))[[1]][[1]]
}