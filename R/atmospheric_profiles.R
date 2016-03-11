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

## pressure levels
p <- c(5, 10, 20, 30, 50, 70, 100, 150, 200, 250, 
       300, 400, 500, 620, 700, 780, 850, 920, 950, 1000)


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
rst_kili <- stack("data/kili_aerial_ll.tif")
spt_kili <- data.frame(t(sp::bbox(rst_kili)))
coordinates(spt_kili) <- ~ s1 + s2
projection(spt_kili) <- "+init=epsg:4326"
ext_kili <- extent(rst_kili)
spy_kili <- as(ext_kili, "SpatialPolygons")

## required parameters
prm <- c("Solar_Zenith", "Skin_Temperature", "Surface_Pressure", 
         "Retrieved_Temperature_Profile", "Retrieved_Moisture_Profile", 
         "Water_Vapor")

## loop over products
for (product in c("MOD07_L2.006", "MYD07_L2.006")) {
  
  ## discard scenes not covering reference extent
  fls <- list.files(paste0("data/", product), pattern = ".hdf$", 
                    full.names = TRUE)
  
  # lst_ext <- list()
  # for (i in 1:length(fls)) {
  #   cat("Processing file", i, "...\n")
  #   lst_ext[[i]] <- try(reset::getAtmosProfBbox(fls[i]), silent = TRUE)
  # }
  # saveRDS(lst_ext, paste0("data/", product, "/inside.rds"))
  
  lst_ext <- readRDS(paste0("data/", product, "/inside.rds"))
  
  inside <- sapply(lst_ext, function(i) {
    if (class(i) != "try-error") {
      rgeos::gCovers(as(i, "SpatialPolygons"), spy_kili)
    } else {
      FALSE
    }
  })
  
  fls <- fls[inside]
  ext <- lst_ext[inside]
  
  lst_out <- getAtmosProfParam(fls[1:10], prm = "Surface_Elevation", ext = rst_kili, 
                               dsn = "data/MYD07_L2.006", verbose = TRUE, 
                               cores = 3L)
  
  ## loop over files
  lst_out <- foreach(h = as.list(fls), g = ext) %do% {
    
    ## status message
    cat(paste0("Commencing with file #", which(h == fls), ":"), h, "\n")
    
    ### process parameters -----------------------------------------------------
    
    ## read metadata
    info <- suppressWarnings(GDALinfo(h, returnScaleOffset = FALSE))
    meta <- attr(info, "mdata")
    sds <- attr(info, "subdsmdata")
    
    ## loop over parameters
    lst_prm <- foreach(i = prm, 
                       .packages = c("Orcs", "rgdal", "raster")) %dopar% {
                         
      # create destination folder and file
      dir_prm <- paste0("data/", product, "/", i)
      if (!dir.exists(dir_prm)) dir.create(dir_prm)
      
      fls_prm <- paste0(dir_prm, Orcs::pureBasename(h, slash = TRUE), "_", i)
      
      # identify relevant sds
      sds_prm <- sds[grep(i, sds)[1]]
      sds_prm <- sapply(strsplit(sds_prm, "="), "[[", 2)
      
      # retrieve scale and offset of current parameter
      info_prm <- suppressWarnings(rgdal::GDALinfo(sds_prm))
      meta_prm <- attr(info_prm, "mdata")
      scale_offset_prm <- sapply(c("^scale_factor", "^add_offset"), function(j) {
        prm <- meta_prm[grep(j, meta_prm)]
        as.numeric(strsplit(prm, "=")[[1]][2])
      })
      
      # single-layer bands
      if (info_prm[["bands"]] == 1) {
        # rasterize band
        rst_prm <- suppressWarnings(
          raster::raster(rgdal::readGDAL(sds_prm, as.is = TRUE, silent = TRUE))
        )
        
        # multi-layer bands
      } else {
        # rasterize band
        rst_prm <- suppressWarnings(
          raster::stack(rgdal::readGDAL(sds_prm, as.is = TRUE, silent = TRUE))
        )
      }
      
      # apply offset and scale factor
      rst_prm <- raster::calc(rst_prm, fun = function(x) {
        (x - scale_offset_prm[2]) * scale_offset_prm[1]
      })
      
      # set extent and crs
      raster::extent(rst_prm) <- g
      raster::projection(rst_prm) <- "+init=epsg:4326"
      
      # crop by reference extent
      rst_prm <- raster::crop(rst_prm, rst_kili, snap = "out")
      
      # write to file and return
      raster::writeRaster(rst_prm, filename = fls_prm, 
                          format = "GTiff", overwrite = TRUE)
    }
  }
}