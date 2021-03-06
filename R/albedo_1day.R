### environment -----

## clear workspace
rm(list = ls(all = TRUE))

## source functions
source("R/uniformExtent.R")

## set working directory
Orcs::setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "H:/",
              path_ext = "kilimanjaro/evapotranspiration/")

## load packages
lib <- c("Rsenal", "MODIS", "doParallel")
Orcs::loadPkgs(lib)

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

## modis options
MODISoptions(localArcPath = "MODIS_ARC", outDirPath = "MODIS_ARC/PROCESSED",
             outProj = "+init=epsg:21037")


### data download -----

## combinded terra/aqua product incl. quality layer
runGdal(product = "MCD43A3", collection = "006", job = "MCD43A3.006", 
        begin = "2015011", end = "2015365", tileH = 21, tileV = 9,
        SDSstring = paste(c(rep(0, 9), 1, rep(0, 19), 1), collapse = ""))


### crop layers -----

## target folders
dir_prd <- "data/MCD43A3.006"
if (!dir.exists(dir_prd)) dir.create(dir_prd)

dir_crp <- "data/MCD43A3.006/crp"
if (!dir.exists(dir_crp)) dir.create(dir_crp)

## retrieve study extent
ext_crp <- uniformExtent(verbose = FALSE)

lst_lyr <- lapply(c("WSA_shortwave", "Quality_shortwave"), function(layer) {
          
  ## list and import available files
  fls <- list.files(paste0(getOption("MODIS_outDirPath"), "/MCD43A3.006"),
                    pattern = layer, full.names = TRUE)
  rst <- stack(fls)
  
  ## perform crop if required
  fls_crp <- paste(dir_crp, basename(fls), sep = "/")
  
  lst_crp <- foreach(i = fls_crp, j = 1:length(fls_crp), .packages = "raster", 
                     .export = ls(envir = globalenv())) %dopar% {
    
    if (file.exists(i)) {
      raster(i)
    } else {
      rst_crp <- crop(rst[[j]], ext_crp, snap = "out")
      dataType(rst_crp) <- ifelse(layer == "WSA_shortwave", "INT2U", "INT1U")
      if (layer == "WSA_shortwave") rst_crp <- rst_crp * 0.001
      writeRaster(rst_crp, filename = i, format = "GTiff", overwrite = TRUE)
    }
  }
  
  stack(lst_crp)
})


### gap-filling -----

mat_alb <- as.matrix(lst_lyr[[1]])

mat_fll <- foreach(i = 1:nrow(mat_alb), .combine = "rbind") %do% {
  val <- mat_alb[i, ]
  
  if (!all(is.na(val))) {
    
    id <- is.na(val)
    id_vld <- which(!id); id_inv <- which(id)
    
    while(length(id_inv) > 0) {
      val[id_inv] <- kza(val, m = 3)$kz[id_inv]
      
      id <- is.na(val)
      id_vld <- which(!id); id_inv <- which(id)
    }
  }
  
  return(val)
}

rst_fll <- setValues(lst_lyr[[1]], mat_fll)

drs_gf <- paste0("data/MCD43A3.006/gf")
if (!dir.exists(drs_gf)) dir.create(drs_gf)

rst_gf <- stack(foreach(i = 1:nlayers(rst_fll), .packages = "raster") %dopar% {
  writeRaster(rst_fll[[i]], 
              filename = paste0(drs_gf, "/", names(rst_fll[[i]]), ".tif"), 
              format = "GTiff", overwrite = TRUE)
})

## close parallel backend
stopCluster(cl)
