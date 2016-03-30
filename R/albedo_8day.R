### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## source functions
source("R/uniformExtent.R")

## set working directory
Orcs::setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/",
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


### data download --------------------------------------------------------------

## combinded terra/aqua product
runGdal(product = "MCD43A3", collection = "005", 
        begin = "2013001", end = "2015365", tileH = 21, tileV = 9,
        SDSstring = paste(c(rep(0, 19), 1), collapse = ""), job = "MCD43A3.005")

## quality information
runGdal(product = "MCD43A2", collection = "005", 
        begin = "2013001", end = "2015365", tileH = 21, tileV = 9,
        SDSstring = "1000", job = "MCD43A2.005")


### crop layers ----------------------------------------------------------------

## retrieve study extent
ext_crp <- uniformExtent(verbose = FALSE)

lst_crp <- foreach(product = c("MCD43A3.005", "MCD43A2.005"), 
                   .packages = c("raster", "MODIS")) %dopar% {
          
  ## which product is currently being processed                   
  quality <- product == "MCD43A2.005"          
  
  ## list and import available files
  fls <- list.files(paste0(getOption("MODIS_outDirPath"), "/", product),
                    pattern = ifelse(!quality, "WSA_shortwave.tif$", "Quality.tif$"), 
                                     full.names = TRUE)
  rst <- raster::stack(fls)
  
  ## perform crop if required
  dir_out <- paste0("data/", product, "/crp")
  if (!dir.exists(dir_out)) dir.create(dir_out)
  
  fls_out <- paste(dir_out, basename(fls), sep = "/")
  
  if (all(file.exists(fls_out))) {
    raster::stack(fls_out)
  } else {
    rst_out <- raster::crop(rst, ext_crp, snap = "out")
    
    # assign data type
    raster::dataType(rst_out) <- ifelse(!quality, "INT2U", "INT1U")

    # apply scale factor
    if (!quality) rst_out <- rst_out * 0.001
    
    # save and return cropped layers
    lst_out <- lapply(1:nlayers(rst_out), function(j)
      raster::writeRaster(rst_out[[j]], filename = fls_out[j],
                          format = "GTiff", overwrite = TRUE)
    )
    
    rst_crp <- raster::stack(lst_out)
  }
}


# ### gap-filling ----------------------------------------------------------------
# ### kolmogorov-zurbenko adaptive
# 
# mat_adj <- as.matrix(rst_adj)
# 
# mat_gf <- t(apply(mat_adj, 1, FUN = function(x) kza(x, m = 5)$kz))
# rst_gf <- setValues(rst_adj, mat_gf)
# 
# dir_gf <- "data/MCD15A2H.006/gf"
# if (!dir.exists(dir_gf)) dir.create(dir_gf)
# 
# lst_gf <- foreach(i = unstack(rst_gf), j = unstack(rst_adj)) %do% 
#   raster::writeRaster(i, filename = paste(dir_gf, names(j), sep = "/"), 
#                       format = "GTiff", overwrite = TRUE)
# 
# rst_gf <- stack(lst_gf)

