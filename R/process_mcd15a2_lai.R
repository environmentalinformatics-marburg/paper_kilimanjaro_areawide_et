## working directory
library(Orcs)
setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/", 
        path_ext = "kilimanjaro/evapotranspiration/")

lib <- c("Rsenal", "rgdal", "doParallel")
jnk <- sapply(lib, function(x) require(x, character.only = TRUE, quiet = TRUE))

# # Required functions
# source("src/number2binary.R")

## parallelize
supcl <- makeCluster(3)
registerDoParallel(supcl)

## reference extent
rst <- kiliAerial(rasterize = TRUE, minNumTiles = 20L)

################################################################################
### import and crop data

lapply(c("MOD15A2H", "MYD15A2H"), function(h) {
  
lst_lai <- 
  foreach(i = c("Lai_500m.tif$", "Lai_QC.tif$", "Extra_QC.tif$"), 
          j = c(.1, 1, 1), .export = ls(envir = globalenv()), 
          .packages = lib) %dopar% {
            
    # available files        
    fls.lai <- list.files("MODIS_ARC/PROCESSED/LAI500m", 
                          pattern = paste(h, i, sep = ".*"), 
                          full.names = TRUE)
    
    # crop
    lst_rst <- foreach(k = fls.lai) %do% {
      filename <- paste0("data/", h, "/crp/CRP_", basename(k))
      crop(raster(k) * j, rst, 
           filename = filename, format = "GTiff", overwrite = TRUE)
    }
    
    # return stacked images
    stack(lst_rst)
  }

################################################################################
### base quality check ('Lai_QC')

lst_lai_qc <- 
  foreach(i = unstack(lst_lai[[1]]), 
          j = unstack(lst_lai[[2]]), .packages = lib) %dopar% {
    overlay(i, j, 
            fun = function(x, y) {
              index <- sapply(y[], function(i) {
                if (!is.na(i)) {
                  # 8-bit string
                  bit <- number2binary(i, 8)
                  # cloud state
                  state <- paste(bit[c(4, 5)], 
                                 collapse = "") %in% c("00", "11")
                  
                  return(state)
                } else {
                  return(FALSE)
                }
              })
              x[!index] <- NA
              return(x)
            }, filename = paste0("data/", h, "/qcb/QCB_", names(i)), 
            overwrite = TRUE, format = "GTiff")
  }

################################################################################
### extra quality check ('Extra_QC')

lst_lai_qcx <- 
  foreach(i = lst_lai_qc, j = unstack(lst_lai[[3]]), .packages = lib) %dopar% {
    overlay(i, j, 
            fun = function(x, y) {
              index <- sapply(y[], function(i) {
                if (!is.na(i)) {
                  # 8-bit string
                  bit <- number2binary(i, 8)
                  
                  # Snow/ice detected
                  snow <- bit[6] == 0
                  # Cirrus detected
                  cirrus <- bit[4] == 0
                  # Internal cloud mask
                  clouds <- bit[3] == 0
                  # Cloud shadow
                  shadow <- bit[2] == 0
                  
                  return(all(snow, cirrus, clouds, shadow))
                } else {
                  return(FALSE)
                }
              })
              x[!index] <- NA
              return(x)
            }, filename = paste0("data/", h, "/qcx/QCX_", names(i)), 
            overwrite = TRUE, format = "GTiff")
  }
})