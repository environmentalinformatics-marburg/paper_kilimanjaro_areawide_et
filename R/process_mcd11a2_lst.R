## working directory
library(Orcs)
setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/", 
        path_ext = "kilimanjaro/evapotranspiration/")

## packages
lib <- c("Rsenal", "rgdal", "doParallel", "MODIS")
jnk <- sapply(lib, function(x) require(x, character.only = TRUE, quietly = TRUE))

## parallelize
supcl <- makeCluster(3)
registerDoParallel(supcl)

## reference extent
rst_ref <- kiliAerial(rasterize = TRUE, minNumTiles = 20L)

################################################################################
### quality control
################################################################################

lst_lst <- lapply(c("MOD11A2", "MYD11A2"), function(h) {

  ## import and crop lst images
  fls <- list.files(paste0(options()$MODIS_outDirPath, "LST1km"), 
                    pattern = paste0("^", h, ".*Day"), full.names = TRUE)
  
  fls_lst <- fls[seq(1, length(fls), 2)]
  rst_lst <- stack(fls_lst)
  lst_lst_crp <- foreach(i = unstack(rst_lst), .packages = lib, 
                         .export = ls(envir = globalenv())) %dopar% {
    file_out <- paste0("data/", h, "/crp/CRP_", names(i), ".tif")
    crop(i, rst_ref, filename = file_out, format = "GTiff", overwrite = TRUE) * 0.02
  }
  rst_lst_crp <- stack(lst_lst_crp)
  
  ## import and crop qa images
  fls_qa <- fls[seq(2, length(fls), 2)]
  rst_qa <- stack(fls_qa)
  lst_qa_crp <- foreach(i = unstack(rst_qa), .packages = lib, 
                        .export = ls(envir = globalenv())) %dopar% {
    file_out <- paste0("data/", h, "/crp/CRP_", names(i), ".tif")
    crop(i, rst_ref, filename = file_out, format = "GTiff", overwrite = TRUE)
  }
  rst_qa_crp <- stack(lst_qa_crp)

  ## overlay
  lst_lst_qc <- foreach(i = 1:nlayers(rst_lst_crp), .packages = lib, 
                        .export = ls(envir = globalenv())) %dopar% {
    
    file_out <- paste0("data/", h, "/qc/QC_", names(rst_lst_crp[[i]]), ".tif")
    
    overlay(rst_lst_crp[[i]], rst_qa_crp[[i]], fun = function(x, y) {
      
      bits <- sapply(y, function(j) number2binary(j, no_bits = 8, to_char = TRUE))
      
      out <- numeric(length(y))
      int_keep <- which(substr(bits, 5, 6) == "00" |                     # overall
                          (substr(bits, 5, 6) == "00" &                  # data quality  
                             substr(bits, 3, 4) %in% c("00", "01") &     # emissivity
                             substr(bits, 1, 2) %in% c("00", "01")))     # lst
      
      out[int_keep] <- x[int_keep]
      out[-int_keep] <- NA
      return(out)
    }, filename = file_out, format = "GTiff", overwrite = TRUE)
  }
  
  return(stack(lst_lst_qc))
})


################################################################################
### merge terra and aqua
################################################################################

# Calculate daily maximum daytime temperatures from Terra and Aqua (usually 
# Aqua higher than Terra, but Terra will be taken in case of missing Aqua cells)
lst_mrg <- 
  foreach(i = unstack(lst_lst[[1]]), j = unstack(lst_lst[[2]]), 
          .packages = lib) %dopar% {
            file_out <- paste0("data/MCD11A2/mrg/MRG_", names(i), ".tif")
            file_out <- gsub("MOD11A2", "MCD11A2", file_out)
            
            overlay(i, j, fun = function(...) max(..., na.rm = TRUE), 
                    filename = file_out, format = "GTiff", overwrite = TRUE)
          }

rst_mrg <- stack(lst_mrg)


################################################################################
### resample to 500 m spatial resolution
################################################################################

# Resample 1km LST raster do 250m (SAVI, NDVI, LAI, etc.) resolution
rst_lai <- raster("data/MOD15A2H/crp/CRP_MOD15A2H.A2012001.Lai_500m.tif")

lst_rsmpl <- foreach(i = unstack(rst_mrg)) %do% {
  file_out <- paste0("data/MCD11A2/rsmpl/RSMPL_", names(i), ".tif")
  resample(i, rst_lai, 
           filename = file_out, format = "GTiff", overwrite = TRUE)
}

rst_rsmpl <- stack(lst_rsmpl)

# reimport
fls_rsmpl <- list.files("data/MCD11A2/rsmpl", pattern = "^RSMPL.*.tif$", 
                        full.names = TRUE)
rst_rsmpl <- stack(fls_rsmpl)

# Deregister parallel backend
stopCluster(supcl)
