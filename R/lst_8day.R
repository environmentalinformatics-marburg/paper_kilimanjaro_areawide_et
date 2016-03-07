### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## set working directory
Orcs::setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/",
              path_ext = "kilimanjaro/evapotranspiration/")

## load packages
lib <- c("Rsenal", "MODIS", "doParallel")
Orcs::loadPkgs(lib)

## source functions
source("src/uniformExtent.R")

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

## modis options
MODISoptions(localArcPath = "MODIS_ARC", outDirPath = "MODIS_ARC/PROCESSED",
             outProj = "+init=epsg:21037")


### data download --------------------------------------------------------------

## loop over single product
lst <- lapply(c("MOD11A2", "MYD11A2"), function(product) {
  
  ## status message
  cat("Commencing with the processing of", product, "...\n")
  
  # ## download data
  # runGdal(product = product,
  #         collection = "005",
  #         begin = "2013001", end = "2015365", tileH = 21, tileV = 9,
  #         SDSstring = "111111110000", job = paste0(product, ".006"))
  

  ### crop layers ----------------------------------------------------------------
  
  ## retrieve study extent
  ext_crp <- uniformExtent(verbose = FALSE)
  
  ## setup output folder
  dir_out <- paste0("data/", product, ".005")
  if (!dir.exists(dir_out)) dir.create(dir_out)
  
  ## perform crop
  pattern <- c("Day_1km", "QC_Day", "Day_view_time", "Day_view_angl", 
               "Night_1km", "QC_Night", "Night_view_time", "Night_view_angl")
  
  # rst_crp <- foreach(i = pattern, .packages = "MODIS") %dopar% {
  # 
  #   # list and import available files
  #   fls <- list.files(paste0(getOption("MODIS_outDirPath"), "/", product, ".005"),
  #                     pattern = paste0(i, ".tif$"), full.names = TRUE)
  #   rst <- raster::stack(fls)
  #   
  #   # crop
  #   dir_crp <- paste0(dir_out, "/crp")
  #   if (!dir.exists(dir_crp)) dir.create(dir_crp)
  #   
  #   fls_crp <- paste(dir_crp, basename(fls), sep = "/")
  #   rst_crp <- raster::crop(rst, ext_crp, snap = "out")
  #   
  #   # if dealing with (day or night) lst bands, convert to 16-bit unsigned 
  #   # integer and apply scale factor of 0.02
  #   if (i %in% c("Day_1km", "Night_1km")) {
  #     raster::dataType(rst_crp) <- "INT2U"
  #     rst_crp <- rst_crp * 0.02
  #     
  #   # else convert to 8-bit unsigned integer
  #   } else {
  #     raster::dataType(rst_crp) <- "INT1U"
  #     
  #     # if dealing with (day or night) view time, apply scale factor of 0.1
  #     if (i %in% c("Day_view_time", "Night_view_time")) {
  #       rst_crp <- rst_crp * 0.1
  #     
  #     # if dealing with (day or night) view angle, apply offset of -65    
  #     } else if (i %in% c("Day_view_angl", "Night_view_angl")) {
  #       rst_crp <- rst_crp - 65
  #     }
  #   }
  #   
  #   # save and return cropped layers
  #   lst_crp <- lapply(1:nlayers(rst_crp), function(j)
  #     raster::writeRaster(rst_crp[[j]], filename = fls_crp[j],
  #                         format = "GTiff", overwrite = TRUE)
  #   )
  #   
  #   raster::stack(lst_crp)
  # }
  
  ## reimport cropped files
  dir_crp <- paste0(dir_out, "/crp")
  
  rst_crp <- foreach(i = pattern, .packages = "raster") %dopar% {
    
    # list and import available files
    fls_crp <- list.files(dir_crp, pattern = paste0(i, ".tif$"), full.names = TRUE)
    raster::stack(fls_crp)
  }
  
  
  ### quality control ----------------------------------------------------------
  ### discard cloudy pixels based on companion quality information ('QC_Day', 
  ### 'QC_Night')
  
  dir_qc <- paste0(dir_out, "/qc")
  if (!dir.exists(dir_qc)) dir.create(dir_qc)
  
  # ## perform quality check for day and night separately
  # lst_qc <- foreach(i = rst_crp[c(1, 5)], j = rst_crp[c(2, 6)]) %do% {
  # 
  #   ## loop over layers
  #   lst_out <- foreach(k = 1:nlayers(i), .packages = lib) %dopar%
  #     overlay(i[[k]], j[[k]], fun = function(x, y) {
  #       id <- sapply(y[], function(l) {
  #         bin <- number2binary(l, 8, TRUE)
  #         mandatory_qa <- substr(bin, 7, 8) == "00"
  #         data_quality <- substr(bin, 5, 6) == "00"
  #         emis_error <- substr(bin, 3, 4) == "00"
  #         lst_error <- substr(bin, 1, 2) == "00"
  # 
  #         all(mandatory_qa, data_quality, emis_error, lst_error)
  #       })
  # 
  #       x[!id] <- NA
  #       return(x)
  #     }, filename = paste(dir_qc, names(i[[k]]), sep = "/"),
  #     overwrite = TRUE, format = "GTiff")
  # 
  #   raster::stack(lst_out)
  # }
  
  ## reimport quality-controlled files
  lst_qc <- lapply(pattern[c(1, 5)], function(i) {
    fls_qc <- list.files(dir_qc, pattern = paste0(i, ".tif$"), full.names = TRUE)
    raster::stack(fls_qc)
  })
  
  return(lst_qc)
})


### combined product -----------------------------------------------------------

lst_comb <- foreach(i = lst[[1]], j = lst[[2]]) %do% {
  
  mat_terra <- raster::as.matrix(i)
  mat_aqua <- raster::as.matrix(j)

  mat_comb <- foreach(k = 1:nrow(mat_terra), .combine = "rbind") %do% {
    dat <- data.frame(x = mat_terra[k, ], y = mat_aqua[k, ])
    
    if (sum(complete.cases(dat)) >= (.2 * nrow(dat))) {
      mod <- lm(y ~ x, data = dat)
      
      id <- which(!is.na(dat$x) & is.na(dat$y))
      newdata <- data.frame(x = dat$x[id])
      dat$y[id] <- predict(mod, newdata)
    }
      
    return(dat$y)
  }
  
  rst_comb <- raster::setValues(i, mat_comb)
  
  dir_comb <- "data/MCD11A2.005/"
  if (!dir.exists(dir_comb)) dir.create(dir_comb)
  
  fls_comb <- paste0(dir_comb, gsub("^MYD", "MCD", names(j)))
  lst_comb <- foreach(i = 1:ncol(mat_aqua), .packages = "raster") %dopar%
    raster::writeRaster(rst_comb[[i]], filename = fls_comb[i], 
                        format = "GTiff", overwrite = TRUE)
  
  raster::stack(lst_comb)
}

## reimport combined files
dir_comb <- "data/MCD11A2.005/"

lst_comb <- foreach(i = pattern[c(1, 5)], .packages = "raster") %dopar% {
  fls_comb <- list.files(dir_comb, pattern = i, full.names = TRUE)
  raster::stack(fls_comb)
}
  