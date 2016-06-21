### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## source functions
source("R/uniformExtent.R")

## set working directory
Orcs::setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "H:/",
              path_ext = "kilimanjaro/evapotranspiration/")

## load packages
lib <- c("Rsenal", "MODIS", "doParallel", "kza")
Orcs::loadPkgs(lib)

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

## modis options
MODISoptions(localArcPath = "MODIS_ARC", outDirPath = "MODIS_ARC/PROCESSED",
             outProj = "+init=epsg:21037")


### data download --------------------------------------------------------------

# ## combinded terra/aqua product
# runGdal(product = "MCD15A3H", collection = "006",
#         begin = "2013001", end = "2015365", tileH = 21, tileV = 9,
#         SDSstring = "011100", job = "MCD15A3H.006")


### crop layers ----------------------------------------------------------------

## retrieve study extent
ext_crp <- uniformExtent(verbose = FALSE)

# ## perform crop
# rst_crp <- foreach(i = c("Lai_500m", "Lai_QC", "Extra_QC"),
#                    .packages = "MODIS") %dopar% {
# 
#     # list and import available files
#     fls <- list.files(paste0(getOption("MODIS_outDirPath"), "/MCD15A3H.006"),
#                       pattern = paste0(i, ".tif$"), full.names = TRUE)
#     rst <- raster::stack(fls)
# 
#     # crop
#     dir_prd <- "data/MCD15A3H.006"
#     if (!dir.exists(dir_prd)) dir.create(dir_prd)
#     
#     dir_out <- paste0(dir_prd, "/crp")
#     if (!dir.exists(dir_out)) dir.create(dir_out)
# 
#     fls_out <- paste(dir_out, basename(fls), sep = "/")
#     rst_out <- raster::crop(rst, ext_crp, snap = "out")
# 
#     # convert to 8-bit unsigned integer
#     raster::dataType(rst_out) <- "INT1U"
# 
#     # apply scale factor
#     if (i == "Lai_500m")
#       rst_out <- rst_out * 0.1
# 
#     # save and return cropped layers
#     lst_out <- lapply(1:nlayers(rst_out), function(j)
#       raster::writeRaster(rst_out[[j]], filename = fls_out[j],
#                           format = "GTiff", overwrite = TRUE)
#     )
# 
#     raster::stack(lst_out)
#   }

## reimport cropped files
rst_crp <- foreach(i = c("Lai_500m", "Lai_QC", "Extra_QC"), 
                   .packages = "raster") %dopar% {

  # list and import available files
  dir_crp <- "data/MCD15A3H.006/crp"
  fls_crp <- list.files(dir_crp, pattern = paste0(i, ".tif$"), full.names = TRUE)
  raster::stack(fls_crp)
}


### quality control, step #1: --------------------------------------------------
### discard cloudy pixels based on companion quality information ('Lai_QC')

dir_qc1 <- "data/MCD15A3H.006/qc1"
if (!dir.exists(dir_qc1)) dir.create(dir_qc1)

# ## perform quality check #1
# lst_qc1 <- foreach(i = 1:nlayers(rst_crp[[1]]), .packages = lib) %dopar% {
#     overlay(rst_crp[[1]][[i]], rst_crp[[2]][[i]], fun = function(x, y) {
#       id <- sapply(y[], function(k) {
#         bin <- number2binary(k, 8, TRUE)
#         modland <- substr(bin, 8, 8) == "0"
#         cloud_state <- substr(bin, 4, 5) == "00"
#         confidence <- substr(bin, 1, 3) %in% c("000", "001")
# 
#         all(modland, cloud_state, confidence)
#       })
# 
#       x[!id] <- NA
#       return(x)
#     }, filename = paste(dir_qc1, names(rst_crp[[1]][[i]]), sep = "/"),
#     overwrite = TRUE, format = "GTiff")
# }
# 
# rst_qc1 <- stack(lst_qc1)

## reimport step #1 quality-controlled files
fls_qc1 <- list.files(dir_qc1, pattern = "Lai_500m.tif$", full.names = TRUE)
rst_qc1 <- raster::stack(fls_qc1)


### quality control, step #2: --------------------------------------------------
### discard cloudy pixels based on 'state_250m' flags

dir_qc2 <- "data/MCD15A3H.006/qc2"
if (!dir.exists(dir_qc2)) dir.create(dir_qc2)

# ## perform quality check #2
# lst_qc2 <- foreach(i = 1:nlayers(rst_qc1), .packages = lib) %dopar% {
#     overlay(rst_qc1[[i]], rst_crp[[2]][[i]], fun = function(x, y) {
#       id <- sapply(y[], function(k) {
#         bin <- number2binary(k, 8, TRUE)
#         snow <- substr(bin, 6, 6) == "0"
#         cirrus <- substr(bin, 4, 4) == "0"
#         intern_cloud <- substr(bin, 3, 3) == "0"
#         cloud_shadow <- substr(bin, 2, 2) == "0"
# 
#         all(snow, cirrus, intern_cloud, cloud_shadow)
#       })
# 
#       x[!id] <- NA
#       return(x)
#     }, filename = paste(dir_qc2, names(rst_qc1[[i]]), sep = "/"),
#     overwrite = TRUE, format = "GTiff")
# }
# 
# rst_qc2 <- raster::stack(lst_qc2)

## reimport step #2 quality-controlled files
fls_qc2 <- list.files(dir_qc2, pattern = "Lai_500m.tif$", full.names = TRUE)
rst_qc2 <- raster::stack(fls_qc2)


### quality control, step #3: --------------------------------------------------
### discard pixels adjacent to clouds

dir_adj <- "data/MCD15A3H.006/adj"
if (!dir.exists(dir_adj)) dir.create(dir_adj)

# lst_adj <- foreach(i = unstack(rst_qc2), .packages = "raster") %dopar% {
#   msk <- raster::focal(i, w = matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 1), 3, 3), 
#                        fun = function(...) any(is.na(...)))
#   
#   raster::overlay(i, msk, fun = function(x, y) {
#     x[y[] == 1] <- NA
#     return(x)
#   }, filename = paste(dir_adj, names(i), sep = "/"), 
#   format = "GTiff", overwrite = TRUE)
# }
# 
# rst_adj <- stack(lst_adj)

## reimport step #3 quality-controlled files
fls_adj <- list.files(dir_adj, pattern = "Lai_500m.tif$", full.names = TRUE)
rst_adj <- stack(fls_adj)


### gap-filling ----------------------------------------------------------------
### kolmogorov-zurbenko adaptive

mat_adj <- as.matrix(rst_adj)

mat_fll <- foreach(i = 1:nrow(mat_adj), .combine = "rbind") %do% {
  val <- mat_adj[i, ]
  
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

rst_fll <- setValues(rst_adj, mat_fll)

dir_gf <- "data/MCD15A3H.006/gf"
if (!dir.exists(dir_gf)) dir.create(dir_gf)

rst_fll <- stack(foreach(i = 1:nlayers(rst_fll), .packages = lib) %dopar% {
  writeRaster(rst_fll[[i]], 
              filename = paste0(dir_gf, "/", names(rst_fll[[i]]), ".tif"), 
              format = "GTiff", overwrite = TRUE)
})


### replicate missing dates -----

dts <- seq(as.Date("2013-01-01"), as.Date("2015-12-31"), 1)
dts <- strftime(dts, format = "%Y%j")

dts_gf <- extractDate(names(rst_fll), 11, 17)$inputLayerDates

dat_rpl <- merge(data.frame(date = dts), 
                 data.frame(date = dts_gf, file = names(rst_fll)), 
                 by = "date", all.x = TRUE)

## replicate available files
for (i in 1:nrow(dat_rpl)) {
  if (is.na(dat_rpl$file[i])) {
    dat_rpl$file[i] <- dat_rpl$file[i-1]
  }
}

## target folder and files
dir_rpl <- "data/MCD15A3H.006/rpl"
if (!dir.exists(dir_rpl)) dir.create(dir_rpl)

fls_rpl <- paste0(dir_gf, "/", dat_rpl$file, ".tif")
rst_rpl <- stack(fls_rpl)

fls_rpl_out <- sapply(1:nrow(dat_rpl), function(i) {
  gsub(extractDate(as.character(dat_rpl$file[i]), 11, 17)$inputLayerDates, 
       as.character(dat_rpl$date[i]), as.character(dat_rpl$file[i]))
})
fls_rpl_out <- paste0(dir_rpl, "/", fls_rpl_out, ".tif")

## write to files
rst_rpl <- stack(foreach(i = 1:nlayers(rst_rpl), j = fls_rpl_out, 
                   .packages = "raster") %dopar% {
  writeRaster(rst_rpl[[i]], filename = j, format = "GTiff", overwrite = TRUE)
})

## close parallel backend
stopCluster(cl)
