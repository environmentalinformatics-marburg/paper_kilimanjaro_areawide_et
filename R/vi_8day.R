### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## set working directory
Orcs::setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/",
              path_ext = "kilimanjaro/evapotranspiration")

## load packages
lib <- c("MODIS", "Rsenal", "rgdal", "doParallel")
Orcs::loadPkgs(lib)

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


### crop layers ----------------------------------------------------------------

## import plot coordinates
suppressWarnings(
  spt_plots <- readOGR("data/shp/", "PlotPoles_ARC1960_mod_20140807_final",
                       p4s = "+init=epsg:21037")
)

chr_plots <- c("sav5", "sav0", "mai4", "mai0", # colline zone
               "cof3", "cof2", "gra3", "gra2", # submontane zone
               "fed1", "hel1", "fer0")         # subalpine and lower alpine zone

spt_plots <- subset(spt_plots, PlotID %in% chr_plots & PoleType == "AMP")

## extend bounding box
ext_plots <- extent(spt_plots)

num_xmin <- xmin(ext_plots) - 1e3
num_xmax <- xmax(ext_plots) + 1e3
num_ymin <- ymin(ext_plots) - 1e3
num_ymax <- ymax(ext_plots) + 1e3

ext_crop <- extent(c(num_xmin, num_xmax, num_ymin, num_ymax))

# ## perform crop
# rst_crp <- foreach(i = c("b01", "b02", "state_250m", "qc_250m"),
#                    .packages = "MODIS") %dopar% {
#
#     # list and import available files
#     fls <- list.files(paste0(getOption("MODIS_outDirPath"), "/MYD09Q1.006"),
#                       pattern = paste0(i, ".tif$"), full.names = TRUE)
#     rst <- raster::stack(fls)
#
#     # crop
#     dir_out <- "data/MYD09Q1.006/crp"
#     if (!dir.exists(dir_out)) dir.create(dir_out)
#
#     fls_out <- paste(dir_out, basename(fls), sep = "/")
#     rst_out <- raster::crop(rst, ext_crop, snap = "out")
#
#     # if dealing with band 1 or band 2, apply scale factor
#     if (i %in% c("b01", "b02")) {
#       raster::dataType(rst_out) <- "INT2S"
#       rst_out <- rst_out * 0.0001
#
#     # else convert to 16-bit unsigned integer
#     } else {
#       raster::dataType(rst_out) <- "INT2U"
#     }
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
rst_crp <- foreach(i = c("b01", "b02", "state_250m", "qc_250m")) %do% {

  # list and import available files
  dir_crp <- "data/MYD09Q1.006/crp"
  fls_crp <- list.files(dir_crp, pattern = paste0(i, ".tif$"), full.names = TRUE)
  raster::stack(fls_crp)
}


### quality control, step #1: --------------------------------------------------
### discard cloudy pixels based on companion quality information ('qc_250m')

dir_qc1 <- "data/MYD09Q1.006/qc1"
if (!dir.exists(dir_qc1)) dir.create(dir_qc1)

# ## perform quality check #1 for each band separately
# lst_qc1 <- foreach(i = rst_crp[1:2]) %do% {
#
#   ## loop over layers
#   lst_out <- foreach(j = 1:nlayers(i), .packages = lib) %dopar%
#     overlay(i[[j]], rst_crp[[4]][[j]], fun = function(x, y) {
#       id <- sapply(y[], function(k) {
#         bin <- number2binary(k, 16, TRUE)
#         modland <- substr(bin, 15, 16) == "00"
#         qual_b1 <- substr(bin, 9, 12) == "0000"
#         qual_b2 <- substr(bin, 5, 8) == "0000"
#
#         all(modland, qual_b1, qual_b2)
#       })
#
#       x[!id] <- NA
#       return(x)
#     }, filename = paste(dir_qc1, names(i[[j]]), sep = "/"),
#     overwrite = TRUE, format = "GTiff")
#
#   raster::stack(lst_out)
# }

## reimport step #1 quality-controlled files
lst_qc1 <- foreach(i = c("b01", "b02")) %do% {
  fls_qc1 <- list.files(dir_qc1, pattern = paste0(i, ".tif$"), full.names = TRUE)
  raster::stack(fls_qc1)
}


### quality control, step #2: --------------------------------------------------
### discard cloudy pixels based on 'state_250m' flags

dir_qc2 <- "data/MYD09Q1.006/qc2"
if (!dir.exists(dir_qc2)) dir.create(dir_qc2)

# ## perform quality check #2 for each band separately
# lst_qc2 <- foreach(i = lst_qc1) %do% {
#
#   ## loop over layers
#   lst_out <- foreach(j = 1:nlayers(i), .packages = lib) %dopar%
#     overlay(i[[j]], rst_crp[[3]][[j]], fun = function(x, y) {
#       id <- sapply(y[], function(k) {
#         bin <- number2binary(k, 16, TRUE)
#         cloud_state <- substr(bin, 15, 16) == "00"
#         cloud_shadow <- substr(bin, 14, 14) == "0"
#         cirrus <- substr(bin, 7, 8) == "00"
#         intern_cloud <- substr(bin, 6, 6) == "0"
#         fire <- substr(bin, 5, 5) == "0"
#         snow <- substr(bin, 4, 4) == "0"
#         cloud_adj <- substr(bin, 3, 3) == "0"
#         intern_snow <- substr(bin, 1, 1) == "0"
#
#         all(cloud_state, cloud_shadow, cirrus, intern_cloud,
#             fire, snow, cloud_adj, intern_snow)
#       })
#
#       x[!id] <- NA
#       return(x)
#     }, filename = paste(dir_qc2, names(i[[j]]), sep = "/"),
#     overwrite = TRUE, format = "GTiff")
#
#   raster::stack(lst_out)
# }

## reimport step #2 quality-controlled files
lst_qc2 <- foreach(i = c("b01", "b02")) %do% {
  fls_qc2 <- list.files(dir_qc2, pattern = paste0(i, ".tif$"), full.names = TRUE)
  raster::stack(fls_qc2)
}


# ### monthly value composites ---------------------------------------------------
#
# indices <- as.numeric(as.factor(as.yearmon(seq(as.Date("2013-01-01"),
#                                                as.Date("2013-02-28"), 8))))
#
# rst.b1.b2.agg <-
#   foreach(i = c(1, 2), .packages = lib) %dopar% {
#     foreach(j = unique(indices), .combine = "stack") %do% {
#       sub <- rst.b1.b2.cc[[i]][[grep(j, indices)]]
#       calc(sub, fun = function(x) {
#         if (all(is.na(x))) return(NA) else return(round(median(x, na.rm = TRUE)))
#       }, filename = paste0("myd09q1/processed/AGG_", names(sub)[1]),
#       format = "GTiff", overwrite = TRUE)
#     }
#   }


### normalized difference vegetation index -------------------------------------

dir_ndvi <- "data/MYD09Q1.006/ndvi"
if (!dir.exists(dir_ndvi)) dir.create(dir_ndvi)

fls_ndvi <- gsub("_b01", "", names(lst_qc1[[1]]))

rst_ndvi <- ndvi(lst_qc2[[1]], lst_qc2[[2]],
                 filename = paste0(dir_ndvi, "/NDVI"),
                 format = "GTiff", overwrite = TRUE,
                 bylayer = TRUE, suffix = fls_ndvi)


### soil-adjusted vegetation index ---------------------------------------------

dir_savi <- "data/MYD09Q1.006/savi"
if (!dir.exists(dir_savi)) dir.create(dir_savi)

fls_savi <- gsub("_b01", "", names(lst_qc1[[1]]))

rst_savi <- savi(lst_qc2[[1]], lst_qc2[[2]],
                 filename = paste0(dir_savi, "/SAVI"),
                 format = "GTiff", overwrite = TRUE,
                 bylayer = TRUE, suffix = fls_savi)

## deregister parallel backend
stopCluster(cl)
