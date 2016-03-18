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

## combinded terra/aqua product
runGdal(product = "MCD15A2H",
        collection = getCollection("MCD15A2H", forceCheck = TRUE),
        begin = "2013001", end = "2015365", tileH = 21, tileV = 9,
        SDSstring = "011100", job = "MCD15A2H.006")


### crop layers ----------------------------------------------------------------

## retrieve study extent
ext_crp <- uniformExtent(verbose = FALSE)

## perform crop
rst_crp <- foreach(i = c("Lai_500m", "Lai_QC", "Extra_QC"),
                   .packages = "MODIS") %dopar% {

    # list and import available files
    fls <- list.files(paste0(getOption("MODIS_outDirPath"), "/MCD15A2H.006"),
                      pattern = paste0(i, ".tif$"), full.names = TRUE)
    rst <- raster::stack(fls)

    # crop
    dir_out <- "data/MCD15A2H.006/crp"
    if (!dir.exists(dir_out)) dir.create(dir_out)

    fls_out <- paste(dir_out, basename(fls), sep = "/")
    rst_out <- raster::crop(rst, ext_crp, snap = "out")

    # convert to 8-bit unsigned integer
    raster::dataType(rst_out) <- "INT1U"

    # apply scale factor
    if (i == "Lai_500m")
      rst_out <- rst_out * 0.1

    # save and return cropped layers
    lst_out <- lapply(1:nlayers(rst_out), function(j)
      raster::writeRaster(rst_out[[j]], filename = fls_out[j],
                          format = "GTiff", overwrite = TRUE)
    )

    raster::stack(lst_out)
  }

## reimport cropped files
rst_crp <- foreach(i = c("Lai_500m", "Lai_QC", "Extra_QC")) %do% {

  # list and import available files
  dir_crp <- "data/MCD15A2H.006/crp"
  fls_crp <- list.files(dir_crp, pattern = paste0(i, ".tif$"), full.names = TRUE)
  raster::stack(fls_crp)
}


### quality control, step #1: --------------------------------------------------
### discard cloudy pixels based on companion quality information ('Lai_QC')

dir_qc1 <- "data/MCD15A2H.006/qc1"
if (!dir.exists(dir_qc1)) dir.create(dir_qc1)

## perform quality check #1
lst_qc1 <- foreach(i = 1:nlayers(rst_crp[[1]]), .packages = lib) %dopar%
    overlay(rst_crp[[1]][[i]], rst_crp[[2]][[i]], fun = function(x, y) {
      id <- sapply(y[], function(k) {
        bin <- number2binary(k, 8, TRUE)
        modland <- substr(bin, 8, 8) == "0"
        cloud_state <- substr(bin, 4, 5) == "00"
        confidence <- substr(bin, 1, 3) %in% c("000", "001")

        all(modland, cloud_state, confidence)
      })

      x[!id] <- NA
      return(x)
    }, filename = paste(dir_qc1, names(rst_crp[[1]][[i]]), sep = "/"),
    overwrite = TRUE, format = "GTiff")

rst_qc1 <- stack(lst_qc1)

## reimport step #1 quality-controlled files
fls_qc1 <- list.files(dir_qc1, pattern = "Lai_500m.tif$", full.names = TRUE)
rst_qc1 <- raster::stack(fls_qc1)


### quality control, step #2: --------------------------------------------------
### discard cloudy pixels based on 'state_250m' flags

dir_qc2 <- "data/MCD15A2H.006/qc2"
if (!dir.exists(dir_qc2)) dir.create(dir_qc2)

## perform quality check #2
lst_qc2 <- foreach(i = 1:nlayers(rst_qc1), .packages = lib) %dopar%
    overlay(rst_qc1[[i]], rst_crp[[2]][[i]], fun = function(x, y) {
      id <- sapply(y[], function(k) {
        bin <- number2binary(k, 8, TRUE)
        snow <- substr(bin, 6, 6) == "0"
        cirrus <- substr(bin, 4, 4) == "0"
        intern_cloud <- substr(bin, 3, 3) == "0"
        cloud_shadow <- substr(bin, 2, 2) == "0"

        all(snow, cirrus, intern_cloud, cloud_shadow)
      })

      x[!id] <- NA
      return(x)
    }, filename = paste(dir_qc2, names(rst_qc1[[i]]), sep = "/"),
    overwrite = TRUE, format = "GTiff")

rst_qc2 <- raster::stack(lst_qc2)

## reimport step #2 quality-controlled files
fls_qc2 <- list.files(dir_qc2, pattern = "Lai_500m.tif$", full.names = TRUE)
rst_qc2 <- raster::stack(fls_qc2)

