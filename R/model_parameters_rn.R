### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("reset", "satellite", "doParallel")
Orcs::loadPkgs(lib)

## set working directory
Orcs::setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/",
              path_ext = "kilimanjaro/evapotranspiration/")

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

## reference extent
rst_ref_utm <- raster("data/MYD09Q1.006/ndvi/NDVI_MYD09Q1.A2013001.sur_refl.tif")
rst_ref <- trim(projectRaster(rst_ref_utm, crs = "+init=epsg:4326"))


### atmospheric transmissivity -----

## solar zenith angle
fls_theta <- list.files("data/MYD05_L2.006/Solar_Zenith", full.names = TRUE, 
                        pattern = "Solar_Zenith.tif$")

dts_theta <- substr(basename(fls_theta), 11, 22)
dat_theta <- data.frame(datetime = dts_theta, theta = fls_theta, 
                        stringsAsFactors = FALSE)

## precipitable water content
fls_ptw <- list.files("data/MYD05_L2.006/Water_Vapor_Near_Infrared/res", 
                      full.names = TRUE, pattern = "Near_Infrared.tif$")

# # merge replicate daily scenes
# dts_ptw <- MODIS::extractDate(fls_ptw, 23, 29, asDate = TRUE)$inputLayerDates
# fct_ptw <- as.numeric(as.factor(dts_ptw))
# rst_ptw <- stackApply(rst_ptw, indices = fct_ptw, fun = mean)
# 
# dir_dly <- paste0(unique(dirname(fls_ptw)), "/daily/")
# if (!dir.exists(dir_dly)) dir.create(dir_dly)
# fls_dly <- paste0(dir_dly, "Surface_Pressure_", strftime(dts_ptw, "%Y%j"), ".tif")
# lst_ptw <- foreach(i = 1:nlayers(rst_ptw), .packages = "raster") %dopar% {
#   if (file.exists(fls_dly[i])) {
#     raster(fls_dly[i])
#   } else {
#     writeRaster(rst_ptw[[i]], filename = fls_dly[i], 
#                 format = "GTiff", overwrite = TRUE)
#   }
# }
# rst_ptw <- stack(lst_ptw)

dts_ptw <- substr(basename(fls_ptw), 23, 34)
dat_ptw <- data.frame(datetime = dts_ptw, precipitable = fls_ptw, 
                      stringsAsFactors = FALSE)

## air pressure (modis)

# import available files
fls_sp <- list.files("data/MYD07_L2.006/Surface_Pressure", full.names = TRUE, 
                     pattern = "Surface_Pressure.tif$")

dts_sp <- substr(basename(fls_sp), 11, 22)
dat_sp <- data.frame(datetime = dts_sp, pressure = fls_sp, 
                     stringsAsFactors = FALSE)

## air pressure (ecmwf)

# import available files
fls_sp_ecmwf <- list.files("data/ecmwf/surface_pressure", full.names = TRUE, 
                           pattern = "surface_pressure.tif$")

dts_sp_ecmwf <- substr(basename(fls_sp_ecmwf), 1, 7)
dat_sp_ecmwf <- data.frame(date = dts_sp_ecmwf, pressure = fls_sp_ecmwf, 
                           stringsAsFactors = FALSE)

## compute atmospheric transmissivity for each complete scene
dat_mrg <- merge(dat_theta, dat_ptw, all = TRUE, by = "datetime")
dat_mrg$date <- substr(dat_mrg$datetime, 1, 7)
dat_tau <- merge(dat_mrg, dat_sp_ecmwf, all = TRUE, by = "date")

# dat_tau <- Reduce(function(...) merge(..., all = TRUE, by = "datetime"), 
#                   list(dat_sp, dat_theta, dat_ptw))

dat_tau <- dat_tau[complete.cases(dat_tau), ]

fls_tau <- gsub("Solar_Zenith", "Atmospheric_Transmissivity", dat_tau$theta)
dir_tau <- unique(dirname(fls_tau))
if (!dir.exists(dir_tau)) dir.create(dir_tau)
lst_tau <- foreach(i = 1:nrow(dat_tau), .packages = "reset") %dopar% {
  if (file.exists(fls_tau[i])) {
    raster::raster(fls_tau[i])
  } else {
    rst <- reset::atmosTrans(Pa = raster::raster(dat_tau$pressure[i]),
                             theta = raster::raster(dat_tau$theta[i]) * pi / 180,
                             W = raster::raster(dat_tau$precipitable[i]))

    if (!(all(is.na(raster::maxValue(rst))))) {
      raster::writeRaster(rst, filename = fls_tau[i], format = "GTiff",
                          overwrite = TRUE)
    } else {
      return(NULL)
    }
  }
}

fls_tau <- list.files(dir_tau, pattern = "Transmissivity.tif", full.names = TRUE)
rst_tau <- stack(fls_tau)
dat_tau <- data.frame(datetime = substr(basename(fls_tau), 11, 22), 
                      transmissivity = fls_tau, stringsAsFactors = FALSE)


### inverse squared earth-sun distance ---

esd <- calcEarthSunDist(strptime(dts_theta, "%Y%j.%H%S"), formula = "Duffie")
dat_esd <- data.frame(datetime = dts_theta, distance = esd, 
                      stringsAsFactors = FALSE)


### solar incidence angle -----

## digital elevation model in epsg:4326
rst_dem <- raster("data/dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- trim(projectRaster(resample(rst_dem, rst_ref_utm), crs = "+init=epsg:4326"))

## scan start time
fls_sst <- list.files("data/MYD05_L2.006/Scan_Start_Time", full.names = TRUE, 
                      pattern = "Scan_Start_Time.tif$")

dts_sst <- substr(basename(fls_sst), 11, 22)
dat_sst <- data.frame(datetime = dts_sst, start_time = fls_sst, 
                      stringsAsFactors = FALSE)

## time from solar zenith
rst <- raster(dat_theta$theta[1])

yrs <- raster::res(rst)[2]
lts <- seq(raster::ymax(rst) - .5 * yrs, raster::ymin(rst), -yrs)

mat <- matrix(ncol = raster::ncol(rst), nrow = raster::nrow(rst))
for (i in 1:nrow(mat)) mat[i, ] <- lts[i]

phi <- raster::raster(mat, template = rst)

fls_opt <- gsub("Scan_Start_Time", "Overpass_Time", fls_sst)
dir_opt <- unique(dirname(fls_opt))
if (!dir.exists(dir_opt)) dir.create(dir_opt)
dts_theta <- as.Date(dat_theta[, 1], "%Y%j.%H%M")

lst_opt <- foreach(i = 1:length(fls_theta), .packages = "satellite") %dopar% {
  if (file.exists(fls_opt[i])) {
    raster::raster(fls_opt[i])
  } else {
    satellite::fromSolarZenith(raster::stack(fls_theta[i]), phi, n = dts_theta[i], 
                               formula = "Spencer", daytime = "PM", 
                               filename = fls_opt[i], format = "GTiff", overwrite = TRUE)
  }
}

## solar incidence angle
dir_sia <- "data/MYD05_L2.006/Solar_Incidence_Angle"
if (!dir.exists(dir_sia)) dir.create(dir_sia)

lst_sia <- foreach(i = lst_opt, .packages = "satellite") %dopar% {
  fls <- gsub("Overpass_Time", "Solar_Incidence_Angle", names(i))
  fls <- paste0(dir_sia, "/", fls, ".tif")
  if (file.exists(fls)) {
    raster(fls)
  } else {
    rst <- satellite::solarIncidenceAngle(i, dem = rst_dem, formula = "Spencer")
    raster::writeRaster(rst, filename = fls, format = "GTiff", overwrite = TRUE)
  }
}

fls_sia <- list.files("data/MYD05_L2.006/Solar_Incidence_Angle", 
                      full.names = TRUE, pattern = "Incidence_Angle.tif")

dts_sia <- substr(basename(fls_sia), 11, 22)
dat_sia <- data.frame(datetime = dts_sia, incidence = fls_sia, 
                      stringsAsFactors = FALSE)


### shortwave downward radiation -----

dat_rsd <- Reduce(function(...) merge(..., all = TRUE, by = "datetime"), 
                  list(dat_esd, dat_sia, dat_tau))
dat_rsd <- dat_rsd[complete.cases(dat_rsd), ]

## target folder and files
fls_rsd <- gsub("Solar_Incidence_Angle", "Shortwave_Downward_Radiation", 
                dat_rsd$incidence)
dir_rsd <- unique(dirname(fls_rsd))
if (!dir.exists(dir_rsd)) dir.create(dir_rsd)

## compute per-scene shortwave downward radiation
lst_rsd <- foreach(i = 1:nrow(dat_rsd), .packages = "reset") %dopar% {
  if (file.exists(fls_rsd[i])) {
    raster::raster(fls_rsd[i])
  } else {
    rst <- reset::swdr(d = dat_rsd$distance[i], 
                       tau = raster::raster(dat_rsd$transmissivity[i]), 
                       theta = raster::raster(dat_rsd$incidence[i]))
  
    raster::writeRaster(rst, filename = fls_rsd[i], 
                        format = "GTiff", overwrite = TRUE)  
  }
}
rst_rsd <- stack(lst_rsd)

## 8-day aggregation
fls_rsd <- list.files("data/MYD05_L2.006/Shortwave_Downward_Radiation", 
                      full.names = TRUE, pattern = "Shortwave_Downward_Radiation.tif")

dts_rsd <- substr(basename(fls_rsd), 11, 22)
dat_rsd <- data.frame(datetime = dts_rsd, swdr = fls_rsd, 
                      stringsAsFactors = FALSE)

dir_agg <- "data/MYD05_L2.006/Shortwave_Downward_Radiation/agg"
if (!dir.exists(dir_agg)) dir.create(dir_agg)

lst_rsd_agg <- lapply(2013:2015, function(i) {
  dat <- dat_rsd[grep(i, dts_rsd), ]
  id <- as.numeric(cut(as.Date(dat$datetime, "%Y%j.%H%M"), "8 days"))
  rst <- raster::stackApply(raster::stack(dat$swdr), id, fun = mean)
  
  fls_agg <- paste0(dir_agg, "/Shortwave_Downward_Radiation_", 
                    unique(strftime(cut(as.Date(dat$datetime, "%Y%j.%H%M"), "8 days"), 
                             "%Y%j")))
  
  lst_agg <- foreach(j = 1:nlayers(rst)) %do% 
    writeRaster(rst[[j]], filename = fls_agg[j], 
                format = "GTiff", overwrite = TRUE)
  
  raster::stack(lst_agg)
})

rst_rsd_agg <- stack(lst_rsd_agg)

fls_rsd_agg <- list.files("data/MYD05_L2.006/Shortwave_Downward_Radiation/agg", 
                          pattern = "^Shortwave.*.tif$", full.names = TRUE)

dts_rsd_agg <- substr(sapply(strsplit(basename(fls_rsd_agg), "_"), "[[", 4), 1, 7)
dat_rsd_agg <- data.frame(datetime = dts_rsd_agg, swdr_agg = fls_rsd_agg, 
                          stringsAsFactors = FALSE)


### shortwave upward radiation -----

## shortwave white-sky albedo
fls_alpha <- list.files("data/MCD43A3.005/crp", full.names = TRUE, 
                        pattern = "^MCD43A3.*shortwave.tif$")
rst_alpha <- stack(fls_alpha)
rst_alpha <- resample(rst_alpha, rst_ref_utm)

dts_alpha <- substr(basename(fls_alpha), 10, 16)
dat_alpha <- data.frame(datetime = dts_alpha, albedo = fls_alpha, 
                        stringsAsFactors = FALSE)

## shortwave upward radiation
dat_rsu <- merge(dat_alpha, dat_rsd_agg, by = "datetime", all = TRUE)
dat_rsu <- dat_rsu[complete.cases(dat_rsu), ]

dir_rsu <- "data/MCD43A3.005/Shortwave_Upward_Radiation"
if (!dir.exists(dir_rsu)) dir.create(dir_rsu)

fls_rsu <- gsub("Albedo_WSA_shortwave", "Shortwave_Upward_Radiation", dat_rsu$albedo)
fls_rsu <- paste0(dir_rsu, "/", basename(fls_rsu))

lst_rsu <- foreach(i = 1:nrow(dat_rsu), .packages = "raster") %dopar% {
  swdr <- raster::resample(
    raster::projectRaster(
      raster::raster(dat_rsu$swdr[i]), crs = "+init=epsg:21037")
  , rst_ref_utm)
  
  albedo <- raster::resample(
    raster::raster(dat_rsu$albedo[i])
  , rst_ref_utm)
  
  raster::overlay(swdr, albedo, fun = function(x, y) x * (1 - y), 
                  filename = fls_rsu[i], format = "GTiff", overwrite = TRUE)
}

rst_rsu <- stack(lst_rsu)


### land surface temperature -----

fls_lst <- list.files("data/MYD11A2.005/qc", full.names = TRUE, 
                        pattern = "^MYD11A2.*Day_1km.tif$")
rst_lst <- stack(fls_lst)
rst_lst <- resample(rst_lst, rst_ref_utm)

dts_lst <- substr(basename(fls_lst), 10, 16)
dat_lst <- data.frame(datetime = dts_lst, lst = fls_lst, 
                      stringsAsFactors = FALSE)


### broadband surface emissivity -----

## ndvi
fls_ndvi <- list.files("data/MYD09Q1.006/ndvi", full.names = TRUE, 
                      pattern = "^NDVI_MYD09Q1.*.tif$")
rst_ndvi <- stack(fls_ndvi)

dts_ndvi <- substr(basename(fls_ndvi), 15, 21)
dat_ndvi <- data.frame(datetime = dts_ndvi, ndvi = fls_ndvi, 
                      stringsAsFactors = FALSE)

## lai
fls_lai <- list.files("data/MCD15A2H.006/qc2", full.names = TRUE, 
                       pattern = "^MCD15A2H.*.tif$")
rst_lai <- stack(fls_lai)
rst_lai <- resample(rst_lai, rst_ndvi[[1]])

dts_lai <- substr(basename(fls_lai), 11, 17)
dat_lai <- data.frame(datetime = dts_lai, lai = fls_lai, 
                       stringsAsFactors = FALSE)

## surface emissivity
dir_bsa <- "data/MCD15A2H.006/Surface_Emissivity"
if (!dir.exists(dir_bsa)) dir.create(dir_bsa)

fls_bsa <- gsub("Lai", "Surface_Emissivity", names(rst_lai))
fls_bsa <- paste0(dir_bsa, "/", fls_bsa)

lst_bsa <- foreach(i = 1:nlayers(rst_ndvi), j = fls_bsa, 
                   .packages = "reset") %dopar% {
  if (file.exists(j)) {
    raster::raster(j)
  } else {
    reset::surfaceEmissivity(rst_ndvi[[i]], rst_lai[[i]], filename = j, 
                           format = "GTiff", overwrite = TRUE)
  }
}

rst_bsa <- stack(lst_bsa)


### longwave upward radiation -----

dir_lwur <- "data/MCD15A2H.006/Longwave_Upward_Radiation"
if (!dir.exists(dir_lwur)) dir.create(dir_lwur)

fls_lwur <- gsub("Surface_Emissivity", "Longwave_Upward_Radiation", names(rst_bsa))
fls_lwur <- paste0(dir_lwur, "/", fls_lwur)

lst_lwur <- foreach(i = 1:nlayers(rst_lst), j = fls_lwur, 
                    .packages = "reset") %dopar% {
  if (file.exists(j)) {
    raster::raster(j)
  } else {
    reset::lwur(lst = rst_lst[[i]], bse = rst_bsa[[i]], filename = j, 
                format = "GTiff", overwrite = TRUE)
  }
}

rst_lwur <- stack(lst_lwur)


## R_ld (longwave downward radiation)

# ta (air temperature)
ta.fls <- list.files("md07/processed", 
                     pattern = "AGGMLY_.*Retrieved_Temperature_Profile.tif$", 
                     full.names = TRUE)
ta <- stack(ta.fls)

# dpt (dewpoint temperature)
dpt.fls <- list.files("md07/processed", 
                      pattern = "^AGGMLY_.*Retrieved_Moisture_Profile.tif$", 
                      full.names = TRUE)
dpt <- stack(dpt.fls)

# e_a (actual vapor pressure, taken from Ryu 2008)
e.a <- 6.11 * exp(19.59 * (dpt - 273.15) / dpt)

# # (taken from Bisht 2005)
# l.v <- 2.5 * 10^6
# r.v <- 461
# t0 <- 273.15
# e.a2 <- 6.11 * exp(l.v / r.v * (1 / t0 - 1 / dpt))

r.ld <- 1.24 * (e.a / ta)^0.14 * sigma * ta^4


## R_n (net radiation)

r.n <- r.sd

stopCluster(cl)