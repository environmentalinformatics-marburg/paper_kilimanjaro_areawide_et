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
fls_theta <- list.files("data/MOD05_L2.006/Solar_Zenith", full.names = TRUE, 
                        pattern = "Solar_Zenith.tif$")

dts_theta <- substr(basename(fls_theta), 11, 22)
dat_theta <- data.frame(datetime = dts_theta, theta = fls_theta, 
                        stringsAsFactors = FALSE)

## precipitable water content
fls_ptw <- list.files("data/MOD05_L2.006/Water_Vapor_Near_Infrared/res", 
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
fls_sp <- list.files("data/MOD07_L2.006/Surface_Pressure", full.names = TRUE, 
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
    
    raster::writeRaster(rst, filename = fls_tau[i], format = "GTiff", 
                        overwrite = TRUE)
  }
}
rst_tau <- stack(lst_tau)

dat_tau <- data.frame(datetime = dat_tau$datetime, transmissivity = fls_tau, 
                      stringsAsFactors = FALSE)


### inverse squared earth-sun distance ---

esd <- calcEarthSunDist(strptime(dts_theta, "%Y%j.%H%S"), formula = "Duffie")
dat_esd <- data.frame(datetime = dts_theta, distance = esd, 
                      stringsAsFactors = FALSE)


### solar incidence angle -----

## digital elevation model in epsg:4326
rst_dem <- raster("data/dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- trim(projectRaster(resample(rst_dem, rst_ref_utm), crs = "+init=epsg:4326"))

## scan start time
fls_sst <- list.files("data/MOD05_L2.006/Scan_Start_Time", full.names = TRUE, 
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

foreach(i = 1:length(fls_theta), .packages = "satellite") %dopar%
  fromSolarZenith(raster::stack(fls_theta[1:3]), phi, n = dts_theta[1:3], 
                  formula = "Spencer", daytime = "AM")

, 
filename = fls_opt[1], format = "GTiff", overwrite = TRUE



## solar incidence angle
dir_sia <- "data/MOD05_L2.006/Solar_Incidence_Angle"
if (!dir.exists(dir_sia)) dir.create(dir_sia)

lst_sia <- foreach(i = unstack(rst_sst), .packages = "satellite") %dopar% {
  fls <- gsub("Scan_Start_Time", "Solar_Incidence_Angle", names(i))
  fls <- paste0(dir_sia, "/", fls, ".tif")
  if (file.exists(fls)) {
    raster(fls)
  } else {
    rst <- solarIncidenceAngle(i, dem = rst_dem, formula = "Spencer")
    writeRaster(rst, filename = fls, format = "GTiff", overwrite = TRUE)
  }
}

fls_sia <- list.files("data/MOD05_L2.006/Solar_Incidence_Angle", 
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


## R_lu (longwave upward radiation)

st <- "2012"
nd <- "2012"

# NDVI
ndvi.fls.terra <- list.files("data/MOD13Q1.006/whittaker", pattern = ".tif$", 
                       full.names = TRUE)

ndvi.fls.terra <- ndvi.fls.terra[grep(st, ndvi.fls.terra)[1]:grep(nd, ndvi.fls.terra)[length(grep(nd, ndvi.fls.terra))]]

ndvi.fls.aqua <- list.files("data/MYD13Q1.006/whittaker", pattern = ".tif$", 
                             full.names = TRUE)

ndvi.fls.aqua <- ndvi.fls.aqua[grep(st, ndvi.fls.aqua)[1]:grep(nd, ndvi.fls.aqua)[length(grep(nd, ndvi.fls.aqua))]]

ndvi.fls <- c(ndvi.fls.terra, ndvi.fls.aqua)
ndvi.fls <- ndvi.fls[order(substr(basename(ndvi.fls), 5, 11))]

ndvi <- stack(ndvi.fls)

ndvi_rsmpl <- resample(ndvi, lai[[1]], filename = "data/MCD13Q1.006/rsmpl/RSMPL", 
                       bylayer = TRUE, suffix = names(ndvi), format = "GTiff", 
                       overwrite = TRUE)

# LAI
lai.fls <- list.files("data/lai_rsmpl", pattern = "^COMB.*.tif", 
                      full.names = TRUE)
lai.fls <- lai.fls[grep(st, lai.fls)[1]:grep(nd, lai.fls)[length(grep(nd, lai.fls))]]
lai <- stack(lai.fls)

# sigma (Stefan-Boltzmann constant, 5.67 * 10^(-8) Wm^(-2)K^(-4))
sigma <- 5.67 * 10^(-8)

# lst (land surface temperature)
lst.fls <- list.files("data/MCD11A2/rsmpl/", pattern = "^RSMPL.*.tif$", 
                      full.names = TRUE)
lst.fls <- lst.fls[grep(st, lst.fls)[1]:grep(nd, lst.fls)[length(grep(nd, lst.fls))]]
lst <- stack(lst.fls)

# epsilon_s (broadband surface emissivity)
epsilon.s <- stack(lapply(seq_len(nlayers(ndvi_rsmpl)), function(i) {
  x <- ndvi_rsmpl[[i]] / 10000
  y <- lai[[i]]

  x[x[] > 0 & y[] <= 3] <- 0.95 + 0.01 * y[x[] > 0 & y[] <= 3]
  x[x[] > 0 & y[] > 3] <- 0.98
  x[x[] <= 0] <- 0.985
  
  return(x)
}))

# R_lu
r.lu <- sigma * epsilon.s * lst^4
r.lu <- writeRaster(r.lu, filename = "model_input/R_lu_13", format = "GTiff", 
                    overwrite = TRUE, bylayer = FALSE)


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