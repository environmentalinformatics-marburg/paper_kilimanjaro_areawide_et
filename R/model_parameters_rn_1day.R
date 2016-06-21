### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("raster", "reset", "satellite", "doParallel", "MODIS")
Orcs::loadPkgs(lib)

## set working directory
Orcs::setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "H:/",
              path_ext = "kilimanjaro/evapotranspiration/")

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

## reference extent
rst_ref_utm <- raster("data/MYD09Q1.006/ndvi/NDVI_MYD09Q1.A2013001.sur_refl.tif")
rst_ref <- trim(projectRaster(rst_ref_utm, crs = "+init=epsg:4326"))

## sensor under investigation
sensor <- "MYD"
dir_sensor <- paste0("data/radiation/", sensor)


### atmospheric transmissivity -----

## solar zenith angle
fls_theta <- list.files(paste0("data/", sensor, "05_L2.006/Solar_Zenith"), 
                        full.names = TRUE, pattern = "Solar_Zenith.tif$")

dts_theta <- substr(basename(fls_theta), 11, 22)
dat_theta <- data.frame(datetime = dts_theta, theta = fls_theta, 
                        stringsAsFactors = FALSE)

## precipitable water content (modis; cm, i.e. multiply with 10 to get mm)
fls_ptw <- list.files(paste0("data/", sensor, "05_L2.006/Water_Vapor_Near_Infrared/res"), 
                      full.names = TRUE, pattern = "Near_Infrared.tif$")

## merge duplicated daily scenes
fls_sez <- list.files(paste0("data/", sensor, "05_L2.006/Sensor_Zenith"), 
                      full.names = TRUE, pattern = "Sensor_Zenith.tif$")

dts_ptw <- extractDate(fls_ptw, 23, 29, asDate = TRUE)$inputLayerDates

dts_dpl <- which(duplicated(dts_ptw))

if (length(dts_dpl) > 0) {
  lst_ptw <- lapply(dts_dpl, function(i) {
    fls <- fls_ptw[which(dts_ptw == dts_ptw[i])]
    rst <- lapply(fls, raster)
    
    vld <- suppressWarnings(sapply(rst, maxValue))
    if (all(is.na(vld))) {
      rm(rst)
      file.remove(fls)
    } else if (sum(is.na(vld)) == 1) {
      rst <- rst[[which(!is.na(vld))]]
      file.remove(fls[which(is.na(vld))])
    } else {
      rst_sez <- stack(fls_sez[which(dts_ptw == dts_ptw[i])])
      rst <- overlay(rst[[1]], rst[[2]], rst_sez[[1]], rst_sez[[2]], 
                     fun = function(w, x, y, z) {
                       val <- sapply(1:ncell(w), function(j) {
                         if (is.na(w[j]) & !is.na(x[j])) {
                           return(x[j])
                         } else if (!is.na(w[j]) & is.na(x[j])) {
                           return(w[j])
                         } else if (!is.na(w[j]) & !is.na(x[j])) {
                           if (y[j] < z[j]) return(w[j]) else return(x[j])
                         } else {
                           return(NA)
                         }
                       })
                     }, filename = fls[which.min(maxValue(rst_sez))], 
                     format = "GTiff", overwrite = TRUE)
      file.remove(fls[which.max(maxValue(rst_sez))])
    }
    
    return(rst)
  })
}

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

# ## air pressure (modis; hPa, i.e. multiply 0.1 to get kPa)
# 
# # import available files
# fls_sp <- list.files(paste0("data/", sensor, "07_L2.006/Surface_Pressure"), 
#                      full.names = TRUE, pattern = "Surface_Pressure.tif$")
# 
# dts_sp <- substr(basename(fls_sp), 11, 22)
# dat_sp <- data.frame(datetime = dts_sp, pressure = fls_sp, 
#                      stringsAsFactors = FALSE)

## air pressure (ecmwf; hPa, i.e. multiply 0.1 to get kPa)

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
dat_tau <- dat_tau[complete.cases(dat_tau), ]

# dat_tau <- Reduce(function(...) merge(..., all = TRUE, by = "datetime"), 
#                   list(dat_sp, dat_theta, dat_ptw))

## target folder and files
fls_tau <- gsub("Solar_Zenith", "Atmospheric_Transmissivity", dat_tau$theta)
dir_tau <- paste0(dir_sensor, "/Atmospheric_Transmissivity")
if (!dir.exists(dir_tau)) dir.create(dir_tau)
fls_tau <- paste0(dir_tau, "/", basename(fls_tau))

lst_tau <- foreach(i = 1:nrow(dat_tau), .packages = lib) %dopar% {
  if (file.exists(fls_tau[i])) {
    raster(fls_tau[i])
  } else {
    rst <- atmosTrans(Pa = raster(dat_tau$pressure[i]) * 0.1,
                      theta = raster(dat_tau$theta[i]) * pi / 180,
                      W = raster(dat_tau$precipitable[i]) * 10)
    
    if (!(all(is.na(maxValue(rst))))) {
      writeRaster(rst, filename = fls_tau[i], format = "GTiff",
                  overwrite = TRUE)
    } else {
      return(NULL)
    }
  }
}

fls_tau <- list.files(dir_tau, pattern = "Transmissivity.tif", full.names = TRUE)
dts_tau <- substr(basename(fls_tau), 11, 22)
rst_tau <- stack(fls_tau); rm(lst_tau)
dat_tau <- data.frame(datetime = substr(basename(fls_tau), 11, 22), 
                      transmissivity = fls_tau, stringsAsFactors = FALSE)


### inverse squared earth-sun distance ---

esd <- calcEarthSunDist(strptime(dts_tau, "%Y%j.%H%S"), formula = "Duffie")
dat_esd <- data.frame(datetime = dts_tau, distance = esd, 
                      stringsAsFactors = FALSE)


### solar incidence angle -----

## digital elevation model in epsg:4326
rst_dem <- raster("data/dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- trim(projectRaster(resample(rst_dem, rst_ref_utm), crs = "+init=epsg:4326"))

## scan start time
fls_sst <- list.files(paste0("data/", sensor, "05_L2.006/Scan_Start_Time"), 
                      full.names = TRUE, pattern = "Scan_Start_Time.tif$")

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
dir_opt <- paste0(dir_sensor, "/Overpass_Time")
if (!dir.exists(dir_opt)) dir.create(dir_opt)
fls_opt <- paste0(dir_opt, "/", basename(fls_opt))

dts_theta <- as.Date(dat_theta[, 1], "%Y%j.%H%M")

lst_opt <- foreach(i = 1:length(fls_theta), .packages = lib) %dopar% {
  if (file.exists(fls_opt[i])) {
    raster(fls_opt[i])
  } else {
    fromSolarZenith(stack(fls_theta[i]), phi, n = dts_theta[i], 
                    formula = "Spencer", daytime = "PM", 
                    filename = fls_opt[i], format = "GTiff", overwrite = TRUE)
  }
}

## solar incidence angle
dir_sia <- paste0(dir_sensor, "/Solar_Incidence_Angle")
if (!dir.exists(dir_sia)) dir.create(dir_sia)

lst_sia <- foreach(i = lst_opt, .packages = lib) %dopar% {
  fls <- gsub("Overpass_Time", "Solar_Incidence_Angle", names(i))
  fls <- paste0(dir_sia, "/", fls, ".tif")
  if (file.exists(fls)) {
    raster(fls)
  } else {
    rst <- solarIncidenceAngle(i, dem = rst_dem, formula = "Spencer")
    writeRaster(rst, filename = fls, format = "GTiff", overwrite = TRUE)
  }
}

fls_sia <- list.files(dir_sia, full.names = TRUE, pattern = "Incidence_Angle.tif")
dts_sia <- substr(basename(fls_sia), 11, 22)
dat_sia <- data.frame(datetime = dts_sia, incidence = fls_sia, 
                      stringsAsFactors = FALSE)


### shortwave downward radiation -----

dat_rsd <- Reduce(function(...) merge(..., by = "datetime"), 
                  list(dat_esd, dat_sia, dat_tau))
dat_rsd$datetime <- sapply(strsplit(dat_rsd$datetime, "\\."), "[[", 1)

## target folder and files
fls_rsd <- gsub("Solar_Incidence_Angle", "Shortwave_Downward_Radiation", 
                dat_rsd$incidence)
dir_rsd <- paste0(dir_sensor, "/Shortwave_Downward_Radiation")
if (!dir.exists(dir_rsd)) dir.create(dir_rsd)
fls_rsd <- paste0(dir_rsd, "/", basename(fls_rsd))

## compute per-scene shortwave downward radiation
rst_rsd <- stack(foreach(i = 1:nrow(dat_rsd), .packages = lib) %dopar% {
  if (file.exists(fls_rsd[i])) {
    raster::raster(fls_rsd[i])
  } else {
    rst <- swdr(d = dat_rsd$distance[i], 
                tau = raster(dat_rsd$transmissivity[i]), 
                theta = raster(dat_rsd$incidence[i]))
  
    writeRaster(rst, filename = fls_rsd[i], 
                        format = "GTiff", overwrite = TRUE)  
  }
})

dat_rsd$swdr <- fls_rsd


### shortwave upward radiation -----

## shortwave white-sky albedo
fls_alpha <- list.files("data/MCD43A3.006/crp", full.names = TRUE, 
                        pattern = "^MCD43A3.*Albedo_WSA_shortwave.tif$")

fls_alpha_res <- gsub("/crp/", "/res/", fls_alpha)
drs_alpha_res <- unique(dirname(fls_alpha_res))
if (!dir.exists(drs_alpha_res)) dir.create(drs_alpha_res)

lst_alpha <- foreach(i = 1:length(fls_alpha), .packages = lib) %dopar% {
  if (file.exists(fls_alpha_res[i])) {
    raster(fls_alpha_res[i])
  } else {
    resample(raster(fls_alpha[i]), rst_ref_utm, filename = fls_alpha_res[i], 
             format = "GTiff", overwrite = TRUE)
  }
}

rst_alpha <- stack(lst_alpha); rm(lst_alpha)
dts_alpha <- substr(basename(fls_alpha), 10, 16)
dat_alpha <- data.frame(datetime = dts_alpha, albedo = fls_alpha_res, 
                        stringsAsFactors = FALSE)

## shortwave upward radiation
dat_rsu <- merge(dat_alpha, dat_rsd, by = "datetime")

dir_rsu <- paste0(dir_sensor, "/Shortwave_Upward_Radiation")
if (!dir.exists(dir_rsu)) dir.create(dir_rsu)

fls_rsu <- gsub("Shortwave_Downward_Radiation", "Shortwave_Upward_Radiation", 
                dat_rsu$swdr)
fls_rsu <- paste0(dir_rsu, "/", basename(fls_rsu))

lst_rsu <- foreach(i = 1:nrow(dat_rsu), .packages = lib) %dopar% {
  if (file.exists(fls_rsu[i])) {
    raster(fls_rsu[i])
  } else {
    swdr <- resample(projectRaster(raster(dat_rsu$swdr[i]), 
                                   crs = "+init=epsg:21037"), rst_ref_utm)
    
    albedo <- raster(dat_rsu$albedo[i])

    overlay(swdr, albedo, fun = function(x, y) y * x, 
            filename = fls_rsu[i], format = "GTiff", overwrite = TRUE)
  }
}
  
rst_rsu <- stack(lst_rsu); rm(lst_rsu)

dat_rsu <- data.frame(datetime = dat_rsu$date, swur_agg = fls_rsu, 
                      stringsAsFactors = FALSE)


### land surface temperature -----

fls_lst <- list.files(paste0("data/", sensor, "11A1.005/qc"), full.names = TRUE, 
                      pattern = "Day_1km.tif$")

fls_lst_res <- gsub("/qc/", "/res/", fls_lst)
dir_lst_res <- unique(dirname(fls_lst_res))
if (!dir.exists(dir_lst_res)) dir.create(dir_lst_res)

rst_lst <- stack(foreach(i = 1:length(fls_lst), .packages = lib) %dopar% {
  if (file.exists(fls_lst_res[i])) {
    raster(fls_lst_res[i])
  } else {
    resample(raster(fls_lst[i]), rst_ref_utm, filename = fls_lst_res[i], 
             format = "GTiff", overwrite = TRUE)
  }
})

dts_lst <- substr(basename(fls_lst), 10, 16)
dat_lst <- data.frame(datetime = dts_lst, lst = fls_lst_res, 
                      stringsAsFactors = FALSE)


### broadband surface emissivity -----

## ndvi
fls_ndvi <- list.files("data/MCD09GQ.006/ndvi/gf", full.names = TRUE, 
                       pattern = "^NDVI_MCD09GQ.*.tif$")

dts_ndvi <- substr(basename(fls_ndvi), 15, 21)
dat_ndvi <- data.frame(datetime = dts_ndvi, ndvi = fls_ndvi, 
                       stringsAsFactors = FALSE)

## lai
fls_lai <- list.files("data/MCD15A3H.006/rpl", full.names = TRUE, 
                      pattern = "^MCD15A3H.*.tif$")

dts_lai <- extractDate(fls_lai, 11, 17)$inputLayerDates
id <- which(!dts_lai %in% dts_ndvi)
fls_lai <- fls_lai[-id]

fls_lai_res <- gsub("/rpl/", "/res/", fls_lai)
dir_lai_res <- unique(dirname(fls_lai_res))
if (!dir.exists(dir_lai_res)) dir.create(dir_lai_res)

rst_lai <- stack(foreach(i = 1:length(fls_lai), .packages = lib) %dopar% {
  if (file.exists(fls_lai_res[i])) {
    raster(fls_lai_res[i])
  } else {
    resample(raster(fls_lai[i]), rst_ref_utm, filename = fls_lai_res[i], 
             format = "GTiff", overwrite = TRUE)
  }
})

dat_lai <- data.frame(datetime = dts_lai[-id], lai = fls_lai_res, 
                      stringsAsFactors = FALSE)

## surface emissivity
dat_vi <- merge(dat_ndvi, dat_lai, by = "datetime")

dir_es <- paste0(dir_sensor, "/Surface_Emissivity")
if (!dir.exists(dir_es)) dir.create(dir_es)

fls_es <- gsub("Lai", "Surface_Emissivity", basename(dat_vi$lai))
fls_es <- paste0(dir_es, "/", fls_es)

rst_es <- stack(foreach(i = 1:length(fls_es), .packages = lib) %dopar% {
  if (file.exists(fls_es[i])) {
    raster(fls_es[i])
  } else {
    surfaceEmissivity(raster(dat_vi$ndvi[i]), raster(dat_vi$lai[i]), 
                      filename = fls_es[i], format = "GTiff", overwrite = TRUE)
  }
})

dts_es <- substr(sapply(strsplit(basename(fls_es), "\\."), "[[", 2), 2, 8)
dat_es <- data.frame(datetime = dts_es, emis_surf = fls_es, 
                     stringsAsFactors = FALSE)


### longwave upward radiation -----

## remove non-required lst layers
id <- which(!dts_lst %in% dts_ndvi)
rst_lst <- rst_lst[[-id]]

dir_lwur <- paste0(dir_sensor, "/Longwave_Upward_Radiation")
if (!dir.exists(dir_lwur)) dir.create(dir_lwur)

fls_lwur <- gsub("Surface_Emissivity", "Longwave_Upward_Radiation", names(rst_es))
fls_lwur <- paste0(dir_lwur, "/", fls_lwur, ".tif")

rst_lwur <- stack(foreach(i = 1:length(fls_lwur), .packages = lib) %dopar% {
  if (file.exists(fls_lwur[i])) {
    raster(fls_lwur[i])
  } else {
    lwur(lst = rst_lst[[i]], bse = rst_es[[i]], filename = fls_lwur[i], 
                format = "GTiff", overwrite = TRUE)
  }
})

fls_rlu <- list.files(dir_lwur, pattern = "Longwave.*.tif$", 
                      full.names = TRUE)

dts_rlu <- substr(basename(fls_rlu), 11, 17)
dat_rlu <- data.frame(datetime = dts_rlu, lwur_agg = fls_rlu, 
                      stringsAsFactors = FALSE)


### longwave downward radiation

## effective atmospheric emissivity
fls_ea <- gsub("Transmissivity", "Emissivity", fls_tau)
dir_ea <- unique(dirname(fls_ea))
if (!dir.exists(dir_ea)) dir.create(dir_ea)

rst_ea <- stack(foreach(i = 1:length(fls_ea), .packages = lib) %dopar% {
  if (file.exists(fls_ea[i])) {
    raster(fls_ea[i])
  } else {
    atmosphericEmissivity(rst_tau[[i]], filename = fls_ea[i], 
                          format = "GTiff", overwrite = TRUE)
  }
})

dts_ea <- substr(sapply(strsplit(basename(fls_ea), "_"), "[[", 2), 5, 11)
dat_ea <- data.frame(datetime = dts_ea, emissivity = fls_ea, 
                     stringsAsFactors = FALSE)

## longwave downward radiation
dat_rld <- merge(dat_ea, dat_lst, by = "datetime")
fls_rld <- gsub("Atmospheric_Emissivity", "Longwave_Downward_Radiation", 
                dat_ea$emissivity)
fls_rld <- gsub("/agg", "", fls_rld)
dir_rld <- unique(dirname(fls_rld))
if (!dir.exists(dir_rld)) dir.create(dir_rld)

rst_rld <- stack(foreach(i = 1:nrow(dat_rld), .packages = lib) %dopar% {
  if (file.exists(fls_rld[i])) {
    raster(fls_rld[i])
  } else {
    ea <- raster(dat_rld$emissivity[i])
    ea <- projectRaster(ea, crs = "+init=epsg:21037")
    ea <- resample(ea, rst_ref_utm)
    
    lst <- raster(dat_rld$lst[i])
    lst <- resample(lst, rst_ref_utm)
    
    lwdr(ea, ta = lst, filename = fls_rld[i], 
         format = "GTiff", overwrite = TRUE)
  }
})

dts_rld <- substr(sapply(strsplit(basename(fls_rld), "_"), "[[", 2), 5, 11)
dat_rld <- data.frame(datetime = dts_rld, lwdr_agg = fls_rld, 
                      stringsAsFactors = FALSE)


# # ta (air temperature)
# ta.fls <- list.files("md07/processed", 
#                      pattern = "AGGMLY_.*Retrieved_Temperature_Profile.tif$", 
#                      full.names = TRUE)
# ta <- stack(ta.fls)
# 
# # dpt (dewpoint temperature)
# dpt.fls <- list.files("md07/processed", 
#                       pattern = "^AGGMLY_.*Retrieved_Moisture_Profile.tif$", 
#                       full.names = TRUE)
# dpt <- stack(dpt.fls)
# 
# # e_a (actual vapor pressure, taken from Ryu 2008)
# e.a <- 6.11 * exp(19.59 * (dpt - 273.15) / dpt)
# 
# # # (taken from Bisht 2005)
# # l.v <- 2.5 * 10^6
# # r.v <- 461
# # t0 <- 273.15
# # e.a2 <- 6.11 * exp(l.v / r.v * (1 / t0 - 1 / dpt))
# 
# r.ld <- 1.24 * (e.a / ta)^0.14 * sigma * ta^4


## R_n (net radiation)

dat_rn <- Reduce(function(...) merge(..., by = "datetime"), 
                 list(dat_rsd, dat_rsu, dat_rld, dat_es, dat_rlu))

fls_rn <- gsub("Longwave_Upward", "Net", dat_rn$lwur)
dir_rn <- unique(dirname(fls_rn))
if (!dir.exists(dir_rn)) dir.create(dir_rn)

lst_rn <- foreach(i = 1:nrow(dat_rn), .packages = lib) %dopar% {
  if (file.exists(fls_rn[i])) {
    raster(fls_rn[i])
  } else {
    rsd <- raster(dat_rn$swdr[i])
    rsd <- projectRaster(rsd, crs = "+init=epsg:21037")
    rsd <- resample(rsd, rst_ref_utm)
    
    rsu <- raster(dat_rn$swur[i])
    
    es <- raster(dat_rn$emis_surf[i])
    rld <- raster(dat_rn$lwdr[i])
    
    rlu <- raster(dat_rn$lwur[i])
    
    netRadiation(rsd, rsu, rld, es, rlu, filename = fls_rn[i], 
                 format = "GTiff", overwrite = TRUE)
  }
}

rst_rn <- stack(lst_rn)

stopCluster(cl)
