switch(Sys.info()[["sysname"]], 
       "Linux" = setwd("/media/fdetsch/XChange/kilimanjaro/evapotranspiration/"), 
       "Windows" = setwd("D:/kilimanjaro/evapotranspiration/"))

lib <- c("raster", "rgdal", "doParallel")
sapply(lib, function(x) library(x, character.only = TRUE))

registerDoParallel(cl <- makeCluster(3))


## R_sd (shortwave downward radiation)

# G_sc (solar constant, 1367 Wm^(-2))
g.sc  <- 1367

# theta (solar incidence angle, i.e. solar zenith angle)
theta.fls <- list.files("md07/processed/", pattern = "^CRP_.*Solar_Zenith.tif$", 
                        full.names = TRUE)

theta.fls.split <- strsplit(basename(theta.fls), "\\.")
theta.fls.date <- as.Date(gsub("A", "", sapply(theta.fls.split, "[[", 2)), "%Y%j")
theta.fls.time <- sapply(theta.fls.split, "[[", 3)

theta.index.day <- as.numeric(theta.fls.time) >= 600 &
  as.numeric(theta.fls.time) <= 1800
# theta.day <- which(theta.index.day)
# theta.night <- which(!theta.index.day)

theta.fls.day <- theta.fls[theta.index.day]

df.theta.fls.day <- data.frame(date = theta.fls.date[theta.index.day], 
                               theta.fls.day, stringsAsFactors = FALSE)

# p (surface pressure from MODIS)
p.fls <- list.files("md07/processed/", pattern = "^CRP_.*Surface_Pressure.tif$", 
                    full.names = TRUE)

p.fls.split <- strsplit(basename(p.fls), "\\.")
p.fls.date <- as.Date(gsub("A", "", sapply(p.fls.split, "[[", 2)), "%Y%j")
p.fls.time <- sapply(p.fls.split, "[[", 3)

p.index.day <- as.numeric(p.fls.time) >= 600 &
  as.numeric(p.fls.time) <= 1800
# p.day <- which(p.index.day)
# p.night <- which(!p.index.day)

p.fls.day <- p.fls[p.index.day]

p <- stack(p.fls.day)
p <- p * 0.1 # hPa -> kPa

# p (surface pressure from ECMWF)
p.fls <- list.files("ecmwf/", pattern = "^DLY.*surface_pressure.*.tif$", 
                    full.names = TRUE)

p.fls.split <- strsplit(basename(p.fls), "_")
p.fls.date <- as.Date(sapply(p.fls.split, "[[", 6), "%Y%j")

df.p.fls <- data.frame(date = p.fls.date, p.fls, stringsAsFactors = FALSE)

k.t <- 1

# W (atmospheric water content)
w.fls <- list.files("md07/processed/", pattern = "^CRP_.*Water_Vapor.tif$", 
                    full.names = TRUE)

w.fls.split <- strsplit(basename(w.fls), "\\.")
w.fls.date <- as.Date(gsub("A", "", sapply(w.fls.split, "[[", 2)), "%Y%j")
w.fls.time <- sapply(strsplit(basename(w.fls), "\\."), "[[", 3)
w.index.day <- as.numeric(w.fls.time) >= 600 & as.numeric(w.fls.time) <= 1800
w.day <- which(w.index.day)
w.night <- which(!w.index.day)

w.fls.day <- w.fls[w.day]

df.w.fls.day <- data.frame(date = w.fls.date[w.day], 
                           w.fls.day, stringsAsFactors = FALSE)

# Continuous data.frame with corresponding theta, w and p files
df.seq <- data.frame(date = seq(as.Date("2013-01-01"), as.Date("2013-12-31"), 1))
ls.seq.theta.w.p.fls.day <- list(df.seq, df.theta.fls.day, df.w.fls.day, df.p.fls)
df.seq.theta.w.p.fls.day <- Reduce(function(...) merge(..., by = "date", all = TRUE), 
                                   ls.seq.theta.w.p.fls.day)
df.seq.theta.w.p.fls.day.cc <- df.seq.theta.w.p.fls.day[complete.cases(df.seq.theta.w.p.fls.day), ]

# Import corresponding raster files
rst.theta.w.p <- foreach(i = 2:4, j = c("theta", "w", "p"), .packages = lib) %dopar% {
  tmp.rst <- stack(df.seq.theta.w.p.fls.day.cc[, i])
  
  if (j == "theta") {
    return(cos(tmp.rst))
  } else if (j == "w") {
    return(tmp.rst * 10) # cm -> mm
  } else {
    return(tmp.rst * 0.1) # hPa -> kPa
  }
}

rst.theta.w.p[[3]] <- resample(rst.theta.w.p[[3]], rst.theta.w.p[[1]])

theta.cos <- rst.theta.w.p[[1]]
w <- rst.theta.w.p[[2]]
p <- rst.theta.w.p[[3]]

# tau_sw (broadband atmospheric transmissivity)
tau.sw <- 0.35 + 0.627 * exp((-0.00146 * p / (k.t * theta.cos)) - (0.075 * (w / theta.cos)^0.4))

# d_r (inverse squared relative earth-sun distance)
doy.fls <- substr(basename(df.seq.theta.w.p.fls.day.cc[, 2]), 19, 21)
doy <- as.numeric(doy.fls)

d.r <- 1 + 0.033 * cos(doy * 2 * pi / 365)

# R_sd
r.sd <- g.sc * d.r * tau.sw * theta.cos
r.sd <- writeRaster(r.sd, filename = "model_input/R_sd_13", format = "GTiff", 
                    overwrite = TRUE, bylayer = FALSE)

indices <- as.numeric(as.factor(substr(df.seq.theta.w.p.fls.day.cc[, 1], 6, 7)))
tmp <- stackApply(r.sd, indices, fun = median)
plot(tmp)


## R_lu (longwave upward radiation)

# NDVI
ndvi.fls <- list.files("myd09gq/processed/", pattern = "^NDVI", 
                       full.names = TRUE)
ndvi <- stack(ndvi.fls)

# LAI
lai <- stack("model_input/lai_13.tif")

# sigma (Stefan-Boltzmann constant, 5.67 * 10^(-8) Wm^(-2)K^(-4))
sigma <- 5.67 * 10^(-8)

# lst (land surface temperature)
lst.fls <- list.files("md11a1/processed/", 
                      pattern = "^RSMPL_AGG1MTH_.*LST_Day_1km.tif$", 
                      full.names = TRUE)
lst <- stack(lst.fls)

# epsilon_s (broadband surface emissivity)
epsilon.s <- stack(lapply(seq_len(nlayers(ndvi)), function(i) {
  x <- ndvi[[i]]
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