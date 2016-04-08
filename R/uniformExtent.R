uniformExtent <- function(f = 1e3, ...) {

  ## load packages
  library(raster)
  library(rgdal)

  ## import plot coordinates
  suppressWarnings(
    spt_plots <- readOGR("data/shp", "PlotPoles_ARC1960_mod_20140807_final",
                         p4s = "+init=epsg:21037", ...)
  )

  ## relevant plots
  chr_plots <- c("sav5", "sav0", "mai4", "mai0", # colline zone
                 "cof3", "cof2", "gra3", "gra2", # submontane zone
                 "fed1", "hel1", "fer0")         # subalpine and lower alpine zone

  spt_plots <- subset(spt_plots, PlotID %in% chr_plots & PoleType == "AMP")

  ## extend bounding box
  ext_plots <- extent(spt_plots)

  num_xmin <- xmin(ext_plots) - f
  num_xmax <- xmax(ext_plots) + f
  num_ymin <- ymin(ext_plots) - f
  num_ymax <- ymax(ext_plots) + f

  extent(c(num_xmin, num_xmax, num_ymin, num_ymax))
}
