### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## set working directory
Orcs::setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/",
              path_ext = "kilimanjaro/evapotranspiration/")

## pressure levels
p <- c(5, 10, 20, 30, 50, 70, 100, 150, 200, 250, 
       300, 400, 500, 620, 700, 780, 850, 920, 950, 1000)

## retrieve metadata
fls <- list.files("md07/raw006", full.names = TRUE, recursive = TRUE, 
                  pattern = paste("^MOD07_L2", "006", ".hdf$", sep = ".*"))

## remove duplicate entries
if (anyDuplicated(basename(fls)))
  fls <- fls[-which(duplicated(basename(fls)))]

## timestamps
timestamps <- sapply(strsplit(basename(fls), "\\."), "[[", 3)
timestamps <- as.integer(timestamps)
id_day <- which(timestamps >= 600 & timestamps <= 2000)
fls <- fls[id_day]

## dates
dates <- sapply(strsplit(basename(fls), "\\."), "[[", 2)
dates <- as.Date(dates, format = "A%Y%j")
