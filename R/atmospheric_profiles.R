### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## set working directory
Orcs::setwdOS(path_lin = "/media/fdetsch/XChange/", path_win = "D:/",
              path_ext = "kilimanjaro/evapotranspiration/")

## load packages
library(doParallel)

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

## pressure levels
p <- c(5, 10, 20, 30, 50, 70, 100, 150, 200, 250, 
       300, 400, 500, 620, 700, 780, 850, 920, 950, 1000)


### rearrange data -------------------------------------------------------------

## loop over products
for (product in c("MOD07_L2", "MYD07_L2")) {
  
  ## status message
  cat("Commencing with processing of", product, "...\n")
  
  ## list available files and remove duplicates
  fls <- list.files("md07/raw006", full.names = TRUE, recursive = TRUE, 
                    pattern = paste(product, "006", ".hdf$", sep = ".*"))

  if (anyDuplicated(basename(fls)) != 0)
    fls <- fls[-which(duplicated(basename(fls)))]

  ## sort by datetime and stop if any date is missing
  meta <- Orcs::list2df(strsplit(basename(fls), "\\."), stringsAsFactors = FALSE)
  names(meta) <- c("product", "date", "time", "collection", "processing", "filetype")
  
  meta$datetime <- strptime(paste(meta$date, meta$time), format = "A%Y%j %H%M")
  fls <- fls[order(meta$datetime)]; # meta <- meta[order(meta$datetime), ]
  
  if (length(unique(meta$date)) != 
        length(seq(as.Date("2013-01-01"), as.Date("2015-12-31"), "day")))
    stop("There are dates missing. Please check available files.\n")
  
  ## copy files to destination folder
  foreach(i = fls) %dopar% {
    fls_out <- paste0("data/", product, ".006/", basename(i))
    if (!file.exists(fls_out)) file.copy(i, fls_out)
  }
}

