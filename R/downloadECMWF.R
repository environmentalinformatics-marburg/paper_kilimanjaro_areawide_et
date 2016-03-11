downloadECMWF <- function(x, y, last, prm = c("geopotential", "temperature"), 
                          dsn = getwd(), overwrite = FALSE) {
  
  ## convert x, y, last to character
  if (is.character(x)) x <- as.Date(x)
  if (is.character(y)) y <- as.Date(y)
  if (is.character(last)) last <- as.Date(last)
  
  ## import python template
  fls <- system.file(paste0("extdata/download_", prm, ".py"), package = "reset")
  lns <- readLines(fls)
  
  ## loop over temporal range in 1-month steps
  st <- seq(x, y, "month")
  nd <- c(st[2:length(st)] - 1, last)
  
  out <- character(length(st))
  
  for (i in 1:length(st)) {
    
    # first run: create new dates and target file
    if (i == 1) {
      dates_old <- unlist(strsplit(lns[grep("date", lns)], "\""))[[4]]
      dates_new <- paste0(st[i], "/to/", nd[i])
      lns <- gsub(dates_old, dates_new, lns)
      
      target_old <- unlist(strsplit(lns[grep("target", lns)], "\""))[[4]]
      target_new <- paste(strftime(st[i], "%Y%j"), strftime(nd[i], "%Y%j"), sep = "_")
      target_new <- paste0(target_new, "__", prm, ".nc")
      lns <- gsub(target_old, target_new, lns)
    
    # subsequent runs: update dates and target file    
    } else {
      # replace date string
      lns <- gsub(st[i-1], st[i], lns)
      lns <- gsub(nd[i-1], nd[i], lns)
      # replace target file
      lns <- gsub(strftime(st[i-1], "%Y%j"), strftime(st[i], "%Y%j"), lns)
      lns <- gsub(strftime(nd[i-1], "%Y%j"), strftime(nd[i], "%Y%j"), lns)
    }
    
    ## write to file if file does not exist or `overwrite = TRUE`
    target <- unlist(strsplit(lns[grep("target", lns)], "\""))[[4]]
    target <- paste(dsn, target, sep = "/")
    
    ## write to file
    if (!file.exists(target) | overwrite) {
      fls_out <- paste(dsn, basename(fls), sep = "/")
      writeLines(lns, fls_out)
      system(paste0("cd ", dsn, "; python ", basename(fls_out)))
    }
    
    out[i] <- target
  }
  
  return(out)
}