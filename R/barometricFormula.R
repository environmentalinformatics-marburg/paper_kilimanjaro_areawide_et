barometricFormula <- function(z, gp, ta, p.levels, ...) {

  diff.abs <- abs(z - gp)
  id.min <- which.min(diff.abs)

  h0 <- gp[id.min]
  h1 <- z
  t <- ta[id.min]
  p0 <- p.levels[id.min]

  ## apply barometric formula 
  ## (taken from https://en.wikipedia.org/wiki/Barometric_formula)
  dh <- h1 - h0
  g <- 9.80665
  M <- 0.0289644
  R <- 8.31432

  p0 * exp((-1) * (g * M * dh) / (R * (t + 273.15)))
}
