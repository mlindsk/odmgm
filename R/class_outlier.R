new_outlier_model <- function(A,
                              sims,
                              sims_mixed,
                              sims_disc_raw,
                              mu,
                              sigma,
                              cdf,
                              cms,
                              sms,
                              pa,
                              tv,
                              vs,
                              hom,
                              approx) {
  structure(
    list(
      A             = A,
      sims          = sims,
      sims_mixed    = sims_mixed,
      sims_disc_raw = sims_disc_raw,
      mu            = mu,
      sigma         = sigma,
      cdf           = cdf,
      cms           = cms,
      sms           = sms,
      pa            = pa,
      vs            = vs,
      hom           = hom,
      approx        = approx
    ),
    class = c("outlier_model", tv, "list")
  )
}

new_outlier <- function(m, dev, pv, cv, a) {
  # m : outlier_model object
  m$dev    <- dev
  m$pval   <- pv
  m$cv     <- cv
  m$alpha  <- a
  class(m) <- c("outlier", class(m))
  return(m)
}
