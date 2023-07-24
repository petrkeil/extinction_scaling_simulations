
# function that extracts number of extinctions at a site
EXTract <- function(MAT, time1, time2)
{
  MAT <- MAT[c(time1, time2),]
  MAT.pres <- MAT
  MAT.pres[MAT > 0] <- 1 # convert abundances to presences
  # indicate extinction events
  MAT.ext  <- 1 * (MAT.pres[-1,] == 0 & MAT.pres[-nrow(MAT.pres),] > 0)
  Ex <- sum(MAT.ext)
  return(Ex)
}

# function that extracts per-species extinction rate at a site
PXTract <- function(MAT, time1, time2)
{
  MAT <- MAT[c(time1, time2),]
  MAT.pres <- MAT
  MAT.pres[MAT > 0] <- 1 # convert abundances to presences
  # indicate extinction events
  MAT.ext  <- 1 * (MAT.pres[-1,] == 0 & MAT.pres[-nrow(MAT.pres),] > 0)
  S1 <- sum(MAT.pres[1,])
  Px <- sum(MAT.ext) / S1 
  return(Px)
}
