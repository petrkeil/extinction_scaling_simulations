# 6/12/2023 by Petr Keil, keil@fzp.czu.cz

# Function that extracts number of extinctions at a site
# Arguments:
# MAT - matrix containing abundances (or incidences) with species as columns, 
#       and time steps as rows
# time1 - first (reference) time step
# time2 - a second time step, in which the extinctions will be counted

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

# Function that extracts per-species extinction rate at a site
# Arguments:
# MAT - matrix containing abundances (or incidences) with species as columns, 
#       and time steps as rows
# time1 - first (reference) time step
# time2 - a second time step, in which the extinctions will be counted

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
