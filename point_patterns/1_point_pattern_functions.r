# FUNCTIONS FOR SIMULATING AND ANALYZING EXTINCTION SCALING IN POINT-PATTERN
# COMMUNITIES

# Author: Petr Keil
# Email: keil@fzp.czu.cz

# ------------------------------------------------------------------------------

## Vojta Bartak's function 

# Parameters:
# alpha  ... interval  <0, 1>
# beta ... interval <-1, 1>
# N ... vector of population abundances
# plot ... should the curve be plotted?

# Value:
# A vector of per-individual probabilities of death

effect.f <- function(alpha=0.5, beta=0.5, N = seq(1,100), plot = TRUE)
{
  # the main calculation:
  P <- exp( N*sign(beta)*exp(-1/abs(beta)) ) /
    ( (1-alpha)/alpha + exp(N*sign(beta)*exp(-1/abs(beta))) )
  
  if(plot)
  {
    plot(N, P, type="l", ylim=c(0,1), main = "Per-individual death rate")
    abline(v=1, col = "grey", lty = 2)
    abline(h=alpha, col = "grey", lty = 2)
  }
  
  return(P)
}

# Test:
# effect.f(alpha=0.5, beta=0.5, N = seq(1,100), plot = TRUE)


# ------------------------------------------------------------------------------

## Create a lognormal SAD, given N, S, and coefficient of variation of abundance

# Parameters:
# N ... total number of individuals in the community
# S ... total number of species
# cv_abund.vect ... vector of coefficients of variation of an SAD

# Value: data frame with parameters of the SAD, 
#        and probability mass (PM) values for a given abundance

SAD.prob.mass <- function(N, S, cv_abund.vect = c(0.1, 1, 10))
{
  dat <- list()
  for(i in 1:length(cv_abund.vect))
  {
    mean_abund = N/S
    sd_abund <- mean_abund * cv_abund.vect[i]
    sigma1 <- sqrt(log(sd_abund^2/mean_abund^2 + 1))
    mu1 <- log(mean_abund) - sigma1^2/2
    PM <- dlnorm(1:40, meanlog = mu1, sdlog=sigma1)
    x.max <- 40
    dat[[i]] <- data.frame(N = 1:x.max, PM = PM, cv_abund = cv_abund.vect[i])
  }
  #print(dat)
  dat <- do.call("rbind", dat)
  return(dat)
}


# ------------------------------------------------------------------------------

# Function that takes a community object (com) subjects its individuals to 
# death according to the Vojta Bartak's function, and then overlays grids of 
# different resolutions over the community and calculates the extinction scaling

# Parameters:
# com ... community produced by the mobsim's sim_thomas_community function
# side.divisions ... vector of grid resolutions, i.e. lengths of the side of a grid cell
# alpha ... first parameters of Bartak's function
# alpha.const ... is alpha given, or ignored and drawn from a uniform 
#                 distribution between 0.01 and 0.99
# beta ... second parameter of Bartak's function
# plot.beta ... should the communities and extiction scaling be plotted?

# Value:
# data frame with parameter values, side divisions, extinction numbers,
# and probabilities

com.extinction <- function(com, 
                           side.divisions,
                           alpha,
                           alpha.const = TRUE,
                           beta,
                           plot.res = TRUE)
{
  # Species-Abundance Distribution - SAD
  SAD <- data.frame(table(com$census$species))
  names(SAD) <- c("species", "N.tot")
  SAD$species <- as.character(SAD$species)
  
  if (alpha.const) { SAD$alpha <- alpha }
  else { SAD$alpha <- runif(nrow(SAD), 0.01, 0.99) }
  
  census <- left_join(com$census, SAD, by="species")

  # per-capita probability of death
  P.minus.vect <- effect.f(alpha = census$alpha, beta, N = census$N.tot, plot = FALSE)

  # sample individuals that will be killed
  deaths <- rbinom(n = length(P.minus.vect), size = 1, prob = P.minus.vect)
 
  # kill the individuals
  com2 <- com
  com2$census <- com2$census[deaths == 0, ]
  
  # sample the grids at given resolutions
  
  res <- list()
  i = 1
  for(div in side.divisions)
  {
    res[[i]] <- sample.grids(com, com2, div)
    i = i + 1
  }
  res <- do.call("rbind", res)
  res <- data.frame(res, 
                    mean.P.minus = mean(P.minus.vect),
                    mean.log.abund = mean(log(community_to_sad(com))),
                    P.x = res$E.x / res$S.t1 )
  
  if(plot.res)
  {
    # plot the communities
    plot(com, main = "Time 1")
    plot(com2, main = "Time 2")
    plot(log10(res$A), res$P.x, type = "b", xlab="log10 Area", 
         ylab="P.x", ylim=c(0,1), main = "Per-species extinction rate")  
    plot(log10(res$A), log10(res$E), type = "b", 
         xlab="log10 Area", ylab="log10 E",
         main = "Number of extinct species")  
  }
  return(res)
}

# Testing:
#com <- sim_thomas_community(n_sim  = 100,
#                            s_pool = 10,
#                            sad_coef = list(cv_abund = 1),
#                            mother_points = 1,
#                            sigma = 0.1)

#com.extinction(com, 
#               side.divisions = c(2,4,8,16), 
#               beta = -0.4, alpha = 0.5,
#               alpha.const = TRUE)

# ------------------------------------------------------------------------------

# Function that takes two communities and overlays them with a grid of a given
# resolution, and calculates species richness, extinctions, etc.

# Arguments:
# com1 ... community produced by the mobsim's sim_thomas_community function, time 1
# com2 ... community produced by the mobsim's sim_thomas_community function, time 2
# side.division ... length of the side of a grid cell

# Value:
# data frame with area, mean species richness, extinction numbers,
# and probabilities

sample.grids <- function(com1, com2, side.division)
{
  # convert the division into area
  A = (1/side.division)^2

  # sample the community at two time points and at two resolutions
  grid1 <- sample_quadrats(com1,
                           plot = FALSE,
                           method = "grid",
                           n_quadrats = 1/A,
                           quadrat_area = A,
                           delta_x = sqrt(A),
                           delta_y = sqrt(A))
  
  grid2 <- sample_quadrats(com2,
                           plot = FALSE,
                           method = "grid",
                           n_quadrats = 1/A,
                           quadrat_area = A,
                           delta_x = sqrt(A),
                           delta_y = sqrt(A))
  
  # create community matrix (site x species)
  abund1 <- grid1$spec_dat
  abund2 <- grid2$spec_dat
  
  # convert abundances to incidences
  incid1 <- abund1
  incid2 <- abund2
  incid1[abund1 > 0] <- 1
  incid2[abund2 > 0] <- 1
  
  # species richness
  S.t1 <- rowSums(incid1)
  S.t2 <- rowSums(incid2)
  
  # numbers of extinct species
  E.x <- rowSums( incid2 - incid1 == -1 )
  
  # return the mean values
  return(data.frame(A = A,
                    S.t1 = mean(S.t1, na.rm=TRUE),
                    S.t2 = mean(S.t2, na.rm=TRUE),
                    E.x = mean(E.x, na.rm=TRUE) ))
  
}


# ------------------------------------------------------------------------------

# Function that takes a vector of simulation parameters, simulates communities
# from them (using mobsim's sim_thomas_community), overlays the grids, 
# and calculates extinction scaling

# Parameters:
# params ... a data.frame with simulation parameters, this is an example:
# params <- expand.grid(alpha = c(0.1, 0.3, 0.5, 0.7, 0.9),
#                      beta = c(-4, -1, -0.4, 0, 0.4, 1, 4),
#                      N.tot = c(100, 1000) , 
#                      S.frac = c(0.05, 0.1, 0.2), 
#                      sigma = c(0.01, 0.1, 1),
#                      cv.abund = c(0.1, 1, 10))
# params <- data.frame(params, S.tot = params$N.tot*params$S.frac)


main.loop <- function(params)
{
  # empty container for results
  results <- list()
  for(i in 1:nrow(params))
  {
    # simulation progress info
    print(paste("Simulation", i, "out of", nrow(params)))
    
    # simulate the community
    # and if realized # of species is 1, try again
    repeat
    {
      # simulate point-pattern community
      com <- mobsim::sim_thomas_community(n_sim  = params$N.tot[i],
                                          s_pool = params$S.tot[i],
                                          sigma  = params$sigma[i],
                                          mother_points = 1,
                                          sad_coef = list(cv_abund = params$cv.abund[i]))
      # actual number of species
      S <- length(unique(com$census$species))
      if(S > 1) { break }
    }
    
    # which grid resolutions should we put over the simulated community?
    side.divisions = c(2,4,8,16)
    
    # calculate "deltas" for the set of  grid resolutions and simulation 
    deltas <- com.extinction(com, 
                             side.divisions, 
                             alpha = params$alpha[i],
                             alpha.const = params$alpha.const[i],
                             beta = params$beta[i],
                             plot.res = FALSE)

    # slopes of the SAR relationship in time 1 and 2
    slope.S1 <- glm(S.t1 ~ log(A), family = "poisson", data = deltas)$coefficients[2]

    # calculate the slopes only if there were any extinctions at all
    if(max(deltas$E.x) > 0)
    {
      # slope of the P.x-Area relationship
      slope.P <- lm(P.x ~ log(A), data = deltas)$coefficients[2]
      # slope of the E-Area relationship
      slope.E <- glm(E.x ~ log(A), family = "poisson", data = deltas)$coefficients[2]
      # take average of the per-individual death probability, and abundance
    }
    if(max(deltas$E.x) == 0) # if no extinctions happen, return NAs
    {
      slope.P <- NA
      slope.E <- NA
    }

    mean.P.minus <- deltas$mean.P.minus[1]
    mean.log.abund <- deltas$mean.log.abund[1]    
    E.happened <- max(deltas$E.x) > 0 # have any extinctions happened at all?
      
    res <- data.frame(params[i,], 
                      slope.P = slope.P, 
                      slope.E = slope.E,
                      slope.S1 = slope.S1,
                      mean.P.minus = mean.P.minus,
                      mean.log.abund = mean.log.abund,
                      E.happened)
    
    # store the results in the container
    results[[i]] <- res
  }
  # merge the results to a data frame
  results <- data.frame(do.call(what = rbind, args = results))
  return(results)
}

# ------------------------------------------------------------------------------
# Function that can be used in the interactive plots using the 'manipulate'
# package.

run.all <- function(alpha, 
                    beta, 
                    side.divisions,
                    N,
                    S,
                    sigma,
                    cv.abund)
{
  
  com <- sim_thomas_community(n_sim  = N,
                              s_pool = S,
                              sad_coef = list(cv_abund = cv.abund),
                              mother_points = 1,
                              sigma  = sigma)
  
  hist(community_to_sad(com), 
       nclass = 20,
       xlab = "Abundance", 
       main = "Species Abundance Distribution")
  
  ef <- effect.f(alpha, beta, plot = TRUE)
  
  deltas <- com.extinction(com, 
                           side.divisions, 
                           alpha = alpha, 
                           beta = beta,
                           plot.res = TRUE)
  #coef <- glm(E.x ~ log(A), family = "poisson", data = deltas)$coefficients[2]
  #print(coef)
  #return(deltas)
}

#par(mfrow=c(3,2))
#del <- run.all(alpha = 0.9, beta = -0.4, side.divisions= c(2,4,8,16), N = 100, 
#        S = 100*0.05, sigma = 0.01, cv.abund = 0.1)

