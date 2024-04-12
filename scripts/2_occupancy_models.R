#________________________________________________#
##### Load libraries #############################
#________________________________________________#
library(unmarked)
library(tidyverse)

#________________________________________________#
##### Load processed SABAP2 cards ################
#________________________________________________#
files <- list.files('data/sabap2_processed_cards/') %>% str_remove('.csv')

for(ii in files){
  
  # load processed SABAP2 cards
  df <- read_csv(paste0('data/sabap2_processed_cards/',ii, '.csv'))

  #________________________________________________#
  ##### Convert data to matrices ###################
  #________________________________________________#

  # Set up parameters
  M <- length(unique(df$Pentad)) # number of site = M
  J <- 4 # number of secondary sample periods (seasons)
  T <- 10 # number of primary sample periods (years)

  # Presence matrix
  yy <- matrix(as.numeric(unlist(c(df[,36:75]))), M, J*T)

  # yearly site covariates
  year <- matrix(as.factor(seq(1,10,1)), M, T, byrow = TRUE)

  # Creating matrix for wc_by_year
  wc_mat <- as.data.frame(year)
  df[,c('Pentad','wc','trend')] %>% drop_na() %>% distinct() -> wc_by_pentad

  wc_mat$Pentad <- wc_by_pentad$Pentad
  wc_mat$trend <- wc_by_pentad$trend
  wc_mat$V10 <- wc_by_pentad$wc *100 # convert to %
  wc_mat$V9 <- wc_mat$V10 - (wc_mat$trend * 1) # hindcast values based on trends
  wc_mat$V8 <- wc_mat$V10 - (wc_mat$trend * 2)
  wc_mat$V7 <- wc_mat$V10 - (wc_mat$trend * 3)
  wc_mat$V6 <- wc_mat$V10 - (wc_mat$trend * 4)
  wc_mat$V5 <- wc_mat$V10 - (wc_mat$trend * 5)
  wc_mat$V4 <- wc_mat$V10 - (wc_mat$trend * 6)
  wc_mat$V3 <- wc_mat$V10 - (wc_mat$trend * 7)
  wc_mat$V2 <- wc_mat$V10 - (wc_mat$trend * 8)
  wc_mat$V1 <- wc_mat$V10 - (wc_mat$trend * 9)

  # Scale the wc data with the mean and sd across all values
  wc_list <- c(wc_mat[c(1:10)])
  wc_mean <- mean(unlist(wc_list))
  wc_sd <- sd(unlist(wc_list))
  scale_wc <- function(x) (x - wc_mean) / wc_sd

  # apply the function and convert to a matrix
  wc_mat <- as.matrix(apply(wc_mat[1:10], 2, scale_wc))

  # extract the first years wc data
  wc_year1 <- wc_mat[,1]

  # Read in site data, select columns and rename
  siteCovs <- read_csv('data/siteCovs/siteCovs.csv')
  df %>% left_join(siteCovs, by = "Pentad") -> df
  minT <-  matrix(c(scale(df$minT)), M, 1, byrow = TRUE)
  ap_rat <- matrix(c(scale(df$ap_rat)), M, 1, byrow = TRUE)

  # Obs covariates: includes jday, total hours, road density
  jday <- matrix(as.numeric(unlist(c(df[,76:115]))), M, J*T)
  tHrs <- matrix(as.numeric(unlist(c(df[,156:195]))), M, J*T)
  r_dens <- matrix(as.numeric(unlist(c(df[,116:155]))), M, J*T)

  #________________________________________________#
  ##### Assemble unmarked df #######################
  #________________________________________________#
  umf <- unmarkedMultFrame(y = yy,
                            yearlySiteCovs = list(year = year, wc = wc_mat),
                            siteCovs = data.frame(wc_start = wc_year1, minT = minT, ap_rat = ap_rat),
                            obsCovs = list(jday = jday, tHrs = tHrs, r_dens = r_dens),
                            numPrimary = 10)

  #________________________________________________#
  ##### Run colext models ##########################
  #________________________________________________#
  sel_mod <- colext(psiformula= ~ wc_start + minT +  ap_rat, 
                    gammaformula = ~ wc, epsilonformula = ~ wc, 
                    pformula = ~ year + wc + jday + tHrs + r_dens, 
                    data = umf)
  # summary(sel_mod)

  save(sel_mod, file = paste0("processed/models/",ii,".RData"))
}

# END