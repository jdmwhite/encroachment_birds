#________________________________________________#
##### Load libraries #############################
#________________________________________________#
library(unmarked)
library(tidyverse)

#________________________________________________#
##### Load processed SABAP2 cards ################
#________________________________________________#

files <- list.files('data/sabap2_processed_cards/') %>% str_remove('.csv') %>% as.numeric()
all_spp_est <- data.frame()

for(ii in files){

# load processed SABAP2 dataset
df <- read_csv(paste0('data/sabap2_processed_cards/',ii, '.csv'))

#________________________________________________#
########## Load colext models ####################
#________________________________________________#

#### load colext model
load(paste0("processed/models/",ii,".RData"))

#________________________________________________#
############ Extract lambda ######################
#________________________________________________#
# see https://assets.peregrinefund.org/docs/project-data/book-applied-raptor-ecology/ch-10-applied-raptor-ecology.pdf
# see McKenzie et al. 2003 for reference: https://doi.org/10.1890/02-3090
smoothed_occ <- smoothed(sel_mod)[2,] # extract occupied vals

lambda <- rep(NA, times = (length(smoothed_occ)-1))
for(i in 1:(length(smoothed_occ)-1)){
  lambda[i] <- smoothed_occ[i+1]/smoothed_occ[i]
}

# Noting that a mean λ<1 indicates population contraction and λ>1 indicates occupancy expansion
mean <- mean(lambda)

# calculate 90% confidence intervals
SD <- sd(lambda)
n <- length(lambda)
margin <- qt(0.95, df = n - 1)*SD/sqrt(n)
lower_interval_90 <- mean - margin
upper_interval_90 <- mean + margin

# combine values together
lambda_vals <- c(ii, mean, lower_interval_90, upper_interval_90)

# format
lambda_df <- data.frame(t(lambda_vals))
names(lambda_df) <- c('species','lambda_mean', 'lambda_0.05', 'lambda_0.95')

#________________________________________________#
############ Extract WC estimates ################
#________________________________________________#
sum_m1 <- summary(sel_mod)
# occ estimate for woody cover
occ_est <- unmarked::coef(sel_mod)[2] # est
occ_p_val <- sum_m1$psi$`P(>|z|)`[2]
occ_90CI <- confint(sel_mod, type = 'psi', level = 0.9)[2,] #90% CIs
occ_vals <- as.data.frame(cbind(ii, occ_est, occ_p_val, occ_90CI[1], occ_90CI[2]))
names(occ_vals) <- c('species', 'occ_est', 'occ_p_val', 'occ_CI_0.05', 'occ_CI_0.95')

# col estimate for woody cover
col_est <- unmarked::coef(sel_mod)[6]
col_p_val <- sum_m1$col$`P(>|z|)`[2]
col_90CI <- confint(sel_mod, type = 'col', level = 0.9)[2,] #90% CIs
col_vals <- as.data.frame(cbind(ii, col_est, col_p_val, col_90CI[1], col_90CI[2]))
names(col_vals) <- c('species', 'col_est', 'col_p_val', 'col_CI_0.05', 'col_CI_0.95')

# ext estimate for woody cover
ext_est <- unmarked::coef(sel_mod)[8]
ext_p_val <- sum_m1$ext$`P(>|z|)`[2]
ext_90CI <- confint(sel_mod, type = 'ext', level = 0.9)[2,] #90% CIs
ext_vals <- as.data.frame(cbind(ii, ext_est, ext_p_val, ext_90CI[1], ext_90CI[2]))
names(ext_vals) <- c('species', 'ext_est', 'ext_p_val', 'ext_CI_0.05', 'ext_CI_0.95')

# join parameter estimates
param_vals <- occ_vals %>%
  left_join(col_vals, by = 'species') %>%
  left_join(ext_vals, by = 'species')

#________________________________________________#
############ Calculate Z-scores ##################
#________________________________________________#

# extract presence/absence data and convert to long format
df %>% 
  dplyr::select(Pentad, year.1.1:year.10.4) %>% 
  pivot_longer(cols = year.1.1:year.10.4, names_to = 'year', values_to = 'pr_ab') %>%
  mutate(year = substr(year, 1, nchar(year)-2)) %>%
  mutate(year = sapply(str_split(year, '[.]'), '[', 2)) %>%
  filter(!is.na(pr_ab)) -> df_long
  
# calculate total cards and reporting rate per epoch
df_long %>%
  mutate(year = as.numeric(year)) %>% 
  filter(year > 1) %>%
  mutate(epoch = case_when(
      year <= 4 ~ 'first_3', 
      year > 7 ~ 'last_3'))  %>% 
  filter(!is.na(epoch)) %>%
  group_by(epoch, pr_ab) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::reframe(rep_rate = n / sum(n),
                     n_total = sum(n)) %>%
  rowid_to_column() %>%
  filter(rowid %in% c(2,4)) %>%
  dplyr::select(-rowid) -> rep_rate_epoch_3s
  
# extract total cards and reporting rate for each epoch
first_epoch_total <- rep_rate_epoch_3s$n_total[1]
first_epoch_rep_rate <- rep_rate_epoch_3s$rep_rate[1]
last_epoch_total <- rep_rate_epoch_3s$n_total[2]
last_epoch_rep_rate <- rep_rate_epoch_3s$rep_rate[2]
  
# calculate P = pooled reporting rate across both epochs
P <- (first_epoch_total * first_epoch_rep_rate + last_epoch_total * last_epoch_rep_rate)/
    (first_epoch_total + last_epoch_total)
  
# calculate Z-scores
z_score <- (last_epoch_rep_rate - first_epoch_rep_rate) /
    sqrt(P * (1 - P)*((1/first_epoch_total) + (1/last_epoch_total)))
  
# format
z_output <- data.frame(ii, z_score)
names(z_output) <- c('species', 'z_score')

# Combine together lambda, occ model estimates and z-scores
all_est <- lambda_df %>%
  left_join(param_vals, by = 'species') %>%
  left_join(z_output, by = 'species')

#________________________________________________#
############ Determine categories ################
#________________________________________________#
all_est %>% mutate(lambda_category =
                   case_when(lambda_0.95 > 1 & lambda_0.05 > 1 ~ 'increasing',
                             lambda_0.95 >= 1 & lambda_0.05 <= 1 ~ 'stable',
                             lambda_0.95 >= 1 & lambda_0.05 <= 1 ~ 'stable',
                             lambda_0.95 < 1 & lambda_0.05 < 1 ~ 'decreasing')) -> all_est

all_est %>% mutate(wc_response =
                   case_when(occ_CI_0.95 > 0 & occ_CI_0.05 > 0 & occ_p_val < 0.05 ~ 'positive',
                             occ_CI_0.95 < 0 & occ_CI_0.05 < 0 & occ_p_val < 0.05 ~ 'negative',
                             occ_CI_0.95 > 0 & occ_CI_0.05 <= 0 & occ_p_val < 0.05 ~ 'neutral',
                             occ_CI_0.95 >= 0 & occ_CI_0.05 <= 0 & occ_p_val < 0.05 ~ 'neutral',
                             occ_p_val >= 0.05 ~ 'neutral')) -> all_est

all_est %>% mutate(outcome = 
                   case_when(lambda_category == 'increasing' & wc_response == 'positive' ~ 'winner',
                             lambda_category == 'decreasing' & wc_response == 'negative' ~ 'loser',
                             lambda_category == 'increasing' & wc_response == 'negative' ~ 'indifferent',
                             lambda_category == 'decreasing' & wc_response == 'positive' ~ 'indifferent',
                             lambda_category == 'stable' ~ 'indifferent',                             
                             wc_response == 'neutral' ~ 'indifferent')) -> all_est

#________________________________________________#
####### Reformat and join on main df #############
#________________________________________________#

all_est %>% dplyr::select(species:lambda_0.95, lambda_category, z_score, occ_est, 
                          occ_CI_0.05, occ_CI_0.95, wc_response, outcome, col_est, 
                          col_CI_0.05:ext_est, ext_CI_0.05, ext_CI_0.95) -> all_est
all_spp_est <- rbind(all_spp_est, all_est)

}

write_csv(all_spp_est, 'processed/all_spp_estimates.csv')

#________________________________________________#
############ Export outputs ######################
#________________________________________________#


# END