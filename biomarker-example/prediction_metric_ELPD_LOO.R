rm(list = ls())

###### Generated data ########

# Load necessary libraries
library(rstan)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(bayesplot)
library(cowplot)
library(shinystan)
library(gridExtra)


# Parameters to load the wanted stan-fit. Default values are the only available on github
date <- '2025-04-17'
experiment <- 'all_but_last_cycle'
number_of_patients <- 10



######## create data_in ############
inputdatafile_PSA <- sprintf('%s_%d_PSA_df.csv', experiment, number_of_patients)
inputdatafile_day <- sprintf('%s_%d_day_df.csv', experiment, number_of_patients)
inputdatafile_Tx <- sprintf('%s_%d_Tx_df.csv', experiment, number_of_patients)

# Read data from files:
data_PSA <- read.table(inputdatafile_PSA, header = TRUE, sep = ",", stringsAsFactors = FALSE)
data_day <- read.table(inputdatafile_day, header = TRUE, sep = ",", stringsAsFactors = FALSE)
data_Tx <- read.table(inputdatafile_Tx, header = TRUE, sep = ",", stringsAsFactors = FALSE)


patient_IDs <- c(1, 4, 13, 14, 15, 17, 24, 26, 28, 29)

# Get the number of sampling times & experimental replicates:
no_reps = ncol(data_PSA);
no_ts_max = nrow(data_PSA)-1; #exclude initial row
ts_lengths = colSums(!is.na(data_day %>% slice(-1)))

#remove NA data
data_day[is.na(data_day)] <- 0
data_PSA[is.na(data_PSA)] <- 0
data_Tx[is.na(data_Tx)] <- 0

# Get initial & sampling times from data:
t0_data = data_day %>% slice(1)
ts_data = data_day %>% slice(-1)

t0_data <- as.numeric(unlist(t0_data))

# Get initial & sampling time population-sizes from data:
y1_data_0 = data_PSA %>% slice(1)
y1_data = data_PSA %>% slice(-1)

y1_data_0 <- as.numeric(unlist(y1_data_0))

# Tx from data:
# Tx_data = data_Tx %>% slice(-1)

Tx_data = data_Tx

# Transform data for later plotting:
inputdata_y1 <- data.frame(populationsize = unlist(y1_data[,1]),replicate = as.vector(matrix(rep(1,length(ts_data[1])),nrow=length(ts_data[1]),byrow=TRUE)),time = rep(ts_data[1],1))

# for plotting whole series

inputdatafile_PSA_all <- sprintf('%d_PSA_df.csv', number_of_patients)
inputdatafile_day_all <- sprintf('%d_day_df.csv', number_of_patients)
inputdatafile_Tx_all <- sprintf('%d_Tx_df.csv', number_of_patients)

# inputdatafile_Params_data <- sprintf('%d_Params_df.csv', number_of_patients)

data_PSA_all <- read.table(inputdatafile_PSA_all, header = TRUE, sep = ",", stringsAsFactors = FALSE)
data_day_all <- read.table(inputdatafile_day_all, header = TRUE, sep = ",", stringsAsFactors = FALSE)
data_Tx_all <- read.table(inputdatafile_Tx_all, header = TRUE, sep = ",", stringsAsFactors = FALSE)

# data_Params <- read.table(inputdatafile_Params_data, header = TRUE, sep = ",", stringsAsFactors = FALSE)


# Get the number of sampling times & experimental replicates:
no_ts_max_all = nrow(data_PSA_all)-1; #exclude initial row
ts_lengths_all = colSums(!is.na(data_day_all %>% slice(-1)))

#keep only no_reps columns
data_day_all = data_day_all[1:no_reps]
data_PSA_all = data_PSA_all[1:no_reps]
data_Tx_all = data_Tx_all[1:no_reps]

#remove NA data
data_day_all[is.na(data_day_all)] <- 0
data_PSA_all[is.na(data_PSA_all)] <- 0
data_Tx_all[is.na(data_Tx_all)] <- 0


# Get initial & sampling times from data:
ts_data_all = data_day_all %>% slice(-1)

Tx_data_all = data_Tx_all

#patient to plot
patient_col <- 1 #not used

patient_ID <- patient_IDs[patient_col]

# Put data in df for plotting:
plot_df<- data.frame(
  x = data_day_all[0:ts_lengths_all[patient_col]+1,patient_col],
  y = data_PSA_all[0:ts_lengths_all[patient_col]+1,patient_col]
)


nr_plot_points = nrow(plot_df)-1 # dont count 0

max_days <- max(data_day_all[, patient_col])
t_gd <- 1:max_days

pool_asgmt <- 1:no_reps


# Data to send to Stan:
data_in = list(
  no_reps = no_reps,
  no_ts_max = no_ts_max,
  
  t0_data = t0_data,
  ts_lengths = ts_lengths,
  ts_data = transpose(ts_data),
  
  y1_data_0 = y1_data_0,
  y1_data = transpose(y1_data),
  Tx_data = transpose(Tx_data),
  
  #for plotting
  no_ts_max_all = no_ts_max_all, 
  ts_data_all = transpose(ts_data_all),
  Tx_data_all = transpose(Tx_data_all),
  
  no_t_gd = length(t_gd),
  t_gd=t_gd,
  
  patient_col=patient_col,
  nr_plot_points = nr_plot_points
  
)

#read fit

fit1 <- readRDS(sprintf("fit_%s_%d_comp_pool_%s.rds",experiment, number_of_patients, date))
mod1 <- stan_model("model_comp_pool_gq.stan")
posterior_draws1 <- as.matrix(fit1, pars = fit1@model_pars)

fit2 <- readRDS(sprintf("fit_%s_%d_no_pool_%s.rds",experiment, number_of_patients, date))
mod2 <- stan_model("model_no_pool_gq.stan")
posterior_draws2 <- as.matrix(fit2, pars = fit2@model_pars)

fit3 <- readRDS(sprintf("fit_%s_%d_part_pool_%s.rds",experiment, number_of_patients, date))
mod3 <- stan_model("model_part_pool_gq.stan")
posterior_draws3 <- as.matrix(fit3, pars = fit3@model_pars)



#fit model to generate quantities

gq_fit_comp_pool <- gqs(mod1, draws = posterior_draws1, data = data_in)
gq_fit_no_pool <- gqs(mod2, draws = posterior_draws2, data = data_in)
gq_fit_part_pool <- gqs(mod3, draws = posterior_draws3, data = data_in)


y_fit_draws_comp_pool <- as.array(gq_fit_comp_pool, pars = "y_fit_tot") #for loo
y_fit_draws_no_pool <- as.array(gq_fit_no_pool, pars = "y_fit_tot") #for loo#
y_fit_draws_part_pool <- as.array(gq_fit_part_pool, pars = "y_fit_tot") #for loo



#### test for specific patient ####
patient_col <- 2

# All observed data for that patient
obs_df <- data.frame(
  x = data_day_all[0:ts_lengths_all[patient_col] + 1, patient_col],
  y = data_PSA_all[0:ts_lengths_all[patient_col] + 1, patient_col]
)

y_obs <- obs_df %>% pull(y)
x_obs <- obs_df %>% pull(x)


# Observed data for interference for that patient
obs_inference_df <- data.frame(
  x = data_day[0:ts_lengths[patient_col] + 1, patient_col],
  y = data_PSA[0:ts_lengths[patient_col] + 1, patient_col]
)

y_obs_inference <- obs_inference_df %>% pull(y)
x_obs_inference <- obs_inference_df %>% pull(x)


# only the unseen observations
obs_unseen_df <- obs_df %>% slice(-(1:length(y_obs_inference)))

y_obs_unseen <- obs_unseen_df %>% pull(y)
x_obs_unseen <- obs_unseen_df %>% pull(x)


y_fit_draws_subset_comp_pool <- y_fit_draws_comp_pool[, , max_days*(patient_col-1) + x_obs_unseen] #test slice
y_fit_draws_subset_no_pool <- y_fit_draws_no_pool[, , max_days*(patient_col-1) + x_obs_unseen] #test slice
y_fit_draws_subset_part_pool <- y_fit_draws_part_pool[, , max_days*(patient_col-1) + x_obs_unseen] #test slice

# extract sigmas
sigma_1_draws <- posterior_draws1[, "sigma_1"]
sigma_2_draws <- posterior_draws1[, "sigma_2"]



# for (s in 1:s_stop){ #or 1:no_ts_max
#   log_lik[s][patient] = lognormal_lpdf(y1_data[patient, s] | log(y_fit[s][3]), sigma_1 + y_fit[s][3]*sigma_2) do this in R


#test on seen data
y_fit_draws_test_subset <- y_fit_draws_no_pool[, , max_days*(patient_col-1) + x_obs_unseen] #test slice
log_lik <- dlnorm(
  x = y_obs_unseen,
  meanlog = log(y_fit_draws_test_subset), #log or not "NOT"
  sdlog = sigma_1_draws + y_fit_draws_test_subset* sigma_2_draws,
  log = TRUE
)

pred_medians <- apply(y_fit_draws_test_subset, 2, median)
print(pred_medians)

# print(mean(log_lik[,1]))

log_lik[is.nan(log_lik)] <- -10
# 
log_lik[log_lik == -Inf] <- -10
# 
log_lik[log_lik < -10] <- -10

loo_estimate_test <- loo(log_lik, is_method="tis")

print(loo_estimate_test)

# loo_results_no_pool <- data.frame(loo_estimate_no_pool$estimate[1,])



### now loop ###

df_row_counter <-1

# creat empty data frames
df <- data.frame(matrix(ncol = 5, nrow = 0))
df <- setNames(df, c("x","y","ymin","ymax", "patient_ID"))


for (patient_col in 1:number_of_patients){
  
  print(sprintf("Patient: %d", patient_IDs[patient_col]))
  
  
  # All observed data for that patient
  obs_df <- data.frame(
    x = data_day_all[0:ts_lengths_all[patient_col] + 1, patient_col],
    y = data_PSA_all[0:ts_lengths_all[patient_col] + 1, patient_col]
  )
  
  y_obs <- obs_df %>% pull(y)
  x_obs <- obs_df %>% pull(x)
  
  
  # Observed data for interference for that patient
  obs_inference_df <- data.frame(
    x = data_day[0:ts_lengths[patient_col] + 1, patient_col],
    y = data_PSA[0:ts_lengths[patient_col] + 1, patient_col]
  )
  
  y_obs_inference <- obs_inference_df %>% pull(y)
  x_obs_inference <- obs_inference_df %>% pull(x)
  
  
  # only the unseen observations
  obs_unseen_df <- obs_df %>% slice(-(1:length(y_obs_inference)))
  
  y_obs_unseen <- obs_unseen_df %>% pull(y)
  x_obs_unseen <- obs_unseen_df %>% pull(x)
  
  
  y_fit_draws_subset_comp_pool <- y_fit_draws_comp_pool[, , max_days*(patient_col-1) + x_obs_unseen] #test slice
  y_fit_draws_subset_no_pool <- y_fit_draws_no_pool[, , max_days*(patient_col-1) + x_obs_unseen] #test slice
  y_fit_draws_subset_part_pool <- y_fit_draws_part_pool[, , max_days*(patient_col-1) + x_obs_unseen] #test slice
  
  # extract sigmas
  sigma_1_draws <- posterior_draws1[, "sigma_1"]
  sigma_2_draws <- posterior_draws1[, "sigma_2"]
  
  
  
  # for (s in 1:s_stop){ #or 1:no_ts_max
  #   log_lik[s][patient] = lognormal_lpdf(y1_data[patient, s] | log(y_fit[s][3]), sigma_1 + y_fit[s][3]*sigma_2) do this in R
  
  #comp pool
  
  print("comp pool")
  
  log_lik <- dlnorm(
    x = y_obs_unseen,
    meanlog = log(y_fit_draws_subset_comp_pool), #log or not
    sdlog = sigma_1_draws + y_fit_draws_subset_comp_pool * sigma_2_draws,
    log = TRUE
  )
  
  log_lik[is.nan(log_lik)] <- -10
  
  log_lik[log_lik == -Inf] <- -10
  
  log_lik[log_lik < -10] <- -10
  
  loo_estimate_comp_pool <- loo(log_lik, is_method="tis")
  print(loo_estimate_comp_pool)
  
  loo_results_comp_pool <- data.frame(loo_estimate_comp_pool$estimate[1,])
  
  
  
  
  df[df_row_counter,] <- data.frame(
    x = 1,  # Assign a single x value for both points
    y = loo_results_comp_pool[1,1],  # Mean estimate
    ymin = loo_results_comp_pool[1,1] - loo_results_comp_pool[2,1],  # Lower bound (mean - SE)
    ymax = loo_results_comp_pool[1,1] + loo_results_comp_pool[2,1],   # Upper bound (mean + SE)
    patient_ID = patient_IDs[patient_col]
  )
  df_row_counter <- df_row_counter + 1
  
  
  #no pool
  
  print("no pool")
  
  log_lik <- dlnorm(
    x = y_obs_unseen,
    meanlog = log(y_fit_draws_subset_no_pool), #log or not
    sdlog = sigma_1_draws + y_fit_draws_subset_no_pool * sigma_2_draws,
    log = TRUE
  )
  
  log_lik[is.nan(log_lik)] <- -10
  
  log_lik[log_lik == -Inf] <- -10
  
  log_lik[log_lik < -10] <- -10
  
  loo_estimate_no_pool <- loo(log_lik, is_method="tis")
  print(loo_estimate_no_pool)
  
  loo_results_no_pool <- data.frame(loo_estimate_no_pool$estimate[1,])
  
  
  df[df_row_counter,] <- data.frame(
    x = 2,  # Assign a single x value for both points
    y = loo_results_no_pool[1,1],  # Mean estimate
    ymin = loo_results_no_pool[1,1] - loo_results_no_pool[2,1],  # Lower bound (mean - SE)
    ymax = loo_results_no_pool[1,1] + loo_results_no_pool[2,1],   # Upper bound (mean + SE)
    patient_ID = patient_IDs[patient_col]
  )
  
  df_row_counter <- df_row_counter + 1
  
  #part pool
  
  print("part pool")
  
  log_lik <- dlnorm(
    x = y_obs_unseen,
    meanlog = log(y_fit_draws_subset_part_pool), #log or not "NOT"
    sdlog = sigma_1_draws + y_fit_draws_subset_part_pool * sigma_2_draws,
    log = TRUE
  )
  
  
  
  log_lik[is.nan(log_lik)] <- -10
  
  log_lik[log_lik == -Inf] <- -10
  
  log_lik[log_lik < -10] <- -10
  
  loo_estimate_part_pool <- loo(log_lik, is_method="tis")
  print(loo_estimate_part_pool)
  
  loo_results_part_pool <- data.frame(loo_estimate_part_pool$estimate[1,])
  
  df[df_row_counter,] <- data.frame(
    x = 3,  # Assign a single x value for both points
    y = loo_results_part_pool[1,1],  # Mean estimate
    ymin = loo_results_part_pool[1,1] - loo_results_part_pool[2,1],  # Lower bound (mean - SE)
    ymax = loo_results_part_pool[1,1] + loo_results_part_pool[2,1],   # Upper bound (mean + SE)
    patient_ID = patient_IDs[patient_col]
  )
  
  
  # update counters
  df_row_counter <- df_row_counter + 1
  
  patient_col <- patient_col + 1
  
  print("####################")
  
}

df$x <- factor(df$x, levels = c(1, 2, 3), labels = c("Complete", "No", "Partial"))


ggplot(df[df$patient_ID == 14, ], aes(x = x, y = y, ymin = ymin, ymax = ymax)) +
  geom_errorbar(width = 0.1) +
  geom_point(size = 3) +
  geom_pointrange() +
  geom_text(aes(label = patient_ID), vjust = -0, hjust = -2) +
  theme_minimal() +
  labs(x = "Pooling Type", y = "ELPD") +  # Optional: Update x-axis label
  theme(axis.text.x = element_text(angle = 0, hjust = 1))  # Rotate for readability if needed

# Save the data frame to a CSV file
write.csv(df, sprintf("ELPD_LOO_tis_scores_%s_%d.csv", experiment, number_of_patients), row.names = TRUE)

print(df[, c(1, 2,5)])  # prints the 1 2 5 columns


# ppc_dens_overlay(y = y_obs_unseen, yrep = y_fit_draws_subset_part_pool)
