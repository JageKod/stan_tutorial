
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


date <- '2025-07-21'

experiment <- 'first_cycle'

number_of_patients <- 10
no_reps <- number_of_patients


######## create data_in ############
inputdatafile_PSA <- sprintf('%s_%d_PSA_df.csv', experiment, number_of_patients)
inputdatafile_day <- sprintf('%s_%d_day_df.csv', experiment, number_of_patients)
inputdatafile_Tx <- sprintf('%s_%d_Tx_df.csv', experiment, number_of_patients)

# Read data from files:
data_PSA <- read.table(inputdatafile_PSA, header = TRUE, sep = ",", stringsAsFactors = FALSE)
data_day <- read.table(inputdatafile_day, header = TRUE, sep = ",", stringsAsFactors = FALSE)
data_Tx <- read.table(inputdatafile_Tx, header = TRUE, sep = ",", stringsAsFactors = FALSE)


patient_IDs <- c(1, 4, 13, 14, 15, 17, 24, 26, 28, 29, 30, 36, 37,39, 44, 50, 55, 58, 60, 61)

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
posterior_draws <- as.matrix(fit1, pars = fit1@model_pars)

gq_fit <- gqs(mod1, draws = posterior_draws, data = data_in)



#### Plot all patients in one PDF #######

if (no_reps>10){

  patient_IDs = c(1, 4, 6, 7, 8, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 24, 25, 26, 28, 29,
    30, 31, 33, 36, 37, 39, 40, 41, 42, 44,
    48, 49, 50, 52, 54, 55, 58, 60, 61, 62,
    63, 66, 71, 73, 75, 77, 78, 79, 80, 81,
    84, 86, 87, 88, 90, 91, 93, 94, 95, 96,
    97, 99, 100, 101, 102, 104, 105, 106, 108, 109)
    # Open PDF
    pdf(sprintf("output_grid_%s_%d_comp_pool_%s.pdf", experiment, number_of_patients, date), width = 21, height = 84)

    plot_list <- list()  # List to store plots

    for (patient_col in 1:no_reps){

      # Get Stan-generated predictions:
      pred1 <- as.data.frame(gq_fit,
                             pars = grep(sprintf("y_fit_tot\\[.*,%d\\]", patient_col),
                                         gq_fit@sim$fnames_oi,
                                         value = TRUE)) %>%
        gather(factor_key = TRUE) %>%  # Reshape data frame
        group_by(key) %>%
        summarize(lb = quantile(value, probs = 0.1),
                  lb_25 = quantile(value, probs = 0.25),
                  median = quantile(value, probs = 0.5),
                  ub_75 = quantile(value, probs = 0.75),
                  ub = quantile(value, probs = 0.9))
      
      pred1_subset <- pred1[1:ts_data_all[ts_lengths_all[patient_col], patient_col], ]
      pred1_subset$key <- droplevels(pred1_subset$key)
      

      tsplot_cp <- ggplot(pred1_subset, aes(x = as.numeric(key), y = median)) +
        geom_vline(xintercept = ts_data[ts_lengths[patient_col],patient_col], linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_ribbon(aes(ymin = lb, ymax = ub), fill = "blue", alpha = 0.25) +
        geom_line(color = "blue") + #median line
        geom_line(aes(y = lb), color = "blue") +  # Lower bound line
        geom_line(aes(y = ub), color = "blue")  # Upper bound line

      # Put data in df for plotting:
      plot_df <- data.frame(
        x = data_day_all[0:ts_lengths_all[patient_col] + 1, patient_col],
        y = data_PSA_all[0:ts_lengths_all[patient_col] + 1, patient_col]
      )


      # Plot data as points:
      tsplot_cp <- tsplot_cp +
        geom_point(data = plot_df, aes(x = x, y = y), show.legend = FALSE, size = 1, color = "black")

      # Add labels:
      tsplot_cp <- tsplot_cp +
        labs(title = sprintf("PSA prediction for patient %d", patient_IDs[patient_col]),
             subtitle = sprintf("Complete pooling %s_%d", experiment, number_of_patients),
             y = "PSA (ug/L)",
             x = "Time (days)")


      # Store the plot in the list
      plot_list[[patient_col]] <- tsplot_cp
    }

    # Arrange all plots in a 2-row, 5-column grid
    grid.arrange(grobs = plot_list, nrow = 14, ncol = 5)

    # Close PDF
    dev.off()
} else {
  # Open PDF
  pdf(sprintf("output_grid_%s_%d_comp_pool_%s.pdf", experiment, number_of_patients, date), width = 21, height = 12)
  
  plot_list <- list()  # List to store plots
  
  for (patient_col in 1:no_reps){
  
    # Get Stan-generated predictions:
    pred1 <- as.data.frame(gq_fit,
                           pars = grep(sprintf("y_fit_tot\\[.*,%d\\]", patient_col),
                                       gq_fit@sim$fnames_oi,
                                       value = TRUE)) %>%
      gather(factor_key = TRUE) %>%  # Reshape data frame
      group_by(key) %>%
      summarize(lb = quantile(value, probs = 0.1),
                median = quantile(value, probs = 0.5),
                ub = quantile(value, probs = 0.9))
  
    # Plot Stan-generated predictions as lines (median) and bands (confidence intervals):
    tsplot_cp <- ggplot(pred1) +
      # geom_rect(aes(xmin = 0, xmax = ts_data[ts_lengths[patient_col],patient_col], ymin = -Inf, ymax = Inf),
      # fill = "gray", alpha = 0.2) +
      geom_vline(xintercept = ts_data[ts_lengths[patient_col],patient_col], linetype = "dashed", color = "black", linewidth = 0.5) +
      geom_ribbon(aes(t_gd, ymin = lb, ymax = ub), color = "blue", fill = "blue", alpha = 0.25) +
      geom_line(aes(t_gd, median), color = "blue", linewidth = 1)
  
    # Put data in df for plotting:
    plot_df <- data.frame(
      x = data_day_all[0:ts_lengths_all[patient_col] + 1, patient_col],
      y = data_PSA_all[0:ts_lengths_all[patient_col] + 1, patient_col]
    )
  
  
    # Plot data as points:
    tsplot_cp <- tsplot_cp +
      geom_point(data = plot_df, aes(x = x, y = y), show.legend = FALSE, size = 1, color = "black")
  
    # Add labels:
    tsplot_cp <- tsplot_cp +
      labs(title = sprintf("PSA prediction for patient %d", patient_IDs[patient_col]),
           subtitle = sprintf("Complete pooling %s_%d", experiment, number_of_patients),
           y = "PSA (ug/L)",
           x = "Time (days)")
  
  
    # Store the plot in the list
    plot_list[[patient_col]] <- tsplot_cp
  }
  
  # Arrange all plots in a 2-row, 5-column grid
  grid.arrange(grobs = plot_list, nrow = 2, ncol = 5)
  
  # Close PDF
  dev.off()
}



############## plot specific patient ###############
patient_col <- 1 # column corresponding to specific patient to plot

patient_ID <- patient_IDs[patient_col] #Patient ID in this column

pdf(sprintf("output_patient_%d_%s_%d_comp_pool_%s.pdf", patient_ID, experiment, number_of_patients, date))

# Get Stan-generated predictions:
pred1 <- as.data.frame(gq_fit, 
                       pars = grep(sprintf("y_fit_tot\\[.*,%d\\]", patient_col), 
                                   gq_fit@sim$fnames_oi, 
                                   value = TRUE)) %>%
  gather(factor_key = TRUE) %>%  # Reshape data frame
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.1),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.9)) 

# Plot Stan-generated predictions as lines (median) and bands (confidence intervals):
tsplot_cp <- ggplot(pred1) +
  geom_vline(xintercept = ts_data[ts_lengths[patient_col],patient_col], linetype = "dashed", color = "black", linewidth = 0.5) + 
  geom_ribbon(aes(t_gd, ymin = lb, ymax = ub), color = "blue", fill="blue", alpha = 0.25,) +
  geom_line(aes(t_gd, median), color = "blue", linewidth=1)

# Put data in df for plotting:
plot_df <- data.frame(
  x = data_day_all[0:ts_lengths_all[patient_col] + 1, patient_col],
  y = data_PSA_all[0:ts_lengths_all[patient_col] + 1, patient_col]
)

# Plot data as points:
tsplot_cp <- tsplot_cp +
  geom_point(data=plot_df, aes(x=x, y=y),show.legend = FALSE, size=3,color = "black")
# geom_point(data=inputdata_y1, aes(x=Day_1, y=populationsize),show.legend = FALSE, size=3,color = "blue")

# Add labels:
tsplot_cp <- tsplot_cp +
  labs(title=sprintf("PSA prediction for patient %d", patient_ID), subtitle=sprintf("output_%s_%d_comp_pool_%s", experiment, number_of_patients, date), y="PSA (ug/L)", x="time (days)", caption="")


plot(tsplot_cp)

dev.off()

############## plot specific patient ###############
patient_col <- 6 # column corresponding to specific patient to plot

patient_ID <- patient_IDs[patient_col] #Patient ID in this column

pdf(sprintf("output_patient_%d_%s_%d_comp_pool_%s.pdf", patient_ID, experiment, number_of_patients, date))

# Get Stan-generated predictions:
pred1 <- as.data.frame(gq_fit, 
                       pars = grep(sprintf("y_fit_tot\\[.*,%d\\]", patient_col), 
                                   gq_fit@sim$fnames_oi, 
                                   value = TRUE)) %>%
  gather(factor_key = TRUE) %>%  # Reshape data frame
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.1),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.9)) 

# Plot Stan-generated predictions as lines (median) and bands (confidence intervals):
tsplot_cp <- ggplot(pred1) +
  geom_vline(xintercept = ts_data[ts_lengths[patient_col],patient_col], linetype = "dashed", color = "black", linewidth = 0.5) + 
  geom_ribbon(aes(t_gd, ymin = lb, ymax = ub), color = "blue", fill="blue", alpha = 0.25,) +
  geom_line(aes(t_gd, median), color = "blue", linewidth=1)

# Put data in df for plotting:
plot_df <- data.frame(
  x = data_day_all[0:ts_lengths_all[patient_col] + 1, patient_col],
  y = data_PSA_all[0:ts_lengths_all[patient_col] + 1, patient_col]
)

# Plot data as points:
tsplot_cp <- tsplot_cp +
  geom_point(data=plot_df, aes(x=x, y=y),show.legend = FALSE, size=3,color = "black")
# geom_point(data=inputdata_y1, aes(x=Day_1, y=populationsize),show.legend = FALSE, size=3,color = "blue")

# Add labels:
tsplot_cp <- tsplot_cp +
  labs(title=sprintf("PSA prediction for patient %d", patient_ID), subtitle=sprintf("output_%s_%d_comp_pool_%s", experiment, number_of_patients, date), y="PSA (ug/L)", x="time (days)", caption="")


plot(tsplot_cp)

dev.off()

