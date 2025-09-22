# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstan)


# Parameters to load the wanted stan-fit
date <- '2025-04-17'

experiment <- 'all_but_last_cycle'

number_of_patients <- 10

# Maximum x-values, edit to make nice plots
ps_max <- 0.075
lambda_max <- 1
alpha_max <- 0.2
rho_max <- 0.003
varphi_max <- 0.5


#read fit
fit1 <- readRDS(sprintf("fit_%s_%d_part_pool_%s.rds",experiment, number_of_patients, date))



# Extract parameter names
param_names <- names(fit1)

# Identify parameters to exclude
exclude_params <- grep("^(ts_array|y_fit|y_fit_tot|log_lik)", param_names, value = TRUE)


# Subset fit1 before converting to data frame
posterior_samples <- as.data.frame(fit1, pars = setdiff(param_names, exclude_params))



############# group posteriors #################

# Reshape to long format and exclude unwanted parameters
mean_data <- posterior_samples %>%
  pivot_longer(cols = starts_with("p"), 
               names_to = "parameter", 
               values_to = "value") %>%
  filter(grepl("_mean$", parameter)) %>%  # Keep only parameters ending in "_mean" or "_std"
  mutate(group = sub("\\[.*", "", parameter))  # Extract group names (e.g., "p1")
# 

# Filter data to include only parameters that end with "_mean" or "sigma"
mean_data <- mean_data %>%
  filter(grepl("_mean$", parameter)) %>%  # Include sigma
  mutate(parameter = recode(parameter,
                            "p1_mean" = "p_s_mean",
                            "p2_mean" = "lambda_mean",
                            "p3_mean" = "alpha_mean",
                            "p4_mean" = "rho_mean",
                            "p5_mean" = "varphi_mean"))


# Manually set the order of the parameters
mean_data$parameter <- factor(mean_data$parameter, levels = c("p_s_mean", "lambda_mean", "alpha_mean", "rho_mean", "varphi_mean"))





# Define prior distributions
prior_params <- list(
  "p_s_mean" = list(mean = 0.0278, sd = 0.1, min_x = 0, max_x = ps_max),
  "lambda_mean" = list(mean = 0.69, sd = 1, min_x = 0, max_x = lambda_max),
  "alpha_mean" = list(mean = 0.036, sd = 0.1, min_x = 0, max_x = alpha_max),
  "rho_mean" = list(mean = 0.000187, sd = 0.001, min_x = 0, max_x = rho_max),
  "varphi_mean" = list(mean = 0.0856, sd = 1, min_x = 0, max_x = varphi_max)

)

# Create a dataset for prior distributions
prior_data <- do.call(rbind, lapply(names(prior_params), function(g) {
  x_vals <- seq(prior_params[[g]]$min_x, prior_params[[g]]$max_x, length.out = 100)
  data.frame(
    parameter = factor(g, levels = c("p_s_mean", "lambda_mean", "alpha_mean", "rho_mean", "varphi_mean")),  # Ensure correct order
    x = x_vals,
    y = dnorm(x_vals, mean = prior_params[[g]]$mean, sd = prior_params[[g]]$sd)
  )
}))



#plot vertical lines at point values
line_positions <- data.frame(
  parameter = unique(mean_data$parameter),
  x_position = c(0.0278, 0.69, 0.036, 0.000187, 0.0856)  # Specify the x-axis positions for the vertical lines in each plot
)



# Plot the posterior density for each parameter alongside prior distributions
ggplot(mean_data, aes(x = value, fill = parameter, color = parameter)) +
  geom_density(alpha = 0.1) +
  facet_wrap(~ parameter, scales = "free", ncol = 3) +  # Adjust ncol for desired number of columns
  geom_line(data = prior_data, aes(x = x, y = y*5), color = "black", linetype = "dashed", size = 0.75) +  # Add prior distributions as black dotted lines
  geom_vline(data = line_positions, aes(xintercept = x_position), linetype = "dotted", color = "black", linewidth = 1) +
  theme_minimal() +
  labs(title = sprintf("Posterior Distributions, Partial Pooling, Mean, %s %d patients", experiment, number_of_patients), x = "Value", y = "Density") +
  theme(legend.position = "none", axis.text.y = element_blank()) +  # Optional: Remove the legend
  theme(strip.text = element_text(size = 10)) +  # Adjust the facet label size
  coord_cartesian()  # Allow axes to adjust for all parameters, or use specific limits if desired




####### Patient specific posteriors #########


# Gather the samples into a long format for ggplot
patient_specific <- posterior_samples %>%
  select(starts_with("p1["), starts_with("p2["), starts_with("p3["), starts_with("p4["), starts_with("p5[")) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")%>%
  mutate(group = sub("\\[.*", "", parameter))# Extract group names (e.g., "p1")


# Define correct group order
group_order <- c("p_s", "lambda", "alpha", "rho", "varphi")

# Map group names to desired labels
patient_specific <- patient_specific %>%
  mutate(group = recode(group,
                        "p1" = "p_s",
                        "p2" = "lambda",
                        "p3" = "alpha",
                        "p4" = "rho",
                        "p5" = "varphi")) %>%
  mutate(group = factor(group, levels = group_order)) %>% # Force correct order
  filter(!(group == "varphi" & value > 0.2)) %>%
  filter(!(group == "rho" & value > 0.003)) #%>%
  # filter(!(group == "p_s" & value > 0.05))




ggplot(patient_specific, aes(x = value, color = parameter, fill = parameter)) +
  geom_density(alpha = 0.1) +
  facet_wrap(~ group, scales = "free") +  
  coord_cartesian() + 
  theme_minimal() +
  labs(title = sprintf("Posterior Distributions, Partial Pooling, Patient Specific, %s %d patients", experiment, number_of_patients), x = "Value", y = "Density") +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12),
    strip.text = element_text(size = 12),
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12),  
    legend.position = "none"
  )



#####################################################



###############################
#--------------- plot sigma ----
##############################


# Add sigma_1 and sigma_2 to the long_data
sigma_data <- posterior_samples %>%
  select(sigma_1, sigma_2) %>%
  pivot_longer(cols = c(sigma_1, sigma_2), names_to = "parameter", values_to = "value") %>%
  mutate(group = "sigma")  # Assign a group name


# Maximum x-values, edit to make nice plots
sigma_1_max <- 1
sigma_2_max <- 1


# Define prior distributions
sigma_prior_params <- list(
  "sigma_1" = list(mean = 0, sd = 1, min_x = 0, max_x = 2),
  "sigma_2" = list(mean = 0, sd = 0.1, min_x = 0, max_x = 1)
)

# Create a dataset for prior distributions
sigma_prior_data <- do.call(rbind, lapply(names(sigma_prior_params), function(g) {
  x_vals <- seq(sigma_prior_params[[g]]$min_x, sigma_prior_params[[g]]$max_x, length.out = 100)
  data.frame(
    parameter = factor(g, levels = c("sigma_1", "sigma_2")),  # Ensure correct order
    x = x_vals,
    y = dnorm(x_vals, mean = sigma_prior_params[[g]]$mean, sd = sigma_prior_params[[g]]$sd)
  )
}))

#plot vertical lines at point values
sigma_line_positions <- data.frame(
  parameter = unique(sigma_data$parameter),
  # x_position = c(1,0.01)  # if simualted data
  x_position = c(0,0) 
)

# Plot the posterior density for each parameter alongside prior distributions
ggplot(sigma_data, aes(x = value, fill = parameter, color = parameter)) +
  geom_density(alpha = 0.1) +
  facet_wrap(~ parameter, scales = "free", ncol = 2) +  # Adjust ncol for desired number of columns
  geom_line(data = sigma_prior_data, aes(x = x, y = y*5), color = "black", linetype = "dashed", size = 0.75) +  # Add prior distributions as black dotted lines
  geom_vline(data = sigma_line_positions, aes(xintercept = x_position), linetype = "dotted", color = "black", linewidth = 1) +
  theme_minimal() +
  labs(title = sprintf("Posterior Distributions, Partial Pooling, sigma, %s %d patients", experiment, number_of_patients), x = "Value", y = "Density") +
  theme(legend.position = "none", axis.text.y = element_blank()) +  # Optional: Remove the legend
  theme(strip.text = element_text(size = 10)) +  # Adjust the facet label size
  coord_cartesian()  # Allow axes to adjust for all parameters, or use specific limits if desired




#################################
#----------------plot std--------
################################


# Reshape to long format and exclude unwanted parameters
std_data <- posterior_samples %>%
  pivot_longer(cols = starts_with("p"), 
               names_to = "parameter", 
               values_to = "value") %>%
  filter(grepl("_std$", parameter)) %>%  # Keep only parameters ending in "_std"
  mutate(parameter = recode(parameter,
                            "p1_std" = "p_s_std",
                            "p2_std" = "lambda_std",
                            "p3_std" = "alpha_std",
                            "p4_std" = "rho_std",
                            "p5_std" = "varphi_std"))



# Manually set the order of the parameters
std_data$parameter <- factor(std_data$parameter, levels = c("p_s_std", "lambda_std", "alpha_std", "rho_std", "varphi_std"))


# Maximum x-values, edit to make nice plots
ps_max_std <- 0.1
lambda_max_std <- 1
alpha_max_std <- 0.15
rho_max_std <- 0.0015
varphi_max_std <- 0.4

# Define prior distributions
prior_params <- list(
  "p_s_std" = list(mean = 0, sd = 0.1, min_x = 0, max_x = ps_max_std),
  "lambda_std" = list(mean = 0, sd = 1, min_x = 0, max_x = lambda_max_std),
  "alpha_std" = list(mean = 0, sd = 0.1, min_x = 0, max_x = alpha_max_std),
  "rho_std" = list(mean = 0, sd = 0.001, min_x = 0, max_x = rho_max_std),
  "varphi_std" = list(mean = 0, sd = 1, min_x = 0, max_x = varphi_max_std)
)

# Create a dataset for prior distributions
prior_data <- do.call(rbind, lapply(names(prior_params), function(g) {
  x_vals <- seq(prior_params[[g]]$min_x, prior_params[[g]]$max_x, length.out = 100)
  data.frame(
    parameter = factor(g, levels = c("p_s_std", "lambda_std", "alpha_std", "rho_std", "varphi_std")),  # Ensure correct order
    x = x_vals,
    y = dnorm(x_vals, mean = prior_params[[g]]$mean, sd = prior_params[[g]]$sd)
  )
}))

#plot vertical lines at point values
std_line_positions <- data.frame(
  parameter = unique(std_data$parameter),
  x_position = c(0.01, 0.1, 0.01, 0.0001, 0.01)  # Specify the x-axis positions for the vertical lines in each plot
)

# Plot the posterior density for each parameter alongside prior distributions
ggplot(std_data, aes(x = value, fill = parameter, color = parameter)) +
  geom_density(alpha = 0.4, color = NA) +
  facet_wrap(~ parameter, scales = "free", ncol = 3) +  # Adjust ncol for desired number of columns
  geom_line(data = prior_data, aes(x = x, y = y*5), color = "black", linetype = "dashed", size = 0.75) +  # Add prior distributions as black dotted lines
  # geom_vline(data = std_line_positions, aes(xintercept = x_position), linetype = "dotted", color = "red", linewidth = 1) + #uncomment if simulated data
  theme_minimal() +
  labs(title = sprintf("Posterior and Prior Distributions, Partial Pooling, Standard deviation, %s %d patients", experiment, number_of_patients), x = "Value", y = "Density") +
  theme(legend.position = "none", axis.text.y = element_blank()) +  # Optional: Remove the legend
  theme(strip.text = element_text(size = 10)) +  # Adjust the facet label size
  coord_cartesian()  # Allow axes to adjust for all parameters, or use specific limits if desired


########################################
########## plot mean and sigma in same figure ####
########################################


long_data <- posterior_samples %>%
  pivot_longer(
    cols = c(p1_mean, p2_mean, p3_mean, p4_mean, p5_mean),
    names_to = "parameter",
    values_to = "value")


long_data <- long_data %>%
  select(parameter, value)

long_data <- long_data %>%
  mutate(parameter = case_when(
    parameter == "p1_mean" ~ "p_s_mean",
    parameter == "p2_mean" ~ "lambda_mean",
    parameter == "p3_mean" ~ "alpha_mean",
    parameter == "p4_mean" ~ "rho_mean",
    parameter == "p5_mean" ~ "varphi_mean"))

long_data <- long_data %>% mutate(group=parameter)


# Combine sigma_1 and sigma_2 under one group for shared subplot
sigma_long_data <- posterior_samples %>%
  select(sigma_1, sigma_2) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(group = "sigma_1 (left), sigma_2 (right)")


# Combine all posterior
all_long_data <- bind_rows(long_data, sigma_long_data)

# Updated group order
all_group_order <- c("p_s_mean", "lambda_mean", "alpha_mean", "rho_mean", "varphi_mean", "sigma_1 (left), sigma_2 (right)")
all_long_data$group <- factor(all_long_data$group, levels = all_group_order)



# Define prior distributions
prior_params <- list(
  "p_s_mean" = list(mean = 0.0278, sd = 0.1, min_x = 0, max_x = ps_max),  
  "lambda_mean" = list(mean = 0.693, sd = 1, min_x = 0, max_x = lambda_max),
  "alpha_mean" = list(mean = 0.036, sd = 0.1, min_x = 0, max_x = alpha_max),
  "rho_mean" = list(mean = 0.000187, sd = 0.001, min_x = 0, max_x = rho_max),
  "varphi_mean" = list(mean = 0.0856, sd = 1, min_x = 0, max_x = varphi_max),
  "sigma_1 (left), sigma_2 (right)" = list(mean = 0, sd = 1, min_x = 0, max_x = sigma_1_max)
)

prior_means <- data.frame(
  group = factor(names(prior_params), levels = all_group_order),  # Explicitly set levels
  mean_value = sapply(prior_params, function(x) x$mean)
)


# Create a dataset for prior distributions
prior_data <- do.call(rbind, lapply(names(prior_params), function(g) {
  x_vals <- seq(prior_params[[g]]$min_x, prior_params[[g]]$max_x, length.out = 100)
  data.frame(
    group = factor(g, levels = all_group_order),  # Ensure correct order
    x = x_vals,
    y = dnorm(x_vals, mean = prior_params[[g]]$mean, sd = prior_params[[g]]$sd)
  )
}))

ggplot(all_long_data, aes(x = value, color = parameter, fill = parameter)) +
  geom_density(alpha = 0.4, color = NA) +
  facet_wrap(~ group, scales = "free", ncol = 3) +
  theme_minimal() +
  labs(title = sprintf("Posterior and Prior Distributions, Partial Pooling, Group Parameters, %s %d patients", experiment, number_of_patients),
       x = "Value", y = "Density") +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12),
    strip.text = element_text(size = 12),
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12),  
    legend.position = "none"
  ) +
  geom_line(data = prior_data,
            aes(x = x, y = y * 3),  # Keep different lines for sigma_1 and sigma_2
            color = "black", linetype = "dashed", linewidth = 0.75, inherit.aes = FALSE) +
  geom_vline(data = prior_means,
             aes(xintercept = mean_value),
             color = "black", linetype = "dotted", linewidth = 1)

