# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)


# Parameters to load the wanted stan-fit
date <- '2025-04-17'

experiment <- 'all_but_last_cycle'

number_of_patients <- 10

# Maximum x-values, edit to make nice plots
ps_max <- 0.2
lambda_max <- 2
alpha_max <- 0.2
rho_max <- 0.003
varphi_max <- 0.2


#read fit
fit1 <- readRDS(sprintf("fit_%s_%d_no_pool_%s.rds",experiment, number_of_patients, date))



############## Posterior distribiutions ##########

# Extract posterior samples
posterior_samples <- as.data.frame(fit1)


# long format for ggplot
long_data <- posterior_samples %>%
  select(starts_with("p1["), starts_with("p2["), starts_with("p3["), starts_with("p4["), starts_with("p5[")) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")%>%
    mutate(group = sub("\\[.*", "", parameter))# Extract group names (e.g., "p1")



# Define correct group order
group_order <- c("p_s", "lambda", "alpha", "rho", "varphi")

# Map group names to desired labels
long_data <- long_data %>%
  mutate(group = recode(group,
                        "p1" = "p_s",
                        "p2" = "lambda",
                        "p3" = "alpha",
                        "p4" = "rho",
                        "p5" = "varphi")) %>%
  mutate(group = factor(group, levels = group_order))  # Force correct order




# Define prior distributions
prior_params <- list(
  "p_s" = list(mean = 0.0278, sd = 0.1, min_x = 0, max_x = ps_max),  
  "lambda" = list(mean = 0.693, sd = 1, min_x = 0, max_x = lambda_max),
  "alpha" = list(mean = 0.036, sd = 0.1, min_x = 0, max_x = alpha_max),
  "rho" = list(mean = 0.000187, sd = 0.001, min_x = 0, max_x = rho_max),
  "varphi" = list(mean = 0.0856, sd = 1, min_x = 0, max_x = varphi_max)
)

long_data <- long_data %>% filter(!(group == "varphi" & value > varphi_max))

prior_means <- data.frame(
  group = factor(names(prior_params), levels = group_order),  # Explicitly set levels
  mean_value = sapply(prior_params, function(x) x$mean)
)



# Create a dataset for prior distributions
prior_params <- do.call(rbind, lapply(names(prior_params), function(g) {
  x_vals <- seq(prior_params[[g]]$min_x, prior_params[[g]]$max_x, length.out = 100)
  data.frame(
    group = factor(g, levels = group_order),  # Ensure correct order
    x = x_vals,
    y = dnorm(x_vals, mean = prior_params[[g]]$mean, sd = prior_params[[g]]$sd)
  )
}))



# Create the plot
ggplot(long_data, aes(x = value, color = parameter, fill = parameter)) +
  geom_density(alpha = 0.1) +

  facet_wrap(~ group, scales = "free") +  
  coord_cartesian() + 
  theme_minimal() +
  labs(title = sprintf("Posterior Distributions, No Pooling, %s %d patients", experiment, number_of_patients), x = "Value", y = "Density") +
  theme(legend.position = "none", axis.text.y = element_blank()) +
  geom_line(data = prior_params, aes(x = x, y = y * 3, group = group),  # Multiply y by a factor to adjust scale
            color = "black", linetype = "dashed", linewidth=0.75, inherit.aes = FALSE) +
  geom_vline(data = prior_means, aes(xintercept = mean_value),
             color = "black", linetype = "dotted", linewidth = 1)


############################
########## plot sigma ####
#########################



# Add sigma_1 and sigma_2 to sigma_data
sigma_data <- posterior_samples %>%
  select(sigma_1, sigma_2) %>%
  pivot_longer(cols = c(sigma_1, sigma_2), names_to = "parameter", values_to = "value") %>%
  mutate(group = "sigma")  # Assign a group name


# Maximum x-values, edit to make nice plots
sigma_1_max <- 2
sigma_2_max <- 2


# Define prior distributions
sigma_prior_params <- list(
  "sigma_1" = list(mean = 0, sd = 1, min_x = 0, max_x = sigma_1_max),
  "sigma_2" = list(mean = 0, sd = 1, min_x = 0, max_x = sigma_2_max)
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
  # x_position = c(1,0.01)  # use this if simulated data
  x_position = c(0,0)
)

# Plot the posterior density for each parameter alongside prior distributions
ggplot(sigma_data, aes(x = value, fill = parameter, color = parameter)) +
  geom_density(alpha = 0.1) +
  facet_wrap(~ parameter, scales = "free", ncol = 2) +  # Adjust ncol for desired number of columns
  geom_line(data = sigma_prior_data, aes(x = x, y = y*5), color = "black", linetype = "dashed", size = 0.75) +  # Add prior distributions as black dotted lines
  geom_vline(data = sigma_line_positions, aes(xintercept = x_position), linetype = "dotted", color = "black", linewidth = 1) +
  theme_minimal() +
  labs(title = sprintf("Posterior Distributions, No Pooling, sigma, %s %d patients", experiment, number_of_patients), x = "Value", y = "Density") +
  theme(legend.position = "none", axis.text.y = element_blank()) +  # Optional: Remove the legend
  theme(strip.text = element_text(size = 10)) +  # Adjust the facet label size
  coord_cartesian()  # Allow axes to adjust for all parameters, or use specific limits if desired


########################################
########## plot both in same figure ####
########################################

# Combine sigma_1 and sigma_2 under one group for shared subplot
sigma_long_data <- posterior_samples %>%
  select(sigma_1, sigma_2) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(group = "sigma_1(left), sigma_2 (right)")

sigma_prior_data <- do.call(rbind, lapply(names(sigma_prior_params), function(param) {
  x_vals <- seq(sigma_prior_params[[param]]$min_x, sigma_prior_params[[param]]$max_x, length.out = 100)
  data.frame(
    group = "sigma_1(left), sigma_2 (right)",  # group both sigmas together
    parameter = param,
    x = x_vals,
    y = dnorm(x_vals, mean = sigma_prior_params[[param]]$mean, sd = sigma_prior_params[[param]]$sd)
  )
}))

sigma_prior_means <- data.frame(
  group = "sigma_1(left), sigma_2 (right)",
  parameter = names(sigma_prior_params),
  mean_value = sapply(sigma_prior_params, function(x) x$mean)
)


# Combine all posterior
all_long_data <- bind_rows(long_data, sigma_long_data)

# Updated group order
all_group_order <- c("p_s", "lambda", "alpha", "rho", "varphi", "sigma_1(left), sigma_2 (right)")
all_long_data$group <- factor(all_long_data$group, levels = all_group_order)

# Combine prior data
all_prior_data <- bind_rows(prior_params, sigma_prior_data)
all_prior_data$group <- factor(all_prior_data$group, levels = all_group_order)

# Combine prior means
all_prior_means <- bind_rows(prior_means, sigma_prior_means)
all_prior_means$group <- factor(all_prior_means$group, levels = all_group_order)

ggplot(all_long_data, aes(x = value, color = parameter, fill = parameter)) +
  geom_density(alpha = 0.1) +
  facet_wrap(~ group, scales = "free", ncol = 3) +
  theme_minimal() +
  labs(title = sprintf("Posterior and Prior Distributions, No pooling,  %s %d patients", experiment, number_of_patients),
       x = "Value", y = "Density") +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12),
    strip.text = element_text(size = 12),
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12),  
    legend.position = "none"
    ) +
  geom_line(data = all_prior_data,
            aes(x = x, y = y * 3, group = parameter),  # Keep different lines for sigma_1 and sigma_2
            color = "black", linetype = "dashed", linewidth = 0.75, inherit.aes = FALSE) +
  geom_vline(data = all_prior_means,
             aes(xintercept = mean_value),
             color = "black", linetype = "dotted", linewidth = 1)