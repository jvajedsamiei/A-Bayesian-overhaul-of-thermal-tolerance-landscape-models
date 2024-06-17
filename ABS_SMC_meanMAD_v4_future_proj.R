##### This script conduct ABC-SMC and the following analyses. Written by Jahangir Vajedsamiei.

## libraries ----
#library(tidyr)
#library(survival)
#library(transport)
#library(openxlsx)
# Load the required libraries
# Install necessary packages
#install.packages("tidyverse")
#install.packages("survival")
#install.packages("zoo")

# Load required libraries
library(tidyverse)
library(dplyr)
library(lubridate)
library(survival)
library(zoo)
library(gridExtra)
library(scales)
library(ggridges)
library(ggpubr)
library(cowplot)
library(modelr)
library(tidybayes)
library(readxl)
library(readr)
library(openxlsx)
library(ggplot2)
library(writexl)
library(brms)
library(officer)
library(flextable)
library(patchwork) # New library for plot_layout

##### ABC-SMC ----
## Addresses for row experimental data and output ----
main_address = paste0("~/Library/Mobile Documents/com~apple~CloudDocs/Jahangir Air/Postdoc_Position:Propos:RelatedStuff/PROPOSAL_postdoc/Selection_docs/Niklas_Mytilus_experiment/")
# Output for Constant Heatwave (CHW) experiment
dir.create(file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "CHW"), recursive = TRUE, showWarnings = FALSE)
address_outputs_CHW <- file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "CHW") 
# Output for Dynamic Heatwave (DHW) experiment
dir.create(file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "DHW"), recursive = TRUE, showWarnings = FALSE)
address_outputs_DHW <- file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "DHW")

## Load your survival data collected in the experiment(s) ----
observed_surv_df <- read_xlsx(file.path(main_address, "obsSurvival_Temperature_dfs/observed_surv_df.xlsx"))

## Generate data (vectorised partly) ----
generate_data <- function(sample_names, proposal, temperature_df, T_ref, observed_surv_df, LTT) {
  # Initialize an empty vector to store Mean Absolute Deviations (MADs) for each sample
  MADs <- sapply(sample_names, function(s_n) {
    # Extracting proposal parameters for the current iteration
    prop_meanlog10LBR_refT <- proposal[1]
    prop_sdlog10LBR_refT <- proposal[2]
    prop_k <- proposal[3]
    
    # Filter temperature data for the current sample
    temperature_sample_df <- subset(temperature_df, sample_name == s_n) # subset temperature data for the sample (treatment level)
    
    # Number of individuals in the population for simulation
    N_ <- 100
    
    # Simulate Lethality Buildup rate (LBR) for each individual
    T_LBR_pop_df <- lapply(1:N_, function(ind) {
      repeat {
        # Generate LBR at reference temperature, retrying in case of NA values
        LBR_T_ref <- rlnorm(1, meanlog = log(10^prop_meanlog10LBR_refT), sdlog = log(10^prop_sdlog10LBR_refT))
        if (!is.na(LBR_T_ref)) {
          break
        }
      }
      
      # Assigning temperature readings for simulation
      Temp <- temperature_sample_df$Temperature
      
      # Calculate LBR at different temperatures, setting LBR to 0 below lethal temperature threshold (LTT)
      LBR_T <- ifelse(Temp <= LTT, 0, LBR_T_ref * (10^(prop_k * (Temp - T_ref))))
      counter <- seq_along(Temp)
      # Cumulative LBR, capped at 100
      LB <- cumsum(LBR_T)
      LB[LB > 100] <- 100
      
      # Create a data frame for each individual with their respective LBR and LB values
      T_LBR_ind_df <- data.frame(t_h = counter, Temp = Temp, ind = ind, LBR = LBR_T, LB = LB)
      return(T_LBR_ind_df)
    })
    
    # Combine individual data frames into a single (population) data frame
    T_LBR_pop_df <- do.call(rbind, T_LBR_pop_df)
    T_LBR_pop_df <- as.data.frame(T_LBR_pop_df)
    colnames(T_LBR_pop_df) <- c('t_h', 'Temp', 'ind', 'LBR', 'LB')
    
    # Dropping the LBR column as it's not needed beyond this point
    T_LBR_pop_df_sub <- subset(T_LBR_pop_df, select = -c(LBR))
    
    # Transform the data from long to wide format, with individuals as columns
    T_LBR_pop_df_wide <- T_LBR_pop_df_sub %>% group_by(t_h, Temp) %>% tidyr::pivot_wider(names_from = ind, values_from = LB)
    T_LBR_surv_pop_df_wide <- as.data.frame(T_LBR_pop_df_wide)
    T_LBR_surv_pop_df_wide$cum_mortality_sim <- rowSums(T_LBR_pop_df_wide[-(1:2)] == 100)
    T_LBR_surv_pop_df_wide$surv_p_sim <- 1 - (T_LBR_surv_pop_df_wide$cum_mortality_sim / N_)
    
    # Attach datetime to the wide data frame
    date_h <- temperature_sample_df[['date_h']]
    T_LBR_surv_pop_df_wide <- cbind(T_LBR_surv_pop_df_wide, date_h)
    T_LBR_surv_pop_df_wide$date_h <- as.POSIXct(T_LBR_surv_pop_df_wide$date_h, format = "%Y-%m-%d %H:%M")
    
    # Perform Kaplan-Meier survival analysis based on observed data (obs_KM_surv_df will have hourly data)
    obs_KM_surv_df <- observed_survival_analysis(temperature_sample_df, observed_surv_df, s_n)
    
    # Merge simulated LB and survival data with observed and KM predicted survival data
    sim_obs_full_df <- full_join(T_LBR_surv_pop_df_wide, obs_KM_surv_df, by = c("date_h"), multiple = "all")
    
    # Filter for mid-day data to match the timing of daily observations
    sim_obs_mid_day_df <- sim_obs_full_df %>% filter(format(date_h, format = "%H:%M:%S") == "12:00:00")
    
    # Calculate MAD distance
    MAD <- mean(abs(sim_obs_mid_day_df$surv_p_sim - sim_obs_mid_day_df$surv_p_fit_interpol))
    MAD
  })
  
  # Calculate mean MAD
  mean_MAD <- mean(MADs)
  return(list(mean_MAD = mean_MAD))
}

## Function for Kaplan-Meier prediction based on observations ----
observed_survival_analysis <- function(temperature_sample_df, observed_surv_df, s_n) {
  # extracting datetime with h (date_h) and cumulative hours (t_h) from temperature_sample_df (this df includes the whole experimental time)
  n <- nrow(temperature_sample_df)
  temperature_sample_df$t_h <- 1:n
  df_time <- temperature_sample_df[, c("date_h", "t_h")]  # experimental time df
  # observed survival data for the sample
  obs_surv_sample <- observed_surv_df[observed_surv_df$sample_name == s_n, ]
  obs_surv_sample1 <- merge(df_time, obs_surv_sample, by = "date_h", all = TRUE)
  obs_surv_sample2 <- obs_surv_sample1[complete.cases(obs_surv_sample1[, c("mortality_status", "cum_mortality", "surv_p")]), ]
  
  # Kaplan-Meier survival prediction
  kmf <- survfit(Surv(t_h, mortality_status) ~ 1, data = obs_surv_sample2, se.fit = F)
  KM_fit_df <- data.frame(t_h=kmf$time, surv_p_fit = kmf$surv)
  
  # extract minimum observed surv_p at 12:00 (only for days when mortality was recorded)
  obs_surv_sample2 = as.data.frame(obs_surv_sample2)
  obs_surv_sample2$day <- as.Date(obs_surv_sample2$date_h)
  obs_surv_sample3 <- aggregate(x = obs_surv_sample2[c("surv_p")], by = obs_surv_sample2[c("date_h", "t_h")],
                              FUN = function(x) {min_val <- min(x, na.rm = TRUE)})
  
  obs_KM_surv_df = left_join(obs_surv_sample3, KM_fit_df, by=c("t_h")) # merge observed and KM predicted survival rates
  
  obs_KM_surv_df <- left_join(df_time, obs_KM_surv_df, by = c("date_h", "t_h"), multiple="all")
  obs_KM_surv_df[1, c("surv_p", "surv_p_fit")] <- 1
  obs_KM_surv_df$surv_p_fit_interpol <- na.locf(obs_KM_surv_df$surv_p_fit, na.rm = FALSE)
  obs_KM_surv_df = as.data.frame(obs_KM_surv_df)
  return(obs_KM_surv_df)
}

## ABC-SMC function (version with automatic save of best posterior)----
# Define the Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC) function
abc_smc <- function(sample_names, n_particles, n_iterations, epsilon_schedule_mean, tolerance) {
  # Initialize lists to store particles and weights for each iteration
  particles <- list()
  weights <- list()
  epsilon_mean <- epsilon_schedule_mean[1]
  
  # Initialize particles: generate initial set of particles based on proposal distribution
  for (i in 1:n_particles) {
    repeat {
      repeat {
        # Generate a proposal for each particle
        proposal <- rnorm(3, startvalue, tolerance)
        # Ensure the second proposal parameter is positive
        if (proposal[2] > 0) {
          break
        }
      }
      # Generate data based on the proposal and calculate the Mean Absolute Deviation (MAD)
      MADs_list <- generate_data(sample_names, proposal, temperature_df, T_ref, observed_surv_df, LTT)
      # Accept the particle if the MAD is less than or equal to the threshold
      mean_MAD <- MADs_list$mean_MAD
      if (mean_MAD <= epsilon_mean) {
        particles[[i]] <- proposal
        weights[[i]] <- 1 / n_particles
        break
      }
    }
  }
  
  # Initialize variables to store the best results achieved during the iterations
  best_particles <- NULL
  best_weights <- NULL
  best_mean_MAD_list <- NULL
  best_iter <- NULL
  mean_MAD_list <- NULL
  
  # Iterate through the SMC stages to sequentially decrease the epsilon tolerance and refine particles
  for (iter in 2:n_iterations) {
    epsilon_mean <- epsilon_schedule_mean[iter]
    
    new_particles <- list()
    new_weights <- list()
    mean_MADs = list()
    
    # Define a threshold for the maximum number of tries per particle
    max_tries_per_particle <- 10 # This means that the particle with success probability < 0.1 is removed
    success_counter <- 0
    tries_per_iter <- 0
    
    # Print the progress of the current iteration
    print(paste0("Stage ", iter, " started for epsilon_mean = ", epsilon_mean))
    cat("Length of particles:", length(particles), "\n")
    cat("Length of weights:", length(weights), "\n")
    cat("Sum of weights:", sum(unlist(weights)), "\n")
    
    # Iterate over each particle to refine them based on the new epsilon tolerance
    for (i in 1:n_particles) {
      valid_particle_found <- FALSE
      try_counter <- 0
      while (!valid_particle_found & try_counter < max_tries_per_particle) {
        # Select a particle based on the weighted distribution of existing particles
        idx <- sample(1:length(particles), 1, prob = weights, replace=T)
        
        # Generate a new proposal based on the selected particle
        repeat {
          proposal <- rnorm(3, particles[[idx]], tolerance)
          if (proposal[2] > 0) {
            break
          }
        }
        
        # Calculate the MAD for the new proposal
        MADs_list <- generate_data(sample_names, proposal, temperature_df, T_ref, observed_surv_df, LTT)
        mean_MAD <- MADs_list$mean_MAD
        
        # Accept the new particle if its MAD is within the current epsilon tolerance
        if (mean_MAD <= epsilon_mean) {
          new_particles[[i]] <- proposal
          mean_MADs[[i]] = mean_MAD
          
          # Calculate the weight for the new particle
          # This weight is based on how close the new proposal is to the existing particles, adjusted by the specified tolerance.
          new_weight <- 0
          for (j in 1:length(particles)) {
            # The weight is calculated as the sum of the exponential of the negative half squared Euclidean distance 
            # between the new proposal and each of the existing particles, normalized by the tolerance.
            # This is essentially a measure of similarity to existing particles.
            new_weight <- new_weight + exp(-0.5 * (sum((proposal - particles[[j]])^2) / sum(tolerance^2))) * weights[[j]]
          }
          
          # Add the calculated weight to the list of new weights for the new particles.
          new_weights[[i]] <- new_weight
          # Increment the counter of successful particles in this iteration.
          success_counter <- success_counter + 1
          # Mark the current particle as successfully found and valid.
          valid_particle_found <- TRUE
        }
        new_particles <- new_particles[!sapply(new_particles, is.null)]
        
        try_counter <- try_counter + 1
        tries_per_iter = tries_per_iter + 1
      }
    }
    
    # Check the stopping criterion for the ABC-SMC algorithm
    success_rate <- success_counter / n_particles
    
    # The success rate is calculated as the ratio of successful particles (those that met the criteria) to the total number of particles.
    # It's a measure of how many particles in this iteration were deemed acceptable based on the defined criteria (e.g., mean_MAD <= epsilon_mean).
    # If the success rate (proportion of valid particles found) is below 20%, stop the iteration.
    # The threshold of 20% is a parameter choice that can be adjusted based on the specific requirements of the model and data.
    if (success_rate < 0.2) {
      # Printing the stage details for diagnostic purposes.
      print(paste0("Stage ", iter, "   epsilon_mean = ", epsilon_mean))
      print(paste0("Overall tries = ", tries_per_iter))
      print(paste0("Success rate = ", success_rate))
      cat("Stopping criterion reached at iteration", iter, "\n")

      # Capture the best set of particles and their weights before termination.
      # These are the parameters that will be returned as the result of the ABC-SMC process.
      best_particles <- particles
      best_weights <- weights
      best_mean_MAD_list <- mean_MAD_list
      best_iter <- iter
      break
    } else {
      # If the stopping criterion is not met, the algorithm proceeds to the next iteration.
      print(paste0("Stage ", iter, " finished for epsilon_mean = ", epsilon_mean))
      print(paste0("Overall tries = ", tries_per_iter))
      print(paste0("Success rate = ", success_rate))
      # Normalize the weights of the new particles.
      # Normalization ensures that the sum of all weights equals 1, which is essential for the proper functioning of the weighted sampling in the next iteration.
      new_weights <- unlist(new_weights)
      new_weights <- new_weights / sum(new_weights)
      
      # Update the particles and weights with the new values for the next iteration.
      particles <- new_particles
      weights <- new_weights
      mean_MAD_list = mean_MADs
    }
  }
  # Return the results of the ABC-SMC process.
  # This includes the best set of particles, their associated weights, the mean MAD values, and the iteration number where the algorithm terminated.
  return(list(particles = best_particles, weights = best_weights, mean_MAD_list = best_mean_MAD_list, iteration = best_iter))
}
                

## Main script CHW----
temperature_df <- read_xlsx(file.path(main_address, "obsSurvival_Temperature_dfs/temperature_df.xlsx")) # for us 
sample_names <- c('26', '27', '28', '29') # You can adjust these sample names as needed
n_particles <- 5000
n_iterations <- 50
epsilon_schedule_mean <- round(seq(0.09, 0.04, length.out = n_iterations), 4)
# SD of proposal distributions
sd_meanlog10LBR_refT <- 0.016
sd_sdlog10LBR <- 0.012
sd_k <- 0.015
tolerance <- c(3*sd_meanlog10LBR_refT, 3*sd_sdlog10LBR, 5*sd_k) # You can adjust the tolerance values as needed
# Start values for the means of proposal distributions
E_meanlog10LBR_refT <- -0.45
E_sdlog10LBR <- 0.15
E_k <- 0.6
startvalue = c(E_meanlog10LBR_refT, E_sdlog10LBR, E_k)
# reference temperature
T_ref = 28
LTT = 22 # NOTE: In our CHW, this has no effect on the outcome (all temperatures were > 22°C)
# Run the ABC-SMC algorithm & save approximate posterior parameters with weights and MADs
system.time({
  results <- abc_smc(sample_names, n_particles, n_iterations, epsilon_schedule_mean, tolerance)
})
setwd(paste0(address_outputs_CHW))
saveRDS(results, "results_CHW.rds")

## Main script DHW----
temperature_df <- read_excel(paste0(main_address, "Temp_data/KOB2022_summer/NEW/df_GHL_long.xlsx")) # def_all_long was made via another script: Temp_data/KOB2022_summer/Temp_all.R
sample_names <- c('E2', 'A1', 'B1', 'F2') # DHW
n_particles <- 5000
n_iterations <- 50
epsilon_schedule_mean <- round(seq(0.09, 0.01, length.out = n_iterations), 4)
# SDs of proposal distributions
sd_meanlog10LBR_refT <- 0.016
sd_sdlog10LBR <- 0.012
sd_k <- 0.015
tolerance <- c(3*sd_meanlog10LBR_refT, 3*sd_sdlog10LBR, 5*sd_k) # You can adjust the tolerance values as needed
# Start values for the means of proposal distributions
E_meanlog10LBR_refT <- -0.45
E_sdlog10LBR <- 0.15
E_k <- 0.6
startvalue = c(E_meanlog10LBR_refT, E_sdlog10LBR, E_k)
# reference temperature
T_ref = 28
LTTs =  c(20) #c(20, 22, 23, 24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5) 
# Run the ABC-SMC algorithm & save approximate posterior parameters with weights and MADs
for(LTT in LTTs){
  print(LTT)
  system.time({
    results <- abc_smc(sample_names, n_particles, n_iterations, epsilon_schedule_mean, tolerance)
  })
  cat(paste0("!!!!! RESULT NOTE !!!!!    For LTT ", LTT, ":  min = ", round(min(unlist(results$mean_MAD_list)),3), ", mean = ", round(mean(unlist(results$mean_MAD_list)),3), " & SD = ", round(sd(unlist(results$mean_MAD_list)), 3)))
  setwd(paste0(address_outputs_DHW))
  saveRDS(results, paste0(LTT, "_results_DHW.rds")) }


## Define a theme and size for plots ----
theme_publication <- function(base_size = 11, base_family = "") {
  theme_light(base_size = base_size, base_family = base_family) +
    theme(text = element_text(size = base_size),
          plot.title = element_blank(),
          axis.title = element_text(size = rel(1.1)),
          axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)), # adjust top margin for x axis title
          axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)), # adjust right margin for y axis title
          axis.text = element_text(size = rel(1)),
          strip.text = element_text(size = rel(1.1)),
          legend.text = element_text(size = rel(1)),
          legend.title = element_text(size = rel(1.1)),
          strip.background = element_blank(),
          panel.spacing = unit(1, "lines"), # existing theme elements here...
          #legend.spacing = unit(0.2, 'lines'),
          legend.box.margin = margin(t = 0, r = 0, b = -7, l = 0, unit = "pt"),
          plot.margin = margin(t = 2, r = 10, b = 2, l = 8, unit = "mm") ) } # adjust the space between legend entries
plot_width_SingleColumn = 8.2 #cm
plot_width_TwoThirdPage = 11
plot_width_FullPage = 17.3

## Posterior distributions (dfs and plots) ----
# form particles_weights_MADs_df from the results.rds file for CHW and DHW
address_outputs <- file.path(main_address, "Outputs", "ABC_SMC_MAD", "CHW") 
setwd(paste0(address_outputs))
results <- readRDS("results_CHW.rds")
particles <- results$particles
weights <- results$weights
mean_MAD_list <- results$mean_MAD_list
mean_MAD = do.call(rbind, mean_MAD_list)
iteration <- results$iteration
particles_weights_MADs_df_CHW <- as.data.frame(do.call(rbind, particles))
colnames(particles_weights_MADs_df_CHW) <- c("meanlog10LBR_refT", "sdlog10LBR_refT", "k")
particles_weights_MADs_df_CHW$weight <- weights
particles_weights_MADs_df_CHW$mean_MAD <- mean_MAD
particles_weights_MADs_df_CHW$iteration <- iteration
length(particles_weights_MADs_df_CHW$meanlog10LBR_refT)

LTT=20
LTTs =  c(20) # 20 °C was the lowest temperature in the experiment. Therefore, in this case, LTT=20 means has a similar effect as not setting any LTT.
for(LTT in LTTs){
  result_LTT = paste0(LTT, "_results_DHW.rds")
  address_outputs <- file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "DHW")
  setwd(paste0(address_outputs))
  results <- readRDS(result_LTT)
  particles <- results$particles
  weights <- results$weights
  mean_MAD_list <- results$mean_MAD_list 
  mean_MAD = do.call(rbind, mean_MAD_list)
  iteration <- results$iteration
  particles_weights_MADs_df_DHW <- as.data.frame(do.call(rbind, particles))
  colnames(particles_weights_MADs_df_DHW) <- c("meanlog10LBR_refT", "sdlog10LBR_refT", "k")
  particles_weights_MADs_df_DHW$weight <- weights
  particles_weights_MADs_df_DHW$mean_MAD <- mean_MAD
  particles_weights_MADs_df_DHW$iteration <- iteration
  
  # Form posterior parameter dfs and compare posteriors between CHW & DHW
  sample_size = 10000
  set.seed(123)  # setting a seed for reproducibility
  particles_weights_MADs_df_CHW_rep <- particles_weights_MADs_df_CHW %>% 
    sample_n(size = sample_size, replace = TRUE, weight = weight)
  set.seed(123) 
  particles_weights_MADs_df_DHW_rep <- particles_weights_MADs_df_DHW %>% 
    sample_n(size = sample_size, replace = TRUE, weight = weight)
  

  
  # Combine the dfs and reshape them and plot
  particles_weights_MADs_df_DHW_rep$experiment <- "DHW"
  particles_weights_MADs_df_CHW_rep$experiment <- "CHW"
  combined_df <- rbind(particles_weights_MADs_df_DHW_rep, particles_weights_MADs_df_CHW_rep)
  combined_df_long <- pivot_longer(combined_df, cols = c(meanlog10LBR_refT, sdlog10LBR_refT, k), names_to = "Parameter", values_to = "para_value")
  
  write_xlsx(as.data.frame(combined_df_long), file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "PostParam_ABCSMC_CHW_DHW_df.xlsx"))
  
  head(combined_df_long)
  
  ## Convert TTL param from log10 to original scale
  # Function to convert log10 mean to original scale
  convert_mean_to_original <- function(log_mean) {
    return(10^log_mean) }
  
  # Function to estimate original scale standard deviation from log scale parameters
  convert_sd_to_original <- function(log_mean, log_sd) {
    mean_original <- convert_mean_to_original(log_mean)
    var_original <- (10^(2 * log_mean)) * ((10^(2 * log_sd)) - 1)
    sd_original <- sqrt(var_original)
    return(sd_original) }

  # Apply conversions row-wise
  combined_df_long$para_value_original_scale <- NA_real_ # Initialize column with NA
  
  # Loop through each row
  for (i in 1:nrow(combined_df_long)) {
    if (combined_df_long$Parameter[i] == "meanlog10LBR_refT") {
      combined_df_long$para_value_original_scale[i] <- convert_mean_to_original(combined_df_long$para_value[i])
    } else if (combined_df_long$Parameter[i] == "sdlog10LBR_refT") {
      # Find meanlog10LBR_refT value from the same iteration and experiment
      same_group <- combined_df_long$iteration[i] == combined_df_long$iteration & 
        combined_df_long$experiment[i] == combined_df_long$experiment &
        combined_df_long$Parameter == "meanlog10LBR_refT"
      log_mean <- mean(combined_df_long$para_value[same_group])
      combined_df_long$para_value_original_scale[i] <- convert_sd_to_original(log_mean, combined_df_long$para_value[i])
    } else if (combined_df_long$Parameter[i] == "k") {
      combined_df_long$para_value_original_scale[i] <- 10^combined_df_long$para_value[i]
    }
  }
  
  # View the modified dataframe
  combined_df_long = combined_df_long[,3:7]
  head(combined_df_long)
  
  ## Creating the summary stat table (with experiment and iteration kept) for posterior params of both log10 and original scale
  summary_stat_table <- combined_df_long %>%
    group_by(iteration, experiment, Parameter) %>%
    summarise(
      median_Log10_scale = median(para_value, na.rm = TRUE),
      Q05_Log10_scale = quantile(para_value, 0.05, na.rm = TRUE),
      Q95_Log10_scale = quantile(para_value, 0.95, na.rm = TRUE),
      median_original_scale = median(para_value_original_scale, na.rm = TRUE),
      Q05_original_scale = quantile(para_value_original_scale, 0.05, na.rm = TRUE),
      Q95_original_scale = quantile(para_value_original_scale, 0.95, na.rm = TRUE)
    ) %>%
    ungroup() # Remove grouping
  summary_stat_table[,4:9] = round(summary_stat_table[,4:9], 3)
  
  # Create a flextable object from the summary_stat_table
  my_flextable <- flextable(summary_stat_table)
  # Create a Word document object
  doc <- read_docx()
  # Add the flextable to the Word document
  doc <- doc %>% 
    body_add_flextable(my_flextable)
  # Define the path and save the Word document
  print(doc, target = file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "Summary_stat_tab_PostParam_bothScales_DHWttl20_CHW.docx"))

  ## Create the comparative posterior distribution figure
  # Create a new column with labels as character strings for the plot
  combined_df_long1 = combined_df_long
  combined_df_long1$CustomParameter <- ifelse(combined_df_long1$Parameter == "k", "k~or~log[10](dot(L))/dT",
                                              ifelse(combined_df_long1$Parameter == "meanlog10LBR_refT", "mean(log[10](dot(L)[r]))",
                                                     ifelse(combined_df_long1$Parameter == "sdlog10LBR_refT", "sd(log[10](dot(L)[r]))", combined_df_long1$Parameter)))
  p = ggplot(combined_df_long1, aes(x = para_value, y = experiment, fill = factor(after_stat(quantile)))) +
    stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, quantiles = c(0.025, 0.975)) +
    scale_fill_manual(name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
                      labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")) +
    labs(x = "Parameter value", y = "Experiment") +
    facet_grid(. ~ CustomParameter, scales = "free_x", labeller = label_parsed) +
    scale_x_continuous(breaks = seq(from = floor(min(combined_df_long1$para_value)), 
                                    to = ceiling(max(combined_df_long1$para_value)), 
                                    by = 0.1)) +
    theme_publication() + theme(strip.text = element_text(colour = "black"),
                                strip.background = element_rect(fill = NA, colour = NA),
                                legend.position = "top", legend.title = element_text(size = 11))
  dir.create(file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "Compare_posteriors_DHWvsCHW"), showWarnings = FALSE)
  address_outputs <- file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "Compare_posteriors_DHWvsCHW")
  setwd(paste0(address_outputs))
  ggsave(paste0("Compare_SMCABCposteriors_DHWLTT", LTT, "vsCHW1.pdf"), plot = p, width = plot_width_FullPage-2, height = 7, units = "cm") 
  
  # plot post param dist with original scale
  combined_df_long1$CustomParameter_original <- ifelse(combined_df_long1$Parameter == "k", "k~or~dot(L)/dT",
                                              ifelse(combined_df_long1$Parameter == "meanlog10LBR_refT", "mean(dot(L)[r])",
                                                     ifelse(combined_df_long1$Parameter == "sdlog10LBR_refT", "sd(dot(L)[r])", combined_df_long1$Parameter)))
  p = ggplot(combined_df_long1, aes(x = para_value_original_scale, y = experiment, fill = factor(after_stat(quantile)))) +
    stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, quantiles = c(0.025, 0.975)) +
    scale_fill_manual(name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
                      labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")) +
    labs(x = "Parameter value", y = "Experiment") +
    facet_grid(. ~ CustomParameter_original, scales = "free_x", labeller = label_parsed) +
    theme_publication() + #theme_light() +
    theme(strip.text = element_text(colour = "black"), strip.background = element_rect(fill = NA, colour = NA),
          legend.position = "top", legend.title = element_text(size = 11))
  dir.create(file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "Compare_posteriors_DHWvsCHW"), showWarnings = FALSE)
  address_outputs <- file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "Compare_posteriors_DHWvsCHW")
  setwd(paste0(address_outputs))
  ggsave(paste0("Compare_SMCABCposteriors_OriginalScale_DHWLTT", LTT, "vsCHW1.pdf"), plot = p, width = plot_width_FullPage-2, height = 7, units = "cm") 
}


### Compare Posteriors between LTTs (only log10 scale) ----
# In this section, we show that, in our case, setting LTTs only as high as >26°C shift the post distributions.But this need to be tested at each study. 
LTTs =  c(22, 23, 24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5) 
all_post_df = NULL
ix=1
for(LTT in LTTs){
  result_LTT = paste0(LTT, "_results_DHW.rds")
  address_outputs <- file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "DHW")
  setwd(paste0(address_outputs))
  results <- readRDS(result_LTT)
  particles <- results$particles
  weights <- results$weights
  mean_MAD_list <- results$mean_MAD_list 
  mean_MAD = do.call(rbind, mean_MAD_list)
  iteration <- results$iteration
  post_df_DHW <- as.data.frame(do.call(rbind, particles))
  colnames(post_df_DHW) <- c("meanlog10LBR_refT", "sdlog10LBR_refT", "k")
  post_df_DHW$weight <- weights
  post_df_DHW$mean_MAD <- mean_MAD
  post_df_DHW$iteration <- iteration
  
  # Form posterior parameter dfs and compare posteriors between CHW & DHW
  sample_size = 1000
  set.seed(123) 
  post_df_DHW_rep <- post_df_DHW %>% 
    sample_n(size = sample_size, replace = TRUE, weight = weight)
  print(head(post_df_DHW_rep))
  
  post_df_DHW_rep$Lethal_TT <- LTT
  
  all_post_df[[ix]] = post_df_DHW_rep
  ix = ix+1
}

all_post_df1 = do.call(rbind, all_post_df)
head(all_post_df1)

## Table
# Ensure that Lethal_TT is a factor
all_post_df1$Lethal_TT <- as.factor(all_post_df1$Lethal_TT)
# Calculate the summary statistics
summary_stats_table <- all_post_df1 %>%
  group_by(Lethal_TT) %>%
  summarize(Median_MAD = median(mean_MAD, na.rm = TRUE),
            Q5 = quantile(mean_MAD, probs = 0.05, na.rm = TRUE),
            Q95 = quantile(mean_MAD, probs = 0.95, na.rm = TRUE),
            .groups = 'drop')
# Create a flextable object from the summary_stat_table
my_flextable <- flextable(summary_stat_table)
# Create a Word document object
doc <- read_docx()
# Add the flextable to the Word document
doc <- doc %>% 
  body_add_flextable(my_flextable)
# Define the path and save the Word document
print(doc, target = file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "Summay_table_posteriors_DHWallLTTs.docx"))

## Plots Posteriors for LTTs
all_post_df_long <- pivot_longer(all_post_df1, cols = c(meanlog10LBR_refT, sdlog10LBR_refT, k), names_to = "Parameter", values_to = "para_value")
head(all_post_df_long)

write_xlsx(as.data.frame(all_post_df_long), file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "PostParam_ABCSMC_DHWallLTTs_df.xlsx"))

# Create a new column with labels as character strings
all_post_df_long1 = all_post_df_long
all_post_df_long1$Lethal_TT = as.factor(all_post_df_long1$Lethal_TT)
all_post_df_long1$CustomParameter <- ifelse(all_post_df_long1$Parameter == "k", "k~or~log[10](dot(L))/dT",
                                            ifelse(all_post_df_long1$Parameter == "meanlog10LBR_refT", "mean(log[10](dot(L)[r]))",
                                                   ifelse(all_post_df_long1$Parameter == "sdlog10LBR_refT", "sd(log[10](dot(L)[r]))", all_post_df_long1$Parameter)))
# Calculate the median para_value for Experiment 20 for each parameter
medians_exp20 <- all_post_df_long1 %>%
  filter(Lethal_TT == 20) %>%
  group_by(CustomParameter) %>%
  summarize(median_para_value = median(para_value))

# Create the plot with the median lines for Experiment 20 in each facet
p = ggplot(all_post_df_long1, aes(x = para_value, y = Lethal_TT, fill = factor(after_stat(quantile)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, quantiles = c(0.025, 0.975)) +
  scale_fill_manual(name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
                    labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")) +
  labs(x = "Parameter value", y = "Lethal temperature thresholds") +
  facet_grid(. ~ CustomParameter, scales = "free_x", labeller = label_parsed) +
  theme_publication() + # theme_light() +
  theme(strip.text = element_text(colour = "black"), strip.background = element_rect(fill = NA, colour = NA),
        legend.position = "top", legend.title = element_text(size = 11)) +
  geom_vline(data = medians_exp20, aes(xintercept = median_para_value), linetype = "dashed", color = "black")
ggsave(paste0("Compare_SMCABCposteriors_DHWallLTTs_vs_CHW1.pdf"), plot = p, width = plot_width_FullPage, height = 14, units = "cm") 


##### Posterior prediction plot of simulations versus observations ----
'Simulating population function (simulate_T_LBR_pop): The function simulates a population of organisms under given temperature conditions. The function takes a dataframe with the temperature data, parameters for a reference temperature, mean and standard deviation for LBR at the reference temperature, and a constant k. It computes LBR at different temperatures and accumulates the sum of LBR. It then generates a dataframe with columns for time, temperature, individual number, LBR and accumulated LBR. This dataframe is then processed to compute the simulated survival probability at different times.
Observed survival analysis function (observed_survival_analysis): This function takes the temperature data, observed survival data, and a sample name as input, and generates Kaplan-Meier survival probabilities at different times. The function also interpolates missing values in the survival probabilities and returns the processed dataframe.
Plotting the simulated vs observed data: The script then enters a loop where it generates simulated survival curves using the function defined above, for different parameter sets. Each set of parameters has a weight associated with it, stored in particles_weights_df. For each sample, the script simulates the survival curves for each set of parameters and plots all these curves together on a ggplot graph. On the same plot, the observed survival probabilities are plotted as points. The simulated and observed curves are colored differently for visual distinction.
Saving the plots to a PDF: Each plot generated for a sample is stored in a list. After all samples are processed, the list of plots is arranged in a grid and saved to a PDF file. The script uses the grid.arrange function from the gridExtra package to do this, specifying one column (i.e., the plots are arranged vertically). The height of each subplot is adjusted for optimal viewing.'


## Load your survival data collected in the experiment(s) ----
observed_surv_df <- read_excel(file.path(main_address, "obsSurvival_Temperature_dfs/observed_surv_df.xlsx"))
## Function to simulate a population of LBRs and LBs (MC) ----
simulate_T_LBR_pop <- function(temperature_sample_df, meanlog10LBR_refT, sdlog10LBR_refT, k, T_ref, LTT) {
  N_ <- 100
  T_LBR_pop_df <- lapply(1:N_, function(ind) {
    repeat {
      LBR_T_ref <- rlnorm(1, meanlog = log(10^meanlog10LBR_refT), sdlog = log(10^sdlog10LBR_refT))
      if (!is.na(LBR_T_ref)) {
        break
      }
    }
    Temp <- temperature_sample_df$Temperature
    # Use only if you prefer to set a Lethal Temperature Threshold below which LBR = 0.
    LBR_T <- ifelse(Temp <= LTT, 0, LBR_T_ref * (10^(k * (Temp - T_ref))))
    counter <- seq_along(Temp)
    LB <- cumsum(LBR_T)
    LB[LB > 100] <- 100
    T_LBR_ind_df <- data.frame(t_h = counter, Temp = Temp, ind = ind, LBR = LBR_T, LB = LB)
    return(T_LBR_ind_df)
  })
  
  T_LBR_pop_df <- do.call(rbind, T_LBR_pop_df)
  T_LBR_pop_df <- as.data.frame(T_LBR_pop_df)
  colnames(T_LBR_pop_df) <- c('t_h', 'Temp', 'ind', 'LBR', 'LB')
  T_LBR_pop_df_sub <- subset(T_LBR_pop_df, select = -c(LBR))
  
  T_LBR_pop_df_wide <- T_LBR_pop_df_sub %>% group_by(t_h, Temp) %>% tidyr::pivot_wider(names_from = ind, values_from = LB)
  T_LBR_surv_pop_df_wide <- as.data.frame(T_LBR_pop_df_wide) ###
  
  T_LBR_pop_df_wide$cum_mortality_sim <- rowSums(T_LBR_pop_df_wide[-(1:2)] == 100)
  T_LBR_pop_df_wide$surv_p_sim <- 1 - (T_LBR_pop_df_wide$cum_mortality_sim / N_)
  T_LBR_pop_df_wide$date_h <- temperature_sample_df[['date_h']]
  T_LBR_pop_df_wide = as.data.frame(T_LBR_pop_df_wide)
  
  return(T_LBR_pop_df_wide)
}

## Function for Kaplan-Meier prediction based on observations ----
observed_survival_analysis <- function(temperature_sample_df, observed_surv_df, s_n) {
  df_time <- temperature_sample_df[, c("date_h", "t_h")]  # experimental time df
  # observed survival data for the sample
  obs_surv_sample <- observed_surv_df[observed_surv_df$sample_name == s_n, ]
  obs_surv_sample1 <- merge(df_time, obs_surv_sample, by = "date_h", all = TRUE)
  obs_surv_sample2 <- obs_surv_sample1[complete.cases(obs_surv_sample1[, c("mortality_status", "cum_mortality", "surv_p")]), ]
  
  # Kaplan-Meier survival prediction
  kmf <- survfit(Surv(t_h, mortality_status) ~ 1, data = obs_surv_sample2, conf.int = 0.95, conf.type = "log-log")
  KM_fit_df <- data.frame(t_h=kmf$time, surv_p_fit = kmf$surv, surv_p_upper = kmf$lower, surv_p_lower = kmf$upper)
  
  # extract minimum observed surv_p at 12:00 (only for days when mortality was recorded)
  obs_surv_sample2 = as.data.frame(obs_surv_sample2)
  obs_surv_sample2$day <- as.Date(obs_surv_sample2$date_h)
  obs_surv_sample3 <- aggregate(x = obs_surv_sample2[c("surv_p")], by = obs_surv_sample2[c("date_h", "t_h")],
                                FUN = function(x) {min_val <- min(x, na.rm = TRUE)})
  
  obs_KM_surv_df = left_join(obs_surv_sample3, KM_fit_df, by=c("t_h")) # merge observed and KM predicted survival rates
  
  obs_KM_surv_df <- left_join(df_time, obs_KM_surv_df, by = c("date_h", "t_h"), multiple="all")
  obs_KM_surv_df[1, c("surv_p", "surv_p_fit", "surv_p_upper", "surv_p_lower")] <- 1
  obs_KM_surv_df$surv_p_fit_interpol <- na.locf(obs_KM_surv_df$surv_p_fit, na.rm = FALSE)
  obs_KM_surv_df$surv_p_upper_interpol <- na.locf(obs_KM_surv_df$surv_p_upper, na.rm = FALSE)
  obs_KM_surv_df$surv_p_lower_interpol <- na.locf(obs_KM_surv_df$surv_p_lower, na.rm = FALSE)
  obs_KM_surv_df = as.data.frame(obs_KM_surv_df)
  obs_KM_surv_df$Date <- as.Date(obs_KM_surv_df$date_h)
  obs_KM_surv_df$day = as.POSIXct(obs_KM_surv_df$Date)
  return(obs_KM_surv_df)
}

## Plot simulations versus observation for CHW or DHW (also calculate MAD to later compare with ABC-SMC) ----
plot_width_SingleColumn = 8.2 #cm
plot_width_TwoThirdPage = 11
plot_width_FullPage = 17.3

# for CHW----
theme_publication <- function(base_size = 11, base_family = "") {
  theme_light(base_size = base_size, base_family = base_family) +
    theme(text = element_text(size = base_size),
          axis.text.y = element_text(size = rel(1), color = "#0b3c5d"),
          axis.text.y.right = element_text(size = rel(1), color = "#CC3311"),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.y.right = element_blank(),
          axis.text.x = element_text(size = rel(1.1), angle = 45, hjust = 1),
          panel.spacing = unit(0.2, "lines"), # existing theme elements here...
          plot.margin = margin(t = 2, r = 6, b = 2, l = 6, unit = "mm") ) } # adjust the space between legend entries
temperature_df <- read_excel(file.path(main_address, "obsSurvival_Temperature_dfs/temperature_df.xlsx")) # for us 
T_ref = 28
LTT=20
particles_weights_df <- particles_weights_MADs_df_CHW
sample_names <- c('26', '27', '28', '29') # You can adjust these sample names as needed
system.time({
  setwd(address_outputs_CHW)
  plot_list_CHW <- list()
  MAD_list_CHW <- NULL
  counter <- 0
  ix = 1
  for (s_n in sample_names) {
    counter <- counter + 1
    temperature_sample_df = subset(temperature_df, sample_name == s_n)
    n <- nrow(temperature_sample_df)
    temperature_sample_df$t_h <- 1:n
    weekly_breaks <- seq(min(temperature_sample_df$date_h), max(temperature_sample_df$date_h), by = "week")
    obs_KM_surv_df <- observed_survival_analysis(temperature_sample_df, observed_surv_df, s_n)
    all_sim_surv_curves <- list()
    for (row_n in 1:1000) {
      sampled_row <- particles_weights_df[sample(nrow(particles_weights_df), 1, prob = particles_weights_df$weight), ]
      meanlog10LBR_refT <- sampled_row[1,1]
      sdlog10LBR_refT <- sampled_row[1, 2]
      k <- sampled_row[1, 3]
      T_LBR_pop_df_wide <- simulate_T_LBR_pop(temperature_sample_df, meanlog10LBR_refT, sdlog10LBR_refT, k, T_ref, LTT)
      T_LBR_pop_df_wide = as.data.frame(T_LBR_pop_df_wide)
      all_sim_surv_curves[[paste("meanlog10LBR_refT", meanlog10LBR_refT, "sdlog10LBR_refT", sdlog10LBR_refT, "k", k, sep = "_")]] <- T_LBR_pop_df_wide
      
      # Merge simulated LB and survival data with observed and KM predicted survival data
      sim_obs_full_df <- full_join(T_LBR_pop_df_wide, obs_KM_surv_df, by = c("date_h"), multiple = "all")
      head(T_LBR_pop_df_wide)
      head(obs_KM_surv_df)
      # Filter for mid-day data to match the timing of daily observations
      sim_obs_mid_day_df <- sim_obs_full_df %>% filter(format(date_h, format = "%H:%M:%S") == "12:00:00")
      MAD_df = data.frame(Parametrization = "ABC-SMC", Experiment = "CHW", sample_name = s_n,
                          meanlog10LBR_refT=meanlog10LBR_refT, sdlog10LBR_refT=sdlog10LBR_refT, k = k,
                          MAD = mean(abs(sim_obs_mid_day_df$surv_p_sim - sim_obs_mid_day_df$surv_p_fit_interpol))) 
      MAD_list_CHW[[ix]] = MAD_df
      ix = ix+1
    }
    p1 <- ggplot()
    for (curve_name in names(all_sim_surv_curves)) {
      p1 <- p1 + geom_line(data = all_sim_surv_curves[[curve_name]], aes(x = date_h, y = surv_p_sim), color = "#33BBEE", size = 0.1)
    }
    p1 <- p1 +
      geom_line(data = temperature_sample_df, aes(x = date_h, y = (Temperature - 20) / 10), size = 0.5, color = "#CC3311") +
      geom_point(data = obs_KM_surv_df, aes(x = day, y = surv_p), shape = 1, size = 1, color = "#0b3c5d") +
      geom_line(data = obs_KM_surv_df, aes(x = day, y = surv_p_fit_interpol), size = 0.7, color = "#0b3c5d") +
      geom_line(data = obs_KM_surv_df, aes(x = day, y = surv_p_upper_interpol), size = 0.7, linetype = "dotted", color = "#0b3c5d") +
      geom_line(data = obs_KM_surv_df, aes(x = day, y = surv_p_lower_interpol), size = 0.7, linetype = "dotted", color = "#0b3c5d") +
      
      scale_x_datetime(labels = if (counter == length(sample_names)) scales::date_format("%b %d") else NULL, breaks = scales::breaks_pretty(n = 6)) +
      scale_y_continuous(name = "Survival probability", limits = c(0, 1), sec.axis = sec_axis(~ . * 10 + 20, name = "Temperature [°C]",
                                                                                              labels = scales::number_format(accuracy = 1),
                                                                                              breaks = scales::pretty_breaks(n = 5)))
    plot_list_CHW[[counter]] <- p1
  }
  setwd(address_outputs_CHW)
  combined_plot_CHW <- plot_list_CHW[[1]] / plot_list_CHW[[2]] / plot_list_CHW[[3]] / plot_list_CHW[[4]] + 
    plot_layout(guides = 'collect') & theme_publication()
  ggsave("combined_plot_CHW.pdf", combined_plot_CHW, width = plot_width_SingleColumn, height = 18, units = "cm")
})
MAD_list_CHW = do.call(rbind, MAD_list_CHW)


# for DHW----
theme_publication <- function(base_size = 11, base_family = "") {
  theme_light(base_size = base_size, base_family = base_family) +
    theme(text = element_text(size = base_size),
          axis.text.y = element_text(size = rel(1), color = "#0b3c5d"),
          axis.text.y.right = element_text(size = rel(1), color = "#CC3311"),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.y.right = element_blank(),
          axis.text.x = element_text(size = rel(1.1), angle = 45, hjust = 1),
          panel.spacing = unit(0.2, "lines"), # existing theme elements here...
          plot.margin = margin(t = 2, r = 6, b = 2, l = 6, unit = "mm") ) } # adjust the space between legend entries
temperature_df <- read_excel(paste0(main_address, "/Temp_data/KOB2022_summer/NEW/df_GHL_long.xlsx")) # def_all_long was made via another script: Temp_data/KOB2022_summer/Temp_all.R
T_ref = 28
particles_weights_df <- particles_weights_MADs_df_DHW
sample_names <- c('E2', 'A1', 'B1', 'F2') # DHW
LTT=20
LTTs = c(20) #c(20, 22, 23, 24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5) 
for(LTT in LTTs){
  system.time({
    setwd(address_outputs_DHW)
    plot_list_DHW <- list()
    MAD_list_DHW <- NULL
    counter <- 0
    ix = 1
    for (s_n in sample_names) {
      counter <- counter + 1
      temperature_sample_df = subset(temperature_df, sample_name == s_n)
      n <- nrow(temperature_sample_df)
      temperature_sample_df$t_h <- 1:n
      weekly_breaks <- seq(min(temperature_sample_df$date_h), max(temperature_sample_df$date_h), by = "week")
      obs_KM_surv_df <- observed_survival_analysis(temperature_sample_df, observed_surv_df, s_n)
      all_sim_surv_curves <- list()
      for (row_n in 1:1000) {
        sampled_row <- particles_weights_df[sample(nrow(particles_weights_df), 1, prob = particles_weights_df$weight), ]
        meanlog10LBR_refT <- sampled_row[1,1]
        sdlog10LBR_refT <- sampled_row[1, 2]
        k <- sampled_row[1, 3]
        T_LBR_pop_df_wide <- simulate_T_LBR_pop(temperature_sample_df, meanlog10LBR_refT, sdlog10LBR_refT, k, T_ref, LTT)
        T_LBR_pop_df_wide = as.data.frame(T_LBR_pop_df_wide)
        all_sim_surv_curves[[paste("meanlog10LBR_refT", meanlog10LBR_refT, "sdlog10LBR_refT", sdlog10LBR_refT, "k", k, sep = "_")]] <- T_LBR_pop_df_wide
        
        # Merge simulated LB and survival data with observed and KM predicted survival data
        sim_obs_full_df <- full_join(T_LBR_pop_df_wide, obs_KM_surv_df, by = c("date_h"), multiple = "all")
        # Filter for mid-day data to match the timing of daily observations
        sim_obs_mid_day_df <- sim_obs_full_df %>% filter(format(date_h, format = "%H:%M:%S") == "12:00:00")
        MAD_df = data.frame(Parametrization = "ABC-SMC", Experiment = "DHW", sample_name = s_n,
                            meanlog10LBR_refT=meanlog10LBR_refT, sdlog10LBR_refT=sdlog10LBR_refT, k = k,
                            MAD = mean(abs(sim_obs_mid_day_df$surv_p_sim - sim_obs_mid_day_df$surv_p_fit_interpol))) 
        MAD_list_DHW[[ix]] = MAD_df
        ix = ix+1
      }
      p1 <- ggplot()
      for (curve_name in names(all_sim_surv_curves)) {
        p1 <- p1 + geom_line(data = all_sim_surv_curves[[curve_name]], aes(x = date_h, y = surv_p_sim), color = "#33BBEE", size = 0.1)
      }
      p1 <- p1 +
        geom_line(data = temperature_sample_df, aes(x = date_h, y = (Temperature - 20) / 10), size = 0.5, color = "#CC3311") +
        geom_point(data = obs_KM_surv_df, aes(x = day, y = surv_p), shape = 1, size = 1, color = "#0b3c5d") +
        geom_line(data = obs_KM_surv_df, aes(x = day, y = surv_p_fit_interpol), size = 0.7, color = "#0b3c5d") +
        geom_line(data = obs_KM_surv_df, aes(x = day, y = surv_p_upper_interpol), size = 0.7, linetype = "dotted", color = "#0b3c5d") +
        geom_line(data = obs_KM_surv_df, aes(x = day, y = surv_p_lower_interpol), size = 0.7, linetype = "dotted", color = "#0b3c5d") +
        
        scale_x_datetime(labels = if (counter == length(sample_names)) scales::date_format("%b %d") else NULL, breaks = scales::breaks_pretty(n = 8)) +
        scale_y_continuous(name = "Survival probability", limits = c(0, 1), sec.axis = sec_axis(~ . * 10 + 20, name = "Temperature [°C]",
                                                                                                labels = scales::number_format(accuracy = 1),
                                                                                                breaks = scales::pretty_breaks(n = 5)))
      plot_list_DHW[[counter]] <- p1
    }
    setwd(address_outputs_DHW)
    combined_plot_DHW <- plot_list_DHW[[1]] / plot_list_DHW[[2]] / plot_list_DHW[[3]] / plot_list_DHW[[4]] + 
      plot_layout(guides = 'collect') & theme_publication()
    ggsave("combined_plot_DHW.pdf", combined_plot_DHW, width = plot_width_SingleColumn, height = 18, units = "cm")
  })
  MAD_list_DHW = do.call(rbind, MAD_list_DHW)
}

# MAD df to use in the other script for comparison with for MAD from BR posteriors
MAD_ABCSMC_df = rbind(MAD_list_CHW, MAD_list_DHW)
write_xlsx(as.data.frame(MAD_ABCSMC_df), file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "MAD_ABCSMC_df.xlsx"))

## Combine posterior prediction plots
# Adjust theme setting for CHW and DHW plots with reduced left margin
theme_publication_CHW <- function(base_size = 11, base_family = "") {
  theme_light(base_size = base_size, base_family = base_family) +
    theme(text = element_text(size = base_size),
          axis.text.y = element_text(size = rel(1), color = "#0b3c5d"),
          axis.text.y.right = element_text(size = rel(1), color = "#CC3311"),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.y.right = element_blank(),
          axis.text.x = element_text(size = rel(1.1), angle = 45, hjust = 1),
          panel.spacing = unit(0.2, "lines"), # existing theme elements here...
          plot.margin = margin(t = 1, r = 3, b = 1, l = 11.8, unit = "mm")) # right margin reduced for CHW
}
theme_publication_DHW <- function(base_size = 11, base_family = "") {
  theme_light(base_size = base_size, base_family = base_family) +
    theme(text = element_text(size = base_size),
          axis.text.y = element_text(size = rel(1), color = "#0b3c5d"),
          axis.text.y.right = element_text(size = rel(1), color = "#CC3311"),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.y.right = element_blank(),
          axis.text.x = element_text(size = rel(1.1), angle = 45, hjust = 1),
          panel.spacing = unit(0.2, "lines"), # existing theme elements here...
          plot.margin = margin(t = 1, r = 11.8, b = 1, l = 3, unit = "mm")) # left margin reduced for DHW
}

# Apply the adjusted theme to both CHW and DHW plot lists before creating the combined plots
for (i in seq_along(plot_list_CHW)) {plot_list_CHW[[i]] <- plot_list_CHW[[i]] + theme_publication_CHW() }
for (i in seq_along(plot_list_DHW)) {plot_list_DHW[[i]] <- plot_list_DHW[[i]] + theme_publication_DHW() }

# Combine the plots, possibly with adjusted rel_widths to compensate for the x-axes text
combined_plot_CHWandDHW <- plot_grid(
  plot_grid(plotlist = plot_list_CHW, ncol = 1, rel_heights = c(1, 1, 1, 1.2725)), 
  plot_grid(plotlist = plot_list_DHW, ncol = 1, rel_heights = c(1, 1, 1, 1.2725)), 
  ncol = 2, rel_widths = c(1, 1) )
combined_plot_CHWandDHW + theme_publication()

# Add global y-axis titles using ggdraw and draw_label
final_plot <- ggdraw() +
  draw_plot(combined_plot_CHWandDHW) +
  draw_label("Survival probability", x = 0.04, y = 0.5, angle = 90, hjust = 0.5, vjust = 0.5, size = 12, color = "#0b3c5d") +
  draw_label("Temperature [°C]", x = 0.96, y = 0.5, angle = 90, hjust = 0.5, vjust = 0.5, size = 12, color = "#CC3311")
ggsave("combined_plot_CHWandDHW.pdf", final_plot, width = plot_width_FullPage-1, height = 18, units = "cm")






##### DHW experiment simulations (dfs) for LB, average_LB, and surv_p based on Posterior means ----
## preparations----
# Again load the posterior results from the .rds file (if not already loaded!)
theme_publication <- function(base_size = 11, base_family = "") {
  theme_light(base_size = base_size, base_family = base_family) +
    theme(text = element_text(size = base_size),
          #axis.text.y = element_text(size = rel(1), color = "#0b3c5d"),
          axis.text.y.right = element_text(size = rel(1), color = "#CC3311"),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.y.right = element_blank(),
          axis.text.x = element_text(size = rel(1.1), angle = 45, hjust = 1),
          panel.spacing = unit(0.2, "lines"), # existing theme elements here...
          plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm") ) } # adjust the space between legend entries
temperature_df <- read_excel(paste0(main_address, "/Temp_data/KOB2022_summer/NEW/df_GHL_long.xlsx")) # def_all_long was made via another script: Temp_data/KOB2022_summer/Temp_all.R
temperature_df <- temperature_df %>%
  arrange(Treatment, date_h)
head(temperature_df)
T_ref = 28
LTT = 15 # arbitrary chosen as very low, with no effect on the outcome

# Extract the most likely posterior parameter set (the one with the highest weight)
posteriors <- particles_weights_MADs_df_DHW[which.max(particles_weights_MADs_df_DHW$weight), ]

# address for saving outputs
dir.create(file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "DHW", "LB_surv_post_mean"), showWarnings = F, recursive = T)
address_outputs <- file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "DHW", "LB_surv_post_mean")
setwd(paste0(address_outputs))


## Function to simulate a population of LBRs and LBs (MC)----
simulate_T_LBR_pop <- function(temperature_sample_df, meanlog10LBR_refT, sdlog10LBR_refT, k, T_ref) {
  N_ <- 100
  T_LBR_pop_df <- lapply(1:N_, function(ind) {
    repeat {
      LBR_T_ref <- rlnorm(1, meanlog = log(10^meanlog10LBR_refT), sdlog = log(10^sdlog10LBR_refT))
      if (!is.na(LBR_T_ref)) {
        break
      }
    }
    Temp <- temperature_sample_df$Temperature
    # Use only if you prefer to set a Lethal Temperature Threshold below which LBR = 0.
    LBR_T <- ifelse(Temp <= LTT, 0, LBR_T_ref * (10^(k * (Temp - T_ref))))
    counter <- seq_along(Temp)
    LB <- cumsum(LBR_T /100)
    LB[LB > 1] <- 1
    T_LBR_ind_df <- data.frame(t_h = counter, Temp = Temp, ind = ind, LBR = LBR_T, LB = LB)
    return(T_LBR_ind_df)
  })
  T_LBR_pop_df <- do.call(rbind, T_LBR_pop_df)
  T_LBR_pop_df <- as.data.frame(T_LBR_pop_df)
  colnames(T_LBR_pop_df) <- c('t_h', 'Temp', 'ind', 'LBR', 'LB')
  T_LBR_pop_df_sub <- subset(T_LBR_pop_df, select = -c(LBR))
  T_LBR_pop_df_wide <- T_LBR_pop_df_sub %>% group_by(t_h, Temp) %>% tidyr::pivot_wider(names_from = ind, values_from = LB)
  T_LBR_surv_pop_df_wide <- as.data.frame(T_LBR_pop_df_wide) ###
  T_LBR_pop_df_wide$cum_mortality_sim <- rowSums(T_LBR_pop_df_wide[-(1:2)] == 1)
  T_LBR_pop_df_wide$surv_p_sim <- 1 - (T_LBR_pop_df_wide$cum_mortality_sim / N_)
  T_LBR_pop_df_wide$average_LB <- rowMeans(T_LBR_pop_df_wide[, 3:(2 + N_)])
  T_LBR_pop_df_wide$date_h <- temperature_sample_df[['date_h']]
  T_LBR_pop_df_wide = as.data.frame(T_LBR_pop_df_wide)
  return(T_LBR_pop_df_wide)
}

## Save excel sheet including posterior predictions of individual LB series and the population's average LB and surv and save the ggplot ----
system.time({
  setwd(paste0(address_outputs))
  plot_list <- list()
  counter <- 0
  sample_names = unique(temperature_df$sample_name)
  for (s_n in sample_names) {
    counter <- counter + 1
    temperature_sample_df = subset(temperature_df, sample_name == s_n)
    n <- nrow(temperature_sample_df)
    temperature_sample_df$t_h <- 1:n
    weekly_breaks <- seq(min(temperature_sample_df$date_h), max(temperature_sample_df$date_h), by = "2 weeks")
    all_sim_surv_curves <- list()
    for (row_n in 1:nrow(posteriors)) {
      meanlog10LBR_refT <- posteriors[row_n, 1]
      sdlog10LBR_refT <- posteriors[row_n, 2]
      k <- posteriors[row_n, 3]
      T_LBR_pop_df_wide <- simulate_T_LBR_pop(temperature_sample_df, meanlog10LBR_refT, sdlog10LBR_refT, k, T_ref)
      T_LBR_pop_df_wide = as.data.frame(T_LBR_pop_df_wide)
      # Add parameter values and sample name as new columns
      T_LBR_pop_df_wide$meanlog10LBR_refT <- meanlog10LBR_refT
      T_LBR_pop_df_wide$sdlog10LBR_refT <- sdlog10LBR_refT
      T_LBR_pop_df_wide$k <- k
      T_LBR_pop_df_wide$sample_name <- s_n
      all_sim_surv_curves[[paste("meanlog10LBR_refT", meanlog10LBR_refT, "sdlog10LBR_refT", sdlog10LBR_refT, "k", k, sep = "_")]] <- T_LBR_pop_df_wide
    }
    # Combine all data frames in the all_sim_surv_curves list into one data frame
    combined_df <- bind_rows(all_sim_surv_curves)
    # Save combined_df to a CSV file
    write.csv(combined_df, file = paste0("all_sim_LB_surv_df_", s_n, ".csv"), row.names = FALSE)

    p1 <- ggplot()
    N_ <- 100
    for (curve_name in names(all_sim_surv_curves)) {
      p1 <- p1 + geom_line(data = all_sim_surv_curves[[curve_name]], aes(x = date_h, y = surv_p_sim), color = "#33BBEE", size = 0.7)
      for (LB_col in 1:N_) {
        p1 <- p1 + geom_line(data = all_sim_surv_curves[[curve_name]], 
                             aes(x = date_h, y = !!sym(as.character(LB_col))), 
                             color = "#009E73", 
                             alpha = 0.5, 
                             size = 0.2) } }
    p1 <- p1 +
      geom_line(data = temperature_sample_df, aes(x = date_h, y = (Temperature - 15) / 15), size = 0.3, color = "#CC3311") +
      
      scale_x_datetime(labels = if (counter %in% c(10, 11, 12)) scales::date_format("%b %d") else NULL,
                       breaks = scales::breaks_pretty(n = 4) ) +
      scale_y_continuous(name = "", limits = c(0, 1),
                         sec.axis = sec_axis(~ . * 15 + 15, name = "",
                                             labels = scales::number_format(accuracy = 1),
                                             breaks = scales::pretty_breaks(n = 5)))
    # Add the plot to the list
    plot_list[[counter]] <- p1
  }
})

for (i in seq_along(plot_list)) {plot_list[[i]] <- plot_list[[i]] + theme_publication() }

combined_plot_LB_Surv_DHW <- plot_grid(plotlist = plot_list, ncol = 3,
                                       rel_heights = c(1, 1, 1, 1.2725,
                                                       1, 1, 1, 1.2725,
                                                       1, 1, 1, 1.2725))
# Add global y-axis titles using ggdraw and draw_label
final_plot <- ggdraw() +
  draw_plot(combined_plot_LB_Surv_DHW) + 
  draw_label("Survival probability", x = -0.012, y = 0.39, angle = 90, hjust = 0.5, vjust = 0.5, size = 12, color = "#33BBEE") +
    draw_label("Lethality buildup (scaled)", x = -0.012, y = 0.67, angle = 90, hjust = 0.5, vjust = 0.5, size = 12, color = "#009E73") +
  draw_label("Temperature [°C]", x = 1.012, y = 0.5, angle = 90, hjust = 0.5, vjust = 0.5, size = 12, color = "#CC3311") +
  theme(plot.margin = margin(t = 2, r = 9, b = 2, l = 9, unit = "mm") ) 
ggsave("posterior_means_LB_combined_ggplot.pdf", final_plot, width = plot_width_FullPage, height = 17, units = "cm")





## post-heatwave average population LBs vs. recruitment probability ----
# Set your working directory
setwd(file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "DHW", "LB_surv_post_mean")) 

# Get the list of all_sim_LB_surv_df (csv) files in the directory
file_list <- list.files(pattern = "*.csv")
post_HW_meanLBs_df <- data.frame()
# Loop over all csv files
for(file in file_list){
  # Read the csv file into a data frame
  df <- read_csv(file, show_col_types = FALSE)
  # Subset the columns you need
  df <- df %>% dplyr::select(average_LB, date_h)
  # Add a new column "sample_name" using the last two characters of the file name (without extension)
  df$sample_name <- substr(file, nchar(file)-5, nchar(file)-4)
  # Keep only the last row to have post-heatwave lethality buildup (the population's average value)
  df <- tail(df, 1)
  # Append the row to the final data frame
  post_HW_meanLBs_df <- rbind(post_HW_meanLBs_df, df)
}
tail(post_HW_meanLBs_df)

# Open recruitment abundance df
recruit_abun_df <- read_excel("~/Library/Mobile Documents/com~apple~CloudDocs/Jahangir Air/Postdoc_Position:Propos:RelatedStuff/PROPOSAL_postdoc/Selection_docs/Niklas_Mytilus_experiment/KOB22_Settlement Panels/Recruitment.xlsx")
recruit_abun_df = recruit_abun_df[, c("Tank", "Abundance")]
colnames(recruit_abun_df) = c("sample_name", "recruit_abu")

# Merge the dfs by "sample_name"
merged_df <- merge(post_HW_meanLBs_df, recruit_abun_df, by = "sample_name")
merged_df$log_average_LB = log(merged_df$average_LB)
merged_df$log_recruit_abu = log(merged_df$recruit_abu)

# Create the scatter plot
p <- ggplot(merged_df, aes(x = log_average_LB, y = log_recruit_abu)) +
  geom_point() +
  theme_light() +
  labs(x = "Average LB [log scale]", y = "Recruit abundance [log scale]") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11, angle = 0, hjust = 1))
ggsave("scatter_logScale_recruitment_vs_meanLB.pdf", plot = p, width = plot_width_SingleColumn, height = 6, units = "cm") 

# Remove rows where Treatment is 1 or 4.5
merged_df <- merged_df %>% filter(sample_name != "A1", sample_name != "B2")

# Fit the model
fit_nb3 <- brm(recruit_abu ~ s(log_average_LB, k=3), data = merged_df, family = negbinomial(link = "log"),
           chains = 4, cores = 4, iter = 10000, warmup = 5000,
           control = list(adapt_delta = 0.99999, stepsize = 0.001, max_treedepth = 20),
           seed = 123)
fit_nb3 <- add_criterion(fit_nb3, c("loo"))
fit_nb4 <- brm(recruit_abu ~ s(log_average_LB, k=4), data = merged_df, family = negbinomial(link = "log"),
              chains = 4, cores = 4, iter = 10000, warmup = 5000,
              control = list(adapt_delta = 0.99999, stepsize = 0.001, max_treedepth = 20),
              seed = 123)
fit_nb4 <- add_criterion(fit_nb4, c("loo"))
# choose best model
loo_df = as.data.frame(loo_compare(fit_nb3, fit_nb4))
sink(paste0("loo_df.txt"))
print(loo_df)
sink()

# Print a summary of the model fit
summary(fit_nb4)

# Plot posterior distributions
pdf(file = "Post_dist_check.pdf",  width = 6, height = 8, onefile=FALSE)
par(mfrow=c(1,1), mar = c(2,4,1,3.8), cex=1, mgp = c(1, 1, 0.2)) 
plot(fit_nb4)
dev.off()

# Plot posterior predictions
theme_publication <- function(base_size = 11, base_family = "") {
  theme_light(base_size = base_size, base_family = base_family) +
    theme(text = element_text(size = base_size),
          plot.title = element_blank(),
          axis.title = element_text(size = rel(1.1)),
          axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)), # adjust top margin for x axis title
          axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)), # adjust right margin for y axis title
          axis.text = element_text(size = rel(1)),
          strip.text = element_text(size = rel(1.1)),
          legend.text = element_text(size = rel(0.9)),
          legend.title = element_text(size = rel(1)),
          legend.position = "top",
          strip.background = element_blank(),
          #panel.spacing = unit(0.5, "lines"), # existing theme elements here...
          #legend.spacing = unit(0.2, 'lines'),
          #legend.box.margin = margin(t = 0, r = 0, b = -2, l = 0, unit = "pt"),
          plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm") ) }
plot_width_FullPage = 17.3
p_logX = merged_df %>%
  data_grid(log_average_LB = seq_range(log_average_LB, n = 50)) %>%
  add_predicted_draws(fit_nb4) %>%
  ggplot(aes(x = log_average_LB, y = recruit_abu)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95,.68)) + #, color = brewer.pal(5, "Blues")[[5]]) +
  geom_point(data = merged_df) +
  scale_fill_brewer() +
  labs(x = "Log(scaled average L)", y = "Recruit abundance") +
  scale_y_continuous(breaks = seq(from = 0, 
                                  to = 350, 
                                  by = 50))  +
  theme_publication()
ggsave("brms_recruitment_vs_logmeanLB.pdf", plot = p_logX, width = plot_width_FullPage/2.8, height = 7, units = "cm") 


p_originalX = merged_df %>%
  data_grid(log_average_LB = seq_range(log_average_LB, n = 50)) %>%
  add_predicted_draws(fit_nb4) %>%
  mutate(average_LB = exp(log_average_LB)) %>%  # Transform back to original scale
  ggplot(aes(x = average_LB, y = recruit_abu)) +  # Use the transformed scale for plotting
  stat_lineribbon(aes(y = .prediction), .width = c(.95,.68)) +
  geom_point(data = merged_df, aes(x = exp(log_average_LB))) +  # Transform points for plotting
  scale_fill_brewer() +
  labs(x = "Scaled average L", y = "Recruit abundance") +
  scale_y_continuous(breaks = seq(from = 0, 
                                  to = 350, 
                                  by = 50))  +
  theme_publication()
ggsave("brms_recruitment_vs_meanLB.pdf", plot = p_originalX, width = plot_width_FullPage/2.8, height = 7, units = "cm") 


'just to test how it works
# Generate data grid for log_average_LB corresponding to average_LB values from 0 to 0.25
log_values_range <- seq(log(0.001), log(0.006), length.out = 3)  # Adjusted for log scale

merged_df %>%
  data_grid(log_average_LB = log_values_range) %>%
  add_predicted_draws(fit_nb4) %>%
  mutate(average_LB = exp(log_average_LB)) %>% 
  ggplot(aes(x = average_LB, y = recruit_abu)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95,.68, .5)) +
  #geom_point(data = merged_df) +
  scale_fill_brewer() +
  labs(x = "Log(scaled average L)", y = "Recruit abundance") +
  theme_publication()


dfb = merged_df %>%
  data_grid(log_average_LB = log_values_range) %>%
  add_predicted_draws(fit_nb4) %>%
  mutate(average_LB = exp(log_average_LB)) %>%
  group_by(average_LB) %>%  # Group by average_LB
  summarise(median_prediction = median(.prediction, na.rm = TRUE)) '


## Posterior predictions of recruitment as a function of average_LB
Ave_LB_df <- tibble(average_LB = seq(0.001, 1, by = 0.001))
predic_df = Ave_LB_df %>%
  mutate(log_average_LB = log(average_LB)) %>%  # Log-transform for prediction
  add_predicted_draws(fit_nb)
tail(predic_df)

# Compute the quantiles of recruitment predicted for each average_LB
quantile_df <- predic_df %>% group_by(average_LB) %>%
  summarise(q_50 = quantile(.prediction, 0.5), q_05 = quantile(.prediction, 0.05), q_95 = quantile(.prediction, 0.95))

# Normalize the quantiles with respect to largest value of q_50 and create a new df
scaledrecruit_vs_meanLB <- quantile_df %>%
  mutate(average_LB = average_LB,
    scaled_recruit_05q = formatC(round((q_50 - min(q_50)) / (max(q_50) - min(q_50)), 3), format = "f", digits = 3),
    scaled_recruit_005q = formatC(round((q_05 - min(q_05)) / (max(q_50) - min(q_50)), 3), format = "f", digits = 3),
    scaled_recruit_095q = formatC(round((q_95 - min(q_95)) / (max(q_50) - min(q_50)), 3), format = "f", digits = 3)) %>%
  dplyr::select(average_LB, scaled_recruit_05q, scaled_recruit_005q, scaled_recruit_095q)

interpol_scaRecruit_vs_meanLB_DHW_df = scaledrecruit_vs_meanLB

# Save combined_df to a CSV file
write_xlsx(interpol_scaRecruit_vs_meanLB_DHW_df, paste0("interpol_scaRecruit_vs_meanLB_DHW_df", ".xlsx"))

ggplot(interpol_scaRecruit_vs_meanLB_DHW_df,aes(x = average_LB, y = scaled_recruit_05q))
ggplot(interpol_scaRecruit_vs_meanLB_DHW_df, aes(x = average_LB, y = as.numeric(scaled_recruit_05q))) +
  geom_point() +
  theme_light()
##### future projections, simulations (dfs) for LB, average_LB, and surv_p based on Posteriors ----
## preparations ----
theme_publication <- function(base_size = 11, base_family = "") {
  theme_light(base_size = base_size, base_family = base_family) +
    theme(text = element_text(size = base_size),
          #axis.text.y = element_text(size = rel(1), color = "#0b3c5d"),
          axis.text.y.right = element_text(size = rel(1), color = "#CC3311"),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.y.right = element_blank(),
          axis.text.x = element_text(size = rel(1.1), angle = 45, hjust = 1),
          panel.spacing = unit(0.2, "lines"), # existing theme elements here...
          plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm") ) } # adjust the space between legend entries
temperature_df <- read_excel(paste0(main_address, "/Temp_data/KOB2022_summer/NEW/df_GHL_long.xlsx")) # def_all_long was made via another script: Temp_data/KOB2022_summer/Temp_all.R
temperature_df <- temperature_df %>%
  arrange(Treatment, date_h)
head(temperature_df)
T_ref = 28
LTT = 15 # arbitrary chosen as very low, with no effect on the outcome

### extract parameter set with the highest weight
posteriors <- particles_weights_MADs_df_DHW[which.max(particles_weights_MADs_df_DHW$weight), ]

### !!! Choose projected temperature data of 5 warmest summers (by 2099) and the respective address for saving plots
main_address = paste0("~/Library/Mobile Documents/com~apple~CloudDocs/Jahangir Air/Postdoc_Position:Propos:RelatedStuff/PROPOSAL_postdoc/Selection_docs/Niklas_Mytilus_experiment/")

#### (1) without daily fluctuations ----
temperature_df <- read_csv("~/Library/Mobile Documents/com~apple~CloudDocs/Jahangir Air/Postdoc_Position:Propos:RelatedStuff/PROPOSAL_postdoc/Selection_docs/Niklas_Mytilus_experiment/Temp_data/Future_Temp_projections/kiel_multiple_grids/kiel_rcp85_A006/R_output/five_warmestsummer_RCP85_projection_1005grid_hourly.csv")
temperature_df$sample_name <- lubridate::year(temperature_df$Datetime)
colnames(temperature_df)[1] = "date_h"
dir.create(file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "projections", "noFluc"), showWarnings = F, recursive = T)
address_outputs <- file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "projections", "noFluc")
setwd(paste0(address_outputs))
sample_names = unique(temperature_df$sample_name)

## Function to simulate a population of LBRs (MC)----
simulate_T_LBR_pop <- function(temperature_sample_df, meanlog10LBR_refT, sdlog10LBR_refT, k, T_ref) {
  N_ <- 100
  T_LBR_pop_df <- lapply(1:N_, function(ind) {
    repeat {
      LBR_T_ref <- rlnorm(1, meanlog = log(10^meanlog10LBR_refT), sdlog = log(10^sdlog10LBR_refT))
      if (!is.na(LBR_T_ref)) {
        break
      }
    }
    Temp <- temperature_sample_df$Temperature
    # Use only if you prefer to set a Lethal Temperature Threshold below which LBR = 0.
    LBR_T <- ifelse(Temp <= LTT, 0, LBR_T_ref * (10^(k * (Temp - T_ref))))
    counter <- seq_along(Temp)
    LB <- cumsum(LBR_T /100)
    LB[LB > 1] <- 1
    T_LBR_ind_df <- data.frame(t_h = counter, Temp = Temp, ind = ind, LBR = LBR_T, LB = LB)
    return(T_LBR_ind_df)
  })
  T_LBR_pop_df <- do.call(rbind, T_LBR_pop_df)
  T_LBR_pop_df <- as.data.frame(T_LBR_pop_df)
  colnames(T_LBR_pop_df) <- c('t_h', 'Temp', 'ind', 'LBR', 'LB')
  T_LBR_pop_df_sub <- subset(T_LBR_pop_df, select = -c(LBR))
  T_LBR_pop_df_wide <- T_LBR_pop_df_sub %>% group_by(t_h, Temp) %>% tidyr::pivot_wider(names_from = ind, values_from = LB)
  T_LBR_surv_pop_df_wide <- as.data.frame(T_LBR_pop_df_wide) ###
  T_LBR_pop_df_wide$cum_mortality_sim <- rowSums(T_LBR_pop_df_wide[-(1:2)] == 1)
  T_LBR_pop_df_wide$surv_p_sim <- 1 - (T_LBR_pop_df_wide$cum_mortality_sim / N_)
  T_LBR_pop_df_wide$average_LB <- rowMeans(T_LBR_pop_df_wide[, 3:(2 + N_)])
  T_LBR_pop_df_wide$date_h <- temperature_sample_df[['date_h']]
  T_LBR_pop_df_wide = as.data.frame(T_LBR_pop_df_wide)
  return(T_LBR_pop_df_wide)
}

## ggplot combined for LB and surv_p_sim based on the posterior_means----
system.time({
  setwd(paste0(address_outputs))
  plot_list <- list()
  counter <- 0
  for (s_n in sample_names) {
    counter <- counter + 1
    temperature_sample_df = subset(temperature_df, sample_name == s_n)
    n <- nrow(temperature_sample_df)
    temperature_sample_df$t_h <- 1:n
    weekly_breaks <- seq(min(temperature_sample_df$date_h), max(temperature_sample_df$date_h), by = "month")
    all_sim_surv_curves <- list()
    
    for (row_n in 1:nrow(posteriors)) {
      meanlog10LBR_refT <- posteriors[row_n, 1]
      sdlog10LBR_refT <- posteriors[row_n, 2]
      k <- posteriors[row_n, 3]
      T_LBR_pop_df_wide <- simulate_T_LBR_pop(temperature_sample_df, meanlog10LBR_refT, sdlog10LBR_refT, k, T_ref)
      T_LBR_pop_df_wide = as.data.frame(T_LBR_pop_df_wide)
      
      # Add parameter values and sample name as new columns
      T_LBR_pop_df_wide$meanlog10LBR_refT <- meanlog10LBR_refT
      T_LBR_pop_df_wide$sdlog10LBR_refT <- sdlog10LBR_refT
      T_LBR_pop_df_wide$k <- k
      T_LBR_pop_df_wide$sample_name <- s_n
      
      all_sim_surv_curves[[paste("meanlog10LBR_refT", meanlog10LBR_refT, "sdlog10LBR_refT", sdlog10LBR_refT, "k", k, sep = "_")]] <- T_LBR_pop_df_wide
    }
    
    # Combine all data frames in the all_sim_surv_curves list into one data frame
    combined_df <- bind_rows(all_sim_surv_curves)
    
    # Save combined_df to a CSV file
    write.csv(combined_df, file = paste0("all_sim_LB_surv_df_", s_n, ".csv"), row.names = FALSE)
    
    
    p1 <- ggplot()
    N_ <- 100
    for (curve_name in names(all_sim_surv_curves)) {
      p1 <- p1 + geom_line(data = all_sim_surv_curves[[curve_name]], aes(x = date_h, y = surv_p_sim), color = "#33BBEE", size = 0.7)
      
      for (LB_col in 1:N_) {
        p1 <- p1 + geom_line(data = all_sim_surv_curves[[curve_name]], 
                             aes(x = date_h, y = !!sym(as.character(LB_col))), 
                             color = "#009E73", 
                             alpha = 0.3, 
                             size = 0.1)
      }
    }
    
    p1 <- p1 +
      geom_line(data = temperature_sample_df, aes(x = date_h, y = (Temperature - 16) / 10), size = 0.4, color = "#CC3311") +
      
      scale_x_datetime(labels = if (counter == length(sample_names)) scales::date_format("%b %d") else NULL,
                       breaks = scales::breaks_pretty(n = 4) ) +
      
      scale_y_continuous(name = if (counter == 3) "Surv. prob. / Leth. buil." else "", limits = c(0, 1),
                         sec.axis = sec_axis(~ . * 10 + 16, name = if (counter == 3) "Temperature [°C]" else "",
                                             labels = scales::number_format(accuracy = 1),
                                             breaks = scales::pretty_breaks(n = 5)))
    
    plot_list[[counter]] <- p1
  }

})

for (i in seq_along(plot_list)) {plot_list[[i]] <- plot_list[[i]] + theme_publication() }

combined_plot_LB_Surv_DHW <- plot_grid(plotlist = plot_list, ncol = 1,
                                       rel_heights = c(1, 1, 1, 1, 1.2725))
# Add global y-axis titles using ggdraw and draw_label
final_plot <- ggdraw() +
  draw_plot(combined_plot_LB_Surv_DHW) + 
  draw_label("Survival probability", x = -0.025, y = 0.39, angle = 90, hjust = 0.5, vjust = 0.5, size = 12, color = "#33BBEE") +
  draw_label("Lethality buildup (scaled)", x = -0.025, y = 0.67, angle = 90, hjust = 0.5, vjust = 0.5, size = 12, color = "#009E73") +
  draw_label("Temperature [°C]", x = 1.025, y = 0.5, angle = 90, hjust = 0.5, vjust = 0.5, size = 12, color = "#CC3311") +
  theme(plot.margin = margin(t = 2, r = 12, b = 2, l = 12, unit = "mm") ) 
ggsave("Projections_NoFluc_ggplot.pdf", final_plot, width = plot_width_FullPage/2.3, height = 19, units = "cm")


## Future projection of recruitment ----
# Set your working directory
setwd(file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "projections", "noFluc")) 

# Get the list of all_sim_LB_surv_df (csv) files in the directory
file_list <- list.files(pattern = "*.csv")
post_HW_proj_meanLBs_df <- data.frame()
# Loop over all csv files
for(file in file_list){
  # Read the csv file into a data frame
  df <- read_csv(file)
  # Subset the columns you need
  df <- df %>% select(average_LB, date_h)
  # Add a new column "sample_name" using the last two characters of the file name (without extension)
  df$sample_name <- paste0("Fluc_", substr(file, nchar(file)-8, nchar(file)-4))
  # Keep only the last row
  df <- tail(df, 1)
  # Append the row to the final data frame
  post_HW_proj_meanLBs_df <- rbind(post_HW_proj_meanLBs_df, df)
}
post_HW_proj_meanLBs_df$average_LB = formatC(round(post_HW_proj_meanLBs_df$average_LB, 3), format = "f", digits = 3)
post_HW_proj_meanLBs_df$average_LB <- as.numeric(post_HW_proj_meanLBs_df$average_LB)
head(post_HW_proj_meanLBs_df)

interpol_scaRecruit_vs_meanLB_DHW_df <- read_excel(paste0(file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "DHW", "LB_surv_post_mean"), "/interpol_scaRecruit_vs_meanLB_DHW_df.xlsx"))
head(interpol_scaRecruit_vs_meanLB_DHW_df)

interpol_scaRecruit_vs_meanLB_DHW_df$average_LB <- round(interpol_scaRecruit_vs_meanLB_DHW_df$average_LB, 3)
post_HW_proj_meanLBs_df$average_LB <- round(post_HW_proj_meanLBs_df$average_LB, 3)

# Perform the inner join
proj_recruit_df <- post_HW_proj_meanLBs_df %>%
  left_join(interpol_scaRecruit_vs_meanLB_DHW_df, by = "average_LB")

# Save combined_df to a CSV file
write_xlsx(proj_recruit_df, paste0("noFluc_proj_recruit_df", ".xlsx"))




#### (2) with daily fluctuations----
temperature_df <- read_csv("~/Library/Mobile Documents/com~apple~CloudDocs/Jahangir Air/Postdoc_Position:Propos:RelatedStuff/PROPOSAL_postdoc/Selection_docs/Niklas_Mytilus_experiment/Temp_data/Future_Temp_projections/kiel_multiple_grids/kiel_rcp85_A006/R_output/five_warmestsummer_RCP85_projection_1005grid_hourly_fluc.csv")
temperature_df$sample_name <- lubridate::year(temperature_df$Datetime)
colnames(temperature_df)[1] = "date_h"
dir.create(file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "projections", "Fluc"), showWarnings = F, recursive = T)
address_outputs <- file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "projections", "Fluc")
setwd(paste0(address_outputs))
sample_names = unique(temperature_df$sample_name)

## Function to simulate a population of LBRs (MC)----
simulate_T_LBR_pop <- function(temperature_sample_df, meanlog10LBR_refT, sdlog10LBR_refT, k, T_ref) {
  N_ <- 100
  T_LBR_pop_df <- lapply(1:N_, function(ind) {
    repeat {
      LBR_T_ref <- rlnorm(1, meanlog = log(10^meanlog10LBR_refT), sdlog = log(10^sdlog10LBR_refT))
      if (!is.na(LBR_T_ref)) {
        break
      }
    }
    Temp <- temperature_sample_df$Temperature
    # Use only if you prefer to set a Lethal Temperature Threshold below which LBR = 0.
    LBR_T <- ifelse(Temp <= LTT, 0, LBR_T_ref * (10^(k * (Temp - T_ref))))
    counter <- seq_along(Temp)
    LB <- cumsum(LBR_T /100)
    LB[LB > 1] <- 1
    T_LBR_ind_df <- data.frame(t_h = counter, Temp = Temp, ind = ind, LBR = LBR_T, LB = LB)
    return(T_LBR_ind_df)
  })
  T_LBR_pop_df <- do.call(rbind, T_LBR_pop_df)
  T_LBR_pop_df <- as.data.frame(T_LBR_pop_df)
  colnames(T_LBR_pop_df) <- c('t_h', 'Temp', 'ind', 'LBR', 'LB')
  T_LBR_pop_df_sub <- subset(T_LBR_pop_df, select = -c(LBR))
  T_LBR_pop_df_wide <- T_LBR_pop_df_sub %>% group_by(t_h, Temp) %>% tidyr::pivot_wider(names_from = ind, values_from = LB)
  T_LBR_surv_pop_df_wide <- as.data.frame(T_LBR_pop_df_wide) ###
  T_LBR_pop_df_wide$cum_mortality_sim <- rowSums(T_LBR_pop_df_wide[-(1:2)] == 1)
  T_LBR_pop_df_wide$surv_p_sim <- 1 - (T_LBR_pop_df_wide$cum_mortality_sim / N_)
  T_LBR_pop_df_wide$average_LB <- rowMeans(T_LBR_pop_df_wide[, 3:(2 + N_)])
  T_LBR_pop_df_wide$date_h <- temperature_sample_df[['date_h']]
  T_LBR_pop_df_wide = as.data.frame(T_LBR_pop_df_wide)
  return(T_LBR_pop_df_wide)
}

## ggplot combined for LB and surv_p_sim based on the posterior_means----
system.time({
  setwd(paste0(address_outputs))
  plot_list <- list()
  counter <- 0
  for (s_n in sample_names) {
    counter <- counter + 1
    temperature_sample_df = subset(temperature_df, sample_name == s_n)
    n <- nrow(temperature_sample_df)
    temperature_sample_df$t_h <- 1:n
    weekly_breaks <- seq(min(temperature_sample_df$date_h), max(temperature_sample_df$date_h), by = "month")
    all_sim_surv_curves <- list()
    
    for (row_n in 1:nrow(posteriors)) {
      meanlog10LBR_refT <- posteriors[row_n, 1]
      sdlog10LBR_refT <- posteriors[row_n, 2]
      k <- posteriors[row_n, 3]
      T_LBR_pop_df_wide <- simulate_T_LBR_pop(temperature_sample_df, meanlog10LBR_refT, sdlog10LBR_refT, k, T_ref)
      T_LBR_pop_df_wide = as.data.frame(T_LBR_pop_df_wide)
      
      # Add parameter values and sample name as new columns
      T_LBR_pop_df_wide$meanlog10LBR_refT <- meanlog10LBR_refT
      T_LBR_pop_df_wide$sdlog10LBR_refT <- sdlog10LBR_refT
      T_LBR_pop_df_wide$k <- k
      T_LBR_pop_df_wide$sample_name <- s_n
      
      all_sim_surv_curves[[paste("meanlog10LBR_refT", meanlog10LBR_refT, "sdlog10LBR_refT", sdlog10LBR_refT, "k", k, sep = "_")]] <- T_LBR_pop_df_wide
    }
    
    # Combine all data frames in the all_sim_surv_curves list into one data frame
    combined_df <- bind_rows(all_sim_surv_curves)
    
    # Save combined_df to a CSV file
    write.csv(combined_df, file = paste0("all_sim_LB_surv_df_", s_n, ".csv"), row.names = FALSE)
    
    
    p1 <- ggplot()
    N_ <- 100
    for (curve_name in names(all_sim_surv_curves)) {
      p1 <- p1 + geom_line(data = all_sim_surv_curves[[curve_name]], aes(x = date_h, y = surv_p_sim), color = "#33BBEE", size = 0.7)
      
      for (LB_col in 1:N_) {
        p1 <- p1 + geom_line(data = all_sim_surv_curves[[curve_name]], 
                             aes(x = date_h, y = !!sym(as.character(LB_col))), 
                             color = "#009E73", 
                             alpha = 0.3, 
                             size = 0.1)
      }
    }
    
    p1 <- p1 +
      geom_line(data = temperature_sample_df, aes(x = date_h, y = (Temperature - 16) / 10), size = 0.1, color = "#CC3311") +
      
      scale_x_datetime(labels = if (counter == length(sample_names)) scales::date_format("%b %d") else NULL,
                       breaks = scales::breaks_pretty(n = 4) ) +
      
      scale_y_continuous(name = if (counter == 3) "Surv. prob. / Leth. buil." else "", limits = c(0, 1),
                         sec.axis = sec_axis(~ . * 10 + 16, name = if (counter == 3) "Temperature [°C]" else "",
                                             labels = scales::number_format(accuracy = 1),
                                             breaks = scales::pretty_breaks(n = 5)))
    
    plot_list[[counter]] <- p1
  }
  
})

for (i in seq_along(plot_list)) {plot_list[[i]] <- plot_list[[i]] + theme_publication() }

combined_plot_LB_Surv_DHW <- plot_grid(plotlist = plot_list, ncol = 1,
                                       rel_heights = c(1, 1, 1, 1, 1.2725))
# Add global y-axis titles using ggdraw and draw_label
final_plot <- ggdraw() +
  draw_plot(combined_plot_LB_Surv_DHW) + 
  draw_label("Survival probability", x = -0.025, y = 0.39, angle = 90, hjust = 0.5, vjust = 0.5, size = 12, color = "#33BBEE") +
  draw_label("Lethality buildup (scaled)", x = -0.025, y = 0.67, angle = 90, hjust = 0.5, vjust = 0.5, size = 12, color = "#009E73") +
  draw_label("Temperature [°C]", x = 1.025, y = 0.5, angle = 90, hjust = 0.5, vjust = 0.5, size = 12, color = "#CC3311") +
  theme(plot.margin = margin(t = 2, r = 12, b = 2, l = 12, unit = "mm") ) 
ggsave("Projections_Fluc_ggplot.pdf", final_plot, width = plot_width_FullPage/2.3, height = 19, units = "cm")


## Future projection of recruitment----
# Set your working directory
setwd(file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "projections", "Fluc")) 

# Get the list of all_sim_LB_surv_df (csv) files in the directory
file_list <- list.files(pattern = "*.csv")
post_HW_proj_meanLBs_df <- data.frame()
# Loop over all csv files
for(file in file_list){
  # Read the csv file into a data frame
  df <- read_csv(file)
  # Subset the columns you need
  df <- df %>% select(average_LB, date_h)
  # Add a new column "sample_name" using the last two characters of the file name (without extension)
  df$sample_name <- paste0("Fluc_", substr(file, nchar(file)-8, nchar(file)-4))
  # Keep only the last row
  df <- tail(df, 1)
  # Append the row to the final data frame
  post_HW_proj_meanLBs_df <- rbind(post_HW_proj_meanLBs_df, df)
}
post_HW_proj_meanLBs_df$average_LB = formatC(round(post_HW_proj_meanLBs_df$average_LB, 3), format = "f", digits = 3)
post_HW_proj_meanLBs_df$average_LB <- as.numeric(post_HW_proj_meanLBs_df$average_LB)
head(post_HW_proj_meanLBs_df)

interpol_scaRecruit_vs_meanLB_DHW_df <- read_excel(paste0(file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "DHW", "LB_surv_post_mean"), "/interpol_scaRecruit_vs_meanLB_DHW_df.xlsx"))
head(interpol_scaRecruit_vs_meanLB_DHW_df)

interpol_scaRecruit_vs_meanLB_DHW_df$average_LB <- round(interpol_scaRecruit_vs_meanLB_DHW_df$average_LB, 3)
post_HW_proj_meanLBs_df$average_LB <- round(post_HW_proj_meanLBs_df$average_LB, 3)

# Perform the inner join
proj_recruit_df <- post_HW_proj_meanLBs_df %>%
  left_join(interpol_scaRecruit_vs_meanLB_DHW_df, by = "average_LB")

# Save combined_df to a CSV file
write_xlsx(proj_recruit_df, paste0("Fluc_proj_recruit_df", ".xlsx"))














































######################## end
## ggplot combined for surv_p_sim (Note: choose proper name for the plot)----
system.time({
  setwd(paste0(address_outputs))
  
  plot_list <- list()
  counter <- 0
  
  for (s_n in sample_names) {
    counter <- counter + 1
    
    temperature_sample_df = subset(temperature_df, sample_name == s_n)
    
    n <- nrow(temperature_sample_df)
    temperature_sample_df$t_h <- 1:n
    weekly_breaks <- seq(min(temperature_sample_df$date_h), max(temperature_sample_df$date_h), by = "month")
    
    all_sim_surv_curves <- list()
    
    for (row_n in 1:nrow(posteriors)) {
      meanlog10LBR_refT <- posteriors[row_n, 1]
      sdlog10LBR_refT <- posteriors[row_n, 2]
      k <- posteriors[row_n, 3]
      
      T_LBR_pop_df_wide <- simulate_T_LBR_pop(temperature_sample_df, meanlog10LBR_refT, sdlog10LBR_refT, k, T_ref)
      T_LBR_pop_df_wide = as.data.frame(T_LBR_pop_df_wide)
      
      all_sim_surv_curves[[paste("meanlog10LBR_refT", meanlog10LBR_refT, "sdlog10LBR_refT", sdlog10LBR_refT, "k", k, sep = "_")]] <- T_LBR_pop_df_wide
    }
    
    p1 <- ggplot()
    
    for (curve_name in names(all_sim_surv_curves)) {
      p1 <- p1 + geom_line(data = all_sim_surv_curves[[curve_name]], aes(x = date_h, y = surv_p_sim), color = "#33BBEE", size = 0.1)
    }
    
    p1 <- p1 +
      geom_line(data = temperature_sample_df, aes(x = date_h, y = (Temperature - 16) / 10), size = 0.2, color = "#CC3311") +
      
      scale_x_datetime(limits = range(weekly_breaks),
                       labels = if (counter == length(sample_names)) scales::date_format("%b %d") else NULL,
                       breaks = weekly_breaks) +
      scale_y_continuous(name = "Survival probability", limits = c(0, 1),
                         sec.axis = sec_axis(~ . * 10 + 16, name = "Temperature [°C]",
                                             labels = scales::number_format(accuracy = 1),
                                             breaks = scales::pretty_breaks(n = 5))) +
      theme_light() +
      theme(plot.margin=unit(c(0.2,0.6,0.2,0.6), "cm"),
            plot.title = element_text(size = 14, face = "bold"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 12, face = "bold", vjust = 2),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            axis.title.y.right = element_text(vjust = 2), # Increase the left margin of the secondary y-axis tick labels
            legend.position = "top")
    
    plot_list[[counter]] <- p1
  }
  
  # Create a vector of relative heights for each subplot
  subplot_heights <- rep(1, length(sample_names))
  subplot_heights[length(sample_names)] <- 1.21 # Adjust the height of the last plot
  
  # Save combined plot to PDF
  pdf(file = "posterior_means_combined_ggplot.pdf", width = 3.5, height = 10)
  grid.arrange(grobs = plot_list, ncol = 1, heights = subplot_heights) # Add the heights parameter
  dev.off()
})

##### Convert ALL PDFs to PNGs----
# This script should find all PDF files in the main directory and its subdirectories and save the extracted images as PNG files in the same directory as the corresponding PDF.
if (!requireNamespace("pdftools", quietly = TRUE)) {
  install.packages("pdftools")
}

library(pdftools)

install.packages("extrafont")
library(extrafont)
font_import()


# Set your main directory path
main_address <- "~/Library/Mobile Documents/com~apple~CloudDocs/Jahangir Air/Postdoc_Position:Propos:RelatedStuff/PROPOSAL_postdoc/Selection_docs/Niklas_Mytilus_experiment"
main_directory <- file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT")

# Find all PDF files in the main directory and its subdirectories
pdf_files <- list.files(path = main_directory, pattern = "\\.pdf$", ignore.case = TRUE, full.names = TRUE, recursive = TRUE)

# Function to extract and save images from a PDF file
process_pdf <- function(pdf_file) {
  # Get the directory of the PDF file
  pdf_directory <- dirname(pdf_file)
  
  # Convert the PDF pages to images
  pages <- pdftools::pdf_convert(pdf_file, dpi = 300)
  
  for (i in seq_along(pages)) {
    output_path <- file.path(pdf_directory, paste0(tools::file_path_sans_ext(basename(pdf_file)), "_page_", i, ".png"))
    file.copy(pages[i], output_path)
  }
}

# Process all PDF files
for (pdf_file in pdf_files) {
  process_pdf(pdf_file)
}

### second approach
install.packages("magick")
library(magick)

# Set your main directory path
main_address <- "~/Library/Mobile Documents/com~apple~CloudDocs/Jahangir Air/Postdoc_Position:Propos:RelatedStuff/PROPOSAL_postdoc/Selection_docs/Niklas_Mytilus_experiment"
main_directory <- file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT")

# Find all PDF files in the main directory and its subdirectories
pdf_files <- list.files(path = main_directory, pattern = "\\.pdf$", ignore.case = TRUE, full.names = TRUE, recursive = TRUE)

# Function to extract and save images from a PDF file
process_pdf <- function(pdf_file) {
  # Get the directory of the PDF file
  pdf_directory <- dirname(pdf_file)
  
  # Convert the PDF pages to images
  pages <- magick::image_read_pdf(pdf_file, density = 300)
  
  for (i in seq_along(pages)) {
    output_path <- file.path(pdf_directory, paste0(tools::file_path_sans_ext(basename(pdf_file)), "_page_", i, ".png"))
    magick::image_write(pages[i], path = output_path, format = "png")
  }
}

# Process all PDF files
for (pdf_file in pdf_files) {
  process_pdf(pdf_file)
}


