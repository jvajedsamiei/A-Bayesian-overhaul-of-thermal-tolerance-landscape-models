' 
The script belongs to the Supplementary Text 3: Bayesian regression analysis of data from constant heatwave experiment
Written by J.Vajedsamiei
'

################# Packages and Libraries----
install.packages("readxl")
install.packages("data.table")
install.packages('ggplot2')
install.packages("plotly")
install.packages("htmlTable")
install.packages('lme4')
install.packages("pscl")
install.packages("MASS")
install.packages("tidyverse")
install.packages("mgcViz")
install.packages('brms')
install.packages('merTools')
install.packages("tidybayes")
install.packages("modelr")
install.packages('survival')
install.packages('survminer')
install.packages('ggsurvfit')
install.packages('zoo')
install.packages('cowplot')
install.packages('writexl')
install.packages("bayesplot")
install.packages("ggplot2")
install.packages("gridExtra")

library(ggridges)
library(writexl)
library(cowplot)
library(zoo)
#library(ggsurvfit)
library(survival)
#library(survminer)
library(modelr)
library(tidybayes)
library(merTools)
library(brms)
library(tidyverse)
library(pscl)
library(readxl)  
library(writexl)
library(data.table)
library(ggplot2)
library(plotly)
library(htmlTable)
library(MASS)
library(ggpubr)
library(texreg)
library(dplyr)
library(scales)
library(gridExtra)
library(bayesplot)
library(ggplot2)
library(brms)  # Assuming you have already installed brms for your model
library(patchwork)
library(ggplot2)
library(patchwork)
library(dplyr)
library(modelr)  # Assuming data_grid and other functions come from this package or similar
library(brms)    # Assuming add_epred_draws and add_predicted_draws come from brms or similar


##### Bayesian linear regression analysis of censored lethality buildup rates over constant temperatures ----
### Address & directory definition ----
main_address = paste0("~/Library/Mobile Documents/com~apple~CloudDocs/Jahangir Air/Postdoc_Position:Propos:RelatedStuff/PROPOSAL_postdoc/Selection_docs/Niklas_Mytilus_experiment")
dir.create(file.path(main_address, "Outputs", "Bayes_regres"), showWarnings = FALSE) #stops warnings if folder already exists

# Address to save regression results (posterior distributions)
dir.create(file.path(main_address, "Outputs", "Bayes_regres", "Bayes_regres_results"), showWarnings = FALSE) #stops warnings if folder already exists
outputs_Bayes_regres_address = file.path(main_address, "Outputs", "Bayes_regres", "Bayes_regres_results")

# Address to save posterior predictions (simulations)
dir.create(file.path(main_address, "Outputs", "Bayes_regres", "simulation"), showWarnings = FALSE) #stops warnings if folder already exists
outputs_simulation_address = file.path(main_address, "Outputs", "Bayes_regres", "simulation")

# Import CHW survival time df
obs_surv_time_CHW_df = read_excel(file.path(main_address, "HeatSurvival_data/Surv_time_df_formatted_for_brm_PART1.xlsx"))

### brm modeling ----
mod <- brm(log10_LBR | cens(cen1) ~ Temp_C, data = obs_surv_time_CHW_df, chains = 3, iter = 20000, warmup = 2000, thin = 4)
summary(mod)
post_summ_df = posterior_summary(mod)

# Check posterior distributions and chains
par(mfrow=c(1,1), mar = c(2,4,1,3.8), cex=1, mgp = c(1, 1, 0.2)) 
plot(mod)

### Plotting theme ----
theme_publication <- function(base_size = 11, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(text = element_text(size = base_size),
          axis.title = element_text(size = rel(1.1)),
          axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)), # adjust top margin for x axis title
          axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)), # adjust right margin for y axis title
          axis.text = element_text(size = rel(1)),
          strip.text = element_text(size = rel(1.1)),
          legend.text = element_text(size = rel(1)),
          legend.title = element_text(size = rel(1.1)),
          strip.background = element_blank(),
          panel.spacing = unit(0.3, "lines"), # existing theme elements here...
          legend.spacing = unit(0.8, 'lines'),
          legend.box.margin = margin(t = 0, r = 0, b = -5, l = 0, unit = "pt"),
          plot.margin = margin(t = 1, r = 5, b = 1, l = 5, unit = "mm") ) } # adjust the space between legend entries
plot_width_FullPage = 17.3

### Plot Posterior Predictions ----
'The first visualization displays a series of predictive lines representing 100 draws from the posterior distribution 
of a Bayesian regression model, overlaid on a scatter plot of observed data. The second visualization shows 
predictive intervals around the median of the model predictions, again overlaid on the same scatter plot of observed data. 
Both visualizations plot temperature against the logarithm (base 10) of a variable'

# Create the first plot and store it into a variable
plot1 <- obs_surv_time_CHW_df %>%
  data_grid(Temp_C = seq_range(Temp_C, n = 50)) %>%
  add_epred_draws(mod, ndraws = 100) %>%
  ggplot(aes(x = Temp_C, y = log10_LBR)) +
  geom_line(aes(y = .epred, group = paste(.draw)), alpha = .1) +
  geom_point(data = obs_surv_time_CHW_df) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = 'Temperature [°C]', y = expression(log[10](dot(L[r])))) + theme_publication()
# Create the second plot and store it into a variable
plot2 <- obs_surv_time_CHW_df %>%
  data_grid(Temp_C = seq_range(Temp_C, n = 50)) %>%
  add_predicted_draws(mod) %>%
  ggplot(aes(x = Temp_C, y = log10_LBR)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.99, .95, .68, .5)) +
  geom_point(data = obs_surv_time_CHW_df) +
  scale_fill_brewer() +
  labs(x = 'Temperature [°C]', y = expression(log[10](dot(L[r])))) + theme_publication()

# Combine the plots using patchwork
combined_PostPred_plot <- plot1 + plot2 + theme_publication()
# Define the output file name and path
output_file_name <- paste0(outputs_Bayes_regres_address, "/brm_combined_PostPred_plots", ".pdf")
# Save the combined plot to a PDF file
ggsave(output_file_name, combined_PostPred_plot, width = plot_width_FullPage, height = 6, units = "cm")


### Posterior Parameter distributions----
# Draw posterior sets of TTL parameters: mean(log10(LBR_refT)), sd(log10LBR), and k [also for log10(LBR_refT)]
set.seed(12345)
T_ref = 28
postpred_manual <- mod |> 
  spread_draws(b_Intercept, b_Temp_C, sigma, ndraws = 10000, seed = NULL) |> 
  mutate(mu = b_Intercept + (b_Temp_C * T_ref),  # This is posterior_linpred()
         pred = rnorm(n(), mean = mu, sd = sigma))  # This is posterior_predict()
post_df = as.data.frame(postpred_manual |> select(.draw:pred))
colnames(post_df) = c('draw', 'Intercept', 'k', 'sdlog10LBR_refT', 'meanlog10LBR_refT', 'log10LBR_refT')
# Long dataframe
PostParam_BR_CHW_df <- post_df[,c(3:5)] %>%                    
  pivot_longer(colnames(post_df[,c(3:5)])) %>% 
  as.data.frame()

# Filter data for each variable
k_data <- PostParam_BR_CHW_df %>% filter(name == "k")
sd_data <- PostParam_BR_CHW_df %>% filter(name == "sdlog10LBR_refT")
mean_data <- PostParam_BR_CHW_df %>% filter(name == "meanlog10LBR_refT")

# Create the plots
p1 <- ggplot(k_data, aes(x = value)) +
  geom_histogram(bins = 30, alpha = 0.7, fill = "gray30") +
  scale_x_continuous(breaks = pretty(k_data$value, n = 3)) +
  labs(title = "k") + theme_publication()
p2 <- ggplot(sd_data, aes(x = value)) +
  geom_histogram(bins = 30, alpha = 0.7, fill = "gray30") +
  scale_x_continuous(breaks = pretty(sd_data$value, n = 3)) +
  labs(title = expression(sd(log[10](dot(L)[r])))) + theme_publication()
p3 <- ggplot(mean_data, aes(x = value)) +
  geom_histogram(bins = 30, alpha = 0.7, fill = "gray30") +
  scale_x_continuous(breaks = pretty(mean_data$value, n = 3)) +
  labs(title = expression(mean(log[10](dot(L)[r])))) + theme_publication()

# Combine the plots
combined_PostParamDist_plot <- p1 + p2 + p3 + plot_layout(ncol = 3) + theme_publication()
output_file_name <- paste0(outputs_Bayes_regres_address, "/PostParamDist_combined_plots", ".pdf")
ggsave(output_file_name, combined_PostParamDist_plot, width = plot_width_FullPage, height = 6, units = "cm")


### Compare ABC_SMC with Bayes_Regres (NOTE: combined_df_long is made by ABS_SMC_meanMAD_v4_future_proj.r)----
PostParam_ABCSMC_CHW_DHW_df <- read_excel(file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "PostParam_ABCSMC_CHW_DHW_df.xlsx")) # for us 
head(PostParam_BR_CHW_df)

PostParam_ABCSMC_CHW_DHW_df_filtered <- PostParam_ABCSMC_CHW_DHW_df[,4:6] %>%
  mutate(experiment = case_when(
    experiment == "DHW" ~ "DHW (ABC-SMC)",
    experiment == "CHW" ~ "CHW (ABC-SMC)",
    TRUE ~ experiment
  ))

PostParam_BR_CHW_df$experiment = "CHW (BR)"
colnames(PostParam_BR_CHW_df)[1:2] = c("Parameter", "para_value")        

# Assuming PostParam_BR_CHW_df and combined_df_filtered are already defined
Post_ABCSMC_BR <- rbind(PostParam_BR_CHW_df, PostParam_ABCSMC_CHW_DHW_df_filtered)

# Specify the desired order of the experiment levels
desired_order <- c("CHW (BR)", "CHW (ABC-SMC)", "DHW (ABC-SMC)")  # Adjust this to your specific desired order

# Reorder the levels of the 'experiment' factor
Post_ABCSMC_BR$experiment <- factor(Post_ABCSMC_BR$experiment, levels = desired_order)

Post_ABCSMC_BR_plot = ggplot(Post_ABCSMC_BR, aes(x = para_value, y = experiment, fill = factor(after_stat(quantile)))) +
  stat_density_ridges(geom = "density_ridges_gradient", bandwidth = 0.02, calc_ecdf = TRUE, quantiles = c(0.025, 0.975)) +
  scale_fill_manual(name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
                    labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")) +
  labs(x = "Parameter value", y = "Experiment") +
  facet_grid(. ~ Parameter, scales = "free_x") +
  theme_publication() + theme_light() +
  theme(plot.margin = unit(c(0.2, 0.6, 0.2, 0.6), "cm"),
        plot.title = element_blank(), strip.text.x = element_blank(),
        legend.position = "top", legend.title = element_text(size = 11))
output_file_name <- paste0(outputs_Bayes_regres_address, "/Post_ABCSMC_BR_plot", ".pdf")
ggsave(output_file_name, Post_ABCSMC_BR_plot, width = plot_width_FullPage, height = 6, units = "cm")


# parameters of the distributions (was used before, but now the distributions are directly used)
# Calculate summary statistics for each parameter within each experiment
summary_stats <- Post_ABCSMC_BR %>%
  group_by(Parameter, experiment) %>%
  summarise(
    Median = median(para_value),
    Quantile_0.05 = quantile(para_value, probs = 0.05),
    Quantile_0.95 = quantile(para_value, probs = 0.95)
  ) %>%
  ungroup()

# View the summary statistics
print(summary_stats)

# Install and load the officer package if you haven't already
if (!requireNamespace("officer", quietly = TRUE)) {
  install.packages("officer")
}
library(officer)

# Create a new Word document
doc <- read_docx()

# Add the summary statistics table to the Word document
doc <- doc %>%
  body_add_table(value = summary_stats, style = "table_template")

# Save the Word document
print(doc, target = paste0(outputs_Bayes_regres_address, "Summary_Statistics.docx"))



### Posterior prediction plot of simulations versus observations ----
'Simulating population function (simulate_T_LBR_pop): The function simulates a population of organisms under given temperature conditions. The function takes a dataframe with the temperature data, parameters for a reference temperature, mean and standard deviation for LBR at the reference temperature, and a constant k. It computes LBR at different temperatures and accumulates the sum of LBR. It then generates a dataframe with columns for time, temperature, individual number, LBR and accumulated LBR. This dataframe is then processed to compute the simulated survival probability at different times.
Observed survival analysis function (observed_survival_analysis): This function takes the temperature data, observed survival data, and a sample name as input, and generates Kaplan-Meier survival probabilities at different times. The function also interpolates missing values in the survival probabilities and returns the processed dataframe.
Plotting the simulated vs observed data: The script then enters a loop where it generates simulated survival curves using the function defined above, for different parameter sets. Each set of parameters has a weight associated with it, stored in particles_weights_df. For each sample, the script simulates the survival curves for each set of parameters and plots all these curves together on a ggplot graph. On the same plot, the observed survival probabilities are plotted as points. The simulated and observed curves are colored differently for visual distinction.
Saving the plots to a PDF: Each plot generated for a sample is stored in a list. After all samples are processed, the list of plots is arranged in a grid and saved to a PDF file. The script uses the grid.arrange function from the gridExtra package to do this, specifying one column (i.e., the plots are arranged vertically). The height of each subplot is adjusted for optimal viewing.'

## Load your survival data collected in the experiment(s) ----
observed_surv_df <- read_excel(file.path(main_address, "obsSurvival_Temperature_dfs/observed_surv_df.xlsx"))
## Function to simulate a population of LBRs and LBs (MC) ----
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
sample_names <- c('26', '27', '28', '29') # You can adjust these sample names as needed
LTTs =  c(20) 
BR_post_df = post_df[1:1001, c(3:5)]
system.time({
  setwd(outputs_simulation_address)
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
      meanlog10LBR_refT <- BR_post_df[row_n, 3]
      sdlog10LBR_refT <- BR_post_df[row_n, 2]
      k <- BR_post_df[row_n, 1]
      T_LBR_pop_df_wide <- simulate_T_LBR_pop(temperature_sample_df, meanlog10LBR_refT, sdlog10LBR_refT, k, T_ref)
      T_LBR_pop_df_wide = as.data.frame(T_LBR_pop_df_wide)
      all_sim_surv_curves[[paste("meanlog10LBR_refT", meanlog10LBR_refT, "sdlog10LBR_refT", sdlog10LBR_refT, "k", k, sep = "_")]] <- T_LBR_pop_df_wide
      
      # Merge simulated LB and survival data with observed and KM predicted survival data
      sim_obs_full_df <- full_join(T_LBR_pop_df_wide, obs_KM_surv_df, by = c("date_h"), multiple = "all")
      head(T_LBR_pop_df_wide)
      head(obs_KM_surv_df)
      # Filter for mid-day data to match the timing of daily observations
      sim_obs_mid_day_df <- sim_obs_full_df %>% filter(format(date_h, format = "%H:%M:%S") == "12:00:00")
      MAD_df = data.frame(Parametrization = "BR", Experiment = "CHW", sample_name = s_n,
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
  setwd(outputs_simulation_address)
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
sample_names <- c('E2', 'A1', 'B1', 'F2') # DHW
LTTs = c(20) #c(20, 22, 23, 24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5) 
BR_post_df = post_df[1:1001, c(3:5)]
nrow(BR_post_df)
for(LTT in LTTs){
  system.time({
    setwd(outputs_simulation_address)
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
        meanlog10LBR_refT <- BR_post_df[row_n, 3]
        sdlog10LBR_refT <- BR_post_df[row_n, 2]
        k <- BR_post_df[row_n, 1]
        T_LBR_pop_df_wide <- simulate_T_LBR_pop(temperature_sample_df, meanlog10LBR_refT, sdlog10LBR_refT, k, T_ref)
        T_LBR_pop_df_wide = as.data.frame(T_LBR_pop_df_wide)
        all_sim_surv_curves[[paste("meanlog10LBR_refT", meanlog10LBR_refT, "sdlog10LBR_refT", sdlog10LBR_refT, "k", k, sep = "_")]] <- T_LBR_pop_df_wide
        
        # Merge simulated LB and survival data with observed and KM predicted survival data
        sim_obs_full_df <- full_join(T_LBR_pop_df_wide, obs_KM_surv_df, by = c("date_h"), multiple = "all")
        # Filter for mid-day data to match the timing of daily observations
        sim_obs_mid_day_df <- sim_obs_full_df %>% filter(format(date_h, format = "%H:%M:%S") == "12:00:00")
        MAD_df = data.frame(Parametrization = "BR", Experiment = "DHW", sample_name = s_n,
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
    setwd(outputs_simulation_address)
    combined_plot_DHW <- plot_list_DHW[[1]] / plot_list_DHW[[2]] / plot_list_DHW[[3]] / plot_list_DHW[[4]] + 
      plot_layout(guides = 'collect') & theme_publication()
    ggsave("combined_plot_DHW.pdf", combined_plot_DHW, width = plot_width_SingleColumn, height = 18, units = "cm")
  })
  MAD_list_DHW = do.call(rbind, MAD_list_DHW)
}

# MAD df for BR
MAD_BR_df = rbind(MAD_list_CHW, MAD_list_DHW)
write_xlsx(as.data.frame(MAD_BR_df), file.path(main_address, "Outputs", "Bayes_regres", "simulation", "MAD_BR_df.xlsx"))

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



### Compare MAD between BR and ABC-SMC----
library(dplyr)
library(ggplot2)
library(readxl)
library(broom) # for tidy()
library(knitr)
library(patchwork)
library(RColorBrewer)

theme_publication <- function(base_size = 11, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(text = element_text(size = base_size),
          axis.title = element_text(size = rel(1.1)),
          axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)), # adjust top margin for x axis title
          axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)), # adjust right margin for y axis title
          axis.text = element_text(size = rel(1)),
          strip.text = element_text(size = rel(1.1)),
          strip.background = element_blank(),
          panel.spacing = unit(0.3, "lines"), # existing theme elements here...
          legend.position = "right",
          legend.spacing = unit(0.8, 'lines'),
          legend.text = element_text(size = rel(1)),
          legend.title = element_text(size = rel(1.1)),
          legend.box.margin = margin(t = 0, r = 0, b = -5, l = 0, unit = "pt"),
          plot.margin = margin(t = 1, r = 5, b = 1, l = 5, unit = "mm") ) } # adjust the space between legend entries
plot_width_FullPage = 17.3

head(MAD_BR_df)

MAD_ABCSMC_df <- read_excel(file.path(main_address, "Outputs", "ABC_SMC_MAD_withLTT", "MAD_ABCSMC_df.xlsx")) # for us 
head(MAD_ABCSMC_df)

# Assuming you have already read the dataframes MAD_BR_df and MAD_ABCSMC_df

# Combine the dataframes
combined_df <- rbind(MAD_BR_df, MAD_ABCSMC_df)
unique(combined_df$sample_name)

# Replace values in the sample_name column using gsub
combined_df$sample_name <- gsub("A1", "4.5", combined_df$sample_name)
combined_df$sample_name <- gsub("B1", "5", combined_df$sample_name)
combined_df$sample_name <- gsub("E2", "4", combined_df$sample_name)
combined_df$sample_name <- gsub("F2", "5.5", combined_df$sample_name)

# Convert the sample_name column back to a factor if it was one originally
combined_df$sample_name <- factor(combined_df$sample_name)

# Check the first few rows to confirm changes
head(combined_df)


# Check the first few rows to confirm changes
tail(combined_df)


# Define a function to calculate median and quantiles for stat_summary
sum_fun <- function(x) {data.frame(y = median(x),
                                   ymin = quantile(x, probs = 0.05),
                                   ymax = quantile(x, probs = 0.95)) }

# Create a separate ggplot for each experiment
plot_list <- list()
for (exp in unique(combined_df$Experiment)) {
  # Filter data for the current experiment
  data_exp <- combined_df %>% filter(Experiment == exp)
  
  # Find the common y-range for all samples within the experiment
  y_range <- range(data_exp$MAD, na.rm = TRUE)
  
  p <- ggplot(data_exp, aes(x = Parametrization, y = MAD, fill = Parametrization)) +
    geom_violin(trim = FALSE) +
    facet_wrap(~ sample_name) + # Facet by sample_name within each experiment
    scale_fill_brewer(palette = "Set2") + # Color-blind friendly palette
    stat_summary(fun.y = median, geom = "point", size = 1.5, color = "blue") +
    stat_summary(fun.data = sum_fun, geom = "errorbar", width = 0.2, color = "blue") +
    labs(title = exp,
         x = "Parametrization Method",
         y = "MAD") + theme_publication() +
    theme(legend.position = "none") +
    guides(fill = guide_legend(title = "Parametrization")) +
    ylim(y_range) # Set common y-range for all plots within the experiment
  
  plot_list[[exp]] <- p
}

# Combine the plots using patchwork with shared y-axis
combined_plot <- patchwork::wrap_plots(plot_list, ncol = 1)

# Add the legend back to the combined plot, on the right
combined_plot <- combined_plot + 
  plot_layout(guides = "collect") & theme_publication()
combined_plot

file = file.path(main_address, "Outputs", "Bayes_regres", "simulation", "MAD_comparison_BRvsABCSMC_plot.pdf")
ggsave(file, combined_plot, width = plot_width_FullPage, height = 14, units = "cm")


# Calculate the summary statistics
summary_stats_table <- combined_df %>%
  group_by(Parametrization, Experiment) %>%
  summarize(Median_MAD = median(MAD, na.rm = TRUE),
    Q5 = quantile(MAD, probs = 0.05, na.rm = TRUE),
    Q95 = quantile(MAD, probs = 0.95, na.rm = TRUE),
    .groups = 'drop')

# View the table in the console
print(summary_stats_table)

# Create a neat table for reporting
neat_summary_table <- kable(summary_stats_table, format = "html", 
                            caption = "Summary Statistics for MAD by Parametrization and Experiment",
                            align = 'c', digits = 3)

# You can display this neat table in R Markdown directly or save it to an HTML file
neat_summary_table # Displays the table in the R console or R Markdown

# To save the table as an HTML file, you can use:
fileConn <- file("summary_stats_table.html")
writeLines(c('<html><body>', neat_summary_table, '</body></html>'), fileConn)
close(fileConn)
