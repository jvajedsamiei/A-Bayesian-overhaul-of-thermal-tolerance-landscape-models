# Load necessary libraries
library(dplyr)


MAD_BR_df <- read_excel("Library/Mobile Documents/com~apple~CloudDocs/Jahangir Air/Postdoc_Position:Propos:RelatedStuff/PROPOSAL_postdoc/Selection_docs/Niklas_Mytilus_experiment/Outputs/Bayes_regres/simulation/MAD_BR_df.xlsx")
View(MAD_BR_df)

MAD_ABCSMC_df <- read_excel("Library/Mobile Documents/com~apple~CloudDocs/Jahangir Air/Postdoc_Position:Propos:RelatedStuff/PROPOSAL_postdoc/Selection_docs/Niklas_Mytilus_experiment/Outputs/ABC_SMC_MAD_withLTT/MAD_ABCSMC_df.xlsx")
View(MAD_ABCSMC_df)


data = MAD_BR_df
data = MAD_ABCSMC_df

# Group by experiment, meanlog10LBR_refT, sdlog10LBR_refT, and k and calculate mean MAD
grouped_data <- data %>%
  group_by(Experiment, meanlog10LBR_refT, sdlog10LBR_refT, k) %>%
  summarise(mean_MAD = mean(MAD, na.rm = TRUE))

# Now, for each experiment, calculate the median, 5th and 95th percentiles of mean MAD
percentile_data <- grouped_data %>%
  group_by(Experiment) %>%
  summarise(
    median_MAD = median(mean_MAD, na.rm = TRUE),
    perc_5 = quantile(mean_MAD, 0.05, na.rm = TRUE),
    perc_95 = quantile(mean_MAD, 0.95, na.rm = TRUE)
  )

# Display the results
print(grouped_data)
print(percentile_data)


