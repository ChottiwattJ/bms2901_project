library(haven)
library(sas7bdat)
library(tidyverse)
library(haven)
library(lubridate)
library(broom)
library(car)
library(dplyr)
library(nortest)
library(ROSE)
library(dplyr)
library(logistf)
library(ggplot2)
library(knitr)

####################
# Data Preparation #
####################

# Read NHANES data
demographics <- read_xpt("Dataset/DEMO_J.XPT")
hep <- read_xpt("Dataset/HEPC_J.XPT")
liver <- read_xpt("Dataset/MCQ_J.XPT")

# Merge datasets
merged_data <- demographics %>%
  left_join(hep, by = "SEQN") %>%
  left_join(liver, by = "SEQN")

# Extract and recode variables 
clean_data <- merged_data %>%
  select(SEQN, age = RIDAGEYR, sex = RIAGENDR, race = RIDRETH1,
         hepatitis_c = LBXHCR, 
         liver_condition = MCQ160L, 
         genotype = LBXHCG) %>% # R
  mutate(hepatitis_c = ifelse(hepatitis_c == 1, 1, 0),
         liver_condition = ifelse(liver_condition == 1, 1, 0))

clean_data <- clean_data %>%
  mutate(liver_condition = ifelse(is.na(liver_condition), 0, liver_condition),
         genotype = ifelse(hepatitis_c == 1 & is.na(genotype), 0, genotype), hepatitis_c = ifelse(is.na(hepatitis_c), 0, hepatitis_c))

####################
# Analysis Section #
####################

# Anderson-Darling Goodness-of-Fit Test

hep_c_pos <- clean_data %>% filter(hepatitis_c == 1)
hep_c_neg <- clean_data %>% filter(hepatitis_c == 0)

ad_test_pos <- ad.test(hep_c_pos$liver_condition)
ad_test_neg <- ad.test(hep_c_neg$liver_condition)

clean_data$hepatitis_c_factor <- as.factor(clean_data$hepatitis_c)

# Levene Equality of Variances Test
levene_test <- car::leveneTest(clean_data$liver_condition, clean_data$hepatitis_c_factor, center = median)

# Kruskal-Wallis One-way Analysis of Variance
kruskal_test <- kruskal.test(liver_condition ~ hepatitis_c, data = clean_data)

# Data Preparation for Oversampling
input_data <- clean_data %>%
  select(hepatitis_c, liver_condition, age, race, sex) %>%
  mutate(age_group = cut(age, breaks = seq(0, 100, by = 20), right = FALSE, include.lowest = TRUE))

# SMOTE Oversampling Parameters
oversampling_ratio <- 0.05  # Amount of oversampling
K <- 5  # Number of nearest neighbors

# Data Oversampling
oversampled_data <- ovun.sample(hepatitis_c ~ ., data = input_data, method = "over", N = nrow(input_data) * (1 + oversampling_ratio), seed = 50295022)$data

# Convert sex to factor
oversampled_data$sex <- as.factor(oversampled_data$sex)

# Convert age_group to factor
oversampled_data$age_group <- as.factor(oversampled_data$age_group)

# Convert race to factor
oversampled_data$race <- as.factor(oversampled_data$race)

# Fitting Logistic Regression
logistic_model <- glm(liver_condition ~ hepatitis_c + sex + age_group + race, data = oversampled_data, family = binomial(link = 'logit'))

# Summarize the results
summary(logistic_model)

# Validation by Firth's Bias-Reduced Logistic Regression
firth_model <- logistf(liver_condition ~ hepatitis_c + sex + age_group + race, data = oversampled_data, family = binomial(link = 'logit'))
summary(firth_model)

#################
# Visualization #
#################

# Calculate predicted probabilities
predicted_probs <- predict(logistic_model, type = "response")

# Add predicted probabilities to the dataset
oversampled_data$predicted_probs <- predicted_probs

ggplot(oversampled_data, aes(x = hepatitis_c, y =  predicted_probs)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, linetype = "solid", color = "blue") +
  labs(
    title = "Logistic Regression Plot",
    y = "Hepatitis C Status (0 = No, 1 = Yes)",
    x = "Predicted Probability of Liver Condition"
  ) +
  theme_minimal()
 
# Original data
ggplot(clean_data, aes(x = factor(1), fill = factor(liver_condition))) +
  geom_bar(position = "fill") +
  coord_flip() +
  labs(title = "Distribution of Liver Condition in Original Data",
       x = NULL, y = "Proportion") +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_discrete(guide = FALSE)

ggplot(clean_data, aes(x = factor(1), fill = factor(hepatitis_c))) +
  geom_bar(position = "fill") +
  coord_flip() +
  labs(title = "Distribution of Hepatitis C in Original Data",
       x = NULL, y = "Proportion") +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_discrete(guide = FALSE)

ggplot(oversampled_data, aes(x = factor(1), fill = factor(liver_condition))) +
  geom_bar(position = "fill") +
  coord_flip() +
  labs(title = "Distribution of Liver Condition in Oversampled Data",
       x = NULL, y = "Proportion") +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_discrete(guide = FALSE)

ggplot(oversampled_data, aes(x = factor(1), fill = factor(hepatitis_c))) +
  geom_bar(position = "fill") +
  coord_flip() +
  labs(title = "Distribution of Hepatitis C in Oversampled Data",
       x = NULL, y = "Proportion") +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_discrete(guide = FALSE)


# Calculate summary statistics
sex_summary <- clean_data %>%
  group_by(sex) %>%
  summarize(
    n = n(),
    n_pct = n() / nrow(clean_data) * 100,
    hep_c_count = sum(hepatitis_c),
    hep_c_pct = hep_c_count / nrow(clean_data) * 100,
    liver_condition_count = sum(liver_condition),
    liver_condition_pct = liver_condition_count / nrow(clean_data) * 100
  )

race_summary <- clean_data %>%
  group_by(race) %>%
  summarize(
    n = n(),
    n_pct = n() / nrow(clean_data) * 100,
    hep_c_count = sum(hepatitis_c),
    hep_c_pct = hep_c_count / nrow(clean_data) * 100,
    liver_condition_count = sum(liver_condition),
    liver_condition_pct = liver_condition_count / nrow(clean_data) * 100
  )

total_mean_age <- mean(clean_data$age, na.rm = TRUE)
total_age_sd <- sd(clean_data$age, na.rm = TRUE)

# Print summary statistics
cat("Sex Summary:\n")
kable(sex_summary)

cat("\nRace Summary:\n")
kable(race_summary)

cat("\nTotal Mean Age:", total_mean_age)
cat("\nTotal Age SD:", total_age_sd)

