# Probability and Statistics (2023-2024) [Semester 2]
# Exam Number 1222962
# March 2024

# Assignment: Biomarkers and pain

# References
#
# Thulin, M. (2021). Modern Statistics with R. Eos Chasma Press. ISBN 9789152701515
# Crawley, Michael J. The R Book. Second edition.. Chichester, West Sussex: Wiley, A John Wiley & Sons Ltd., Publicaton, 2013.

################################################################################
#
# Data Preparation
#
################################################################################

# Import the provided biomarker and covariate dataset.
install.packages(c("openxlsx", "data.table"))
library(openxlsx)
library(data.table)

biomarkers <- as.data.table(read.xlsx("biomarkers.xlsx"))
covariates <- as.data.table(read.xlsx("covariates.xlsx"))

# Are there any duplicate observations or covariate patient IDs in the dataset?
# There are no duplicates in the provided dataset, but if there were, they would need to be resolved.
stopifnot(length(unique(biomarkers$Biomarker)) == length(biomarkers$Biomarker))
stopifnot(length(unique(covariates$PatientID)) == length(covariates$PatientID))

# Are there any invalid (negative) observations in the biomarker dataset?
# There are no invalid observations in the provided dataset, but if there were, they would need to be resolved.
stopifnot(FALSE == any(apply(biomarkers[, -1], 2, function(col) (any(col < 0)))))

# From the assignment instructions, columns other than the first column represent individual biomarkers.
# Get the biomarker names.
stopifnot(length(colnames(biomarkers)) > 0)
biomarker_names <- (colnames(biomarkers))[-1]

# From the assignment instructions, the Biomarker column in the biomarkers table is in the format PatientID.timepoint.
# As we checked for duplicates above, we can assume PatientID is a unique numeric per-patient key with which we can index all the other data.
# Split the Biomarker column into PatientID and timepoint columns.
biomarkers[, c("PatientID", "timepoint") := tstrsplit(Biomarker, "-", fixed = TRUE)]

# Reformat PatientID as a numeric and sort the table by PatientID.
biomarkers[, PatientID := as.numeric(PatientID)]
biomarkers[order(PatientID), ]

# Reformat the table with PatientID as the primary key.
biomarkers_patientid <- dcast(biomarkers, PatientID ~ timepoint, value.var = biomarker_names)

# Merge the biomarker data with the covariates, assuming that PatientID can be used as a common key.
# We use a left join such that covariate data for patients without any observations is discarded.
# We will check for missing and invalid observation and covariate data for individual patients as appropriate during subsequent analysis.
biomarkers_covariates_patientid <- merge(biomarkers_patientid, covariates, all.x = TRUE, by = "PatientID")

# Print the number of patients for which we have data.
cat("Data prepared for", nrow(biomarkers_covariates_patientid), "patients.\n")

################################################################################
#
# Statistical Hypothesis Testing
#
################################################################################

# Build single-sex 0 week biomarker tables.
biomarkers_covariates_patientid_male <- biomarkers_covariates_patientid[`Sex.(1=male,.2=female)` == 1, colnames(biomarkers_covariates_patientid)[grepl("0weeks|PatientID", colnames(biomarkers_covariates_patientid))], with = FALSE]
biomarkers_covariates_patientid_female <- biomarkers_covariates_patientid[`Sex.(1=male,.2=female)` == 2, colnames(biomarkers_covariates_patientid)[grepl("0weeks|PatientID", colnames(biomarkers_covariates_patientid))], with = FALSE]

# Remove patients with any missing observations.
biomarkers_covariates_patientid_male <- biomarkers_covariates_patientid_male[rowSums(is.na(biomarkers_covariates_patientid_male)) == 0, ]
biomarkers_covariates_patientid_female <- biomarkers_covariates_patientid_female[rowSums(is.na(biomarkers_covariates_patientid_female)) == 0, ]

# Print the number of patients for which we have data.
cat("Data prepared for", nrow(biomarkers_covariates_patientid_male), "males.\n")
cat("Data prepared for", nrow(biomarkers_covariates_patientid_female), "females.\n")

histlab = 'Biomarker Level'

# Plot histograms for each male biomarker.
par(mfrow = c(3,3))
hist(biomarkers_covariates_patientid_male$`IL-8_0weeks`, main = "IL-8 Levels At Inclusion (Male)", xlab = histlab)
hist(biomarkers_covariates_patientid_male$`VEGF-A_0weeks`, main = "VEGF Levels At Inclusion (Male)", xlab = histlab)
hist(biomarkers_covariates_patientid_male$`OPG_0weeks`, main = "OPG Levels At Inclusion (Male)", xlab = histlab)
hist(biomarkers_covariates_patientid_male$`TGF-beta-1_0weeks`, main = "TGF-beta-1 Levels At Inclusion (Male)", xlab = histlab)
hist(biomarkers_covariates_patientid_male$`IL-6_0weeks`, main = "IL-6 Levels At Inclusion (Male)", xlab = histlab)
hist(biomarkers_covariates_patientid_male$`CXCL9_0weeks`, main = "CXCL9 Levels At Inclusion (Male)", xlab = histlab)
hist(biomarkers_covariates_patientid_male$`CXCL1_0weeks`, main = "CXCL1 Levels At Inclusion (Male)", xlab = histlab)
hist(biomarkers_covariates_patientid_male$`IL-18_0weeks`, main = "IL-18 Levels At Inclusion (Male)", xlab = histlab)
hist(biomarkers_covariates_patientid_male$`CSF-1_0weeks`, main = "CSF-1 Levels At Inclusion (Male)", xlab = histlab)

# Plot histograms for each female biomarker.
par(mfrow = c(3,3))
hist(biomarkers_covariates_patientid_female$`IL-8_0weeks`, main = "IL-8 Levels At Inclusion (Female)", xlab = histlab)
hist(biomarkers_covariates_patientid_female$`VEGF-A_0weeks`, main = "VEGF Levels At Inclusion (Female)", xlab = histlab)
hist(biomarkers_covariates_patientid_female$`OPG_0weeks`, main = "OPG Levels At Inclusion (Female)", xlab = histlab)
hist(biomarkers_covariates_patientid_female$`TGF-beta-1_0weeks`, main = "TGF-beta-1 Levels At Inclusion (Female)", xlab = histlab)
hist(biomarkers_covariates_patientid_female$`IL-6_0weeks`, main = "IL-6 Levels At Inclusion (Female)", xlab = histlab)
hist(biomarkers_covariates_patientid_female$`CXCL9_0weeks`, main = "CXCL9 Levels At Inclusion (Female)", xlab = histlab)
hist(biomarkers_covariates_patientid_female$`CXCL1_0weeks`, main = "CXCL1 Levels At Inclusion (Female)", xlab = histlab)
hist(biomarkers_covariates_patientid_female$`IL-18_0weeks`, main = "IL-18 Levels At Inclusion (Female)", xlab = histlab)
hist(biomarkers_covariates_patientid_female$`CSF-1_0weeks`, main = "CSF-1 Levels At Inclusion (Female)", xlab = histlab)

# Create a results table.
hyptest_results <- matrix(c(1:54), ncol = 6)
colnames(hyptest_results) <- c('Sample Mean (Male)', 'Sample Mean (Female)', 'p-value H1', 'p-value H2', 'H1 rejects H0 (5%)', 'H2 rejects H0 (5%)')
rownames(hyptest_results) <- biomarker_names

# Compute the sample means for the male and female levels.
hyptest_results[1, 1] <- mean(biomarkers_covariates_patientid_male$`IL-8_0weeks`)
hyptest_results[2, 1] <- mean(biomarkers_covariates_patientid_male$`VEGF-A_0weeks`)
hyptest_results[3, 1] <- mean(biomarkers_covariates_patientid_male$`OPG_0weeks`)
hyptest_results[4, 1] <- mean(biomarkers_covariates_patientid_male$`TGF-beta-1_0weeks`)
hyptest_results[5, 1] <- mean(biomarkers_covariates_patientid_male$`IL-6_0weeks`)
hyptest_results[6, 1] <- mean(biomarkers_covariates_patientid_male$`CXCL9_0weeks`)
hyptest_results[7, 1] <- mean(biomarkers_covariates_patientid_male$`CXCL1_0weeks`)
hyptest_results[8, 1] <- mean(biomarkers_covariates_patientid_male$`IL-18_0weeks`)
hyptest_results[9, 1] <- mean(biomarkers_covariates_patientid_male$`CSF-1_0weeks`)

hyptest_results[1, 2] <- mean(biomarkers_covariates_patientid_female$`IL-8_0weeks`)
hyptest_results[2, 2] <- mean(biomarkers_covariates_patientid_female$`VEGF-A_0weeks`)
hyptest_results[3, 2] <- mean(biomarkers_covariates_patientid_female$`OPG_0weeks`)
hyptest_results[4, 2] <- mean(biomarkers_covariates_patientid_female$`TGF-beta-1_0weeks`)
hyptest_results[5, 2] <- mean(biomarkers_covariates_patientid_female$`IL-6_0weeks`)
hyptest_results[6, 2] <- mean(biomarkers_covariates_patientid_female$`CXCL9_0weeks`)
hyptest_results[7, 2] <- mean(biomarkers_covariates_patientid_female$`CXCL1_0weeks`)
hyptest_results[8, 2] <- mean(biomarkers_covariates_patientid_female$`IL-18_0weeks`)
hyptest_results[9, 2] <- mean(biomarkers_covariates_patientid_female$`CSF-1_0weeks`)

# Compute the H1 p-value for each biomarker and populate the results table.
hyptest_results[1, 3] <- t.test(biomarkers_covariates_patientid_male$`IL-8_0weeks`, biomarkers_covariates_patientid_female$`IL-8_0weeks`, alternative="greater")$p.value
hyptest_results[2, 3] <- t.test(biomarkers_covariates_patientid_male$`VEGF-A_0weeks`, biomarkers_covariates_patientid_female$`VEGF-A_0weeks`, alternative="greater")$p.value
hyptest_results[3, 3] <- t.test(biomarkers_covariates_patientid_male$`OPG_0weeks`, biomarkers_covariates_patientid_female$`OPG_0weeks`, alternative="greater")$p.value
hyptest_results[4, 3] <- t.test(biomarkers_covariates_patientid_male$`TGF-beta-1_0weeks`, biomarkers_covariates_patientid_female$`TGF-beta-1_0weeks`, alternative="greater")$p.value
hyptest_results[5, 3] <- t.test(biomarkers_covariates_patientid_male$`IL-6_0weeks`, biomarkers_covariates_patientid_female$`IL-6_0weeks`, alternative="greater")$p.value
hyptest_results[6, 3] <- t.test(biomarkers_covariates_patientid_male$`CXCL9_0weeks`, biomarkers_covariates_patientid_female$`CXCL9_0weeks`, alternative="greater")$p.value
hyptest_results[7, 3] <- t.test(biomarkers_covariates_patientid_male$`CXCL1_0weeks`, biomarkers_covariates_patientid_female$`CXCL1_0weeks`, alternative="greater")$p.value
hyptest_results[8, 3] <- t.test(biomarkers_covariates_patientid_male$`IL-18_0weeks`, biomarkers_covariates_patientid_female$`IL-18_0weeks`, alternative="greater")$p.value
hyptest_results[9, 3] <- t.test(biomarkers_covariates_patientid_male$`CSF-1_0weeks`, biomarkers_covariates_patientid_female$`CSF-1_0weeks`, alternative="greater")$p.value

# Compute the H2 p-value for each biomarker and populate the results table.
hyptest_results[1, 4] <- t.test(biomarkers_covariates_patientid_male$`IL-8_0weeks`, biomarkers_covariates_patientid_female$`IL-8_0weeks`, alternative="less")$p.value
hyptest_results[2, 4] <- t.test(biomarkers_covariates_patientid_male$`VEGF-A_0weeks`, biomarkers_covariates_patientid_female$`VEGF-A_0weeks`, alternative="less")$p.value
hyptest_results[3, 4] <- t.test(biomarkers_covariates_patientid_male$`OPG_0weeks`, biomarkers_covariates_patientid_female$`OPG_0weeks`, alternative="less")$p.value
hyptest_results[4, 4] <- t.test(biomarkers_covariates_patientid_male$`TGF-beta-1_0weeks`, biomarkers_covariates_patientid_female$`TGF-beta-1_0weeks`, alternative="less")$p.value
hyptest_results[5, 4] <- t.test(biomarkers_covariates_patientid_male$`IL-6_0weeks`, biomarkers_covariates_patientid_female$`IL-6_0weeks`, alternative="less")$p.value
hyptest_results[6, 4] <- t.test(biomarkers_covariates_patientid_male$`CXCL9_0weeks`, biomarkers_covariates_patientid_female$`CXCL9_0weeks`, alternative="less")$p.value
hyptest_results[7, 4] <- t.test(biomarkers_covariates_patientid_male$`CXCL1_0weeks`, biomarkers_covariates_patientid_female$`CXCL1_0weeks`, alternative="less")$p.value
hyptest_results[8, 4] <- t.test(biomarkers_covariates_patientid_male$`IL-18_0weeks`, biomarkers_covariates_patientid_female$`IL-18_0weeks`, alternative="less")$p.value
hyptest_results[9, 4] <- t.test(biomarkers_covariates_patientid_male$`CSF-1_0weeks`, biomarkers_covariates_patientid_female$`CSF-1_0weeks`, alternative="less")$p.value

# Calculate whether H0 was rejected by each alternative hypothesis at the 5% significance level.
for (i in 1:2)
{
  for (j in 1:9)
  {
    hyptest_results[j, i + 4] = 'No'

    if (hyptest_results[j, i + 2] < 0.05)
    {
      hyptest_results[j, i + 4] = 'Yes'
    }
  }
}

View(hyptest_results)

# Print the family-wise error rate assuming 18 independent tests at the 5% significance level.
cat("FWER:", 1 - (1 - 0.05)^18, "\n")

# Duplicate the results table and recalculate H0 rejection at the 5% significance level with Bonferroni correction.
hyptest_results_bonferroni <- hyptest_results
colnames(hyptest_results_bonferroni) <- c('Sample Mean (Male)', 'Sample Mean (Female)', 'p-value H1', 'p-value H2', 'H1 rejects H0 (5% / m)', 'H2 rejects H0 (5% / m)')

for (i in 1:2)
{
  for (j in 1:9)
  {
    hyptest_results_bonferroni[j, i + 4] = 'No'
    
    if (hyptest_results_bonferroni[j, i + 2] < (0.05 / 18))
    {
      hyptest_results_bonferroni[j, i + 4] = 'Yes'
    }
  }
}

View(hyptest_results_bonferroni)

# Print the Bonferroni-corrected family-wise error rate assuming 18 independent tests at the 5% significance level.
cat("Bonferroni-corrected FWER:", 1 - (1 - (0.05 / 18))^18, "\n")

################################################################################
#
# Regression Modelling
#
################################################################################

# Build regression (PatientID, levels at inclusion, covariates) table.
biomarkers_covariates_patientid_regression <- biomarkers_covariates_patientid[, colnames(biomarkers_covariates_patientid)[grepl("0weeks|PatientID|Age|Sex|Smoker|VAS|Vas", colnames(biomarkers_covariates_patientid))], with = FALSE]

# Remove patients with any missing observations.
biomarkers_covariates_patientid_regression <- biomarkers_covariates_patientid_regression[rowSums(is.na(biomarkers_covariates_patientid_regression)) == 0, ]

# Plot the 12-month VAS dependence on biomarker levels and covariates.
attach(biomarkers_covariates_patientid_regression)
par(mfrow=c(3, 3))
plot(`IL-8_0weeks`, `Vas-12months`)
m <- lm(`Vas-12months` ~ `IL-8_0weeks`)
abline(coef(m), col = 2)

plot(`VEGF-A_0weeks`, `Vas-12months`)
m <- lm(`Vas-12months` ~ `VEGF-A_0weeks`)
abline(coef(m), col = 2)

plot(`OPG_0weeks`, `Vas-12months`)
m <- lm(`Vas-12months` ~ `OPG_0weeks`)
abline(coef(m), col = 2)

plot(`TGF-beta-1_0weeks`, `Vas-12months`)
m <- lm(`Vas-12months` ~ `TGF-beta-1_0weeks`)
abline(coef(m), col = 2)

plot(`IL-6_0weeks`,`Vas-12months`)
m <- lm(`Vas-12months` ~ `IL-6_0weeks`)
abline(coef(m), col = 2)

plot(`CXCL9_0weeks`,`Vas-12months`)
m <- lm(`Vas-12months` ~ `CXCL9_0weeks`)
abline(coef(m), col = 2)

plot(`CXCL1_0weeks`,`Vas-12months`)
m <- lm(`Vas-12months` ~ `CXCL1_0weeks`)
abline(coef(m), col = 2)

plot(`IL-18_0weeks`,`Vas-12months`)
m <- lm(`Vas-12months` ~ `IL-18_0weeks`)
abline(coef(m), col = 2)

plot(`CSF-1_0weeks`,`Vas-12months`)
m <- lm(`Vas-12months` ~ `CSF-1_0weeks`)
abline(coef(m), col = 2)

plot(`Age`,`Vas-12months`)
m <- lm(`Vas-12months` ~ `Age`)
abline(coef(m), col = 2)

plot(`Sex.(1=male,.2=female)`,`Vas-12months`)
m <- lm(`Vas-12months` ~ `Sex.(1=male,.2=female)`)
abline(coef(m), col = 2)

plot(`Smoker.(1=yes,.2=no)`,`Vas-12months`)
m <- lm(`Vas-12months` ~ `Smoker.(1=yes,.2=no)`)
abline(coef(m), col = 2)

plot(`VAS-at-inclusion`,`Vas-12months`)
m <- lm(`Vas-12months` ~ `VAS-at-inclusion`)
abline(coef(m), col = 2)
detach()

# Split the table 80/20.
biomarkers_covariates_patientid_regression_80 <- biomarkers_covariates_patientid_regression[c(1:(nrow(biomarkers_covariates_patientid_regression) * 0.8)), colnames(biomarkers_covariates_patientid_regression), with = FALSE]
biomarkers_covariates_patientid_regression_20 <- biomarkers_covariates_patientid_regression[c(((nrow(biomarkers_covariates_patientid_regression) * 0.8) + 1):(nrow(biomarkers_covariates_patientid_regression))), colnames(biomarkers_covariates_patientid_regression), with = FALSE]

# Print the number of patients for which we have data.
cat("80% regression data prepared for", nrow(biomarkers_covariates_patientid_regression_80), "patients.\n")
cat("20% regression data prepared for", nrow(biomarkers_covariates_patientid_regression_20), "patients.\n")

# Using 80% of the patient data, fit a multiple linear regression model with 12-month VAS as the response variable and the biomarker levels at inclusion and covariates as explanatory variables.
multiple_regression_model <- lm(`Vas-12months` ~ `IL-8_0weeks` + `VEGF-A_0weeks` + `OPG_0weeks` + `TGF-beta-1_0weeks` + `IL-6_0weeks` + `CXCL9_0weeks` + `CXCL1_0weeks` + `IL-18_0weeks` + `CSF-1_0weeks` + `Age` + `Sex.(1=male,.2=female)` + `Smoker.(1=yes,.2=no)` + `VAS-at-inclusion`, biomarkers_covariates_patientid_regression_80)
summary(multiple_regression_model)$r.squared

# Add the fitted values and residuals as new columns in the 80% data table.
biomarkers_covariates_patientid_regression_80[, `Vas-12months-fitted` := predict(multiple_regression_model)]
biomarkers_covariates_patientid_regression_80[, `Vas-12months-residuals` := residuals(multiple_regression_model)]

# Plot the fitted values and residuals against actual 12-month VAS.
attach(biomarkers_covariates_patientid_regression_80)
par(mfrow=c(1, 2))
plot(`Vas-12months`, `Vas-12months-fitted`)
plot(`Vas-12months`, `Vas-12months-residuals`)
detach()

# Predict 12 month VAS for the remaining 20% of the patient data and add them as a new column.
biomarkers_covariates_patientid_regression_20[, `Vas-12months-predicted` := predict(multiple_regression_model, newdata=biomarkers_covariates_patientid_regression_20)]

# Calculate residuals for the predicted 12 month VAS and add them as a new column.
biomarkers_covariates_patientid_regression_20[, `Vas-12months-residuals`:= biomarkers_covariates_patientid_regression_20$`Vas-12months` - biomarkers_covariates_patientid_regression_20$`Vas-12months-predicted`]

# Plot the fitted values and residuals against actual 12-month VAS.
attach(biomarkers_covariates_patientid_regression_20)
par(mfrow=c(1, 2))
plot(`Vas-12months`, `Vas-12months-predicted`)
plot(`Vas-12months`, `Vas-12months-residuals`)
detach()

# Calculate the prediction and confidence intervals.
predicted_vas_prediction <- predict(multiple_regression_model, int="prediction", newdata = biomarkers_covariates_patientid_regression_20)
predicted_vas_confidence <- predict(multiple_regression_model, int="confidence", newdata = biomarkers_covariates_patientid_regression_20)

View(predicted_vas_confidence)
View(predicted_vas_prediction)
