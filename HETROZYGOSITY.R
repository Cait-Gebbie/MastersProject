library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(readxl)

# 1️⃣ Load PLINK .het file
het <- read.table("C:/Users/caitg/OneDrive/Documents/Work/Masters/Results/OutputFiles/NoWhiteLion/hetGVCF_unrelated.het", header = TRUE)

# Calculate HO, HE, Delta
het <- het %>%
  mutate(
    HO = 1 - (O.HOM. / N.NM.),
    HE = 1 - (E.HOM. / N.NM.),
    Delta_HO_HE = HO - HE
  )
# Print just the Delta column
het$Delta_HO_HE

# Print with individual ID alongside it (assuming an ID column, e.g. INDV)
het[, c("HO", "HE", "Delta_HO_HE")]

# Or with dplyr, more readable
het %>% select( HO, HE, Delta_HO_HE)

# View the whole updated data frame
print(het)
wilcox.test(het$HO, het$HE, paired = TRUE, alternative = "greater")
ggplot(het_long, aes(x = Heterozygosity, fill = Type)) +
  geom_density(alpha = 0.5) +
  scale_color_brewer(palette = "Greys") +
  theme_classic() +
  labs(title = "Density of Observed and Expected Heterozygosity",
       x = "Heterozygosity",
       y = "Density")
# Combine HO and HE manually

# Basic stats
n        <- length(het$Delta_HO_HE)
mean_val <- mean(het$Delta_HO_HE)
sd_val   <- sd(het$Delta_HO_HE)
se_val   <- sd_val / sqrt(n)

# 95% CI using the t-distribution (better than normal approx for smaller n)
t_crit <- qt(0.975, df = n - 1)
ci_lower <- mean_val - t_crit * se_val
ci_upper <- mean_val + t_crit * se_val

cat("n =", n, "\n")
cat("Mean =", mean_val, "\n")
cat("SD =", sd_val, "\n")
cat("SE =", se_val, "\n")
cat("95% CI = [", ci_lower, ",", ci_upper, "]\n")
het_long <- data.frame(
  Heterozygosity = c(het$HO, het$HE),
  Type = rep(c("HO", "HE"), each = nrow(het))
)
ggplot(het_long, aes(x = Heterozygosity, fill = Type)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("HO" = "#94D2BD", "HE" = "#0A9396")) +
  theme_classic() +
  labs(title = "Density of Observed and Expected Heterozygosity",
       x = "Heterozygosity",
       y = "Density")
het_obs <- het_long %>%
  filter(Type == "Observed")

ggplot(het_obs, aes(x = het_obs, y = Heterozygosity)) +
  geom_boxplot(outlier.shape = NA, fill = "grey80") +
  geom_jitter(width = 0, alpha = 0.5, size = 1) +
  theme_classic() +
  labs(title = "Observed Heterozygosity per Individual",
       x = "",
       y = "Observed Heterozygosity")
ggplot(het_long, aes(x = Type, y = Heterozygosity, fill = Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  theme(legend.position = "none") +
  labs(title = "Observed vs Expected Heterozygosity per Individual",
       x = "",
       y = "Heterozygosity")

ggplot(het, aes(x = Delta_HO_HE)) +
  geom_density(fill = "#005F73", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Set2") +
  theme_classic() +
  labs(title = "Density of Δ(HO - HE)",
       x = expression(Delta~(H[O] - H[E])),
       y = "Density")

ggplot(het, aes(x = "", y = Delta_HO_HE)) +
  geom_boxplot(fill = "grey70") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "",
       y = expression(Delta~(H[O] - H[E])),
       title = "Deviation from Expected Heterozygosity")


sex_data <- read_excel("C:/Users/caitg/OneDrive/Documents/Work/Masters/Lion_Sequences/ROH/Sexdata.xlsx")
colnames(het_obs)    

colnames(het)
head(het)
het_sex <- het_obs %>%
  left_join(sex_data, by = c("IID" = "ID"))

ggplot(het_sex, aes(x = SEX, y = Heterozygosity, fill = SEX)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = c("Male" = "grey40", "Female" = "grey80")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Observed Heterozygosity by Sex",
       x = "Sex",
       y = "Observed Heterozygosity")

sex_data <- read_excel("C:/Users/caitg/OneDrive/Documents/Work/Masters/Lion_Sequences/ROH/Sexdata.xlsx")

# Merge observed heterozygosity with sex data
het_sex <- het %>%
  left_join(sex_data, by = c("IID" = "ID"))

ggplot(het_sex, aes(x = SEX, y = HE, fill = SEX)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = c("Male" = "grey40", "Female" = "grey80")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Expected Heterozygosity by Sex",
       x = "Sex",
       y = "Observed Heterozygosity (HO)")

shapiro.test(het_sex$HO[het_sex$SEX == "Male"])
shapiro.test(het_sex$HO[het_sex$SEX == "Female"])
t.test(HO ~ SEX, data = het_sex)
wilcox.test(HO ~ SEX, data = het_sex)


