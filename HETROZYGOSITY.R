library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(readxl)

het <- read.table("C:/Users/caitg/OneDrive/Documents/Work/Masters/Results/OutputFiles/NoWhiteLion/hetGVCF_unrelated.het", header = TRUE)


het <- het %>%
  mutate(
    HO = 1 - (O.HOM. / N.NM.),
    HE = 1 - (E.HOM. / N.NM.),
    Delta_HO_HE = HO - HE
  )

het$Delta_HO_HE


het[, c("HO", "HE", "Delta_HO_HE")]
het %>% select( HO, HE, Delta_HO_HE)

print(het)
wilcox.test(het$HO, het$HE, paired = TRUE, alternative = "greater")
ggplot(het_long, aes(x = Heterozygosity, fill = Type)) +
  geom_density(alpha = 0.5) +
  scale_color_brewer(palette = "Greys") +
  theme_classic() +
  labs(title = "Density of Observed and Expected Heterozygosity",
       x = "Heterozygosity",
       y = "Density")

n        <- length(het$Delta_HO_HE)
mean_val <- mean(het$Delta_HO_HE)
sd_val   <- sd(het$Delta_HO_HE)
se_val   <- sd_val / sqrt(n)


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





