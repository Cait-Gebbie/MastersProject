library(ggplot2)
library(scales)
library(readr)

# ---- 1. Read the .het file ----
# vcftools --het output is whitespace-delimited with a header:
# FID  IID  O(HOM)  E(HOM)  N(NM)  F
het <- read.table("C:/Users/caitg/OneDrive/Documents/Work/Masters/Lion_Sequences/ROH/Final/hetGVCF_urelated.het", header = TRUE, stringsAsFactors = FALSE)

# Some vcftools versions name the F column "F", others may quote it —
# this keeps things safe regardless of exact column name formatting
names(het)[names(het) == "F"] <- "FHOM"

# ---- 2. Plot F_HOM per individual, same structure as the ROH plot ----
p <- ggplot(het, aes(x = "", y = FHOM)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0, alpha = 0.3, size = 2) +
  theme(plot.title = element_text(size = 10),
        aspect.ratio = 1.5,
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey95")) +
  scale_y_continuous(labels = comma) +
  labs(x = "",
       y = expression(F[HOM]),
       title = expression(paste("Distribution of ", F[HOM], " per Individual")))

print(p)

# ---- 3. Save ----
ggsave("FHOM_distribution.png", plot = p, width = 5, height = 6, dpi = 300)
