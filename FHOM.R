library(ggplot2)
library(scales)
library(readr)



het <- read.table("C:/Users/caitg/OneDrive/Documents/Work/Masters/Lion_Sequences/ROH/Final/hetGVCF_urelated.het", header = TRUE, stringsAsFactors = FALSE)


names(het)[names(het) == "F"] <- "FHOM"


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
