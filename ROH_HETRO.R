
library(ggplot2)
library(dplyr)

roh <- read.table("C:/Users/caitg/OneDrive/Documents/Work/Masters/Results/OutputFiles/NoWhiteLion/filteredvcf_results.hom", header=TRUE)

roh_indiv <- read.table("C:/Users/caitg/OneDrive/Documents/Work/Masters/Results/OutputFiles/NoWhiteLion/filteredvcf_results.hom.indiv", header=TRUE)


roh$LengthKb <- roh$KB


ggplot(roh, aes(x = LengthKb)) +
  geom_histogram(binwidth = 250, colour = "white", linewidth = 0.3,fill = "#005F73") +
  theme_classic() +
  labs(x = "ROH Length (Kb)",
       y = "Count",
       title = "Distribution of ROH Lengths")
roh_sex <- roh %>%
  left_join(sex_data, by = c("IID" = "ID"))

ggplot(roh_sex, aes(x = LengthKb)) +
  geom_histogram(binwidth = 500, colour = "white", linewidth = 0.3, fill = "#005F73") +
  facet_wrap(~ SEX) +
  theme_classic() +
  labs(x = "ROH Length (Kb)",
       y = "Count",
       title = "Distribution of ROH Lengths by Sex")


ggplot(roh, aes(x = IID, y = LengthKb)) +
  geom_boxplot() +
  theme_classic() +
  scale_y_continuous(labels = comma) +
  theme(axis.text.x = element_blank()) +
  labs(x = "Individuals",
       y = "ROH Length (Kb)",
       title = "Distribution of ROH Length per Individual")

ggplot(roh_sex, aes(x = LengthKb, fill = SEX)) +
  geom_histogram(binwidth = 500, colour = "white", linewidth = 0.3, 
                 position = "identity", alpha = 0.5) +
  scale_fill_manual(values = c("Male" = "blue", "Female" = "pink")) +
  theme_classic() +
  labs(x = "ROH Length (Kb)",
       y = "Count",
       title = "Distribution of ROH Lengths by Sex")
# Total ROH per individual
total_roh <- roh %>%
  group_by(IID) %>%
  summarise(TotalROH_Kb = sum(LengthKb),
            N_ROH = n())


ggplot(total_roh, aes(x = reorder(IID, TotalROH_Kb), y = TotalROH_Kb)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_y_continuous(labels = comma) +
  theme(axis.text.x = element_blank()) +
  labs(x = "Individuals",
       y = "Total ROH (Kb)",
       title = "Individual Length")


ggplot(total_roh, aes(x = "", y = TotalROH_Kb)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0, alpha = 0.3, size = 2) +
  scale_y_continuous(labels = comma) +
  theme(plot.title = element_text(size = 10),
        aspect.ratio = 1.5,
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey95")) +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  labs(x = "",
       y = "Total ROH (Kb)",
       title = "Distribution of Total ROH per Individual")


ggplot(total_roh, aes(x = reorder(IID, TotalROH_Kb), y = TotalROH_Kb)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_y_continuous(labels = comma) +
  theme(axis.text.x = element_blank()) +
  labs(x = "Individuals",
       y = "Total ROH (Kb)",
       title = "Individual Length")
avg_roh <- roh %>%
  group_by(IID) %>%
  summarise(MeanROH_Kb = mean(LengthKb))


ggplot(avg_roh, aes(x = "", y = MeanROH_Kb)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(labels = comma) +
  geom_jitter(width = 0, alpha = 0.3, size= 2) +
  theme(plot.title = element_text(size = 10),
        aspect.ratio = 1.5,
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey95")) +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  theme(plot.title = element_text(size = 12),
        aspect.ratio = 1.5) +
  labs(x = "",
       y = "Mean ROH Length (Kb)",
       title = "Average ROH Length per Individual")


genome_length_kb <- 2.4e6   


total_roh$FROH <- total_roh$TotalROH_Kb / genome_length_kb


p <- ggplot(total_roh, aes(x = "", y = FROH)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0, alpha = 0.3, size = 2) +
  theme(plot.title = element_text(size = 10),
        aspect.ratio = 1.5,
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey95")) +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  labs(x = "",
       y = expression(F[ROH]),
       title = expression(paste("Distribution of ", F[ROH], " per Individual")))

print(p)

# ---- 5. Save ----
ggsave("FROH_distribution.png", plot = p, width = 5, height = 6, dpi = 300)
library(dplyr)
library(ggplot2)
library(scales)

# Merge average ROH with sex data

avg_roh_sex <- avg_roh %>%
  left_join(sex_data, by = c("IID" = "ID"))

ggplot(avg_roh_sex, aes(x = SEX, y = MeanROH_Kb, fill = SEX)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 2) +
  scale_fill_manual(values = c("Male" = "grey40", "Female" = "grey80")) +
  scale_y_continuous(labels = comma, limits = c(0, NA)) +
  theme_classic() +
  theme(plot.title = element_text(size = 12),
        aspect.ratio = 1.5,
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey95"),
        legend.position = "none") +
  labs(x = "Sex",
       y = "Mean ROH Length (Kb)",
       title = "Average ROH Length per Individual by Sex")
# Check distribution shape within each group
shapiro.test(avg_roh_sex$MeanROH_Kb[avg_roh_sex$SEX == "Male"])
shapiro.test(avg_roh_sex$MeanROH_Kb[avg_roh_sex$SEX == "Female"])
t.test(MeanROH_Kb ~ SEX, data = avg_roh_sex)



froh <- read.table("C:/Users/caitg/OneDrive/Documents/Work/Masters/Lion_Sequences/ROH/filteredvcf_results_froh.txt", header = TRUE)

ggplot(froh, aes(x = FROH)) +
  geom_histogram(binwidth = 0.02, colour = "white") +
  theme_classic() +
  labs(x = "Genomic Inbreeding Coefficient (FROH)",
       y = "Individual Lion Count",
       title = "Distribution of Genomic Inbreeding Coefficients")


