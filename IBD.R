library(ggplot2)
library(dplyr)
library(ltc)



pal <- ltc("heatmap0")   
plot(pal)                

ibd <- read.table("C:/Users/caitg/OneDrive/Documents/Work/Masters/Results/OutputFiles/NoWhiteLion/IBD_pairs_v2.genome", header = TRUE)


ggplot(ibd, aes(x = PI_HAT)) +
  geom_histogram(binwidth = 0.025,
                 fill = pal[2],
                 color = "white") +
  geom_vline(xintercept = c(0.0625, 0.125, 0.25, 0.5),
             linetype = "dashed",
             color = pal[1],
             linewidth = 0.7) +
  annotate("text", x = c(0.0625, 0.125, 0.25, 0.5),
           y = Inf, vjust = 2,
           label = c("3rd degree", "2nd degree", "1st degree", "Full siblings/Parent-offspring"),
           color = pal[1],
           size = 3, angle = 90, hjust = 1) +
  labs(title = "Distribution of Pairwise IBD (PI_HAT)",
       x = "PI_HAT",
       y = "Number of pairs") +
  theme_classic()



ggplot(ibd, aes(x = Z0, y = Z1, color = PI_HAT)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_colour_gradientn(colours = pal, name = "PI_HAT") +
  labs(title = "IBD Coefficients: Z0 vs Z1",
       x = "Z0 (P(IBD=0))", y = "Z1 (P(IBD=1))") +
  theme_classic()



ggplot(ibd, aes(x = IID1, y = IID2, fill = PI_HAT)) +
  geom_tile() +
  scale_fill_gradientn(colours = pal, name = "PI_HAT") +
  labs(title = "Pairwise IBD Heatmap", x = "", y = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6))



ibd <- ibd %>%
  mutate(Relationship = case_when(
    PI_HAT >= 0.35 & PI_HAT < 0.85 ~ "1st degree",
    PI_HAT >= 0.1875 & PI_HAT < 0.35 ~ "2nd degree",
    PI_HAT >= 0.09375 & PI_HAT < 0.1875 ~ "3rd degree",
    TRUE ~ "Unrelated"
  ))

relationship_summary <- ibd %>%
  group_by(Relationship) %>%
  summarise(
    n_pairs     = n(),
    mean_PI_HAT = mean(PI_HAT),
    sd_PI_HAT   = sd(PI_HAT),
    min_PI_HAT  = min(PI_HAT),
    max_PI_HAT  = max(PI_HAT)
  ) %>%
  arrange(desc(mean_PI_HAT))

print(relationship_summary)

# Bar chart of relationship classes
ggplot(ibd, aes(x = Relationship, fill = Relationship)) +
  geom_bar() +
  scale_fill_manual(values = pal) +
  labs(title = "Inferred Relationship Classes",
       x = "Relationship", y = "Number of pairs") +
  theme_classic() +
  theme(legend.position = "none")

ggplot(ibd, aes(x = Relationship, fill = Relationship)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.5) +
  scale_fill_manual(values = pal) +
  labs(title = "Inferred Relationship Classes",
       x = "Relationship", y = "Number of pairs") +
  theme_classic() +
  theme(legend.position = "none")

ibd %>%
  count(Relationship) %>%
  arrange(desc(n)) %>%
  print()

