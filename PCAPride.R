
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(readxl)
library(RColorBrewer)


PCA_DIR <- "C:/Users/caitg/OneDrive/Documents/Work/Masters/Results/OutputFiles/NoWhiteLion/"
OUT_DIR <- PCA_DIR

PRIDE_FILE <- paste0(PCA_DIR, "Lion - WGS Study_2026.08.21_Pride Estimates.xlsx")


eigenval <- read.table(paste0(PCA_DIR, "PCAgvcf_results.eigenval"), header = FALSE)
cat("eigenval rows:", nrow(eigenval), "\n")   # should be 89
colnames(eigenval) <- "Eigenvalue"
eigenval$PC <- paste0("PC", 1:nrow(eigenval))
eigenval <- eigenval %>%
  mutate(
    Var_explained  = Eigenvalue / sum(Eigenvalue) * 100,
    Cumulative_var = cumsum(Var_explained)
  )
cat("=== Variance Explained ===\n")
print(eigenval %>% select(PC, Eigenvalue, Var_explained, Cumulative_var) %>% head(10))

pc1_var <- round(eigenval$Var_explained[1], 2)
pc2_var <- round(eigenval$Var_explained[2], 2)
pc3_var <- round(eigenval$Var_explained[3], 2)


readLines(paste0(PCA_DIR, "PCAgvcf_results.eigenvec"), n = 3)
eigenvec <- read.table(paste0(PCA_DIR, "PCAgvcf_results.eigenvec"), header = FALSE, sep = "")
cat("eigenvec rows:", nrow(eigenvec), "\n")   # should also be 89
colnames(eigenvec) <- c("FID", "IID", paste0("PC", 1:20))
eigenvec$IID <- as.character(eigenvec$IID)


pride_data <- read_excel(PRIDE_FILE)

cat("Pride file columns:\n")
print(colnames(pride_data))

pride_data <- pride_data %>%
  rename(IID   = lion_id,
         Pride = pride_affiliation_estimated) %>%
  mutate(IID = as.character(IID))

cat("pride_data rows (before dedup):", nrow(pride_data), "\n")

# --- Duplicate check: list every duplicated ID, not just flag that one exists ---
dup_ids <- pride_data %>% count(IID) %>% filter(n > 1)
if (nrow(dup_ids) > 0) {
  cat(nrow(dup_ids), "duplicated lion_id(s) found in pride file:\n")
  print(dup_ids)
  cat("Rows for duplicated ID(s):\n")
  print(pride_data %>% filter(IID %in% dup_ids$IID))
} else {
  cat("No duplicate lion_ids found in pride file.\n")
}

pride_data <- pride_data %>%
  distinct(IID, .keep_all = TRUE)

cat("pride_data rows (after dedup):", nrow(pride_data), "\n")
cat("pride_data duplicated IIDs (should be 0):", sum(duplicated(pride_data$IID)), "\n")

pca_data <- eigenvec %>%
  left_join(pride_data %>% select(IID, Pride), by = "IID")

cat("pca_data rows after join (should match eigenvec rows):", nrow(pca_data), "\n")

n_unmatched <- sum(is.na(pca_data$Pride))
if (n_unmatched > 0) {
  cat(n_unmatched, "sample(s) had no matching pride record - labelled 'Unknown'.\n")
  cat("Unmatched IIDs:\n")
  print(pca_data$IID[is.na(pca_data$Pride)])
}
pca_data$Pride[is.na(pca_data$Pride)] <- "Unknown"

cat("Number of prides:", length(unique(pca_data$Pride)), "\n")


pride_order <- pca_data %>%
  count(Pride) %>%
  arrange(Pride) %>%
  pull(Pride)

pca_data <- pca_data %>%
  mutate(Pride = factor(Pride, levels = pride_order)) %>%
  group_by(Pride) %>%
  mutate(Pride_n = paste0(Pride, " (n = ", n(), ")")) %>%
  ungroup() %>%
  mutate(Pride_n = factor(Pride_n,
                          levels = unique(Pride_n[order(Pride)])))


pal <- readRDS(paste0(PCA_DIR, "pride_palette.rds"))

pride_raw      <- levels(pca_data$Pride)     # raw pride names, no "(n = X)" suffix
pride_n_labels <- levels(pca_data$Pride_n)

palette_colours <- setNames(pal$colours[pride_raw], pride_n_labels)
palette_shapes  <- setNames(pal$shapes[pride_raw],  pride_n_labels)


p <- ggplot(pca_data, aes(x = PC1, y = PC3,
                          colour = Pride_n,
                          shape  = Pride_n)) +
  geom_point(size = 3, stroke = 0.7) +
  scale_colour_manual(values = palette_colours, name = "KNP_pride") +
  scale_shape_manual(values = palette_shapes, name = "KNP_pride") +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),   # was pc2_var - now matches the PC1 axis
    y = paste0("PC3 (", pc3_var, "%)")    # was pc3_var - now matches the PC2 axis
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_blank(),
    axis.title = element_text(face = "bold")
  )

print(p)

