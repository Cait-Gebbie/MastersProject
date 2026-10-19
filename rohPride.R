library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(readxl)
library(RColorBrewer)

##############################
# File paths
##############################
ROH_DIR <- "C:/Users/caitg/OneDrive/Documents/Work/Masters/Results/OutputFiles/NoWhiteLion/"
OUT_DIR <- ROH_DIR

# PLINK --homozyg output (per-individual segment summary)
HOM_INDIV_FILE <- paste0(ROH_DIR, "filteredvcf_results.hom.indiv")   # <-- update filename if different

# Pride/pack affiliation file (same one used for the PCA plot)
PRIDE_FILE <- paste0(ROH_DIR, "Lion - WGS Study_2026.08.21_Pride Estimates.xlsx")

##############################
# 1. Load per-individual ROH summary
##############################
hom_indiv <- read.table(HOM_INDIV_FILE, header = TRUE, stringsAsFactors = FALSE)
hom_indiv$IID <- as.character(hom_indiv$IID)
hom_indiv$Total_length_Mb <- hom_indiv$KB / 1000   # convert Kb -> Mb, to match the other ROH plot

cat("=== ROH per-individual summary ===\n")
print(head(hom_indiv))
cat("hom_indiv rows:", nrow(hom_indiv), "\n")

##############################
# 2. Load pride/pack affiliation, check + remove duplicates, then merge
##############################
pride_data <- read_excel(PRIDE_FILE)
View(PRIDE_FILE)
cat("Pride file columns:\n")
print(colnames(pride_data))

pride_data <- pride_data %>%
  rename(IID  = lion_id,
         Pack = pride_affiliation_estimated) %>%
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

# --- Remove duplicates, keeping the first occurrence per IID ---
# Review the printed duplicate rows above first - if they disagree on Pack
# (not just incidental fields like date/GPS/age class), decide which record
# should win before trusting this default.
pride_data <- pride_data %>%
  distinct(IID, .keep_all = TRUE)

cat("pride_data rows (after dedup):", nrow(pride_data), "\n")
cat("pride_data duplicated IIDs (should be 0):", sum(duplicated(pride_data$IID)), "\n")

# --- Join and confirm row count didn't inflate ---
roh_data <- hom_indiv %>%
  left_join(pride_data %>% select(IID, Pack), by = "IID")

cat("roh_data rows after join (should match hom_indiv rows):", nrow(roh_data), "\n")

n_unmatched <- sum(is.na(roh_data$Pack))
if (n_unmatched > 0) {
  cat(n_unmatched, "sample(s) had no matching pack record - labelled 'Unknown'.\n")
  cat("Unmatched IIDs:\n")
  print(roh_data$IID[is.na(roh_data$Pack)])
}
roh_data$Pack[is.na(roh_data$Pack)] <- "Unknown"

cat("Number of prides:", length(unique(roh_data$Pack)), "\n")

##############################
# 3. Build legend labels with sample sizes, e.g. "A (n = 3)"
##############################
pack_order <- roh_data %>%
  count(Pack) %>%
  arrange(Pack) %>%
  pull(Pack)

roh_data <- roh_data %>%
  mutate(Pack = factor(Pack, levels = pack_order)) %>%
  group_by(Pack) %>%
  mutate(Pack_n = paste0(Pack, " (n = ", n(), ")")) %>%
  ungroup() %>%
  mutate(Pack_n = factor(Pack_n,
                         levels = unique(Pack_n[order(Pack)])))

##############################
# 4. Colour + shape palette - shared master palette
##############################
# Sourced from pride_palette.R so every plot (PCA, ROH, IBD, etc.) uses the
# EXACT same colour/shape per pride/pack. Run pride_palette.R once beforehand
# (or whenever the pride list changes) to (re)generate the .rds.
pal <- readRDS(paste0(ROH_DIR, "pride_palette.rds"))

pack_raw      <- levels(roh_data$Pack)     # raw pack names, no "(n = X)" suffix
pack_n_labels <- levels(roh_data$Pack_n)

palette_colours <- setNames(pal$colours[pack_raw], pack_n_labels)
palette_shapes  <- setNames(pal$shapes[pack_raw],  pack_n_labels)

##############################
# 5. Plot: number of ROH segments vs total ROH length, by pack
##############################
p <- ggplot(roh_data, aes(x = Total_length_Mb, y = NSEG,
                          colour = Pack_n,
                          shape  = Pack_n)) +
  geom_point(size = 3, stroke = 0.7) +
  scale_colour_manual(values = palette_colours, name = "KNP_pride") +
  scale_shape_manual(values = palette_shapes, name = "KNP_pride") +
  scale_x_continuous(labels = comma) +
  labs(
    x = "Sum total length of ROH (Mb)",
    y = "Number of ROH segments"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_blank(),
    axis.title = element_text(face = "bold")
  ) +
  guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))

print(p)

##############################
# 6. Save
##############################
ggsave(paste0(OUT_DIR, "ROH_segments_by_pack.png"), plot = p,
       width = 9.5, height = 6, dpi = 300)
##############################
# 6. Save
##############################
ggsave(paste0(OUT_DIR, "ROH_segments_by_pack.png"), plot = p,
       width = 9.5, height = 6, dpi = 300)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(readxl)

##############################
# File paths
##############################
ROH_DIR <- "C:/Users/caitg/OneDrive/Documents/Work/Masters/Results/OutputFiles/NoWhiteLion/"
OUT_DIR <- ROH_DIR

# PLINK --homozyg segment-level output (one row per ROH segment, not the
# per-individual .hom.indiv summary used in the other ROH script)
HOM_FILE <- paste0(ROH_DIR, "filteredvcf_results.hom")   # <-- update filename if different

# Pride/pack affiliation file (same one used for the PCA/ROH plots)
PRIDE_FILE <- paste0(ROH_DIR, "Lion - WGS Study_2026.08.21_Pride Estimates.xlsx")

# Shared colour/shape palette (run pride_palette.R first to generate this)
PALETTE_RDS <- paste0(ROH_DIR, "pride_palette.rds")

##############################
# 1. Load segment-level ROH calls
##############################
# plink --homozyg's .hom file has one row per ROH segment, with KB = that
# segment's length. This is different from .hom.indiv (already summed per
# individual) used in the other ROH script - here we need per-segment
# lengths so we can bin them into length categories.
hom <- read.table(HOM_FILE, header = TRUE, stringsAsFactors = FALSE)
hom$IID <- as.character(hom$IID)

cat("=== ROH segment-level data ===\n")
print(head(hom))

##############################
# 2. Bin segments into length categories (Mb)
##############################
hom <- hom %>%
  mutate(
    Length_Mb = KB / 1000,
    Length_category = cut(
      Length_Mb,
      breaks = c(0.5, 1.0, 1.5, 2.0, 2.5, 3, 3.5, 4, 4.5, 5,Inf),
      labels = c("(0.5-1.0]", "(1.0-1.5]", "(1.5-2.0]", "(2.0-2.5]", "(2.5-3]", "(3-3.5]", "(3.5-4]", "(4-4.5]", "(4.5-5]", ">5"),
      right  = TRUE
    )
  ) %>%
  filter(!is.na(Length_category))   # drops segments below 0.5 Mb, if any

##############################
# 3. Sum total ROH length per individual, per length category
##############################
roh_by_category <- hom %>%
  group_by(IID, Length_category) %>%
  summarise(Total_length_Mb = sum(Length_Mb), .groups = "drop")

##############################
# 4. Load pack affiliation and merge
##############################
pride_data <- read_excel(PRIDE_FILE) %>%
  rename(IID  = lion_id,
         Pack = pride_affiliation_estimated) %>%
  mutate(IID = as.character(IID))

roh_by_category <- roh_by_category %>%
  left_join(pride_data %>% select(IID, Pack), by = "IID")


n_unmatched <- sum(is.na(roh_by_category$Pack))
if (n_unmatched > 0) {
  cat(n_unmatched, "row(s) had no matching pack record - labelled 'Unknown'.\n")
}
roh_by_category$Pack[is.na(roh_by_category$Pack)] <- "Unknown"

##############################
# 5. Apply shared colour + shape palette
##############################
pal <- readRDS(PALETTE_RDS)

roh_by_category <- roh_by_category %>%
  mutate(Pack = factor(Pack, levels = pal$levels))

palette_colours <- pal$colours
palette_shapes  <- pal$shapes

# Build "A (n = 3)" style legend labels from distinct individuals per pack
# (not rows in roh_by_category, since each individual has several rows -
# one per length category)
pack_n <- roh_by_category %>%
  distinct(IID, Pack) %>%
  count(Pack, name = "n")

pack_labels <- setNames(
  paste0(pack_n$Pack, " (n = ", pack_n$n, ")"),
  pack_n$Pack
)
# Preserve factor order/levels so labels line up with palette_colours/shapes
pack_labels <- pack_labels[levels(roh_by_category$Pack)]

##############################
# 6. Plot: total ROH length per individual across length categories
##############################

p <- ggplot(roh_by_category, aes(x = Length_category, y = Total_length_Mb,
                                 group = IID,
                                 colour = Pack,
                                 shape  = Pack)) +
  geom_line(linewidth = 0.4, alpha = 0.7) +
  geom_point(size = 2.5, alpha = 0.9, stroke = 0.6) +
  scale_colour_manual(values = palette_colours,labels=pack_labels, name = "KNP_pride") +
  scale_shape_manual(values = palette_shapes,labels=pack_labels, name = "KNP_pride") +
  scale_y_continuous(labels = comma) +
  labs(
    x = "ROH length category (Mb)",
    y = "Total length of ROH per individual (Mb)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey70"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_blank(),
    axis.title = element_text(face = "bold")
  ) +
  guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))

print(p)

##############################
# 7. Save
##############################
ggsave(paste0(OUT_DIR, "ROH_length_by_category_pack.png"), plot = p,
       width = 9.5, height = 6, dpi = 300)


