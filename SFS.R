
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)

# Read vcftools file
freq <- read.table(
  "C:/Users/caitg/OneDrive/Documents/Work/Masters/Results/OutputFiles/NoWhiteLion/SFS.frq",
  header = FALSE,
  skip = 1,
  stringsAsFactors = FALSE
)

colnames(freq) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "ALLELE1", "ALLELE2")

# Total chromosomes
n_chr <- freq$N_CHR

# Extract allele frequencies
allele1 <- as.numeric(sub(".*:", "", freq$ALLELE1))
allele2 <- as.numeric(sub(".*:", "", freq$ALLELE2))

# Minor allele frequency
maf <- pmin(allele1, allele2)

# Convert to minor allele count
mac <- round(maf * n_chr)

# Remove monomorphic
mac <- mac[mac > 0]

# Create SFS
sfs <- as.data.frame(table(mac))
colnames(sfs) <- c("MAC", "Count")
sfs$MAC <- as.numeric(as.character(sfs$MAC))
library(scales)

ggplot(sfs, aes(x = MAC, y = Count)) +
  geom_bar(stat = "identity",
           fill = "#0A9396",
           colour = "white") +
  scale_y_continuous(labels = label_comma()) +
  theme_classic() +
  labs(
    x = "Minor Allele Count",
    y = "Number of SNPs"
  )
View(sfs)
