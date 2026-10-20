#!/usr/bin/env bash
set -euo pipefail

echo "Starting PCA pipeline (v2 - using depth-refiltered IBD VCF)"

############################################################
# Load software
############################################################

source ~/spack/share/spack/setup-env.sh

spack load plink@1.9
spack load plink2
spack load bcftools

export LD_LIBRARY_PATH=/opt/spack/opt/spack/linux-ubuntu22.04-icelake/gcc-12.3.0/openblas-0.3.27-5vlpkj6rte7ndbth25yp7f45g6ih6hyz/lib:$LD_LIBRARY_PATH

############################################################
# Directories / input
############################################################

BASE_DIR="/home/caitling/Outputs"
OUT_DIR="${BASE_DIR}/Output_notrealated/PCA/related"
mkdir -p "$OUT_DIR"

# Built directly from the v2 depth-refiltered, biallelic-SNP VCF from
# the IBD pipeline.
SHARED_VCF="${BASE_DIR}/IBD/IBD_merged_biallelic_v2.vcf.gz"

if [[ ! -s "$SHARED_VCF" ]]; then
    echo "ERROR: v2 filtered VCF from the IBD pipeline not found at ${SHARED_VCF}"
    echo "Run rerun_depth_filtering.sh first, or update SHARED_VCF to the correct path."
    exit 1
fi

############################################################
# Convert shared VCF to PLINK
############################################################

echo ""
echo "Converting VCF to PLINK format..."

RAW_PLINK="${OUT_DIR}/PCAgvcf_raw"

plink \
    --vcf "$SHARED_VCF" \
    --allow-extra-chr \
    --double-id \
    --make-bed \
    --out "$RAW_PLINK"

############################################################

############################################################

echo ""
echo "Removing related individuals..."

RELATED_SAMPLES="${OUT_DIR}/PCAgvcf_related_samples_to_remove.txt"

cat > "${RELATED_SAMPLES}" << 'SAMPLES'
C165 C165
C152 C152
C90 C90
C143 C143
C96 C96
C125 C125
C168 C168
SAMPLES

UNRELATED="${OUT_DIR}/PCAgvcf_unrelated"

plink \
    --bfile "$RAW_PLINK" \
    --remove "$RELATED_SAMPLES" \
    --allow-extra-chr \
    --make-bed \
    --out "$UNRELATED"

N_KEPT=$(wc -l < "${UNRELATED}.fam")
echo "Samples retained after relatedness removal: ${N_KEPT}"

############################################################
# Quality filtering
############################################################

echo ""
echo "Applying quality filters..."

FILTERED="${OUT_DIR}/PCAgvcf_filtered"

plink \
    --bfile "$UNRELATED" \
    --allow-extra-chr \
    --maf 0.05 \
    --geno 0.05 \
    --mind 0.10 \
    --make-bed \
    --out "$FILTERED"

############################################################
# Create unique SNP IDs
############################################################

echo ""
echo "Creating unique variant IDs..."

UNIQUE="${OUT_DIR}/PCAgvcf_unique"

plink2 \
    --bfile "$FILTERED" \
    --allow-extra-chr \
    --set-all-var-ids @:#:\$r:\$a \
    --make-bed \
    --out "$UNIQUE"

############################################################
# LD pruning
############################################################

echo ""
echo "Performing LD pruning..."

PRUNE="${OUT_DIR}/PCAgvcf_pruned"

plink \
    --bfile "$UNIQUE" \
    --allow-extra-chr \
    --indep-pairwise 50 10 0.1 \
    --out "$PRUNE"

############################################################
# Extract independent SNPs
############################################################

FINAL="${OUT_DIR}/PCAgvcf_final"

plink \
    --bfile "$UNIQUE" \
    --allow-extra-chr \
    --extract "${PRUNE}.prune.in" \
    --make-bed \
    --out "$FINAL"

############################################################
# Summary after filtering
############################################################

NSAMP=$(wc -l < "${FINAL}.fam")
NSNP=$(wc -l < "${FINAL}.bim")

echo ""
echo "Samples retained : ${NSAMP}"
echo "Independent SNPs : ${NSNP}"

############################################################
# Run PCA
############################################################

echo ""
echo "Running PCA..."

plink2 \
    --bfile "$FINAL" \
    --allow-extra-chr \
    --pca 20 \
    --out "${OUT_DIR}/PCAgvcf_results"

############################################################
# Missingness report
############################################################

plink \
    --bfile "$FINAL" \
    --allow-extra-chr \
    --missing \
    --out "${OUT_DIR}/PCAgvcf_missing"

############################################################
# Variance explained
############################################################

echo ""
echo "Variance explained"

awk '
{
    total += $1
    eig[NR]=$1
}
END{
    for(i=1;i<=NR;i++)
        printf("PC%-2d  Eigenvalue = %.4f   %.2f%%\n",
               i,eig[i],100*eig[i]/total)
}' "${OUT_DIR}/PCAgvcf_results.eigenval"

############################################################
# Finished
############################################################

echo ""
echo "========================================="
echo "PCA COMPLETE"
echo "========================================="

echo "Samples : $(wc -l < ${FINAL}.fam)"
echo "SNPs    : $(wc -l < ${FINAL}.bim)"

echo ""
echo "Outputs"

echo "Eigenvectors : ${OUT_DIR}/PCAgvcf_results.eigenvec"
echo "Eigenvalues  : ${OUT_DIR}/PCAgvcf_results.eigenval"
echo "Missingness  : ${OUT_DIR}/PCAgvcf_missing.imiss"
