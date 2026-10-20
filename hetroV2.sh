#!/usr/bin/env bash
set -euo pipefail
echo "Starting multi-batch heterozygosity analysis (v2 - using depth-refiltered IBD VCF)"

##############################
# 1. Load environment
##############################
source ~/spack/share/spack/setup-env.sh
spack load bcftools
spack load plink@1.9
export LD_LIBRARY_PATH=/opt/spack/opt/spack/linux-ubuntu22.04-icelake/gcc-12.3.0/openblas-0.3.27-5vlpkj6rte7ndbth25yp7f45g6ih6hyz/lib:$LD_LIBRARY_PATH

##############################
# 2. Directories
# SHARED_VCF is the v2 output of the IBD pipeline (rerun_depth_filtering.sh) -
# joint genotyped, GATK hard-filtered, depth-filtered (without the flat
# MAX_GENO_DP cap), and already restricted to biallelic SNPs.
##############################
BASE_DIR="/home/caitling/Outputs"
IBD_OUT_DIR="${BASE_DIR}/IBD"
SHARED_VCF="${IBD_OUT_DIR}/IBD_merged_biallelic_v2.vcf.gz"

OUT_DIR="${BASE_DIR}/Output_notrealated/Hetro/related"
mkdir -p "${OUT_DIR}"

if [[ ! -s "$SHARED_VCF" ]]; then
    echo "ERROR: v2 filtered VCF from the IBD pipeline not found at ${SHARED_VCF}"
    echo "Run rerun_depth_filtering.sh first, or update SHARED_VCF to the correct path."
    exit 1
fi

##############################
# 4. Convert QC'd VCF to PLINK (FULL DATASET)
##############################
echo "Converting full dataset to PLINK..."
plink --vcf "$SHARED_VCF" \
      --allow-extra-chr \
      --double-id \
      --make-bed \
      --out "${OUT_DIR}/hetGVCF_full"

##############################
# 4b. Remove related samples (PI_HAT > 0.25)
# Updated list (44 samples) derived via minimum vertex cover of the
# v2 IBD relatedness graph.
##############################
echo "Removing 44 related samples (PI_HAT > 0.25)..."
RELATED_SAMPLES="${OUT_DIR}/hetGVCF_related_samples_to_remove.txt"

cat > "${RELATED_SAMPLES}" << 'SAMPLES'
C165 C165
C152 C152
C90 C90
C143 C143
C96 C96
C125 C125
C168 C168
SAMPLES

plink --bfile "${OUT_DIR}/hetGVCF_full" \
      --remove "${RELATED_SAMPLES}" \
      --make-bed \
      --allow-extra-chr \
      --out "${OUT_DIR}/hetGVCF_unrelated"

N_KEPT=$(wc -l < "${OUT_DIR}/hetGVCF_unrelated.fam")
echo "Samples retained after relatedness removal: ${N_KEPT}"

##############################
# 5. Heterozygosity - FULL unrelated dataset (QC'd, biallelic)
##############################
echo "Calculating heterozygosity (full unrelated dataset)..."
plink --bfile "${OUT_DIR}/hetGVCF_unrelated" \
      --allow-extra-chr \
      --het \
      --out "${OUT_DIR}/hetGVCF_unrelated"

##############################
# 6. Apply MAF + HWE filtering
##############################
echo "Applying MAF + HWE filtering..."
plink --bfile "${OUT_DIR}/hetGVCF_unrelated" \
      --allow-extra-chr \
      --maf 0.05 \
      --hwe 1e-5 \
      --mind 0.1 \
      --geno 0.05\
      --make-bed \
      --out "${OUT_DIR}/hetGVCF_filtered"

##############################
# 7. LD pruning
##############################
echo "Performing LD pruning..."
plink --bfile "${OUT_DIR}/hetGVCF_filtered" \
      --allow-extra-chr \
      --indep-pairwise 200 20 0.2 \
      --out "${OUT_DIR}/hetGVCF_prune"

plink --bfile "${OUT_DIR}/hetGVCF_filtered" \
      --allow-extra-chr \
      --extract "${OUT_DIR}/hetGVCF_prune.prune.in" \
      --make-bed \
      --out "${OUT_DIR}/hetGVCF_filtered_LD"

##############################
# 8. Heterozygosity - Filtered + LD pruned dataset (QC'd, biallelic)
##############################
echo "Calculating heterozygosity (filtered + LD pruned)..."
plink --bfile "${OUT_DIR}/hetGVCF_filtered_LD" \
      --allow-extra-chr \
      --het \
      --out "${OUT_DIR}/hetGVCF_filtered_LD"

echo "Heterozygosity analysis complete!"
echo "Full unrelated het results:      ${OUT_DIR}/hetGVCF_unrelated.het"
echo "Filtered + LD-pruned het results: ${OUT_DIR}/hetGVCF_filtered_LD.het"
