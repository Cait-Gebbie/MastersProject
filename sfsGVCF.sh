#!/usr/bin/env bash
set -euo pipefail
echo "Starting Site Frequency Spectrum analysis (v2 - using depth-refiltered IBD VCF)"

##############################
# 1. Load environment
##############################
source ~/spack/share/spack/setup-env.sh
spack load bcftools
spack load vcftools

##############################
# 2. Directories / input
# SHARED_VCF is the v2 output of the IBD pipeline (rerun_depth_filtering.sh) -
# joint genotyped, GATK hard-filtered, depth-filtered (without the flat
# MAX_GENO_DP cap), and already restricted to biallelic SNPs.
##############################
BASE_DIR="/home/caitling/Outputs"
SHARED_VCF="${BASE_DIR}/IBD/IBD_merged_biallelic_v2.vcf.gz"

OUTDIR="${BASE_DIR}/Output_norelated/SFS/related"
mkdir -p "${OUTDIR}"

if [[ ! -s "$SHARED_VCF" ]]; then
    echo "ERROR: v2 filtered VCF from the IBD pipeline not found at ${SHARED_VCF}"
    echo "Run rerun_depth_filtering.sh first, or update SHARED_VCF to the correct path."
    exit 1
fi

##############################
# 3. Remove related samples (PI_HAT > 0.25)
# Updated list (44 samples) derived via minimum vertex cover of the
# v2 IBD relatedness graph (post-GVCF, depth-filter-corrected run) -
# removing this set resolves all related pairs while keeping the
# maximum number of samples (47 retained from 91). Related
# individuals share alleles by descent rather than independently,
# which biases the SFS toward the rare-variant classes if left in.
##############################
echo "Removing 44 related samples (PI_HAT > 0.25)..."
EXCLUDE_SAMPLES="${OUTDIR}/related_samples_to_exclude.txt"

cat > "${EXCLUDE_SAMPLES}" << 'SAMPLES'
C165
C152
C90
C143
C96
C125
C168
SAMPLES

bcftools view \
    -S ^"${EXCLUDE_SAMPLES}" \
    "$SHARED_VCF" \
    -Oz -o "${OUTDIR}/SFS_unrelated.vcf.gz"
bcftools index "${OUTDIR}/SFS_unrelated.vcf.gz"

N_KEPT=$(bcftools query -l "${OUTDIR}/SFS_unrelated.vcf.gz" | wc -l)
echo "Samples retained after relatedness removal: ${N_KEPT}"

##############################
echo "Applying site-level missingness filter (<10%)..."
bcftools view \
    -i 'F_MISSING<0.1' \
    "${OUTDIR}/SFS_unrelated.vcf.gz" \
    -Oz -o "${OUTDIR}/SFS_filtered.vcf.gz"
bcftools index "${OUTDIR}/SFS_filtered.vcf.gz"

##############################
# 5. Calculate allele frequencies
##############################
echo "Calculating allele frequencies..."
vcftools \
  --gzvcf "${OUTDIR}/SFS_filtered.vcf.gz" \
  --freq \
  --out "${OUTDIR}/SFS"

echo "SFS frequency file created:"
echo "${OUTDIR}/SFS.frq"
echo "SFS processing complete!"
