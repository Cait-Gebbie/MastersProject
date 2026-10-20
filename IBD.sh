#!/usr/bin/env bash
set -euo pipefail
 
echo "🧬 Starting Pairwise IBD Analysis from VCF files"
 
##############################
# 1️⃣ Load environment
##############################
source ~/spack/share/spack/setup-env.sh
spack load bcftools
spack load plink@1.9
export LD_LIBRARY_PATH=/opt/spack/opt/spack/linux-ubuntu22.04-icelake/gcc-12.3.0/openblas-0.3.27-5vlpkj6rte7ndbth25yp7f45g6ih6hyz/lib:$LD_LIBRARY_PATH
 
##############################
# 2️⃣ Define directories
##############################
BASE_DIR="/home/caitling/Outputs"
BATCHES=("Batch1" "Batch2" "Batch3")
 
OUT_DIR="${BASE_DIR}/Output_notrealated/IBD"
mkdir -p "${OUT_DIR}"
 
VCF_LIST="${OUT_DIR}/IBD_vcfs.txt"
> "${VCF_LIST}"
 
##############################
# 3️⃣ Collect VCFs
##############################
echo "📂 Collecting VCF files..."
 
for BATCH in "${BATCHES[@]}"; do
    for SAMPLE_DIR in "${BASE_DIR}/${BATCH}"/*/ ; do
        SAMPLE=$(basename "${SAMPLE_DIR}")
        VCF="${SAMPLE_DIR}/${SAMPLE}_variantsGATK.vcf.gz"
 
        if [[ -f "$VCF" ]]; then
            if [[ -s "$VCF" ]]; then
                echo "$VCF" >> "${VCF_LIST}"
            else
                echo "⚠️ Empty VCF: ${BATCH}/${SAMPLE}"
            fi
        else
            echo "⚠️ Missing VCF: ${BATCH}/${SAMPLE}"
        fi
    done
done
 
TOTAL_VCFS=$(wc -l < "${VCF_LIST}")
echo "Total valid VCFs collected: $TOTAL_VCFS"
if [[ $TOTAL_VCFS -eq 0 ]]; then
    echo "❌ No valid VCFs found. Exiting."
    exit 1
fi
 
##############################
# 4️⃣ Merge VCFs
##############################
echo "🔗 Merging VCFs..."
MERGED_VCF="${OUT_DIR}/IBD_merged.vcf.gz"
 
bcftools merge \
    -l "${VCF_LIST}" \
    --missing-to-ref \
    -Oz \
    -o "$MERGED_VCF"
 
bcftools index "$MERGED_VCF"
 
##############################
# 5️⃣ Filter biallelic SNPs
##############################
echo "🧹 Filtering to biallelic SNPs..."
FILTERED_VCF="${OUT_DIR}/IBD_merged_biallelic.vcf.gz"
 
bcftools view \
    -m2 -M2 -v snps \
    "$MERGED_VCF" \
    -Oz -o "$FILTERED_VCF"
 
tabix -p vcf "$FILTERED_VCF"
 
##############################
# 6️⃣ Convert to PLINK format
##############################
echo "🔄 Converting VCF to PLINK format..."
RAW_PLINK="${OUT_DIR}/IBD_raw"
 
plink --vcf "$FILTERED_VCF" \
      --make-bed \
      --allow-extra-chr \
      --double-id \
      --out "$RAW_PLINK"
 
##############################
# 7️⃣ SNP filtering
##############################
echo "🧹 Filtering SNPs (MAF, missingness, HWE)..."
FILTERED_PLINK="${OUT_DIR}/IBD_filtered"
 
plink --bfile "$RAW_PLINK" \
      --maf 0.05 \
      --geno 0.05 \
      --mind 0.1 \
      --hwe 1e-5 \
      --make-bed \
      --allow-extra-chr \
      --out "$FILTERED_PLINK"
 
NSNP=$(wc -l < "${FILTERED_PLINK}.bim")
echo "SNPs after filtering: $NSNP"
if [[ $NSNP -lt 10000 ]]; then
    echo "⚠️ Warning: Very few SNPs ($NSNP). IBD estimates may be unreliable."
fi
 
##############################
# 8️⃣ Assign unique variant IDs
# FIX: Replace '.' IDs with CHR:POS:REF:ALT to avoid
#      "Duplicate ID '.'" error in --genome
##############################
echo "🔧 Assigning unique variant IDs..."
NAMED_PLINK="${OUT_DIR}/IBD_named"
 
plink --bfile "$FILTERED_PLINK" \
      --allow-extra-chr \
      --set-missing-var-ids @:#[b37]\$1\$2 \
      --make-bed \
      --out "$NAMED_PLINK"
 
# Also check for and remove any remaining duplicates after renaming
# (can occur if multiallelic sites weren't fully resolved)
echo "🔍 Checking for remaining duplicate IDs..."
DEDUP_PLINK="${OUT_DIR}/IBD_dedup"
 
cut -f2 "${NAMED_PLINK}.bim" | sort | uniq -d > "${OUT_DIR}/duplicate_ids.txt"
DUP_COUNT=$(wc -l < "${OUT_DIR}/duplicate_ids.txt")
 
if [[ $DUP_COUNT -gt 0 ]]; then
    echo "⚠️ Found $DUP_COUNT duplicate IDs after renaming — removing them..."
    plink --bfile "$NAMED_PLINK" \
          --allow-extra-chr \
          --exclude "${OUT_DIR}/duplicate_ids.txt" \
          --make-bed \
          --out "$DEDUP_PLINK"
else
    echo "✅ No duplicate IDs found."
    cp "${NAMED_PLINK}.bed" "${DEDUP_PLINK}.bed"
    cp "${NAMED_PLINK}.bim" "${DEDUP_PLINK}.bim"
    cp "${NAMED_PLINK}.fam" "${DEDUP_PLINK}.fam"
fi
 
##############################
# 9️⃣ LD pruning
##############################
echo "✂️ Performing LD pruning..."
PRUNE_BASE="${OUT_DIR}/IBD_prune"
 
plink --bfile "$DEDUP_PLINK" \
      --allow-extra-chr \
      --indep-pairwise 50 5 0.1 \
      --out "$PRUNE_BASE"
 
PRUNED_PLINK="${OUT_DIR}/IBD_pruned"
 
plink --bfile "$DEDUP_PLINK" \
      --extract "${PRUNE_BASE}.prune.in" \
      --make-bed \
      --allow-extra-chr \
      --out "$PRUNED_PLINK"
 
NSNP_PRUNED=$(wc -l < "${PRUNED_PLINK}.bim")
echo "SNPs after LD pruning: $NSNP_PRUNED"
 
##############################
# 🔟 Calculate allele frequencies
##############################
echo "📊 Calculating allele frequencies..."
plink --bfile "$PRUNED_PLINK" \
      --allow-extra-chr \
      --freq \
      --out "$PRUNED_PLINK"
 
##############################
# 1️⃣1️⃣ Pairwise IBD calculation
##############################
echo "🧬 Calculating pairwise IBD..."
plink --bfile "$PRUNED_PLINK" \
      --allow-extra-chr \
      --read-freq "${PRUNED_PLINK}.frq" \
      --genome \
      --out "${OUT_DIR}/IBD_pairs"
 
echo ""
echo "📊 Quick summary of PI_HAT distribution:"
awk 'NR>1 {print $10}' "${OUT_DIR}/IBD_pairs.genome" | sort -n | \
  awk 'BEGIN{n=0; s=0} {a[n++]=$1; s+=$1} END{
    print "N pairs: " n;
    print "Mean PI_HAT: " s/n;
    print "Min PI_HAT: " a[0];
    print "Max PI_HAT: " a[n-1];
  }'
 
##############################
# 1️⃣2️⃣ Identify & remove related samples
##############################
echo ""
echo "🔗 Identifying related samples (PI_HAT > 0.25)..."

N_RELATED=$(awk 'NR>1 && $10 > 0.25' "${OUT_DIR}/IBD_pairs.genome" | wc -l)
echo "Related pairs above threshold: $N_RELATED"

if [[ $N_RELATED -eq 0 ]]; then
    echo "✅ No related pairs found. No pruning needed."
    cp "${PRUNED_PLINK}.bed" "${OUT_DIR}/IBD_unrelated.bed"
    cp "${PRUNED_PLINK}.bim" "${OUT_DIR}/IBD_unrelated.bim"
    cp "${PRUNED_PLINK}.fam" "${OUT_DIR}/IBD_unrelated.fam"
else
    echo "🧹 Running greedy relatedness pruning..."
    # --rel-cutoff runs its own internal IBD calculation from the .bed file
    # It does NOT accept --read-genome; pass the allele freq file instead
    plink --bfile "$PRUNED_PLINK" \
          --allow-extra-chr \
          --read-freq "${PRUNED_PLINK}.frq" \
          --rel-cutoff 0.25 \
          --out "${OUT_DIR}/IBD_related_pruned"

    # Report removals
    N_BEFORE=$(wc -l < "${PRUNED_PLINK}.fam")
    N_KEEP=$(wc -l < "${OUT_DIR}/IBD_related_pruned.rel.id")
    N_REMOVED=$(( N_BEFORE - N_KEEP ))
    echo "Samples before : $N_BEFORE"
    echo "Samples kept   : $N_KEEP"
    echo "Samples removed: $N_REMOVED"

    # Save list of removed samples for QC records
    awk 'NR==FNR{keep[$1"_"$2]=1; next}
         !keep[$1"_"$2]{print $1, $2}' \
        "${OUT_DIR}/IBD_related_pruned.rel.id" \
        "${PRUNED_PLINK}.fam" \
        > "${OUT_DIR}/IBD_removed_samples.txt"

    ##############################
    # 1️⃣3️⃣ Apply filter
    ##############################
    echo "💾 Creating final unrelated dataset..."
    plink --bfile "$PRUNED_PLINK" \
          --allow-extra-chr \
          --keep "${OUT_DIR}/IBD_related_pruned.rel.id" \
          --make-bed \
          --out "${OUT_DIR}/IBD_unrelated"

    echo "✅ Done! Removed samples: ${OUT_DIR}/IBD_removed_samples.txt"
    echo "📄 Final dataset: ${OUT_DIR}/IBD_unrelated.bed/bim/fam"
fi
