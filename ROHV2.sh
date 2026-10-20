#!/usr/bin/env bash
set -euo pipefail

echo "Starting corrected multi-batch ROH pipeline (v2 - reusing depth-refiltered IBD VCF)"

##############################
# 1. Load environment
##############################
source ~/spack/share/spack/setup-env.sh
spack load bcftools
spack load plink@1.9
export LD_LIBRARY_PATH=/opt/spack/opt/spack/linux-ubuntu22.04-icelake/gcc-12.3.0/openblas-0.3.27-5vlpkj6rte7ndbth25yp7f45g6ih6hyz/lib:$LD_LIBRARY_PATH

##############################
# 2. Directories
##############################
BASE_DIR="/home/caitling/Outputs"
IBD_OUT_DIR="${BASE_DIR}/IBD"
SHARED_VCF="${IBD_OUT_DIR}/IBD_merged_biallelic_v2.vcf.gz"

OUT_DIR="${BASE_DIR}/Output_notrealated/ROH/related"
mkdir -p "${OUT_DIR}"

if [[ ! -s "$SHARED_VCF" ]]; then
    echo "ERROR: v2 filtered VCF from the IBD pipeline not found at ${SHARED_VCF}"
    echo "Run rerun_depth_filtering.sh first, or update SHARED_VCF to the correct path."
    exit 1
fi

##############################
# 3. Convert to PLINK
##############################
echo "Converting to PLINK..."

RAW_PLINK="${OUT_DIR}/filteredvcf_raw"

plink --vcf "$SHARED_VCF" \
      --allow-extra-chr \
      --double-id \
      --make-bed \
      --out "$RAW_PLINK"

##############################
echo "Removing samples ..."
RELATED_SAMPLES="${OUT_DIR}/related_samples_to_remove.txt"

cat > "${RELATED_SAMPLES}" << 'SAMPLES'
C165 C165
C152 C152
C90 C90
C143 C143
C96 C96
C125 C125
C168 C168
SAMPLES

UNRELATED_PLINK="${OUT_DIR}/filteredvcf_unrelated"

plink --bfile "$RAW_PLINK" \
      --remove "${RELATED_SAMPLES}" \
      --make-bed \
      --allow-extra-chr \
      --out "$UNRELATED_PLINK"

N_KEPT=$(wc -l < "${UNRELATED_PLINK}.fam")
echo "Samples retained after relatedness removal: ${N_KEPT}"

##############################
# 4. Assign unique variant IDs
##############################
echo "Assigning unique variant IDs..."

cp "${UNRELATED_PLINK}.bed" "${OUT_DIR}/filteredvcf_named.bed"
cp "${UNRELATED_PLINK}.fam" "${OUT_DIR}/filteredvcf_named.fam"

awk '{print $1, $1"_"$4"_"$5"_"$6, $3, $4, $5, $6}' \
    "${UNRELATED_PLINK}.bim" > "${OUT_DIR}/filteredvcf_named.bim"

NAMED_PLINK="${OUT_DIR}/filteredvcf_named"
DEDUP_PLINK="${OUT_DIR}/filteredvcf_dedup"

echo "Checking for duplicate IDs..."

cut -f2 "${NAMED_PLINK}.bim" | sort | uniq -d > "${OUT_DIR}/filteredvcfduplicate_ids.txt"
DUP_COUNT=$(wc -l < "${OUT_DIR}/filteredvcfduplicate_ids.txt")

if [[ $DUP_COUNT -gt 0 ]]; then
    echo "Removing $DUP_COUNT duplicate variants..."

    awk 'FNR==NR {dup[$1]++; next}
         dup[$2] > 1 {dup[$2]--; next}
         {print $2}' \
        "${OUT_DIR}/filteredvcfduplicate_ids.txt" "${NAMED_PLINK}.bim" \
        > "${OUT_DIR}/keep_ids.txt"

    plink --bfile "$NAMED_PLINK" \
          --allow-extra-chr \
          --extract "${OUT_DIR}/keep_ids.txt" \
          --make-bed \
          --out "$DEDUP_PLINK"
else
    echo "No duplicates found."
    cp "${NAMED_PLINK}.bed" "${DEDUP_PLINK}.bed"
    cp "${NAMED_PLINK}.bim" "${DEDUP_PLINK}.bim"
    cp "${NAMED_PLINK}.fam" "${DEDUP_PLINK}.fam"
fi

##############################
# 5. PLINK-level filtering
##############################
echo "Applying PLINK-level filtering..."

FILTERED_PLINK="${OUT_DIR}/filteredvcf_filtered"

plink \
    --bfile "$DEDUP_PLINK" \
    --allow-extra-chr \
    --maf 0.05 \
    --mind 0.1 \
    --geno 0.05\
    --make-bed \
    --out "$FILTERED_PLINK"

##############################
# 6. Missingness check
##############################
echo "Checking missingness..."

plink --bfile "$FILTERED_PLINK" \
      --allow-extra-chr \
      --missing \
      --out "${OUT_DIR}/filteredvcf_missingness"

##############################
# 7. Run ROH detection
##############################
echo "Running ROH detection..."

plink --bfile "$FILTERED_PLINK" \
      --allow-extra-chr \
      --homozyg \
      --homozyg-window-snp 25 \
      --homozyg-snp 25 \
      --homozyg-kb 500 \
      --homozyg-gap 1000 \
      --homozyg-density 5000 \
      --homozyg-window-missing 5 \
      --homozyg-window-het 1 \
      --out "${OUT_DIR}/filteredvcf_results"

##############################
# Summary
##############################
echo ""
echo "ROH analysis complete!"
echo ""
echo "Quick ROH summary:"

awk '
NR>1 {
    nroh[$2]++
    total[$2]+=$9
}
END{
    printf "%-8s %-8s %-15s\n","Sample","N_ROH","Total_ROH_Mb"
    for(i in total)
        printf "%-8s %-8d %-15.2f\n", i, nroh[i], total[i]/1000
}' "${OUT_DIR}/filteredvcf_results.hom" | sort
