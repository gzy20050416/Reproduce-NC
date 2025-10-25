#!/usr/bin/env bash
set -euo pipefail
source ./00_config.env

if ! command -v cellranger >/dev/null 2>&1; then
  echo "未检测到 cellranger，自动安装..."
  cd "$SOFT"
  wget -c -O cellranger-9.0.1.tar.gz https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz
  tar -xzf cellranger-9.0.1.tar.gz
  ln -sf "$SOFT/cellranger-9.0.1/cellranger" /usr/local/bin/cellranger
  cd "$ROOT"
fi

mv $FASTQ/RNA/SRR13847285_1.fastq.gz $FASTQ/RNA/RNA_S1_L001_R1_001.fastq.gz
mv $FASTQ/RNA/SRR13847285_2.fastq.gz $FASTQ/RNA/RNA_S1_L001_R2_001.fastq.gz

mv $FASTQ/ADTSRR13847286_1.fastq.gz $FASTQ/ADTADT/ADT_S1_L001_R1_001.fastq.gz
mv $FASTQ/ADTSRR13847286_2.fastq.gz $FASTQ/ADTADT_S1_L001_R2_001.fastq.gz

mv $FASTQ/HTO/SRR13847287_1.fastq.gz $FASTQ/HTO/HTO_S1_L001_R1_001.fastq.gz
mv $FASTQ/HTO/SRR13847287_2.fastq.gz $FASTQ/HTO/HTO_S1_L001_R2_001.fastq.gz

echo "===== 运行 cellranger count ====="
OUT="$DATA/output"
mkdir -p "$OUT"

cellranger count \
  --id=Bcell_count \
  --libraries="$META/libraries.csv" \
  --feature-ref="$META/feature_reference.csv" \
  --transcriptome="$REF_PATH" \
  --localcores=$THREADS \
  --localmem=64

echo "Cellranger count 完成"
