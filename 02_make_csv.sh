#!/usr/bin/env bash
set -euo pipefail
source ./00_config.env

mkdir -p "$META"

# === 生成 libraries.csv ===
LIB_CSV="$META/libraries.csv"
cat > "$LIB_CSV" <<'EOF'
library_id,fastqs,feature_types
RNA,$FASTQ/RNA,Gene Expression
ADT,$FASTQ/ADT,Antibody Capture
HTO,$FASTQ/HTO,Antibody Capture
EOF
echo "写入: $LIB_CSV"

# === 生成 feature_reference.csv ===
FEAT_CSV="$META/feature_reference.csv"
cat > "$FEAT_CSV" <<'EOF'
id,name,feature_type,read,pattern,sequence
HTO1,Hashtag1,Antibody Capture,R2,CCCCCCCCNNNNNNNN(BC),ACCCACCAGTAAGAC
HTO2,Hashtag2,Antibody Capture,R2,CCCCCCCCNNNNNNNN(BC),GGTCGAGAGCATTCA
ADT1,B220,Antibody Capture,R2,CCCCCCCCNNNNNNNN(BC),CCTACACCTCATAAT
ADT2,CD19,Antibody Capture,R2,CCCCCCCCNNNNNNNN(BC),ATCAGCCATGTCAGT
ADT3,CD93,Antibody Capture,R2,CCCCCCCCNNNNNNNN(BC),GGTATTTCCTGTGGT
ADT4,CD25,Antibody Capture,R2,CCCCCCCCNNNNNNNN(BC),ACCATGAGACACAGT
ADT5,IgM,Antibody Capture,R2,CCCCCCCCNNNNNNNN(BC),AGCTACGCATTCAAT
ADT6,Streptavidin,Antibody Capture,R2,CCCCCCCCNNNNNNNN(BC),AACCTTTGCCACTGC
EOF
echo "写入: $FEAT_CSV"

echo "CSV 文件生成完成"
