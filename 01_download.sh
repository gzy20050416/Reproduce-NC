#!/usr/bin/env bash
set -euo pipefail
source ./00_config.env

command -v fastp >/dev/null 2>&1 || { echo "fastp 未安装"; exit 1; }

echo "fastp: $(fastp --version)"

### === 数据下载 ===
echo "===== 下载 FASTQ 文件 ====="

cd "$FASTQ/RNA"
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/085/SRR13847285/SRR13847285_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/085/SRR13847285/SRR13847285_2.fastq.gz

cd "$FASTQ/ADT"
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/086/SRR13847286/SRR13847286_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/086/SRR13847286/SRR13847286_2.fastq.gz

cd "$FASTQ/HTO"
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/087/SRR13847287/SRR13847287_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/087/SRR13847287/SRR13847287_2.fastq.gz

echo "FASTQ 数据下载完成"

### === 参考基因组下载 ===
echo "===== [4/6] 下载并解压官方参考基因组 ====="
cd "$REF"
wget -c http://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz

REF_PATH="$REF/refdata-gex-mm10-2020-A"

echo "参考基因组路径: $REF_PATH"
