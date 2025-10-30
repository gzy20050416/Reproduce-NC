文献链接

https://www.nature.com/articles/s41467-021-27232-5

注意下载的fastq.gz文件的完整性，与ENA网站的md5校验之进行比对是否一致
```text
md5sum "yourdata"
```
cellranger运行时需要规范命名，已进行改名

主要工作目录如下，已在env文件写好相对路径
```text
NC/
├── Bcell/
│   ├── fastq/
│   │   ├── ADT/
│   │   ├── HTO/
│   │   └── RNA/
│   ├── meta/
│   └── refdata/
├── Bcell_count
├── Ranalyze/
│   └── raw/
└── software/
    └── cellranger-9.0.1/
```

如不支持进行上游复现，control_a是对作者发布在geo数据库上数据进行处理，也可运行
访问https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168158，将下载作者的以下文件
```text
GSM5130034_WT_1_singlecell_gex_raw_counts.txt.gz	
GSM5130034_WT_2_singlecell_gex_raw_counts.txt.gz	
GSM5130035_WT_1_singlecell_adt_raw_counts.txt.gz	
GSM5130035_WT_2_singlecell_adt_raw_counts.txt.gz	
GSM5130036_WT_1_singlecell_hto_raw_counts.txt.gz	
GSM5130036_WT_2_singlecell_hto_raw_counts.txt.gz	
```
将以上结果保存到Ranalyze/raw

在终端界面按以下运行代码即可
```text
bash 01_download.sh
bash 02_make_csv.sh
bash 03_run_cellranger.sh
Rscript Control_a.R
#或者运行下面这个若无法进行上面的分析
#Rscript Control.R
Rscript clear.R
```
