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
│   ├── refdata/
├── Ranalyze/
└── software/
    └── cellranger-9.0.1/
```

如不支持进行上游复现，control_a是对作者发布在geo数据库上数据进行处理，也可运行

在终端界面按以下运行代码即可
```text
bash 01_download.sh
bash 02_make_csv.sh
bash 03_run_cellranger.sh
Rscript Control_a.R
#或者运行Rscript Control.R
Rscript clear.R
```
