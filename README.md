主要工作目录如下，已在env文件写好相对路径
注意下载的fastq.gz文件的完整性，与ENA网站的md5校验之进行比对是否一致
```text
md5sum "yourdata"
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
