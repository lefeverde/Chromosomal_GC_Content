# Introduction

This is the script used in performing the GC content of the human genome in [CITATION]. The output of the script is two folders "unique_csvs" and "raw_csvs", and a file "concat_file.csv". The files in the raw_csvs folder are the raw results for every CDS feature listed in the assembeld genome. The files in the unique_csvs folder are the CDS features which contain the CDS entries with unique flanking regions, and the concat_file is just a concatenation of these.

The data from the concat_file was used to generate the histogram shown in [FIGURE]. It was created using the ggplot2 package in R.

# Dependencies

- Python 2.7
- Biopython

# Usage

```
python chromo_gc_content.py  -i GCF_000001405.34_GRCh38.p8_genomic.gbff
```

# Input File Information

The input file was download from the NCBI ftp server. The file can be obtained by running the following comand:
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.34_GRCh38.p8/GCF_000001405.34_GRCh38.p8_genomic.gbff.gz
```

