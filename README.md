# Create Admixture model files from WBBC project of Westlake University
Create admixture model files (*.F, *.alleles) based on WBBC project of Westlake University. Microarray chip files (TSV format) from popular genetic testing companies are used as referenced SNPs. Linkage disequilibrium SNPs could be eliminated optionally (https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)). The output admixture model files could be used for model based admixture calculator such as https://geneu.xyz/user-profile

This program features multiple threads acording to user's CPU cores, default is 4 threads.

The required VCF files of WBBC could be download at: https://wbbc.westlake.edu.cn/downloads.html

### run

```
python main.py -tp tsv_tmpl -hld high_ld/high_ld_hg19.txt`

-tf: file name of TSV raw data
-tp: path of TSV raw data
-mp: path of admix model files
-af: file name of allele and frequency without extension
-ad: decimal digits of alleles frequency
-sd: threshold of standard deviation of alleles frequency
-hld: filename of regions of High Linkage Disequilibrium. Refer to: https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)
-th: threads for concurrency
```

### Citation

https://doi.org/10.1038/s41467-022-30526-x
