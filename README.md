# Create Admixture model files from WBBC project of West Lake University
Create admixture model files (*.F, *.alleles) based on WBBC project of West Lake University. Microarray chip files (TSV format) from popular genetic testing companies are used as referenced SNPs. Linkage disequilibrium SNPs could be eliminated optionally. The output admixture model files could be used for model based admixture calculator such as https://geneu.xyz/user-profile

This program features multiple threads acording to user's CPU cores, default is 4 threads.

The required VCF files of WBBC could be download at: https://wbbc.westlake.edu.cn/downloads.html

run: python wbbc.py -tp tsv_tmpl -hld high_ld/high_ld_hg19.txt

ref: https://doi.org/10.1038/s41467-022-30526-x
