## SSC_liftover_imputation
#Basic scripts used for imputing the SSC genotyped datasets


These are basic steps to conduct imputation in the SSC dataset. Note, this can be modified and used for any dataset. 

In general, the steps are as follows:

#1. Liftover to hg19 from hg18 (all SSC files are in hg18)
#2. Quality control
#3. Generation of genetic PCs and checking for and removing outliers against the HapMap dataset
#4. Imputation (we will do this on the imputation servers).









##Resources:

1. liftOver: https://genome.sph.umich.edu/wiki/LiftOver
2. liftOverPlink wrapper: https://github.com/sritchie73/liftOverPlink
3. Plink 1.9 :https://www.cog-genomics.org/plink2/
4. Plink 2.0: https://www.cog-genomics.org/plink/2.0/
5. HapMap files: ftp://ftp.ncbi.nlm.nih.gov/hapmap/
6. Imputation servers:https://imputationserver.sph.umich.edu/index.html and https://imputation.sanger.ac.uk/
7. General information on imputation and QC: https://sites.google.com/a/broadinstitute.org/ricopili/



