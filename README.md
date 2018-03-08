## SSC_liftover_imputation
# Basic scripts used for imputing the SSC genotyped datasets


These are basic steps to conduct imputation in the SSC dataset. Note, this can be modified and used for any dataset. 

In general, the steps are as follows:

 1. LiftOver to hg19 from hg18 (all SSC files are in hg18)
 2. Quality control
 3. Generation of genetic PCs and checking for and removing outliers against the HapMap dataset
 4. Imputation (we will do this on the imputation servers)
 
 
 # Step 1: LiftOver
 
The SSC files are on hg18, which makes it tricky to do any downstream analysis. Imputation servers require hg19 build. So the first step is to convert it from hg18 to hg19.

The files provided are in Plink format (map and ped), and LiftOver works on the UCSC bed format. Luckily, tools have been developed to get around this (LiftOverPlink, see resources section). 

The LiftOverPlink is a really nifty tool. It converts the map file to a UCSC bed file, runs liftOver, and then updates the ped and map files. The most time-consuming step is the final one, updating the ped file. 

```bash
#install liftOver
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver

#get the chain file
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
gunzip hg18ToHg19.over.chain.gz

#get the LiftOverPlink
git clone https://github.com/sritchie73/liftOverPlink.git

#run the programme
python2 liftOverPlink.py -m path/to/map/file -p path/to/ped/file -e path/to/liftover -o outputname -c path/to/chainfile

```









## Resources:

1. liftOver: https://genome.sph.umich.edu/wiki/LiftOver
2. liftOverPlink wrapper: https://github.com/sritchie73/liftOverPlink
3. Plink 1.9 :https://www.cog-genomics.org/plink2/
4. Plink 2.0: https://www.cog-genomics.org/plink/2.0/
5. HapMap files: ftp://ftp.ncbi.nlm.nih.gov/hapmap/
6. Imputation servers:https://imputationserver.sph.umich.edu/index.html and https://imputation.sanger.ac.uk/
7. General information on imputation and QC: https://sites.google.com/a/broadinstitute.org/ricopili/



