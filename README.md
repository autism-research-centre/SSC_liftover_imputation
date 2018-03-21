## SSC_liftover_imputation
# Basic scripts used for imputing the SSC genotyped datasets


These are basic steps to conduct imputation in the SSC dataset. Note, this can be modified and used for any dataset. 

There are 3 SSC files. Steps 1 - 3 must be done for each of the three files seperately. 

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

 # Step 2: Quality control
 
 The QC steps are largely done using a combination of Plink/ Plink2 and R. The dataset is in trios so that's convenient to check for mendelian errors.

The next chunk of code deals with removing 
 - individuals who have low genotyping specifically, if greater than 5% of the SNPs in each individual is missing (--mind 0.05)
 - families with mendelian errors > 5% (--me **0.05** 0.1) (Families where more than 5% of the variants have evidence for non-mendelian transmission. Likely reflects genotyping errors. )
 - SNPs that significantly deviate from hardy-weinberg equilibrium (--hwe 0.000001)
 - SNPs that low genotype rate, specifically, if they are not genotyped in more than 10% of the individuals (--geno 0.01)
 - SNPs with mendelian errors > 10% (--me 0.5 **0.1**)

We take the ped/map file, run these commands, and then write the results out as a bedfile called QC1output

After this, we check for individuals with excessively high or low heterozygosity, and individuals with discordant sex information, and write this information out as QC1checkhet and QC1checksex respectively. 

The sexcheck essentially compares the X chromosome dosage to the reported sex. 
 
 ```bash
 ./plink --file path/to/inputfile --geno 0.1 --mind 0.05 --hwe 0.000001 --me 0.05 0.1 --make-bed --out QC1output
 
 ./plink --bfile QC1output --het --out QC1checkhet
 
 ./plink --bfile QC1output --check-sex --out QC1checksex
 
```
Note, sometimes liftover may throw in some odd chromosomes/snps. This will get Plink to give the following error message:
" Invalid chromosome code '11_gl000202_random' on line 1563099 of .map"

As far as I am aware there is no quick way around it. Read the map file in R. Check for invalid chromosomes (outside 1:22, X, Y), and exclude the SNPs in those using the --exclude command when creating the bfile. Alternatively, use the --allow-extra-chr command to ignore this in Plink. 

The next step is to create a file with problematic individuals and remove it from subsequent analysis. This is done in R.


```R
sex = read.delim("QC1checksex.sexcheck", sep = "")
head(sex)
sexproblem = subset(sex, STATUS == "PROBLEM")


het = read.delim("QC1checkhet.het", sep = "")
het$HET = (het$N.NM. - het$O.HOM.)/het$N.NM. #create heterozygosity stats

mean = mean(het$HET)
sd = sd(het$HET)
het$Z = (het$HET - mean)/sd #create Z scores of heterozygosity

hetoutlier = subset(het, abs(Z) > 3)

het2 = hetoutlier[,c(1:2)]
sex2 = sexproblem[,c(1:2)]
failedsample = rbind(het2, sex2)

write.table(failedsample, file = "failedsample.txt", row.names = F, col.names = T, quote = F)
```

Next, we remove the failed samples.

```bash

./plink --bfile QC1output --remove failedsample.txt --make-bed --out QC2output

```


./plink --bfile SSC_Omni2.5_binary_QC2 --filter-founders --make-bed --out SSC_Omni2.5_binary_QC2foundersonly

./plink --bfile SSC_1Mv3_binary_QC2 --filter-founders --make-bed --out SSC_1Mv3_binary_QC2foundersonly

./plink --bfile SSC_1Mv1_binary_QC2 --filter-founders --make-bed --out SSC_1Mv1_binary_QC2foundersonly

./plink --bfile merged --bmerge  SSC_1Mv1_binary_QC1 --out merged2

./plink --bfile merged2 --filter-founders --make-bed --out merged2founders 







## Resources:

1. liftOver: https://genome.sph.umich.edu/wiki/LiftOver
2. liftOverPlink wrapper: https://github.com/sritchie73/liftOverPlink
3. Plink 1.9 :https://www.cog-genomics.org/plink2/
4. Plink 2.0: https://www.cog-genomics.org/plink/2.0/
5. HapMap files: ftp://ftp.ncbi.nlm.nih.gov/hapmap/
6. Imputation servers:https://imputationserver.sph.umich.edu/index.html and https://imputation.sanger.ac.uk/
7. General information on imputation and QC: https://sites.google.com/a/broadinstitute.org/ricopili/
8. A good primer on data QC for GWAS: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3025522/



