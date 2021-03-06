
## Author: Varun Warrier
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

# Step 3: Removing ancestry outliers
## Step 3a: Downloading and converting the Hapmap3 files

The next step is to retain only individuals who are primarily of European ancestry (CEU and TSI from HapMap 3). 

There are a few things we need to do for this:


We first need to download the hapmap3 files. This can be downloaded from here:
ftp://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/

```bash

wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/hapmap3_r1_b36_fwd.qc.poly.tar.bz2
bunzip2 hapmap3_r1_b36_fwd.qc.poly.tar.bz2
tar -xvf hapmap3_r1_b36_fwd.qc.poly.tar

```
This will create multiple files based on ancestry. For our current analyses, we need only the CEU and TSI files. We need to merge them using Plink. We use the --recode command to keep it in Ped/Map format as opposed to the default bed/bim/fam format as we need this in the next step for LiftOver.

```bash

./plink --file hapmap3_r1_b36_fwd.CEU.qc.poly.recode --merge hapmap3_r1_b36_fwd.TSI.qc.poly.recode --recode --out hapmap3_hg18_eur
```

The next step is to perform LiftOver on the merged file. The HapMap3 files are in NCBI build 36 (hg18). So we need to LiftOver to Hg19 for any subsequent analysis. 

```bash
python2 liftOverPlink.py -m path/to/map/file -p path/to/ped/file -e path/to/liftover -o outputname -c path/to/chainfile

# python2 liftOverPlink.py -m ./hapmap3_hg18_eur.map -p ./hapmap3_hg18_eur.ped -e ./liftOver -o hapmap3_hg19_eur -c ./hg18ToHg19.over.chain

```

Finally, we need to convert it to a binary bed file to merge it with the SSC files

```bash
./plink --file hapmap3_hg19_eur --make-bed --out hapmap3_hg19_eur
```

## Step 3b: Merging with the SSC files, running PCA, and removing ancestry outliers

To generate PCs ideally we need not closely individuals. Both the HapMap files and the SSC files are family based. So, we next create seperate files that has only founders.The filter-founders option excludes individuals who have atleast 1 parental ID. With this option, we are retaining only parents who we assume are not closely related. 

```bash
./plink --bfile hapmap3_hg19_eur --filter-founders --make-bed --out hapmap3_hg19_eurfoundersonly

./plink --bfile QC2output --filter-founders --make-bed --out QC2outputfoundersonly

```

After this, merge the HapMap3 file with the SSC file. This will be need to be done in a few steps usually.
The first merge command will produce a .missnp file. This is a list of SNPs that are either multiallelic or need to be flipped. 
We will first try to flip them, and then merge again.
This will additionally produce another list .missnp
We will then exclude this from both the files, and then merge again. 
We will recycle file names to save space. 

```bash
./plink --bfile hapmap3_hg19_eurfoundersonly --bmerge QC2outputfoundersonly --make-bed --out HapMap3SSCfileforPC

./plink --bfile QC2outputfoundersonly --flip HapMap3SSCfileforPC-merge.missnp --make-bed --out SSC_flippedfile

./plink --bfile hapmap3_hg19_eurfoundersonly --bmerge SSC_flippedfile --make-bed --out HapMap3SSCfileforPC

./plink --bfile hapmap3_hg19_eurfoundersonly --exclude HapMap3SSCfileforPC-merge.missnp --make-bed --out Hapmap3_formerging

./plink --bfile SSC_flippedfile --exclude HapMap3SSCfileforPC-merge.missnp --make-bed --out SSC_flippedfileformerging

./plink --bfile Hapmap3_formerging --bmerge SSC_flippedfileformerging --make_bed --out HapMap3SSCfileforPC

```


The next step is to generate PCs. To do this, we need another round of quality control, prune SNPs, and then generate PCs.

```bash
./plink --bfile HapMap3SSCfileforPC --geno 0.1 --hwe 0.000001 --make-bed --out HapMap3SSCfileforPCQC1

./plink --bfile HapMap3SSCfileforPCQC1 --maf 0.05 --indep-pairwise 100 50 0.2 --out HapMap3SSCfileforPCQC1pruned 

./plink --bfile HapMap3SSCfileforPCQC1 --exclude HapMap3SSCfileforPCQC1pruned.prune.out --pca --out SSC_pcaall
```



Phew! Now we get to plot the PCs and remove outliers. Nearly there...

Now, let's import the PCs into R and calculate things.

```R
library(data.table)
library(plyr)
pc = fread("SSC_pcaall.eigenvec")
selectedRows <- pc[grep("NA", pc$V2), ] #all the hapmapsamples have FID with NA

#Calculate the mean and SD of PC1 and PC2 based on the hapmapsamples
meanV1 = mean(selectedRows$V3)
sdV1 = sd(selectedRows$V3)
meanV2 = mean(selectedRows$V4)
sdV2 = sd(selectedRows$V4)

pc$ZPC1 = (abs(pc$V3 - meanV1))/sdV1
pc$ZPC2 = (abs(pc$V4 - meanV2))/sdV2

selectedRows2 <- pc[!grep("NA", pc$V2), ] #Now restrict it to the SSC samples
PCOK = subset(selectedRows2, ZPC1 < 5 & ZPC2 < 5) #include only samples that are less than 5 SDs away from the mean

count = count(PCOK$V1)
setnames(count, 1, "V1")
PCOK2 = merge(PCOK, count, by = "V1")
PCOKparents = subset(PCOK2, freq == "2")

famfile = fread("QC2output.fam")
keepfile = famfile[famfile$V1 %in% PCOKparents$V1,]
write.table(keepfile, file = "keepfile.txt", row.names = F, col.names = F, quote = F)

```

We've now got a list of individuals to keep in the SSC files for imputation in the keepfile. We can now include only these individuals and create a new bedfile for imputation. We then seperate it into the 22 autosomes and then create VCFs. Finally, we sort it and compress it.

```bash
./plink --bfile QC2output.fam --keep keepfile.txt --make-bed --out SSCimputationfile 
```

So far so good. But we now need to check if the files are valid. LiftOver can be a bit problematic, and luckily, there is a solution.

```bash
http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.9.zip
unzip HRC-1000G-check-bim-v4.2.9.zip

http://www.well.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz
gunzip 1000GP_Phase3_combined.legend.gz

./plink --bfile SSCimputationfile --freq --out SSCimputationfilefreq

perl HRC-1000G-check-bim.pl -b ./SSC_1Mv3/SSCimputationfile.bim -f SSCimputationfilefreq.frq -r 1000GP_Phase3_combined.legend -g -p EUR

chmod +x Finalstep.sh

run Finalstep.sh

```



Et, viola! You are now done. Upload the files onto your favourite imputation servers and pray to the gods of the interenet that it works!

I used the Michigan Imputation Server. Imputation was conducted using 1000G Phase 3 V5. Phasing was done using Eagle v2.3. Population was restricted to EUR. The mode was Quality Control and imputation.

## Post-imputation QC and data management

Now, it's not all done from here, and there are a few additional steps. First, you need to unzip the files. This is password protected. 

Next, once you have the VCFs, convert it into plink binary, and unzip the info file. It's best to unzip the info file and create a list of files to exclude before converting it to plink binary, as this way, we can reduce an extra step.

```bash

```

Let's now do some QC. Here, we are removing SNPs with 0.5 < ALT_Freq < 0.95. This will ensure that all our SNPs have a minor allele freq > 0.05. We also want to remove poorly imputed SNPs. We do this by removing SNPs with Rsq < 0.6. Finally, we create a file to that can be passed on to the --extract command in Plink. 

```R
for (i in 1:22){
  a = read.table(paste0("chr", i, ".info"), header = T)
  a$Rsq = as.numeric(as.character(a$Rsq))
  b = subset(a, Rsq > 0.6)
  b = subset(b, ALT_Frq > 0.05 & ALT_Frq < 0.95)
  write.table(b[,1], file = paste0("chr", i, "exclude.txt"), row.names = F, col.names = T, quote = F)
}
```

Now we create the binary Plink file. This does 
```bash

for i in {1..22}; do ./plink --vcf ./Imputed/chr${i}.dose.vcf.gz --make-bed --out ./imputed_plinkfile/1Mv1_chr${i}  --extract ./Imputed/chr${i}exclude.txt --const-fid 0; done

```

Next, let's combine all the files, and recode the SNP IDs. 
To recode the SNP IDs, you need to download the VCF file that matches your build (GRCh37), and then manipulate it to get the file you need to pass that onto plink. We will restrict it to only the common variants. 

```bash
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/common_all_20170710.vcf.gz
zgrep -v "^##" common_all_20170710.vcf.gz | cut -f1-3 > fileforrecoding.txt
awk '{print $1":"$2"\t"$3}' < fileforrecoding.txt > plinkrecodingfile.txt
```


Merge the files, update SNP name, and do QC
```bash
i in {1..22}; do ./plink --bfile ./imputed_plinkfile/1Mv1imputed_chr${i} --extract ./Imputed_1000G/chr${i}exclude.txt --make-bed --maf 0.05 --geno 0.05 --hwe 0.000001 --out ./imputed_plinkfile/1Mv1_chr${i}; done
      
./plink --bfile 1Mv1_chr22 --merge-list 1Mv1mergelist.txt --make-bed -biallelic-only --out 1Mv1_merged

for i in {1..22}; do ./plink --bfile 1Mv1_chr${i} --exclude 1Mv1_merged-merge.missnp --make-bed --out 1Mv1_chr${i}_2; done

./plink --bfile 1Mv1_chr22_2 --merge-list 1Mv1mergelist2.txt --make-bed -biallelic-only --out 1Mv1_merged

./plink --bfile 1Mv1_merged --maf 0.05 --update-name ~/SFARI/liftOverPlink/plinkrecodingfile.txt --hwe 0.000001 --geno 0.05 --mind 0.05 --make-bed --out 1Mv1_mergedQC

```

Finally, update the Fam file, as this gets messed up in the whole process. To do this, open the fam file of the imputed merged version, and the fam file of the non-imputed version in R.

```R
library(data.table)
library(tidyr)
setwd()

fileimputed = fread("")
filenonimputed = fread("")

fileimputed = fileimputed %>% separate(V2, into = c('FID', 'IID'), sep = 6)
setnames(filenonimputed, 2, "IID")

merged = merge(fileimputed, filenonimputed, by = "IID")

merged$oldIID = paste0(merged$FID, merged$IID)

setnames(merged, "IID", "NewIID")
setnames(merged, "V1.x", "oldFID")
setnames(merged, "V1.y", "newFID")

write.table(merged[,c("oldFID", "oldIID", "newFID", "NewIID")], file = "1Mv1updatenames.txt", row.names = F, col.names = F, quote = F)
write.table(filenonimputed[,1:4], file = "1Mv1updateparents.txt", row.names = F, col.names = F, quote = F)
write.table(filenonimputed[,c(1,2,5)], file = "1Mv1updatesex.txt", row.names = F, col.names = F, quote = F)
write.table(filenonimputed[,c(1,2,6)], file = "1Mv1updatepheno.txt", row.names = F, col.names = F, quote = F)

```

```bash
./plink --bfile 1Mv1_mergedQC --update-ids 1Mv1updatenames.txt --make-bed  --out 1Mv1_mergedQC2

./plink --bfile 1Mv1_mergedQC2 --update-parents 1Mv1updateparents.txt --update-sex 1Mv1updatesex.txt --pheno 1Mv1updatepheno.txt --make-bed  --out 1Mv1_mergedQC2
```


## Resources:

1. liftOver: https://genome.sph.umich.edu/wiki/LiftOver
2. liftOverPlink wrapper: https://github.com/sritchie73/liftOverPlink
3. Plink 1.9 :https://www.cog-genomics.org/plink2/
4. Plink 2.0: https://www.cog-genomics.org/plink/2.0/
5. HapMap files: ftp://ftp.ncbi.nlm.nih.gov/hapmap/
6. Imputation servers:https://imputationserver.sph.umich.edu/index.html and https://imputation.sanger.ac.uk/
7. General information on imputation and QC: https://sites.google.com/a/broadinstitute.org/ricopili/
8. A good primer on data QC for GWAS: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3025522/



