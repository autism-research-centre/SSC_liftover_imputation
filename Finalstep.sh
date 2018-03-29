./plink --bfile SSCimputationfile --exclude Exclude-SSCimputationfile-1000G.txt --make-bed --out TEMP1
./plink --bfile TEMP1 --update-map Chromosome-SSCimputationfile-1000G.txt --update-chr --make-bed --out TEMP2
./plink --bfile TEMP2 --update-map Position-SSCimputationfile-1000G.txt --make-bed --out TEMP3
./plink --bfile TEMP3 --flip Strand-Flip-SSCimputationfile-1000G.txt --make-bed --out TEMP4
./plink --bfile TEMP4 --reference-allele Force-Allele1-SSCimputationfile-1000G.txt --make-bed --out SSCimputationfile-updated
rm TEMP*
for i in {1..22}; do ./plink --bfile SSCimputationfile-updated --reference-allele Force-Allele1-SSCimputationfile-1000G.txt --chr ${i} --recode-vcf --out SSC_file_chr${i}; done
for i in {1..22}; do vcf-sort SSC_file_chr${i}.vcf | bgzip -c > SSC_file_chr${i}.vcf.gz; done
for i in {1..22}; do rm SSC_file_chr${i}.vcf | rm SSC_file_chr${i}.log; done
