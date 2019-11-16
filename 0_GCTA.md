


###FOR EACH SEPARATE RANGE, miss50

filt_maf05.miss50.het0.7.good.vcf

module load java/jdk1.7.0_21
module load gatk/3.5.0
java -jar /opt/sw/gatk/3.5.0/GenomeAnalysisTK.jar \
   -T SelectVariants \
-R /nfs/shares/hodgins-grp/ragweed/HiRise_v2/ragweed_10Feb2016_2ABsE_uppercase70.fasta  \
   -V filt_maf05.miss50.het0.7.good.vcf \
   -o filt_maf05.miss50.het0.7.good.AU.vcf \
--ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES \
   -sf au.list

java -jar /opt/sw/gatk/3.5.0/GenomeAnalysisTK.jar \
   -T SelectVariants \
-R /nfs/shares/hodgins-grp/ragweed/HiRise_v2/ragweed_10Feb2016_2ABsE_uppercase70.fasta  \
   -V filt_maf05.miss50.het0.7.good.vcf \
   -o filt_maf05.miss50.het0.7.good.AU.vcf \
--ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES \
   -sf au.list

java -jar /opt/sw/gatk/3.5.0/GenomeAnalysisTK.jar \
   -T SelectVariants \
-R /nfs/shares/hodgins-grp/ragweed/HiRise_v2/ragweed_10Feb2016_2ABsE_uppercase70.fasta  \
   -V filt_maf05.miss50.het0.7.good.vcf \
   -o filt_maf05.miss50.het0.7.good.EU.vcf \
--ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES \
   -sf eu.list

java -jar /opt/sw/gatk/3.5.0/GenomeAnalysisTK.jar \
   -T SelectVariants \
-R /nfs/shares/hodgins-grp/ragweed/HiRise_v2/ragweed_10Feb2016_2ABsE_uppercase70.fasta  \
   -V filt_maf05.miss50.het0.7.good.vcf \
   -o filt_maf05.miss50.het0.7.good.NA.vcf \
--ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES \
   -sf na.list



cat filt_maf05.miss50.het0.7.good.NA.vcf | cut -f1 | grep -v '#' | uniq > tmp
awk '{print $1, 1}' tmp > keys
cat filt_maf05.miss50.het0.7.good.NA.vcf | grep -v '#' > tmp2
awk '{print $1}' tmp2 > tmp3
awk 'NR==FNR{a[$1]=$2} NR>FNR{$1=a[$1];print}' keys tmp3 > numbercolumnNA
cut -f2- tmp2 > tmp4
paste -d"\t" numbercolumnNA tmp4 > tmp5
cat filt_maf05.miss50.het0.7.good.NA.vcf | grep '#' > filt_maf05.miss50.het0.7.numscafgood.NA.vcf
cat tmp5 >> filt_maf05.miss50.het0.7.numscafgood.NA.vcf


#get the loci names
cat filt_maf05.miss50.het0.7.good.NA.vcf  | cut -f1-2 | grep -v '#' > locinamesNA

####rename SampleID
cat filt_maf05.miss50.het0.7.numscafgood.NA.vcf | grep '##' > filt_maf05.miss50.het0.7.numscafgood_RNID.NA.vcf
cat filt_maf05.miss50.het0.7.numscafgood.NA.vcf| grep '#' > tmp
cat tmp | grep -v '##' > tmp2

-edit tmp2 (rename to tmp3)

cat tmp3 >> filt_maf05.miss50.het0.7.numscafgood_RNID.vcf
cat filt_maf05.miss50.het0.7.numscafgood.vcf | grep -v '#' >> filt_maf05.miss50.het0.7.numscafgood_RNID.vcf


#### for full GRM
module load vcftools
vcftools --vcf filt_maf05.miss50.het0.7.numscafgood_RNID.vcf --plink --out filt_maf05.miss50.het0.7
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --file filt_maf05.miss50.het0.7 --make-bed --out filt_maf05.miss50.het0.7 --noweb

/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --noweb --bfile filt_maf05.miss50.het0.7 --indep-pairwise 3 10 .1
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --noweb --bfile filt_maf05.miss50.het0.7 --extract plink.prune.in 

/nfs/shares/hodgins-grp/bin/gcta/gcta64  --bfile filt_maf05.miss50.het0.7 --make-grm --extract plink.prune.in --out filt_maf05.miss50.het0.7_sub
/nfs/shares/hodgins-grp/bin/gcta/gcta64 --grm filt_maf05.miss50.het0.7_sub --grm-cutoff 0.025 --make-grm --out grm_adj
/nfs/shares/hodgins-grp/bin/gcta/gcta64  --grm grm_adj  --pca 10 --out filt_maf05.miss50.het0.7.PCA
/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss50.het0.7 --grm grm_adj --qcovar filt_maf05.miss50.het0.7.PCA.eigenvec --pheno filt_maf05.miss50.het0.7.pheno --out filt_maf05.miss50.het0.7_filters


for i in *pheno; do /nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss50.het0.7NA --grm grm_adjNA --qcovar filt_maf05.miss50.het0.7.NA.PCA.eigenvec --pheno $i --out ${i/.pheno/} > ${i/.pheno/}.log 2>&1; done

- find if they all converged (print list of files not containing the words'
grep -L "Log-likelihood ratio converged" *log




###



#####FOR 60% MISSING DATA###
vcftools --remove bad.list --vcf filt_maf05.miss90.het0.7.vcf --recode --out filt_maf05.miss90.het0.7
mv filt_maf05.miss90.het0.7.recode.vcf filt_maf05.miss90.het0.7.good.vcf
vcftools --vcf filt_maf05.miss90.het0.7.good.vcf --max-missing 0.4 --recode --out filt_maf05.miss60.het0.7.vcf

mv filt_maf05.miss60.het0.7.vcf.recode.vcf filt_maf05.miss60.het0.7.good.vcf

cat filt_maf05.miss60.het0.7.good.vcf  | cut -f1 | grep -v '#' | uniq > tmp
awk '{print $1, 1}' tmp > keys
cat filt_maf05.miss60.het0.7.good.vcf | grep -v '#' > tmp2
awk '{print $1}' tmp2 > tmp3
awk 'NR==FNR{a[$1]=$2} NR>FNR{$1=a[$1];print}' keys tmp3 > numbercolumn
cut -f2- tmp2 > tmp4
paste -d"\t" numbercolumn tmp4 > tmp5
cat filt_maf05.miss60.het0.7.good.vcf | grep '#' > filt_maf05.miss60.het0.7.numscafgood.vcf
cat tmp5 >> filt_maf05.miss60.het0.7.numscaf.vcf

#get the loci names
cat filt_maf05.miss60.het0.7.good.vcf  | cut -f1-2 | grep -v '#' > locinames

####rename SampleID
cat filt_maf05.miss60.het0.7.numscafgood.vcf | grep '##' > filt_maf05.miss60.het0.7.numscafgood_RNID.vcf
cat filt_maf05.miss60.het0.7.numscafgood.vcf| grep '#' > tmp
cat tmp | grep -v '##' > tmp2

-edit tmp2 (rename to tmp3) --> download and change all names to i1, i2 etc

cat tmp3 >> filt_maf05.miss60.het0.7.numscafgood_RNID.vcf
cat filt_maf05.miss60.het0.7.numscafgood.vcf | grep -v '#' >> filt_maf05.miss60.het0.7.numscafgood_RNID.vcf


#### for full GRM
module load vcftools
vcftools --vcf filt_maf05.miss60.het0.7.numscafgood_RNID.vcf --plink --out filt_maf05.miss60.het0.7
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --file filt_maf05.miss60.het0.7 --make-bed --out filt_maf05.miss60.het0.7 --noweb

/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --noweb --bfile filt_maf05.miss60.het0.7 --indep-pairwise 3 10 .1
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --noweb --bfile filt_maf05.miss60.het0.7 --extract plink.prune.in 

/nfs/shares/hodgins-grp/bin/gcta/gcta64  --bfile filt_maf05.miss60.het0.7 --make-grm --extract plink.prune.in --out filt_maf05.miss60.het0.7_sub
/nfs/shares/hodgins-grp/bin/gcta/gcta64 --grm filt_maf05.miss60.het0.7_sub --grm-cutoff 0.025 --make-grm --out grm_adj
/nfs/shares/hodgins-grp/bin/gcta/gcta64  --grm grm_adj  --pca 10 --out filt_maf05.miss60.het0.7.PCA
/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss60.het0.7 --grm grm_adj --qcovar filt_maf05.miss60.het0.7.PCA.eigenvec --pheno filt_maf05.miss60.het0.7.pheno --out filt_maf05.miss60.het0.7_filters

# without additional grm linkage and relatedness filtering...
/nfs/shares/hodgins-grp/bin/gcta/gcta64  --bfile filt_maf05.miss60.het0.7 --make-grm --out filt_maf05.miss60.het0.7_sub
/nfs/shares/hodgins-grp/bin/gcta/gcta64  --grm filt_maf05.miss60.het0.7_sub --pca 10 --out filt_maf05.miss60.het0.7.PCA
/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss60.het0.7 --grm filt_maf05.miss60.het0.7_sub --qcovar filt_maf05.miss60.het0.7.PCA.eigenvec --pheno filt_maf05.miss60.het0.7.pheno --out filt_maf05.miss60.het0.7

/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss60.het0.7 --grm filt_maf05.miss60.het0.7_sub --qcovar filt_maf05.miss60.het0.7.PCA.eigenvec --pheno ${f} --out ${f}

To print the GRM:

ReadGRMBin=function(prefix, AllN=F, size=4){
 sum_i=function(i){
 return((1+i)*i/2)
 }

  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}

GRM <- ReadGRMBin("filt_maf05.miss50.het0.7_N_sub")
write.table(GRM,"GRM.txt")





####FOR 50% MISSING DATA, with explanations#####
vcftools --remove bad.list --vcf filt_maf05.miss50.het0.7.vcf --recode --out filt_maf05.miss50.het0.7
mv filt_maf05.miss50.het0.7.recode.vcf filt_maf05.miss50.het0.7.good.vcf

#### for full GRM
module load vcftools
vcftools --vcf filt_maf05.miss50.het0.7.numscafgood_RNID.vcf --plink --out filt_maf05.miss50.het0.7
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --file filt_maf05.miss50.het0.7 --make-bed --out filt_maf05.miss50.het0.7 --noweb

/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --noweb --bfile filt_maf05.miss50.het0.7 --indep-pairwise 3 10 .1
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --noweb --bfile filt_maf05.miss50.het0.7 --extract plink.prune.in 

/nfs/shares/hodgins-grp/bin/gcta/gcta64  --bfile filt_maf05.miss50.het0.7 --make-grm --extract plink.prune.in --out filt_maf05.miss50.het0.7_sub
/nfs/shares/hodgins-grp/bin/gcta/gcta64 --grm filt_maf05.miss50.het0.7_sub --grm-cutoff 0.025 --make-grm --out grm_adj
/nfs/shares/hodgins-grp/bin/gcta/gcta64  --grm grm_adj  --pca 10 --out filt_maf05.miss50.het0.7.PCA
/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss50.het0.7 --grm grm_adj --qcovar filt_maf05.miss50.het0.7.PCA.eigenvec --pheno filt_maf05.miss50.het0.7.pheno --out filt_maf05.miss50.het0.7_filters

# without additional grm linkage and relatedness filtering...
/nfs/shares/hodgins-grp/bin/gcta/gcta64  --bfile filt_maf05.miss50.het0.7 --make-grm --out filt_maf05.miss50.het0.7_sub
/nfs/shares/hodgins-grp/bin/gcta/gcta64  --grm filt_maf05.miss50.het0.7_sub --pca 10 --out filt_maf05.miss50.het0.7.PCA
/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss50.het0.7 --grm filt_maf05.miss50.het0.7_sub --qcovar filt_maf05.miss50.het0.7.PCA.eigenvec --pheno filt_maf05.miss50.het0.7.pheno --out filt_maf05.miss50.het0.7



####FOR 90% MISSING DATA, with explanations#####

-copy the vcf
cp filt_maf05.miss90.het0.7.randunl.vcf  vcf

-Plink isn't able to deal with scaffold names, so create a key file with the numbers
cat filt_maf05.miss90.het0.7.randunl.vcf  | cut -f1 | grep -v '##' | uniq > tmp

- in vi, delete the first line (#CHROM)
- add the number 1 to this file
awk '{print $1, 1}' tmp > keys2

- remove the ## header from the vcf file
cat filt_maf05.miss90.het0.7.randunl.vcf | grep -v '#' > tmp5

-Get the first column of the file
awk '{print $1}' tmp5 > tmp6

- change only the first column of the vcf: in the 'keys' file (with scaffold names and corresponding numbers), find the scaffold name
in column 1, make an array a with corresponding column 2. In the 'tmp6' file, find this and replace scaffold name with number
awk 'NR==FNR{a[$1]=$2} NR>FNR{$1=a[$1];print}' keys2 tmp6 > numbercolumn

-remove first column from tmp5 and add the number column to this
cut -f2- tmp5 > tmp7
paste -d"\t" numbercolumn tmp7 > tmp8

- Cat only the vcf header to the new file
cat filt_maf05.miss90.het0.7.randunl.vcf | grep '#' > filt_maf05.miss90.het0.7.randunl_numscaf.vcf

- cat tmp8  to new file
cat tmp8 >> filt_maf05.miss90.het0.7.randunl_numscaf.vcf

-Plink isn't able to deal with scaffold names, so create a key file with the numbers
cat filt_maf05.miss90.het0.7.vcf | cut -f1 | grep -v '##' | uniq > tmp

- in vi, delete the first line (#CHROM)
- add the number 1 to this file
awk '{print $1, 1}' tmp > keys2

- remove the ## header from the vcf file
cat vcf | grep -v '#' > tmp5

-Get the first column of the file
awk '{print $1}' tmp5 > tmp6

- change only the first column of the vcf: in the 'keys' file (with scaffold names and corresponding numbers), find the scaffold name
in column 1, make an array a with corresponding column 2. In the 'tmp6' file, find this and replace scaffold name with number
awk 'NR==FNR{a[$1]=$2} NR>FNR{$1=a[$1];print}' keys2 tmp6 > numbercolumn

-remove first column from tmp5 and add the number column to this
cut -f2- tmp5 > tmp7
paste -d"\t" numbercolumn tmp7 > tmp8

- Cat only the vcf header to the new file
cat filt_maf05.miss90.het0.7.randunl.vcf | grep '#' > filt_maf05.miss90.het0.7.numscaf.vcf

- cat tmp8  to new file
cat tmp8 >> filt_maf05.miss90.het0.7.numscaf.vcf

- remove individuals with >90% missing genotypes
vcftools --remove bad.list --vcf filt_maf05.miss90.het0.7.randunl_numscaf.vcf --recode --out filt_maf05.miss90.het0.7.randunl_numscaf2.vcf 
vcftools --remove bad.list --vcf filt_maf05.miss90.het0.7.numscaf.vcf --recode --out filt_maf05.miss90.het0.7.numscaf2.vcf

mv filt_maf05.miss90.het0.7.randunl_numscaf2.vcf.recode.vcf filt_maf05.miss90.het0.7.randomunl_numscafgood.vcf
mv filt_maf05.miss90.het0.7.numscaf2.vcf.recode.vcf  filt_maf05.miss90.het0.7.numscafgood.vcf 

#rename SampleID
cat filt_maf05.miss90.het0.7.randomunl_numscafgood.vcf | grep '##' > filt_maf05.miss90.het0.7.randomunl_numscafgood_RNID.vcf
cat filt_maf05.miss90.het0.7.randomunl_numscafgood.vcf | grep '#' > tmp
cat tmp | grep -v '##' > tmp2
-edit tmp2 (rename to tmp3)
cat tmp3 >> filt_maf05.miss90.het0.7.randomunl_numscafgood_RNID.vcf
cat filt_maf05.miss90.het0.7.randomunl_numscafgood.vcf | grep -v '#' >> filt_maf05.miss90.het0.7.randomunl_numscafgood_RNID.vcf

cat filt_maf05.miss90.het0.7.numscafgood.vcf | grep '##' > filt_maf05.miss90.het0.7.numscafgood_RNID.vcf
cat tmp3 >> filt_maf05.miss90.het0.7.numscafgood_RNID.vcf
cat filt_maf05.miss90.het0.7.numscafgood.vcf | grep -v '#' >> filt_maf05.miss90.het0.7.numscafgood_RNID.vcf




module load vcftools
vcftools --vcf filt_maf05.miss90.het0.7.randomunl_numscafgood_RNID.vcf --plink --out filt_maf05.miss90.het0.7_N
vcftools --vcf filt_maf05.miss90.het0.7.numscafgood_RNID.vcf --plink --out filt_maf05.miss90.het0.7

/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --file filt_maf05.miss90.het0.7 --make-bed --out filt_maf05.miss90.het0.7 --noweb
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --file filt_maf05.miss90.het0.7_N --make-bed --out filt_maf05.miss90.het0.7_N --noweb
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --noweb --bfile filt_maf05.miss90.het0.7_N --indep-pairwise 3 10 .1
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --noweb --bfile filt_maf05.miss90.het0.7_N --extract plink.prune.in 

#grm based on additional linkage filtering...
/nfs/shares/hodgins-grp/bin/gcta/gcta64  --bfile filt_maf05.miss90.het0.7_N --make-grm --extract plink.prune.in --out filt_maf05.miss90.het0.7_N_sub

/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss90.het0.7 --grm filt_maf05.miss90.het0.7_N_sub --pheno filt_maf05.miss90.het0.7.pheno --out filt_maf05.miss90.het0.7

/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss90.het0.7 --pheno filt_maf05.miss90.het0.7.pheno --out filt_maf05.miss90.het0.7


#add a pca component
/nfs/shares/hodgins-grp/bin/gcta/gcta64  --grm filt_maf05.miss90.het0.7_N_sub --pca 10 --out filt_maf05.miss90.het0.7.PCA
#run again
/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss90.het0.7 --grm filt_maf05.miss90.het0.7_N_sub --qcovar filt_maf05.miss90.het0.7.PCA.eigenvec --pheno filt_maf05.miss90.het0.7.pheno --out filt_maf05.miss90.het0.7

#Increase # iterations
--reml-maxit 5000

/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss90.het0.7 --reml-maxit 5000 --grm filt_maf05.miss90.het0.7_N_sub --qcovar filt_maf05.miss90.het0.7.PCA.eigenvec --pheno filt_maf05.miss90.het0.7.pheno --out filt_maf05.miss90.het0.7
/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss90.het0.7 --reml-maxit 5000 --grm filt_maf05.miss90.het0.7_N_sub --pheno filt_maf05.miss90.het0.7.pheno --out filt_maf05.miss90.het0.7


# for full GRM
module load vcftools
vcftools --vcf filt_maf05.miss90.het0.7.numscafgood_RNID.vcf --plink --out filt_maf05.miss90.het0.7
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --file filt_maf05.miss90.het0.7 --make-bed --out filt_maf05.miss90.het0.7 --noweb

#grm based on additional linkage filtering...
/nfs/shares/hodgins-grp/bin/gcta/gcta64  --bfile filt_maf05.miss90.het0.7 --make-grm --out filt_maf05.miss90.het0.7_sub
/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss90.het0.7 --grm filt_maf05.miss90.het0.7_sub --pheno filt_maf05.miss90.het0.7.pheno --out filt_maf05.miss90.het0.7


#add a pca component
/nfs/shares/hodgins-grp/bin/gcta/gcta64  --grm filt_maf05.miss90.het0.7_sub --pca 10 --out filt_maf05.miss90.het0.7.PCA
#run again
/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss90.het0.7 --grm filt_maf05.miss90.het0.7_sub --qcovar filt_maf05.miss90.het0.7.PCA.eigenvec --pheno filt_maf05.miss90.het0.7.pheno --out filt_maf05.miss90.het0.7

/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss90.het0.7 --reml-maxit 5000 --reml-bendV --grm filt_maf05.miss90.het0.7_sub --qcovar filt_maf05.miss90.het0.7.PCA.eigenvec --pheno filt_maf05.miss90.het0.7.pheno --out filt_maf05.miss90.het0.7





#### WITH PRIOR LD FILTERING ####

vcftools --vcf filt_maf05.miss90.het0.7.vcf --geno-r2 --min-r2 0.5 --ld-window 5 --ld-window-bp 50 --out filt_maf05.miss90.het0.7.vcf_ldresults5
cat filt_maf05.miss90.het0.7.vcf_ldresults5.geno.ld| cut -f1-2 | grep 'SC' > snps.0.5.90miss_hetero0.7.LD505.5
grep -vwF -f snps.0.5.90miss_hetero0.7.LD505.5 filt_maf05.miss90.het0.7.vcf > filt_maf05.miss90.het0.7.LD505.5.vcf
cat filt_maf05.miss90.het0.7.LD505.5.vcf | grep -v '#' | wc -l
68391

vcftools --vcf filt_maf05.miss90.het0.7.vcf --geno-r2 --min-r2 0.1 --ld-window 5 --ld-window-bp 50 --out filt_maf05.miss90.het0.7.vcf_ldresults5
cat filt_maf05.miss90.het0.7.vcf_ldresults5.geno.ld| cut -f1-2 | grep 'SC' > snps.0.5.90miss_hetero0.7.LD505.1
grep -vwF -f snps.0.5.90miss_hetero0.7.LD505.1 filt_maf05.miss90.het0.7.vcf > filt_maf05.miss90.het0.7.LD505.1.vcf
cat filt_maf05.miss90.het0.7.LD505.1.vcf | grep -v '#' | wc -l
54012


copy the vcf
cp filt_maf05.miss90.het0.7.LD505.1.vcf vcf

-Plink isn't able to deal with scaffold names or more than 95 scaffolds, so create a key file with the numbers
cat vcf | cut -f1 | grep -v '##' | uniq > tmp

- in vi, delete the first line (#CHROM)

- add the number 1 to this file
awk '{print $1, 1}' tmp > keys2

- remove the ## header from the vcf file
cat vcf | grep -v '#' > tmp5

-Get the first column of the file
awk '{print $1}' tmp5 > tmp6

- change only the first column of the vcf: in the 'keys' file (with scaffold names and corresponding numbers), find the scaffold name
in column 1, make an array a with corresponding column 2. In the 'tmp6' file, find this and replace scaffold name with number
awk 'NR==FNR{a[$1]=$2} NR>FNR{$1=a[$1];print}' keys2 tmp6 > numbercolumn

-remove first column from tmp5 and add the number column to this
cut -f2- tmp5 > tmp7
paste -d"\t" numbercolumn tmp7 > tmp8

- Cat only the vcf header to the new file
cat vcf | grep '#' > filt_maf05.miss90.het0.7.LD505.1._scaf1.vcf

- cat tmp8  to new file
cat tmp8 >> filt_maf05.miss90.het0.7.LD505.1._scaf1.vcf

module load vcftools
vcftools --vcf filt_maf05.miss90.het0.7.LD505.1._scaf1.vcf --plink --out filt_maf05.miss90.het0.7.LD505.1
vcftools --vcf filt_maf05.miss90.het0.7.LD505.1._scaf1.vcf --plink --out filt_maf05.miss90.het0.7.LD505.1_N

/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --file filt_maf05.miss90.het0.7.LD505.1_N --make-bed --out filt_maf05.miss90.het0.7.LD505.1_N --noweb
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --noweb --bfile filt_maf05.miss90.het0.7.LD505.1_N --indep-pairwise 3 10 .1
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --noweb --bfile filt_maf05.miss90.het0.7.LD505.1_N --extract plink.prune.in 
/nfs/shares/hodgins-grp/bin/gcta/gcta64  --bfile filt_maf05.miss90.het0.7.LD505.1_N --make-grm --extract plink.prune.in --out filt_maf05.miss90.het0.7.LD505.1_N_sub
/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss90.het0.7.LD505.1 --grm filt_maf05.miss90.het0.7.LD505.1_N_sub --pheno pheno --out filt_maf05.miss90.het0.7.LD505.1


#####





/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --file filt_maf05.miss90.het0.7 --allow-no-sex --make-bed --out filt_maf05.miss90.het0.7 --noweb

-filter based on linkage: command  that specifies 50 5 0.5 would a) consider a window of 50 SNPs, b) calculate LD between each pair of SNPs in the window, b) remove one of a pair of SNPs if the LD is greater than 0.5, c) shift the window 5 SNPs forward and repeat the procedure.
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --noweb --bfile filt_maf05.miss90.het0.7 --indep-pairwise 50 5 .5
- Note that filtering with r2=0.2, there were 473 SNPs remaining, with 0.5 there are 687


/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --noweb --allow-no-sex  --bfile filt_maf05.miss90.het0.7 --extract plink.prune.in 


#check if more than 22 chromosomes are allowed
/nfs/shares/hodgins-grp/bin/gcta/gcta64  --bfile filt_maf05.miss90.het0.7 --make-grm --extract plink.prune.in  --autosome --out filt_maf05.miss90.het0.7_sub


/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile filt_maf05.miss90.het0.7 --grm filt_maf05.miss90.het0.7_sub --pheno pheno --out filt_maf05.miss90.het0.7





--rename-chrs file
rename chromosomes according to the map in file, with "old_name new_name\n" pairs separated by whitespaces, each on a separate line.



for f in `cat prefix2.list`
do
vcftools --vcf $f\_VCFallsim1_N.vcf --plink --out $f\_N
done
```


Genotype-Phenotype Associations (GPA) with control for population structure
was conducted using GCTA using the mixed linear model based association
analysis (MLMA) module. Population structure was adjusted by incorporating both the
genetic relationship matrix (estimated with GCTA) as a random effect and the first 10
MDS components as covariates.

for f in `cat prefix2.list`
do
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --file $f\_N --make-bed --out $f\_N --noweb
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --noweb --bfile $f\_N --indep-pairwise 3 10 .1
/nfs/shares/hodgins-grp/bin/plink-1.07-x86_64/plink --noweb --bfile $f\_N --extract plink.prune.in 
/nfs/shares/hodgins-grp/bin/gcta/gcta64  --bfile $f\_N --make-grm --extract plink.prune.in --out $f\_N_sub
paste <(awk '{print $1"\t"$2}' $f.pheno) <(grep -v ^V $f\_MDS)  > $f\_MDS.mds



/nfs/shares/hodgins-grp/bin/gcta/gcta64 --mlma --bfile $f --grm $f\_N_sub --qcovar $f\_MDS.mds --pheno $f.pheno --out $f
done
```

and resulting pheno like:
 i0      i0      0.826778
i1      i1      0.504232
i2      i2      0.194862
i3      i3      -0.0772802
i4      i4      0.0257363
i5      i5      0.317192
i6      i6      -0.845423
i7      i7      0.194289
i8      i8      -0.276146
i9      i9      -0.797334

_MDS is a mulitdimensional scaling file direct from R, 
paste <(awk '{print $1"\t"$2}' $f.pheno) <(grep -v ^V $f\_MDS)  > $f\_MDS.mds
is pasting the population(?) data with the MDS data, line numbers in each file correspond
I guess the MDS file is optional, or maybe this is the kinship matrix?




