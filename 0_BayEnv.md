BAYENV 
===


**This script requires climate data. Produce this data with  
QGIS_coordinate_climatevalue_Walkthrough.md**

BAYENV:  
1) Variance covariance matrix with subset, multiple runs  
2) Check convergence MCMC runs/correlation matrices (cov2cor in R) # not in this walkthrough 
3) XTX on subset —> identify outliers  
4) remove outliers original subset and repeat 1&2  
5) Env/xtx  on complete dataset (using corrected variance/covariance matrix)  

-------------

**A) Transform SNPtable to Bayenv format**
===

To transform the custom-filtered SNPtable to Bayenv(2) format, my choice was to use PGDspider to change the input file for Structure.
http://www.cmpg.unibe.ch/software/PGDSpider/. 

One thing to note though is that you might want different input files for Structure and Bayenv. For Structure, I applied very
stringent filtering, with <50% missing data per locus. For Bayenv though I allowed for
missing data up to 90% per locus, as I need a dataset more representative of the entire
dataset (cite, Coop?). Make sure to put in the correct format for Structure in PGDspider.
The script in GBS6_STRUCTURE will output a file with diploid data on 2 consecutive rows,
no phase information, missing value code of -9, SNP format, inclusion of locus, individual
names and population identifier. There is no additional information in the file.

Additionally, there is a functionality in PGDSpider to transform vcf directly to Bayenv

-------------

**B)    Estimate covariance matrix**
===


First, we need to make the covariance matrix. This will take a long time. For 90k SNPs over 30 pops it ran ~100k simulations per day

To run Bayenv2, type in command line below command (-i is input file, -p is number of populations, -k is
number of iterations, -r is the root, this can be any number). This will print to "matrices" the
updated covariance matrix for every 5000 iterations:

``` bayenv2 -i SNPSFILE -p 89 -k 100000 -r 63479 > matrices ```

This CANNON be run as an array job, each sub job would simply overwrite output of another job.


```
#!/bin/sh
### job shell -S
#$ -S /bin/sh
### keep yours files together in the same directory
#$ -cwd
### descriptive name of job, avoid spaces
#$ -N covau2
### walltime hh:mm:ss per subjob
#$ -l h_rt=100:00:00
### memory requirements
#$ -l h_vmem=10G
#$ -o covau2.log 
#$ -e covau2.err
#$ -pe smp 4
#$ -M lotte.van.boheemen@monash.edu
#$ -m abes

. /etc/profile
module load bayenv2
bayenv2 -i au42_48.bay.txt  -p 7  -k 500000 -r 52681 > matrix_aurun2.out
```

--------------

**C)    Run Bayenv with env file**
===

This step will require 3 input files, ENV, SNPFILE and MATRIX_FILE. The latter 2 look the same as
files produced/used above, but they are not! The code to run this is as follows (don't run yet, I am
going to explain every step):

``` bayenv2 -i SNPFILE -m matrix_file -e env -p 4 -k 500000 -n 8 -t -r 429 ```

-p is the number of populations, -k is the number if iterations (needs to be larger then 1000 or it
will divide by 0 in the bayenv source code), -n is the number of environmental variables, -t
indicates 'rest' mode, and will make sure the Bayes Factor is calculated for every SNP, -r is the
random root.

1]   Environmental variables are standardised by substracting the mean and dividing through the
standard deviation. Even though I will be analysing subsets for each of the continents and all of
the samples, I won't subset the environmental variables and standardize using data within each
range. If my goal is to look at parellelism between introduced ranges, I need to hold the gradients
similar. I cannot compare if standardisation is different.

The format of the environmental data for 4 populations and 8 environmental factors is shown below.
Each column is a population, in the same order as the SNPSFILE and input file.

-2.053517756	-1.905325197	-1.828770599	-1.836090653
1.563096989	1.594375406	1.600992044	1.574604643
-1.181446926	-1.171381351	-0.980135416	0.378717279
-1.465018067	-1.47330868	-1.550562123	-0.793855224
0.614591847	0.56667764	1.892304032	-0.455492108
-0.572885276	0.313220562	0.060047465	0.249927288
0.482329177	1.521168194	2.759319348	-0.571609245
0.872405436	1.162237909	2.359372036	-0.627162575


2]  The matrix_file is a summary or excerpt from above "matrices". Either average among all computed
covariance matrices (not done here, needs evaluation of convergence), or use the last one. 
Note that bayenv is VERY sensitive in the use of whitespaces vs tabs, and a simple copy 
paste might change tabs into whitespaces. To replace
multiple spaces with a single tab, use the following:

``` cat spaced-file | sed 's/ \+/\t/g' > tabbed-file ```

Another option is to yank the lines out if the matrix file. Go to the first line of the last matrix
<linenumber>gg and type the number of lines you want to copy (should be equal to number of populations in your file
<numberoflines>yy. Go to a new file and type p.    
If there are more then 50 lines, we need to adjust vim to store more in the registers. Type
```:set viminfo='20,<1000,s1000```
Followed by the yanking.

Open the new file and simply type <p>

Now, you want to average over multiple covariance matrices from diffferent runs. Paste all the
final matrices of each run into a new file and load R to get a mean matrix:

```
module load R
A<-read.table("matrix_sel2_allr")
B<-read.table("matrix_sel_allr")
C<-read.table("matrix_sel3_allr")
D<-read.table("matrix_sel4_allr")
E<-read.table("matrix_sel5_allr")

my.list<-list(A,B,C,D,E)
meanm<-Reduce("+", my.list) / length(my.list)
write.table(meanm,"mean_matrixallr",sep="\t",col.names=FALSE,row.names=FALSE)
```


3] SNPfiles contain the allele counts (on 2 lines for diploids) for ONE SNP ONLY. This means that
SNPSFILE needs to be cut into little pieces. Below script will do this. Sidenote on this script:
authorship is unknown, copied from shared folder, edits by Simon Michnowicz & Philip Chan

```	
#bayenv.pl
#!/usr/bin/perl

my $in = $ARGV[0]; #merged bsnp table loci (rows) pop (columns)
my $which = $ARGV[1];  # int that identifies which block of two lines to process
my $tmpdir = $ARGV[2];  # temporary directory to stash the lociXXXX files

my @array=();
print "Input file is $in\n";
open IN, $in;
while (<IN>){
   chomp;
   push (@array, $_);
#  print "$_\n";           
}
close IN;

my $len=$#array+1;
my $counter=1;
#print "Size of input file is $len lines\n";

my $fileIndex = $which;
my $which = ($which-1) * 2;

for(my $f=$which; $f < $len; $f++){
   print "Index is $f\n";
   if ($counter==1) {
      $counter=2;
      $fileIndex=$which;

      open OUT, ">$tmpdir/loci$fileIndex";
      my @temp=split(" ", $array[$f]);
      foreach (@temp){
          print OUT "$_\t";
      }
                print OUT "\n";
    }
    elsif ($counter == 2) {
      $counter=1;

      my @temp=split(" ", $array[$f]);
      foreach (@temp){
         print OUT "$_\t";
      }
      print OUT "\n";
      close OUT;

      #print "Fileanme=loci$fileIndex\n";

      my $command="bayenv2 -i $tmpdir/loci$fileIndex -m mean_matrixau -e env -p 7 -k 500000 -n 23 -t -X -r 666 -o bf_$fileIndex";
      print   "$command\n";
      system ("$command");
      #system "rm loci*";
      exit 0;
    }
}

```

Run this script as below. allr473_allsnp48.bay.txt is the SNPSFILE

``` perl bayenv.pl allr473_allsnp48.bay.txt  ```

This will produce files containing information for every locus, the SNPfiles. Note that the above
script is producing the right number of files in the right order, but the naming is off, as it skips
every even number.

Now, this script can become really terrible with many SNPs as it will spit out a file for every 2 lines in your infile. Therefore, the above script should be used together with below script and run on the cluster. Below job script is written by Philip Chan, and will create parent and child directories for the SNPs

```
#!/bin/sh

### job shell -S
#$ -S /bin/sh

### keep yours files together in the same directory
#$ -cwd

### descriptive name of job, avoid spaces
#$ -N BayEau

### walltime hh:mm:ss per subjob
#$ -l h_rt=200:00:00

### memory requirements
#$ -l h_vmem=12G

#$ -o bayEau.log 
#$ -e bayEau.err

#$ -l dpod=1
#$ -l passwd=lotteanv

#$ -pe smp 4

#$ -m n

#$ -t 1-30420

. /etc/profile
module load bayenv2

ID=$SGE_TASK_ID

P=`expr $ID / 1000`
C=`expr $ID % 1000`

PDIR=`printf "%02d" $P`
CDIR=`printf "%03d" $C`

echo "Working on $PDIR / $CDIR"

if [ ! -d $PDIR ]
then
  mkdir $PDIR
fi

cd $PDIR

if [ ! -d $CDIR ]
then
  mkdir $CDIR
  fi

cd $CDIR

IND=`expr $ID - 1`
INDEX=`expr $IND \* 2`

if [ -f bf_$INDEX.bf ]
then
  echo "Output file for $PDIR / $CDIR already present!"
  exit 0
fi

pwd

cp -p ../../env .
cp -p ../../*matrix* .

ulimit -s 100000
perl ../../bayenv.pl ../../bay_au $ID $TMPDIR > out 2> err
```

--------------

**D)	Running XtX matrix**
===

We can use the same script as above to run the XtX matrices. This also requires an environment input
file, but won't use it. The -t flag is left out, the -X flag is added. For this, we can also create a dummy environment file, filled with 0 and the amount of columns as populations

Note that this script will also run the environment test if a correct environment file is given
```
bayenv2 -i loci1a -m matrix_sel_allcra -e env -p 86 -k 100000 -n 1 -t -X -r 429
```

Or change in bayenv_split:
```
my $command="bayenv2 -i loci$fileIndex -m matrix -e env -p 7 -k 100000 -n 23 -t -X -r 429";
```

-------------

**E)	Pasting the results together**
===

Now all the results are in subfolders with one file per directory. We want to extract all this information and put it in the same file.
```
cat */*/*.bf > bf_result_allr
cat */*/*.xtx > xtx_result_allr
```

This needs to be done in every replicate run (with different roots), over which we can then average the values

```
module load R
A<-read.table("bf_result_allr_run1")
B<-read.table("bf_result_allr_run2")
C<-read.table("bf_result_allr_run3")

my.list<-list(A,B,C)
meanm<-Reduce("+", my.list) / length(my.list)
write.table(meanm,"bf_result_allr_mean",sep="\t",col.names=FALSE,row.names=FALSE)
```


-------

**Residual scripting**

In below command (for BAYENV), the first '0' indicates we want to estimate the covariance matrix,
'89' is the number of populations, '-23457' is the negative seed, '10000' is the number of
iterations (100,000 is standard for covariance matrix), 'snp48_bayenv' is the snpstable we want to
estimate the covariance matrix from. Note this is a random subset of non-linked SNPs, as described
in the structure walkthrough. According to the Bayenv2 manual, increasing the number of SNPs (e.g.
thousands) doesn't result in large differences in estimation of the covariance matrix. This SNPtable
is different from the SNP table we will use to do further calculations in bayenv, which will be more
inclusive.

``` bayenv 0 89 -23457 100000 snp48_bayenv > matrixall.out ```