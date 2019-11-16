R-scripts for phenotype-environment-allele associations

Sequence data are available at the National Center for Biotechnology Information (NCBI) Sequence Read Archive under Bioproject PRJNA449949.Processed genomic data is available on https://doi.org/10.6084/m9.figshare.9783311.v1. Phenotypic, environmental and population structure data is available on https://doi.org/10.6084/m9.figshare.7217744.v1

Citation: under review


Files in this repo

0_BayEnv.md: BayEnv2 walkthrough

0_GCTA.md: GCTA walkthrough

0_GCTA_CnP.txt: GCTA code only

1_dataprepare.Rmd: Data reformatting from BayEnv2 output

2a_xtx_wind.Rmd: window analysis XTX

2b_bf_windowtest.Rmd: window analysis BF

2c_gcta.Rmd: window analysis GCTA

3a_xtxwind_conv.Rmd: binomial outlier window analysis xtx

3b_bfwind_conv.Rmd: binomial outlier window analysis bf

4a_null-wxtx.Rmd: null-W analysis xtx

4b_null-w_bf.Rmd: null-W analysis BF

4c_GCTA_nullw.Rmd: null-W analysis GCTA

5a_sum_parallel.Rmd: summarise the parallel results among ranges for xtx

5b_sum_parallelbf.Rmd: summarise the parallel results among ranges for BF

The outputs of 5a and 5b have been engineered in an excel book, where significant results (p<0.005) from the null-W were marked as 'yes'

6_overlap_gctabfxtx.Rmd: detect overlap among the xtx, bf and gcta windows

7_LDconvertSNPtable2LD.Rmd: convert a SNP table input to an LD matrix

8_LD_rand.Rmd: select random SNPs to test in windows

9_LD_process.Rmd: analyse the LD results