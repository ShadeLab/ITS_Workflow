#Analysis of ITS_apple_replant_45samples 

```
########### GLBRC/Bonito's Lab ITS Workflow #######################
```
```
# all the analyses results are stored in /mnt/research/ShadeLab/WorkingSpace/Bintarti/ITS_apple_replant/ITS_45_aprep/Workflow1_ITS1_denovoClust/
# ITS raw sequence data stored on HPCC 
# the analyses were conducted in new centos7 HPCC system (ssh dev-intel18)
# quality filtering was conducted using USEARCH10
# primer/adapter removal was performed using cutadapt 1.17 with Python 2.7.15.
# trimming was conducted with length of 200 bp
# ZOTU table was built using denoising 
# OTU table was built using de novo clustering with 97% of similarity
# Classifying was conducted using "consensus_taxonomy.txt" file produced by CONSTAX
# There is no need to filter mitochondria, chloroplast, or protist because CONSTAX has been modified
# further analyses were conducted using QIIME 1.9.1 pipeline (you have to install it in your home directory because it cannot be loaded either in old centos6 or new centos7 HPCC system)
# unfortunately, taxonomy file (consensus_taxonomy.txt) produced from CONSTAX is biom unfriendly so it cannot be added into OTU table biom file using command: biom add-metadata
# add taxonomy to OTU table can be conducted/modified using R software, Excel, and Phyloseq 

apple replant soil ITS:  45 total files
Copied apple replant ITS sequences to working space for analysis
ShadeLab/WorkingSpace/Bintarti/ITS_apple_replant/ITS_45_aprep/raw_reads/
```
## QUALITY CHECKING AND PRE-FILTERING

a) count read numbers
```
for fastq in ../rawreads/*.fastq
do wc -l $fastq
done > reads_raw.counts
```
b) produce reads quality graphs using FastQC
```
mkdir stats

cat ../rawreads/*R1_001.fastq > raw_reads_R1.fastq; cat ../rawreads/*R2_001.fastq > raw_reads_R2.fastq

/mnt/home/bintarti/FastQC/fastqc raw_reads_R1.fastq raw_reads_R2.fastq -o stats && rm -rf raw_reads_R1.fastq raw_reads_R2.fastq
```

##MERGE PAIRED READS

```
mkdir mergedfastq

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs ../rawreads/*R1*.fastq -reverse ../rawreads/*R2*.fastq -fastq_maxdiffs 10 -fastq_minmergelen 50 -relabel @ -fastqout mergedfastq/merged2.fastq
```
###output
```
6818052  Pairs (6.8M)
   4601240  Merged (4.6M, 67.49%)
   2563437  Alignments with zero diffs (37.60%)
   1854285  Too many diffs (> 10) (27.20%)
    234054  No alignment found (3.43%)
         0  Alignment too short (< 16) (0.00%)
    128473  Merged too short (< 50)
    578828  Staggered pairs (8.49%) merged & trimmed
    202.81  Mean alignment length
    281.03  Mean merged length
      1.04  Mean fwd expected errors
      0.85  Mean rev expected errors
      0.23  Mean merged expected errors
```
##CHECK SEQUENCE QUALITY OF MERGED SEQS USING USEARCH AND VSEARCH
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 mergedfastq/merged2.fastq -output stats_eestats2_USEARCH.txt

/mnt/home/bintarti/vsearch-2.8.4-linux-x86_64/bin/vsearch -fastq_stats mergedfastq/merged2.fastq -log stats_results_VSEARCH.txt
```
###output
```
4601240 reads, max len 484, avg 281.0

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    4555392( 99.0%)    4597024( 99.9%)    4601216(100.0%)
   100    4272714( 92.9%)    4421917( 96.1%)    4440687( 96.5%)
   150    4117291( 89.5%)    4371651( 95.0%)    4438364( 96.5%)
   200    3944832( 85.7%)    4245254( 92.3%)    4360838( 94.8%)
   250    3770636( 81.9%)    4091828( 88.9%)    4238271( 92.1%)
   300     677838( 14.7%)     794362( 17.3%)     867429( 18.9%)
   350     173782(  3.8%)     218912(  4.8%)     256250(  5.6%)
   400      15972(  0.3%)      22710(  0.5%)      29544(  0.6%)
   450       1319(  0.0%)       2299(  0.0%)       3547(  0.1%)
```
## REMOVE PRIMER AND ADAPTERS WITH CUTADAPT, CHECK THE STATS
```
#upgrade pip
pip install --upgrade pip

#install latest version of cutadapt 1.17
pip install --user cutadapt
```
```
CS1-ITS1 (fwd): 5’- CTTGGTCATTTAGAGGAAGTAA – 3’ (EMP/Smith and Peay 2014)
CS2-ITS2 (rev): 5’- GCTGCGTTCTTCATCGATGC – 3’ (EMP/Smith and Peay 2014)
Reverse complement of adapter-reverse primer: GCATCGATGAAGAACGCAGC

/mnt/home/bintarti/.local/bin/cutadapt -g CTTGGTCATTTAGAGGAAGTAA -a GCATCGATGAAGAACGCAGC -f fastq -n 2 --discard-untrimmed --match-read-wildcards -o cut_merged.fastq mergedfastq/merged2.fastq > cut_adpt_results.txt

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 cut_merged.fastq -output cutdapt_eestats2_USEARCH.txt

#output: 4601240 reads, max len 466, avg 239.1
Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    4385216( 95.3%)    4436883( 96.4%)    4441885( 96.5%)
   100    4244089( 92.2%)    4412022( 95.9%)    4439955( 96.5%)
   150    4047081( 88.0%)    4290590( 93.2%)    4368298( 94.9%)
   200    3862899( 84.0%)    4139471( 90.0%)    4257802( 92.5%)
   250    1450615( 31.5%)    1607192( 34.9%)    1694306( 36.8%)
   300     203609(  4.4%)     249182(  5.4%)     286902(  6.2%)
   350      39578(  0.9%)      52640(  1.1%)      65240(  1.4%)
   400       3164(  0.1%)       4947(  0.1%)       7113(  0.2%)
   450         52(  0.0%)         78(  0.0%)        102(  0.0%)
```

## FILTERING, TRIMMING, QUALITY CHECK
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_filter cut_merged.fastq -fastq_maxee 1 -fastq_trunclen 200 -fastq_maxns 0 -fastaout filtered_cut_merged.fa -fastqout filtered_cut_merged.fastq

###output:
   4601240  Reads (4.6M)                    
    328219  Discarded reads length < 200
    133550  Discarded reads with expected errs > 1.00
   4139471  Filtered reads (4.1M, 90.0%)

#quality check
/mnt/home/bintarti/FastQC/fastqc filtered_cut_merged.fastq
```
##CLUSTERING AND DENOISING (OTUs vs. ESV)

###generating Exact Sequence Variants (ESV) also called 0-radius OTUs
```
1) DEREPLICATE

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques filtered_cut_merged.fastq -fastaout derep_filtered_cut_merged.fasta -sizeout

#output: 4139471 seqs, 548090 uniques, 373069 singletons (68.1%)
00:11 3.0Gb  Min size 1, median 1, max 401391, avg 7.55

2) DENOISING-MAKING ZOTU

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -unoise3 derep_filtered_cut_merged.fasta -tabbedout unoise_zotus.txt -zotus zotus.fasta

#output:
00:02 190Mb   100.0% 3018 amplicons, 1097779 bad (size >= 8)
00:39 202Mb   100.0% 3013 good, 5 chimeras                  
00:39 202Mb   100.0% Writing zotus

/mnt/home/bintarti/python_scripts-master/fasta_number.py zotus.fasta OTU_ > zotus_numbered.fasta

3) CONSTRUCT OTU TABLE-MAP READS BACK TO zotus_numbered.fasta AS DATABASE FILE

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global mergedfastq/merged2.fastq -db zotus_numbered.fasta -strand plus -id 0.97 -otutabout otu_table_ITS_UNOISE.txt

#OUTPUT:4250658 / 4601240 mapped to OTUs (92.4%) 
```
### CLUSTERING-MAKING OTU TABLE at 97% sequence similarity (de novo otu-picking)-UPARSE
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_otus derep_filtered_cut_merged.fasta -minsize 2 -otus otus.fasta -uparseout uparse_otus.txt -relabel OTU_

#output: 100.0% 3262 OTUs, 576 chimeras

## CONSTRUCT OTU TABLE-MAP READS BACK TO otus.fasta AS DATABASE FILE

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global mergedfastq/merged2.fastq -db otus.fasta -strand plus -id 0.97 -otutabout otu_table_ITS_UPARSE.txt

#output: 93.1% matched
4282789 / 4601240 mapped to OTUs (93.1%)
```
## TAXONOMIC CLASSIFICATION USING CONSTAX (or eventually RDP)
```
please refer to how "Running CONSTAX on the MSU HPCC on lab guru : https://my.labguru.com/knowledge/documents/330

# output directories: outputs, taxonomy_assignments, training_files
# move the output directories to: /mnt/research/ShadeLab/WorkingSpace/Bintarti/ITS_apple_replant/ITS_45_aprep/BonitoLab_workflow/
# Use the taxonomy file : outputs/consensus_taxonomy.txt
```

## CONVERT .TXT FILE INTO .BIOM FILE - CONSTAX output is not BIOM friendly, thus I can't add using: biom add_metadata
```
1.) FOR: otu_table_ITS_UPARSE.tx

biom convert -i otu_table_ITS_UPARSE.txt -o otu_table_ITS_UPARSE.biom --table-type="OTU table" --to-json

#output file: otu_table_ITS_UPARSE.biom

2.) FOR: otu_table_ITS_UNOISE.txt

biom convert -i otu_table_ITS_UNOISE.txt -o otu_table_ITS_UNOISE.biom --table-type="OTU table" --to-json

#output file: otu_table_ITS_UNOISE.biom
```

## ALIGN SEQUENCES USING MUSCLE
```
module load MUSCLE/3.8.31

1.) otus.fasta

muscle -in otus.fasta -out OTUS_aligned.fasta -maxiters 2 -diags1

2.) zotus_numbered.fasta

muscle -in zotus_numbered.fasta -out ZOTUS_aligned.fasta -maxiters 2 -diags1
```

# START THE ANALYSES USING QIIME 1.9.1 PIPELINE

```
# should install QIIME 1.9.1 using miniconda 2 (Python2.7) environment
# download miniconda 2 and put it in home directory
```
##install miniconda 
```
bash Miniconda2-latest-Linux-x86_64
```
##create qiime1 environment and install qiime
```
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
```
##activate qiime1 environment and test qiime installation
```
source activate qiime1
print_qiime_config.py -t
```
## FILTER EXCESS GAPS FROM ALIGNMENT USING QIIME 1.9.1
```
filter_alignment.py -i OTUS_aligned.fasta -o OTUS_filtered_alignment

filter_alignment.py -i  ZOTUS_aligned.fasta -o ZOTUS_filtered_alignment
```
## MAKE PHYLOGENY WITH FASTREE
```
make_phylogeny.py -i OTUS_filtered_alignment/OTUS_aligned_pfiltered.fasta -o OTUS_rep_set.tre

make_phylogeny.py -i ZOTUS_filtered_alignment/ZOTUS_aligned_pfiltered.fasta -o ZOTUS_rep_set.tre
```
## Rarefy OTU table to lowest sequencing depth
```
# summarize OTU and ZOTU table

biom summarize-table -i otu_table_ITS_UPARSE.biom -o OTU_table_sum.txt
biom summarize-table -i otu_table_ITS_UNOISE.biom -o ZOTU_table_sum.txt

# rarefaction

single_rarefaction.py -d 57144 -o   rarefied/OTU_single_rare.biom -i otu_table_ITS_UPARSE.biom
single_rarefaction.py -d 56581 -o   rarefied/ZOTU_single_rare.biom -i otu_table_ITS_UNOISE.biom

# summarize rarefied OTU and ZOTU table

biom summarize-table -i OTU_single_rare.biom -o OTU_single_rare_sum.txt
biom summarize-table -i ZOTU_single_rare.biom -o ZOTU_single_rare_sum.txt


# Convert OTU table file: .biom into .txt

biom convert -i OTU_single_rare.biom -o OTU_single_rare.txt --to-tsv
biom convert -i ZOTU_single_rare.biom -o ZOTU_single_rare.txt --to-tsv
```

