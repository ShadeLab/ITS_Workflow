# Analysis of ITS_apple_replant_45samples

```
######## Modified from Dr. Patrick Kearns Workflow for ITS Sequences###########
```
```
# all the analyses results are stored in /mnt/research/ShadeLab/WorkingSpace/Bintarti/ITS_apple_replant/ITS_45_aprep/Workflow3_ITS1_OpenRef/
# ITS raw sequence data stored on HPCC 
# the analyses were conducted in old centos7 HPCC system
# quality filtering was conducted using USEARCH10 and FASTX 0.0.14
# primer/adapter removal was performed using cutadapt 1.17 with Python 2.7.15
# no trimming was performed
# OTU table was built using both open reference of UNITE V.7.2 and de novo clustering with 97% of similarity
# Classifying was conducted using "consensus_taxonomy.txt" output file produced by CONSTAX
# There is no need to filter mitochondria, chloroplast, or protist because CONSTAX has been modified
# further analyses were conducted using QIIME 1.9.1 pipeline (you have to download and install it (using miniconda2 environment) in your home directory because this version does not exist either in old centos6 or new centos7 HPCC system)
# unfortunately, taxonomy file (consensus_taxonomy.txt) produced from CONSTAX is biom unfriendly so it cannot be added into OTU table biom file using command: biom add-metadata
# add taxonomy to OTU table can be conducted/modified using R software, Excel, and Phyloseq 

apple replant soil ITS:  45 total files
Copied apple replant ITS sequences to working space for analysis
ShadeLab/WorkingSpace/Bintarti/ITS_apple_replant/ITS_45_aprep/rawreads/
```

## MERGE PAIRED END READS
```
#decompress the reads
gunzip *.gz

mkdir mergedfastq

#will use a more variable sequence size 200-500bp due to the variability of ITS amplicons.
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs ../rawreads/*R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastqout mergedfastq/merged.fq -fastq_merge_maxee 1.0 -fastq_minmergelen 200 -fastq_maxmergelen 500

###output:

6818052  Pairs (6.8M)
   4373305  Merged (4.4M, 64.14%)
   2563437  Alignments with zero diffs (37.60%)
   1854285  Too many diffs (> 10) (27.20%)
    234054  No alignment found (3.43%)
         0  Alignment too short (< 16) (0.00%)
    356408  Merged too short (< 200)
         0  Merged too long (> 500)
         0  Exp.errs. too high (max=1.0) (0.00%)
    578828  Staggered pairs (8.49%) merged & trimmed
    208.28  Mean alignment length
    290.58  Mean merged length
      0.90  Mean fwd expected errors
      0.69  Mean rev expected errors
      0.24  Mean merged expected errors
```
## CHECK SEQUENCE QUALITY OF MERGED SEQS
```
mkdir fastq_info

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 mergedfastq/merged.fq -output fastq_info/eestats.txt

##output:
Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    4334401( 99.1%)    4369982( 99.9%)    4373284(100.0%)
   100    4205653( 96.2%)    4353143( 99.5%)    4371680(100.0%)
   150    4052944( 92.7%)    4304023( 98.4%)    4369440( 99.9%)
   200    3944832( 90.2%)    4245254( 97.1%)    4360838( 99.7%)
   250    3770636( 86.2%)    4091828( 93.6%)    4238271( 96.9%)
   300     677838( 15.5%)     794362( 18.2%)     867429( 19.8%)
   350     173782(  4.0%)     218912(  5.0%)     256250(  5.9%)
   400      15972(  0.4%)      22710(  0.5%)      29544(  0.7%)
   450       1319(  0.0%)       2299(  0.1%)       3547(  0.1%)
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

/mnt/home/bintarti/.local/bin/cutadapt -g CTTGGTCATTTAGAGGAAGTAA -a GCATCGATGAAGAACGCAGC -f fastq -n 2 -m 20 --discard-untrimmed --match-read-wildcards -o cut_merged.fastq mergedfastq/merged.fq > cut_adpt_results.txt

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 cut_merged.fastq -output cutdapt_eestats2_USEARCH.txt

#output: 4372177 reads, max len 463, avg 248.6

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    4315745( 98.7%)    4366816( 99.9%)    4371786(100.0%)
   100    4178074( 95.6%)    4342928( 99.3%)    4369905( 99.9%)
   150    4045464( 92.5%)    4288545( 98.1%)    4365967( 99.9%)
   200    3862516( 88.3%)    4138727( 94.7%)    4256775( 97.4%)
   250    1450277( 33.2%)    1606508( 36.7%)    1693344( 38.7%)
   300     203478(  4.7%)     248866(  5.7%)     286097(  6.5%)
   350      39477(  0.9%)      52476(  1.2%)      64854(  1.5%)
   400       3066(  0.1%)       4801(  0.1%)       6918(  0.2%)
   450          0(  0.0%)          1(  0.0%)          1(  0.0%)
```
## QUALITY FILTER FASTQ SEQUENCES USING FASTX toolkit (v. 0.0.14)-ONLY CAN BE LOADED IN CENTOS6 HPCC SYSTEM
```
# module load FASTX/0.0.14
fastq_quality_filter -q 30 -i cut_merged.fastq -o filtered_cut_merged.fq -p 50
```
## DEREPLICATE
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques filtered_cut_merged.fq -fastaout derep_filtered_cut_merged.fasta -sizeout

#output: 4372031 seqs, 877374 uniques, 645307 singletons (73.5%)
```
## Remove Singeltons
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sortbysize derep_filtered_cut_merged.fasta -fastaout nosigs_derep_filtered_cut_merged.fasta -minsize 2

#output: Sorting 232067 sequences
```
## Precluster Sequences
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_fast nosigs_derep_filtered_cut_merged.fasta -centroids denoised_nosigs_derep_filtered_cut_merged.fasta -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size

#output: Seqs  232067 (232.1k)
  Clusters  9697
  Max size  626378 (626.4k)
  Avg size  384.3
  Min size  2
Singletons  0, 0.0% of seqs, 0.0% of clusters
```
## Reference-based OTU picking against UNITE fungal ITS database (v. 7.2) at 97% sequence similarity
```
 /mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global denoised_nosigs_derep_filtered_cut_merged.fasta -id 0.97 -db /mnt/research/ShadeLab/UNITE_v7.2/sh_refs_qiime_ver7_97_s_01.12.2017.fasta -strand plus -uc ref_seqs.uc -dbmatched UNITE_reference.fasta -notmatched UNITE_failed_reference.fq
```
## Sort by size and then de novo OTU picking on sequences that failed to hit UNITE
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64  -sortbysize UNITE_failed_reference.fq -fastaout sorted_UNITE_failed_reference.fq

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64  -cluster_otus sorted_UNITE_failed_reference.fq -minsize 2 -otus denovo_otus.fasta -relabel OTU_dn_ -uparseout denovo_out.up

#output: 2740 OTUs, 515 chimeras
```
## Combine the rep sets between de novo and reference-based OTU picking
```
cat UNITE_reference.fasta denovo_otus.fasta > FULL_REP_SET.fna

/mnt/home/bintarti/python_scripts-master/fasta_number.py FULL_REP_SET.fna OTU_ > NUMBERED_FULL_REP_SET.fasta
```
## Map rep_set back to pre-dereplicated sequences and make OTU tables
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global mergedfastq/merged.fq -db NUMBERED_FULL_REP_SET.fasta  -strand plus -id 0.97 -uc OTU_map.uc -otutabout OPEN_REF_OTU_TABLE_ITS.txt

#output: 4101281 / 4373305 mapped to OTUs (93.8%)
#count OTUs
grep -c '>' NUMBERED_FULL_REP_SET.fasta
#OTUs: 3719
```
## TAXONOMIC CLASSIFICATION USING CONSTAX
```
please refer to how "Running CONSTAX on the MSU HPCC on lab guru : https://my.labguru.com/knowledge/documents/330

# YOU HAVE TO RUN CONSTAX WITHIN THE OLD CENTOS6 HPCC SYSTEM!
# add "NUMBERED_FULL_REP_SET.fasta" into: mnt/home/bintarti/CONSTAX_hpcc/otus
# change the PATH in: mnt/home/bintarti/CONSTAX_hpcc/config
# output directories: outputs, taxonomy_assignments, training_files
# move the output directories to: /mnt/research/ShadeLab/WorkingSpace/Bintarti/ITS_apple_replant/ITS_45_aprep/Workflow3_ITS1_OpenRef/CONSTAX_OTU
# Use the taxonomy file : outputs/consensus_taxonomy.txt
```
## CONVERT .TXT FILE INTO .BIOM FILE - CONSTAX output is not BIOM friendly, thus do not use this command: biom add_metadata
```
biom convert -i OPEN_REF_OTU_TABLE_ITS.txt -o OPEN_REF_OTU_TABLE_ITS.biom --table-type="OTU table" --to-json

#output file: OPEN_REF_OTU_TABLE_ITS.biom
```
## ALIGN SEQUENCES USING MUSCLE (in CENTOS 6)
```
module load MUSCLE/3.8.31
muscle -in NUMBERED_FULL_REP_SET.fasta -out OTUS_aligned.fasta -maxiters 2 -diags1
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
```
## MAKE PHYLOGENY WITH FASTREE
```
make_phylogeny.py -i OTUS_filtered_alignment/OTUS_aligned_pfiltered.fasta -o OTUS_rep_set.tre
```
## Rarefy OTU table to lowest sequencing depth
```
# summarize OTU table
biom summarize-table -i OPEN_REF_OTU_TABLE_ITS.biom -o OTU_table_sum.txt

# rarefaction
single_rarefaction.py -d 54825 -o rarefied/OTU_single_rare.biom -i OPEN_REF_OTU_TABLE_ITS.biom

# summarize rarefied OTU table
biom summarize-table -i rarefied/OTU_single_rare.biom -o rarefied/OTU_single_rare_sum.txt

# Convert rarefied OTU table file: .biom into .txt for further analyses on R
biom convert -i rarefied/OTU_single_rare.biom -o rarefied/OTU_single_rare.txt --to-tsv
```