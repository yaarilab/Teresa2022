The pipeline process TCR that were produce in the same fashion as those in [Teresa Rubio et al. 2020](https://www.sciencedirect.com/science/article/pii/S2667119022000040)

Library preperation and sequencing method:

Whole CD4 T cell RNA was extracted using the miRNeasy Micro Kit, then sequenced on an Illumina HiSeq2500 with HiSeq Sequencing v4 Chemistry. RNA-seq employed ultra-low input library prep with polyA selection, utilizing strand-specificity and 125 bp paired-end reads at a depth of 50 million for short-read sequencing.


Input files:

1. The read 
2. csv file with information of the sample

To test the pipeline:

you can download the read from [github](https://github.com/ConesaLab/TCR_nextflow/tree/main/data/reads) and upload directly to dolphinnext.


Output files:
1. report for each steps
2. global report from the last process


Pipeline container:

* Docker: ssnnm/mhecd4tcr:0.1.0
	
* For the first step(mixcr) - milaboratory/mixcr:latest (the lisence in need to given in the RunOptions)


Sequence processing steps:

1. mixcr_analyze 
1.1 mixcr_qc 
2. data_filtering
3. dataset_overview
4. correlations 
5. overlap
6. diversity
7. kmers
8. network 
9. ddbb
10. report


Files used:

* [sample data file](https://github.com/ConesaLab/TCR_nextflow/blob/main/sampleslist.csv)

