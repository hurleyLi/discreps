# Code for Manuscript:
> Exome Variant Discrepancies due to Reference Genome Differences (*AJHG*; in press)

In this study, we determined variant call discrepancies between GRCh37 and GRCh38 using exome sequencing data (n = 1,572) from the Baylor-Hopkins Center for Mendelian Genomics. The code generated from this study are provided here. 

## Variant filtering and lift-over
`filter_and_liftover.submit.sh`

* This script is optimized for submitting jobs on the computing clusters of Human Genome Sequencing Center (HGSC) at the Baylor College of Medicine. However, the code should allow general usage on single machines with minimum modifications.  

* The GRCh37 and GRCh38 reference files are too big to be included on GitHub. We provide download links here:
    * [hs37d5](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz)
    * [GRCh38](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa)

* The input files should be obtained from [dbGaP](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000711.v7.p2) (phs000711.v7.p2)

## Identification of DISCordant REference Patches (DISCREPs)
`discreps.py`

* The analyses were performed on GRCh37 and GRCh38 separately. 
* The genome was divided into 10kb windows. In each window, the total number of distinct variants across all samples was counted. We only kept windows with >10 distinct variants for analysis. 
* In each genomic window, one-sided Fisher’s exact test was used to compare concordant vs. discordant variant counts again the baseline level across the whole exome.  

## Enrichment of genomic features in DISCREPs
`LOLA.R`

* The following genomic elements were tested for enrichment within DISCREPs using LOLA:
`simple tandem repeats`, `microsatellite`, `segmental duplications`, `interrupted repeats`, `known assembly problems`, `loci with fix patches`, `loci with alternate haplotypes`, `loci with known genome assembly differences`, and `gaps in the assembly`
* For DISCREPs regions that overlapped between GRCh37 and GRCh38, counts of variants within each genomic windows were combined from both GRCh37 and GRCh38.
* Location of the genomic features were downloaded from UCSC genome browser.

## Identification of genes influenced by the reference assembly
`discordant_genes.py`

* Genes with GENCODE annotations on both the GRCh37 and GRCh38 references were analyzed. 
* Variants were found in a total of 19,003 genes and filtered to keep genes with at least one discordant variant call 
