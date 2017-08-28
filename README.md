# rnaseq_pipeline

Mark Ziemann 2017

mark.ziemann@gmail.com

This script is an attempt to authomate some of the steps we routinely use in RNA-seq analysis including:
- Download fastq.gz files from service provider
- Perform quality trimming
- Download and index the genome for analysis with STAR mapper
- Map with STAR and generate the counts file - no SAM or BAM files are generated
- Enbedded R script to read in the counts files and perform differential expression analysis
- Pathway analysis using GSEA with any number of GMT files

The script requires the user to specify several things to run successfully
- PROJ_NAME this is the name of the experiment, ie experiment_123. it will be created by the script
- USER this is the username for downloading the fastq.gz files from the service provider
- PW this is the password for downloading the fastq.gz files from the service provider 
- URL this is the location of the data FOLDER on service providers FTP server
- GENOME_URL this is the URL to download Ensembl genome to use as a reference
- GTF_URL this is the URL to download the Ensembl gene annotation set
- STRAND this is the strandedness of the library preparation 0=unstranded, 1=pos strand, 2=neg strand
- GSEA_UTILS this is the folder which contains the GSEA JAR file and the GMT files to use in pathway analysis

Other things to remember before kicking it off
- Dependencies: Skewer, STAR, GNU parallel, R (limma/edgeR/parallel/reshape/plyr/statmod/locfit/
- Sample name specifications
   - No spaces or non alphanumeric characters
   - The first string before the underscore "_" is the sample group
   - The second string is the replicate number
   - Text after the second underscore doesn't matter
   - Follow this schema: TreatmentGroup_replicate1_ANYTHING_XX43437841



