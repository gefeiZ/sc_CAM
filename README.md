# sc_CAM

A statistical method for single cell Chromatin Accessibility QTL Mapping based on scATAC-seq data.


Introduction
------------

A python package for **discovering caQTL with single cell ATAC-seq data, without WGS or WES data**. 


Installation
------------

pip install sc_CAM

------------

### Dependencies

-   `samtools`
-   `featureCounts`
-   `GATK`
-   `picard`
-   `java`
-   `bam-readcount`
-   `Python3`
-   `Pandas`
  
------------

## Usage
------------
### Step0: 


#### Build index for mapping software STAR

#### Building Peak matrix, export barcodes for each cell types.

#### From Fastq files to Bam files; SNV calling from Bam files

Check  -   `GATK_usage.ipynb`


### Step1: Build SNV metadata pairs


Check  -   `usage.ipynb`

### Step 2: Conduct caQTL analysis 


Check  -   `usage.ipynb`

------------
Format of caQTL Result,

-   `SNVid`: Id of SNV, represented as CHR\_\_POS
-   `Peakid`: Peak 
-   `sample_size_Ref`, `sample_size_Alt`: number of cells in REF and ALT
    group
-   `theta_Ref`, `theta_Alt`, `mu_Ref`, `mu_Alt`, `size_Ref`,
    `size_Alt`, `prob_Ref`, `prob_Alt`: estimated parameters of ZINB for
    REF group and ALT group
-   `total_mean_Ref`, `total_mean_Alt`: mean of the REF group and ALT
    group
-   `foldChange`: fold change of mean of REF group (`total_mean_Ref`)
    with respect to mean of ALT group (`total_mean_Alt`)
-   `chi2LR1`: chi square statistic for the test
-   `pvalue`, `adjusted_pvalue`: pvalue and adjusted pvalue of the test.
    If adjusted p-value is smaller than some threshold, this SNV shows
    significant eQTL effect on the target gene


## Running sc_CAM on example data:

1.Download the sc_CAM from Pypi

2.Download the toy data from Zenodo

3.Follow the usage instruction


The structure should be like:
------------
    output_folder
    │   
    │
    |─── bam
    |   │─── cell_type1
    |   |   │  cell_1.bam
    |   |   │  cell_1.bam.bai
    |   |   │  cell_2.bam
    |   |   │  cell_2.bam.bai
    |   |   | ...
    |   │─── cell_type2
    |   |   │  cell_1.bam
    |   |   │  cell_1.bam.bai
    |   |   │  cell_2.bam
    |   |   │  cell_2.bam.bai
    |   |   | ...
    |─── snv
    |   |─── cell_level_snv
    |   |   |─── cell_type1
    |   |   │      cell_1_filtered_pass.vcf
    |   |   |      cell_1_filtered_pass.vcf.idx
    |   |   │      cell_2_filtered_pass.vcf
    |   |   |      cell_2_filtered_pass.vcf.idx
    |   |   | ...
    |   |   |─── cell_type2
    |   |   │      cell_1_filtered_pass.vcf
    |   |   |      cell_1_filtered_pass.vcf.idx
    |   |   │      cell_2_filtered_pass.vcf
    |   |   |      cell_2_filtered_pass.vcf.idx
    |   |   |  ...
    |   |─── snv_matrix
    |   |   │  celltype1_SNV_matrix.csv
    |   |   │  celltype2_SNV_matrix.csv
    |   |   |  ...
    |   |─── snv_meta
    |   |   |  celltype1_SNV_meta.csv
    |   |   |  ...
    |─── QTLs
