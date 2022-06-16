# A workflow to extract reads from a region of interest
This snakemake workflow chains bcmap commands and unix commands to extract all reads from user defined regions of interest. 

## Prerequisites
bcmap installed

Snakemake ([available via conda](https://anaconda.org/bioconda/snakemake))

## Data requirements
- Paired-end Linked-reads
- Barcodes are stored in BX:Z: flag of read Ids
- Sorted by barcode (use i.e. [bcctools](https://github.com/kehrlab/bcctools) or [samtools](https://github.com/samtools/samtools))

- A reference genome in fasta or fastq format

## Usage
Update set_up.yml to define your input files and name your output files.

Run the workflow specifying the available number of cpus:

    snakemake --cores {#cpus}
