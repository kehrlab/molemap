# bcmap
Maps barcodes to a reference genome and returns genomic windows from which the barcoded reads most likely originate. Each window is assessed with a quality score.

# Installation
- Install [SeqAn2](https://seqan.readthedocs.io/en/master/Infrastructure/Use/Install.html#infra-use-install)
- Update line XX in the makefile to link and include the SeqAn directory
- run make

## Prerequisites
- SeqAn2
- gcc version XXX

## Data requirements
- 10XGenomics paired-end Linked-reads
- Sorted by barcode (use [bcctools](https://github.com/kehrlab/bcctools))
- Barcodes are stored in BX:Z: flag of read Ids

# Example 

    # building the index for chr21.fa using 45000000 buckets (based on the size of chromosome 21)
    ./bcmap index example/chr21.fa example/Index -b 45000000
    
    # mapping the reads of readfile 1 and 2 to chromosome 21
    ./bcmap map example/readfile.1.fq example/readfile.2.fq example/Index example/BarcodeIndex -o example/results.bed
    
    # extracting the first barcode from the results
    awk 'BEGIN {print($4)}' results.bed > FirstBarcode.txt
    
    # extracting all reads belonging to this barcode
    ./bcmap get example/readfile.1.fq examle/readfile.2.fq example/BarcodeIndex barcodes.txt -o readsOfFirstBarcode.fq


# Commands
For detailed information on Arguments and parameters use:

    ./bcmap [command] --help

## index
Builds an open addressing k-mer index of the reference genome. The index is required to run "map".

    ./bcmap index reference.fa RefIndexName

## map
Maps the barcodes of the provided readfiles to the reference and creates a barcode index of the readfiles to quickly retrieve all reads of a given barcode.

    ./bcmap map readfile1.fastq readfile2.fastq RefIndexName BarcodeIndexName

## get
Returns all reads of the given barcodes. Barcodes can be provided directly as argument or in a file.

    ./bcmap get readfile1.fastq readfile2.fastq BarcodeIndexName Barcodes
 
