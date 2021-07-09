# bcmap
Maps barcodes to a reference genome and returns genomic windows from which the barcoded reads most likely originate. Each window is assessed with a quality score.

## Prerequisites
- SeqAn2
- gcc version XXX

# Installation
- Install [SeqAn2](https://seqan.readthedocs.io/en/master/Infrastructure/Use/Install.html#infra-use-install)
- Update line 14 and 15 in the makefile to point to the SeqAn /lib and /include directory
- run make

# Data requirements
- 10XGenomics paired-end Linked-reads
- Sorted by barcode (use i.e. [bcctools](https://github.com/kehrlab/bcctools) or [samtools](https://github.com/samtools/samtools))
- Barcodes are stored in BX:Z: flag of read Ids

# Commands
For detailed information on Arguments and Options:

    ./bcmap [command] --help

## index
Builds an open addressing k-mer index of the reference genome. The index is required to run "map". The bucket count should be set to a value close to the size of the reference using the -b option. The default bucket number is optimized for hg38 and k=31.

    ./bcmap index reference.fa RefIndexName

## map
Maps the barcodes of the provided readfiles to the reference and creates a barcode index of the readfiles to quickly retrieve all reads of a given barcode.

    ./bcmap map readfile1.fastq readfile2.fastq RefIndexName BarcodeIndexName

## get
Returns all reads of the given barcodes. Barcodes can be provided directly as argument or in a file.

    ./bcmap get readfile1.fastq readfile2.fastq BarcodeIndexName Barcodes
 
# Example 
This small example demonstrates how to use bcmap and allows you to check if it is properly installed. Navigate do the bcmap folder and run the commands listed below. 

    # building the index for chr21.fa using 45000000 buckets (based on the size of chromosome 21)
    ./bcmap index example/chr21.fa example/Index -b 45000000
    
    # mapping the reads of readfile 1 and 2 to chromosome 21
    ./bcmap map example/readfile.1.fq example/readfile.2.fq example/Index example/BarcodeIndex -o example/results.bed
    
    # extracting the first barcode from the results
    awk '{if(NR==1) print($4)}' example/results.bed > example/FirstBarcode.txt
    
    # extracting all reads belonging to the first barcode
    ./bcmap get example/readfile.1.fq example/readfile.2.fq example/BarcodeIndex example/FirstBarcode.txt -o example/readsOfFirstBarcode.fq
    
    # extracting reads of barcode AACATCGCAAACAGTA
    ./bcmap get example/readfile.1.fq example/readfile.2.fq example/BarcodeIndex AACATCGCAAACAGTA -o example/readsOfAACATCGCAAACAGTA.fq

