# Barcode_assignment

# Installation

## Prerequisites
- Seqan
- gcc version XXX
- Other libs
Needed directories etc.


# Usage 
1. build an index of the reference genome using countK
2. Map barcodes to reference using bcmap

## Data requirements
- 10XGenomics paired end Linked reads
- sorted by barcode (use [bcctools](https://github.com/kehrlab/bcctools))
- barcodes are stored in -RX: flag of read Ids
- Barcode whitelist required (can be infered using [bcctools](https://github.com/kehrlab/bcctools))

## countK
- ./countK reference.fa IndexName

## bcmap
- ./bcmap readfile1.fastq readfile2.fastq IndexName BCI(holder) -o outputFile
 
