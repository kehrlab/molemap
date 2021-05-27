# Barcode_assignment
Maps barcodes to a reference genome and returns genomic windows from which the barcoded reads most likely originate. Each window is assessed with a quality score.

# Installation
- Install [SeqAn2](https://seqan.readthedocs.io/en/master/Infrastructure/Use/Install.html#infra-use-install)
- Update line XX in the makefile to link and include the SeqAn directory
- run make

## Prerequisites
- SeqAn2
- gcc version XXX

# Usage 
1. Build an index of the reference genome using countK
2. Map barcodes to reference using bcmap

- arguments and options are explained using --help for every command

## Data requirements
- 10XGenomics paired-end Linked-reads
- Sorted by barcode (use [bcctools](https://github.com/kehrlab/bcctools))
- Barcodes are stored in BX:Z: flag of read Ids

## countK
- ./countK reference.fa IndexName

## bcmap
- ./bcmap readfile1.fastq readfile2.fastq IndexName BarcodeIndexName

## getreads
- ./getreads readfile1.fastq readfile2.fastq BarcodeIndexName Barcodes
 
