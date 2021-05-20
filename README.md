# Barcode_assignment
Maps barcodes to a reference genome and returns genomic windows from which the barcoded reads most likely originate. Each window is assessed with a quality score.

# Installation

## Prerequisites
- Seqan
- gcc version XXX
- Other libs
Needed directories etc.

# Usage 
1. Build an index of the reference genome using countK
2. Map barcodes to reference using bcmap

## Data requirements
- 10XGenomics paired-end Linked-reads
- Sorted by barcode (use [bcctools](https://github.com/kehrlab/bcctools))
- Barcodes are stored in RX:Z: flag of read Ids
- Barcode whitelist required (can be infered using [bcctools](https://github.com/kehrlab/bcctools))

## countK
- ./countK reference.fa IndexName

## bcmap
- ./bcmap readfile1.fastq readfile2.fastq IndexName BarcodeIndexName
 
