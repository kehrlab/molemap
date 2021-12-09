# bcmap
Maps barcodes to a reference genome and returns genomic windows from which the barcoded reads most likely originate. Each window is assessed with a quality score representing the trustworthines of the mapping. A barcode index is constructed alongside the mapping and can be used to quickly retrieve all reads belonging to a barcode.

## Prerequisites
- gcc version 7.2.0

# Installation
    git clone https://github.com/kehrlab/bcmap.git
    cd bcmap
    make

# Data requirements
- 10XGenomics paired-end Linked-reads
- Sorted by barcode (use i.e. [bcctools](https://github.com/kehrlab/bcctools) or [samtools](https://github.com/samtools/samtools))

    ./script/run_bcctools -f fastq first.fq.gz second.fq.gz

- Barcodes are stored in BX:Z: flag of read Ids

# Commands
For detailed information on Arguments and Options:

    ./bcmap [command] --help

## index
Builds an minimized open addressing k-mer index of the reference genome. The index is required to run "map".

    ./bcmap index reference.fa [options]

## map
Maps the barcodes of the provided readfiles to the reference and creates a barcode index of the readfiles to quickly retrieve all reads of a given barcode.

    ./bcmap map readfile1.fastq readfile2.fastq [options]

Content of output bed-file:
* *chromosome  startposition  endposition  barcode  mapping_score*

![plot](./example/Example_score_hist.png)

Bcmap returns a _output.hist_ file that can be ploted using _plot_score_histogram.py_ resulting in a plot like the one above. To create a set of mappings with very high precision (at the cost of some recall), the local minimum inbetween the two peaks should be set as the score threshold. A lower theshold yields better recall at the cost of precision, a higher threshold is not recomended.

## get
Returns all reads of the given barcodes. Barcodes can be provided directly as argument or in a file.

    ./bcmap get readfile1.fastq readfile2.fastq Barcodes [options]
 
# Example 
This small example demonstrates how to use bcmap and allows you to check if it is properly installed. Navigate to the bcmap folder and run the commands listed below.  

    # building the index for chr21.fa
    ./bcmap index example/chr21.fa -o example/Index
    
    # mapping the reads of readfile 1 and 2 to chromosome 21
    ./bcmap map example/readfile.1.fq example/readfile.2.fq -i example/Index -b example/BarcodeIndex -o example/results.bed
    
    # extracting the first barcode from the results
    awk '{if(NR==1) print($4)}' example/results.bed > example/FirstBarcode.txt
    
    # extracting all reads belonging to the first barcode
    ./bcmap get example/readfile.1.fq example/readfile.2.fq example/FirstBarcode.txt -b example/BarcodeIndex -o example/readsOfFirstBarcode.fq
    
    # extracting reads of barcode AACATCGCAAACAGTA
    ./bcmap get example/readfile.1.fq example/readfile.2.fq AACATCGCAAACAGTA -b example/BarcodeIndex -o example/readsOfAACATCGCAAACAGTA.fq

