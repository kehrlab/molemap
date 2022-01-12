# Local assembly of a 766 bp structural variant

## Linked-read extraction

Task: Construct kmer index of GRCh38 \
Input: [GRCh38.fq] human reference genome \
Output: [Index_GRCh38] folder containing the bcmap kmer index

    ./bcmap index GRCh38.fq -o Index_GRCh38

Task: Map barcodes of NA12878 \
Input: [NA12878_linked_reads_1.fq, NA12878_linked_reads_2.fq], [Index_GRCh38] \
Output: [NA12878_mapped.bed] barcode index produced by bcmap \
        [NA12878_read_index] read index for NA12878 read files
  
    ./bcmap map NA12878_linked_reads_1.fq NA12878_linked_reads_2.fq -i Index_GRCh38 -o NA12878_mapped.bed -b NA12878_read_index

Task: Extract Barcodes of interest from barcode index \
Input: [NA12878_mapped.bed] \
Output: [Barcodes_NA12878_50K_region.txt] barcodes mapping to region of interest

    awk '{if($1=="chr17" && (int($2)<=17831079+25000 && int($3)>=17831079-25000)) print($0)}' NA12878_mapped.bed | awk '{if(int($5)>8) print($4)}' > Barcodes_NA12878_50K_region.txt

Task: Extract reads from read index \
Input:  [NA12878_linked_reads_1.fq, NA12878_linked_reads_2.fq], [Barcodes_NA12878_50K_region.txt], [NA12878_read_index] \
Output: [NA12878_50K_region.fq]

    ./bcmap get NA12878_linked_reads_1.fq NA12878_linked_reads_2.fq Barcodes_NA12878_50K_region.txt -b NA12878_read_index -o NA12878_50K_region.fq 

## Linked-read assembly

Task: Assembly of PE linked-reads.
Input: [NA12878_50K_region.fq] An input file of interleaved pairs of linked-reads.
Output: [assembly_k121.unitigs.fa] A set of unitig from the assemblers final iteration.

  gatb --12 NA12878_50K_region.fa --no-scaffolding --nb-cores 8 > logs/assembly.log 2>&1


Task: Simplifying the set of unitigs in their de Bruijn Graph representation, i.e. removing tips and singletons.
Input: [assembly_k121.unitigs.fa] The set of unitigs from the previous step.
Output: [assembly_k121.unitigs.fa.gfa] A simplified de Bruijn Graph in GFA format.

  Bifrost build -r assembly_k121.unitigs.fa -t 8 -k 121 -m 81 --clip-tips --del-isolated -o assembly_k121.unitigs.fa

Task: Extract unitigs from simplified de Bruijn Graph.
Input: [assembly_k121.unitigs.fa.gfa] The simplified de Bruijn Graph from the previous step.
Output: [assembly_k121.unitigs.bifrost.fa] The unitigs of the simplified de Bruijn Graph in FASTA format.

  awk '$0 ~ /^S/ {print ">"$2"\n"$3}' assembly_k121.unitigs.fa.gfa > assembly_k121.unitigs.bifrost.fa

## Assembly validation.

![plot](./766bp-NRS.png)
Alignment of the unitigs from the local assembly to the HG38 reference using the UCSC web application of BLAT \cite{kent_blatblast-like_2002

## Sequence validation of the non-reference sequence variant.
Local sequence alignment of the continuous subsequence of unitig 264 that does not align to the reference genome and the insertion sequence at breakpoint chr17:17,831,079 reported in the GIAB callset.
The two sequences, both of 766 bp length, align with 100\% identity and similarity. 
The alignment was computed using the EMBOSS Water web application \thomas{https://academic.oup.com/nar/article/47/W1/W636/5446251} of the Smith-Waterman algorithm for local sequence alignment.

    ########################################
    # Program: water
    # Rundate: Thu  6 Jan 2022 16:01:26
    # Commandline: water
    #    -auto
    #    -stdout
    #    -asequence emboss_water-I20220106-160436-0340-65113364-p2m.aupfile
    #    -bsequence emboss_water-I20220106-160436-0340-65113364-p2m.bupfile
    #    -datafile EDNAFULL
    #    -gapopen 10.0
    #    -gapextend 0.5
    #    -aformat3 pair
    #    -snucleotide1
    #    -snucleotide2
    # Align_format: pair
    # Report_file: stdout
    ########################################

    #=======================================
    #
    # Aligned_sequences: 2
    # 1: chr17_17831079_pbsv.INS.40992
    # 2: 264
    # Matrix: EDNAFULL
    # Gap_penalty: 10.0
    # Extend_penalty: 0.5
    #
    # Length: 766
    # Identity:     766/766 (100.0%)
    # Similarity:   766/766 (100.0%)
    # Gaps:           0/766 ( 0.0%)
    # Score: 3830.0
    # 
    #
    #=======================================

## References
[Kent, 2002]Kent,W.J.(2002)BLAT—The BLAST-Like Alignment Tool.Genome Research,12(4), 656–664. Company: Cold Spring Harbor Laboratory Press Distributor: Cold SpringHarbor Laboratory Press Institution: Cold Spring Harbor Laboratory Press Label: ColdSpring Harbor Laboratory Press Publisher: Cold Spring Harbor Lab.
