

def getbarcode(line):
    barcode=line.split('\t')[3]
    print(barcode)
    return barcode

bcmap_res=open('resall.bed','r')
# bwa_res=open('','r')
# readfile=open('159916111600 25. Feb 14:36 NA12878_WGS_v2_S1_L001_all_corrected.1.fastq','r')

old_barcode=''

for line in bcmap_res:
    barcode=getbarcode(line)
    # if barcode==old_barcode
        #append mappings

    # else
        #evaluate and create new mappings
