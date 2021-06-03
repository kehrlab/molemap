

def getbarcode(line):
    barcode=line.split('\t')[3]
    print(barcode)
    return barcode

def getmapping(line):
    mapping=line.split('\t')[0:2]
    pring(mapping)
    return mapping


bcmap_res=open('resall.bed','r')
# bwa_res=open('','r')
# readfile=open('159916111600 25. Feb 14:36 NA12878_WGS_v2_S1_L001_all_corrected.1.fastq','r')

old_barcode=''
mappings=[[]]
for line in bcmap_res:
    barcode=getbarcode(line)
    if barcode==old_barcode:
        mappings+=getmapping(line)
        #append mappings

    else:
        mappings=[[]]
        #evaluate and create new mappings
