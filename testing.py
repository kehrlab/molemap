

def getbarcode(line):
    barcode=line.split('\t')[3]
    print(barcode)
    return barcode

def getmapping(line):
    mapping=line.split('\t')[0:2]
    print(mapping)
    return mapping

# def evaluate(bwa_line, mappings):
#     for mapping in mappings:
#         if #ref correct
#             if #pos correct
#                 return 1
# parameters:
tp_per=0.9 #fraction of reads that have to be bwa_mapped to an bcmap_identified position to count BC as TP

#files:
bcmap_res=open('res5bc.bed','r')
bwa_res=open('./bwa/resall32.sam','r')
bwa_line=bwa_res.readline()
while bwa_line[0]=='@':
    bwa_line=bwa_res.readline()
# print('bwa_line#1: ' , bwa_line, "\n")
readfile=open('./testdata/NA12878_WGS_v2_S1_L001_all_corrected.1.fastq','r')
readcount=0
barcodecount=0
tp=0
old_barcode=''
mappings=[[]]
for line in bcmap_res:
    readcount+=1
    barcode=getbarcode(line)
    if barcode==old_barcode:
        #append mappings
        mappings+=getmapping(line)

    else:
        barcodecount+=1
        # #evaluate
        # while getbarcode(readfile1.readline())<old_barcode:
        #     readfile1.readline()
        #     readfile1.readline()
        #     readfile1.readline()
        #     bwa_res.readline()
        # bwa_line=bwa_res.readline()
        # tp+=evaluate(bwa_line, mappings)
        # #evaluate bwa_line
        # while getbarcode(readfile1.readline())==old_barcode:
        #     bwa_line=bwa_res.readline()


        #create new mappings
        old_barcode=barcode
        mappings=[[]]
        mappings+=getmapping(line)

print("\n")
print("readcount:     " , readcount , "\n")
print("barcodecount:  " , barcodecount , "\n")
print("true positives:" , tp, "\n")
