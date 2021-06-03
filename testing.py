

def getbarcode(line):
    barcode=line.split('\t')[3]
    # print(barcode)
    return barcode

def get10xbarcode(line):
    print("line: ",line,"\n")
    barcode=line.split(' ')[1][5:]
    print(barcode, " ")
    return barcode

def getmapping(line):
    mapping=line.split('\t')[0:3]
    # print(mapping)
    return mapping

def evaluate(bwa_line, mappings):
    # bwamap=bwa_line.split('\t')
    # print("bwamap: ",bwamap,"\n")
    # for mapping in mappings:
    #     if mapping[0]==bwamap[2]:
    #         if (int(bwamap[3])>int(mapping[1]) and int(bwamap[3])<int(mapping[2]) and int(bwamap[7])>int(mapping[1]) and int(bwamap[7])<int(mapping[2])):
    #             print("mapping: ",mapping,"\n")
    #             return 1
    return 0
# parameters:
tp_per=0.9 #fraction of reads that have to be bwa_mapped to an bcmap_identified position to count BC as TP

#files:
bcmap_res=open('res5bc.bed','r')
bwa_res=open('res5bc.sam','r')
bwa_line=bwa_res.readline()
while bwa_line[0]=='@':
    bwa_line=bwa_res.readline()
# print('bwa_line#1: ' , bwa_line, "\n")
readfile=open('./testdata/NA12878_WGS_v2_S1_L001_all_corrected.1.fastq','r')
# readcount=0
barcodecount=0
tp=0
old_barcode=''
mappings=[[]]
for line in bcmap_res:
    barcode=getbarcode(line)
    if barcode==old_barcode:
        #append mappings
        mappings+=[getmapping(line)]

    else:
        barcodecount+=1
        #evaluate
        correct=0
        tenXbc=get10xbarcode(readfile.readline())
        while tenXbc<old_barcode and tenXbc!='*':
            print(readfile1.readline(),"\n")
            print(readfile1.readline(),"\n")
            print(readfile1.readline(),"\n")
            tenXbc=get10xbarcode(readfile.readline())
            bwa_res.readline()
            bwa_res.readline()
        bwa_line=bwa_res.readline()
        correct+=evaluate(bwa_line, mappings)
        bwa_line=bwa_res.readline()
        correct+=evaluate(bwa_line, mappings)
        #evaluate bwa_line
        while get10xbarcode(readfile.readline())==old_barcode:
            readfile1.readline()
            readfile1.readline()
            readfile1.readline()
            bwa_line=bwa_res.readline()
            correct+=evaluate(bwa_line, mappings)
            bwa_line=bwa_res.readline()
            correct+=evaluate(bwa_line, mappings)
        #create new mappings
        print("bc: ",old_barcode,"\n")
        print(mappings,"\n\n")
        old_barcode=barcode
        mappings=[[]]
        mappings[0]=getmapping(line)

print("\n")
# print("readcount:     " , readcount , "\n")
print("barcodecount:  " , barcodecount , "\n")
print("true positives:" , tp, "\n")
