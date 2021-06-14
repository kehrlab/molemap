

def getbarcode(line):
    barcode=line.split('\t')[3]
    # print(barcode)
    return barcode

def get10xbarcode(line):
    # print("line: ",line,"\n")
    barcode=line.split(' ')[1][5:21]
    # print(barcode, " ")
    return barcode

def get10xID(line):
    id=line.split(' ')[0][1:]
    return id

def getbwaID(line):
    id=line.split('\t')[0]
    return id

def getmapping(line):
    mapping=line.split('\t')[0:3]
    # print(mapping)
    return mapping

def evaluate(bwa_line, mappings):
    bwamap=bwa_line.split('\t')
    give=20000
    for mapping in mappings:
        if mapping[0]==bwamap[2]:
            if (int(bwamap[3])>int(mapping[1])-give and int(bwamap[3])<int(mapping[2])+give and int(bwamap[7])>int(mapping[1])-give and int(bwamap[7])<int(mapping[2])+give):
                # print("mapping: ",mapping,"\n")
                return 1

    # print("bwamap: ",bwamap,"\n")
    return 0
# parameters:
tp_per=0.9 #fraction of reads that have to be bwa_mapped to an bcmap_identified position to count BC as TP

#files:
bcmap_res=open('resallsorted.bed','r')
bwa_res=open('resallbwa.sam','r')
readfile=open('./testdata/new_and_corrected.1.fastq','r')
bwa_line=bwa_res.readline()
while bwa_line[0]=='@':
    bwa_line=bwa_res.readline()
readfileline=readfile.readline()
tenXbc=get10xbarcode(readfileline)
tenXid=get10xID(readfileline)
readfile.readline()
readfile.readline()
readfile.readline()
# while getbwaID(bwa_line)==tenXid:
#     bwa_line=bwa_res.readline()

while tenXbc=='*':
    while getbwaID(bwa_line)==tenXid:
        bwa_line=bwa_res.readline()
    readfileline=readfile.readline()
    tenXbc=get10xbarcode(readfileline)
    tenXid=get10xID(readfileline)
    readfile.readline()
    readfile.readline()
    readfile.readline()

# print("tenXbc: ", tenXbc)
# print("readfile: ", readfileline)
# print("bwa_res: ", bwa_line)
# print(readfile.readline())

# readcount=0
barcodecount=0
tp=0
old_barcode=getbarcode(bcmap_res.readline())
mappings=[]
for line in bcmap_res:
    barcode=getbarcode(line)
    if barcode==old_barcode:
        #append mappings
        mappings+=[getmapping(line)]

    else:
        barcodecount+=1
        #evaluate
        correct=0
        reads=0
        # tenXbc=get10xbarcode(readfile.readline())
        # readfile.readline()
        # readfile.readline()
        # readfile.readline()
        print("tenXbc: ", tenXbc, " old_barcode: ", old_barcode, " comparison: " , tenXbc<old_barcode)
        while tenXbc<old_barcode:# or tenXbc=='*':      # skipping till barcodes match
            while getbwaID(bwa_line)==tenXid:
                bwa_line=bwa_res.readline()
            readfileline=readfile.readline()
            tenXbc=get10xbarcode(readfileline)
            tenXid=get10xID(readfileline)
            readfile.readline()
            readfile.readline()
            readfile.readline()
            # print("tenXbc: ", tenXbc, " old_barcode: ", old_barcode, " comparison: " , tenXbc<old_barcode, "\n")
            # bwa_res.readline()
            # bwa_res.readline()
        #evaluate bwa_line
        while tenXbc==old_barcode:
            # bwa_line=bwa_res.readline()
            while getbwaID(bwa_line)==tenXid:
                correct+=evaluate(bwa_line, mappings)
                reads+=1
                bwa_line=bwa_res.readline()
            readfileline=readfile.readline()
            if readfileline!="":
                tenXbc=get10xbarcode(readfileline)
                tenXid=get10xID(readfileline)
                readfile.readline()
                readfile.readline()
                readfile.readline()
            else:
                break;

        # bwa_line=bwa_res.readline()
        # reads+=1
        # correct+=evaluate(bwa_line, mappings)
        # bwa_line=bwa_res.readline()
        # reads+=1
        # correct+=evaluate(bwa_line, mappings)
        #create new mappings
        print("bc: ",old_barcode,"\n")
        print(mappings,"\n")
        print("reads:   ", reads, "\n")
        print("correct: ", correct, "\n\n")
        old_barcode=barcode
        mappings=[[]]
        mappings[0]=getmapping(line)

print("\n")
# print("readcount:     " , readcount , "\n")
print("barcodecount:  " , barcodecount , "\n")
print("true positives:" , tp, "\n")
