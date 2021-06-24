

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
    res=-1
    for mapping in mappings:
        res+=1
        if mapping[0]==bwamap[2]:
            if (int(bwamap[3])>int(mapping[1])-give and int(bwamap[3])<int(mapping[2])+give):
                # print("mapping: ",mapping,"\n")
                return res

    # print("bwamap: ",bwamap,"\n")
    return -1


def cluster(unmapped):
    unmapped=[map for map in unmapped if int(map[2])>0 and map[0]!='*']
    unmapped.sort()
    # print("\n\n" , unmapped , "\n\n")
    cluster=0
    ref=0
    pos=0
    FN=0
    for map in unmapped:
        if ref==map[0]:
            if map[1]-pos<300000:
                cluster+=1
            else:
                if cluster>5:
                    print("\n\n" , unmapped , "\n\n")
                    FN+=1

                pos=map[1]:
                cluster=1

        else:
            if cluster>5:
                print("\n\n" , unmapped , "\n\n")
                FN+=1
            ref=map[0]
            pos=map[1]
            cluster=1
    if cluster>5:
        print("\n\n" , unmapped , "\n\n")
        FN+=1
    return FN;






# parameters:
tp_per=0.5 #fraction of reads that have to be bwa_mapped to an bcmap_identified position to count BC as TP

#files:
bcmap_res=open('resallsorted.bed','r')
bwa_res=open('/fast/users/luepkenr_c/scratch/BIH_TRASH/2021-06-18/resallbwa.sam','r')
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
TP=0
FN=0
FP=0
old_barcode=getbarcode(bcmap_res.readline())
mappings=[]
unmapped=[]
for line in bcmap_res:
    barcode=getbarcode(line)
    if barcode==old_barcode:
        #append mappings
        mappings+=[getmapping(line)]

    else:
        mappinglist=[0]*len(mappings)
        barcodecount+=1
        #evaluate
        correct=0
        reads=0
        # tenXbc=get10xbarcode(readfile.readline())
        # readfile.readline()
        # readfile.readline()
        # readfile.readline()
        # print("tenXbc: ", tenXbc, " old_barcode: ", old_barcode, " comparison: " , tenXbc<old_barcode)
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
                res=evaluate(bwa_line, mappings)
                if res!=-1:
                    mappinglist[res]+=1
                else:
                    unmapped+=[bwa_line.split('\t')[2:5]];
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
        # print("bc: ",old_barcode,"\n")
        # print(mappings,"\n")
        # print("reads:   ", reads, "\n")
        # print("correct: ", correct, "\n\n")
        for mapping in mappinglist:
            if mapping==0:
                FP+=1
            else:
                TP+=1

        FN+=cluster(unmapped)
        if sum(mappinglist)/reads<tp_per:
            FN+=1
        print("FN: " , round(FN/barcodecount*100,2),"% ","FP: ", round(FP/(FP+TP),2)*100,"% TP: ", round(TP/(FP+TP),2)*100,"%")

        # print("\n\n" , unmapped , "\n" , mappings, "\n\n")
        old_barcode=barcode
        mappings=[[]]
        mappings[0]=getmapping(line)
        unmapped=[]

print("\n")
# print("readcount:     " , readcount , "\n")
print("barcodecount:  " , barcodecount , "\n")
