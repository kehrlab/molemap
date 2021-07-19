import sys

def getbarcode(line):
    barcode=line.split('\t')[22][5:21]
    # print("barcode: ",barcode)
    return barcode

def get10xbarcode(line):
    # print("line: ",line,"\n")
    barcode=line.split(' ')
    if len(barcode)>1:
        barcode=barcode[1][5:21]
    else:
        barcode="BAD_BARCODE_____"
    # print(barcode, " ")
    return barcode

def get10xID(line):
    id=line.split(' ')[0][1:]
    # print("10xid: " , id)
    return id

def getbwaID(line):
    id=line.split('\t')[0]
    # print("id:    " , id)
    return id

def getmapping(line):
    mapping=line.split('\t')[0:3]
    # print(mapping)
    return mapping

def evaluate(lr_line, mappings):
    bwamap=lr_line.split('\t')
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
    ref="*"
    pos=0
    FN=0
    for map in unmapped:
        if ref==map[0]:
            if int(map[1])-pos<300000:
                if int(map[1])-pos>200:
                    cluster+=1
                    pos=int(map[1])
            else:
                if cluster>10:
                    # print("\n\n" , unmapped , "\n\n")
                    FN+=1

                pos=int(map[1])
                cluster=1

        else:
            if cluster>10:
                # print("\n\n" , unmapped , "\n\n")
                FN+=1
            ref=map[0]
            pos=int(map[1])
            cluster=1
    if cluster>10:
        # print("\n\n" , unmapped , "\n\n")
        FN+=1
    return FN;






# parameters:
tp_per=0.5 #fraction of reads that have to be bwa_mapped to an bcmap_identified position to count BC as TP
# rfile="resallq" + str(sys.argv[1]) + "sorted.bed"
rfile="/fast/users/luepkenr_c/work/P03/Barcode_assignment/reslongrangersorted.bed"
#files:
bcmap_res=open(rfile,'r')
# lr_res=open('/fast/users/luepkenr_c/scratch/BIH_TRASH/2021-06-18/resallbwa.sam','r')
lr_res=open('/fast/users/luepkenr_c/scratch/longranger/lariat_bc_sorted_valid.sam','r')
# readfile=open('/fast/users/luepkenr_c/scratch/longranger/longranger.1.fq','r')
lr_line=lr_res.readline()
lr_line=lr_res.readline()
print(lr_line)
print(getbarcode(lr_line))


while lr_line[0]=='@':
    lr_line=lr_res.readline()
readfileline=readfile.readline()
tenXbc=get10xbarcode(readfileline)
tenXid=get10xID(readfileline)
readfile.readline()
readfile.readline()
readfile.readline()
# while getbwaID(lr_line)==tenXid:
#     lr_line=lr_res.readline()

while tenXbc=='*':
    while getbwaID(lr_line)==tenXid:
        lr_line=lr_res.readline()
    readfileline=readfile.readline()
    tenXbc=get10xbarcode(readfileline)
    tenXid=get10xID(readfileline)
    readfile.readline()
    readfile.readline()
    readfile.readline()

# print("tenXbc: ", tenXbc)
# print("readfile: ", readfileline)
# print("lr_res: ", lr_line)
# print(readfile.readline())

# readcount=0
barcodecount=0
TP=1
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
            while getbwaID(lr_line)==tenXid:
                lr_line=lr_res.readline()
            readfileline=readfile.readline()
            tenXbc=get10xbarcode(readfileline)
            tenXid=get10xID(readfileline)
            readfile.readline()
            readfile.readline()
            readfile.readline()
            # print("tenXbc: ", tenXbc, " old_barcode: ", old_barcode, " comparison: " , tenXbc<old_barcode, "\n")
            # lr_res.readline()
            # lr_res.readline()
        #evaluate lr_line
        while tenXbc==old_barcode:
            # lr_line=lr_res.readline()
            while getbwaID(lr_line)==tenXid:
                res=evaluate(lr_line, mappings)
                if res!=-1:
                    mappinglist[res]+=1
                else:
                    unmapped+=[lr_line.split('\t')[2:5]];
                reads+=1
                lr_line=lr_res.readline()
            readfileline=readfile.readline()
            if readfileline!="":
                tenXbc=get10xbarcode(readfileline)
                tenXid=get10xID(readfileline)
                readfile.readline()
                readfile.readline()
                readfile.readline()
            else:
                break;

        # lr_line=lr_res.readline()
        # reads+=1
        # correct+=evaluate(lr_line, mappings)
        # lr_line=lr_res.readline()
        # reads+=1
        # correct+=evaluate(lr_line, mappings)
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
        # if sum(mappinglist)/reads<tp_per:
        #     FN+=1
        # print("barcodecount: ",barcodecount)
        # print("FN: " , FN ," FP: ", FP , " TP: ",TP)
        # print("Precision: " , round(TP/(TP+FP)*100,2), "Recall: " ,round(TP/(TP+FN)*100,2))

        # print("\n\n" , unmapped , "\n" , mappings, "\n\n")
        old_barcode=barcode
        mappings=[[]]
        mappings[0]=getmapping(line)
        unmapped=[]

print("\n\n")
# print("readcount:     " , readcount , "\n")
print("score_threshold: ", str(sys.argv)[1])
print("barcodecount:  " , barcodecount , "\n")
print("TP: " , TP)
print("FP: " , FP)
print("FN: " , FN , "\n\n")
print("Precision: " , round(TP/(TP+FP)*100,2), "Recall: " ,round(TP/(TP+FN)*100,2))