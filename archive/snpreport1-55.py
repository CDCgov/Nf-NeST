import argparse
import pandas as pd
import subprocess
import csv
import itertools
def install(name):
    subprocess.call(['pip', 'install', name])

install('xlrd==1.2.0')
import xlrd

parser = argparse.ArgumentParser(description='name')
parser.add_argument('-v1', dest='unfiltered', type=str, help="name of unfilterd merged vcf file")
parser.add_argument('-v2', dest='filtered', type=str, help="name of filtered merged vcf file")
parser.add_argument('-b1', dest='bam', type=str, help="name of bam file")
parser.add_argument('-b2', dest='bed', type=str, help="name of bed file")
parser.add_argument('-o1', dest='name', type=str, help="name of output")
parser.add_argument('-e1', dest='candidates', type=str, help="name of candidates")
parser.add_argument('-e2', dest='voi', type=str, help="name of variants of interest")
parser.add_argument('-f1', dest='fasta', type=str, help="name of fasta file")
parser.add_argument('-c', dest='coveragefile',type=str, help="name of coverage files")

args = parser.parse_args()

AAdic= {"Ala":"A", "Arg":"R", "Asn":"N", "Asp":"D", "Asx":"B", "Cys":"C", "Glu":"E", "Gln":"Q", "Glx":"Z", "Gly":"G", "His":"H", "Ile":"I", "Leu":"L", "Lys":"K", "Met":"M", "Phe":"F", "Pro":"P", "Ser":"S", "Thr":"T", "Trp":"W", "Tyr":"Y", "Val":"V", "STOP":"STOP"}

currentcodon=""
keywords = [''.join(i) for i in itertools.product(["T","G","A","C"], repeat = 3)]
codondic={}
for x in keywords:
    currentcodon=""
    for y in x:
        currentcodon+=y
    #print(currentcodon)
    if currentcodon == "TTT" or currentcodon == "TTC":
        #print("True")
        codondic[currentcodon]="Phe"
    if currentcodon == "TTA" or currentcodon == "TTG": 
        codondic[currentcodon]="Leu"
    if currentcodon == "CTT" or currentcodon =="CTC" or currentcodon =="CTA" or currentcodon =="CTG":
        codondic[currentcodon]="Leu"
    if currentcodon == "ATT" or currentcodon =="ATC" or currentcodon =="ATA":
        codondic[currentcodon]="Ile"
    if currentcodon == "ATG":
        codondic[currentcodon]="Met"
    if currentcodon == "GTT" or currentcodon =="GTC" or currentcodon =="GTA" or currentcodon =="GTG":
        codondic[currentcodon]="Val"
    if currentcodon == "TCT" or currentcodon =="TCC" or currentcodon =="TCA" or currentcodon =="TCG":
        codondic[currentcodon]="Ser"
    if currentcodon == "CCT" or currentcodon =="CCC" or currentcodon =="CCA" or currentcodon =="CCG": 
        codondic[currentcodon]="Pro"
    if currentcodon == "ACT" or currentcodon =="ACC" or currentcodon =="ACA" or currentcodon =="ACG":
        codondic[currentcodon]="Thr"
    if currentcodon == "GCT" or currentcodon =="GCC" or currentcodon =="GCA" or currentcodon =="GCG":
        codondic[currentcodon]="Ala"
    if currentcodon == "TAT" or currentcodon =="TAC":
        codondic[currentcodon]="Tyr"
    if currentcodon == "CAT" or currentcodon =="CAC":
        codondic[currentcodon]="His"
    if currentcodon == "CAA" or currentcodon =="CAG":
        codondic[currentcodon]="Gln"
    if currentcodon == "AAT" or currentcodon =="AAC":
        codondic[currentcodon]="Asn"
    if currentcodon == "AAA" or currentcodon =="AAG":
        codondic[currentcodon]="Lys"
    if currentcodon == "GAT" or currentcodon =="GAC":
        codondic[currentcodon]="Asp"
    if currentcodon == "GAA" or currentcodon =="GAG":
        codondic[currentcodon]="Glu"
    if currentcodon == "TGT" or currentcodon== "TGC":
        #print("I am correct")
        codondic[currentcodon]="Cys"
    if currentcodon == "TGG":
        codondic[currentcodon]="Trp"
    if currentcodon == "AGA" or currentcodon =="AGG" or currentcodon =="CGT" or currentcodon =="CGC" or currentcodon =="CGA" or currentcodon =="CGG":
        codondic[currentcodon]="Arg"
    if currentcodon == "AGT" or currentcodon == "AGC":
        codondic[currentcodon]="Ser"
    if currentcodon == "GGT" or currentcodon =="GGC" or currentcodon =="GGA" or currentcodon =="GGG":
        codondic[currentcodon]="Gly"
    if currentcodon == "TGA" or currentcodon == "TAA" or currentcodon == "TAG":
        codondic[currentcodon]="STOP"
    #print(currentcodon)
    #print(currentcodon)
    #if currentcodon=="TGT" or currentcodon=="TGC":
    #    print("I am correct")
    #print(codondic[currentcodon])

#print(codondic)


#PfCRTStartlist1=[96,364,812,1157,1443,1638,1810,2020,2208,2413,2699,2891,3115]
#PfCRTEndlist1=[x+y for x,y in zip([96,364,812,1157,1443,1638,1810,2020,2208,2413,2699,2891,3115],[90,268,172,132,71,75,82,50,56,92,44,54,76])]
#DHFRStartlist1=[1]
#DHFREndlist1=[x+y for x,y in zip([1],[1827])]
#K13Startlist1=[1]
#K13Endlist1=[x+y for x,y in zip([1],[2181])]
#PfMDR1Startlist1=[1]
#PfMDR1Endlist1=[x+y for x,y in zip([1],[4260])]
#PfDHPSStartlist1=[1,313,2302]
#PfDHPSEndlist1=[x+y for x,y in zip([1,313,2302],[135,1868,115])]
#mitochondrialgenomeStartlist1=[734,1933,3492]
#mitochondrialgenomeEndlist1=[x+y for x,y in zip([734,1933,3492],[839,1538,1130])]

with open(args.filtered, "r") as f1:
    #with open(args.fasta, "r") as f2:
    Genelist1=[]
    POSlist1=[]
    Reflist1=[]
    Altlist1=[]
    #AAlist1=[]
    #AAlist2=[]
    #AAPoslist1=[]
    currentgene1=""
    currentpos1=""
    currentAminoacid=""
    currentposlen1=0
    for line in f1:
        if "p." in line:
            count=0
            for word in line.split():
                #print(line)
                if count==0:
                    Genelist1+=[word]
                    currentgene1=word
                    #print(word)
                if count==1:
                    POSlist1+=[word]
                    #currentposlen1=len(POSlist1)
                    #if currentgene1 == "PfCRT": 
                    #    for z in range(len(PfCRTStartlist1)):
                    #        if PfCRTStartlist1[z]<int(currentpos1)<PfCRTEndlist1[z]:
                    #            POSlist1+=[currentpos1]
                    #if currentgene1 == "DHFR":
                    #    for z in range(len(DHFRStartlist1)):
                    #        #print(type(DHFRStartlist1[z]))
                            #print((depthlist1[y-1][2]))
                    #        if DHFRStartlist1[z]<int(currentpos1)<DHFREndlist1[z]:
                    #            POSlist1+=[currentpos1]
                    #if currentgene1 == "K13": 
                    #    for z in range(len(K13Startlist1)):
                    #        if K13Startlist1[z]<int(currentpos1)<K13Endlist1[z]:
                    #            POSlist1+=[currentpos1]
                    #if currentgene1 == "PfMDR1":
                    #    for z in range(len(PfMDR1Startlist1)):
                    #        if PfMDR1Startlist1[z]<int(currentpos1)<PfMDR1Endlist1[z]:
                    #            POSlist1+=[currentpos1]
                    #if currentgene1 == "PfDHPS":
                    #    for z in range(len(PfDHPSStartlist1)):
                    #        if PfDHPSStartlist1[z]<int(currentpos1)<PfDHPSEndlist1[z]:
                    #            POSlist1+=[currentpos1]
                    #if currentgene1 == "mitochondrial_genome":
                    #    for z in range(len(mitochondrialgenomeStartlist1)):
                    #        if mitochondrialgenomeStartlist1[z]<int(currentpos1)<mitochondrialgenomeEndlist1[z]:
                    #            POSlist1+=[currentpos1]
                    #if len(POSlist1)==currentposlen1:
                        #print("test")
                    #   Genelist1=Genelist1[:-1]
                    #    break
                    #print(word)
                if count==2:
                    Reflist1+=[word]
                if word.startswith("c.") or word.startswith("n."):
                    Altlist1+=[word[-1]]
                if word.startswith("p."):
                    temppos=""
                    for x in word:
                        #print(word)
                        if x.isdigit():
                            #print(x)
                            temppos+=x
                    if word[2:5] not in AAdic or word[5+len(temppos):5+len(temppos)+3] not in AAdic:
                        Genelist1=Genelist1[:-1]
                        POSlist1=POSlist1[:-1]
                        Reflist1=Reflist1[:-1]
                        Altlist1=Altlist1[:-1]
                    #print(word[5+len(temppos):5+len(temppos)+3])
                count+=1
        elif ("p." not in line and "c." in line) or ("p." not in line and "n." in line):
            count=0
            currentAminoacid=""
            for word in line.split():
                #print(line)
                if count==0:
                    Genelist1+=[word]
                    currentgene1=word
                    #print(word)
                if count==1:
                    POSlist1+=[word]
                    currentpos1=word
                    #print(currentpos1)
                if count==2:
                    Reflist1+=[word]
                if word.startswith("c.") or word.startswith("n."):
                    Altlist1+=[word[-1]]
                count+=1
            #linecount=0
            #print(currentgene1)
            #print(currentpos1)
            #print(currentAminoacid)
            #print(currentAminoacid)
            #if currentAminoacid in codondic:
                #print("TRUE")
                #if codondic[currentAminoacid] in AAdic:
                    #print("TRUE")
                    #print(AAdic[codondic[line[int(currentpos1)-2:int(currentpos1)+1]]])
                    #AAlist1+=[AAdic[codondic[currentAminoacid]]]
                    #AAlist2+=[AAdic[codondic[currentAminoacid]]]
                    #AAPoslist1+=[int(currentpos1)/3]
                    #line[int(currentpos1)-2:int(currentpos1)+1]
            #else:
            #    print(currentAminoacid)
            #    print("false")
            #print(currentAminoacid)
    #print(len(Genelist1))
    #print(len(POSlist1))
    #print(len(Reflist1))
    #print(len(Altlist1))
    #print(len(AAlist1))
    #print(len(AAPoslist1))
    #print(len(AAlist2))
    #print(len(Genelist1))
    #print(len(AAlist2))
    #print(len(Genelist1))
    #print(Genelist1)
    #print(AAlist2)
    #print(POSlist1)
    process=subprocess.Popen(['samtools', 'depth', '-a', args.name+ "_SR.bam"], stdout=subprocess.PIPE)
    stdout, stderr = process.communicate()
    reader = csv.reader(stdout.decode('ascii').splitlines(),
                            delimiter='\t', skipinitialspace=True)
    depthlist1=[]
    for row in reader:
        depthlist1+=[row]
    depthlist1up=[]
    #for x in range(len(Genelist1)):
    #    if Genelist1[x] in depthlist1:
    #        print(Genelist1[x])
    #    if POSlist1[x] in depthlist1:
    #        print(POSlist1[x])
    #print(depthlist1)
    #print(depthlist1)
    #print(Genelist1)
    #print(POSlist1)
    #print(Genelist1)
    depthpair1=[]
    for x in range(len(depthlist1)):
        depthpair1+=[[depthlist1[x][0],depthlist1[x][1]]]

    for x in range(len(Genelist1)):
        if [Genelist1[x], POSlist1[x]] not in depthpair1:
            depthlist1up+=[None]

        else:
            for y in range(len(depthlist1)):
                if Genelist1[x] in depthlist1[y] and POSlist1[x] in depthlist1[y]:
                    #print(depthlist1[y])
                    #print(depthlist1[y][2])
                    depthlist1up+=[depthlist1[y][2]]
                    break
    #print(len(depthlist1up))
    #print(depthpair1)
    #print(depthlist1)
    #print(len(depthlist1up))
    #print(depthlist1up)
    #print(Genelist1)
    #print(POSlist1)
    codondepthlist1=[]
    #sadcount=0
    for x in range(len(Genelist1)):
        if [Genelist1[x],POSlist1[x]] not in depthpair1:
                #print("True")
                codondepthlist1+=[None]
                #break
        else:
            for y in range(len(depthlist1)):
                if Genelist1[x] in depthlist1[y] and POSlist1[x] in depthlist1[y]:
                    #print(Genelist1[x],POSlist1[x])
                    #print(depthlist1[y])
                    #sadcount+=1
                    tempcodondepthlist1=0
                    if Genelist1[x] == "PfCRT": 
                        #for z in range(len(PfCRTStartlist1)):
                            #if PfCRTStartlist1[z]<int(depthlist1[y-1][2])<PfCRTEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y-1][2])
                            #if PfCRTStartlist1[z]<int(depthlist1[y][2])<PfCRTEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y][2])
                            #if PfCRTStartlist1[z]<int(depthlist1[y+1][2])<PfCRTEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y+1][2])
                            #print(tempcodondepthlist1)
                        codondepthlist1+=[int(tempcodondepthlist1/3)]
                        break
                    if Genelist1[x] == "DHFR":
                        #for z in range(len(DHFRStartlist1)):
                            #print(type(DHFRStartlist1[z]))
                            #print((depthlist1[y-1][2]))
                            #if DHFRStartlist1[z]<int(depthlist1[y-1][2])<DHFREndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y-1][2])
                            #if DHFRStartlist1[z]<int(depthlist1[y][2])<DHFREndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y][2])
                            #if DHFRStartlist1[z]<int(depthlist1[y+1][2])<DHFREndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y+1][2])
                        codondepthlist1+=[int(tempcodondepthlist1/3)]
                        break
                    if Genelist1[x] == "K13": 
                        #for z in range(len(K13Startlist1)):
                            #if K13Startlist1[z]<int(depthlist1[y-1][2])<K13Endlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y-1][2])
                            #if K13Startlist1[z]<int(depthlist1[y][2])<K13Endlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y][2])
                            #if K13Startlist1[z]<int(depthlist1[y+1][2])<K13Endlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y+1][2])
                        codondepthlist1+=[int(tempcodondepthlist1/3)]
                        break
                    if Genelist1[x] == "PfMDR1":
                        #for z in range(len(PfMDR1Startlist1)):
                            #if PfMDR1Startlist1[z]<int(depthlist1[y-1][2])<PfMDR1Endlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y-1][2])
                            #if PfMDR1Startlist1[z]<int(depthlist1[y][2])<PfMDR1Endlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y][2])
                            #if PfMDR1Startlist1[z]<int(depthlist1[y+1][2])<PfMDR1Endlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y+1][2]) 
                        codondepthlist1+=[int(tempcodondepthlist1/3)]
                    if Genelist1[x] == "DHPS":
                        #for z in range(len(PfDHPSStartlist1)):
                            #if PfDHPSStartlist1[z]<int(depthlist1[y-1][2])<PfDHPSEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y-1][2])
                            #if PfDHPSStartlist1[z]<int(depthlist1[y][2])<PfDHPSEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y][2])
                            #if PfDHPSStartlist1[z]<int(depthlist1[y+1][2])<PfDHPSEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y+1][2])
                        codondepthlist1+=[int(tempcodondepthlist1/3)]
                        break
                    if Genelist1[x] == "mitochondrial_genome":
                        #for z in range(len(mitochondrialgenomeStartlist1)):
                            #if mitochondrialgenomeStartlist1[z]<int(depthlist1[y-1][2])<mitochondrialgenomeEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y-1][2])
                            #if mitochondrialgenomeStartlist1[z]<int(depthlist1[y][2])<mitochondrialgenomeEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y][2])
                            #if mitochondrialgenomeStartlist1[z]<int(depthlist1[y+1][2])<mitochondrialgenomeEndlist1[z]:
                        tempcodondepthlist1+=int(depthlist1[y+1][2])
                        codondepthlist1+=[int(tempcodondepthlist1/3)]
                        break
                #elif Genelist1[x] == "DHFR": 
                #elif Genelist1[x] == "K13": 
                #elif Genelist1[x] == "PfMDR1": 
                    #break
    #print(len(Genelist1))
    #print(len(codondepthlist1))



#PfCRT   1   3399    PfCRT   .   +   96  3191    0   13  90,268,172,132,71,75,82,50,56,92,44,54,76,  96,364,812,1157,1443,1638,1810,2020,2208,2413,2699,2891,3115,
#MT  1   5967    COXIII  .   -   734 1573    0   1   839,    734,
#MT  1   5967    COL .   +   1933    3471    0   1   1538,   1933,
#MT  1   5967    CYTOb   .   +   3492    4622    0   1   1130,   3492,
#PfDHFR  1   1827    PfDHFR  .   +   1   1827    0   1   1827,   1,
#PfDHPS  1   2417    PfDHPS  .   +   1   2417    0   3   135,1868,115,   1,313,2302,
#PfK13   1   2181    PfK13   .   +   1   2181    0   1   2181,   1,
#PfMDR1  1   4260    PfMDR1  .   +   1   4260    0   1   4260,   1,

    
    #print(len(Genelist1))
    #print(len(POSlist1))
    #print(len(depthlist1up))
    #print(depthlist1up)
    #with open(args.unfiltered, "r") as f2:
    #    for line in f2:
    #        print(line)
   # print(len(Genelist1))
    #print(POSlist1)
    #print(POSlist1)
    #print(Genelist1)
Genelist2=[]
POSlist2=[]
AFlist1=[]
QDlist1=[]
SORlist1=[]
MQlist1=[]
MQRankSumlist1=[]
DP4list1=[]
f2= open(args.unfiltered, "r")
f3= open(args.filtered, "r")
wholefilnumsx=[]
wholefilnums=[]


DepthVal={}
for line in f3:
    count=0
    tempword=[]
    tempsumall=0
    tempsumalt=0
    tempgeneval1=""
    tempposval1=""
    for word in line.split():
        count+=1
        if count==8:
            for x in word.split(","):
                tempword+=[x]
        if count==1:
            tempgeneval1=word
        if count==2:
            tempposval1=word
    if len(tempword)==4:
        for x in range(len(tempword)):
            tempsumall+=int(tempword[x])
            if x>1:
                tempsumalt+=int(tempword[x])
        DepthVal[tempgeneval1,tempposval1]=[tempsumalt/tempsumall]
    if len(tempword)==2:
        for x in range(len(tempword)):
            tempsumall+=int(tempword[x])
            if x>0:
                tempsumalt+=int(tempword[x])
        DepthVal[tempgeneval1,tempposval1]=[tempsumalt/tempsumall]

#print(DepthVal)
#print(tempVal)



    #if len(tempword)==2:
    #    print(tempword)



for line in f3: 
    count=0
    tempnum=""
    #print(line)
    if line.find("p.")!=-1 and line.startswith("#")==False:
        for word in line.split():
            count+=1
            if count==2:
                tempnum=word
            #print(tempnum)
                #print(word)
            #print(tempnum)
            if word.startswith("c."):
                if word.find("*")!=-1:
                    wholefilnumsx+=[tempnum]
                    #print(word[word.find("c.*"):word.find(tempnum)+len(tempnum)])
                    wholefilnums+=[word[word.find("c.*"):word.find(tempnum)+len(tempnum)]]
                elif word.find("*")==-1:
                    #print(word[word.find("c.*"):word.find(tempnum)+len(tempnum)])
                    wholefilnumsx+=[tempnum]
                    wholefilnums+=[word[word.find("c."):word.find(tempnum)+len(tempnum)]]
            if word.startswith("p."):
                #print(word[5+len(temppos):5+len(temppos)+3])
                #print(word[5+len(temppos):5+len(temppos)+3])
                temppos=""
                for x in word:
                    #print(word)
                    if x.isdigit():
                        #print(x)
                        temppos+=x
                if word[2:5] not in AAdic or word[5+len(temppos):5+len(temppos)+3] not in AAdic:
                    #print("working")
                    if word[2:5] in AAdic and word[-1:] == "*":
                        break
                    wholefilnumsx=wholefilnumsx[:-1]
                    wholefilnums=wholefilnums[:-1]
#print(wholefilnumsx)
        #print(tempnum)
#print(wholefilnums)
#print(len(Genelist1))
#print(len(wholefilnums))
#print(wholefilnums)
#print(wholefilnumsx)
#print(wholefilnums)
#print(wholefilnums)
#print(wholefilnumsx)
#print(Genelist1)
#print(POSlist1)
Geneset1=[]
for k in range(len(Genelist1)):
    Geneset1+=[[Genelist1[k],POSlist1[k]]]
GenePosSet1=[]
for line in f2:
    count=0
    tempDP4list=[]
    for word in line.split():
        #print(line)
        #print(Genelist1[x])
        #print(word)
        #if count==0 and word == Genelist1[x]:
            #print(word+"test1111")
        if count==0:
            word1=""
            if word.startswith("#")==False:
                word1=word
                #print(word1)
        #print(word1)
        #print(count)
        #print(word1)
        #temppos3=""
        #tempgene1=""
        if count==1:
            #print(word1)
            #print(word1)
            if word1!="":
                #filtercount=-1
                #print(word1)
                #print(word)
                #print(POSlist1)
                #print(Genelist1)
                #print(word1)
                for x in range(len(Genelist1)):
                    #print(str(int(POSlist1[x])-1))
                    #if word1== Genelist1[x]:print(word1,word)
                    #print(POSlist1)
                    #print(word)
                    #if word.startswith("c."):
                        #temppos2=""
                        #for x in word:
                            #print(word)
                            #if x.isdigit():
                                #print(x)
                                #temppos2+=x
                        #print(temppos2)
                    #print(word1)
                    #print(Genelist1[x])
                    #print(POSlist1[x])
                    #print(word)
                    #print(Genelist1[x],POSlist1[x])
                    #print(word1,word)
                    #print(POSlist1)
                    #print(type(POSlist1[x]))
                    #if word1 == Genelist1[x]:print("no")
                    #print(word)
                    #print(POSlist1[x])
                    #if word == POSlist1[x]: print("yes")
                    if word1== Genelist1[x] and word== POSlist1[x] and line.startswith("#")==False:
                        #print(word)
                        #print(word)
                        #print(word)
                        #filtercount+=1
                        #print(POSlist1[x])
                        #print("test"+POSlist1[x])
                        #print(Genelist1[x],POSlist1[x])
                        Genelist2+=[Genelist1[x]]
                        POSlist2+=[POSlist1[x]]
                        GenePosSet1+=[[Genelist1[x],POSlist1[x]]]
                        #temppos3=POSlist1[x]
                        if line.find(";AF=") != -1:
                            #print(line[int(line.find("AF="))+3:len(line)].find(";"))
                            #print(line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(";")])
                            #print(line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(",")].isdigit())
                            #tempAF=""
                            countAF=0
                            for k in (line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(",")]):
                                if k.isdigit() == False and k !=".":
                                    AFlist1+=[line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(";")]]
                                    break
                                #print(len(line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(",")]))
                                #print("COUNTAF"+str(countAF))
                                if countAF == len(line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(",")])-1:
                                    AFlist1+=[line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(",")]]
                                    break
                                countAF+=1
                            #else:AFlist1+=[line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(";")]]
                            #for k in range(int(line.find("AF="))+3,line[int(line.find("AF="))+3:len(line)].find(";")):
                            #    if line[k]==";":
                            #        break
                            #    if line[line.find("AF=")+3:k].find(",") == -1:
                            #        AFlist1+=[float(line[line.find("AF=")+3:k])]
                            #    elif line[line.find("AF=")+3:k].find(",") != -1:
                            #        AFlist1+=[float(line[line.find("AF=")+3:line.find("AF=")+5])]
                        else: AFlist1+=["None"]
                        #print(AFlist1)
                        if line.find("MQ=") != -1:
                        #if line.find("MQ=")+3 != 2:
                        #    for k in range(int(line.find("MQ=")+3),len(line)):
                        #        if line[k]==";":
                        #            break
                        #        MQlist1+=[float(line[line.find("MQ=")+3:k])]
                            MQlist1+=[line[int(line.find("MQ="))+3:int(line.find("MQ="))+3+line[int(line.find("MQ="))+3:len(line)].find(";")]]
                        else: MQlist1+=["None"]
                        if line.find("MQRankSum=") != -1:
                        #    for k in range(int(line.find("MQRankSum="))+10,len(line)):
                        #        if line[k]==";":
                        #           break
                        #        MQRankSumlist1+=[float(line[line.find("MQRankSum=")+10:k])]
                            #print(int(line.find("MQRankSum="))+10:int(line.find("MQRankSum="))+10+line[int(line.find("MQRankSum="))+10])
                            #print(int(line.find("MQRankSum="))+10+line[int(line.find("MQRankSum="))+10:len(line)].find(";"))
                            #print(line[int(line.find("MQRankSum="))+10:int(line.find("MQRankSum="))+10+line[int(line.find("MQRankSum="))+10:len(line)].find(";")])
                            MQRankSumlist1+=[line[int(line.find("MQRankSum="))+10:int(line.find("MQRankSum="))+10+line[int(line.find("MQRankSum="))+10:len(line)].find(";")]]
                        else:MQRankSumlist1+=["None"]
                        if line.find("SOR=") != -1:
                            #for k in range(int(line.find("SOR="))+4,len(line)):
                            #    if line[k]==";":
                            #        break
                            #    SORlist1+=[float(line[line.find("SOR=")+4:k])]
                            SORlist1+=[line[int(line.find("SOR="))+4:int(line.find("SOR="))+4+line[int(line.find("SOR="))+4:len(line)].find(";")]]
                        else: SORlist1+=["None"]
                        if line.find("QD=") != -1:
                            #for k in range(int(line.find("QD="))+3,len(line)):
                            #    if line[k]==";":
                            #        break
                            #    QDlist1+=[float(line[line.find("QD=")+3:k])]
                            QDlist1+=[line[int(line.find("QD="))+3:int(line.find("QD="))+3+line[int(line.find("QD="))+3:len(line)].find(";")]]
                        else: QDlist1+=["None"]
                        #print(filternumber)
                        DP4list1+=DepthVal[word1,word]
                        break
        #print(len(QDlist1))
        #print(GenePosSet1)
        #print(AFlist1)
        if count==7:
            if word1 != "":
                #print(word)
                if line.startswith("#")==False:
                    if word.startswith("AC=") or word.startswith("AB="):
                        #print(word[word.find("c."):word.find("c.")+10])
                        #print(line)
                        for x in range(len(wholefilnums)):
                            #print(x)
                            #print(wholefilnumsx[count])
                            #print(POSlist2)
                            #print(wholefilnumsx[count2])
                            #print([word1,str(int(wholefilnumsx[x]))])
                            #print(word)
                            #find part of c.138 in c.1381 in different gene.....
                            #print(line.find(wholefilnums[x]))
                            #if word.find(wholefilnums[x]) != -1:
                                #print(wholefilnums[x])
                                #print(word.find(wholefilnums[x]))
                                #print(len(wholefilnums[x]))
                                #print(wholefilnums[x])
                                #print(word[word.find(wholefilnums[x])+len(wholefilnums[x])])
                                #print(word[word.find(wholefilnums[x])+len(wholefilnums[x])-len(wholefilnumsx[x])-1])
                            if word.find(wholefilnums[x])!=-1 and str(word.find(wholefilnums[x])+len(wholefilnums[x])+1).isdigit()==False and [word1,wholefilnumsx[x]] in Geneset1 and [word1,str(int(wholefilnumsx[x]))] not in GenePosSet1:
                                #print(word1)
                                #print(wholefilnumsx[x])
                                #print(wholefilnums[x])
                                #print(wholefilnumsx[x])
                                #print(word1)
                                #print(word1)
                                #print(wholefilnumsx[x])
                                Genelist2+=[word1]
                                #print(wholefilnumsx[x])
                                POSlist2+=[wholefilnumsx[x]]
                                #"QD": QDlist1, "SOR":SORlist1, "MQ":MQlist1, "MQRankSum":MQRankSumlist1,
                                AFlist1+=[AFlist1[-1]]
                                QDlist1+=[QDlist1[-1]]
                                SORlist1+=[SORlist1[-1]]
                                MQlist1+=[MQlist1[-1]]
                                MQRankSumlist1+=[MQRankSumlist1[-1]]
                                GenePosSet1+=[[word1,str(int(wholefilnumsx[x]))]]
                                DP4list1+=[DP4list1[-1]]
                            #print(Genelist1[count])
                            #print(wholefilnumsx[count])
        count+=1
#print(len(wholefilnums))
#print(len(Genelist2))
#print(len(QDlist1))
#print(len(Genelist2))
#print(POSlist1)
#print(POSlist2)
#print(Genelist1)
#print(Genelist2)
#print(DP4list1[2])
Descriptionlist1=[]
filterscorelist1=[]

for x in range(len(QDlist1)):
    filternumber=0
    Description=""
    #print(QDlist1[filtercount])
    if QDlist1[x] != "None":
        filternumber+=float(QDlist1[x])/80
        if float(SORlist1[x])>30:
            Description+="Good,"
        else:
            Description+="Bad,"
    else:Description+="NA,"
    if SORlist1[x] != "None":
        if float(SORlist1[x])<3:
            Description+="Good,"
        else:
            Description+="Bad,"
    else:Description+="NA,"
    if MQlist1[x] != "None":
        if float(MQlist1[x])<70:
            filternumber+=float(MQlist1[x])/140
        else:filternumber+=1/2
        if 40<float(MQlist1[x])<60:
            Description+="Good,"
        else:
            Description+="Bad,"
    else:Description+="NA,"
    #print(MQRankSumlist1)
    if MQRankSumlist1[x] != "None":
        if -5<abs(float(MQRankSumlist1[x]))<5:
            Description+="Good"
        else:
            Description+="Bad"
    else:Description+="NA"
    Descriptionlist1+=[Description]
    #print(Description)
    filterscorelist1+=[filternumber]
#print(filterlist1)
                    #print(word+"test2222")
#for x in range(len(Genelist1)):
    #print(Genelist1[x])
#    for line in f2:
        #count=0
        #print(line)
#        for word in line.split():
            #if count==0:
                #print(word + "test")
#            if Genelist1[x] == word:
#                print(Genelist1[x] +"test222")
                #    print(line)
                #    break
        #    count+=1
            #if count==1:
            #    if POSlist1[x] == word:
            #        print(line)

#print(len(Genelist1))
#print(len(depthlist1up))
#(sadcount)
#print(len(AAPoslist1))
#print(len(codondepthlist1))
Reflist2=[]
Altlist2=[]
AAlist21=[]
AAlist22=[]
depthlist2up=[]
AAPoslist12=[]
codondepthlist2=[]
#print(Genelist2)
#print(POSlist2)
#print(Genelist1)
#print(POSlist1)
#print(len(Genelist2))
#print(len(POSlist2))
#print(len(AAPoslist12))
#print(len(AAlist21))
#print(len(AAlist22))
#print(len(depthlist1up))
#print(len(Genelist1))
#print(len(AAlist2))
#print(len(Reflist1))
#print(len(Altlist1))
#print(len(Reflist1))
#print(len(AAlist1))
#print(len(AAlist2))
#print(Altlist1)

count=0
for x in range(len(Genelist2)):
    for y in range(len(Genelist1)):
        if  Genelist2[x]==Genelist1[y] and POSlist2[x]==POSlist1[y]:
            count+=1
            Reflist2+=[Reflist1[y]]
            Altlist2+=[Altlist1[y]]
            depthlist2up+=[depthlist1up[y]]
            codondepthlist2+=[codondepthlist1[y]]
            break
    #Geneset2=[]
    #for k in range(len(Genelist1)):
    #    Geneset2+=[[Genelist1[k],POSlist1[k]]]
    #for x in range(len(Genelist2)):
    #    for y in range(len(Genelist1)):
    #        if [Genelist2[x],POSlist2[x]] not in Geneset2:
    #            print(Genelist2[x],POSlist2[x])
    #for y in POSlist2:
    #    if y not in POSlist1:
    #        print(y)
    #print(len(Genelist2))
    #print(count)
    #print(len(wholefilnumsx))
    #print(len(wholefilnums))
    #        print("order", x, wholelenlist[x], len(Genelist2))
xls=pd.ExcelFile(args.candidates)
df1=pd.read_excel(xls)
df1Gene=df1['Gene'].tolist()
df1RefAA=df1['RefAA'].tolist()
df1AltAA=df1['AltAA'].tolist()
df1AAPos=df1['AAPos'].tolist()
#df1.to_csv("test1.csv", index=False)
df2 = pd.read_csv(args.voi)
df2Gene=df2['Gene'].tolist()
df2RefAA=df2['RefAA'].tolist()
df2AltAA=df2['AltAA'].tolist()
df2AAPos=df2['AAPos'].tolist()
    #print(df2RefAA)
    #print(df2Gene)
    #print(df2RefAA)
 
filteredGenelist2=[]
filteredPoslist2=[]
filtereddepthlist2up=[]
filteredReflist2=[]
filteredAltlist2=[]
filteredAAlist21=[]
filteredAAlist22=[]
filteredAAPoslist12=[]
filteredcodondepthlist2=[]
filteredAFlist1=[]
filteredmutationlist1=[]
filteredQDlist1=[]
filteredSORlist1=[]
filteredMQlist1=[]
filteredMQRankSumlist1=[]
filteredfilterscorelist1=[]
filteredDescriptionlist1=[]
filteredcandidateslist1=[]
filteredreportablelist1=[]


combineddflist1=[]
combineddflist2=[]
for x in range(len(df1Gene)):
    if df1Gene[x] == "Pfcrt":
        combineddflist1+=[["PfCRT",str(df1AAPos[x])]]
    elif df1Gene[x] == "Pfdhfr":
        combineddflist1+=[["DHFR",str(df1AAPos[x])]]
    elif df1Gene[x] == "PfK13":
        combineddflist1+=[["K13",str(df1AAPos[x])]]
    elif df1Gene[x] == "Pfmdr1":
        combineddflist1+=[["PfMDR1",str(df1AAPos[x])]]
    elif df1Gene[x] == "Pfdhps":
        combineddflist1+=[["DHPS",str(df1AAPos[x])]]
    elif df1Gene[x] == "MT":
        combineddflist1+=[["mitochondrial_genome",str(df1AAPos[x])]]
for x in range(len(df2Gene)):
    if df2Gene[x] == "Pfcrt":
        combineddflist2+=[["PfCRT",str(df2AAPos[x])]]
    elif df2Gene[x] == "Pfdhfr":
        combineddflist2+=[["DHFR",str(df2AAPos[x])]]
    elif df2Gene[x] == "PfK13":
        combineddflist2+=[["K13",str(df2AAPos[x])]]
    elif df2Gene[x] == "Pfmdr1":
        combineddflist2+=[["PfMDR1",str(df2AAPos[x])]]
    elif df2Gene[x] == "Pfdhps":
        combineddflist2+=[["DHPS",str(df2AAPos[x])]]
    elif df2Gene[x] == "MT":
        combineddflist2+=[["mitochondrial_genome",str(df2AAPos[x])]]

#print(combineddflist2)


bedname1=[]
bedstart1=[]
bedend1=[]
tempstart1=[]
tempend1=[]
with open(args.bed, "r") as fb1:
    for line in fb1:
        count=0
        #print(len(line.split()))
        if line.split()[0]==line.split()[3] or line.split()[3]=="CYTOb":
            for word in line.split():
                count+=1
                if count==1:
                    bedname1+=[word]
                if count==len(line.split()):
                    tempstart1=[]
                    for x in word[:-1].split(","):
                        #print(x)
                        tempstart1+=[int(x)]
                    bedstart1+=[tempstart1]
                    #tempstart1=[word[:-1]]
                    #print(word[:-1])
                if count==len(line.split())-1:
                    #print(tempstart1)
                    #print(word[:-1])
                    #print([x+y for x,y in zip(tempstart1,[word[:-1]])])
                    tempend1=[]
                    for x in word[:-1].split(","):
                        #print(x)
                        tempend1+=[int(x)]
            bedend1+=[[x+y for x,y in zip(tempstart1,tempend1)]]
        #print(bedend1)
#print(bedstart1)
#print(bedend1)
#print(bedname1)
dictrange={}
#print(bedstart1)
#print(bedname1)
#print(bedend1)
for x in range(len(bedname1)):
    #print(bedstart1[x])
    exec(bedname1[x]+"Startlist1"+"="+str(bedstart1[x]))
    #print(bedname1[x])
    #print(PfCRTStartlist1)
    exec(bedname1[x]+"Endlist1"+"="+str(bedend1[x]))
    #bedname1[x]+"Startlist1"=[bedstart1[x]]
    #bedname1[x]+"Endlist1"=bedend1[x]
    dictrange[bedname1[x]]=[bedstart1[x],bedend1[x]]

#print(dictrange)
#for x in range(len(bedname1)):
#print(PfDHFREndlist1)''
#print(MTEndlist1)
#print(dictrange)
dicttrans1={}
dicttrans2={}

with open(args.fasta, "r") as fa1:
    for line in fa1:
        if line.find(">")!=-1:
            for x in range(len(bedname1)):
                #print(line)
                if line[1::].strip()==bedname1[x]:
                    tempname=bedname1[x]
        if line.find(">")==-1:
            dicttrans1[tempname]=line.strip()

#print(Genelist2)
#print(POSlist2)
#print(Altlist2)
GenePosSet2=[]
GenePosdict={}
for x in range(len(Genelist2)):
    GenePosSet2+=[[Genelist2[x],POSlist2[x]]]
    GenePosdict[Genelist2[x],POSlist2[x]]=Altlist2[x]
#print(Reflist2)
#print(GenePosSet2)

#############################################
###############no mutated base applied#######
#############################################
testdic={}
for x in range(len(bedname1)):
    tempfa1=""
    for y in dicttrans1[bedname1[x]]:
        tempfa1+=y
    testdic[bedname1[x]]=tempfa1
#print(dictrange)

#print(dicttrans1["PfCRT"])
#print(GenePosdict)
#############################################
###############mutated base applied##########
#############################################
#print(GenePosSet2)
for x in range(len(bedname1)):
    tempfa1=""
    count=0
    for y in dicttrans1[bedname1[x]]:
        count+=1
        #print([bedname1[x],str(count)] in GenePosSet2)
        if ([bedname1[x],str(count)]) in GenePosSet2:
            #print("True")
            #if bedname1[x]=="PfCRT":
                #print(count)
                #print(GenePosdict[bedname1[x],str(count)])
            tempfa1+=GenePosdict[bedname1[x],str(count)]
        else:
            tempfa1+=y
    dicttrans2[bedname1[x]]=tempfa1
#print(dictrange)
dicttrans22={}
testdicttrans22={}
#############################################
###############cds region applied############
#############################################

#print(1%3)
#print(dictrange)
#for x in range(3):
#    print(x)

#print(dictrange)
#print(1 in dictrange["K13"][0])
BASEAAPOSdic1={}
for x in range(len(bedname1)):
    tempfa2=""
    temptestfa2=""
    tempAAPOS=""
    count=0
    AAcount=0
    for y in range(len(dictrange[bedname1[x]][0])):
        #print(dictrange[bedname1[x][0]])
        if len(dictrange[bedname1[x]][0])==1:
            for z in range(dictrange[bedname1[x]][0][y]-1,dictrange[bedname1[x]][1][y]-1):
                #print(z)
                if count%3==1:
                    AAcount+=1
                count+=1
                if ([bedname1[x],str(z)]) in GenePosSet2: 
                    BASEAAPOSdic1[bedname1[x],str(z)]=AAcount
                tempfa2+=dicttrans2[bedname1[x]][z]
                temptestfa2+=testdic[bedname1[x]][z]
        if len(dictrange[bedname1[x]][0])>1:
            for z in range(dictrange[bedname1[x]][0][y]-1,dictrange[bedname1[x]][1][y]):
                #print(z)
                if count%3==1:
                    #print(count)
                    AAcount+=1
                count+=1
                if ([bedname1[x],str(z)]) in GenePosSet2: 
                    #print(z)
                    #print(AAcount)
                    BASEAAPOSdic1[bedname1[x],str(z)]=AAcount
                tempfa2+=dicttrans2[bedname1[x]][z]
                temptestfa2+=testdic[bedname1[x]][z]
    #print(len(tempfa2))
    dicttrans22[bedname1[x]]=tempfa2
    testdicttrans22[bedname1[x]]=temptestfa2

#print(dicttrans2["PfCRT"])
#print(dicttrans22["DHPS"])
#print(testdicttrans22["DHPS"])
#print(BASEAAPOSdic1)
#for x in range(len(dicttrans22["PfCRT"])):
#    if dicttrans22["PfCRT"][x]!=testdicttrans22["PfCRT"][x]:
#        print(x)
#        print(dicttrans22["PfCRT"][x])
#print(BASEBASEAAPOSdic1)
#print(dicttrans22["PfCRT"])
#print(testdicttrans22["DHPS"])
#print((testdicttrans22["DHPS"]))
#############################################
###############translation applied############
#############################################
#print(dicttrans22["PfCRT"])

dicttransAA1={}
testdicttransAA1={}
for x in range(len(bedname1)):
    tempAA1=""
    testtempAA1=""
    count=0
    #print(bedname1[x])
    for z in range(len(dicttrans22[bedname1[x]])-2):
        count+=1
        if count==1:
            #print(line.strip()[z:z+3])
            #print(bedname1[x])
            #print(dicttrans22[bedname1[x]][z:z+3])
            tempAA1+=AAdic[codondic[dicttrans22[bedname1[x]][z:z+3]]]
        if count==3:
            count=0
    dicttransAA1[bedname1[x]]=tempAA1
    count=0
    for k in range(len(testdicttrans22[bedname1[x]])-2):
        count+=1
        if count==1:
            #if bedname1[x]=="PfMDR1":
                #print(testdicttrans22[bedname1[x]][k:k+3])
            testtempAA1+=AAdic[codondic[testdicttrans22[bedname1[x]][k:k+3]]]
        if count==3:
            count=0
    testdicttransAA1[bedname1[x]]=testtempAA1




#print(dicttrans22["PfCRT"])
#print(testdicttrans22["PfCRT"])
#print(dicttransAA1["DHPS"])
#print(testdicttransAA1["DHPS"])
#print(POSlist1)
#print(Genelist2)
#print(POSlist2)
#print(customAA1)
#d3={"Gene":Genelist2, "BasePOS":POSlist2, "customAA":customAA1}
#df3=pd.DataFrame(data=d3)
#df3.to_csv('custom1.csv', index=False,) 
#print(len(Genelist2))
#print(len(Genelist1))
#print(len(AAPoslist12))
#print(len(AAlist21))
#print(len(AAlist22))
#print(Genelist2)
wholepotlist1=[]
for k in range(len(df1Gene)):
    if df1Gene[k] == "Pfcrt":
        wholepotlist1.append(["PfCRT",str(df1AAPos[k])])
    elif df1Gene[k] == "Pfdhfr":
        wholepotlist1.append(["DHFR",str(df1AAPos[k])])
    elif df1Gene[k] == "PfK13":
        wholepotlist1.append(["K13",str(df1AAPos[k])])
    elif df1Gene[k] == "Pfmdr1":
        wholepotlist1.append(["PfMDR1",str(df1AAPos[k])])
    elif df1Gene[k] == "Pfdhps":
        wholepotlist1.append(["DHPS",str(df1AAPos[k])])
    elif df1Gene[k] == "MT":
        wholepotlist1.append(["mitochondrial_genome",str(df1AAPos[k])])
for k in range(len(df2Gene)):
    if df2Gene[k] == "Pfcrt":
        wholepotlist1.append(["PfCRT",str(df2AAPos[k])])
    elif df2Gene[k] == "Pfdhfr":
        wholepotlist1.append(["DHFR",str(df2AAPos[k])])
    elif df2Gene[k] == "PfK13":
        wholepotlist1.append(["K13",str(df2AAPos[k])])
    elif df2Gene[k] == "Pfmdr1":
        wholepotlist1.append(["PfMDR1",str(df2AAPos[k])])
    elif df2Gene[k] == "Pfdhps":
        wholepotlist1.append(["DHPS",str(df2AAPos[k])])
    elif df2Gene[k] == "MT":
        wholepotlist1.append(["mitochondrial_genome",str(df2AAPos[k])])


    #print(x[0])
    #print(x[1])
    #print(AAuntranslatewilddic1[x[0],x[1]])
    #print(AAtranslatewilddic1[x[0],x[1]])  

totalgene=[]
with open(args.fasta, "r") as f1:
    for line in f1:
        if line.startswith(">"):
            totalgene+=[(line.strip(">")).strip("\n")]

#print(totalgene)


#samtools index ${bam_path}
result=subprocess.Popen(['samtools', 'index', args.name+ "_SR.bam"])
process=subprocess.Popen(['samtools', 'coverage', args.name+ "_SR.bam", '-r', 'mitochondrial_genome:3492-4622'], stdout=subprocess.PIPE)
stdout, stderr = process.communicate()
reader2 = csv.reader(stdout.decode('ascii').splitlines(),
                        delimiter='\t', skipinitialspace=True)

tempcov1=[]
tempcovwordlist1=[]
for line in reader2:
    count=0
    tempcovword1=""
    #print(line)
    if line[0].startswith("#")==False:
        for word in line:
            count+=1
            if count==1:
                tempcovword1=word
            #print(word)
            if count==6 and float(word)>75:
                #print(word)
                tempcov1+=[word]
                tempcovwordlist1+=[tempcovword1]

process=subprocess.Popen(['samtools', 'coverage', args.name+ "_SR.bam"], stdout=subprocess.PIPE)
stdout, stderr = process.communicate()
reader2 = csv.reader(stdout.decode('ascii').splitlines(),
                        delimiter='\t', skipinitialspace=True)

for line in reader2:
    count=0
    tempcovword1=""
    #print(line)
    if line[0].startswith("#")==False:
        for word in line:
            count+=1
            if count==1:
                tempcovword1=word
            #print(word)
            if count==6 and float(word)>75 and tempcovword1 not in tempcovwordlist1:
                #print(word)
                tempcov1+=[word]
                tempcovwordlist1+=[tempcovword1]




#print(tempcov1)
#print(tempcovwordlist1)
#print(totalgene2)

#d3={"Gene":wildgene1, "AAPOS":wildgeneAApos1, "previousAA": wildAAprev1, "customAA":wildAAtrans1, "type":wildtype1}
#df3=pd.DataFrame(data=d3)
#df3.to_csv('wildtype.csv', index=False,) 
#print(str(codoncount)ra

filteredGenelist2=[]
filteredPoslist2=[]
filtereddepthlist2up=[]
filteredReflist2=[]
filteredAltlist2=[]
filteredAAlist21=[]
filteredAAlist22=[]
filteredAAPoslist12=[]
filteredcodondepthlist2=[]
filteredAFlist1=[]
filteredmutationlist1=[]
filteredQDlist1=[]
filteredSORlist1=[]
filteredMQlist1=[]
filteredMQRankSumlist1=[]
filteredfilterscorelist1=[]
filteredDescriptionlist1=[]
filteredcandidateslist1=[]
filteredreportablelist1=[]

#print(len(Genelist2))
#print(Genelist2)
#print(POSlist2)
#print(BASEAAPOSdic1)
combinedGenePOSlist1=[]
for x in range(len(Genelist2)):
    if (Genelist2[x],(POSlist2[x])) in BASEAAPOSdic1:
        #print("True")
        #print([Genelist2[x],BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]])
        #print(BASEAAPOSdic1[Genelist2[x],(POSlist2[x])])
        #print(Genelist2[x],POSlist2[x])
        #print([Genelist2[x],str(BASEAAPOSdic1[Genelist2[x],(POSlist2[x])])])
        if [Genelist2[x],str(BASEAAPOSdic1[Genelist2[x],(POSlist2[x])])] in combineddflist1 or [Genelist2[x],str(BASEAAPOSdic1[Genelist2[x],(POSlist2[x])])] in combineddflist2:
            #print("True")
            combinedGenePOSlist1+=[[Genelist2[x],BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]]]

#print(combinedGenePOSlist1)
#print(combinedGenePOSlist1)
for x in range(len(Genelist2)):
    #print(Genelist2[x],(POSlist2[x]))
    if (Genelist2[x],(POSlist2[x])) in BASEAAPOSdic1:
        #print([Genelist2[x],BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]])
        if [Genelist2[x],str(BASEAAPOSdic1[Genelist2[x],(POSlist2[x])])] in combineddflist1 or [Genelist2[x],str(BASEAAPOSdic1[Genelist2[x],(POSlist2[x])])] in combineddflist2:
    #print(AAPoslist12[x])
            #print("True")
            filteredGenelist2+=[Genelist2[x]]
            filteredPoslist2+=[POSlist2[x]]
            filtereddepthlist2up+=[depthlist2up[x]]
            filteredReflist2+=[Reflist2[x]]
            filteredAltlist2+=[Altlist2[x]]
            #print(BASEAAPOSdic1[Genelist2[x],(POSlist2[x])])
            #print(testdicttransAA1[Genelist2[x])
            #print(dicttransAA1[Genelist2[x]])
            filteredAAlist21+=[testdicttransAA1[Genelist2[x]][BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]-1]]
            filteredAAlist22+=[dicttransAA1[Genelist2[x]][BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]-1]]
            #print(AAuntranslatewilddic1[Genelist2[x],AAPoslist12[x]])
            #print(AAtranslatewilddic1[Genelist2[x],AAPoslist12[x]])
            filteredAAPoslist12+=[BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]]
            filteredcodondepthlist2+=[codondepthlist2[x]]
            filteredAFlist1+=[AFlist1[x]]
            if testdicttransAA1[Genelist2[x]][BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]-1]==dicttransAA1[Genelist2[x]][BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]-1]:
                #print("False")
                filteredmutationlist1+=["Wildtype"]
            if DP4list1[x]!="None":
                if float(DP4list1[x])>=0.5:
                    if testdicttransAA1[Genelist2[x]][BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]-1]!=dicttransAA1[Genelist2[x]][BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]-1]:
                        filteredmutationlist1+=["Major"]
                if float(DP4list1[x])<0.5:
                    if testdicttransAA1[Genelist2[x]][BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]-1]!=dicttransAA1[Genelist2[x]][BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]-1]:
                        filteredmutationlist1+=["Minor"]
            if DP4list1[x]=="None":
                #print('TRUE')
                if testdicttransAA1[Genelist2[x]][BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]-1]!=dicttransAA1[Genelist2[x]][BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]-1]:
                    filteredmutationlist1+=["Wildtype"]
            filteredQDlist1+=[QDlist1[x]]
            filteredSORlist1+=[SORlist1[x]]
            filteredMQlist1+=[MQlist1[x]]
            filteredMQRankSumlist1+=[MQRankSumlist1[x]]
            filteredfilterscorelist1+=[filterscorelist1[x]]
            filteredDescriptionlist1+=[Descriptionlist1[x]]



filteredGenelist22=[]
filteredPoslist22=[]
filtereddepthlist2up2=[]
filteredReflist22=[]
filteredAltlist22=[]
filteredAAlist212=[]
filteredAAlist222=[]
filteredAAPoslist122=[]
filteredcodondepthlist22=[]
filteredAFlist12=[]
filteredmutationlist12=[]
filteredQDlist12=[]
filteredSORlist12=[]
filteredMQlist12=[]
filteredMQRankSumlist12=[]
filteredfilterscorelist12=[]
filteredDescriptionlist12=[]
filteredcandidateslist12=[]
filteredreportablelist12=[]


combinedwildgeneAA=[]

#print(wholepotlist1)
#print(combinedGenePOSlist1)
#print(combinedGenePOSlist1)
for y in wholepotlist1:
    #count=0
    for x in range(len(combinedGenePOSlist1)):
        #count+=1
        #print(x)
        #print(len(combinedGenePOSlist1))
        if [y[0],int(y[1])] == combinedGenePOSlist1[x]:
            #print(combinedGenePOSlist1[x])
            #print("True")
            filteredGenelist22+=[filteredGenelist2[x]]
            filteredPoslist22+=[filteredPoslist2[x]]
            filtereddepthlist2up2+=[filtereddepthlist2up[x]]
            filteredReflist22+=[filteredReflist2[x]]
            filteredAltlist22+=[filteredAltlist2[x]]
            filteredAAlist212+=[filteredAAlist21[x]]
            filteredAAlist222+=[filteredAAlist22[x]]
            filteredAAPoslist122+=[filteredAAPoslist12[x]]
            filteredcodondepthlist22+=[filteredcodondepthlist2[x]]
            filteredAFlist12+=[filteredAFlist1[x]]
            filteredmutationlist12+=[filteredmutationlist1[x]]
            filteredQDlist12+=[filteredQDlist1[x]]
            filteredSORlist12+=[filteredSORlist1[x]]
            filteredMQlist12+=[filteredMQlist1[x]]
            filteredMQRankSumlist12+=[filteredMQRankSumlist1[x]]
            filteredfilterscorelist12+=[filteredfilterscorelist1[x]]
            filteredDescriptionlist12+=[filteredDescriptionlist1[x]]
        if x==len(combinedGenePOSlist1)-1 and ([y[0],int(y[1])] not in combinedGenePOSlist1) and ([y[0],int(y[1])] not in combinedwildgeneAA) and (y[0] in tempcovwordlist1):
            combinedwildgeneAA+=[[y[0],int(y[1])]]
            filteredGenelist22+=[y[0]]
            filteredPoslist22+=["NA"]
            filtereddepthlist2up2+=["NA"]
            filteredReflist22+=["NA"]
            filteredAltlist22+=["NA"]
            filteredAAlist212+=[testdicttransAA1[y[0]][int(y[1])-1]]
            filteredAAlist222+=[dicttransAA1[y[0]][int(y[1])-1]]
            filteredAAPoslist122+=[y[1]]
            filteredcodondepthlist22+=["NA"]
            filteredAFlist12+=["NA"]
            filteredmutationlist12+=["Wildtype"]
            filteredQDlist12+=["NA"]
            filteredSORlist12+=["NA"]
            filteredMQlist12+=["NA"]
            filteredMQRankSumlist12+=["NA"]
            filteredfilterscorelist12+=["NA"]
            filteredDescriptionlist12+=["NA,NA,NA,NA"]



candidateslist1=[]
reportablelist1=[]
#print(df1Gene)
#print(df1AAPos)
#print(df1RefAA)
#print(df1AltAA)
wholeAAlist1=[]
wholeCandlist1=[]
wholeReportlist1=[]
#print(len(Genelist2))
#print(len(Genelist1))
#print(len(AAPoslist12))
#print(len(AAlist21))
#print(len(AAlist22))
#print(Genelist2)
for x in range(len(filteredGenelist22)):
    #print(len(Genelist2))
    #print(len(Genelist1))
    #print(len(AAPoslist12))
    #print(len(AAlist21))
    #print(len(AAlist22))
    wholeAAlist1.append([filteredGenelist22[x],filteredAAPoslist122[x],filteredAAlist212[x],filteredAAlist222[x]])
for k in range(len(df1Gene)):
    if df1Gene[k] == "Pfcrt":
        wholeCandlist1.append(["PfCRT",(df1AAPos[k]),df1RefAA[k],df1AltAA[k]])
    elif df1Gene[k] == "Pfdhfr":
        wholeCandlist1.append(["DHFR",(df1AAPos[k]),df1RefAA[k],df1AltAA[k]])
    elif df1Gene[k] == "PfK13":
        wholeCandlist1.append(["K13",(df1AAPos[k]),df1RefAA[k],df1AltAA[k]])
    elif df1Gene[k] == "Pfmdr1":
        wholeCandlist1.append(["PfMDR1",(df1AAPos[k]),df1RefAA[k],df1AltAA[k]])
    elif df1Gene[k] == "Pfdhps":
        wholeCandlist1.append(["DHPS",(df1AAPos[k]),df1RefAA[k],df1AltAA[k]])
    elif df1Gene[k] == "MT":
        wholeCandlist1.append(["mitochondrial_genome",(df1AAPos[k]),df1RefAA[k],df1AltAA[k]])
for k in range(len(df2Gene)):
    if df2Gene[k] == "Pfcrt":
        wholeReportlist1.append(["PfCRT",(df2AAPos[k]),df2RefAA[k],df2AltAA[k]])
    elif df2Gene[k] == "Pfdhfr":
        wholeReportlist1.append(["DHFR",(df2AAPos[k]),df2RefAA[k],df2AltAA[k]])
    elif df2Gene[k] == "PfK13":
        wholeReportlist1.append(["K13",(df2AAPos[k]),df2RefAA[k],df2AltAA[k]])
    elif df2Gene[k] == "Pfmdr1":
        wholeReportlist1.append(["PfMDR1",(df2AAPos[k]),df2RefAA[k],df2AltAA[k]])
    elif df2Gene[k] == "Pfdhps":
        wholeReportlist1.append(["DHPS",(df2AAPos[k]),df2RefAA[k],df2AltAA[k]])
    elif df2Gene[k] == "mitochondrial_genome":
        wholeReportlist1.append([df2Gene[k],(df2AAPos[k]),df2RefAA[k],df2AltAA[k]])
#print(wholeAAlist1)
#print(wholeReportlist1)
#print(['PfCRT', 74, 'M', 'I'])
#candidateslist1=[]
#print(wholeAAlist1)
for x in ((wholeAAlist1)):
    #print(x)
    if x in wholeCandlist1:
        candidateslist1+=["candidates"]
    else: candidateslist1+=["no"]
#    print(x)
    if x in wholeReportlist1:
        reportablelist1+=["reportable"]
    else:reportablelist1+=["no"]
    #print("test")
    #print(Genelist2)
    #print(Genelist2[x])
    #print(AAPoslist12[x])
    #print(AAlist21[x])
    #print(AAlist22[x])
filteredcandidateslist12=candidateslist1
filteredreportablelist12=reportablelist1
#print(len(filteredGenelist2))
#print(len(filteredPoslist2))
#print(len(filtereddepthlist2up))
#print(len(filteredReflist2))
#print(len(filteredAltlist2))
#print(len(filteredAAlist21))
#print(len(filteredAAlist22))
#print(len(filteredAAPoslist12))
#print(len(filteredcodondepthlist2))
#print(len(filteredAFlist1))
#print(len(filteredmutationlist1))
#print("test")
#print(len(filteredQDlist1))
#print(len(filteredSORlist1))
#print(len(filteredMQlist1))
#print(len(filteredMQRankSumlist1))
#print(len(filteredfilterscorelist1))
#print(len(filteredDescriptionlist1))
#print(len(filteredcandidateslist1))
#print(len(filteredreportablelist1))


d3={"Gene":filteredGenelist22, "BasePOS":filteredPoslist22, "BaseDepth":filtereddepthlist2up2, "Ref":filteredReflist22, "Alt":filteredAltlist22, "AAref":filteredAAlist212, "AAalt":filteredAAlist222, "AAPOS":filteredAAPoslist122, "CodonCoverage":filteredcodondepthlist22, "AF": filteredAFlist12, "Mutation": filteredmutationlist12, "QD": filteredQDlist12, "SOR":filteredSORlist12, "MQ":filteredMQlist12, "MQRankSum":filteredMQRankSumlist12, "Filter":filteredfilterscorelist12, "FilterDescription":filteredDescriptionlist12, "Candidates":filteredcandidateslist12, "Reportable":filteredreportablelist12}
df3=pd.DataFrame(data=d3)
#df.to_csv('snpreport4.csv', index=False,) 
df3 = df3.sort_values(df3.columns[0], ascending = False)
df3.to_csv(args.name+'.csv', index=False,) 