class GPRA:
    #Obtain Gene, Position, Alt

    def __init__(self, filtered_path):
        self.filtered_path = filtered_path
        return

    def GPRAPROCESS(self):
        #print(self.filtered_path)

        AAdic = {
            "Ala": "A",
            "Arg": "R",
            "Asn": "N",
            "Asp": "D",
            "Asx": "B",
            "Cys": "C",
            "Glu": "E",
            "Gln": "Q",
            "Glx": "Z",
            "Gly": "G",
            "His": "H",
            "Ile": "I",
            "Leu": "L",
            "Lys": "K",
            "Met": "M",
            "Phe": "F",
            "Pro": "P",
            "Ser": "S",
            "Thr": "T",
            "Trp": "W",
            "Tyr": "Y",
            "Val": "V",
            "STOP": "STOP",
        }
        
        with open(self.filtered_path, "r") as f1:
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
            currentagreement1=""
            currentag1=[]
            for line in f1:
                #print(line)
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
                        if count==6:
                            currentag1+=[word]
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
                                currentag1=currentag1[:-1]
                            #print(word[5+len(temppos):5+len(temppos)+3])
                            #if Altlist1!=[]:
                            #    print(currentagreement1)
                            #    print(Genelist1[-1])
                            #    print(POSlist1[-1])
                            #    print(Reflist1[-1])
                            #    print(Altlist1[-1])
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
                        if count==6:
                            currentag1+=[word]
                        if word.startswith("c.") or word.startswith("n."):
                            Altlist1+=[word[-1]]
                        count+=1
        return Genelist1,POSlist1,Reflist1,Altlist1