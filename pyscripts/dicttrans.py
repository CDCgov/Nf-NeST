class dicttrans:

    def __init__(self, fasta_path):
        self.fasta_path=fasta_path
        return

    def dicttransprocess(self, bedname1, Genelist2, POSlist2, Altlist2):
        dicttrans1 = {}
        dicttrans2 = {}

        with open(self.fasta_path, "r") as fa1:
            for line in fa1:
                if line.find(">") != -1:
                    for x in range(len(bedname1)):
                        # print(line)
                        if line[1::].strip() == bedname1[x]:
                            tempname = bedname1[x]
                if line.find(">") == -1:
                    dicttrans1[tempname] = line.strip()

        # print(dicttrans1["PfMDR1"][3099])
        # print(dicttrans1["PfMDR1"][3100])
        # print(dicttrans1["PfMDR1"][3101])

        # print(Genelist2)
        # print(POSlist2)
        # print(Altlist2)
        GenePosSet2 = []
        GenePosdict = {}
        for x in range(len(Genelist2)):
            GenePosSet2 += [[Genelist2[x], POSlist2[x]]]
            GenePosdict[Genelist2[x], POSlist2[x]] = Altlist2[x]
        # print(Reflist2)
        # print(GenePosSet2)

        #############################################
        ###############no mutated base applied#######
        #############################################
        testdic = {}
        for x in range(len(bedname1)):
            tempfa1 = ""
            for y in dicttrans1[bedname1[x]]:
                tempfa1 += y
            testdic[bedname1[x]] = tempfa1
        # print(dictrange)

        # print(dicttrans1)

        # print(dicttrans1["PfCRT"])
        # print(GenePosdict)
        #############################################
        ###############mutated base applied##########
        #############################################
        # print(GenePosSet2)
        for x in range(len(bedname1)):
            tempfa1 = ""
            count = 0
            for y in dicttrans1[bedname1[x]]:
                count += 1
                # print([bedname1[x],str(count)] in GenePosSet2)
                if ([bedname1[x], str(count)]) in GenePosSet2:
                    # print("True")
                    # if bedname1[x]=="PfCRT":
                    # print(count)
                    # print(GenePosdict[bedname1[x],str(count)])
                    tempfa1 += GenePosdict[bedname1[x], str(count)]
                else:
                    tempfa1 += y
            dicttrans2[bedname1[x]] = tempfa1
        return dicttrans2,testdic,GenePosSet2