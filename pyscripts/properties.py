class properties:

    def __init__(self, unfiltered_path):
        self.unfiltered_path = unfiltered_path
        return

    def propertiesprocess(self,Genelist1,POSlist1,DepthVal,DepthVal3,currentag1,wholefilnums,wholefilnumsx):
        Genelist2 = []
        POSlist2 = []
        currentag2 = []
        AFlist1 = []
        QDlist1 = []
        SORlist1 = []
        MQlist1 = []
        MQRankSumlist1 = []
        DP4list1 = []
        ADlist1 = []
        f2 = open(self.unfiltered_path, "r")
        Geneset1 = []
        for k in range(len(Genelist1)):
            Geneset1 += [[Genelist1[k], POSlist1[k]]]
        GenePosSet1 = []
        for line in f2:
            count = 0
            tempDP4list = []
            #print(line)
            for word in line.split():
                # print(line)
                # print(Genelist1[x])
                # print(word)
                # if count==0 and word == Genelist1[x]:
                # print(word+"test1111")
                if count == 0:
                    word1 = ""
                    if word.startswith("#") == False:
                        word1 = word
                        # print(word)
                        # print(word1)
                # print(word1)
                # print(count)
                # print(word1)
                # temppos3=""
                # tempgene1=""
                # print(line)
                if count == 1:
                    # print(word1)
                    # print(word1)
                    if word1 != "":
                        # filtercount=-1
                        # print(word1)
                        # print(word)
                        # print(POSlist1)
                        # print(Genelist1)
                        # print(word1)
                        for x in range(len(Genelist1)):
                            # print(str(int(POSlist1[x])-1))
                            # if word1== Genelist1[x]:print(word1,word)
                            # print(POSlist1)
                            # print(word)
                            # if word.startswith("c."):
                            # temppos2=""
                            # for x in word:
                            # print(word)
                            # if x.isdigit():
                            # print(x)
                            # temppos2+=x
                            # print(temppos2)
                            # print(word1)
                            # print(Genelist1[x])
                            # print(POSlist1[x])
                            # print(word)
                            # print(Genelist1[x],POSlist1[x])
                            # print(word1,word)
                            # print(POSlist1)
                            # print(type(POSlist1[x]))
                            # if word1 == Genelist1[x]:print("no")
                            # print(word)
                            # print(POSlist1[x])
                            # if word == POSlist1[x]: print("yes")
                            #print(word1,word)
                            #print(Genelist1[x],POSlist1[x])
                            #print(DepthVal)
                            if (
                                word1 == Genelist1[x]
                                and word == POSlist1[x]
                                and line.startswith("#") == False
                                and (word1, word) in DepthVal
                            ):
                                # print(word)
                                # print(word)
                                # print(word)
                                # filtercount+=1
                                # print(POSlist1[x])
                                # print("test"+POSlist1[x])
                                # print(Genelist1[x],POSlist1[x])
                                Genelist2 += [Genelist1[x]]
                                #print(POSlist1[x])
                                POSlist2 += [POSlist1[x]]
                                currentag2 += [currentag1[x]]
                                GenePosSet1 += [[Genelist1[x], POSlist1[x]]]
                                # temppos3=POSlist1[x]
                                if line.find(";AF=") != -1:
                                    # print(line[int(line.find("AF="))+3:len(line)].find(";"))
                                    # print(line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(";")])
                                    # print(line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(",")].isdigit())
                                    # tempAF=""
                                    countAF = 0
                                    for k in line[
                                        int(line.find("AF="))
                                        + 3 : int(line.find("AF="))
                                        + 3
                                        + line[int(line.find("AF=")) + 3 : len(line)].find(",")
                                    ]:
                                        if k.isdigit() == False and k != ".":
                                            AFlist1 += [
                                                line[
                                                    int(line.find("AF="))
                                                    + 3 : int(line.find("AF="))
                                                    + 3
                                                    + line[
                                                        int(line.find("AF=")) + 3 : len(line)
                                                    ].find(";")
                                                ]
                                            ]
                                            break
                                        # print(len(line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(",")]))
                                        # print("COUNTAF"+str(countAF))
                                        if (
                                            countAF
                                            == len(
                                                line[
                                                    int(line.find("AF="))
                                                    + 3 : int(line.find("AF="))
                                                    + 3
                                                    + line[
                                                        int(line.find("AF=")) + 3 : len(line)
                                                    ].find(",")
                                                ]
                                            )
                                            - 1
                                        ):
                                            AFlist1 += [
                                                line[
                                                    int(line.find("AF="))
                                                    + 3 : int(line.find("AF="))
                                                    + 3
                                                    + line[
                                                        int(line.find("AF=")) + 3 : len(line)
                                                    ].find(",")
                                                ]
                                            ]
                                            break
                                        countAF += 1
                                    # else:AFlist1+=[line[int(line.find("AF="))+3:int(line.find("AF="))+3+line[int(line.find("AF="))+3:len(line)].find(";")]]
                                    # for k in range(int(line.find("AF="))+3,line[int(line.find("AF="))+3:len(line)].find(";")):
                                    #    if line[k]==";":
                                    #        break
                                    #    if line[line.find("AF=")+3:k].find(",") == -1:
                                    #        AFlist1+=[float(line[line.find("AF=")+3:k])]
                                    #    elif line[line.find("AF=")+3:k].find(",") != -1:
                                    #        AFlist1+=[float(line[line.find("AF=")+3:line.find("AF=")+5])]
                                else:
                                    AFlist1 += ["None"]
                                # print(AFlist1)
                                if line.find("MQ=") != -1:
                                    # if line.find("MQ=")+3 != 2:
                                    #    for k in range(int(line.find("MQ=")+3),len(line)):
                                    #        if line[k]==";":
                                    #            break
                                    #        MQlist1+=[float(line[line.find("MQ=")+3:k])]
                                    MQlist1 += [
                                        line[
                                            int(line.find("MQ="))
                                            + 3 : int(line.find("MQ="))
                                            + 3
                                            + line[int(line.find("MQ=")) + 3 : len(line)].find(
                                                ";"
                                            )
                                        ]
                                    ]
                                else:
                                    MQlist1 += ["None"]
                                if line.find("MQRankSum=") != -1:
                                    #    for k in range(int(line.find("MQRankSum="))+10,len(line)):
                                    #        if line[k]==";":
                                    #           break
                                    #        MQRankSumlist1+=[float(line[line.find("MQRankSum=")+10:k])]
                                    # print(int(line.find("MQRankSum="))+10:int(line.find("MQRankSum="))+10+line[int(line.find("MQRankSum="))+10])
                                    # print(int(line.find("MQRankSum="))+10+line[int(line.find("MQRankSum="))+10:len(line)].find(";"))
                                    # print(line[int(line.find("MQRankSum="))+10:int(line.find("MQRankSum="))+10+line[int(line.find("MQRankSum="))+10:len(line)].find(";")])
                                    MQRankSumlist1 += [
                                        line[
                                            int(line.find("MQRankSum="))
                                            + 10 : int(line.find("MQRankSum="))
                                            + 10
                                            + line[
                                                int(line.find("MQRankSum=")) + 10 : len(line)
                                            ].find(";")
                                        ]
                                    ]
                                else:
                                    MQRankSumlist1 += ["None"]
                                if line.find("SOR=") != -1:
                                    # for k in range(int(line.find("SOR="))+4,len(line)):
                                    #    if line[k]==";":
                                    #        break
                                    #    SORlist1+=[float(line[line.find("SOR=")+4:k])]
                                    SORlist1 += [
                                        line[
                                            int(line.find("SOR="))
                                            + 4 : int(line.find("SOR="))
                                            + 4
                                            + line[int(line.find("SOR=")) + 4 : len(line)].find(
                                                ";"
                                            )
                                        ]
                                    ]
                                else:
                                    SORlist1 += ["None"]
                                if line.find("QD=") != -1:
                                    # for k in range(int(line.find("QD="))+3,len(line)):
                                    #    if line[k]==";":
                                    #        break
                                    #    QDlist1+=[float(line[line.find("QD=")+3:k])]
                                    QDlist1 += [
                                        line[
                                            int(line.find("QD="))
                                            + 3 : int(line.find("QD="))
                                            + 3
                                            + line[int(line.find("QD=")) + 3 : len(line)].find(
                                                ";"
                                            )
                                        ]
                                    ]
                                else:
                                    QDlist1 += ["None"]
                                # print(filternumber)
                                # print((word1,word) in DepthVal)
                                # print(DepthVal[word1,word])
                                DP4list1 += DepthVal[word1, word]
                                if (word1, word) in DepthVal3:
                                    ADlist1 += DepthVal3[word1, word]
                                    # print(DepthVal3[word1,word])
                                else:
                                    ADlist1 += ["NA"]
                                # print(ADlist1)
                                break
                # print(len(QDlist1))
                # print(GenePosSet1)
                # print(AFlist1)
                # print(line)
                if count == 7:
                    if word1 != "":
                        # print(word)
                        if line.startswith("#") == False:
                            if word.startswith("AC=") or word.startswith("AB="):
                                # print(word[word.find("c."):word.find("c.")+10])
                                # print(line)
                                # print(wholefilnums)
                                for x in range(len(wholefilnums)):
                                    # print(x)
                                    # print(wholefilnumsx[count])
                                    # print(POSlist2)
                                    # print(wholefilnumsx[count2])
                                    # print([word1,str(int(wholefilnumsx[x]))])
                                    # print(word)
                                    # find part of c.138 in c.1381 in different gene.....
                                    # print(line.find(wholefilnums[x]))
                                    # if word.find(wholefilnums[x]) != -1:
                                    # print(wholefilnums[x])
                                    # print(word.find(wholefilnums[x]))
                                    # print(len(wholefilnums[x]))
                                    # print(wholefilnums[x])
                                    # print(word[word.find(wholefilnums[x])+len(wholefilnums[x])])
                                    # print(word[word.find(wholefilnums[x])+len(wholefilnums[x])-len(wholefilnumsx[x])-1])
                                    if (
                                        word.find(wholefilnums[x]) != -1
                                        and str(
                                            word.find(wholefilnums[x])
                                            + len(wholefilnums[x])
                                            + 1
                                        ).isdigit()
                                        == False
                                        and [word1, wholefilnumsx[x]] in Geneset1
                                        and [word1, str(int(wholefilnumsx[x]))]
                                        not in GenePosSet1
                                    ):
                                        # print(word1)
                                        # print(wholefilnumsx[x])
                                        # print(wholefilnums[x])
                                        # print(wholefilnumsx[x])
                                        # print(word1)
                                        # print(word1)
                                        # print(wholefilnumsx[x])
                                        # print(word1)
                                        Genelist2 += [word1]
                                        # print(wholefilnumsx[x])
                                        POSlist2 += [wholefilnumsx[x]]
                                        # "QD": QDlist1, "SOR":SORlist1, "MQ":MQlist1, "MQRankSum":MQRankSumlist1,
                                        AFlist1 += [AFlist1[-1]]
                                        QDlist1 += [QDlist1[-1]]
                                        SORlist1 += [SORlist1[-1]]
                                        MQlist1 += [MQlist1[-1]]
                                        MQRankSumlist1 += [MQRankSumlist1[-1]]
                                        GenePosSet1 += [[word1, str(int(wholefilnumsx[x]))]]
                                        DP4list1 += [DP4list1[-1]]
                                        ADlist1 += [ADlist1[-1]]
                                    # print(Genelist1[count])
                                    # print(wholefilnumsx[count])
                count += 1
        return Genelist2,POSlist2,AFlist1,QDlist1,SORlist1,MQlist1,MQRankSumlist1,GenePosSet1,DP4list1,ADlist1, currentag2