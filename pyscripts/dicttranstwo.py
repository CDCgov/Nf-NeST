class dicttranstwo:

    def __init__(self):
        return

    def dicttranstwoprocess(self, bedname1, dicttrans2, testdic, dictrange, GenePosSet2):
        dicttrans22 = {}
        testdicttrans22 = {}
        #############################################
        ###############cds region applied############
        #############################################
        # print(dicttrans2)
        # print(1%3)
        # print(dictrange)
        # for x in range(3):
        #    print(x)

        # print(dictrange)
        # print(1 in dictrange["K13"][0])
        BASEAAPOSdic1 = {}
        for x in range(len(bedname1)):
            tempfa2 = ""
            temptestfa2 = ""
            tempAAPOS = ""
            count = 0
            AAcount = 0
            for y in range(len(dictrange[bedname1[x]][0])):
                # print(dictrange[bedname1[x][0]])
                if len(dictrange[bedname1[x]][0]) == 1:
                    for z in range(
                        dictrange[bedname1[x]][0][y] - 1, dictrange[bedname1[x]][1][y] - 1
                    ):
                        # print(z)
                        if count % 3 == 1:
                            AAcount += 1
                        count += 1
                        if ([bedname1[x], str(z)]) in GenePosSet2:
                            BASEAAPOSdic1[bedname1[x], str(z)] = AAcount
                        tempfa2 += dicttrans2[bedname1[x]][z]
                        temptestfa2 += testdic[bedname1[x]][z]
                if len(dictrange[bedname1[x]][0]) > 1:
                    for z in range(
                        dictrange[bedname1[x]][0][y] - 1, dictrange[bedname1[x]][1][y]
                    ):
                        # print(z)
                        if count % 3 == 1:
                            # print(count)
                            AAcount += 1
                        count += 1
                        if ([bedname1[x], str(z)]) in GenePosSet2:
                            # print(z)
                            # print(AAcount)
                            BASEAAPOSdic1[bedname1[x], str(z)] = AAcount
                        tempfa2 += dicttrans2[bedname1[x]][z]
                        temptestfa2 += testdic[bedname1[x]][z]
            # print(len(tempfa2))
            dicttrans22[bedname1[x]] = tempfa2
            testdicttrans22[bedname1[x]] = temptestfa2
        return dicttrans22,testdicttrans22, BASEAAPOSdic1