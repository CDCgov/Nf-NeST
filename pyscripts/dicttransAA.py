class dicttransAA:

    def __init__(self):
        return

    def dicttransAAprocess(self, bedname1, dicttrans22, testdicttrans22, AAdic, codondic):
        dicttransAA1 = {}
        testdicttransAA1 = {}
        for x in range(len(bedname1)):
            tempAA1 = ""
            testtempAA1 = ""
            count = 0
            # print(bedname1[x])
            for z in range(len(dicttrans22[bedname1[x]]) - 2):
                count += 1
                if count == 1:
                    # print(line.strip()[z:z+3])
                    # print(bedname1[x])
                    # print(dicttrans22[bedname1[x]][z:z+3])
                    tempAA1 += AAdic[codondic[dicttrans22[bedname1[x]][z : z + 3]]]
                if count == 3:
                    count = 0
            dicttransAA1[bedname1[x]] = tempAA1
            count = 0
            for k in range(len(testdicttrans22[bedname1[x]]) - 2):
                count += 1
                if count == 1:
                    # if bedname1[x]=="PfMDR1":
                    # print(testdicttrans22[bedname1[x]][k:k+3])
                    testtempAA1 += AAdic[codondic[testdicttrans22[bedname1[x]][k : k + 3]]]
                if count == 3:
                    count = 0
            testdicttransAA1[bedname1[x]] = testtempAA1

        return dicttransAA1,testdicttransAA1