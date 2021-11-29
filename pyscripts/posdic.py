class posdicwi:

    def __init__(self):
        return

    def posdicwiprocess(self,bedname1, wholepotlist1, dictrange, depthlist1):
        posdic1 = {}

        # print(bedname1)
        # print(bedname1)
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
                        # print(bedname1[x])
                        # print(AAcount)
                        # print(count)
                        if [bedname1[x], str(AAcount)] in wholepotlist1 and bedname1[
                            x
                        ] == "mitochondrial_genome":
                            # if
                            # print(posdic1[bedname1[x],str(AAcount)])
                            posdic1[bedname1[x], str(AAcount)] = z - 3493
                        if [bedname1[x], str(AAcount)] in wholepotlist1 and bedname1[
                            x
                        ] != "mitochondrial_genome":
                            posdic1[bedname1[x], str(AAcount)] = z - 2
                        count += 1

                if len(dictrange[bedname1[x]][0]) > 1:
                    for z in range(
                        dictrange[bedname1[x]][0][y] - 1, dictrange[bedname1[x]][1][y]
                    ):
                        # print(z)
                        if count % 3 == 1:
                            # print(count)
                            AAcount += 1
                        if [bedname1[x], str(AAcount)] in wholepotlist1:
                            # print(bedname1[x])
                            if bedname1[x] == "PfCRT":
                                posdic1[bedname1[x], str(AAcount)] = z - 2
                            if bedname1[x] == "DHPS":
                                # if
                                # print(posdic1[bedname1[x],str(AAcount)])
                                posdic1[bedname1[x], str(AAcount)] = z - 2
                                # if (bedname1[x],str(AAcount)) in posdic1:
                                # print(AAcount)
                                # print(bedname1[x])
                                # print(posdic1[bedname1[x],str(AAcount)])
                        count += 1
        # print(posdic1)
        # print(posdic1["DHPS",str(540)])
        # print(posdic1['PfCRT', '72'])
        ######################################################################
        # for x in wholepotlist1:
        wildic1 = {}
        for x in wholepotlist1:
            # print(x)
            for y in depthlist1:
                # print(y)
                # print(x[0])
                # print(posdic1[x[0],x[1]])
                # print(y[1])
                # print(y[1])
                # if x[0]=="mitochondrial_genome":
                #    print(x[1])
                #        print("True")
                # print(y[1])
                if x[0] == y[0] and str(posdic1[x[0], x[1]]) == y[1]:
                    # print(posdic1[x[0],x[1]])
                    # if x[0]=="mitochondrial_genome":
                    #    print("True")
                    wildic1[x[0], x[1]] = y[2]
                if x[0] == y[0] == "mitochondrial_genome" and str(posdic1[x[0], x[1]]) == str(
                    int(y[1]) - 3493
                ):
                    wildic1[x[0], x[1]] = y[2]
        return posdic, wildic1