class combinewholepot:

    def __init__(self):
        return

    def combinewholepotprocess(self,filteredADlist1,filteredAGlist2,filteredGenelist2, filteredPoslist2, filtereddepthlist2up, filteredReflist2, filteredAltlist2, filteredAAlist21,
    filteredAAlist22, filteredAAPoslist12, filteredcodondepthlist2, filteredAFlist1, filteredmutationlist1, filteredQDlist1, filteredSORlist1,
    filteredMQlist1, filteredMQRankSumlist1, filteredfilterscorelist1, filteredDescriptionlist1, filteredcandidateslist1, filteredreportablelist1,
    filteredDP4list1, tempcovwordlist1, combinedGenePOSlist1,testdicttransAA1, dicttransAA1, wholepotlist1, wildic1, posdic1):
        filteredGenelist22 = []
        filteredPoslist22 = []
        filtereddepthlist2up2 = []
        filteredReflist22 = []
        filteredAltlist22 = []
        filteredAAlist212 = []
        filteredAAlist222 = []
        filteredAAPoslist122 = []
        filteredcodondepthlist22 = []
        filteredAFlist12 = []
        filteredmutationlist12 = []
        filteredQDlist12 = []
        filteredSORlist12 = []
        filteredMQlist12 = []
        filteredMQRankSumlist12 = []
        filteredfilterscorelist12 = []
        filteredDescriptionlist12 = []
        filteredcandidateslist12 = []
        filteredreportablelist12 = []
        filteredDP4list12 = []
        filteredADlist12 = []
        filteredAGlist12 = []


        combinedwildgeneAA = []
        for y in wholepotlist1:
            # count=0
            for x in range(len(combinedGenePOSlist1)):
                # count+=1
                # print(x)
                # print(len(combinedGenePOSlist1))
                if [y[0], int(y[1])] == combinedGenePOSlist1[x] and y[0] in tempcovwordlist1:
                    # print(combinedGenePOSlist1[x])
                    # print("True")
                    filteredGenelist22 += [filteredGenelist2[x]]
                    if filteredGenelist2[x] == "mitochondrial_genome":
                        filteredPoslist22 += [str(int(filteredPoslist2[x]) - 3491)]
                    else:
                        filteredPoslist22 += [filteredPoslist2[x]]
                    filtereddepthlist2up2 += [filtereddepthlist2up[x]]

                    filteredReflist22 += [filteredReflist2[x]]
                    filteredAltlist22 += [filteredAltlist2[x]]
                    filteredAAlist212 += [filteredAAlist21[x]]
                    # print(filteredAAlist21[x])
                    # print(filteredAAlist22[x])
                    filteredAAlist222 += [filteredAAlist22[x]]
                    filteredAGlist12 += [filteredAGlist2[x]]
                    filteredAAPoslist122 += [filteredAAPoslist12[x]]
                    filteredcodondepthlist22 += [filteredcodondepthlist2[x]]
                    filteredAFlist12 += [filteredAFlist1[x]]
                    filteredmutationlist12 += [filteredmutationlist1[x]]
                    filteredQDlist12 += [filteredQDlist1[x]]
                    filteredSORlist12 += [filteredSORlist1[x]]
                    filteredMQlist12 += [filteredMQlist1[x]]
                    filteredMQRankSumlist12 += [filteredMQRankSumlist1[x]]
                    filteredfilterscorelist12 += [filteredfilterscorelist1[x]]
                    filteredDescriptionlist12 += [filteredDescriptionlist1[x]]
                    filteredDP4list12 += [filteredDP4list1[x]]
                    filteredADlist12 += [filteredADlist1[x]]
                # print(combinedGenePOSlist1)
                if (
                    x == len(combinedGenePOSlist1) - 1
                    and ([y[0], int(y[1])] not in combinedGenePOSlist1)
                    and ([y[0], int(y[1])] not in combinedwildgeneAA)
                    and (y[0], y[1]) in wildic1
                ):
                    combinedwildgeneAA += [[y[0], int(y[1])]]
                    filteredGenelist22 += [y[0]]
                    # print(y[0])
                    filteredPoslist22 += [str(posdic1[y[0], y[1]])]
                    if (y[0], y[1]) in wildic1:
                        filtereddepthlist2up2 += [wildic1[y[0], y[1]]]
                    else:
                        filtereddepthlist2up2 += ["NA"]
                    filteredReflist22 += ["NA"]
                    filteredAltlist22 += ["NA"]
                    # print([testdicttransAA1[y[0]][int(y[1])-1]])
                    # print([dicttransAA1[y[0]][int(y[1])-1]])
                    filteredAAlist212 += [testdicttransAA1[y[0]][int(y[1]) - 1]]
                    filteredAAlist222 += [dicttransAA1[y[0]][int(y[1]) - 1]]
                    #print(y[0],y[1])
                    #print([testdicttransAA1[y[0]][int(y[1]) - 1]])
                    #print([dicttransAA1[y[0]][int(y[1]) - 1]])
                    filteredAGlist12 += ["NA"]
                    filteredAAPoslist122 += [y[1]]
                    filteredcodondepthlist22 += ["NA"]
                    filteredAFlist12 += ["NA"]
                    filteredmutationlist12 += ["Wildtype"]
                    filteredQDlist12 += ["NA"]
                    filteredSORlist12 += ["NA"]
                    filteredMQlist12 += ["NA"]
                    filteredMQRankSumlist12 += ["NA"]
                    filteredfilterscorelist12 += ["NA"]
                    filteredDescriptionlist12 += ["NA,NA,NA,NA"]
                    filteredDP4list12 += ["NA"]
                    filteredADlist12 += ["NA"]

        # print(wildic1)
        # print(wholepotlist1)
        # print(combinedGenePOSlist1)
        if combinedGenePOSlist1 == []:
            for y in wholepotlist1:
                # count=0
                # print("True")
                # print(y[0])
                # if ([y[0],int(y[1])] not in combinedwildgeneAA) and (y[0],y[1]) in wildic1:# and int(wildic1[y[0],y[1]])>0:
                # print("True")
                # print(y[0])
                combinedwildgeneAA += [[y[0], int(y[1])]]
                filteredGenelist22 += [y[0]]
                filteredPoslist22 += [str(posdic1[y[0], y[1]])]
                if (y[0], y[1]) in wildic1:
                    filtereddepthlist2up2 += [wildic1[y[0], y[1]]]
                else:
                    filtereddepthlist2up2 += ["NA"]
                filteredReflist22 += ["NA"]
                filteredAltlist22 += ["NA"]
                filteredAGlist12 += ["NA"]
                filteredAAlist212 += [testdicttransAA1[y[0]][int(y[1]) - 1]]
                filteredAAlist222 += [dicttransAA1[y[0]][int(y[1]) - 1]]
                filteredAAPoslist122 += [y[1]]
                filteredcodondepthlist22 += ["NA"]
                filteredAFlist12 += ["NA"]
                filteredmutationlist12 += ["Wildtype"]
                filteredQDlist12 += ["NA"]
                filteredSORlist12 += ["NA"]
                filteredMQlist12 += ["NA"]
                filteredMQRankSumlist12 += ["NA"]
                filteredfilterscorelist12 += ["NA"]
                filteredDescriptionlist12 += ["NA,NA,NA,NA"]
                filteredDP4list12 += ["NA"]
                filteredADlist12 += ["NA"]
        return (filteredGenelist22,filteredPoslist22,filtereddepthlist2up2,filteredReflist22,filteredAltlist22,filteredAAlist212,
            filteredAAlist222,filteredAAPoslist122,filteredcodondepthlist22,filteredAFlist12,filteredmutationlist12,filteredQDlist12,
            filteredSORlist12,filteredMQlist12,filteredMQRankSumlist12,filteredfilterscorelist12,filteredDescriptionlist12,
            filteredcandidateslist12,filteredreportablelist12,filteredDP4list12,filteredADlist12,filteredAGlist12)