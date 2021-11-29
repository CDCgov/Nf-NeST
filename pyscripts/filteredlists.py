class filteredlists:
    def __init__(self):
        return

    def filteredlistsprocess(self, Genelist2, POSlist2, Reflist2, Altlist2, BASEAAPOSdic1, combineddflist1,combineddflist2,currentag2, depthlist2up, AFlist1, testdicttransAA1, dicttransAA1, DP4list1, ADlist1, QDlist1, SORlist1, MQlist1, MQRankSumlist1, filterscorelist1, Descriptionlist1, codondepthlist2):
        filteredGenelist2 = []
        filteredPoslist2 = []
        filtereddepthlist2up = []
        filteredReflist2 = []
        filteredAltlist2 = []
        filteredAAlist21 = []
        filteredAAlist22 = []
        filteredAAPoslist12 = []
        filteredcodondepthlist2 = []
        filteredAFlist1 = []
        filteredmutationlist1 = []
        filteredQDlist1 = []
        filteredSORlist1 = []
        filteredMQlist1 = []
        filteredMQRankSumlist1 = []
        filteredfilterscorelist1 = []
        filteredDescriptionlist1 = []
        filteredcandidateslist1 = []
        filteredreportablelist1 = []
        filteredDP4list1 = []
        filteredADlist1 = []
        filteredAGlist2 = []
        # print(testdicttransAA1)
        # print(dicttransAA1)

        # print(len(Genelist2))
        # print(Genelist2)
        # print(POSlist2)
        # print(BASEAAPOSdic1)
        combinedGenePOSlist1 = []
        for x in range(len(Genelist2)):
            if (Genelist2[x], (POSlist2[x])) in BASEAAPOSdic1:
                # print("True")
                # print([Genelist2[x],BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]])
                # print(BASEAAPOSdic1[Genelist2[x],(POSlist2[x])])
                # print(Genelist2[x],POSlist2[x])
                # print([Genelist2[x],str(BASEAAPOSdic1[Genelist2[x],(POSlist2[x])])])
                if [
                    Genelist2[x],
                    str(BASEAAPOSdic1[Genelist2[x], (POSlist2[x])]),
                ] in combineddflist1 or [
                    Genelist2[x],
                    str(BASEAAPOSdic1[Genelist2[x], (POSlist2[x])]),
                ] in combineddflist2:
                    # print("True")
                    combinedGenePOSlist1 += [
                        [Genelist2[x], BASEAAPOSdic1[Genelist2[x], (POSlist2[x])]]
                    ]

        # print(combinedGenePOSlist1)
        # print(combinedGenePOSlist1)
        for x in range(len(Genelist2)):
            # print(Genelist2[x],(POSlist2[x]))
            if (Genelist2[x], (POSlist2[x])) in BASEAAPOSdic1:
                # print([Genelist2[x],BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]])
                if [
                    Genelist2[x],
                    str(BASEAAPOSdic1[Genelist2[x], (POSlist2[x])]),
                ] in combineddflist1 or [
                    Genelist2[x],
                    str(BASEAAPOSdic1[Genelist2[x], (POSlist2[x])]),
                ] in combineddflist2:
                    # print(AAPoslist12[x])
                    # print("True")
                    filteredGenelist2 += [Genelist2[x]]
                    filteredPoslist2 += [POSlist2[x]]
                    filteredAGlist2 += [currentag2[x]]
                    filtereddepthlist2up += [depthlist2up[x]]
                    filteredReflist2 += [Reflist2[x]]
                    filteredAltlist2 += [Altlist2[x]]
                    # print(BASEAAPOSdic1[Genelist2[x],(POSlist2[x])])
                    # print(testdicttransAA1[Genelist2[x])
                    # print(dicttransAA1[Genelist2[x]])
                    # print([testdicttransAA1[Genelist2[x]][BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]-1]])
                    # print([dicttransAA1[Genelist2[x]][BASEAAPOSdic1[Genelist2[x],(POSlist2[x])]-1]])
                    filteredAAlist21 += [
                        testdicttransAA1[Genelist2[x]][
                            BASEAAPOSdic1[Genelist2[x], (POSlist2[x])] - 1
                        ]
                    ]
                    filteredAAlist22 += [
                        dicttransAA1[Genelist2[x]][
                            BASEAAPOSdic1[Genelist2[x], (POSlist2[x])] - 1
                        ]
                    ]
                    # print(AAuntranslatewilddic1[Genelist2[x],AAPoslist12[x]])
                    # print(AAtranslatewilddic1[Genelist2[x],AAPoslist12[x]])
                    filteredAAPoslist12 += [BASEAAPOSdic1[Genelist2[x], (POSlist2[x])]]
                    filteredcodondepthlist2 += [codondepthlist2[x]]
                    filteredAFlist1 += [AFlist1[x]]
                    if (
                        testdicttransAA1[Genelist2[x]][
                            BASEAAPOSdic1[Genelist2[x], (POSlist2[x])] - 1
                        ]
                        == dicttransAA1[Genelist2[x]][
                            BASEAAPOSdic1[Genelist2[x], (POSlist2[x])] - 1
                        ]
                    ):
                        # print("False")
                        filteredmutationlist1 += ["Wildtype"]
                    if DP4list1[x] != "None":
                        if float(DP4list1[x]) >= 0.5:
                            if (
                                testdicttransAA1[Genelist2[x]][
                                    BASEAAPOSdic1[Genelist2[x], (POSlist2[x])] - 1
                                ]
                                != dicttransAA1[Genelist2[x]][
                                    BASEAAPOSdic1[Genelist2[x], (POSlist2[x])] - 1
                                ]
                            ):
                                filteredmutationlist1 += ["Major"]
                        if float(DP4list1[x]) < 0.5:
                            if (
                                testdicttransAA1[Genelist2[x]][
                                    BASEAAPOSdic1[Genelist2[x], (POSlist2[x])] - 1
                                ]
                                != dicttransAA1[Genelist2[x]][
                                    BASEAAPOSdic1[Genelist2[x], (POSlist2[x])] - 1
                                ]
                            ):
                                filteredmutationlist1 += ["Minor"]
                    if DP4list1[x] == "None":
                        # print('TRUE')
                        if (
                            testdicttransAA1[Genelist2[x]][
                                BASEAAPOSdic1[Genelist2[x], (POSlist2[x])] - 1
                            ]
                            != dicttransAA1[Genelist2[x]][
                                BASEAAPOSdic1[Genelist2[x], (POSlist2[x])] - 1
                            ]
                        ):
                            filteredmutationlist1 += ["Wildtype"]
                    filteredDP4list1 += [DP4list1[x]]
                    filteredADlist1 += [ADlist1[x]]
                    filteredQDlist1 += [QDlist1[x]]
                    filteredSORlist1 += [SORlist1[x]]
                    filteredMQlist1 += [MQlist1[x]]
                    filteredMQRankSumlist1 += [MQRankSumlist1[x]]
                    filteredfilterscorelist1 += [filterscorelist1[x]]
                    filteredDescriptionlist1 += [Descriptionlist1[x]]
        return (filteredGenelist2, filteredPoslist2, filtereddepthlist2up, filteredReflist2, filteredAltlist2, filteredAltlist2, filteredAAlist21,
            filteredAAlist22, filteredAAPoslist12, filteredcodondepthlist2, filteredAFlist1, filteredmutationlist1, filteredQDlist1, filteredSORlist1,
            filteredMQlist1, filteredMQRankSumlist1, filteredfilterscorelist1, filteredDescriptionlist1, filteredcandidateslist1, filteredreportablelist1,
            filteredDP4list1, combinedGenePOSlist1, filteredAGlist2, filteredADlist1)
