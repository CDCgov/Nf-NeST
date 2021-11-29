class candirepo:

    def __init__(self):
        return

    def candirepoprocess(self,filteredGenelist22,filteredAAPoslist122,filteredAAlist212,filteredAAlist222,df1AAPos,df1RefAA,df1AltAA,
        df2AAPos, df2RefAA, df2AltAA, df2Gene, df1Gene):
        candidateslist1 = []
        reportablelist1 = []
        # print(df1Gene)
        # print(df1AAPos)
        # print(df1RefAA)
        # print(df1AltAA)
        wholeAAlist1 = []
        wholeCandlist1 = []
        wholeReportlist1 = []
        # print(len(Genelist2))
        # print(len(Genelist1))
        # print(len(AAPoslist12))
        # print(len(AAlist21))
        # print(len(AAlist22))
        # print(Genelist2)
        for x in range(len(filteredGenelist22)):
            # print(len(Genelist2))
            # print(len(Genelist1))
            # print(len(AAPoslist12))
            # print(len(AAlist21))
            # print(len(AAlist22))
            wholeAAlist1.append(
                [
                    filteredGenelist22[x],
                    filteredAAPoslist122[x],
                    filteredAAlist212[x],
                    filteredAAlist222[x],
                ]
            )
        for k in range(len(df1Gene)):
            if df1Gene[k] == "Pfcrt":
                wholeCandlist1.append(["PfCRT", (df1AAPos[k]), df1RefAA[k], df1AltAA[k]])
            elif df1Gene[k] == "Pfdhfr":
                wholeCandlist1.append(["DHFR", (df1AAPos[k]), df1RefAA[k], df1AltAA[k]])
            elif df1Gene[k] == "PfK13":
                wholeCandlist1.append(["K13", (df1AAPos[k]), df1RefAA[k], df1AltAA[k]])
            elif df1Gene[k] == "Pfmdr1":
                wholeCandlist1.append(["PfMDR1", (df1AAPos[k]), df1RefAA[k], df1AltAA[k]])
            elif df1Gene[k] == "Pfdhps":
                wholeCandlist1.append(["DHPS", (df1AAPos[k]), df1RefAA[k], df1AltAA[k]])
            elif df1Gene[k] == "MT":
                wholeCandlist1.append(
                    ["mitochondrial_genome", (df1AAPos[k]), df1RefAA[k], df1AltAA[k]]
                )
        for k in range(len(df2Gene)):
            if df2Gene[k] == "Pfcrt":
                wholeReportlist1.append(["PfCRT", (df2AAPos[k]), df2RefAA[k], df2AltAA[k]])
            elif df2Gene[k] == "Pfdhfr":
                wholeReportlist1.append(["DHFR", (df2AAPos[k]), df2RefAA[k], df2AltAA[k]])
            elif df2Gene[k] == "PfK13":
                wholeReportlist1.append(["K13", (df2AAPos[k]), df2RefAA[k], df2AltAA[k]])
            elif df2Gene[k] == "Pfmdr1":
                wholeReportlist1.append(["PfMDR1", (df2AAPos[k]), df2RefAA[k], df2AltAA[k]])
            elif df2Gene[k] == "Pfdhps":
                wholeReportlist1.append(["DHPS", (df2AAPos[k]), df2RefAA[k], df2AltAA[k]])
            elif df2Gene[k] == "mitochondrial_genome":
                wholeReportlist1.append([df2Gene[k], (df2AAPos[k]), df2RefAA[k], df2AltAA[k]])
        # print(wholeAAlist1)
        # print(wholeReportlist1)
        # print(['PfCRT', 74, 'M', 'I'])
        # candidateslist1=[]
        # print(wholeAAlist1)
        for x in wholeAAlist1:
            # print(x)
            if x in wholeCandlist1:
                candidateslist1 += ["candidates"]
            else:
                candidateslist1 += ["no"]
            #    print(x)
            if x in wholeReportlist1:
                reportablelist1 += ["reportable"]
            else:
                reportablelist1 += ["no"]
            # print("test")
            # print(Genelist2)
            # print(Genelist2[x])
            # print(AAPoslist12[x])
            # print(AAlist21[x])
            # print(AAlist22[x])
        filteredcandidateslist12 = candidateslist1
        filteredreportablelist12 = reportablelist1
        return filteredcandidateslist12, filteredreportablelist12