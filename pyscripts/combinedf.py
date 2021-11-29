class combinedf:

    def __init__(self):
        return

    def combinedfprocess(self, df1Gene, df2Gene, df1AAPos, df2AAPos):
        combineddflist1 = []
        combineddflist2 = []
        for x in range(len(df1Gene)):
            if df1Gene[x] == "Pfcrt":
                combineddflist1 += [["PfCRT", str(df1AAPos[x])]]
            elif df1Gene[x] == "Pfdhfr":
                combineddflist1 += [["DHFR", str(df1AAPos[x])]]
            elif df1Gene[x] == "PfK13":
                combineddflist1 += [["K13", str(df1AAPos[x])]]
            elif df1Gene[x] == "Pfmdr1":
                combineddflist1 += [["PfMDR1", str(df1AAPos[x])]]
            elif df1Gene[x] == "Pfdhps":
                combineddflist1 += [["DHPS", str(df1AAPos[x])]]
            elif df1Gene[x] == "MT":
                combineddflist1 += [["mitochondrial_genome", str(df1AAPos[x])]]
        for x in range(len(df2Gene)):
            if df2Gene[x] == "Pfcrt":
                combineddflist2 += [["PfCRT", str(df2AAPos[x])]]
            elif df2Gene[x] == "Pfdhfr":
                combineddflist2 += [["DHFR", str(df2AAPos[x])]]
            elif df2Gene[x] == "Pfk13":
                combineddflist2 += [["K13", str(df2AAPos[x])]]
            elif df2Gene[x] == "Pfmdr1":
                combineddflist2 += [["PfMDR1", str(df2AAPos[x])]]
            elif df2Gene[x] == "Pfdhps":
                combineddflist2 += [["DHPS", str(df2AAPos[x])]]
            elif df2Gene[x] == "MT":
                combineddflist2 += [["mitochondrial_genome", str(df2AAPos[x])]]
        return combineddflist1,combineddflist2