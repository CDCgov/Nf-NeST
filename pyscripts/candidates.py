class candidates:

    def __init__(self,candidates_path,voi_path):
        self.candidates_path=candidates_path
        self.voi_path=voi_path
        return

    def candidatesprocess(self,Genelist1,Genelist2,POSlist1,POSlist2,Reflist1,Altlist1,depthlist1up,codondepthlist1):
        import pandas as pd
        Reflist2 = []
        Altlist2 = []
        AAlist21 = []
        AAlist22 = []
        depthlist2up = []
        AAPoslist12 = []
        codondepthlist2 = []
        # print(Genelist2)
        # print(POSlist2)
        # print(Genelist1)
        # print(POSlist1)
        # print(len(Genelist2))
        # print(len(POSlist2))
        # print(len(AAPoslist12))
        # print(len(AAlist21))
        # print(len(AAlist22))
        # print(len(depthlist1up))
        # print(len(Genelist1))
        # print(len(AAlist2))
        # print(len(Reflist1))
        # print(len(Altlist1))
        # print(len(Reflist1))
        # print(len(AAlist1))
        # print(len(AAlist2))
        # print(Altlist1)

        count = 0
        for x in range(len(Genelist2)):
            for y in range(len(Genelist1)):
                if Genelist2[x] == Genelist1[y] and POSlist2[x] == POSlist1[y]:
                    count += 1
                    Reflist2 += [Reflist1[y]]
                    Altlist2 += [Altlist1[y]]
                    depthlist2up += [depthlist1up[y]]
                    codondepthlist2 += [codondepthlist1[y]]
                    break
            # Geneset2=[]
            # for k in range(len(Genelist1)):
            #    Geneset2+=[[Genelist1[k],POSlist1[k]]]
            # for x in range(len(Genelist2)):
            #    for y in range(len(Genelist1)):
            #        if [Genelist2[x],POSlist2[x]] not in Geneset2:
            #            print(Genelist2[x],POSlist2[x])
            # for y in POSlist2:
            #    if y not in POSlist1:
            #        print(y)
            # print(len(Genelist2))
            # print(count)
            # print(len(wholefilnumsx))
            # print(len(wholefilnums))
            #        print("order", x, wholelenlist[x], len(Genelist2))
        xls = pd.ExcelFile(self.candidates_path)
        df1 = pd.read_excel(xls)
        df1Gene = df1["Gene"].tolist()
        df1RefAA = df1["RefAA"].tolist()
        df1AltAA = df1["AltAA"].tolist()
        df1AAPos = df1["AAPos"].tolist()
        # df1.to_csv("test1.csv", index=False)
        df2 = pd.read_csv(self.voi_path)
        df2Gene = df2["Gene"].tolist()
        df2RefAA = df2["RefAA"].tolist()
        df2AltAA = df2["AltAA"].tolist()
        df2AAPos = df2["AAPos"].tolist()

        return df1Gene,df1RefAA,df1AltAA,df1AAPos,df2Gene,df2RefAA,df2AltAA,df2AAPos, Reflist2, Altlist2, depthlist2up, codondepthlist2