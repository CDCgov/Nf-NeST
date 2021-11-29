class filterdes:

    def __init__(self):
        return

    def filterdesprocess(self,QDlist1,SORlist1,MQlist1,MQRankSumlist1):
        Descriptionlist1 = []
        filterscorelist1 = []

        for x in range(len(QDlist1)):
            filternumber = 0
            Description = ""
            # print(QDlist1[filtercount])
            if QDlist1[x] != "None":
                filternumber += float(QDlist1[x]) / 80
                if float(SORlist1[x]) > 30:
                    Description += "Good,"
                else:
                    Description += "Bad,"
            else:
                Description += "NA,"
            if SORlist1[x] != "None":
                if float(SORlist1[x]) < 3:
                    Description += "Good,"
                else:
                    Description += "Bad,"
            else:
                Description += "NA,"
            if MQlist1[x] != "None":
                if float(MQlist1[x]) < 70:
                    filternumber += float(MQlist1[x]) / 140
                else:
                    filternumber += 1 / 2
                if 40 < float(MQlist1[x]) < 60:
                    Description += "Good,"
                else:
                    Description += "Bad,"
            else:
                Description += "NA,"
            # print(MQRankSumlist1)
            if MQRankSumlist1[x] != "None":
                if -5 < abs(float(MQRankSumlist1[x])) < 5:
                    Description += "Good"
                else:
                    Description += "Bad"
            else:
                Description += "NA"
            Descriptionlist1 += [Description]
            # print(Description)
            filterscorelist1 += [filternumber]
        return Descriptionlist1,filterscorelist1