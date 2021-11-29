class bedrange:

    def __init__(self, bedfile_path):
        self.bedfile_path=bedfile_path
        return

    def bedrangeprocess(self):
        bedname1 = []
        bedstart1 = []
        bedend1 = []
        tempstart1 = []
        tempend1 = []
        with open(self.bedfile_path, "r") as fb1:
            # print(line)
            for line in fb1:
                count = 0
                # print(len(line.split()))
                if line.split()[0] == line.split()[3] or line.split()[3] == "CYTOb":
                    for word in line.split():
                        count += 1
                        if count == 1:
                            bedname1 += [word]
                        if count == len(line.split()):
                            tempstart1 = []
                            for x in word[:-1].split(","):
                                # print(x)
                                tempstart1 += [int(x)]
                            bedstart1 += [tempstart1]
                            # tempstart1=[word[:-1]]
                            # print(word[:-1])
                        if count == len(line.split()) - 1:
                            # print(tempstart1)
                            # print(word[:-1])
                            # print([x+y for x,y in zip(tempstart1,[word[:-1]])])
                            tempend1 = []
                            for x in word[:-1].split(","):
                                # print(x)
                                tempend1 += [int(x)]
                    bedend1 += [[x + y for x, y in zip(tempstart1, tempend1)]]
                # print(bedend1)
        # print(bedstart1)
        # print(bedend1)
        # print(bedname1)
        dictrange = {}
        # print(bedstart1)
        # print(bedname1)
        # print(bedend1)
        for x in range(len(bedname1)):
            # print(bedstart1[x])
            exec(bedname1[x] + "Startlist1" + "=" + str(bedstart1[x]))
            # print(bedname1[x])
            # print(PfCRTStartlist1)
            exec(bedname1[x] + "Endlist1" + "=" + str(bedend1[x]))
            # bedname1[x]+"Startlist1"=[bedstart1[x]]
            # bedname1[x]+"Endlist1"=bedend1[x]
            dictrange[bedname1[x]] = [bedstart1[x], bedend1[x]]
        return dictrange, bedname1