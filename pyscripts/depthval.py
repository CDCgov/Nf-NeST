class depthval:
    #########################
    #### this is where filtered spread file is used 
    #### along with DP and DP4 from to filter out variants
    #### that does not have DP4 or DP information
    ###########

    def __init__(self, filtered_path):
        self.filtered_path = filtered_path
        return

    def depthvalprocess(self):
        f3 = open(self.filtered_path, "r")
        DepthVal = {}
        DepthVal2 = {}
        DepthVal3 = {}
        for line in f3:
            # print(line)
            count = 0
            tempword = []
            tempword2 = []
            tempword22 = []
            tempword3 = []
            tempsumall = 0
            tempsumalt = 0
            tempgeneval1 = ""
            tempposval1 = ""
            currenttool=""
            for word in line.split():
                count += 1
                if count ==7:
                    currenttool=word
                if count == 8:
                    for x in word.split(","):
                        tempword += [x]
                if count == 9 and len(tempword) == 1:
                    for x in word.split(","):
                        tempword2 += [x]
                if count == 9 and len(tempword) == 4:
                    for x in word.split(","):
                        tempword22 += [x]
                if count == 10 and len(tempword22) == 1:
                    for x in word.split(","):
                        tempword3 += [x]
                if count == 1:
                    tempgeneval1 = word
                if count == 2:
                    tempposval1 = word
                #print(tempgeneval1,tempposval1)
            #print(tempword)
            #print(tempword2)
            #print(tempgeneval1,tempposval1)
            #print(tempword)
            #print(tempword2)
            #print(tempgeneval1,tempposval1)
            if len(tempword) == 4:
                #print(tempgeneval1,tempposval1)
                for x in range(len(tempword)):
                    tempsumall += int(tempword[x])
                    if x > 1:
                        tempsumalt += int(tempword[x])
                # if tempsumalt>=tempsumall:
                #    print(tempsumalt)
                #    print(tempsumall)
                #print(tempgeneval1,tempposval1)
                DepthVal[tempgeneval1, tempposval1] = [tempsumalt / tempsumall]
                DepthVal2[tempgeneval1, tempposval1] = [tempsumalt / tempsumall]
                if len(tempword22) == 1 and len(tempword3) == 2:
                    # print(int(tempword3[1])/int(tempword22[0]))
                    # print(tempword3[1])
                    # print(tempword22[0])
                    DepthVal3[tempgeneval1, tempposval1] = [
                        int(tempword3[1]) / int(tempword22[0])
                    ]
            if currenttool=="Freebayes":
                #print(tempword)
                DepthVal[tempgeneval1, tempposval1] = tempword
            if len(tempword) == 1 and len(tempword2) == 2:
                # print(tempword2[1])
                # print(tempword[0])
                #print(tempword)
                #print(tempword2)
                #print(line)
                #print(tempgeneval1,tempposval1)
                DepthVal[tempgeneval1, tempposval1] = [int(tempword2[1]) / int(tempword[0])]
                DepthVal3[tempgeneval1, tempposval1] = [int(tempword2[1]) / int(tempword[0])]
        return DepthVal, DepthVal3