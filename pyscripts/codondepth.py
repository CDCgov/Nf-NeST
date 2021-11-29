class codondepthlist:

    def __init__(self):
    #    self.name = name
        return

    def codondepprocess(self, Genelist1, POSlist1, depthpair1, depthlist1):
        codondepthlist1 = []
        # sadcount=0
        for x in range(len(Genelist1)):
            if [Genelist1[x], POSlist1[x]] not in depthpair1:
                # print("True")
                codondepthlist1 += [None]
                # break
            else:
                for y in range(len(depthlist1)):
                    if Genelist1[x] in depthlist1[y] and POSlist1[x] in depthlist1[y]:
                        # print(Genelist1[x],POSlist1[x])
                        # print(depthlist1[y])
                        # sadcount+=1
                        tempcodondepthlist1 = 0
                        if Genelist1[x] == "PfCRT":
                            # for z in range(len(PfCRTStartlist1)):
                            # if PfCRTStartlist1[z]<int(depthlist1[y-1][2])<PfCRTEndlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y - 1][2])
                            # if PfCRTStartlist1[z]<int(depthlist1[y][2])<PfCRTEndlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y][2])
                            # if PfCRTStartlist1[z]<int(depthlist1[y+1][2])<PfCRTEndlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y + 1][2])
                            # print(tempcodondepthlist1)
                            codondepthlist1 += [int(tempcodondepthlist1 / 3)]
                            break
                        if Genelist1[x] == "DHFR":
                            # for z in range(len(DHFRStartlist1)):
                            # print(type(DHFRStartlist1[z]))
                            # print((depthlist1[y-1][2]))
                            # if DHFRStartlist1[z]<int(depthlist1[y-1][2])<DHFREndlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y - 1][2])
                            # if DHFRStartlist1[z]<int(depthlist1[y][2])<DHFREndlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y][2])
                            # if DHFRStartlist1[z]<int(depthlist1[y+1][2])<DHFREndlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y + 1][2])
                            codondepthlist1 += [int(tempcodondepthlist1 / 3)]
                            break
                        if Genelist1[x] == "K13":
                            # for z in range(len(K13Startlist1)):
                            # if K13Startlist1[z]<int(depthlist1[y-1][2])<K13Endlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y - 1][2])
                            # if K13Startlist1[z]<int(depthlist1[y][2])<K13Endlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y][2])
                            # if K13Startlist1[z]<int(depthlist1[y+1][2])<K13Endlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y + 1][2])
                            codondepthlist1 += [int(tempcodondepthlist1 / 3)]
                            break
                        if Genelist1[x] == "PfMDR1":
                            #print(y)
                            #print(len(depthlist1))
                            #print(depthlist1[y+1])
                            # for z in range(len(PfMDR1Startlist1)):
                            # if PfMDR1Startlist1[z]<int(depthlist1[y-1][2])<PfMDR1Endlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y - 1][2])
                            # if PfMDR1Startlist1[z]<int(depthlist1[y][2])<PfMDR1Endlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y][2])
                            # if PfMDR1Startlist1[z]<int(depthlist1[y+1][2])<PfMDR1Endlist1[z]:
                            if y<len(depthlist1)-1:
                                tempcodondepthlist1 += int(depthlist1[y + 1][2])
                            codondepthlist1 += [int(tempcodondepthlist1 / 3)]
                        if Genelist1[x] == "DHPS":
                            # for z in range(len(PfDHPSStartlist1)):
                            # if PfDHPSStartlist1[z]<int(depthlist1[y-1][2])<PfDHPSEndlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y - 1][2])
                            # if PfDHPSStartlist1[z]<int(depthlist1[y][2])<PfDHPSEndlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y][2])
                            # if PfDHPSStartlist1[z]<int(depthlist1[y+1][2])<PfDHPSEndlist1[z]:
                            if y + 1 < len(depthlist1):
                                tempcodondepthlist1 += int(depthlist1[y + 1][2])
                            codondepthlist1 += [int(tempcodondepthlist1 / 3)]
                            break
                        if Genelist1[x] == "mitochondrial_genome":
                            # for z in range(len(mitochondrialgenomeStartlist1)):
                            # if mitochondrialgenomeStartlist1[z]<int(depthlist1[y-1][2])<mitochondrialgenomeEndlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y - 1][2])
                            # if mitochondrialgenomeStartlist1[z]<int(depthlist1[y][2])<mitochondrialgenomeEndlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y][2])
                            # if mitochondrialgenomeStartlist1[z]<int(depthlist1[y+1][2])<mitochondrialgenomeEndlist1[z]:
                            tempcodondepthlist1 += int(depthlist1[y + 1][2])
                            codondepthlist1 += [int(tempcodondepthlist1 / 3)]
                            break
        return codondepthlist1