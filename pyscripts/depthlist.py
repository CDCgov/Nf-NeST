import subprocess
import csv

class depthlist:
    #Obtain Gene, Position, Alt

    def __init__(self, name):
        self.name = name
        return

    def depthprocess(self, Genelist1, POSlist1):
        process = subprocess.Popen(
                ["samtools", "depth", "-a", self.name + "_SR.bam"], stdout=subprocess.PIPE
            )
        stdout, stderr = process.communicate()
        reader = csv.reader(
            stdout.decode("ascii").splitlines(), delimiter="\t", skipinitialspace=True
        )
        depthlist1 = []
        for row in reader:
            #print(row)
            if row[2] != "0":
                # print(row)
                depthlist1 += [row]
        depthlist1up = []
        # for x in range(len(Genelist1)):
        #    if Genelist1[x] in depthlist1:
        #        print(Genelist1[x])
        #    if POSlist1[x] in depthlist1:
        #        print(POSlist1[x])
        # print(depthlist1)
        # print(depthlist1)
        # print(Genelist1)
        # print(POSlist1)
        # print(Genelist1)
        depthpair1 = []
        #   print(POSlist1)
        for x in range(len(depthlist1)):
            depthpair1 += [[depthlist1[x][0], depthlist1[x][1]]]

        for x in range(len(Genelist1)):
            if [Genelist1[x], POSlist1[x]] not in depthpair1:
                depthlist1up += [None]

            else:
                for y in range(len(depthlist1)):
                    if Genelist1[x] in depthlist1[y] and POSlist1[x] in depthlist1[y]:
                        # print(depthlist1[y])
                        # print(depthlist1[y][2])
                        depthlist1up += [depthlist1[y][2]]
                        break
        return depthlist1up, depthpair1, depthlist1