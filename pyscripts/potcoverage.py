class potcoverage:
    def __init__(self,fasta_path, name):
        self.fasta_path=fasta_path
        self.name=name
        return

    def potcoverageprocess(self, df2Gene, df2AAPos):
        import subprocess
        import csv
        wholepotlist1 = []
        # print(df1Gene)
        for k in range(len(df2Gene)):
            if df2Gene[k] == "Pfcrt":
                wholepotlist1.append(["PfCRT", str(df2AAPos[k])])
            elif df2Gene[k] == "Pfdhfr":
                wholepotlist1.append(["DHFR", str(df2AAPos[k])])
            elif df2Gene[k] == "Pfk13":
                wholepotlist1.append(["K13", str(df2AAPos[k])])
            elif df2Gene[k] == "Pfmdr1":
                wholepotlist1.append(["PfMDR1", str(df2AAPos[k])])
            elif df2Gene[k] == "Pfdhps":
                wholepotlist1.append(["DHPS", str(df2AAPos[k])])
            elif df2Gene[k] == "MT":
                wholepotlist1.append(["mitochondrial_genome", str(df2AAPos[k])])


        # print(wholepotlist1)
        # print(x[0])
        # print(x[1])
        # print(AAuntranslatewilddic1[x[0],x[1]])
        # print(AAtranslatewilddic1[x[0],x[1]])

        # print(wholepotlist1)

        totalgene = []
        with open(self.fasta_path, "r") as f1:
            for line in f1:
                if line.startswith(">"):
                    totalgene += [(line.strip(">")).strip("\n")]

        # print(totalgene)


        # samtools index ${bam_path}
        result = subprocess.Popen(["samtools", "index", self.name + "_SR.bam"])
        process = subprocess.Popen(
            [
                "samtools",
                "coverage",
                self.name + "_SR.bam",
                "-r",
                "mitochondrial_genome:3492-4622",
            ],
            stdout=subprocess.PIPE,
        )
        stdout, stderr = process.communicate()
        reader2 = csv.reader(
            stdout.decode("ascii").splitlines(), delimiter="\t", skipinitialspace=True
        )

        tempcov1 = []
        tempcovwordlist1 = []

        for line in reader2:
            count = 0
            tempcovword1 = ""
            # print(line)
            if line[0].startswith("#") == False:
                for word in line:
                    count += 1
                    if count == 1:
                        tempcovword1 = word
                    # print(word)
                    if count == 6 and float(word) > 0:
                        # print(word)
                        tempcov1 += [word]
                        tempcovwordlist1 += [tempcovword1]

        # print(tempcovwordlist1)

        result = subprocess.Popen(["samtools", "index", self.name + "_SR.bam"])
        process = subprocess.Popen(
            ["samtools", "coverage", self.name + "_SR.bam"], stdout=subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        reader2 = csv.reader(
            stdout.decode("ascii").splitlines(), delimiter="\t", skipinitialspace=True
        )

        for line in reader2:
            count = 0
            tempcovword1 = ""
            # print(line)
            if line[0].startswith("#") == False:
                for word in line:
                    count += 1
                    if count == 1:
                        tempcovword1 = word
                    # print(word)
                    if (
                        count == 6
                        and float(word) > 0
                        and (tempcovword1 not in tempcovwordlist1)
                    ):
                        # print(word)
                        # print(word)
                        # print(tempcovword1)
                        tempcov1 += [word]
                        tempcovwordlist1 += [tempcovword1]
        return tempcov1,tempcovwordlist1, wholepotlist1