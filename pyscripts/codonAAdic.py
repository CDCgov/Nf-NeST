class codonAAdic:

    def __init__(self):
        return

    def codonAAdicprocess(self):
        import itertools
        ###### Assign letters to the amino acid names #########
        AAdic = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Asx": "B",
        "Cys": "C",
        "Glu": "E",
        "Gln": "Q",
        "Glx": "Z",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
        "STOP": "X",
        }
        #######Convert combination of bases to a codon############
        currentcodon = ""
        keywords = ["".join(i) for i in itertools.product(["T", "G", "A", "C"], repeat=3)]
        codondic = {}
        for x in keywords:
            currentcodon = ""
            for y in x:
                currentcodon += y
            # print(currentcodon)
            if currentcodon == "TTT" or currentcodon == "TTC":
                # print("True")
                codondic[currentcodon] = "Phe"
            if currentcodon == "TTA" or currentcodon == "TTG":
                codondic[currentcodon] = "Leu"
            if (
                currentcodon == "CTT"
                or currentcodon == "CTC"
                or currentcodon == "CTA"
                or currentcodon == "CTG"
            ):
                codondic[currentcodon] = "Leu"
            if currentcodon == "ATT" or currentcodon == "ATC" or currentcodon == "ATA":
                codondic[currentcodon] = "Ile"
            if currentcodon == "ATG":
                codondic[currentcodon] = "Met"
            if (
                currentcodon == "GTT"
                or currentcodon == "GTC"
                or currentcodon == "GTA"
                or currentcodon == "GTG"
            ):
                codondic[currentcodon] = "Val"
            if (
                currentcodon == "TCT"
                or currentcodon == "TCC"
                or currentcodon == "TCA"
                or currentcodon == "TCG"
            ):
                codondic[currentcodon] = "Ser"
            if (
                currentcodon == "CCT"
                or currentcodon == "CCC"
                or currentcodon == "CCA"
                or currentcodon == "CCG"
            ):
                codondic[currentcodon] = "Pro"
            if (
                currentcodon == "ACT"
                or currentcodon == "ACC"
                or currentcodon == "ACA"
                or currentcodon == "ACG"
            ):
                codondic[currentcodon] = "Thr"
            if (
                currentcodon == "GCT"
                or currentcodon == "GCC"
                or currentcodon == "GCA"
                or currentcodon == "GCG"
            ):
                codondic[currentcodon] = "Ala"
            if currentcodon == "TAT" or currentcodon == "TAC":
                codondic[currentcodon] = "Tyr"
            if currentcodon == "CAT" or currentcodon == "CAC":
                codondic[currentcodon] = "His"
            if currentcodon == "CAA" or currentcodon == "CAG":
                codondic[currentcodon] = "Gln"
            if currentcodon == "AAT" or currentcodon == "AAC":
                codondic[currentcodon] = "Asn"
            if currentcodon == "AAA" or currentcodon == "AAG":
                codondic[currentcodon] = "Lys"
            if currentcodon == "GAT" or currentcodon == "GAC":
                codondic[currentcodon] = "Asp"
            if currentcodon == "GAA" or currentcodon == "GAG":
                codondic[currentcodon] = "Glu"
            if currentcodon == "TGT" or currentcodon == "TGC":
                # print("I am correct")
                codondic[currentcodon] = "Cys"
            if currentcodon == "TGG":
                codondic[currentcodon] = "Trp"
            if (
                currentcodon == "AGA"
                or currentcodon == "AGG"
                or currentcodon == "CGT"
                or currentcodon == "CGC"
                or currentcodon == "CGA"
                or currentcodon == "CGG"
            ):
                codondic[currentcodon] = "Arg"
            if currentcodon == "AGT" or currentcodon == "AGC":
                codondic[currentcodon] = "Ser"
            if (
                currentcodon == "GGT"
                or currentcodon == "GGC"
                or currentcodon == "GGA"
                or currentcodon == "GGG"
            ):
                codondic[currentcodon] = "Gly"
            if currentcodon == "TGA" or currentcodon == "TAA" or currentcodon == "TAG":
                codondic[currentcodon] = "STOP"
        return AAdic, codondic