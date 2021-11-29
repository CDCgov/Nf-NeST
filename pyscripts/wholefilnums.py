class wholefilnums:

    def __init__(self, filtered_path):
        self.filtered_path = filtered_path
        return

    def wholefilnumsprocess(self, AAdic):
        f3 = open(self.filtered_path, "r")

        wholefilnumsx = []
        wholefilnums = []
        
        for line in f3:
            count = 0
            tempnum = ""
            # print(line)
            if line.find("p.") != -1 and line.startswith("#") == False:
                for word in line.split():
                    count += 1
                    if count == 2:
                        tempnum = word
                    # print(tempnum)
                    # print(word)
                    # print(tempnum)
                    if word.startswith("c."):
                        if word.find("*") != -1:
                            wholefilnumsx += [tempnum]
                            # print(word[word.find("c.*"):word.find(tempnum)+len(tempnum)])
                            wholefilnums += [
                                word[word.find("c.*") : word.find(tempnum) + len(tempnum)]
                            ]
                        elif word.find("*") == -1:
                            # print(word[word.find("c.*"):word.find(tempnum)+len(tempnum)])
                            wholefilnumsx += [tempnum]
                            wholefilnums += [
                                word[word.find("c.") : word.find(tempnum) + len(tempnum)]
                            ]
                    if word.startswith("n."):
                        if word.find("*") != -1:
                            wholefilnumsx += [tempnum]
                            # print(word[word.find("c.*"):word.find(tempnum)+len(tempnum)])
                            wholefilnums += [
                                word[word.find("n.*") : word.find(tempnum) + len(tempnum)]
                            ]
                        elif word.find("*") == -1:
                            # print(word[word.find("c.*"):word.find(tempnum)+len(tempnum)])
                            wholefilnumsx += [tempnum]
                            wholefilnums += [
                                word[word.find("n.") : word.find(tempnum) + len(tempnum)]
                            ]
                    if word.startswith("p."):
                        # print(word[5+len(temppos):5+len(temppos)+3])
                        # print(word[5+len(temppos):5+len(temppos)+3])
                        temppos = ""
                        for x in word:
                            # print(word)
                            if x.isdigit():
                                # print(x)
                                temppos += x
                        if (
                            word[2:5] not in AAdic
                            or word[5 + len(temppos) : 5 + len(temppos) + 3] not in AAdic
                        ):
                            # print("working")
                            if word[2:5] in AAdic and word[-1:] == "*":
                                break
                            wholefilnumsx = wholefilnumsx[:-1]
                            wholefilnums = wholefilnums[:-1]
        return wholefilnumsx, wholefilnums