class finaldoutput:

    def __init__(self, name):
        self.name=name
        return

    def finaldoutputprocess(self,filteredGenelist22,filteredPoslist22,filtereddepthlist2up2,filteredAGlist12,filteredReflist22,filteredAltlist22
        ,filteredAAlist212,filteredAAlist222,filteredAAPoslist122,filteredcodondepthlist22,filteredDP4list12,filteredAFlist12,filteredmutationlist12,
        filteredQDlist12,filteredSORlist12,filteredMQlist12,filteredMQRankSumlist12,filteredfilterscorelist12,filteredDescriptionlist12,filteredcandidateslist12,
        filteredreportablelist12):
        import pandas as pd

        d3 = {
            "Gene": filteredGenelist22,
            "BasePOS": filteredPoslist22,
            "BaseDepth": filtereddepthlist2up2,
            "Agreents": filteredAGlist12,
            "Ref": filteredReflist22,
            "Alt": filteredAltlist22,
            "AAref": filteredAAlist212,
            "AAalt": filteredAAlist222,
            "AAPOS": filteredAAPoslist122,
            "CodonCoverage": filteredcodondepthlist22,
            "VAF(DP4)": filteredDP4list12,
            "AF": filteredAFlist12,
            "Mutation": filteredmutationlist12,
            "QD": filteredQDlist12,
            "SOR": filteredSORlist12,
            "MQ": filteredMQlist12,
            "MQRankSum": filteredMQRankSumlist12,
            "Filter": filteredfilterscorelist12,
            "FilterDescription": filteredDescriptionlist12,
            "Candidates": filteredcandidateslist12,
            "Reportable": filteredreportablelist12,
        }
        df3 = pd.DataFrame(data=d3)
        # df.to_csv('snpreport4.csv', index=False,)
        df3 = df3.sort_values(df3.columns[0], ascending=False)
        df3.to_csv(
            self.name + ".csv",
            index=False,
        )