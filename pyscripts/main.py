####This is the main pyton script for running NfNeST. This automatically imports different steps of analyzing SNP####

###Main steps for the SNP analysis Pipeline ######

#1Generate POS/ALT/GENE/REF (Position of SNP, alternate base, name of Sample, reference base)

#2 SAMTOOLS Depth figure out (Figure out depth of the SNP like how many reads)

#3 DepthVal

#4 Genelist2/CurrentPOS2/GENE2/REF2/ALT2 (Filter Position, alternate base, name of sample, reference base based on "spread" and "final VCF")

#5 AF/MQ/QD/MQRankSUM/SOR (Figure out more information about alleles)

#6 filter number filter score (Add filter for quality of the SNP. Also, add explanation for the filter quality)

#7 dicttrans1 (Translate to obtain the CDS region and amino acid sequence part 1)

#8 dicttrans2 (Translate to obtain the CDS region and amino acid sequence part 2)

#9 coverage (Find the coverage for the genes for each samples)

#10 filterlists (Final lists for wildtypes and mutations)

#11 filteredlists2 (Final lists for wildtypes and mutations)


###### Import python libraries and scripts #####
import argparse
import pandas as pd
import subprocess
import csv
import itertools
###### Import dictionary for combination of bases to amino acid and combination of amino acids to codons###########
from codonAAdic import codonAAdic
###### Import script for Gene name, Position, Reference, Alternate bases #######
from GPRA import GPRA
###### Import script for finding depths for different positions of bases for each genes ######
from depthlist import depthlist
##### Obtain the depth of codon of interest #####
from codondepth import codondepthlist
##### Obtain the depth value #######
from depthval import depthval
from wholefilnums import wholefilnums
##### Obtain AF/MQ/QD/MQRankSUM/SOR ###### 
from properties import properties
from filterdes import filterdes
##### Build a dictionary for candidates #####
from candidates import candidates
from combinedf import combinedf
from bedrange import bedrange
from dicttrans import dicttrans
from dicttranstwo import dicttranstwo
from dicttransAA import dicttransAA
from potcoverage import potcoverage
from filteredlists import filteredlists
from posdicwi import posdicwi
from combinewholepot import combinewholepot
from candirepo import candirepo
from finaldoutput import finaldoutput
def install(name):
    subprocess.call(['pip3', 'install', name])

install('xlrd==1.2.0')
import xlrd


#def main(unfiltered_path,filtered_path,bamfile_path,bedfile_path,candidates_path,voi_path,fasta_path):
def main(arguments):
    #print(arguments)
    #print(arguments)
    #arguements=parse_args(arguements)
    ##########################################################################################
    #########import arguments#################
    unfiltered_path = arguments.unfiltered_path
    filtered_path = arguments.filtered_path
    #print(filtered_path)
    name = arguments.name
    bamfile_path = arguments.bamfile_path
    bedfile_path = arguments.bedfile_path
    candidates_path = arguments.candidates_path
    voi_path = arguments.voi_path
    fasta_path = arguments.fasta_path
    #########################################
    ##########################################################################################
    ######First class (Get Gene, Ref, POS, Alt lists)############
    #Use a filtered file (spread file) to obtain Gene name, Position, Reference base, and Alternate base for a given variant
    ########
    codonAA1=codonAAdic()
    AAdic, codondic = codonAA1.codonAAdicprocess()

    grpa=GPRA(filtered_path)
    Genelist1,POSlist1,Reflist1,Altlist1,currentag1=grpa.GPRAPROCESS(AAdic)
    #print(Genelist1,Reflist1,POSlist1,Altlist1)
    #print(POSlist1)
    ##########################################################################################
    #####Second class (Get depths) ##########################
    #####Obtain the base depths of different position and pair them up with Gene and Position of variants######
    #####basically creating a dictionary
    depthpro=depthlist(name)
    depthlist1up,depthpair1,depthlist1=depthpro.depthprocess(Genelist1, POSlist1)
    #print(depthlist1up,depthpair1)
    #########################################
    #print(depthlist1up,depthpair1)
    ##########################################################################################
    #print(depthlist1)
    #####Third class (Codondepthlist) ############
    ##########################################################################################
    codondepth1=codondepthlist()
    #print(codondepth1)
    #print(depthlist1)
    codondepthlist1=codondepth1.codondepprocess(Genelist1, POSlist1, depthpair1, depthlist1)
    #print(codondepthlist1)
    #codondepprocess(Genelist1, POSlist1, depthpair1, depthlist1)
    ##########################################################################################
    ##########################################################################################

    ######Fourth class Depthval###################
    depthval1=depthval(filtered_path)
    depthvallist1, dpethvallist31=depthval1.depthvalprocess()
    ##################################################


    ######FIFTH class wholefilnum1###################
    wholefilnum1=wholefilnums(filtered_path)
    wholefilnumsxlist1, wholefilnumslist1 = wholefilnum1.wholefilnumsprocess(AAdic)
    #print(wholefilnumsxlist1)
    #print(wholefilnumslist1)
    ##################################################


    ######Sixth class properties (AF,QD,SOR,MQ,MQRankSum,GenePos,DP4,AD)###################
    properties1=properties(unfiltered_path)
    Genelist2,POSlist2,AFlist1,QDlist1,SORlist1,MQlist1,MQRankSumlist1,GenePosSet1,DP4list1,ADlist1,currentag2 = properties1.propertiesprocess(Genelist1,POSlist1,depthvallist1,dpethvallist31,currentag1,wholefilnumslist1,wholefilnumsxlist1)
    #print(POSlist2)
    ##################################################

    ######Seventh class Filterdes with description###################
    filterdes1=filterdes()
    Descriptionlist1,filterscorelist1=filterdes1.filterdesprocess(QDlist1,SORlist1,MQlist1,MQRankSumlist1)
    #print(Descriptionlist1)

    ######Eigth class candidates###################
    candidates1=candidates(candidates_path,voi_path)
    df1Gene,df1RefAA,df1AltAA,df1AAPos,df2Gene,df2RefAA,df2AltAA,df2AAPos,Reflist2, Altlist2, depthlist2up, codondepthlist2=candidates1.candidatesprocess(Genelist1,Genelist2,POSlist1,POSlist2,Reflist1,Altlist1,depthlist1up,codondepthlist1)
    #print(df1Gene)

    ######9th class combinedf###################
    combinedf1=combinedf()
    combineddflist1,combineddflist2=combinedf1.combinedfprocess(df1Gene, df2Gene, df1AAPos, df2AAPos)
    #print(combineddflist2)

    ######10th class bedrange###################
    bedrange1=bedrange(bedfile_path)
    dictrange, bedname1=bedrange1.bedrangeprocess()
    #print(dictrange)
    ######################################

    ######11th class dictrange###################
    dicttrans1=dicttrans(fasta_path)
    dicttrans2,testdic,GenePosSet2=dicttrans1.dicttransprocess(bedname1, Genelist2, POSlist2, Altlist2)
    #print(dicttrans2)
    #print(dictrange)
    ######################################

    ######12th classs dicttranstwo##############
    dicttranstwo1=dicttranstwo()
    dicttrans22,testdicttrans22, BASEAAPOSdic1=dicttranstwo1.dicttranstwoprocess(bedname1, dicttrans2, testdic, dictrange,GenePosSet2)
    #print(dicttrans22)


    ######13th classs dicttransAA##############
    dicttransAAlist1=dicttransAA()
    dicttransAA1,testdicttransAA1=dicttransAAlist1.dicttransAAprocess(bedname1, dicttrans22, testdicttrans22, AAdic, codondic)
    #print(dicttransAA1)


    ######14th classs potcoverage##############
    potcoverage1=potcoverage(fasta_path, name)
    tempcov1,tempcovwordlist1, wholepotlist1 = potcoverage1.potcoverageprocess(df2Gene, df2AAPos)
    #print(tempcov1)

    ######15th classs filteredlists1##############
    filteredlists1=filteredlists()
    (filteredGenelist2, filteredPoslist2, filtereddepthlist2up, filteredReflist2, filteredAltlist2, filteredAltlist2, filteredAAlist21,
    filteredAAlist22, filteredAAPoslist12, filteredcodondepthlist2, filteredAFlist1, filteredmutationlist1, filteredQDlist1, filteredSORlist1,
    filteredMQlist1, filteredMQRankSumlist1, filteredfilterscorelist1, filteredDescriptionlist1, filteredcandidateslist1, filteredreportablelist1,
    filteredDP4list1, combinedGenePOSlist1,filteredAGlist2, filteredADlist1) = filteredlists1.filteredlistsprocess(Genelist2, POSlist2, Reflist2, Altlist2, BASEAAPOSdic1, combineddflist1,combineddflist2, 
    currentag2, depthlist2up, AFlist1, testdicttransAA1, dicttransAA1, DP4list1, ADlist1, QDlist1, SORlist1, MQlist1, 
    MQRankSumlist1, filterscorelist1, Descriptionlist1, codondepthlist2)
    #print(filteredAAlist21,filteredAAlist22)
    #print(filteredAAlist22)
    

    ######16th classs posdic/wildic##############
    posdicwi1 = posdicwi()
    posdiclist1, wildic1 = posdicwi1.posdicprocess(bedname1, wholepotlist1, dictrange, depthlist1)
    #print(posdiclist1)

    ######17th classs posdic/wildic##############
    combinewholepot1 = combinewholepot()
    (filteredGenelist22,filteredPoslist22,filtereddepthlist2up2,filteredReflist22,filteredAltlist22,filteredAAlist212,
    filteredAAlist222,filteredAAPoslist122,filteredcodondepthlist22,filteredAFlist12,filteredmutationlist12,filteredQDlist12,
    filteredSORlist12,filteredMQlist12,filteredMQRankSumlist12,filteredfilterscorelist12,filteredDescriptionlist12,
    filteredcandidateslist12,filteredreportablelist12,filteredDP4list12,filteredADlist12,filteredAGlist12) = combinewholepot1.combinewholepotprocess(filteredADlist1,filteredAGlist2,filteredGenelist2, filteredPoslist2, filtereddepthlist2up, filteredReflist2, filteredAltlist2, filteredAAlist21,
    filteredAAlist22, filteredAAPoslist12, filteredcodondepthlist2, filteredAFlist1, filteredmutationlist1, filteredQDlist1, filteredSORlist1,
    filteredMQlist1, filteredMQRankSumlist1, filteredfilterscorelist1, filteredDescriptionlist1, filteredcandidateslist1, filteredreportablelist1,
    filteredDP4list1, tempcovwordlist1, combinedGenePOSlist1, testdicttransAA1, dicttransAA1, wholepotlist1, wildic1 ,posdiclist1)
    #print(filteredAAlist212,filteredAAlist222)
    #print(filteredAAlist222)
    #print(testdicttransAA1["mitochondrial_genome"])
    #print(dicttransAA1["mitochondrial_genome"])

    ######18th classs candidates/reportables##############
    candirepo1 = candirepo()
    filteredcandidateslist12, filteredreportablelist12 = candirepo1.candirepoprocess(filteredGenelist22,filteredAAPoslist122,filteredAAlist212,filteredAAlist222,df1AAPos,df1RefAA,df1AltAA,df2AAPos, df2RefAA, df2AltAA, df2Gene, df1Gene)
    #print(filteredcandidateslist12)
    #print(Genelist1)
    #print(POSlist1)
    #print(filteredAAPoslist12)
    ######19th classs candidates/reportables##############
    finaldoutput1 = finaldoutput(name)
    finaldoutput1.finaldoutputprocess(filteredGenelist22,filteredPoslist22,filtereddepthlist2up2,filteredAGlist12,filteredReflist22,filteredAltlist22
        ,filteredAAlist212,filteredAAlist222,filteredAAPoslist122,filteredcodondepthlist22,filteredDP4list12,filteredAFlist12,filteredmutationlist12,
        filteredQDlist12,filteredSORlist12,filteredMQlist12,filteredMQRankSumlist12,filteredfilterscorelist12,filteredDescriptionlist12,filteredcandidateslist12,
        filteredreportablelist12)
    return(0)

#########Set up parameters to run the command ##############
def parse_arguments():
    parser = argparse.ArgumentParser(description="name")
    #############Import final_vcf file########
    parser.add_argument(
        "-v1", dest="unfiltered_path", type=str, help="name of unfilterd merged vcf file"
    )
    #############Import spread file#######
    parser.add_argument(
        "-v2", dest="filtered_path", type=str, help="name of filtered merged vcf file"
    )
    parser.add_argument("-b1", dest="bamfile_path", type=str, help="name of bam file")
    parser.add_argument("-b2", dest="bedfile_path", type=str, help="name of bed file")
    parser.add_argument("-e1", dest="candidates_path", type=str, help="name of candidates")
    parser.add_argument("-e2", dest="voi_path", type=str, help="name of variants of interest")
    parser.add_argument("-f1", dest="fasta_path", type=str, help="name of fasta file")
    parser.add_argument("-o1", dest="name", type=str, help="name of output")
    args = parser.parse_args()
    #print(args)
    return args



if __name__ == '__main__':
    arguments = parse_arguments()
    #main(args.unfiltered_path,args.filtered_path,args.bamfile_path,args.bedfile_path,args.candidates_path,args.voi_path,args.fasta_path)
    #print(arguments)
    main(arguments)