import pandas as pd
import os
import subprocess
import glob
import argparse
import matplotlib.pyplot as plt
import dexplot as dxp
import seaborn as sns

parser = argparse.ArgumentParser(description='name')
parser.add_argument('-f', dest='snpfilter', type=str, help="snpfilter output")
parser.add_argument('-voi', dest='voifile',type=str,help="voi files")
parser.add_argument('-d',dest='directory',type=str, help="snpfilter output")
parser.add_argument('-cod', dest='cutoff_dir', type=str, help="cutoff files directory")
parser.add_argument('-o', dest='output', type=str, help="name of output")
args = parser.parse_args()
all_df=pd.DataFrame()

for a in glob.glob(args.snpfilter):
       if not a.startswith("list"):
              name=a.split(".csv")[0]
              a=pd.read_csv(a)
              a['Sample name']=name
              all_df=all_df.append(a)

all_df.AAPOS=all_df.AAPOS.astype(str)
all_df['Reportable mutation with gene']=all_df['Gene']+"_"+all_df['AAref']+all_df['AAPOS']+all_df['AAalt']
all_df['Reportable mutation']=all_df['AAref']+all_df['AAPOS']+all_df['AAalt']
all_df.Reportable=all_df.Reportable.str.replace('no','No')
all_df.Reportable=all_df.Reportable.str.replace('reportable','Yes')
all_df.Candidates=all_df.Candidates.str.replace('no','No')
all_df.Candidates=all_df.Candidates.str.replace('candidates','Yes')
all_df['first']=all_df['AAref']+all_df['AAPOS'].astype('str')

voi=pd.read_csv(args.voifile)
voi['first']=voi["RefAA"]+voi["AAPos"].astype('str')
merge1=pd.merge(all_df,voi,on=['first'])
merge1=merge1.drop_duplicates()

merge1.columns=['Gene','BasePOS', 'BaseDepth', "Agreents", 'Ref', 'Alt', 'AAref', 'AAalt',
       'AAPOS', 'CodonCoverage', 'VAF(DP4)','AF','Mutation', 'QD', 'SOR', 'MQ',
       'MQRankSum', 'Filter', 'FilterDescription', 'Candidates', 'Reportable',
       'Sample name', 'Reportable mutation with gene','Reportable mutation', 'first','Chr','Gene_y', 
       'RefAA', 'AAPos', 'AltAA']
merge1['AAPOS']=merge1['AAPOS'].astype(str)
#print(merge1['AAref'])
#merge1['Reportable mutation']=merge1['AAref']+merge1['AAPOS']+merge1['AAalt']
merge1['Reportable mutation']=''
for i in range(0,len(merge1)):
    if merge1.Mutation[i]=="Wildtype":
      #print(merge1['AAPOS'][i])
      #print(merge1['AAref'][i])
      merge1['Reportable mutation'][i]=merge1['AAPOS'][i]+merge1['AAref'][i]
    else:
        merge1['Reportable mutation'][i]=merge1['AAPOS'][i]+merge1['AAref'][i]+">"+merge1['AAalt'][i]
merge1['Drug Resistance Marker']=merge1['AAPos'].astype(str)+merge1['RefAA'].astype(str)+">"+merge1['AltAA'].astype(str)
merge1['Reportable mutation with gene']=merge1['Gene']+"_"+merge1['AAref']+merge1['AAPOS'].astype("str")+merge1['AAalt']
merge1=merge1[['Gene','BasePOS', 'BaseDepth', "Agreents", 'Ref', 'Alt', 'AAref', 'AAalt',
       'AAPOS', 'CodonCoverage', 'VAF(DP4)','AF','Mutation', 'QD', 'SOR', 'MQ',
       'MQRankSum', 'Filter', 'FilterDescription', 'Candidates', 'Reportable',
       'Sample name', 'Reportable mutation with gene','Reportable mutation','Drug Resistance Marker']]

summary=merge1[['Sample name','Gene','Ref', 'Alt', 'AAref', 'AAalt',
       'AAPOS', 'CodonCoverage', 'VAF(DP4)','Mutation','Drug Resistance Marker','QD', 
       'SOR', 'MQ','MQRankSum', 'Filter', 'FilterDescription']]


merge1_AAonly=merge1[['Gene','AAref', 'AAalt', 'AAPOS',
       'CodonCoverage', 'AF', 'Mutation', 'QD', 'SOR', 'MQ', 'MQRankSum',
       'Filter', 'FilterDescription', 'Candidates', 'Reportable',
       'Sample name','VAF(DP4)','Reportable mutation', 'Drug Resistance Marker',
       ]]
merge1_AAonly.columns=['Gene','AAref', 'AAalt', 'AAPOS',
       'CodonCoverage', 'AF', 'Mutation', 'QD', 'SOR', 'MQ', 'MQRankSum',
       'Filter', 'FilterDescription', 'Candidates', 'Reportable',
       'Sample name','VAF(DP4)','Reportable mutation', 'Drug Resistance Marker',
       ]
merge1_AAonly['VAF']=merge1_AAonly['VAF(DP4)'].astype("float",errors='ignore')
merge1_AAonly['Filter']=merge1_AAonly['Filter'].astype("float",errors='ignore')
merge1_AAonly=merge1_AAonly.drop_duplicates()
#merge1_AAonly.to_csv("BI_input.csv")

haplo=merge1_AAonly[['Sample name','Gene','AAalt','AAPOS','Reportable mutation','VAF']]
pd.to_numeric(haplo.VAF,errors='ignore')
pd.to_numeric(haplo.VAF,errors='coerce')

haplo=haplo[haplo.VAF >= 0.95]
haplo.AAPOS=haplo.AAPOS.astype('int')
haplo=haplo.drop_duplicates()
#print(haplo['Sample name'])
h=haplo.pivot_table(columns=['Reportable mutation'],values='AAalt',index=['Sample name','Gene'],aggfunc='first')
h=h.sort_index()
h.columns=h.columns.sort_values()
#h_list=["C72S","V73V","M74I","N75E","K76T","N86Y","Y184F","S1034C","N1042D","D1246Y","S436A","A437G","K540E","A581G","A613S","N51I","C59R","S108N",
#       "I258M","Y268S","N458Y","Y493H","R539T","I543I","R561H","C580Y"]
#h_list=["72C>S","73V>V","74M>I","75N>E","76K>T","86N>Y","184Y>F","1034S>C","1042N>D","1246D>Y","436S>A","437A>G","540K>E","581A>G",
#"613A>S","51N>I","59C>R","108S>N","258I>M","268Y>S","458N>Y","493Y>H","539R>T","543I>I","561R>H","580C>Y"]
h_list=["72C>S","73V>V","74M>I","75N>E","76K>T","86N>Y","184Y>F","1034S>C","1042N>D","1246D>Y","436S>A","437A>G","540K>E","581A>G","613A>S","51N>I","59C>R","108S>N"]
import numpy as np
for col in h_list:
    if col not in h.columns:
        h[col] = np.nan
#h=h[["C72S","V73V","M74I","N75E","K76T","N86Y","Y184F","S1034C","N1042D","D1246Y","S436A","A437G","K540E","A581G","A613S","N51I","C59R","S108N"]]
h=h[["72C>S","73V>V","74M>I","75N>E","76K>T","86N>Y","184Y>F","1034S>C","1042N>D","1246D>Y","436S>A","437A>G","540K>E","581A>G","613A>S","51N>I","59C>R","108S>N"]]
h['haplotype']=pd.Series(h.fillna('').values.tolist(),index=h.index).map(lambda x: ''.join(map(str,x)))
h_reset=h.reset_index()
haplo_table=pd.merge(summary,h_reset,on=["Sample name","Gene"])
haplo_table=haplo_table[['Sample name','Gene','haplotype']]
haplo_table=haplo_table.drop_duplicates()
haplo_table.to_csv(args.output+"_haplo.csv")