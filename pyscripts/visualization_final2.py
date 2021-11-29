
import pandas as pd
import os
import subprocess
import glob
import argparse
import matplotlib.pyplot as plt


def install(name):
    subprocess.call(['pip3', 'install', name])

install('seaborn')
install('dexplot')

import dexplot as dxp
import seaborn as sns
parser = argparse.ArgumentParser(description='name')
parser.add_argument('-f', dest='snpfilter', type=str, help="snpfilter output")
parser.add_argument('-voi', dest='voifile',type=str,help="voi files")
parser.add_argument('-d',dest='directory',type=str, help="snpfilter output")
parser.add_argument('-cod', dest='cutoff_dir', type=str, help="cutoff files directory")
args = parser.parse_args()
all_df=pd.DataFrame()
dic={'ANZa':(-6.2667,14.2333),'ANLS':(-10.2592,20.6249),'ANBe':(-12.5763,13.4055),'GNMa':(9.546,-13.281),'GNHa':(11.56667,-10.05),'GNDo':(7.75833,-8.80472),'GNLS':(11.03856,-12.04793)}
dic2={'ANBe':'Benguela','ANZa':'Zaire','ANLS':'Lunda Sul','GNLa':'Labé','GNHa':'Dabola','GNDo':'Nzérékoré','GNMa':'Forécariah Centre','ANxx':'control'}
dic3={'A':'AL','B':'CQ','C':'CQ+PQ','D':'AL+PQ','E':'AL+MQ','F':'AS+AQ','G':'AS+SP','H':'AS+MQ','I':'DHA+PQ','J':'AS+SP+PQ','K':'AL+MQ+PQ','L':'Doxxycline','M':'Malarone','N':'Quinine','O':'Quinine+Doxycycline','P':'Atovaquone/Proguanil+AL','Q':'Mefloquine','R':'Vancomycin+Rocephin','S':'Malarone+Doxycycline','x':'None'}


for a in glob.glob("*.csv"):
       if not a.startswith("list"):
           if not a.startswith("21US"):
                name=a.split(".csv")[0]
                a=pd.read_csv(a)
                a['Sample name']=name
                a['region']=dic2[name[2:6]]
                a['Year']='20'+name[0:2]
                a['Treatment']=dic3[name[8:9]]
                a['Day']=name[6:8]
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

merge1.columns=['Gene', 'BasePOS', 'BaseDepth', 'Agreents', 'Ref', 'Alt', 'AAref',
       'AAalt', 'AAPOS', 'CodonCoverage', 'VAF(DP4)', 'AF', 'Mutation', 'QD',
       'SOR', 'MQ', 'MQRankSum', 'Filter', 'FilterDescription', 'Candidates',
       'Reportable', 'Sample name', 'region', 'Year', 'Treatment', 'Day',
       'Reportable mutation with gene', 'Reportable mutation', 'first', 'Chr',
       'Gene_y', 'RefAA', 'AAPos', 'AltAA']
merge1['AAPOS']=merge1['AAPOS'].astype(str)
merge1=merge1.reset_index()
#merge1['Reportable mutation']=merge1['AAref']+merge1['AAPOS']+merge1['AAalt']
merge1['Reportable mutation']=''
for i in range(0,len(merge1)):
    if merge1.Mutation[i]=="Wildtype":
        merge1['Reportable mutation'][i]=str(merge1['AAPOS'][i])+str(merge1['AAref'][i])
    else:
        merge1['Reportable mutation'][i]=str(merge1['AAPOS'][i])+merge1['AAref'][i]+">"+merge1['AAalt'][i]
merge1['Drug Resistance Marker']=merge1['AAPos'].astype(str)+merge1['RefAA'].astype(str)+">"+merge1['AltAA'].astype(str)
merge1['Reportable mutation with gene']=merge1['Gene']+"_"+merge1['AAref']+merge1['AAPOS'].astype("str")+merge1['AAalt']
merge1=merge1[['Gene','BasePOS', 'BaseDepth', 'Ref', 'Alt', 'AAref', 'AAalt',
       'AAPOS', 'CodonCoverage', 'VAF(DP4)','AF','Mutation', 'QD', 'SOR', 'MQ',
       'MQRankSum', 'Filter', 'FilterDescription', 'Candidates', 'Reportable',
       'Sample name', 'Reportable mutation with gene','Reportable mutation','Drug Resistance Marker','region','Treatment','Day']]

summary=merge1[['Sample name','Gene','Ref', 'Alt', 'AAref', 'AAalt',
       'AAPOS', 'CodonCoverage', 'VAF(DP4)','Mutation','Drug Resistance Marker','QD', 
       'SOR', 'MQ','MQRankSum', 'Filter', 'FilterDescription','region','Treatment','Day']]
summary.to_csv("summary.csv")
merge1.to_csv("depth_filter_input.csv")

list0=pd.read_csv(args.cutoff_dir+"list.csv",sep="\t",header=None)
list1=pd.read_csv(args.cutoff_dir+"list1.csv",sep="\t",header=None)
list0.columns=['file name','Gene']
list0['Sample name']=list0['file name'].apply(lambda x:x.split("_final")[0])

merge1_wild=merge1[merge1.Mutation=="Wildtype_N/A"]

list1.columns=['file name','Gene','BasePOS','Depth']
list1['Sample name']=list1['file name'].apply(lambda x:x.split("_depth")[0])

test2=pd.merge(merge1,list1,on=['Sample name','Gene','BasePOS'])
test2=test2.drop_duplicates()
test3=pd.concat([merge1_wild,test2])
test3=test3.drop_duplicates()
test3.to_csv("cutoff_SNP.csv")

all_df1=test3[test3.Reportable=="Yes"]
merge1_AAonly=test3[['Gene','AAref', 'AAalt', 'AAPOS',
       'CodonCoverage', 'AF', 'Mutation', 'QD', 'SOR', 'MQ', 'MQRankSum',
       'Filter', 'FilterDescription', 'Candidates', 'Reportable',
       'Sample name','VAF(DP4)','Reportable mutation', 'Drug Resistance Marker', 'region','Treatment'
       ]]
merge1_AAonly.columns=['Gene','AAref', 'AAalt', 'AAPOS',
       'CodonCoverage', 'AF', 'Mutation', 'QD', 'SOR', 'MQ', 'MQRankSum',
       'Filter', 'FilterDescription', 'Candidates', 'Reportable',
       'Sample name','VAF','Reportable mutation', 'Drug Resistance Marker','region','Treatment'
       ]
merge1_AAonly['VAF']=merge1_AAonly['VAF'].astype("float",errors='ignore')
merge1_AAonly['Filter']=merge1_AAonly['Filter'].astype("float",errors='ignore')
merge1_AAonly=merge1_AAonly.drop_duplicates()
#merge1_AAonly.to_csv("BI_input.csv")

haplo=merge1_AAonly[['Sample name','Gene','AAalt','AAPOS','Reportable mutation','VAF']]
pd.to_numeric(haplo.VAF,errors='ignore')
pd.to_numeric(haplo.VAF,errors='coerce')

haplo=haplo[haplo.VAF >= 0.95]
haplo.AAPOS=haplo.AAPOS.astype('int')
haplo=haplo.drop_duplicates()
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
haplo_table.to_csv("haplotype.csv")

haplotype_count=pd.DataFrame(haplo_table.haplotype.value_counts()).reset_index().drop(0,axis=0)
haplotype_count.columns=['haplotype','count']
h2=pd.merge(haplo_table,haplotype_count,on="haplotype",how="left")
h2['count']=h2['count'].fillna(0)
h3=h2.pivot_table(columns=['haplotype'],values='count',index=['Sample name','Gene'],aggfunc=lambda x:x).reset_index()
h3_list=["CVMNK","CVIEK","SVMNT","CVIET","NYSND","NFSND","YFSND","YYSNY"]
for col in h3_list:
    if col not in h3.columns:
        h3[col] = np.nan
h3=h3.fillna(0)
h3.head()
h3['Year']=h3['Sample name'].apply(lambda x:x[0:2])
h3['region']=h3['Sample name'].apply(lambda x:x[2:6])
h3['region']=h3.region.apply(lambda x:dic2[x])
h3['Treatment']=h3['Sample name'].apply(lambda x:dic3[x[8]])
h3.to_csv("haplotype_count.csv")



merge1_sub=merge1[['Gene', 'AAref', 'AAalt', 'AAPOS', 
       'Candidates', 'Reportable',
       'Sample name', 'VAF(DP4)',
       'Reportable mutation','region','Treatment'
       ]]
merge1_sub.columns=['Gene', 'AAref', 'AAalt', 'AAPOS',
       'Candidates', 'Reportable',
       'Sample name', 'VAF',
        'Reportable mutation','region','Treatment'
       ]
merge1_sub=merge1_sub.drop_duplicates()
b=merge1_sub.groupby(['Gene','Reportable mutation','Treatment']).size()
b.columns=['Gene','Reportable mutation','Treatment','Count']
b.to_csv("summary_reportable.csv")

all_df_sub=merge1[merge1.Mutation!="Wildtype"]
all_df_NCBI=all_df_sub[['Sample name','Gene','Reportable mutation']]
all_df_NCBI=all_df_NCBI.pivot_table(index='Sample name',columns="Gene",values='Reportable mutation',aggfunc=lambda x:', '.join(x.astype(str).unique()))
all_df_NCBI.to_csv("NCBI_feature_column.csv")


merge1_AAonly=merge1_AAonly.reset_index()
test3=test3.reset_index()
test3['combine_key']=''
summary['combine_key']=''
for i in range(0,len(test3)):
    test3['combine_key'][i]=str(test3['AAPOS'][i])+test3['AAalt'][i]+test3['Sample name'][i]
for i in range(0,len(summary)):
    summary['combine_key'][i]=str(summary['AAPOS'][i])+summary['AAalt'][i]+summary['Sample name'][i]
a=pd.merge(test3,summary,on="combine_key")[['combine_key']].drop_duplicates()
bi=a.combine_key.unique()
summary['flagged']=''
for i in range(0,len(summary)):
    if summary.Mutation[i] != "Wildtype":
       if summary.combine_key[i] not in bi:
          summary['flagged'][i]="Q1"
summary=summary.drop_duplicates()
bi_out=summary[['Sample name','Gene', 'AAref', 'AAalt', 'AAPOS', 'VAF(DP4)', 'Mutation',
        'Drug Resistance Marker','QD', 'SOR', 'MQ','MQRankSum', 'Filter', 'FilterDescription',
        'flagged','region','Treatment']]
bi_out.to_csv("PowerBI_input.csv")