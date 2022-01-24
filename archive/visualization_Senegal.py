
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

for a in glob.glob("*.csv"):
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

merge1.columns=['Gene','BasePOS', 'BaseDepth', 'Ref', 'Alt', 'AAref', 'AAalt',
       'AAPOS', 'CodonCoverage', 'VAF(DP4)','VAF(AD/DP)','AF','Mutation', 'QD', 'SOR', 'MQ',
       'MQRankSum', 'Filter', 'FilterDescription', 'Candidates', 'Reportable',
       'Sample name', 'Reportable mutation with gene','Reportable mutation', 'first','Chr','Gene_y', 
       'RefAA', 'AAPos', 'AltAA']
merge1['AAPOS']=merge1['AAPOS'].astype(str)
#merge1['Reportable mutation']=merge1['AAref']+merge1['AAPOS']+merge1['AAalt']
merge1['Reportable mutation']=''
for i in range(0,len(merge1)):
    if merge1.Mutation[i]=="Wildtype":
        merge1['Reportable mutation'][i]=merge1['AAPOS'][i]+merge1['AAref'][i]
    else:
        merge1['Reportable mutation'][i]=merge1['AAPOS'][i]+merge1['AAref'][i]+">"+merge1['AAalt'][i]
merge1['Drug Resistance Marker']=merge1['AAPos'].astype(str)+merge1['RefAA'].astype(str)+">"+merge1['AltAA'].astype(str)
merge1['Reportable mutation with gene']=merge1['Gene']+"_"+merge1['AAref']+merge1['AAPOS'].astype("str")+merge1['AAalt']
merge1=merge1[['Gene','BasePOS', 'BaseDepth', 'Ref', 'Alt', 'AAref', 'AAalt',
       'AAPOS', 'CodonCoverage', 'VAF(DP4)','VAF(AD/DP)','AF','Mutation', 'QD', 'SOR', 'MQ',
       'MQRankSum', 'Filter', 'FilterDescription', 'Candidates', 'Reportable',
       'Sample name', 'Reportable mutation with gene','Reportable mutation','Drug Resistance Marker']]

summary=merge1[['Sample name','Gene','Ref', 'Alt', 'AAref', 'AAalt',
       'AAPOS', 'CodonCoverage', 'VAF(DP4)','Mutation','Drug Resistance Marker','QD', 
       'SOR', 'MQ','MQRankSum', 'Filter', 'FilterDescription']]
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
       'Sample name','VAF(DP4)','Reportable mutation', 'Drug Resistance Marker',
       ]]
merge1_AAonly.columns=['Gene','AAref', 'AAalt', 'AAPOS',
       'CodonCoverage', 'AF', 'Mutation', 'QD', 'SOR', 'MQ', 'MQRankSum',
       'Filter', 'FilterDescription', 'Candidates', 'Reportable',
       'Sample name','VAF','Reportable mutation', 'Drug Resistance Marker',
       ]
merge1_AAonly['VAF']=merge1_AAonly['VAF'].astype("float",errors='ignore')
merge1_AAonly['Filter']=merge1_AAonly['Filter'].astype("float",errors='ignore')
merge1_AAonly=merge1_AAonly.drop_duplicates()
#merge1_AAonly.to_csv("BI_input.csv")

haplo=merge1_AAonly[['Sample name','Gene','AAalt','AAPOS','Reportable mutation']]
haplo.AAPOS=haplo.AAPOS.astype('int')
haplo=haplo.drop_duplicates()
b=haplo.pivot_table(columns=['Reportable mutation'],values='AAalt',index=['Sample name','Gene'],aggfunc='first')
b=b.sort_index()
b.columns=b.columns.sort_values()

merge1_sub=merge1[['Gene', 'AAref', 'AAalt', 'AAPOS', 
       'Candidates', 'Reportable',
       'Sample name', 'VAF(DP4)',
       'Reportable mutation',
       ]]
merge1_sub.columns=['Gene', 'AAref', 'AAalt', 'AAPOS',
       'Candidates', 'Reportable',
       'Sample name', 'VAF',
        'Reportable mutation',
       ]
merge1_sub=merge1_sub.drop_duplicates()
b=merge1_sub.groupby(['Gene','Reportable mutation']).size()
b.columns=['Gene','Reportable mutation','Count']
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
    if summary.combine_key[i] not in bi:
        summary['flagged'][i]="Q1"
summary=summary.drop_duplicates()
bi_out=summary[['Sample name','Gene', 'AAref', 'AAalt', 'AAPOS', 'VAF(DP4)',
        'Drug Resistance Marker','QD', 'SOR', 'MQ','MQRankSum', 'Filter', 'FilterDescription',
        'flagged']]
bi_out.to_csv("PowerBI_input.csv")