
import pandas as pd
import os
import subprocess
import glob
import argparse
import matplotlib.pyplot as plt


def install(name):
    subprocess.call(['pip', 'install', name])

install('seaborn')
install('dexplot')

import dexplot as dxp
import seaborn as sns
parser = argparse.ArgumentParser(description='name')
parser.add_argument('-f', dest='snpfilter', type=str, help="snpfilter output")
parser.add_argument('-voi', dest='voifile',type=str,help="voi files")
parser.add_argument('-d',dest='directory',type=str, help="snpfilter output")
args = parser.parse_args()


all_df=pd.DataFrame()
dic={'ANZa':(-6.2667,14.2333),'ANLS':(-10.2592,20.6249),'ANBe':(-12.5763,13.4055),'GNMa':(9.546,-13.281),'GNHa':(11.56667,-10.05),'GNDo':(7.75833,-8.80472),'GNLS':(11.03856,-12.04793)}
dic2={'ANBe':'Benguela','ANZa':'Zaire','ANLS':'Lunda Sul','GNLS':'Labé','GNHa':'Dabola','GNDo':'Nzérékoré','GNMa':'Forécariah Centre','ANxx':'control'}
dic3={'A':'AL','B':'CQ','C':'CQ+PQ','D':'AL+PQ','E':'AL+MQ','F':'AS+AQ','G':'AS+SP','H':'AS+MQ','I':'DHA+PQ','J':'AS+SP+PQ','K':'AL+MQ+PQ','L':'Doxxycline','M':'Malarone','N':'Quinine','O':'Quinine+Doxycycline','P':'Atovaquone/Proguanil+AL','Q':'Mefloquine','R':'Vancomycin+Rocephin','S':'Malarone+Doxycycline','x':'None'}

for a in glob.glob("*.csv"):
    if not a.startswith("19US"):
#for a in glob.glob(args.directory+"*.csv"):
#for filename in os.listdir(args.directory+".csv"):
    #open(os.path.join(ar
    # gs.directory, filename)
    #with open(os.path.join(args.directory, filename), "r") as a:

        name=a.split(".csv")[0]
        a=pd.read_csv(a)
        a['Sample name']=name
        a['region']=dic2[name[2:6]]
        a['Year']='20'+name[0:2]
        a['Treatment']=dic3[name[8:9]]
    #a['pool_sample']=name[10:14]
        a['Day']=name[6:8]
        a['VAF']=a.Mutation.apply(lambda x:x.split("_")[1])
        a['Mutation_Allele']=a.Mutation.apply(lambda x:x.split("_")[0])
    # for key in dic:
    #    a['Latitude']=dic[name[2:6]][0]
    #    a['Longitude']=dic[name[2:6]][1]
        all_df=all_df.append(a)

# all_df=pd.DataFrame()
# with open(args.snpfilter, "r") as f1:
#     f1=pd.read_csv(f1)
#     all_df=all_df.append(f1)
# all_df.AAPOS=all_df.AAPOS.astype(str)
#all_df['Reportable mutation with gene']=all_df['Gene']+"_"+all_df['AAref']+all_df['AAPOS']+all_df['AAalt']
#all_df['Reportable mutation']=all_df['AAref']+all_df['AAPOS']+all_df['AAalt']
all_df.Reportable=all_df.Reportable.str.replace('no','No')
all_df.Reportable=all_df.Reportable.str.replace('reportable','Yes')
all_df.Candidates=all_df.Candidates.str.replace('no','No')
all_df.Candidates=all_df.Candidates.str.replace('candidates','Yes')
all_df['first']=all_df['AAref']+all_df['AAPOS'].astype('str')
#all_df=all_df[all_df.CodonCoverage >1 ]

voi=pd.read_csv(args.voifile)
voi['first']=voi["RefAA"]+voi["AAPos"].astype('str')
merge1=pd.merge(all_df,voi,on='first')
merge1.columns=['Gene','BasePOS', 'BaseDepth', 'Ref', 'Alt', 'AAref', 'AAalt',
       'AAPOS', 'CodonCoverage', 'AF', 'Mutation', 'QD', 'SOR', 'MQ',
       'MQRankSum', 'Filter', 'FilterDescription', 'Candidates', 'Reportable',
       'Sample name', 'region', 'Year', 'Treatment', 'Day', 'VAF',
       'Mutation_Allele', 'first', 'Chr', 'Gene_y',
       'RefAA', 'AAPos', 'AltAA']
merge1['Reportable mutation']=merge1['RefAA']+merge1['AAPos'].astype("str")+merge1['AltAA']
merge1['Drug Resistance Marker']=merge1['AAPos'].astype("str")+merge1['RefAA']+merge1['AltAA']
merge1['Reportable mutation with gene']=merge1['Gene']+"_"+merge1['RefAA']+merge1['AAPos'].astype("str")+merge1['AltAA']
merge1.to_csv("depth_filter_input.csv")

all_df1=merge1[merge1.Reportable=="Yes"]

merge1_AAonly=merge1[['Gene','AAref', 'AAalt', 'AAPOS',
       'CodonCoverage', 'AF', 'Mutation', 'QD', 'SOR', 'MQ', 'MQRankSum',
       'Filter', 'FilterDescription', 'Candidates', 'Reportable',
       'Sample name', 'region', 'Year', 'Treatment', 'Day', 'VAF',
       'Mutation_Allele', 'first','Reportable mutation', 'Drug Resistance Marker',
       'Reportable mutation with gene']]
merge1_AAonly['VAF']=merge1_AAonly['VAF'].astype("float",errors='ignore')
merge1_AAonly['Filter']=merge1_AAonly['Filter'].astype("float",errors='ignore')
merge1_AAonly=merge1_AAonly.drop_duplicates()
merge1_AAonly.to_csv("BI_input.csv")

# haplo=merge1_AAonly[['Sample name','Gene','AAalt','AAPOS','Reportable mutation']]
# haplo.AAPOS=haplo.AAPOS.astype('int')
# haplo=haplo.drop_duplicates()
# b=haplo.pivot(columns=['Reportable mutation'],values='AAalt',index=['Sample name','Gene'])
# b=b.sort_index()
# b.columns=b.columns.sort_values()


merge1_sub=merge1[['Gene', 'AAref', 'AAalt', 'AAPOS',
       'Candidates', 'Reportable',
       'Sample name', 'region', 'Year', 'Treatment', 'Day', 'VAF',
       'Mutation_Allele', 'first', 'Chr', 'Gene_y',
       'RefAA', 'AAPos', 'AltAA', 'Reportable mutation',
       'Reportable mutation with gene',]]
merge1_sub=merge1_sub.drop_duplicates()
# #merge1_sub.to_csv("BI_input.csv")

# for ax_num, score in zip(range(1,5), ['f1', 'recall', 'accuracy', 'precision']):
#     plt.subplot(2,2,ax_num)
#     sns.barplot(x='model_name', y='score', hue='train_val_test',
#                 data=classification_scores[classification_scores['score_name'] == score])
#     plt.xticks(rotation=15, fontsize=14)
#g=sns.catplot(y="Reportable mutation",data=merge1_sub,hue="Mutation_Allele",col="Gene",kind="count",sharey=False,dodge = False,col_wrap=2)
#g.savefig('summaryplot_reportable.png')

#g2=sns.catplot(x="Reportable mutation",y="CodonCoverage",data=all_df1,col="Gene",kind='violin',sharex=False,dodge = False,col_wrap=2)
#sns.catplot(x="Reportable mutation",y="CodonCoverage",data=all_df1,col="Gene",sharex=False,dodge = False)
# for ax in g2.axes:
#     plt.setp(ax.get_xticklabels(), visible=True, rotation=45)
# plt.gcf().autofmt_xdate()
# g2.savefig('coverageplot_reportable.png')

# sns.set(font_scale=0.5)
# g3=dxp.count(val='Reportable mutation with gene',data=merge1_sub,split='Mutation_Allele',normalize=['Reportable mutation with gene'],stacked=True,sharey=False,sharex=True,orientation='h')
# g3.legend(loc='upper left', bbox_to_anchor=(0.6,0.98))
# g3.savefig('normalization_reportable.png')



b=merge1_sub.groupby(['region','Treatment','Day','Gene','Reportable mutation','Mutation_Allele']).size()
#b=b.drop_duplicates()
# b=pd.DataFrame(b)
# b.columns=['count']
b.to_csv("summary_reportable.csv")

all_df_sub=merge1[merge1.Mutation_Allele!="Wildtype"]
all_df_NCBI=all_df_sub[['Sample name','Gene','Reportable mutation']]
all_df_NCBI=all_df_NCBI.pivot_table(index='Sample name',columns="Gene",values='Reportable mutation',aggfunc=lambda x:', '.join(x.astype(str).unique()))
all_df_NCBI.to_csv("NCBI_feature_column.csv")
##if len(all_df2)!= 0:    
# g=sns.catplot(y="Reportable mutation",data=all_df2,hue="Mutation",col="Gene",kind="count",sharey=False,dodge = False)
# g.savefig('summaryplot_candidates.png')
# g2=sns.catplot(x="Reportable mutation",y="CodonCoverage",data=all_df2,col="Gene",kind='violin',sharex=False,dodge = False)
# g2.savefig('coverageplot_candidates.png')
# sns.set(font_scale=0.5)
# g3=dxp.count(val='Reportable mutation',data=all_df2,split='Mutation',normalize=['Reportable mutation with gene'],stacked=True,sharey=False,sharex=True,orientation='h')
# g3.savefig('normalization_candidates.png')

# b2=all_df1.groupby(['Gene','Reportable mutation','Mutation']).size()
# b2=pd.DataFrame(b2)
# b2.columns=['count']
# b2.to_csv("summary_candidates.csv")
