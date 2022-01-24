
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
parser.add_argument('-cod', dest='cutoff_dir', type=str, help="cutoff files directory")
args = parser.parse_args()


all_df=pd.DataFrame()
dic={'ANZa':(-6.2667,14.2333),'ANLS':(-10.2592,20.6249),'ANBe':(-12.5763,13.4055),'GNMa':(9.546,-13.281),'GNHa':(11.56667,-10.05),'GNDo':(7.75833,-8.80472),'GNLS':(11.03856,-12.04793)}
dic2={'ANBe':'Benguela','ANZa':'Zaire','ANLS':'Lunda Sul','GNLS':'Labé','GNHa':'Dabola','GNDo':'Nzérékoré','GNMa':'Forécariah Centre','ANxx':'control'}
dic3={'A':'AL','B':'CQ','C':'CQ+PQ','D':'AL+PQ','E':'AL+MQ','F':'AS+AQ','G':'AS+SP','H':'AS+MQ','I':'DHA+PQ','J':'AS+SP+PQ','K':'AL+MQ+PQ','L':'Doxxycline','M':'Malarone','N':'Quinine','O':'Quinine+Doxycycline','P':'Atovaquone/Proguanil+AL','Q':'Mefloquine','R':'Vancomycin+Rocephin','S':'Malarone+Doxycycline','x':'None'}

for a in glob.glob("19*.csv"):
    if not a.startswith("19US"):
        name=a.split(".csv")[0]
        a=pd.read_csv(a)
        a['Sample name']=name
        a['region']=dic2[name[2:6]]
        a['Year']='20'+name[0:2]
        a['Treatment']=dic3[name[8:9]]
        a['Day']=name[6:8]
        a['VAF']=a.Mutation.apply(lambda x:x.split("_")[1])
        a['Mutation_Allele']=a.Mutation.apply(lambda x:x.split("_")[0])

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


#os.chdir(".")
# list1=pd.DataFrame()
# cutcov=pd.read_csv(args.cutoff_dir+"first_q_coverage.csv")
# for file1 in glob.glob('*/*coverage.txt'):
#     f=pd.read_csv(file1,sep="\t")
#     for i in range(0,len(f)):
#         for j in range(0,len(cutcov)):
#             if (f['#rname'][i] == cutcov['#rname'][j]) & (f['coverage'][i] >= cutcov['25%'][j]):
#                 list1['Sample name']=file1.split("/")[0]
#                 list1['Gene']=f['#rname'][i]
#                 list1=list1.append(list1)
#                 # with open(os.getcwd()+"list.csv",'a') as ofile:
#                 #     ofile.write(file1.split("/")[0]+"\t"+f['#rname'][i]+"\n")

# #depth=pd.read_csv("../cutoff/first_q_depth.csv")
# depth=pd.read_csv(args.cutoff_dir+"first_q_depth.csv")
# for file1 in glob.glob('*/*depth.txt'):
#     if os.stat(file1).st_size > 0:
#         f=pd.read_csv(file1,sep="\t",header=None)
#         f.columns=['#rname','position','depth']
#         for i in range(0,len(f)):
#             for j in range(0,len(depth)):
#             #print(f.head())
#                 if (f['#rname'][i] == cutcov['#rname'][j]) & (f['depth'][i] >= depth['25%'][j]):
#                     with open(os.getcwd()+"list1.csv",'a') as ofile1:
#                         ofile1.write(file1.split("/")[0]+"\t"+f['#rname'][i]+"\t"+str(f['position'][i])+"\t"+str(f['depth'][i])+"\n")
                    # ofile1.close

#input_=pd.read_csv("../../vis11//visualization/depth_filter_input.csv")

list0=pd.read_csv(args.cutoff_dir+"list.csv",sep="\t",header=None)
list1=pd.read_csv(args.cutoff_dir+"list1.csv",sep="\t",header=None)
list0.columns=['file name','Gene']
list0['Sample name']=list0['file name'].apply(lambda x:x.split("_final")[0])
#list0.columns=['Sample name','Gene']

#test1=pd.merge(merge1,list0,on=['Gene','Sample name'],how='inner')
#test1=test1.drop_duplcicates()
#test1_wild=test1[test1.Mutation=="Wildtype_N/A"]
merge1_wild=merge1[merge1.Mutation=="Wildtype_N/A"]

list1.columns=['file name','Gene','BasePOS','Depth']
list1['Sample name']=list1['file name'].apply(lambda x:x.split("_depth")[0])

test2=pd.merge(merge1,list1,on=['Sample name','Gene','BasePOS'])
test2=test2.drop_duplicates()
#test3=pd.concat([test1_wild,test2])
test3=pd.concat([merge1_wild,test2])
test3=test3.drop_duplicates()
test3.to_csv("cutoff_SNP.csv")


all_df1=test3[test3.Reportable=="Yes"]
merge1_AAonly=test3[['Gene','AAref', 'AAalt', 'AAPOS',
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
