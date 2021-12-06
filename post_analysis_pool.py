import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

parser = argparse.ArgumentParser(description='name')
parser.add_argument('-indi', dest='indi_powerbi', type=str, help="individual PowerBI input")
parser.add_argument('-pool', dest='pool_powerbi',type=str,help="pool PowerBI input")
parser.add_argument('-rr',dest='recru_reinfect',type=str, help="molecular classification file")
parser.add_argument('-poolinfo', dest='poolinfo', type=str, help="pool information file")
parser.add_argument('-f', dest='snpfilter', type=str, help="snpfilter output")
args = parser.parse_args()

indi=pd.read_csv(args.indi_powerbi)
pool=pd.read_csv(args.pool_powerbi)
poolname=pd.read_csv(args.poolinfo)
RR=pd.read_csv(args.recru_reinfect)
all_df=pd.DataFrame()

# for a in glob.glob(snpfilter"*.csv"):
#        if not a.startswith("list"):
#            if not a.startswith("21US"):
#                 name=a.split(".csv")[0]
#                 a=pd.read_csv(a)
#                 a['Sample name']=name
#                 all_df=all_df.append(a)

# agreements=all_df[['Sample name','Agreents','Gene','AAPOS']]    
# agreements.columns=['Sample name','Agreement','Gene','AAPOS']
# pool=pd.read_csv("./output/visualization/PowerBI_input.csv")
# # pool master file from lab
# indi=pd.read_csv("../inputfiles/output2/snpfilter/PowerBI_input.csv")
# poolname=pd.read_csv("/Users/subinpark/Nf-NeST/gn2_pool_master.csv")



dic2={'ANBe':'Benguela','ANZa':'Zaire','ANLS':'Lunda Sul','GNLa':'Labé','GNHa':'Dabola','GNDo':'Nzérékoré','GNMa':'Forécariah Centre','ANxx':'control','USxx':'control'}


pool['VAF(DP4)']=pool['VAF(DP4)'].fillna(0)
pool['Year']=pool['Sample name'].apply(lambda x:x[0:2]).astype(int)
pool['poolname']=pool['Sample name'].apply(lambda x:x[9:13])
pool['Study site']=''
for i in range(0,len(pool)):
    pool['Study site'][i]=dic2[pool['Sample name'][i][2:6]]
poolname['study_site_pool']=''
for i in range(0,len(poolname)):
    poolname['study_site_pool'][i]=poolname['Study site'][i]+"_"+poolname['Final ID'][i]
pool['study_site_pool']=''
for i in range(0,len(pool)):
    pool['study_site_pool'][i]=pool['Study site'][i]+"_"+pool['poolname'][i]
poolname.CSID=poolname.CSID.astype(str)
poolname['poolsize']=poolname['study_site_pool'].map(poolname['study_site_pool'].value_counts())

a=pd.merge(pool,poolname,left_on="study_site_pool",right_on="study_site_pool")
a['Year']=a['Sample name'].apply(lambda x:x[0:2])

a=a[['Sample name', 'Gene', 'AAref', 'AAalt', 'AAPOS', 'VAF(DP4)',
       'Mutation', 'Drug Resistance Marker', 'QD', 'SOR', 'MQ', 'MQRankSum',
       'Filter', 'FilterDescription', 'flagged', 'poolname', 'poolsize', 'CSID',
       'Country Sample ID','Collection date', 'country', 'Day', 'Treatment/Drug arm',
       'Treatment Code ', 'Study site_x', 'Microcopy (species)',
       'type', 'markers', 'repeat', 'AMD Sample ID (concatenate) ','study_site_pool','Year']]
a.columns=['Sample name', 'Gene', 'AAref', 'AAalt', 'AAPOS', 'VAF(DP4)',
       'Mutation', 'Drug Resistance Marker', 'QD', 'SOR', 'MQ', 'MQRankSum',
       'Filter', 'FilterDescription', 'flagged', 'poolname', 'poolsize', 'CSID',
       'Country Sample ID','Collection date', 'country', 'Day', 'Treatment/Drug arm',
       'Treatment Code ', 'Study site', 'Microcopy (species)',
       'type', 'markers', 'repeat', 'AMD Sample ID (concatenate) ','study_site_pool','Year']
a=a.drop_duplicates()
a.to_csv("pool_PowerBI_input.csv")

pool_counts=a[['study_site_pool','poolsize']].drop_duplicates().reset_index()
pool_counts=pool_counts.drop('index',axis=1)
# pool_counts.to_csv("pool_counts.csv")
pool_weight=pd.DataFrame(a.groupby(['Study site','Gene','Drug Resistance Marker','VAF(DP4)','AAPOS','study_site_pool','Year']).size()).reset_index()
pool_weight.columns=['Study site', 'Gene', 'Drug Resistance Marker', 'VAF(DP4)', 'AAPOS','poolname','Year','count_pool']
pool_weight["pool_weight_sum"]=pool_weight['VAF(DP4)']*pool_weight.count_pool
pool_weight=pool_weight[['Study site','Gene','Drug Resistance Marker','pool_weight_sum','count_pool','Year']]
pool_weight_sum=pd.DataFrame(pool_weight.groupby(['Study site','Gene','Drug Resistance Marker','Year'])[['pool_weight_sum','count_pool']].agg('sum')).reset_index()

indi=indi.drop("Unnamed: 0",axis=1)
indi['day']=indi['Sample name'].apply(lambda x:x[6:8])
indi['day']=indi['day'].astype(int)
indi=indi[indi['day']==0]
indi['Year']=indi['Sample name'].apply(lambda x:x[0:2]).astype(int)
indi['VAF(DP4)']=indi['VAF(DP4)'].fillna(0)
indi_weight=pd.DataFrame(indi.groupby(['region','Gene','Drug Resistance Marker','VAF(DP4)','AAPOS','Year']).size()).reset_index()
indi_weight.columns=['region', 'Gene', 'Drug Resistance Marker', 'VAF(DP4)', 'AAPOS','Year','count_indi']
indi_weight["indi_weight_sum"]=indi_weight['VAF(DP4)']*indi_weight.count_indi
indi_weight_sum=pd.DataFrame(indi_weight.groupby(['region','Gene','Drug Resistance Marker','Year'])['indi_weight_sum','count_indi'].sum().reset_index())

pool_weight_sum.columns=['region','Gene','Drug Resistance Marker','Year','pool_weight_sum',
       'count_pool']
pool_weight_sum['Year']=pool_weight_sum['Year'].astype(str)
indi_weight_sum['Year']=indi_weight_sum['Year'].astype(str)
final=pd.merge(pool_weight_sum,indi_weight_sum,on=['region','Gene','Drug Resistance Marker','Year'])
final['final_weight']=(final.pool_weight_sum+final.indi_weight_sum)/(final.count_pool+final.count_indi)

final['type']=''
for i in range(0,len(final)):
    if final.final_weight[i]==0:
        final.type[i]="Wildtype"
    else:
        final.type[i]="Mutant"

final['key']=final.region+final.Gene+final['Drug Resistance Marker']+final['Year']
final['total']=final['count_indi']+final['count_pool']
final.to_csv("weighted.csv")
pool['VAF(DP4)']=pool['VAF(DP4)'].fillna(0)
indi['VAF(DP4)']=indi['VAF(DP4)'].fillna(0)

all_=pd.concat([a,indi])
indi['Country Sample ID']=indi['Sample name']
all_=all_.reset_index()
all_['type2']=''
for i in range(0,len(all_)):
    if all_['VAF(DP4)'][i]==0:
        all_['type2'][i]="Wildtype"
    else:
        all_['type2'][i]="Mutant"
all_['poolsize']=all_['poolsize'].fillna(1)
all_['poolsize']=all_['poolsize'].astype(int)
all_['Year']=all_['Sample name'].apply(lambda x:x[0:2])

#all_.to_csv("pool_indi.csv")

count_table=pd.DataFrame(all_.groupby(['Gene','Drug Resistance Marker','Year','Study site'])['poolsize'].count())
count_table=count_table.reset_index()
count_table.columns=['Gene', 'Drug Resistance Marker', 'Year', 'Study site', 'total_count']
all2_=pd.merge(count_table,all_,on=['Gene','Drug Resistance Marker','Year','Study site'])
all2_['key']=all2_['Study site']+all2_.Gene+all2_['Drug Resistance Marker']+all2_['Year']
# all2_.columns=['Gene', 'Drug Resistance Marker', 'Year', 'region', 'total_count',
#        'index', 'Sample name', 'AAref', 'AAalt', 'AAPOS', 'VAF(DP4)',
#        'Mutation', 'QD', 'SOR', 'MQ', 'MQRankSum', 'Filter',
#        'FilterDescription', 'flagged', 'poolname', 'poolsize', 'CSID',
#        'Country Sample ID', 'Collection date', 'country', 'Day',
#        'Treatment/Drug arm', 'Treatment Code ', 'Microcopy (species)', 'type_sex',
#        'markers', 'repeat', 'AMD Sample ID (concatenate) ', 'study_site_pool',
#        'AA', 'Treatment', 'VAF_final', 'day', 'type', 'key']
#all3=pd.merge(all2_,agreements,on=['Sample name','Agreement','Gene','AAPOS'])
all2_['newname']=''
for i in range(0,len(all2_)):
    all2_['newname'][i]=all2_.Gene[i]+":"+all2_['Drug Resistance Marker'][i]+"(N="+str(all2_['total_count'][i])+")"
all2_['site_year']=all2_['Study site']+all2_['Year'].astype(str)
all2_=all2_[(all2_.Gene == "K13") |(all2_.Gene == "PfCRT") | (all2_.Gene == "PfMDR1")]
all2_.to_csv("pool_indi.csv")

for a in all2_.site_year.unique():
    
    df_sub=all2_[all2_.site_year==a]

#g = sns.FacetGrid(tips, col = 'size',  row = 'smoker', hue = 'day')
    plt.figure(figsize=(80, 32))
    total = df_sub.groupby(['newname','AAPOS','Gene','site_year'])['poolsize'].sum().reset_index()
    wildtype = df_sub[df_sub.type2=="Wildtype"].groupby(['newname','AAPOS','Gene','site_year'])['poolsize'].sum().reset_index()
    wildtype['poolsize'] = [i / j * 100 for i,j in zip(wildtype['poolsize'], total['poolsize'])]
    total['poolsize']=[i / j * 100 for i,j in zip(total['poolsize'],total['poolsize'])]
    wildtype=wildtype.sort_values(by=['Gene','AAPOS'],ascending=True)
    total=total.sort_values(by=['Gene','AAPOS'],ascending=True)
#g=sns.FacetGrid(df, col = 'poolsize',  row = 'newname', hue = 'site_year')
# bar chart 1 -> top bars (group of 'smoker=No')
    bar1 = sns.barplot(x="poolsize",  y='newname', data=total, color='darkblue')
# bar chart 2 -> bottom bars (group of 'smoker=Yes')
    bar2 = sns.barplot(x="poolsize", y='newname', data=wildtype, color='lightblue')
# add legend
    top_bar = mpatches.Patch(color='darkblue', label='Mutant')
    bottom_bar = mpatches.Patch(color='lightblue', label='Wildtype')
    plt.legend(handles=[top_bar, bottom_bar],loc='upper right',fontsize=29)
    plt.xlabel("SNP allele frequency")
    plt.ylabel("")
    
    size=50
    params = {'legend.fontsize': 'large',
          'figure.figsize': (20,8),
          'axes.labelsize': size,
          'axes.titlesize': size,
          'xtick.labelsize': size*0.75,
          'ytick.labelsize': size*0.75,
          'axes.titlepad': 25}
    plt.rcParams.update(params)
    plt.tight_layout()
    
    plt.savefig(a+".pdf")
    

test=pd.merge(all2_,final,on=['Gene','Drug Resistance Marker','Year','region','type','key'],how='left').drop_duplicates()
# test1=test[(test.Gene == "K13") |(test.Gene == "PfCRT") | (test.Gene == "PfMDR1")]
# count1=test1.groupby(["Year","region"])['total_count'].max().reset_index()

indi=indi.reset_index()
indi['AMD ID']=''
for i in range(0,len(indi)):
    indi['AMD ID'][i]=indi['Sample name'][i].split("_")[0]
#all2_['AMD ID']=all2_['Sample name'].str.split("_")[0]
rr_all=pd.merge(indi,RR,on="AMD ID",how="left").drop_duplicates()
rr_all['VAF(DP4)']=rr_all['VAF(DP4)'].fillna(0)
rr_all['day']=rr_all['Sample name'].apply(lambda x:x[6:8])
rr_all['Year']=rr_all['Sample name'].apply(lambda x:x[0:2]).astype(int)
rr_all.to_csv("molecular_classification.csv")

