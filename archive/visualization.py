
import pandas as pd
import os
import subprocess
import glob
import argparse

def install(name):
    subprocess.call(['pip', 'install', name])

install('seaborn')
install('dexplot')

import dexplot as dxp
import seaborn as sns
parser = argparse.ArgumentParser(description='name')
parser.add_argument('-f', dest='snpfilter', type=str, help="snpfilter output")
args = parser.parse_args()

all_df=pd.DataFrame()
with open(args.snpfilter, "r") as f1:
    f1=pd.read_csv(f1)
    all_df=all_df.append(f1)

all_df.AAPOS=all_df.AAPOS.astype(str)
all_df['Reportable mutation']=all_df['AAref']+all_df['AAPOS']+all_df['AAalt']
all_df['Reportable mutation with gene']=all_df['Gene']+"_"+all_df['AAref']+all_df['AAPOS']+all_df['AAalt']
all_df.Reportable=all_df.Reportable.str.replace('no','No')
all_df.Reportable=all_df.Reportable.str.replace('reportable','Yes')
all_df.Candidates=all_df.Candidates.str.replace('no','No')
all_df1=all_df[all_df.Reportable=="Yes"]
all_df1=all_df1[all_df1.CodonCoverage >1 ]
all_df2=all_df[all_df.Candidates=="Yes"]


g=sns.catplot(y="Reportable mutation",data=all_df1,hue="Mutation",col="Gene",kind="count",sharey=False,dodge = False)
g.savefig('summaryplot_reportable.png')




g2=sns.catplot(x="Reportable mutation",y="CodonCoverage",data=all_df1,col="Gene",kind='violin',sharex=False,dodge = False)
g2.savefig('coverageplot_reportable.png')


sns.set(font_scale=0.5)
g3=dxp.count(val='Reportable mutation with gene',data=all_df1,split='Mutation',normalize=['Reportable mutation with gene'],stacked=True,sharey=False,sharex=True,orientation='h')
g3.savefig('normalization_reportable.png')



b=all_df1.groupby(['Gene','Reportable mutation','Mutation']).size()
b=pd.DataFrame(b)
b.columns=['count']
b.to_csv("summary_reportable.csv")

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
