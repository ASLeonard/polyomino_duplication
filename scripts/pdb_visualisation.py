import matplotlib.pyplot as plt
import seaborn as sns
import glob
import pandas as pd

class Dimer(object):
     __slots__ = ('id','homomer','homology','TM','BSA')
     def __init__(self,pdb,mer,hom,tm,sasa):
          self.id=str(pdb)
          self.homomer=bool(mer)
          self.homology=float(hom)
          self.TM=float(tm)
          self.BSA=float(sasa)
          
     def __repr__(self):
          return '{}dimer: {}, with {:.2f}% homology and BSA of {}A2'.format('Homo' if self.homomer else 'Hetero', self.id, self.homology, self.BSA)
          
def scrubInput(path,local_alignment=True,similarity=True):
     df = pd.DataFrame(index=['Hom', 'Het'],columns=['Null', 'Dimer'])
     rows_list = []
     BSAs_raw=[line.rstrip().split() for line in open('/scratch/asl47/PDB/Inputs/{}.txt'.format(path))]
     BSAs={}
     for bsa in BSAs_raw[1:]:
          BSAs[bsa[0]]=(float(bsa[2])-float(bsa[1]))/2
          
     for result in glob.glob('/scratch/asl47/PDB/results/{}/*.results'.format(path)):
          raw=[line.rstrip() for line in open(result)]

          try:
               d={'homomer':bool(raw[0]=='1')}
               tag=result.split('.')[0]
               if '_' in tag:
                    d['id']= tag[-6:]
               else:
                    d['id']= tag[-4:]
                    
               if len(raw)==1:
                    print("Empty ID on {}".format(d['id']))
                    continue
               
               if any('#' in line for line in raw):
                    slice_index=1+2*local_alignment+similarity
                    d['homology']=float(raw[slice_index].split()[-1][:-1])
                    d['H_pred']='yes' if d['homology']>=30 else 'no'

               if any('Total' in line for line in raw):
                    bsa=tuple(filter(lambda element: 'Total' in element,raw))
                    d['BSA']=((float(bsa[2].split()[-1])+float(bsa[1].split()[-1]))-float(bsa[0].split()[-1]))/2
               elif d['id'] in BSAs:
                    d['BSA']=BSAs[d['id']]
                    
               if any('TM-score' in line for line in raw):
                    aligns=tuple(filter(lambda element: 'TM-score' in element,raw))
                    d['TM']=float(aligns[0].split()[1])
                    d['S_pred']= 'no' if d['TM']<=.2 else 'yes' if d['TM']>=.5 else 'maybe'

          except Exception as e:
               print(e)
               print("except on ", d['id'])
               continue

          rows_list.append(d)
     return pd.DataFrame(rows_list)
          
         # rows_list.append
import numpy as np
import scipy.stats as ss
def plotData(data,stat_func=''):
     #g=sns.regplot(x="TM", y="BSA", data=data)
     low=data.loc[data['TM']<=.2]
     high=data.loc[data['TM']>=.5]
     #g = sns.violinplot(y="BSA",data=low)#xlim=(0,1))
     #plt.figure()
     #g3 = sns.violinplot(y="BSA",data=high)#xlim=(0,1))
     #g2 = sns.jointplot(x="TM", y="BSA",data=data,xlim=(0,1))
     
     print("H: ", np.nanmedian(low['BSA'])," and ",np.nanstd(low['BSA']))
     print("HH: ",np.nanmedian(high['BSA'])," +- ",np.nanstd(high['BSA']))
     sns.relplot(x='TM', y='homology', size="BSA",sizes=(40, 400), alpha=.5, height=6, data=data)
     plt.figure()
     ax = sns.violinplot(x="S_pred", y="BSA",hue='H_pred', data=data, palette="muted",split=False,  scale="count",scale_hue=False,inner="quartile",bw=.1)
     
     plt.show(block=False)
     print(ss.mannwhitneyu(high['BSA'],low['BSA'],alternative='greater'))
     if stat_func == 'AD':
          print(ss.anderson_ksamp([data.loc[data['S_pred']=='yes']['BSA'], data.loc[data['S_pred']=='no']['BSA']]))
          print(ss.anderson_ksamp([data.loc[data['S_pred']=='yes']['BSA'], data.loc[data['S_pred']=='maybe']['BSA']]))
          print(ss.anderson_ksamp([data.loc[data['S_pred']=='maybe']['BSA'], data.loc[data['S_pred']=='no']['BSA']]))
     elif stat_func == 'KS':
          print(ss.ks_2samp(data.loc[data['S_pred']=='yes']['BSA'], data.loc[data['S_pred']=='no']['BSA']))
          print(ss.ks_2samp(data.loc[data['S_pred']=='yes']['BSA'], data.loc[data['S_pred']=='maybe']['BSA']))
          print(ss.ks_2samp(data.loc[data['S_pred']=='maybe']['BSA'], data.loc[data['S_pred']=='no']['BSA']))
     
     



import pandas as pd
def loadLevy():
     data=[]
     for line in open('../../../../Downloads/levydata'):
          vs=line.split()
          data.append((vs[0],float(vs[3]),*[int(i) for i in vs[6:9]]))
     return data

def plotLevy(df,plo='BSA'):
     plt.figure()
     if plo == 'BSA':
          plt.hist(df['BSA'],bins=100)   
          plt.show(block=False)

     elif plo == 'Hom':
          ax = sns.violinplot(x="homol", y="BSA", data=df.loc[df['iden']==0], palette="Set2",  scale="count",inner="quartile",bw=.1)
          plt.show(block=False)
     data=df.loc[df['iden']==0]
     #return data
     print(np.median(data.loc[data['homol']==1]['BSA']),np.median(data.loc[data['homol']==0]['BSA']))
     print(stats.mannwhitneyu(data.loc[data['homol']==1]['BSA'],data.loc[data['homol']==0]['BSA'],alternative='greater'))
