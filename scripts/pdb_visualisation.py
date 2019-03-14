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

def getDomains(path='Complex_Heterodimers'):
     data=[line.rstrip().split() for line in open('/scratch/asl47/PDB/Inputs/{}.txt'.format(path))]
     domain_dict={}
     for l in data:
          ds=l[-1].split(':')
          if ds == ['?']:
               domain_dict[l[0]]='?'
               continue
          
          domain_dict[l[0]]={ds[i-1][-1]:ds[i][:(-1 if i<(len(ds)-1) else None)] for i in range(1,len(ds))}
          continue
          
          

     return domain_dict
def domainMatch(data,c1,c2):
     if c1 not in data and c2 not in data:
          return 'NA'
     elif c1 not in data or c2 not in data:
          return 'none'
     elif data[c1] == data[c2]:
          return 'full'
     elif not set(data[c1]).isdisjoint(data[c2]):
          return 'partial'
     else:
          return 'none'

import requests
from collections import defaultdict
def pullSCOP(pdb_id):
     data=requests.get('http://www.ebi.ac.uk/pdbe/api/mappings/scop/{}'.format(pdb_id)).json()
     domain_arch=defaultdict(list)
     for superfam in data[pdb_id]['SCOP'].values():
          scop_id=superfam['superfamily']['sunid']
          for chain_element in superfam['mappings']:
               domain_arch[chain_element['chain_id']].append(scop_id)
     return dict(domain_arch)

def pullCATH(pdb_id):
     data=requests.get('http://www.ebi.ac.uk/pdbe/api/mappings/cath_b/{}'.format(pdb_id)).json()
     domain_arch=defaultdict(list)
     for cath_id, superfam in data[pdb_id]['CATH-B'].items():
          for chain_element in superfam['mappings']:
               domain_arch[chain_element['chain_id']].append(cath_id)
     return dict(domain_arch)
     
import json
def writeDomains(pdb_list,class_type='SCOP'):
     domains={}
     pullDomains=pullSCOP if class_type=='SCOP' else pullCATH
     for pdb in pdb_list:
          domains[pdb]=pullDomains(pdb)

     with open('domain_architectures_{}.json'.format(class_type), 'w') as file_out:
          file_out.write(json.dumps(domains))

def readDomains(file_name='SCOP'):
     with open('domain_architectures_{}.json'.format(file_name)) as f:
          return json.load(f)

from collections import Counter
def getUniqueHomodimerDomains(top_X=None):
     if top_X is None:
          return {tuple(v) for vals in readDomains('HCath').values() for v in vals.values()}
     else:
          counted_domains=Counter([tuple(v) for vals in readDomains('HCath').values() for v in vals.values()])
          return {mc[0] for mc in counted_domains.most_common(top_X)}
          
     
     
def scrubInput(path,local_alignment=True,similarity=True,domain='SCOP'):
     df = pd.DataFrame(index=['Hom', 'Het'],columns=['Null', 'Dimer'])
     rows_list = []
     BSAs_raw=[line.rstrip().split() for line in open('/scratch/asl47/PDB/Inputs/{}.txt'.format(path))]
     domain_dict=readDomains(domain)
     homomer_domains=getUniqueHomodimerDomains(100)
     BSAs={}
     for bsa in BSAs_raw[1:]:
          BSAs[bsa[0]]=(float(bsa[2])-float(bsa[1]))/2
          
     for result in glob.glob('/scratch/asl47/PDB/results/{}/*.results'.format(path)):
          d={'id':result.split('.')[0].split('/')[-1].lower()}
          raw=[line.rstrip() for line in open(result)]
          
          try:
               d['homomer']=bool(raw[0]=='1')
                    
               if len(raw)==1:
                    print("Empty ID on {}".format(d['id']))
                    continue
               
               if any('#' in line for line in raw):
                    homologies=tuple(filter(lambda element: '#' in element,raw))
                    slice_index=2+2*local_alignment+similarity
                    d['homology']=float(homologies[2*local_alignment+similarity].split()[-1][:-1]) 
                    d['H_pred']='yes' if d['homology']>=30 else 'no'

               if any('Total' in line for line in raw):
                    bsa=tuple(filter(lambda element: 'Total' in element,raw))
                    d['BSA']=((float(bsa[2].split()[-1])+float(bsa[1].split()[-1]))-float(bsa[0].split()[-1]))/2
                    
               if d['id'][:4] in domain_dict:
                    chains=raw[1].split()
                    d['domain']=domainMatch(domain_dict[d['id'][:4]],*chains)
                    if d['domain']=='full' or d['domain']=='partial':
                         d['d_group']=any([tuple(domain_dict[d['id'][:4]][c]) in homomer_domains for c in chains])
                    #elif 'Homodimer' in path:
                    #     d['d_group']=None
                    else:
                         d['d_group']='none'
                    
                         
                         
               else:
                    d['domain']='unknown'

               if any('TM-score' in line for line in raw):
                    aligns=tuple(filter(lambda element: 'TM-score' in element,raw))
                    d['TM']=float(aligns[0].split()[1])
                    d['S_pred']= 'no' if d['TM']<=.2 else 'yes' if d['TM']>=.5 else 'maybe'
                    

          except Exception as e:
               print("Exception",e)
               print("except on ", d['id'], "len ",len(raw))
               continue

          rows_list.append(d)
     return pd.DataFrame(rows_list)
          
         # rows_list.append
import numpy as np
import scipy.stats as ss
def plotData(data,stat_func=''):
     #g=sns.regplot(x="TM", y="BSA", data=data)
     
                   
     #g = sns.violinplot(y="BSA",data=low)#xlim=(0,1))
     #plt.figure()
     #g3 = sns.violinplot(y="BSA",data=high)#xlim=(0,1))
     #g2 = sns.jointplot(x="TM", y="BSA",data=data,xlim=(0,1))
     sns.relplot(x='TM', y='BSA', size="BSA",sizes=(40, 400),hue='domain', alpha=.75, height=6, data=data)
     plt.figure()
     ax = sns.violinplot(x="domain", y="BSA",hue='d_group', data=data, palette="muted",split=False,  scale="width",scale_hue=False,inner="quartile",bw=.1)
     
     plt.show(block=False)

     low=data.loc[data['S_pred']=='no']
     high=data.loc[data['S_pred']=='yes']

     for val, dic in zip(('High','Low'),(high,low)):
          print("{}: ({})".format(val,len(dic))," = ", np.nanmedian(dic['BSA']))

     
     
     print("\np-value: {}\n".format(ss.mannwhitneyu(high['BSA'],low['BSA'],alternative='greater')[1]))
     for overlap in ('full','partial','none','NA'):
          print(overlap,"({})".format(len(data.loc[data['domain']==overlap])),np.nanmedian([float(i) for i in data.loc[data['domain']==overlap]['BSA']]))

     domain_yes=x=data.loc[(data['domain']=='full') | (data['domain']=='partial')]
     domain_no=x=data.loc[~(data['domain']=='full') & ~(data['domain']=='partial')]
     print("\np-value: {}".format(ss.mannwhitneyu(domain_yes['BSA'],domain_no['BSA'],alternative='greater')[1]))
     print(np.nanmedian(domain_yes['BSA']),np.nanmedian(domain_no['BSA']))
     



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
