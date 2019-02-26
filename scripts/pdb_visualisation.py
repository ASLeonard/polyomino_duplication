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
     BSAs_raw=[line.rstrip().split() for line in open('/scratch/asl47/PDB/Hets2.txt')]
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
                    #continue
                    d['id']= tag[-4:]
               if len(raw)==1:
                    print("Empty ID on {}".format(d['id']))
                    continue
               if('#' in raw[1]):
                    slice_index=1+2*local_alignment+similarity
                    d['homology']=float(raw[slice_index].split()[-1][:-1])

               if any('Total' in line for line in raw):
                    bsa=tuple(filter(lambda element: 'Total' in element,raw))
                    d['BSA']=((float(bsa[2].split()[-1])+float(bsa[1].split()[-1]))-float(bsa[0].split()[-1]))/2

               else:
                    d['BSA']=BSAs[d['id']]
               if any('TM-score' in line for line in raw):
                    d['TM']=float(raw[-1].split()[1])
          except Exception as e:
               print(e)
               print("except on ", d['id'])

          rows_list.append(d)
     return pd.DataFrame(rows_list)
          
         # rows_list.append
          
def plotData(data):
     #g=sns.regplot(x="TM", y="BSA", data=data)
     low=data.loc[data['TM']<.25]
     high=data.loc[data['TM']>.5]
     g = sns.violinplot(y="BSA",data=low)#xlim=(0,1))
     plt.figure()
     g3 = sns.violinplot(y="BSA",data=high)#xlim=(0,1))
     g2 = sns.jointplot(x="TM", y="BSA",data=high,xlim=(0,1))
     
     print(np.mean(low['BSA'])," and ",np.mean(high['BSA']))
     
     plt.show(block=False)
