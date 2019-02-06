import matplotlib.pyplot as plt
import seaborn as sns
import glob
import pandas as pd

class Dimer(object):
     __slots__ = ('id','homomer','homology','BSA')
     def __init__(self,pdb,mer,hom,sasa):
          self.id=str(pdb)
          self.homomer=bool(mer)
          self.homology=float(hom)
          self.BSA=float(sasa)
          
     def __repr__(self):
          return '{}dimer: {}, with {:.2f}% homology and BSA of {}A2'.format('Homo' if self.homomer else 'Hetero', self.id, self.homology, self.BSA)
          
def scrubInput(path,local_alignment=True,similarity=True):
     df = pd.DataFrame(index=['Hom', 'Het'],columns=['Null', 'Dimer'])
     rows_list = []

     for result in glob.glob('/scratch/asl47/PDB/results/{}/*.results'.format(path)):
          raw=[line.rstrip() for line in open(result)]

          try:
               d={'homomer':bool(raw[0]=='1')}
               tag=result.split('.')[0]
               if '_' in tag:
                    d['id']= tag[-6:-2]
               else:
                    #continue
                    d['id']= tag[-4:]
               if len(raw)==1:
                    print("Empty ID on {}".format(d['id']))
                    continue
               if('#' in raw[1]):
                    slice_index=1+2*local_alignment+similarity
                    d['homology']=float(raw[slice_index].split()[-1][:-1])
                    #if similarity:
                    #     d['homology']=float(raw[2].split()[-1][:-1])
                    #else:
                    #d['homology']=float(raw[1].split()[-1][:-1])
               if len(raw)>5:
                    #print(raw)
                    try:
                         d['BSA']=((float(raw[-1].split()[-1])+float(raw[-2].split()[-1]))-float(raw[-3].split()[-1]))/2
                    except:
                         print(raw)
          except Exception as e:
               print(e)
               print(d['id'])

          rows_list.append(d)
     return pd.DataFrame(rows_list)
          
         # rows_list.append
          
def plotData(data):
     g=sns.regplot(x="homology", y="BSA", data=data)
     #g = sns.jointplot(x="homology", y="BSA", data=data,xlim=(0,100))

     plt.show(block=False)
