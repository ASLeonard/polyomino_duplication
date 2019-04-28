#!/usr/bin/env python3                                          
import sys
import os

from collections import defaultdict
from operator import itemgetter

import subprocess


def needleAlign(pdb_1,chain_1,pdb_2,chain_2,needle_EXEC='./needle'):
    for pdb, chain in ((pdb_1,chain_1),(pdb_2,chain_2)):
        try:
            subprocess.run('grep -i {0}_{1} -A1 all_fasta.txt > FASTA/{0}_{1}.fasta.txt'.format(pdb,chain),shell=True,check=True)
        except CalledProcessError as e:
            print('need to wget and append, then recall')
    subprocess.run('{EXEC} FASTA/{0}_{1}.fasta.txt FASTA/{2}_{3}.fasta.txt -gapopen 10.0 -gapextend 0.5 -outfile NEEDLE/{0}_{1}_{2}_{3}.needle'.format(pdb_1,chain_1,pdb_2,chain_2,EXEC=needle_EXEC),shell=True)

def makeTypes(pdb):
    similaritythresh = 2.0
    type_ = {}
    with open('newINT/{}.int'.format(pdb)) as file_:
        for line in file_:
            l = line.strip().split('\t')
            chain = l[0]
            res = l[1]
            if len(l) == 4:
                ll = l[3].split()
                if float(ll[2]) > 0:
                    type_[chain+res] = tuple(reversed(ll[0].split('_')))
            if len(l) > 4:
                jjj = []
                for j in l[3:]:
                    jj = j.split()
                    if float(jj[2]) > 0:
                        jjj.append((jj,float(jj[2])))
                if jjj != []:
                    if not(len(jjj) == 2 and abs(jjj[0][1]-jjj[1][1]) < similaritythresh and jjj[0][0][0] == jjj[1][0][3] and jjj[0][0][3] == jjj[1][0][0]): #SYMMETRIC INTERFACE: ALWAYS PICK FIRST COLUMN, OTHERWISE SORT
                        jjj.sort(key=itemgetter(1),reverse=True)
                    type_[chain+res] = tuple(reversed(jjj[0][0][0].split('_')))
    return type_

##scrape FASTA sequence for both chains, and the alignment interactions 
def readNeedle(pdbs):
    seqs, align = ['', ''], ''
    read = 0
    with open('NEEDLE/{}_{}_{}_{}.needle'.format(*pdbs)) as file_:
        for line in file_:
            ##can skip empty or commented lines in .needle file
            if line == '\n' or line[0] == '#':
                continue
            
            ll = line.strip().split()
            ##add alternating sequence lines
            if len(ll) > 0 and ('_' in ll[0] or ll[0] == 'SEQUENCE'):
                seqs[read]+=ll[2]
                read = not read
            ##if between chain A and chain B, record the alignment interactions
            elif read == 1:
                align += line[21:71].replace(' ','_').strip()

    return seqs, align
        

def makeInteractions(type_, seq, chain):

    def numalpha(s):
        if int(s) > 9:
            s = str(10+(int(s)-10)%26) # WRAP ROUND NUMBERS LARGER THAN 35 BACK TO a, b, c...
            s = chr(97+int(s)-10)
        return s

    c = 0
    inter_A = ''
    inter_B = ''
    for i in range(len(seq)):
        if seq[i] != '-':
            chain_key = chain+str(c)
            if chain_key in type_:
                inter_A += numalpha(type_[chain_key][0])
                inter_B += type_[chain_key][1]
            else:
                inter_A += ' '
                inter_B += ' '
            c += 1
        else:
            inter_A += ' '
            inter_B += ' '
    return inter_A, inter_B

def writePialign(pdbs,inter1a,inter1b,inter2a,inter2b,align,seq1,seq2):
    iter_chunk=50
    with open('PALIGN/{}_{}_{}_{}.pialign'.format(*pdbs),'w') as pialign_out:
        idx, idx_A, idx_B = 0,0,0
        for idx in range(0,max(len(seq1),len(seq2)),iter_chunk):
            slice_ = slice(idx,idx+iter_chunk)
            idx_A+=iter_chunk-seq1[slice_].count('-')
            idx_B+=iter_chunk-seq2[slice_].count('-')

            a=align[idx:idx+50]
           
            pialign_out.write('{}\n{}\n{}   {}\n{}\n{}   {}\n{}\n{}\n\n'.format(inter1b[slice_],inter1a[slice_],seq1[slice_],idx_A,align[slice_],seq2[slice_],idx_B,inter2a[slice_],inter2b[slice_]))

def writePmatrix(pdbs,inter1a,inter1b,inter2a,inter2b):
    mat = defaultdict(int)
    key1, key2 = set(), set()
    
    for a,b,c,d in zip(inter1a,inter1b,inter2a,inter2b):
        mat[(a+b,c+d)] += 1
        key1.add(a+b)
        key2.add(c+d)

    with open('PALIGN/{}_{}_{}_{}.pmatrix'.format(*pdbs),'w') as pmatrix_out:
        pmatrix_out.write('\t')
        for i in sorted(key2):
            pmatrix_out.write('{}\t'.format(i))
        pmatrix_out.write('\n')
        for i in sorted(key1):
            pmatrix_out.write('{}\t'.format(i))
            for j in sorted(key2):
                pmatrix_out.write('{}\t'.format(mat[(i,j)]))
            pmatrix_out.write('\n')
        pmatrix_out.write('\n')
    
                
def generateAssistiveFiles(pdbs):
    type1=makeTypes(pdbs[0])
    type2=makeTypes(pdbs[2])
    (seq1,seq2),align=readNeedle(pdbs)

    int_1A,int_1B=makeInteractions(type1,seq1,pdbs[1])
    int_2A,int_2B=makeInteractions(type2,seq2,pdbs[3])

    writePialign(pdbs,int_1A,int_1B,int_2A,int_2B,align,seq1,seq2)
    writePmatrix(pdbs,int_1A,int_1B,int_2A,int_2B)
    

if __name__ == "__main__":
    if len(sys.argv) == 5:       
        needleAlign(*sys.argv[1:])
        generateAssistiveFiles(sys.argv[1:])
    else:
        print('Wrong arguments')

