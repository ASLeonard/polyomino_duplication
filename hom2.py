import os

def do():
    for line in [l.split(', ') for l in open('PDB_list.txt')][0]:
        os.system('/u/fs1/asl47/Documents/PolyDev/duplication/dimer_analysis.pl {}_1 0 0'.format(line))
        

def main():
    do()
    print('Finished')
                        
if __name__ == "__main__":
    os.chdir('/scratch/asl47/PDB')
    main()
