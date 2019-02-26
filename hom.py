#!/usr/bin/python3

import sys
import os





def tex(pdb_in,pdb_ID, ch1, ch2):
        with open('{}_{}.pdb'.format(pdb_in,ch1), 'w') as ch1_out:
                with open('{}_{}.pdb'.format(pdb_in,ch2), 'w') as ch2_out:
                        for line in open('{}.pdb{}'.format(pdb_in,pdb_ID)):
                                if 'ATOM' in line:
                                        if ch1 in line.split()[:-1]:
                                                ch1_out.write(line)
                                        elif ch2 in line.split()[:-1]:
                                                ch2_out.write(line)

def main():
        if len(sys.argv)!=5:
	        sys.exit('Incorrect arguments')
        tex(*sys.argv[1:5])
        print('Finished')
                        
if __name__ == "__main__":
	os.chdir('/scratch/asl47/PDB')
	main()
