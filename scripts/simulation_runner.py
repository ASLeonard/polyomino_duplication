import subprocess
import sys
import os
import argparse

def runnerCore(L,S_c,mu,dup_rate,generations,runs):

    default_args=f'-N 3 -P 100 -B 20 -X .33 -F 1 -A 1 -V 0 -T 10 -M {mu} -D {runs} -J {dup_rate} -L {dup_rate} '

    print(f'Attempting simulation for L: {L} S_c: {S_c}, duplication: {dup_rate}')
    full_args = default_args + f'-Y {S_c} -G {generations}'

    os.chdir(os.path.expanduser('~/Documents/PolyDev/duplication'))

    executable_name = f'bin/NEW_DuplicationEvolution_L{L}'

    if not os.path.isfile(executable_name):
        print('Building simulation for arguments')
        subprocess.run(f'make clean && make ILen={L}',shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        print('Compiled successfull\nRunning simulation')

    subprocess.run(f'{executable_name} -E '+full_args,shell=True)
    print('Simulation complete')
    return

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-L','--Length', type=int,required=True)
    parser.add_argument('-S','--Strength', type=float,required=True)
    parser.add_argument('-M','--Mutation', type=float,required=True)
    parser.add_argument('--dup_rates', nargs='+', type=float,required=True)
    parser.add_argument('-G','--Generations', type=int)
    parser.add_argument('-R','--Runs', type=int)

    parser.set_defaults(Generations=0,Runs=0)
    
    args = parser.parse_args()

    for dup_rate in args.dup_rates:
        runnerCore(args.Length,args.Strength,args.Mutation,dup_rate,args.Generations,args.Runs)    

if __name__ == '__main__':
    main()
