import subprocess
import os
import argparse

def compileExecutable(executable_name,L,root_path,fullwrite,recompile):
    print('Checking if compiling is required')
    if not os.path.isfile(executable_name) or recompile:
        print('**Compiling requested**\nBuilding simulation for arguments')
        print(f'Current directory is {os.getcwd()}')
        subprocess.run(f'make clean && make ILen={L} ROOT_FILE_PATH={root_path}{"/" if root_path else ""} FW={int(fullwrite)}',shell=True,check=True)
        print('Compiled successfull\nRunning simulation')
    else:
        print('Executable already exists, proceeding')

    return True

def runnerCore(executable_name,L,S_c,generations,full_args,verbose):
    print(f'Attempting simulation for L: {L} S_c: {S_c}')

    if verbose:
        print(f'about to run following command:\n{executable_name} -E {full_args}')

    subprocess.run(f'{executable_name} -E {full_args}',shell=True)
    print('Simulation complete')
    return True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-L','--Length', type=int,required=True)
    parser.add_argument('-S','--Strength', type=float,required=True)
    parser.add_argument('-M','--Mutation', type=float,required=True)
    parser.add_argument('--dup_rates', nargs='+', type=float,required=True)
    parser.add_argument('-G','--Generations', type=int,required=True)
    parser.add_argument('-R','--Runs', type=int)

    parser.add_argument('--pathway', type=str)
    parser.add_argument('--verbose', dest='verbose', default=False, action='store_true')
    parser.add_argument('--fullwrite', dest='fullwrite', default=False, action='store_true')
    parser.add_argument('--recompile', dest='recompile', default=False, action='store_true')

    parser.set_defaults(Runs=1,pathway='')
    args = parser.parse_args()

    executable_name = f'bin/DuplicationEvolution_L{args.Length}'
    try:
        compileExecutable(executable_name,args.Length,args.pathway,args.fullwrite,args.recompile)
    except subprocess.CalledProcessError:
        print('Error in compiling, try and fix it?')
        return

    default_args = f'-N 3 -P 100 -B 20 -X .33 -F 1 -A 1 -T 10'

    full_args = default_args + f' -M {args.Mutation} -D {args.Runs} -Y {args.Strength} -G {args.Generations}'

    for index, dup_rate in enumerate(args.dup_rates):
        offset = index*args.Runs
        runnerCore(executable_name,args.Length,args.Strength,args.Generations,full_args + f' -J {dup_rate} -L {dup_rate} -V {offset}',args.verbose)

if __name__ == '__main__':
    main()

