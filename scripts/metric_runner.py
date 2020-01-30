from scripts.simulation_runner import compileExecutable
import subprocess
import os
import argparse

def runnerCore(executable_name,run_args,full_args,verbose):
    print(f'Running metrics for chosen data')

    command = f'{executable_name} -M{run_args} {full_args}'
    if verbose:
        print(f'about to run following command:\n\t{command}')

    subprocess.run(command,shell=True)
    print('Simulation complete')
    return True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-L','--Length', type=int,required=True)
    parser.add_argument('-S','--Strength', type=float,required=True)
    parser.add_argument('-G','--Generations', type=int,required=True)
    parser.add_argument('-R','--Runs', type=int)

    parser.add_argument('--formation', dest='formation', default=False, action='store_true')
    parser.add_argument('--decay', dest='decay', default=False, action='store_true')
    parser.add_argument('--gamma', dest='gamma', default=False, action='store_true')

    parser.add_argument('--pathway', type=str)    
    parser.add_argument('--verbose', dest='verbose', default=False, action='store_true')
    parser.add_argument('--recompile', dest='recompile', default=False, action='store_true')

    parser.set_defaults(Runs=100,pathway='')
    args = parser.parse_args()

    executable_name = f'bin/DuplicationEvolution_L{args.Length}'
    try:
        compileExecutable(executable_name,args.Length,args.pathway,False,args.recompile)
    except subprocess.CalledProcessError:
        print('Error in compiling, try and fix it?')
        return

    run_args = f'{int(args.formation)}{int(args.decay)}{int(args.gamma)}'
    default_args = f'-D {args.Runs} -Y {args.Strength} -G {args.Generations}'
    runnerCore(executable_name,run_args,default_args,args.verbose)

if __name__ == '__main__':
    main()

