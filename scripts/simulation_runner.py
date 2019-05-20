import subprocess
import sys
#import argparse



def runnerCore(G_RATE):
    generations={.671875:1250, .6875:2500, .703125:5000, .71875:10000, .734375:21250, .75:50000}
    #1250 2500 5000 10000 21250 50000
    default_args=' -N 3 -P 100 -B 20 -X .25 -F 1 -A 1 -V 0 -T 10 -M .0015 -D 1000 -J {0} -L {0}'.format(G_RATE)


    for S_c, gen in generations.items():
        print('Running simulation at S: {}'.format(S_c))
        full_args = default_args + ' -Y {} -G {}'.format(S_c, gen)
        
        subprocess.run('../bin/DuplicationEvolution -E'+full_args,shell=True)


if __name__ == '__main__':
    if len(sys.argv) == 2:
        print('Running simulations')
        runnerCore(float(sys.argv[1]))
        print('Done')
    else:
        print('Wrong args')
