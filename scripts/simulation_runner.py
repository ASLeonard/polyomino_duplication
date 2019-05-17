import subprocess
import sys
#import argparse



def runnerCore(S_CRITICAL,DU_FLAG='J'):
    print('Running simulation at S: {}'.format(S_CRITICAL))
    generations={.671875:1500,.6875:2500,.703125:4000, .71875:9000, .734375:17500,.75:25000}

    default_args=' -N 3 -P 100 -B 50 -X .25 -F 1 -A 1 -V 0 -T 10 -Y {} -M .0015 -D 250 -G {}'.format(S_CRITICAL, generations[S_CRITICAL])

    dead_flag={'J':'K','K':'J'}

    for G_RATE in (0.01,):
        special_args=' -L {} -{} {} -{} 0'.format(G_RATE, DU_FLAG, G_RATE, dead_flag[DU_FLAG])
        subprocess.run('../bin/DuplicationEvolution -E'+default_args+special_args,shell=True)


if __name__ == '__main__':
    if len(sys.argv) == 2:
        print('Running simulations')
        for S in range(86,98,2):
            runnerCore(S/128)
        print('Done')
    else:
        print('Wrong args')
