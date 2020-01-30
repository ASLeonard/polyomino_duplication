from scripts import simulation_runner
import os

def test_compiler(capsys):
    executable_name = f'bin/DuplicationEvolution_L40'
    with capsys.disabled():
        print('Testing compiler')
        assert simulation_runner.compileExecutable(executable_name,L=40,root_path='tests',fullwrite=True,recompile=True), 'Compiling failed'
        print('Compiling worked!')    

def test_data_generation(capsys):
    executable_name = f'bin/DuplicationEvolution_L40'
    default_args = f'-N 3 -P 100 -B 10 -X .33 -F 1 -A 1 -T 10'
    full_args = default_args + f' -M 0.005 -D 1 -Y .75 -G 250'

    with capsys.disabled():
        print('Running for 2 duplication rates')
        for offset, dup_rate in enumerate((0,.05)):
            assert simulation_runner.runnerCore(executable_name,40,.75,250,full_args + f' -J {dup_rate} -L {dup_rate} -V {offset}',verbose=True)
            os.remove(f'tests/EvoRecord_Mu0.005000_S0.750000_D{dup_rate:.6f}.txt')
            for fname in ('Homology','Interactions','PIDs','Selections','Size','Strengths'):
                os.remove(f'tests/{fname}_Run{offset}.txt') 
        print('All done!')


