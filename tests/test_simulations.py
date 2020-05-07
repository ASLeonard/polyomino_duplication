from scripts import simulation_runner, phylogeny_tracking
import os

def test_data_generation(capsys):
    executable_name = f'bin/{simulation_runner.extractExectutableName()}_L40'
    default_args = '-N 3 -P 100 -B 10 -X .33 -F 1 -T 10'
    full_args = default_args + ' -M 0.005 -D 1 -Y .75 -G 500'

    with capsys.disabled():
        print('Running for 2 duplication rates')
        for offset, dup_rate in enumerate((0,.05)):
            assert simulation_runner.runnerCore(executable_name,40,.75,full_args + f' -J {dup_rate} -L {dup_rate} -V {offset}',verbose=True)

            for fname in ('Homology','Interactions','PIDs','Selections','Size','Strengths'):
                os.remove(f'tests/{fname}_Run{offset}.txt')
        print('All done!')

def test_evolution_record_analysis(capsys):
    dup_rates = (0,0.05)
    assert phylogeny_tracking.loadManyRecords([.75],dup_rates,0.005,40,use_raw=False,FNAME='tests/'), 'making record (without raw generations) failed'
    assert phylogeny_tracking.loadManyRecords([.75],dup_rates,0.005,40,use_raw=True,FNAME='tests/'), 'making record (with raw generations) failed'
    for dup_rate in dup_rates:
        os.remove(f'tests/EvoRecord_Mu0.005000_S0.750000_D{dup_rate:.6f}.txt')


