from scripts import metric_runner, interaction_dynamics, simulation_runner
import os
import warnings

def test_metric_run(capsys):
    executable_name = f'bin/{simulation_runner.extractExectutableName()}_L40'
    with capsys.disabled():
        print('Running metrics')
        metric_runner.runnerCore(executable_name,'111','-D 250 -Y .75 -G 5000',True)

def test_formation_calculation(capsys):
    data = interaction_dynamics.loadFormations(.75,'tests/')
    assert interaction_dynamics.calculateMeanMetric(data,(1,)).any(), 'Calculation failed'
    os.remove('tests/Formation_0.750000.BIN')

def test_decay_calculation(capsys):
    data = interaction_dynamics.loadDecays(40,.75,'tests/')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        assert interaction_dynamics.calculateMeanMetric(data,(2,)).any(), 'Calculation failed'
    os.remove('tests/Decay_0.750000.BIN')

def test_gamma_calculation(capsys):
    assert interaction_dynamics.calculateGammaFactors(40,.75,'tests/').any(), 'Calculation failed'
    os.remove('tests/Gamma_0.750000.BIN')
