from scripts import interaction_dynamics
import numpy as np

def test_steady_state(capsys):
    with capsys.disabled():
        print('testing steady state distribution')
    assert np.isclose(interaction_dynamics.steadyStates(4,.75),[2/3,1/3]).all(), 'steady state not close to expectation'

def test_form_time(capsys):
    with capsys.disabled():
        print('testing form time expectation')
    assert np.isclose(round(interaction_dynamics.formTime(40,.75),2),2944.76), 'form time not close to expectation'

def test_drop_time(capsys):
    with capsys.disabled():
        print('testing drop time expectation')
    assert np.isclose(round(interaction_dynamics.dropTime(40,.75),2),3.88), 'drop time not close to expectation'

def test_mutual_exclusion(capsys):
    with capsys.disabled():
        print('testing mutual interaction exclusion')
    assert np.isclose(round(interaction_dynamics.MutualExclusion(5,.75,40),4),0.8908), 'mutual exclusion not close to expectation'

