from scripts import simulation_runner

def test_compiler(capsys):
    with capsys.disabled():
        print('Testing compiler')
        assert simulation_runner.compileExecutable(None,L=40,root_path='tests',fullwrite=True,recompile=True), 'Compiling failed'
        print('Compiling worked!')
