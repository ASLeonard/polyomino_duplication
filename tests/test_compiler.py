from scripts import simulation_runner

def test_compiler(capsys):
    executable_name = f'bin/DuplicationEvolution_L40'
    with capsys.disabled():
        print('Testing compiler')
        assert simulation_runner.compileExecutable(executable_name,L=40,root_path='tests',fullwrite=True,recompile=True), 'Compiling failed'
        print('Compiling worked!')
