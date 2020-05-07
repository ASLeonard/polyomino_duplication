[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python application](https://github.com/ASLeonard/duplication/workflows/Python%20application/badge.svg)](https://github.com/ASLeonard/duplication/actions?query=workflow%3A%22Python+application%22)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/e64e2359118d4555ae6916c9e5d540a4)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ASLeonard/duplication&amp;utm_campaign=Badge_Grade)

# Duplication and lattice self-assembly

Polyomino model of lattice self-assembly for exploring the role of duplication in the evolution of structural complexity in protein complexes. Further information can be found in this [preprint](https://www.biorxiv.org/content/10.1101/2020.04.22.054783v1).

## Overview
There are two main components to the duplication evolution simulations. There is the evolution framework implemented in c++ (found in the src/, include/, and polyomino_core/ folders), as well as analysis and plotting tools implemented in python (founds in the scripts/ folder). The main methods have been made easily accessible to command line calls (see later), but can also be run or edited at a lower level if desired. Some more details on the polyomino background can be found in the original repository [here](https://github.com/ASLeonard/polyomino_core), with more details on the specific binary binding site model [here](https://github.com/ASLeonard/polyomino_interfaces).

## Install

There are several requirements to using the software
- recent c++ compiler which supports c++2a (tested on g++ 8)
- recent python (developed on 3.6.9)
- several common python packages (detailed in requirements.txt)
  - numpy
  - scipy
  - pandas
  - matplotlib
- pytest (optional)

The installation process is as follows 

```shell
git clone --recurse-submodules -j8 https://github.com/ASLeonard/polyomino_duplication
cd polyomino_duplication
pip install -r requirements.txt
```

### Testing

If pytest is present/installed, the software can be quickly tested using

```shell
python -m pytest tests/
```
Which verifies the compiling works and the simulations/analyses occur properly. Several compiler flags can be added depending on what is available on the user's system if desired, like g++ has support for multi-threading (-fopenmp) and link-time optimization (-flto), but are not required to be used.

## Usage examples

Data generation can be run directly from the command line either via the c++ executable or an easier python interface. An example of the python interface is shown below.

```python
python simulation_runner.py -L 100 -S .75 -M 0.0025 -G 2500 --dup_rates 0 0.10 -R 10 --pathway data/
```
This will generate evolution data for polyominoes with binding site lengths of 100 bits, a strength threshold of 0.75, and a mutation rate of 0.0025 per bit. The simulations run for 2500 generations, once for a duplication rate of 0 and once for a duplication rate 0.10. Each simulation configuration will run for 10 independent evolutions, and the results are written to a folder called data.


Similarly, metrics for interaction formation, decay, and the duplication-heteromerisation epsilon factor can be generated with the following command.
```python
python metric_runner.py -L 100 -S .75 -G 2500 --formation --decay --epsilon
```

Results can then be analysed or plotted within the evolution_plotting.py file
```python
data = loadSimulationResults()
printTransitionTable(data)
plotTimeDistributions(*data)
```
