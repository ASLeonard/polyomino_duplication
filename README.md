[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python application](https://github.com/ASLeonard/duplication/workflows/Python%20application/badge.svg)](https://github.com/ASLeonard/duplication/actions?query=workflow%3A%22Python+application%22)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/e64e2359118d4555ae6916c9e5d540a4)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ASLeonard/duplication&amp;utm_campaign=Badge_Grade)

# Duplication and lattice self-assembly

## Overview

## Install

There are several requirements to using the software
- recent c++ compiler which supports c++2a (tested on g++ 8)
- recent python (developed on 3.6.9)
- several common python package found in requirements.txt
 - numpy
 - scipy
 - etc
- pytest (optional)

The installation process is as follows 

```shell
git clone --recurse-submodules https://github.com/ASLeonard/polyomino_duplication
cd polyomino_duplication
pip install -r requirements.txt
```

If pytest is present/installed, the software can be quickly tested using

```shell
python -m pytest
```
Which verifies the compiling works and the simulations/analyses occur properly. Several compiler flags can be added depending on what is available on the user's system if desired, like g++ has support for multi-threading (-fopenmp) and link-time optimization (-flto), but are not required to be used.

## Usage



### Examples
