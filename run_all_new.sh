#!/bin/bash

python run.py 'params/default_lowRAM.yml' 'results/default_lowRAM.pickle'

python run.py 'params/default.yml' 'results/default.pickle'

python run.py 'params/default_Lutomirski.yml' 'results/default_Lutomirski.pickle'

python run.py 'params/default_vandeHulst.yml' 'results/default_vandeHulst.pickle'

python run.py 'params/error.yml' 'results/error.pickle'

