#!/bin/bash

python run.py 'params/default.yml' 'results/default.pickle'

python run.py 'params/default_Lutomirski.yml' 'results/default_Lutomirski.pickle'

python run.py 'params/default_vandeHulst.yml' 'results/default_vandeHulst.pickle'

python run_error.py 'params/error.yml' 'results/error.pickle'

python run.py 'params/low_vol_radial.yml' 'results/low_vol_radial.pickle'

python run.py 'params/high_vol_radial.yml' 'results/high_vol_radial.pickle'

python run.py 'params/low_vol_z.yml' 'results/low_vol_z.pickle'

python run.py 'params/high_vol_z.yml' 'results/high_vol_z.pickle'

python run.py 'params/normal_vol_z.yml' 'results/normal_vol_z.pickle'
