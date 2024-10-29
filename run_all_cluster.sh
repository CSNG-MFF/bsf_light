#!/bin/bash

python run.py 'params/low_vol_radial.yml' 'results/low_vol_radial_cluster.pickle'

python run.py 'params/default.yml' 'results/default_cluster.pickle'

python run.py 'params/high_vol_radial.yml' 'results/high_vol_radial_cluster.pickle'

python run.py 'params/very_high_vol_radial.yml' 'results/very_high_vol_radial_cluster.pickle'

python run.py 'params/very_very_high_vol_radial.yml' 'results/very_very_high_vol_radial_cluster.pickle'

python run.py 'params/low_vol_z.yml' 'results/low_vol_z_cluster.pickle'

python run.py 'params/normal_vol_z.yml' 'results/normal_vol_z_cluster.pickle'

python run.py 'params/high_vol_z.yml' 'results/high_vol_z_cluster.pickle'

python run.py 'params/very_high_vol_z.yml' 'results/very_high_vol_z_cluster.pickle'

python run.py 'params/very_very_high_vol_z.yml' 'results/very_very_high_vol_z_cluster.pickle'
