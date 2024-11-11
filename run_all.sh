#!/bin/bash

python scripts/run.py 'params/default.yml' 'results/default.pickle'

python scripts/run.py 'params/default_Lutomirski.yml' 'results/default_Lutomirski.pickle'

python scripts/run.py 'params/default_vandeHulst.yml' 'results/default_vandeHulst.pickle'

python scripts/run_error.py 'params/error.yml' 'results/error.pickle'

python scripts/run.py 'params/low_vol_radial.yml' 'results/low_vol_radial.pickle'

python scripts/run.py 'params/high_vol_radial.yml' 'results/high_vol_radial.pickle'

python scripts/run.py 'params/very_high_vol_radial.yml' 'results/very_high_vol_radial.pickle'

python scripts/run.py 'params/very_very_high_vol_radial.yml' 'results/very_very_high_vol_radial.pickle'

python scripts/run.py 'params/low_vol_z.yml' 'results/low_vol_z.pickle'

python scripts/run.py 'params/normal_vol_z.yml' 'results/normal_vol_z.pickle'

python scripts/run.py 'params/high_vol_z.yml' 'results/high_vol_z.pickle'

python scripts/run.py 'params/very_high_vol_z.yml' 'results/very_high_vol_z.pickle'

python scripts/run.py 'params/very_very_high_vol_z.yml' 'results/very_very_high_vol_z.pickle'

python scripts/plot_figure1.py

python scripts/plot_figure2.py

python scripts/plot_figure3.py
