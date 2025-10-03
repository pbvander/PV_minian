See upstream MiniAn repo for full details on MiniAn package.

# Quick Start Guide for Paul's version

1. Clone this repository to your computer
1. Download Anaconda Navigator (https://www.anaconda.com/products/navigator)
1. Install MATLAB https://www.mathworks.com/products/matlab.html and add PV_minian GitHub directory (with subfolders) to MATLAB path
1. Launch Anaconda Prompt
1. Install mamba (much faster than conda for environment management): `conda install mamba`
1. Change directory to cloned repository `cd [path to GitHub folder]/PV_minian`
1. Create new environment from environment.yml file: `mamba env create -n [environment-name] -f environment.yml`
1. Install MATLAB Python Engine (to run NormCorre) `pip install matlabengine`
1. Fire up jupyter: `jupyter notebook` and open the notebook "pipeline.ipynb"
