See upstream MiniAn repo for full details on MiniAn package.

# Quick Start Guide for Paul's version

1. Clone this repository to your computer
1. Install Miniforge (https://github.com/conda-forge/miniforge/tree/main?tab=readme-ov-file)
1. Install MATLAB R2023a (https://www.mathworks.com/downloads/) and add PV_minian GitHub directory (with subfolders) to MATLAB path
1. Download and install K-Lite Codec Pack (Standard) (https://codecguide.com/download_k-lite_codec_pack_standard.htm) (allows MATLAB to read videos using FFV1 codec that Miniscope V4 generates)
1. Launch MiniForge Prompt
1. Change directory to cloned repository `cd [path to GitHub folder]/PV_minian`
1. Create new environment from environment.yml file: `mamba env create -n PV_minian -f environment.yml` (be patient this step may take a while to execute)
1. Activate your new environment `mamba activate PV_minian`
1. Install MATLAB Python Engine (to run NormCorre) `pip install matlabengine`
1. Fire up jupyter: `jupyter notebook` and open the notebook "pipeline.ipynb"

# Re-launching after installation is complete
1. Launch MiniForge Prompt
1. Activate your new environment `mamba activate PV_minian`
1. Fire up jupyter: `jupyter notebook` and open the notebook "pipeline.ipynb"
