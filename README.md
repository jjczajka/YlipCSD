# YlipCSD


Y. lipolytica computational strain design algorithm.

Environmental .yml files have been provided to replicate the environments used for training and predictions.
A .yml file is provided for Linux (condaPY36lin_.yml), Mac (condaPY36mac.yml), and Windows (condaPY36win.yml) operating systems.

The Miniconda or Anaconda python package manager can be downloaded here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/
To install and activate the environment, use the following commands from the conda prompt:

conda env create -f condaPY36xxx.yml
conda activate condaPY36xxx

To export the environment to a Jupyter Kernel, use the following commands:
python-m ipykernel install --user --name "KERNEL NAME: condaPY36..."
