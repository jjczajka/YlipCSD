# YlipCSD


This repository contains a computational strain design (CSD) algorithm train and validated on the oleaginous yeast *Yarrowia lipolytica*. <br>

The CSD alorithm takes environmental and genetic conditions and a target biochemical molecule and outputs titer predictions. The algorithm can be coupled with flux  balancing analysis strain design algorithms to generate titer predictions for strain knock-outs. Tutorials on specific aspects of the CSD algorithm are found within the appropriate directories.  


* **PredictionFiles** directory contains the trained algorithm and necessary files for generating strain predictions.
* **modelTrainingFiles** directory contains the code and initial database that was used for training and validating the model.  
* **environmentalYAMLFiles** directory contains .yml files to replicate the python environments.
  * Linux (condaPY36lin_.yml)
  * Mac (condaPY36mac.yml)
  * Windows (condaPY36windows.yml) operating systems.

The Conda (Miniconda or Anaconda) python package manager was utilized for initiating the python environment. Information on downloading and installing the Conda package manager can be [found here.](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

To install and activate the environment from the yaml file, navigate to the directory in the conda prompt and use the following commands:

>conda env create -f condaPY36xxx.yml  
>conda activate condaPY36xxx

To install Jupyter notebooks use the following lines after activating the conda environment:
>pip3 install --upgrade pip  
>pip3 install jupyter

To export the environment to a Jupyter Kernel, use the following commands:
>python -m ipykernel install --user --name "KERNEL NAME: <condaPY36...>"



UPDATE: Currently the algorithm is only set up to analyze either Overexpression or Knockouts for strain design predictions. In Supplemental Excel File 6- CSD Template.xlsx there is a sheet where you can specify either KO or OE to enable the option.
