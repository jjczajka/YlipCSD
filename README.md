# YlipCSD


This is a computational strain design algorithm train and validated on the oleaginous yeast *Yarrowia lipolytica*. <br>

* **PredictionFiles** directory contains the trained algorithm and necessary files for generating strain predictions.
* **modelTrainingFiles** directory contains the code and initial database that was used for training and validating the model.  
* **environmentalYAMLFiles** directory contains .yml files to replicate the python environments.
  * Linux (condaPY36lin_.yml)
  * Mac (condaPY36mac.yml)
  * Windows (condaPY36win.yml) operating systems.

The Conda (Miniconda or Anaconda) python package manager was utilized for initiating the python environment. Information on downloading and installing the Conda package manager can be [found here.](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

To install and activate the environment from the yaml file, navigate to the directory in the conda prompt and use the following commands:

>conda env create -f condaPY36xxx.yml
>conda activate condaPY36xxx

To export the environment to a Jupyter Kernel, use the following commands:
>python-m ipykernel install --user --name "KERNEL NAME: <condaPY36...>"
