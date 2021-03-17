# YlipCSD 
## modelTrainingFiles

This directory contains the code used for training and analyzing the *Yarrowia lipolytica* machine learning model. The databases are provided as supplemental information in the following manuscript.  <br>

See the ***README*** file on the [main page](https://github.com/jjczajka/YlipCSD) for instructions on installing and setting up the environment necessary to run the files.

The data handling and model training/validation was performed using the following scripts, with the first three being used for pre-processing of data.
1. **ML_pipeline_part1.ipynb** - script for importing the database, creating genetic features, and splitting the data into train/test sets.
2. **ML_pipeline_part2.ipynb** - script for creating FBA features based on genome scale models that replicate the provided genetic background detailed in the strain instance.  
3. **ML_pipeline_part3.ipynb** - script for encoding the features for input into the ML pipeline.


ML training and validation scripts
4. **ML_pipeline_part4.ipynb** - script for creating train/validate datasets and ML pipeline and model training.
5. **FeatureImportanceScores.ipynb** - script for assessing the importance of individual features on model predictions.
6. **ModelEvaluation.ipynb** - script for evaluating model performance on validate or test data.

Other files needed for training are listed below and provided as supplemental information with the accompanying manuscript.
* **Supplemental Excel File 1- Database** - Database containing the original data used to train the model. 
* **Supplemental Excel File 2- DataCharacteristics&Encoding.xlsx** - file containing encoding dictionaries needed in the *ML_pipeline_part2.ipynb & ML_pipeline_part3.ipynb* script.
* **iYLI647_corr.mat** - Genome scale model (.mat) for generating *Y. lipolytica* features in *ML_pipeline_part2.ipynb*. Originally published [10.1186/s12918-018-0542-5](10.1186/s12918-018-0542-5) and then curated [DOI:10.1007/s12257-019-0208-1](https://link.springer.com/article/10.1007%2Fs12257-019-0208-1). The curated model is provided here. 
A step-by-step guide for utilizing the algorithm is provided with the following manuscript.