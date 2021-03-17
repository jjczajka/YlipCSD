# YlipCSD 
## PredictionFiles


This directory contains the necessary files for generating strain design predictions for *Yarrowia lipolytica*. <br>

See the ***README*** file on the [main page](https://github.com/jjczajka/YlipCSD) for instructions on installing and setting up the environment necessary to run the files.

The user needs to only access/run two files to generate predictions.
1. **Supplemental Excel File 6- CSD Template.xlsx** - file for entering the genetic and environmental conditions and potential strain constructs to be examined.
2. **computationalStrainPrediction.ipynb** - jupyter notebook iPython file for running the code.  

Other files in the directory are needed for the algorithm to function and are described here.
* **M21iYL.pickle** - pickle file containing the machine learning trained model described in the accompining manuscript.
* **FBA_function_.py** - function that generates the FBA features from GSM *iYLI647_corr* for input into the prediction algorithm.
* **encodingFunction_.py** - function that encodes the features (environmental, genetic, FBA) for input into the prediction algorithm.
* **iYLI647_corr.mat** - Genome scale model (.mat) for generating *Y. lipolytica* features. Originally published [10.1186/s12918-018-0542-5](10.1186/s12918-018-0542-5) and then curated [here, 10.1007/s12257-019-0208-1](10.1007/s12257-019-0208-1). The curated model is provided here. 
* **Supplemental Excel File 2- DataCharacteristics&Encoding.xlsx** - file containing encoding dictionaries needed in the *encodingFunction.py* script.

A step-by-step guide for utilizing the algorithm is provided with the following manuscript.
