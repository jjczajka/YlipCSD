# YlipCSD 
## PredictionFiles


This directory contains the necessary files for generating strain design predictions for *Yarrowia lipolytica*. <br>

See the ***README*** file on the [main page](https://github.com/jjczajka/YlipCSD) for instructions on installing and setting up the environment necessary to run the files.

The user needs to only access/run two files to generate predictions.
1. **Supplemental Excel File 6- CSD Template.xlsx** - file for entering the genetic and environmental conditions and potential strain constructs to be examined.
2. **computationalStrainPrediction.ipynb** - jupyter notebook iPython file for running the code.  

Other files in the directory are needed for the algorithm to function and are described here.
* **M21iYL.pickle** - pickle file containing the machine learning trained model described in the accompanying manuscript.
* **FBA_function_.py** - KNOCK OUT function that generates the FBA features from GSM *iYLI647_corr* for input into the prediction algorithm.
* **FBA_functionOE_.py** - OVER Expression function that generates the FBA features from GSM *iYLI647_corr* for input into the prediction algorith1m.
* **encodingFunction_.py** - function that encodes the features (environmental, genetic, FBA) for input into the prediction algorithm.
* **iYLI647_corr.mat** - Genome scale model (.mat) for generating *Y. lipolytica* features. Originally published [DOI 10.1186/s12918-018-0542-5](https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-018-0542-5) and then curated [DOI:10.1007/s12257-019-0208-1](https://link.springer.com/article/10.1007%2Fs12257-019-0208-1). The curated model is provided here. 
* **Supplemental Excel File 2- DataCharacteristics&Encoding.xlsx** - file containing encoding dictionaries needed in the *encodingFunction.py* script.


A step-by-step guide for utilizing the algorithm is provided with the following manuscript.

UPDATE: Currently the algorithm is only set up to analyze either Overexpression or Knockouts for strain design predictions. In **Supplemental Excel File 6- CSD Template.xlsx** there is a sheet where you can specify either KO or OE to enable the option.


If you would like to simulate for non-glucose carbon sources, use the **updatedPredictionFiles** Directory. The new directory has a modified version of the model, where the FBA simulation uses the individual carbon sources normalized to 10 mmol/h/gDCW of glucose. For example, if you wanted to screen on glycerol, the FBA model would use ~20 mmol/h/gDCW as an uptake rate for glycerol.

The model that has been added was described in the manuscript and had similar performance metrics. 
