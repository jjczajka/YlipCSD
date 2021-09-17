def perform_MLPrediction(encodedData,output):
    import pandas as pd
    import pickle
    from collections import defaultdict
    import warnings
    import numpy as np
    import os

    from sklearn.preprocessing import StandardScaler, RobustScaler, QuantileTransformer,MinMaxScaler, MaxAbsScaler,Normalizer,PowerTransformer
    from sklearn.pipeline import Pipeline
    import warnings
    from sklearn.linear_model import ElasticNet,Ridge
    from sklearn.neighbors import KNeighborsRegressor
    from sklearn.svm import SVR
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF
    import xgboost as xgb
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.preprocessing import MinMaxScaler,StandardScaler
    from mlxtend.regressor import StackingRegressor
    from sklearn.model_selection import learning_curve

    encodedData['mw_Lipids'] = encodedData['mw']
   
    cols_train__set = [
    'mw_Lipids'
    ,'pH'
    ,'product_deltaGo'
    ,'foldCarbonFed'
    ,'product_name'
    ,'rxt_volume'
    ,'inputThermo(kJ/L)'
    ,'FermentationTime'
    ,'atp_cost'
    ,'precursorsRequiredEncoded'
    ,'nadh_nadph_cost'
    ,'Pathway_enzymatic_steps'
    ,'averageThermBarrier'
    ,'media'
    ,'number_genes_het'
    ,'number_native_genes_overexp'
    ,'ATP_iYLI647'
    ,'NADPH_iYLI647'
    ,'PPP_iYLI647'
    ,'TCA_iYLI647'
    ,'PrdtYield_iYLI647'
    ]

    warnings.simplefilter('ignore')


    useful_cols = []
    useful_cols.extend(cols_train__set)
    data = pd.DataFrame()


    data = encodedData.loc[:,useful_cols]
    for column in data:
        data[column] = data[column].astype(np.float32)
        #obtain features used from the data
    
    
    warnings.simplefilter('ignore')


    useful_cols = []
    useful_cols.extend(cols_train__set)
    data = pd.DataFrame()

    data = encodedData.loc[:,useful_cols]
    for column in data:
        data[column] = data[column].astype(np.float32)

    #open the ML model for prediction
    with open('M21iYL_cs_GlcNormalized_noO2noEMPnoBio.pickle','rb') as f:
        masterGrid = pickle.load(f)

    masterGrid = masterGrid[0]

    #perform prediction on data
    x_testData = data.copy()
    target = 'Product_titer(g/L)'
    x_testData.PrdtYield_iYLI647.fillna(0,inplace=True)

    #prediction
    y_prediction = np.exp(masterGrid[target].predict(x_testData))
    len(y_prediction)
    
    #output dataframe
    MLOutput = pd.DataFrame()
    MLOutput['TiterPrediction(g/L)'] = y_prediction
    MLOutput['% of Original Strain Production'] = y_prediction/y_prediction[0]*100
    MLOutput['FBA predicted Biomass'] = output['Biomass_iYLI647']
    MLOutput['FBA predicted Yield'] = data['PrdtYield_iYLI647']
    MLOutput['Input Reaction Tested'] = output['geneMod']
    # MLOutput.at['rxns'] = optKnockRxns
    MLOutput.index = data.index
    MLOutput.round({'% of Original Strain Production':2})
    return(MLOutput)  
