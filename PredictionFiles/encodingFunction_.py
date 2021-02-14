
def encodeTransform(workingData):

    import os as os
    import pandas as pd

    #categorically encodes the volume to 1-5 (smallest to largest) based on cultivation volume 
    def rxtVolumeEncoding(x):
        if x <= 0.01:
            return 1
        elif x <= 0.075:
            return 2
        elif x <= 0.25:
            return 3
        elif x < 1:
            return 4
        else:
            return 5
    #categorically encodes the reactor vessel type (micro-reactors, shaking flasks, batch fedbatch or continuous vessels)
    #smallest (1) to largest (3)
    def reactorTypeEncoding(x):
        if x == 2:
            return 3
        elif x == 4:
            return 3
        elif x==3:
            return 3
        elif x == 1:
            return 2
        elif x == 5:
            return 1
  
    #corrects the database encoding to categorically encode from lowest oxygen level (1) to highest (3)    
    def oxygenEncodingFix(x):
        if x==1:
            return 3 #3 not oxygen sufficient
        elif x==2:
            return 1 #1 now oxygen insufficeint
        elif x==3:
            return 2 #2 now intermidiate

    #changes the database encoding to categorically encode from lowest nitrogen level (1) to highest (3)                
    def nitrogenEncodingFix(x):
        if x==1:
            return 3 #nitrogen sufficent
        elif x==2:
            return 1 #nitrogen limited
        elif x==3:
            return 2 #intermidiate


    def funcReturnStr(x):
        return str(x)


    #obtain dictionaries on how to encode the data; not all dictionaries are utilized for model features.
    encoding_Data = pd.ExcelFile('Supplemental Excel File 2- DataCharateristics & Encoding.xlsx').parse('Encoding')
    
    strainDict = dict(zip(encoding_Data.strain_background,encoding_Data.strain_class))
    mediaDict = dict(zip(encoding_Data.media,encoding_Data.media_class))
    productDict = dict(zip(encoding_Data.Product, encoding_Data.prdt_class))

    carbonSourceMWDict = dict(zip(encoding_Data.carbonSource,encoding_Data.carbonSourceMW))
    N2sourceDict = dict(zip(encoding_Data.N2Source,encoding_Data.N2source_class))
    promoterDict = dict(zip(encoding_Data.Promoters,encoding_Data.prom_class))
    integrationSiteDict = dict(zip(encoding_Data.integrationSite,encoding_Data.int_class))

    ####################################################################################################
    ####################################################################################################
    ####################################################################################################
    workingData['product_name2'] = workingData['product_name']
    workingData['product_name'] = workingData.product_name.map(productDict).fillna(workingData.product_name)
    workingData['media'] = workingData.media.map(mediaDict).fillna(workingData.media)

    workingData['carbonSourceOneMolecularWeight'] = workingData.cs1.map(carbonSourceMWDict).fillna(workingData.cs1)
    workingData['carbonSourceTwoMolecularWeight'] = workingData.cs2.map(carbonSourceMWDict).fillna(workingData.cs2)


    temp1 = workingData.cs_conc1/workingData.carbonSourceOneMolecularWeight*workingData['cs1_heatCombustion(kJ/mol)']
    temp2 = workingData.cs_conc2/workingData.carbonSourceTwoMolecularWeight*workingData['cs2_heatCombustion(kJ/mol)']
    temp2.fillna(0,inplace=True)

    df = pd.DataFrame()
    df['one'] = temp1
    df['two'] = temp2

    temp3={}

    for y,z in enumerate(df.one):
        if df.two.iloc[y]!=0:
            temp3[df.index[y]]=df.one.iloc[y]+df.two.iloc[y]
        else:
            temp3[df.index[y]]=df.one.iloc[y]

    temp3 = pd.Series(temp3)

    workingData['inputThermo(kJ/L)'] = temp3


    ####################################################################################################
    #precursors required
    ####################################################################################################
    temp2 = pd.DataFrame()
    temp2 = workingData.precursor_required.apply(funcReturnStr).str.split(';',expand=True).fillna(0) #TAG
    temp2 = temp2.apply(pd.to_numeric)
    workingData['precursorsRequiredEncoded'] = temp2.sum(axis=1)

    ########################################   thermoOptions   #########################################
    totalTherm={}
    averageTherm={}

    for dataPoint in workingData.index:


        stoichNADPH=(workingData.nadh_nadph_cost.loc[dataPoint])
        stoichATP=(workingData.atp_cost.loc[dataPoint])
        stoichprecursor={}
        temp1={}
        temp2={}
        
        #The thermo option for lipids utilized the fatty acid compononents due to availability of Gibbs free eneregy 
        if workingData.loc[dataPoint]['product_name2']=='Lipids':
            stoichATP = stoichATP/3
            stoichNADPH = stoichNADPH/3
            temp2=workingData.loc[dataPoint].precursor_required.strip().split(';')
            stoichprecursor[0] = float(temp2[0])/3

            temp1[0] = -3341.2 #deltaGo
            prec = ['Acetyl-CoA']

        # non-lipid products
        else:
            prec = workingData.loc[dataPoint].central_carbon_precursor.strip().split(';')

            if isinstance(workingData.loc[dataPoint].precursor_required,str):
                stoichprecursor=workingData.loc[dataPoint].precursor_required.strip().split(';')
                temp1 = workingData.loc[dataPoint]['ccm_precursor_deltaGo'].strip().split(';')

            else:
                stoichprecursor[0]=workingData.loc[dataPoint].precursor_required
                temp1[0] = workingData.loc[dataPoint]['ccm_precursor_deltaGo']


        thermoTemp = 0
        for i,j in enumerate(prec):
            thermoTemp += float(stoichprecursor[i])*float(temp1[i])
            if j == 'Acetyl-CoA':
                stoichCoA = float(stoichprecursor[i])
            else:
                stoichCoA = 0

        ATP_nadph_tempThermo = stoichATP*-31.8+-28.8*stoichNADPH #deltaGo

        # withCoA
        totalTherm[dataPoint] = round(ATP_nadph_tempThermo + workingData.loc[dataPoint]['product_deltaGo']-thermoTemp+stoichCoA*-3202.2)
        
        averageTherm[dataPoint]=round(totalTherm[dataPoint]/(workingData.loc[dataPoint]['Pathway_enzymatic_steps']))



    workingData['totalThermBarrier'] = pd.Series(totalTherm)
    workingData['averageThermBarrier'] = pd.Series(averageTherm)
    workingData['averageThermBarrier'].fillna(0,inplace=True)
    workingData['averageThermBarrier']

    ####################################################################################################
    ###################################   N2 Source Encoding   #########################################
    ####################################################################################################
    temp1 = pd.DataFrame()

    temp1 = workingData.N2Source.apply(funcReturnStr).str.split(';',expand=True).fillna('NaN') #FA

    for col in temp1.columns:
        temp1[col]=temp1[col].map(N2sourceDict)

    workingData['N2SourceEncoded(max)'] = temp1.max(axis=1).fillna(1)



    temp1={}
    for dataPoint in workingData.index:
        if workingData.loc[dataPoint]['N2SourceEncoded(max)']==1 and workingData.loc[dataPoint]['N2_content']>1:
            temp1[dataPoint]=4 # high N2, organic
        elif workingData.loc[dataPoint]['N2SourceEncoded(max)']==1 and workingData.loc[dataPoint]['N2_content']<2:
            temp1[dataPoint]=2 #low N2, organic
        elif workingData.loc[dataPoint]['N2SourceEncoded(max)']==0 and workingData.loc[dataPoint]['N2_content']<2:
            temp1[dataPoint]=1 #low N2, inorganic
        elif workingData.loc[dataPoint]['N2SourceEncoded(max)']==0 and workingData.loc[dataPoint]['N2_content']>1:
            temp1[dataPoint]=3 #high N2, inorganic

    temp1 = pd.Series(temp1)
    workingData['N2_contentEncoded']=temp1


    ####################################################################################################
    #####################  Integration site & Promoter strength Encoding  ##############################
    ####################################################################################################
    temp1 = pd.DataFrame()
    temp2 = pd.DataFrame()

    temp1 = workingData.integration_site_Filled.apply(funcReturnStr).str.split(';',expand=True).fillna('NaN')

    for col in temp1:
        temp1[col]=temp1[col].map(integrationSiteDict)

    workingData['integrationSiteEncoded(Sum)'] = temp1.sum(axis=1).fillna(0)

    ####################################################################################################
    ####################################################################################################

    workingData.rxt_volume = workingData.rxt_volume.apply(rxtVolumeEncoding)
    workingData.reactor_type = workingData.reactor_type.apply(reactorTypeEncoding)
    workingData.oxygen = workingData.oxygen.apply(oxygenEncodingFix)
    workingData.pH = workingData.pH.fillna(0)

    workingData['csConcTotal'] = workingData['cs_conc1']+ workingData['cs_conc2']
    workingData['dir_evo'].fillna(0,inplace=True)


    workingData['PrdtFlux_iYLI647'].fillna(workingData['PrdtFlux_iYLI647'].min(),inplace=True)
    workingData['O2Uptake_iYLI647'] = abs(workingData['O2Uptake_iYLI647'])
    workingData['GlcUptake_iYLI647'] = abs(workingData['GlcUptake_iYLI647'])
    workingData['PrdtYield_iYLI647'].fillna(0,inplace=True)

    return(workingData)
