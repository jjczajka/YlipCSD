def FBA_FeatureExtraction(FBATrainData,optKnockRxns,FBA_models):


    import cobra
    from cobra.flux_analysis import single_gene_deletion, single_reaction_deletion, double_gene_deletion,double_reaction_deletion
    from cobra import Reaction, Metabolite, Model
    import pandas as pd
    import pickle
    from collections import defaultdict
    import warnings
    import numpy as np
    import os

    ####USER DEFINED OPTIONS.
    #Do you want to perform GENETIC KNOCKOUTS?
    KO_option = 1

    #Do you want to perform OVEREXPRESSIONS
    OE_option = 1

    #Do you want to perform PRODUCT-BASED SIMULATIONS?
    Product_option = 1

    #OVEREXPRESSION and PRODUCT OPTIONS
    epsilon = [0.05,10,1,0.75,0.02,0.5] #percent to increase/decrease gene expression. new bounds if flux is all 0. bounds if flux is 0, but upper/lowerbound is not 0, biomass percent change

    ####################################################################################################
    #######################################      USED FXNS      ########################################
    ###########################################################################################
    #Create common-name to Genome-Scale-Model (GSM) gene name & is gene in GSM dictionary.
    def createGeneDict():
        productInfo = pd.ExcelFile('Supplemental Excel File 2- DataCharateristics & Encoding.xlsx').parse('Encoding')

        df = pd.DataFrame()
        df['bname'] = productInfo.bname
        df['traditionalName'] = productInfo.traditionalName
        df['iYLI647'] = productInfo.in_iYLI647
    #     df['iMK735'] = productInfo.in_iMK735
    #     df['iYali4'] = productInfo.in_iYali4
        # df['iNL895'] = productInfo.in_iNL895
    #     df['iYL_2.0'] = productInfo['in_iYL_2.0']

        df = df.T
        df = df.rename(columns=df.loc['traditionalName'])
        df = df.drop('traditionalName')#,axis=0)
        geneDict = df.to_dict()


        df2 = pd.DataFrame()
        df2['CCM'] = productInfo['Central Carbon']
    #     df2['iYL_2.0'] = productInfo['iYL_2metabolites']
        df2['iYLI647'] = productInfo['iYLI647metabolites']
        # df2['iNL895'] = productInfo['iNL895metabolites']
    #     df2['iMK735'] = productInfo['iMK735metabolites']
    #     df2['iYali4'] = productInfo['iYali4metabolites']

        df2 = df2.T
        df2.rename(columns=df2.loc['CCM'],inplace=True)
        df2.drop('CCM',inplace=True)

        fbaModelMetaboliteDict = df2.to_dict()

        return(geneDict,fbaModelMetaboliteDict)

    #Generate Gene-product-assocated dictionary for used Genome-Scale-Model
    def generateOEGeneGPR(GSM,model):
        GPR_dict=defaultdict(list)
        for x in geneDict.keys():
            if geneDict[x][GSM]==1:
                tempGene = model.genes.get_by_id(geneDict[x]['bname'])
                rxn_list=[]
                for reaction in tempGene.reactions:
                    temp_dict={}
                    temp_dict['mets']=[x.id for x in reaction.metabolites]
                    temp_dict['mets_coefs']=[x for x in reaction.get_coefficients(reaction.metabolites)]
                    temp_dict['lower_bound']=reaction.lower_bound
                    temp_dict['upper_bound']=reaction.upper_bound
                    temp_dict['id']=reaction.id
                    temp_dict['name']=reaction.name
                    temp_dict['subsystem']=reaction.subsystem
                    temp_dict['gpr']=reaction.gene_reaction_rule
                    rxn_list.append(temp_dict)
                GPR_dict[x]=rxn_list

        return(GPR_dict)

        #Simulate default Genome-scale-model flux with biomass as objective function
    def defaultObjFunction(dGSM):
        defaultObj = 'biomass_C'

        model = cobra.io.load_matlab_model(dGSM+'_corr.mat')
        model.objective = model.reactions.get_by_id(defaultObj)
        
        defaultFlux = model.optimize()
        
        
        
        
        
        return(model, defaultFlux.objective_value,defaultFlux)

        #Perform Genetic Knockout for model simulation

    def search(list, g):
        """
        Search genome scale model for engineered gene of interest.

        Parameters
        ----------
        list: list
            List of all genes in the genome-scale-model.
        gene: list
            The engineered gene name entered as either cobra model gene name or as a generic name provided in the gene_dict (see createGeneDict function).

        Returns
        -------
        True or False and the gene name in the list.
        """
        gene = [g]
        for ind,temp in enumerate(list):
            if temp == gene[0]:
                return True,temp
        return False,temp


    def searchGeneDict(gene,GSM,GPR_dict,gene_list):
        """
        Extract the gene-reaction rules from the selected genome-scale-model for each particular gene.
        Parameters
        ----------
        gene: list
        The engineered gene name entered as either cobra model gene name or as a generic name provided in the gene_dict (see createGeneDict function).
        GSM: list
        Genome-scale model name
        GPR_dict: .dict
        Gene-protein-reaction dictionary, see generateOEGeneGPR
        gene_list: list
        List of all genes in the genome-scale-model.
        Returns
        -------
        Gene associated reactions or FALSE if engineered gene of interest is not in the GSM.
        """
        try:
            #Check if there is externally provided generic name-genome scale gene name dict
            #Check if the generic name is the externally provided dictionary
            if (GPR_dict[gene] and geneDict[gene][GSM]==1):
                return(GPR_dict[gene])
        except Exception as e:
            # search for the engineered gene of interest is within the genome-scale model.
            # temp1 = True or False (if gene in model)
            # temp2 = gene in the gene_ist
            temp1,temp2 = search(gene_list,gene)
            if (temp1==True):
                return(GPR_dict[temp2])
            else:
                return(False)

    #Implement gene OVEREXPRESSION for each overexpressed native gene
    def performGeneKOs(modelKO,GSM,genesKO,geneMO):
        """
        Performs GSM model knock-outs.

        Parameters
        ----------
        modelKO: cobra.Model
            Model on which to act
        GSM: list
            Name of GSM that is being acted on
        genesKO: list
            Vector of 1 or 0 corresponding to whether the list of genes are Knocked out
        geneMO:
            Vector of gene names that are being genetically modified

        Returns
        -------
        Modifed GSM with the corresponding genetic knock-outs.
        """
        gene_list=[z.id for z in modelKO.genes]
        for i,KO in enumerate(genesKO):
            try:
                if (KO==1 and geneDict[geneMO[i]][GSM]==1) | (KO=='1' and geneDict[geneMO[i]][GSM]==1):
                    try:
                        cobra.manipulation.delete_model_genes(modelKO,(pd.Series(geneDict[geneMO[i]]['bname'])))
                    except Exception as e:
                        print(repr(e),geneMO[i],GSM)

            except Exception as e:
                temp1,temp2 = search(gene_list,geneMO[i])


                if (KO==1 and temp1==True) | (KO==1 and temp1==True):
                    try:
                        cobra.manipulation.delete_model_genes(modelKO,(pd.Series(temp2)))
                    except Exception as e:
                        print(repr(e),geneMO[i],GSM)
                else:
                    print(geneMO[i],'not in GSM, no KO modification performed')
        return(modelKO)

    def performGeneOE(tempOEModel,GSM,genesOE,genesMO,hetGenes,tempKOSol,GPR_dict,ep0,OE_f,ep1,ep2,ep5,f1a,f2a,f3a,f4a,f5a,f6a):
        """
        Performs GSM model OE.

        Parameters
        ----------
        tempOEModel: cobra.Model
            Model on which to act
        GSM: list
            Name of GSM that is being acted on
        genesOE: list
            Vector of 1 or 0 corresponding to whether the list of genes are overexpressed
        geneMO:
            Vector of gene names that are being genetically modified
        tempKOSol: corba.Solution
            Flux solution for the GSM before overexpression is implemented and optimized for biomass growth
        GPR_dict: defaultdict
            Gene-reaction rules for each gene modified.
        ep0: int
            Percent to increase or decrease overexpressed flux
        OE_f: int
            Count of number of times the overexpression results in infeasible flux solutions
        ep1: int
            new bounds if the non-overexpressed model has flux bounds of 0.

        ep2: int
            new reaction bound if the non-overexpressed solution has a lower or upper flux value of 0.
        ep5: int

        f1a: int
            Count of number of times model fails to overexpress a reaction that had a prior flux solution of 0 with bounds set to 0.
        f2a: int
            Count of number of times model fails to overexpress a reaction that had a prior flux solution of 0 with an upper bound set to 0.
        f3a: int
            Count of number of times model fails to overexpress a reaction that had a prior flux solution of 0 with a lower bound set to 0.
        f4a: int
            Count of number of times model fails to overexpress a reaction that had a prior flux solution that was negative.
        f5a: int
            Count of number of times model fails to overexpress a reaction that had a prior flux solution that was positive.
        f6a: int
            Count of number of times model fails to overexpress a reaction that had a prior flux solution with a 0, and fluxes were reset to original bounds (i.e, no resulting modifications).


        Returns
        -------
        tempOEModel: cobra.Model
            Model with the resulting overexpression implemented.
        OE_f: int
            Count of number of times the overexpression results in infeasible flux solutions
        f1a: int
            Count of number of times model fails to overexpress a reaction that had a prior flux solution of 0 with bounds set to 0.
        f2a: int
            Count of number of times model fails to overexpress a reaction that had a prior flux solution of 0 with an upper bound set to 0.
        f3a: int
            Count of number of times model fails to overexpress a reaction that had a prior flux solution of 0 with a lower bound set to 0.
        f4a: int
            Count of number of times model fails to overexpress a reaction that had a prior flux solution that was negative.
        f5a: int
            Count of number of times model fails to overexpress a reaction that had a prior flux solution that was positive.
        f6a: int
            Count of number of times model fails to overexpress a reaction that had a prior flux solution with a 0, and fluxes were reset to original bounds (i.e, no resulting modifications).

        """
        gene_list=[z.id for z in modelKO.genes]

        for i,OE in enumerate(genesOE):

            GPR_dict_list=None
            try:
                '''
                #types1 = [type(k) for k in GPR_dict.keys()]
                #print(types1)
                print(genesMO[i])
                print(geneDict)
                #print(GPR_dict.keys())
                print(GPR_dict[genesMO[i]])
                print(geneDict[genesMO[i]],'i')
                '''
                if ((GPR_dict[genesMO[i]]) and (geneDict[genesMO[i]][GSM]==1)):
                    GPR_dict_list = GPR_dict[genesMO[i]]
                    
            except Exception as e:
                if(OE=='1' and GPR_dict_list==None) | (OE==1 and GPR_dict_list==None):
                    GPR_dict_list = searchGeneDict(genesMO[i],GSM,GPR_dict,gene_list)# print(i,genesMO,OE)


            if (GPR_dict_list!=None):
                for rxn in GPR_dict_list:
                    lower = tempOEModel.reactions.get_by_id(rxn['id']).lower_bound
                    upper = tempOEModel.reactions.get_by_id(rxn['id']).upper_bound
                    

                    #obtain the previous solution reaction flux value
                    rxnKOFlux = tempKOSol.fluxes[rxn['id']]
                    #if the flux is positive, implement OE to lower bound (i.e., set lower bound higher than default flux value)
                    if rxnKOFlux>0:
                        brk=0
                        tempOEModel.reactions.get_by_id(rxn['id']).lower_bound = (rxnKOFlux+((rxnKOFlux)*ep0))
                        tempeps0 = ep0

                        while (tempOEModel.optimize().status!='optimal') and (brk<1):
                            # print('5fail')
                            tempeps0 = tempeps0/2
                            tempOEModel.reactions.get_by_id(rxn['id']).lower_bound = (rxnKOFlux+((rxnKOFlux)*tempeps0))
                            if tempeps0<1e-5:
                                # print('5 total failure')
                                tempOEModel.reactions.get_by_id(rxn['id']).lower_bound = rxnKOFlux
                                f5a+=1
                                brk=2
                                break
                    #if the flux is negative, implement OE to upper bound (i.e., set upper bound lower than default flux value)
                    elif rxnKOFlux<0:
                        brk=0
                        tempOEModel.reactions.get_by_id(rxn['id']).upper_bound = (rxnKOFlux+((rxnKOFlux)*ep0))
                        tempeps0 = ep0

                        while (tempOEModel.optimize().status!='optimal') and (brk<1):
                            tempeps0 = tempeps0/2
                            tempOEModel.reactions.get_by_id(rxn['id']).upper_bound = (rxnKOFlux+((rxnKOFlux)*tempeps0))
                            if tempeps0<1e-5:
                                f4a+=1
                                brk=2
                                tempOEModel.reactions.get_by_id(rxn['id']).upper_bound = rxnKOFlux
                                break
                    
#                     else:
#                         brk=0
#                         if (tempOEModel.reactions.get_by_id(rxn['id']).lower_bound==0 and tempOEModel.reactions.get_by_id(rxn['id']).upper_bound==0):
#                             tempOEModel.reactions.get_by_id(rxn['id']).bounds = (-ep1,ep1)
#                             if (tempOEModel.optimize().status!='optimal'):
#                                 f1a+=1
#                                 tempOEModel.reactions.get_by_id(rxn['id']).bounds = (0,0)

#                         elif (tempOEModel.reactions.get_by_id(rxn['id']).lower_bound==0):

#                             tempOEModel.reactions.get_by_id(rxn['id']).lower_bound=ep2

#                             tempeps2 = ep2
#                             while ((tempOEModel.optimize().status!='optimal') and (brk<1)):
#                                 tempeps2 = tempeps2/10
#                                 tempOEModel.reactions.get_by_id(rxn['id']).lower_bound=tempeps2


#                                 if tempeps2<1e-5:
#                                     tempOEModel.reactions.get_by_id(rxn['id']).lower_bound=0
#                                     brk = 2
#                                     f2a+=1
#                                     break
#                         elif (tempOEModel.reactions.get_by_id(rxn['id']).upper_bound==0):

#                             tempOEModel.reactions.get_by_id(rxn['id']).upper_bound=-ep2
#                             tempeps2 = ep2
#                             while (tempOEModel.optimize().status!='optimal') and (brk<1):
#                                 tempeps2 = tempeps2/10
#                                 tempOEModel.reactions.get_by_id(rxn['id']).upper_bound=-tempeps2
#                                 if tempeps2<1e-5:
#                                     tempOEModel.reactions.get_by_id(rxn['id']).upper_bound=0
#                                     brk=2
#                                     f3a+=1
#                                     break

#                         else:
#                             lower = tempOEModel.reactions.get_by_id(rxn['id']).lower_bound
#                             upper = tempOEModel.reactions.get_by_id(rxn['id']).upper_bound
#                             tempeps5 = ep5
#                             tempOEModel.reactions.get_by_id(rxn['id']).lower_bound=tempeps5

#                             brk2=0
#                             brk=0
#                             while (tempOEModel.optimize().status!='optimal') and brk<1:
#                                 tempOEModel.reactions.get_by_id(rxn['id']).lower_bound=tempeps5
#                                 tempeps5 = tempeps5/10
#                                 if (tempeps5 < 1e-5) and (tempOEModel.optimize().status!='optimal'): #was-8
#                                     tempOEModel.reactions.get_by_id(rxn['id']).lower_bound=lower
#                                     tempeps5 = ep5
#                                     brk = 1
#                                     while (tempOEModel.optimize().status!='optimal') and brk2<1:
#                                         tempOEModel.reactions.get_by_id(rxn['id']).upper_bound=-tempeps5
#                                         tempeps5 = tempeps5/10
#                                         if (tempeps5 < 1e-5) and (tempOEModel.optimize().status!='optimal'):#was -8
#                                             tempOEModel.reactions.get_by_id(rxn['id']).upper_bound=upper
#                                             brk2 = 1
#                                             f6a+=1
#                                             break


            else:
                print('Gene:',genesMO[i],'not in Genome scale model, OE simulation performed without accounting for gene')
        tempOESol = tempOEModel.optimize()
        #!!!! Does not return KO model... (too slow)
        if tempOESol.status!='optimal':
            OE_f+=1

        return(tempOEModel,OE_f,f1a,f2a,f3a,f4a,f5a,f6a)


    #Product flux
    def maximizeProduct(model,defaultBioObj,ep3,ep4,fbaModelMetaboliteDict,dataPoint,counterProductFailTemp,gsm,prod_f,isRbflvOption):

        """
        Adds the pseudo-reaction simulating product flux to the GSM, sets the biomass to a set value, and optimizes for the pseudo-reaction.

        Parameters
        ----------
        model: cobra.Model
            Model on which to act
        defaultBioObj: float
            Prior model biomass objective function (before addition of pseudoreaction)
        ep3: int

        ep4: int

        fbaModelMetaboliteDict: defaultdict
            Dictionary mapping names of metabolites to the model metabolite names (e.g., ATP to atp[c])
            Generated in "createGeneDict()" function
        dataPoint: int
            Index of database construct
        counterProductFailTemp: int
            Count of number of times the pseudo product reaction results in infeasible flux solution and a resulting decrease in the biomass constraint
        gsm: list
            The name of the GSM being modified
        prod_f: int
            Count of number of times that the product pseudo-reaction resulted in an infeasible flux solution with 0 biomass flux.

        isRbflvOption: int
            0 or 1 indicating whether the product is riboflavin, resulting in the correct application of reactant consumption.
        Returns
        -------
        finalProductFluxSolnTemp: cobra.Solution
            Model flux soluiton from the resulting genetic engineering and product reaction implementation.
        counterProductFailTemp: int
            Count of number of times the pseudo product reaction results in infeasible flux solution and a resulting decrease in the biomass constraint    prod_f int
            Count of number of times that the product pseudo-reaction resulted in an infeasible flux solution with 0 biomass flux.
        """
        modelP=model.copy()
        modelP.reactions.get_by_id('biomass_C').upper_bound = (defaultBioObj*(ep3))
        modelP.reactions.get_by_id('biomass_C').lower_bound = (defaultBioObj*(ep3))


        #FBATrainData global?
        # stoichNADPH=(FBATrainData.nadh_nadph_cost.loc[dataPoint])
        # stoichATP=(FBATrainData.atp_cost.loc[dataPoint])
        stoichNADPH=round(FBATrainData.nadh_nadph_cost.loc[dataPoint])
        stoichATP=round(FBATrainData.atp_cost.loc[dataPoint])
        stoichATP
        stoichNADPH
        prec = FBATrainData.loc[dataPoint].central_carbon_precursor.strip().split(';')

        #create product reaction
        reaction__product = []
        reaction__product = Reaction('Prdt_r')
        reaction__product.name = 'Prdt_r'
        reaction__product.subsystem = 'Exchange'
        reaction__product.lower_bound = 0
        reaction__product.upper_bound = 1000
        prdt_m = Metabolite('prdt_m', formula = '', name = 'Prdt_m', compartment = 'e')

        stoichprecursor={}

        modelP.add_reactions([reaction__product])

        #adds energy and cofactors (NADPH only)
        if isRbflvOption==0:
            reaction__product.add_metabolites({
            prdt_m: 1.0,
            modelP.metabolites.get_by_id(fbaModelMetaboliteDict['ATP'][gsm].strip('\'"')).id: -stoichATP,
            modelP.metabolites.get_by_id(fbaModelMetaboliteDict['NADPH'][gsm].strip('\'"')).id: -stoichNADPH,
            modelP.metabolites.get_by_id(fbaModelMetaboliteDict['NADP'][gsm].strip('\'"')).id : stoichNADPH,
            modelP.metabolites.get_by_id(fbaModelMetaboliteDict['ADP'][gsm].strip('\'"')).id : stoichATP
            })
        else:
            reaction__product.add_metabolites({
            prdt_m: 1.0,
            #
            modelP.metabolites.get_by_id(fbaModelMetaboliteDict['ATP'][gsm].strip('\'"')).id: -stoichATP,
            modelP.metabolites.get_by_id(fbaModelMetaboliteDict['NADPH'][gsm].strip('\'"')).id: -stoichNADPH,
            modelP.metabolites.get_by_id(fbaModelMetaboliteDict['NADP'][gsm].strip('\'"')).id : stoichNADPH,
            #modelP.metabolites.get_by_id(fbaModelMetaboliteDict['ADP'][gsm].strip('\'"')).id : stoichATP
                
            })
        if isinstance(FBATrainData.loc[dataPoint].precursor_required,str):
            stoichprecursor=FBATrainData.loc[dataPoint].precursor_required.strip().split(';')
        else:
            stoichprecursor[0]=FBATrainData.loc[dataPoint].precursor_required
        for i,j in enumerate(prec):
            met = fbaModelMetaboliteDict[j][gsm].strip('\'"')
            reaction__product.add_metabolites({modelP.metabolites.get_by_id(met).id: -round(float(stoichprecursor[i]))})
            if j=='Acetyl-CoA':
                reaction__product.add_metabolites({model.metabolites.get_by_id(fbaModelMetaboliteDict['CoenzymeA'][gsm].strip('\'"')).id: round(float(stoichprecursor[i]))})

        #adds to model
        demand = modelP.add_boundary(modelP.metabolites.prdt_m,type="demand")

        modelP.objective = 'Prdt_r'
        # modelP.objective = 'DM_prdt_m'


        finalProductFluxSolnTemp = modelP.optimize()

        c=1
        while finalProductFluxSolnTemp.status!='optimal':
            #print('Infeasible',c)
            c+=1
            ep3-=.05

            modelP.reactions.get_by_id('biomass_C').upper_bound = (defaultBioObj*(ep3))#-ep3),defaultBioObj*(ep3+ep4)) #lower_bound, upper_bound sets Biomass
            modelP.reactions.get_by_id('biomass_C').lower_bound = (defaultBioObj*(ep3))#

            if ep3<0:
                counterProductFailTemp+=1
                prod_f+=1
                finalProductFluxSolnTemp = modelP.optimize()
                break
        finalProductFluxSolnTemp = modelP.optimize()



        return(finalProductFluxSolnTemp,counterProductFailTemp,prod_f)


    def FBAFeatureExtraction(featureModelSoln,GSM):
        if (GSM=='iYLI647'):
            bio2 = featureModelSoln.fluxes['biomass_C']
            EMP2 = featureModelSoln.fluxes['GAPD']/2 # FBA,PFK
            PPP2 = featureModelSoln.fluxes['GND'] # GND 646
            TCA2 = featureModelSoln.fluxes['CSm'] # FUMm
            NADPH2 = featureModelSoln.fluxes['GND']+featureModelSoln.fluxes['G6PDH2']+featureModelSoln.fluxes['ICDHy']+featureModelSoln.fluxes['ICDHym']
            if (featureModelSoln.fluxes['MTHFDm']>0):
                NADPH2 = NADPH2 + featureModelSoln.fluxes['MTHFDm']
            if (featureModelSoln.fluxes['MTHFD']>0):
                NADPH2 = NADPH2 + featureModelSoln.fluxes['MTHFD']
            NADH2 = featureModelSoln.fluxes['GAPD']+featureModelSoln.fluxes['PDHm']+featureModelSoln.fluxes['PGCD']+featureModelSoln.fluxes['MDHm']+featureModelSoln.fluxes['ICDHxm']+featureModelSoln.fluxes['PDHm']
            ATP2 = featureModelSoln.fluxes['ATPS3m']+featureModelSoln.fluxes['PYK']
            if (featureModelSoln.fluxes['PGK']<0):
                ATP2 = ATP2 - featureModelSoln.fluxes['PGK']
            if (featureModelSoln.fluxes['SUCOASm']<0):
                ATP2 = ATP2 - featureModelSoln.fluxes['SUCOASm']
            if (featureModelSoln.fluxes['FACOAL140']<0):
                ATP2 = ATP2 - featureModelSoln.fluxes['FACOAL140']#correct?
            Precursors2 = {}
            PrdtFlux2 = featureModelSoln.fluxes['Prdt_r']
            O2 = featureModelSoln.fluxes['EX_o2(e)']
            Glc = featureModelSoln.fluxes['EX_glc(e)']


        return(EMP2,PPP2,TCA2,NADPH2,ATP2,PrdtFlux2,bio2,O2,Glc)


    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    workingData2 = pd.DataFrame()
    # workingData2
    #Import product info containing the class of each product along with the in or out of model information
    geneDict,fbaModelMetaboliteDict = createGeneDict()
    output=pd.DataFrame()




    # FBA_models=['iNL895']
    #FBATrainData
    ############################### FBA loop start ###################################################
    counterProductFail=0
    for GSM in FBA_models:
        prod_fail = 0
        OE_fail=0
        defaultModel, defaultObj, defaultFluxSol = defaultObjFunction(GSM)
    #####################
        EMP,PPP,TCA,NADPH,ATP,PrdtFlux,PrdtYield,bio,O2uptake,Glcuptake = {},{},{},{},{},{},{},{},{},{}
        EMP2,PPP2,TCA2,NADPH2,ATP2,PrdtFlux2,PrdtYield2,bio2,O2uptake2,Glcuptake2 = {},{},{},{},{},{},{},{},{},{}
        print(defaultObj,GSM)


        ##generates GPR dict for each gene in model (once instead of repeatedly)
        GPR_dict = generateOEGeneGPR(GSM,defaultModel)


        counter=0
        counterOEFail=0
        counterKOFail=0
        fail1 = 0
        fail2 = 0
        fail3 = 0
        fail4 = 0
        fail5 = 0
        fail6 = 0

    #####################altenrative for modeling######################################
        for dataPoint in FBATrainData.index:

            counter+=1
            modelKO = defaultModel.copy()

    ############### Determine if KO, GE instances, perform model simulation ############################

            MW = FBATrainData.loc[dataPoint].mw/1000
            print(MW)
            #Are there gene Knock-outs?
            if (FBATrainData.loc[dataPoint].number_genes_deleted!=0 and KO_option==1):
                #get gene KO data
                tempGenesModified = FBATrainData.genes_modified_updated[dataPoint].strip().split(';')
                tempKO = FBATrainData.gene_deletion[dataPoint].strip().split(';')

                #perform model KO
                modelKO = performGeneKOs(modelKO,GSM,tempKO,tempGenesModified)
                tempKOSol = modelKO.optimize()

                #Did the model produce an infeasible solution? Yes-revert to default soln
                if tempKOSol.status!='optimal':
                    print('geneKO growth failed')
                    sim_grw_flag=0
                    defaultKOBioObj = defaultFluxSol.objective_value
                    forPrdtModel = defaultModel.copy()
                else:
                    defaultKOBioObj = tempKOSol.objective_value
                    forPrdtModel = modelKO.copy()

                #Are there also gene overexpressions? (after KO)
                if (FBATrainData.loc[dataPoint].number_native_genes_overexp!=0 and OE_option==1):
                    #get gene overexpression data, heterologous genes
                    tempGenesOE = FBATrainData.loc[dataPoint].gene_overexpression.strip().split(';')
                    tempHetGenes = FBATrainData.loc[dataPoint].heterologous_gene.strip().split(';')

                    #perform model overexpression
                    modelKO,OE_fail,fail1,fail2,fail3,fail4,fail5,fail6 = performGeneOE(modelKO,GSM,tempGenesOE,tempGenesModified,tempHetGenes,tempKOSol,GPR_dict,epsilon[0],OE_fail,epsilon[1],epsilon[2],epsilon[5],fail1,fail2,fail3,fail4,fail5,fail6)

                    #perform OE FBA analysis with Biomass as objective
                    tempOESol = modelKO.optimize()

                    #Did the model produce an infeasible solution? Yes-Keep default (KO) soln
                    if tempOESol.status!='optimal':
                        counterKOFail+=1
                        print('KO & OE optmizing failed:',counterKOFail, 'OE genes', tempGenesModified[tempGenesOE==1],dataPoint)
                    else:
                        tempKOSol = tempOESol
                        #if infeasible, keep KO copy only, else take new model
                        forPrdtModel = modelKO.copy()

                    #UNUSED
                    tempOEdefaultBiomass = tempKOSol.objective_value
                    #use tempKO sol... then
                    #add the product/overexpression.....
                centCarbPrecursor = FBATrainData.loc[dataPoint].central_carbon_precursor.strip().split(';')

                cobra.manipulation.undelete_model_genes(defaultModel)

                dataPointFBASol = tempKOSol

            #There are no genetic Knock-outs, but are there gene overexpressions?
            elif (FBATrainData.loc[dataPoint].number_native_genes_overexp!=0 and OE_option==1):
                #Get gene overexpression data
                tempGenesModified = FBATrainData.loc[dataPoint].genes_modified_updated.strip().split(';')

                # tempGenesModified
                tempGenesOE = FBATrainData.loc[dataPoint].gene_overexpression.strip().split(';')
                tempHetGenes = FBATrainData.loc[dataPoint].heterologous_gene.strip().split(';')


                modelOE,OE_fail,fail1,fail2,fail3,fail4,fail5,fail6 = performGeneOE(modelKO,GSM,tempGenesOE,tempGenesModified,tempHetGenes,defaultFluxSol,GPR_dict,epsilon[0],OE_fail,epsilon[1],epsilon[2],epsilon[5],fail1,fail2,fail3,fail4,fail5,fail6)
                tempOESol = modelOE.optimize()

                #Did the model produce an infeasible solution? Yes-revert to default soln
                if tempOESol.status!='optimal':
                    counterOEFail+=1
                    print('OE optmizing failed:',counterOEFail,tempGenesModified[tempGenesOE==1], dataPoint)
                    tempOESol = defaultFluxSol
                    forPrdtModel = defaultModel.copy()
                else:
                    forPrdtModel = modelOE.copy()

                tempOEdefaultBiomass = tempOESol.objective_value

                dataPointFBASol = tempOESol

            #There are no genetic modifications
            else:
                dataPointFBASol = noGeneticMOSol = defaultFluxSol
                forPrdtModel = defaultModel.copy()

            if FBATrainData.loc[dataPoint].product_name == 'Riboflavin':
                isRbflv=1
            else:
                isRbflv=0

            if Product_option == 1:
                finalProdFluxSoln,counterProductFail,prod_fail = maximizeProduct(forPrdtModel,dataPointFBASol.objective_value,epsilon[3],epsilon[4],fbaModelMetaboliteDict,dataPoint,counterProductFail,GSM,prod_fail,isRbflv)
                EMP[dataPoint], PPP[dataPoint], TCA[dataPoint], NADPH[dataPoint], ATP[dataPoint], PrdtFlux[dataPoint],bio[dataPoint],O2uptake[dataPoint],Glcuptake[dataPoint] = FBAFeatureExtraction(finalProdFluxSoln,GSM)
                print(PrdtFlux[dataPoint])
            else:
                EMP[dataPoint], PPP[dataPoint], TCA[dataPoint], NADPH[dataPoint], ATP[dataPoint],PrdtFlux[dataPoint],bio[dataPoint],O2uptake[dataPoint],Glcuptake[dataPoint] = FBAFeatureExtraction(dataPointFBASol,GSM)
            PrdtYield[dataPoint] = PrdtFlux[dataPoint]*MW

            
            #Are there knock-outs to screen?
            tempOptKnock=[]
            optKnockModel = forPrdtModel.copy()
            additionalKnocks=0
            if (optKnockRxns.empty==False):
                for optKnockDataPoint in range(0,len(optKnockRxns)):
                    print(optKnockDataPoint,optKnockRxns.loc[optKnockDataPoint].rxns_deleted_updated_)
                    if Product_option == 1:
                        optKO = optKnockRxns.loc[optKnockDataPoint].rxns_deleted_updated_.strip().split(',')
                        # forPrdtModel2 = forPrdtModel.copy()
                        tempOptKnock=[]
                        tempKO2=[]
                        
                        #for each reaction, grab the associated gene to knock-out.
                        for junk in optKO:
                            try:
                                t = (optKnockModel.reactions.get_by_id(junk).gene_reaction_rule)
                                t2 = t.replace("(", "").replace(")", "").replace(" ", "").split('or')
                                t2=[ot.split('and', 1)[0] for ot in t2]
                                tempOptKnock = tempOptKnock + t2
                                tempKO2 = [1 for i in range(len(tempOptKnock))]
                                additionalKnocks = len(tempOptKnock)
                            except Exception as e:
                                t = junk
                                tempKO2 = t2 = t.replace("(", "").replace(")", "").replace(" ", "").split('or')
                                t2=[ot.split('and', 1)[0] for ot in t2]
                                tempOptKnock = tempOptKnock + t2
                                tempKO2 = [1 for i in range(len(tempOptKnock))]
                                additionalKnocks = len(tempOptKnock)
                    optKnockModel2 = optKnockModel.copy()
                    optKnockModel2 = performGeneKOs(optKnockModel2,GSM,tempKO2,tempOptKnock)
                    tempOptKnockSol = optKnockModel2.optimize()

                    if tempOptKnockSol.status!='optimal':
                        print('gtempOptKnockSol growth failed')
                        sim_grw_flag=0
                        defaultPrdtModelBioObj = tempOptKnockSol.objective_value
                        # optKnockModel2 = forPrdtModel.copy()
                    else:
                        defaultPrdtModelBioObj = tempOptKnockSol.objective_value
                        # forPrdtModel = modelKO.copy()

                    if Product_option == 1:
                        finalProdFluxSoln,counterProductFail,prod_fail = maximizeProduct(optKnockModel2,defaultPrdtModelBioObj,epsilon[3],epsilon[4],fbaModelMetaboliteDict,dataPoint,counterProductFail,GSM,prod_fail,isRbflv)
                        EMP2[optKnockDataPoint], PPP2[optKnockDataPoint], TCA2[optKnockDataPoint], NADPH2[optKnockDataPoint], ATP2[optKnockDataPoint], PrdtFlux2[optKnockDataPoint],bio2[optKnockDataPoint],O2uptake2[optKnockDataPoint],Glcuptake2[optKnockDataPoint] = FBAFeatureExtraction(finalProdFluxSoln,GSM)
                        print(PrdtFlux2[optKnockDataPoint])
                    else:
                        EMP2[optKnockDataPoint], PPP2[optKnockDataPoint], TCA2[optKnockDataPoint], NADPH2[optKnockDataPoint], ATP2[optKnockDataPoint], PrdtFlux2[optKnockDataPoint],bio2[optKnockDataPoint],O2uptake2[optKnockDataPoint],Glcuptake2[optKnockDataPoint] = FBAFeatureExtraction(dataPointFBASol,GSM)
                    PrdtYield2[optKnockDataPoint] = PrdtFlux2[optKnockDataPoint]*MW
                    print(PrdtYield2[optKnockDataPoint],PrdtFlux2[optKnockDataPoint])
                temp1 = FBATrainData.loc[dataPoint,'number_genes_deleted']
                temp1+=additionalKnocks
                FBATrainData.loc[dataPoint,'number_genes_deleted']=temp1
                temp2 = FBATrainData.loc[dataPoint].number_genes_mod+additionalKnocks
                FBATrainData.loc[dataPoint,'number_genes_mod']=temp2

                workingData2['EMP_'+GSM]=pd.concat([pd.Series(EMP),pd.Series(EMP2)],axis=0,ignore_index=True)
                workingData2['PPP_'+GSM]=pd.concat([pd.Series(PPP),pd.Series(PPP2)],axis=0,ignore_index=True)
                workingData2['TCA_'+GSM]=pd.concat([pd.Series(TCA),pd.Series(TCA2)],axis=0,ignore_index=True)
                workingData2['NADPH_'+GSM]=pd.concat([pd.Series(NADPH),pd.Series(NADPH2)],axis=0,ignore_index=True)
                workingData2['ATP_'+GSM]=pd.concat([pd.Series(ATP),pd.Series(ATP2)],axis=0,ignore_index=True)
                # workingData2['NADH_'+GSM]=pd.concat([pd.Series(NADH),pd.Series(NADH2)],axis=0,ignore_index=True)
                workingData2['PrdtFlux_'+GSM]=pd.concat([pd.Series(PrdtFlux),pd.Series(PrdtFlux2)],axis=0,ignore_index=True)
                workingData2['PrdtYield_'+GSM]=pd.concat([pd.Series(PrdtYield),pd.Series(PrdtYield2)],axis=0,ignore_index=True)
                workingData2['Biomass_'+GSM]=pd.concat([pd.Series(bio),pd.Series(bio2)],axis=0,ignore_index=True)
                workingData2['O2Uptake_'+GSM]=pd.concat([pd.Series(O2uptake),pd.Series(O2uptake2)],axis=0,ignore_index=True)
                workingData2['GlcUptake_'+GSM]=pd.concat([pd.Series(Glcuptake),pd.Series(Glcuptake2)],axis=0,ignore_index=True)


                test = pd.DataFrame()
                test = pd.DataFrame(FBATrainData.loc[dataPoint]).transpose()
                test = pd.concat([test]*(len(EMP2)+1), ignore_index=True)
                test = pd.concat([test,workingData2],axis=1)
                output = pd.concat([output,test],axis=0,ignore_index=True)

            else:
                test = pd.DataFrame()
                test = pd.DataFrame(FBATrainData.loc[dataPoint]).transpose()

                test['EMP_'+GSM]=pd.Series(EMP)
                test['PPP_'+GSM]=pd.Series(PPP)
                test['TCA_'+GSM]=pd.Series(TCA)
                test['NADPH_'+GSM]=pd.Series(NADPH)
                test['ATP_'+GSM]=pd.Series(ATP)
                # test['NADH_'+GSM]=pd.Series(NADH)
                test['PrdtFlux_'+GSM]=pd.Series(PrdtFlux)
                test['PrdtYield_'+GSM]=pd.Series(PrdtYield)
                test['Biomass_'+GSM]=pd.Series(bio)
                test['O2Uptake_'+GSM]=pd.Series(O2uptake)
                test['GlcUptake_'+GSM]=pd.Series(Glcuptake)
                # print(workingData2)

                output = pd.concat([output,test],axis=0,ignore_index=True)
                # output = test.copy()
        if (counter%50)==0:
            srint(counter)


    print(OE_fail,'OE failures')
    print(prod_fail,'Prod failures')
    print(fail1,fail2,fail3,fail4,fail5,fail6,'failure cases 1-6')
    return(output)
