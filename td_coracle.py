# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 07:15:57 2022

@author: JohnDoe
"""

import pandas as pd
from hierarchical_filtering import fs
import coracle

directory = "/Users/JohnDoe/Desktop/PhD/10_data/datasets/"

#%%


### CBASS84 dataset
directory = "/Users/JohnDoe/Desktop/data/datasets/"
#cbass84_meta = pd.read_csv("cbass84_metadata.txt", sep= "[ \t \t]+", error_bad_lines=False) #, usecols=(["Sample", "Site", "Reef", "Species"]))
cbass84_ASV = pd.read_csv(directory + "cbass84_ASVs.txt", sep= "\s+")
cbass84_ED50 = pd.read_csv(directory + "cbass84_ED50s.txt", sep= "\s+")
cbass84_tax = cbass84_ASV.iloc[:, -5:] #get taxonomic information
cbass84_ASV = cbass84_ASV.iloc[:, :-5] #get rest
cbass84_ASV = cbass84_ASV.transpose()
cbass84_ASV = cbass84_ASV.loc[:, (cbass84_ASV != 0).any(axis=0)]
cbass84_ASV = cbass84_ASV.transpose()
cbass84_ASV = cbass84_ASV.merge(cbass84_tax, left_index=True, right_index=True)
cbass84_tax = cbass84_ASV.iloc[:, -5:] #get taxonomic information

#change ASV
cbass84_ASV_genus = cbass84_ASV.groupby(["Genus"]).sum()
cbass84_ASV_genus = cbass84_ASV_genus.transpose()
cbass84_ASV_family = cbass84_ASV.groupby(["Family"]).sum()
cbass84_ASV_family = cbass84_ASV_family.transpose()
cbass84_ASV_order = cbass84_ASV.groupby(["Order"]).sum()
cbass84_ASV_order = cbass84_ASV_order.transpose()
cbass84_ASV_class = cbass84_ASV.groupby(["Class"]).sum()
cbass84_ASV_class = cbass84_ASV_class.transpose()


#change ED50
cbass84_ED50["Site"] = None
cbass84_ED50["Site"][cbass84_ED50["Sample"].str.contains("AF")] = "Al Fahal (AF)"
cbass84_ED50["Site"][cbass84_ED50["Sample"].str.contains("ExT")] = "Tahala (ExT)"
cbass84_ED50["Site"][cbass84_ED50["Sample"].str.contains("PrT")] = "Tahala (PrT)"
cbass84_ED50["Site"][cbass84_ED50["Sample"].str.contains("ICN")] = "Interuniversity Institute for Marine Science (IUI)"
cbass84_ED50["Species"] = "Stylophora pistillata"
cbass84_ED50.set_index(["Sample"], inplace = True)
#specific for correlation analysis:
cbass84_ED50 = pd.DataFrame(cbass84_ED50["ED50"])

#get filtered ASV
ASV = cbass84_ASV.iloc[:, :-5]
cbass84_filtered_ASV = fs(ASV, cbass84_ED50, cbass84_tax)



#change ASV
cbass84_ASV_genus = cbass84_ASV.groupby(["Genus"]).sum()
cbass84_ASV_genus = cbass84_ASV_genus.transpose()
cbass84_ASV_family = cbass84_ASV.groupby(["Family"]).sum()
cbass84_ASV_family = cbass84_ASV_family.transpose()
cbass84_ASV_order = cbass84_ASV.groupby(["Order"]).sum()
cbass84_ASV_order = cbass84_ASV_order.transpose()
cbass84_ASV_class = cbass84_ASV.groupby(["Class"]).sum()
cbass84_ASV_class = cbass84_ASV_class.transpose()
cbass84_ASV_phylum = cbass84_ASV.groupby(["Phylum"]).sum()
cbass84_ASV_phylum = cbass84_ASV_phylum.transpose()

cbass84_ASV_genus_f = cbass84_filtered_ASV.groupby(["Genus"]).sum()
cbass84_ASV_genus_f = cbass84_ASV_genus_f.transpose()
cbass84_ASV_family_f = cbass84_filtered_ASV.groupby(["Family"]).sum()
cbass84_ASV_family_f = cbass84_ASV_family_f.transpose()
cbass84_ASV_order_f = cbass84_filtered_ASV.groupby(["Order"]).sum()
cbass84_ASV_order_f = cbass84_ASV_order_f.transpose()
cbass84_ASV_class_f = cbass84_filtered_ASV.groupby(["Class"]).sum()
cbass84_ASV_class_f = cbass84_ASV_class_f.transpose()
cbass84_ASV_phylum_f = cbass84_filtered_ASV.groupby(["Phylum"]).sum()
cbass84_ASV_phylum_f = cbass84_ASV_phylum_f.transpose()



#merge
cbass84_genus = cbass84_ED50.merge(cbass84_ASV_genus, left_index=True, right_index=True)
cbass84_family = cbass84_ED50.merge(cbass84_ASV_family, left_index=True, right_index=True)
cbass84_order = cbass84_ED50.merge(cbass84_ASV_order, left_index=True, right_index=True)
cbass84_class = cbass84_ED50.merge(cbass84_ASV_class, left_index=True, right_index=True)
cbass84_phylum = cbass84_ED50.merge(cbass84_ASV_phylum, left_index=True, right_index=True)

cbass84_genus_f = cbass84_ED50.merge(cbass84_ASV_genus_f, left_index=True, right_index=True)
cbass84_family_f = cbass84_ED50.merge(cbass84_ASV_family_f, left_index=True, right_index=True)
cbass84_order_f = cbass84_ED50.merge(cbass84_ASV_order_f, left_index=True, right_index=True)
cbass84_class_f = cbass84_ED50.merge(cbass84_ASV_class_f, left_index=True, right_index=True)
cbass84_phylum_f = cbass84_ED50.merge(cbass84_ASV_phylum_f, left_index=True, right_index=True)


#%%

"""
x = cbass84_phylum.iloc[:,3:]
y = cbass84_phylum[['ED50']]

start = time.time()
result_p84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = cbass84_phylum_f.iloc[:,3:]
y_f = cbass84_phylum_f[['ED50']]

start = time.time()
result_p84_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)
"""

#%%

def td_coracle(ASV, target, tax, threshold = 0.2, uc = True, uc_threshold = 0.15):
    """
    Algorithm to analyse hierarchical ASV datasets and identify bacteria/ASV that are associated with a continuous target variable
    First it uses the UniCor algorithm to propagate the most unique and valuable features to the next higher level (preserve information). Then a top-down skimming approach uses Coracle to select the [threshold] top partition consecutively at each level. The lowest level is once again anylized with Coracle and given back as the final result.

    Parameters
    ----------
    ASV : pd.dataframe
        ASV dataset
    target_var : pd.dataframe
        Continuous target variable
    tax : pd.dataframe
        Taxonomic hierarchy
    threshold : bool, optional
        Propagation threshold. Has to be between 0 and 1. Gives the proportion of the best performing groups that are propagated to the next lower level. Default is 0.2. In general you want this value to be as small as possible for runtime efficiency reasons.
    uc_threshold : bool, optional
        UniCor metric threshold. Can only be between 0 and 1. Optimal values depend on the dataset. The default is 0.15.

    Returns
    -------
    resultO : Final Coracle result on the lowest level

    """
    
    ### check input types
    ###########################################################################
    #check type
    if not isinstance(ASV, (pd.DataFrame)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for ASV")
    if not isinstance(target, (pd.DataFrame)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for target")
    if not isinstance(tax, (pd.DataFrame, pd.Series)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for tax")
    if not isinstance(threshold, (float)):
        raise TypeError("TypeError exception thrown. Expected float for threshold")
    if not isinstance(uc, (bool)):
        raise TypeError("TypeError exception thrown. Expected bool for uc")
    if not isinstance(uc_threshold, (float)):
        raise TypeError("TypeError exception thrown. Expected float uc_threshold")
    #check values    
    if threshold > 1 or threshold < 0:
        raise TypeError("ValueError exception thrown. threshold is expected to be a float between 0 and 1")
    if uc_threshold > 1 or uc_threshold < 0:
        raise TypeError("ValueError exception thrown. uc_threshold is expected to be a float between 0 and 1")
    #check dimensions
    if ASV.ndim != 2 or tax.ndim != 2:  
        raise ValueError("ValueError exception thrown. Expected ASV and tax to have two dimensions")
    if ASV.shape[0] != tax.shape[0]:
        raise ValueError("ValueError exception thrown. Expected ASV and tax to have the same number of samples")
    if ASV.shape[1] != (target.shape[0] + tax.shape[1]):
        raise ValueError("ValueError exception thrown. Expected ASV and target to have the same number of samples")
    
    
    
    
    ### check order/hierarchy of taxonomy
    ###########################################################################
    tax_levels = list(tax.columns) #get tax level names
    ASV_tax = ASV.merge(tax, how="left", left_index=True, right_index=True) #merge to compare size
    size = {} #dict to check sizes
    last_size = 0 #helper variable to save the last size
    count = 0 #helper variable
    
    for i in tax_levels: #for every taxonomic level
        print(i)
        size[i] = len(ASV_tax[i].unique()) #get number of unique entries in that level
        print(size[i])
        
        if size[i] < last_size: #if number of unique entries in current level is smaller than in the last level
            count += 1 #increase count
        last_size = size[i] #set current size to future last size
    
    #check if hierarchical order is fulfilled
    if count == len(tax_levels)-1: #check if hierarchical order is ascending from left to right(-1 because first comparison is with count=0)
        tax_levels.reverse() #if that is the case, reverse the order
        print("taxonomic order has been reversed") #print a note
    elif count > 0 and count < len(tax_levels)-1: #if order is not unambiguous: (-1 because first comparison is with count=0)
        print("Warning: unclear hierarchy as number of unique entries is neither clearly ascending nor descending! The order will be unchanged and used as if descending from left to right") #print a warning
    ###########################################################################
        
        
    ### UniCor to filter ASV
    ###########################################################################
    if uc == True: #if UniCor is wanted
        filtered_ASV = fs(ASV, target, tax, threshold=uc_threshold)
    else: #if not
        filtered_ASV = ASV.merge(tax, how="left", left_index=True, right_index=True)
    ###########################################################################
    
    ### Top-Down Coracle propagation
    ###########################################################################
    for i in range(len(tax_levels)): #for every taxonomic level
        selected = [] #list of selected groups to propagate
        highest_level = filtered_ASV.groupby([tax_levels[i]]).sum().transpose() #start at the highest level
        highest_level = target.merge(highest_level, left_index=True, right_index=True) #merge with target variable
        x = highest_level.iloc[:,1:] 
        y = highest_level.iloc[:,0].to_frame()
        
        resultO = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4)) #run coracle
        result = resultO.iloc[3:,0] #select only the groups and their respective score
        for j in range(int(result.size * threshold)+1): #for every group within the threshold partition (always round up to prevent zero features)
            if result[j] > 0: #if score is not null
                #print(result.index[j])
                selected.append(result.index[j]) #add to selected-list
            else: #else print a warning that no more features will be added
                print("Warning: as the score reached zero, threshold partition will not be fulfilled but will include the top", ((100*j)/(result.size * threshold)), "percent")
                break
        
        
        filtered_ASV = filtered_ASV[filtered_ASV[tax_levels[i]].isin(selected)] #propagate them to the next level
        
        
        
    ###########################################################################
    
    return resultO #return final coracle result
    

def td_coracle_wo_uc(ASV, target, tax, threshold=0.5):
    
    
    
    tax_levels = list(tax.columns) #get tax level names
    ASV_tax = ASV.merge(tax, how="left", left_index=True, right_index=True) #merge to compare size
    
    
    ### check order/hierarchy of taxonomy
    ###########################################################################
    size = {} #dict to check sizes
    last_size = 0 #helper variable to save the last size
    count = 0 #helper variable
    for i in tax_levels:
        print(i)
        size[i] = len(ASV_tax[i].unique())
        print(size[i])
        if last_size:
            if size[i] < last_size: 
                count += 1
        last_size = size[i]
    
    if count == len(tax_levels)-1: #-1 because first comparison is with count=0
        tax_levels.reverse()
        print("taxonomic order has been reversed")
    elif count > 0 and count < len(tax_levels)-1: #-1 because first comparison is with count=0
        print("Warning: unclear hierarchy as number of unique entries is neither clearly ascending nor descending! The order will be unchanged and used as if descending from left to right")
    ###########################################################################
        
        
    ### UniCor to filter ASV
    ###########################################################################
    filtered_ASV = ASV.merge(tax, how="left", left_index=True, right_index=True)
    ###########################################################################
    
    ### Top-Down Coracle propagation
    ###########################################################################
    for i in range(len(tax_levels)):
        selected = []
        highest_level = filtered_ASV.groupby([tax_levels[i]]).sum().transpose()
        highest_level = target.merge(highest_level, left_index=True, right_index=True)
        x = highest_level.iloc[:,1:]
        y = highest_level.iloc[:,0].to_frame()
        
        resultO = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
        result = resultO.iloc[3:,0]
        for j in range(int(result.size * threshold)+1): #always round up to prevent zero features
            if result[j] > 0:
                #print(result.index[j])
                selected.append(result.index[j])
            else:
                print("Warning: as the score reached zero, threshold partition will not be fulfilled but will include the top", ((100*j)/(result.size * threshold)), "percent")
                break
        
        #if i <= range(len(tax_levels)-1):
        filtered_ASV = filtered_ASV[filtered_ASV[tax_levels[i]].isin(selected)]
        
        
        
    ###########################################################################
    
    return resultO

#%%

#test 
#start = time.time()
test04 = td_coracle(cbass84_ASV, cbass84_ED50, cbass84_tax, threshold=0.4)
#print(time.time()- start)
#%%
#filtered_ASV = filtered_ASV[filtered_ASV[i].isin(selected)]



#test02 = td_coracle(ASV, cbass84_ED50, cbass84_tax[cbass84_tax.columns[::-1]], threshold=0.5)
#test03 = td_coracle(ASV, cbass84_ED50, cbass84_tax[cbass84_tax.columns[::-1]].merge(cbass84_tax, how="left", left_index=True, right_index=True), threshold=0.5)

#%%
"""
#optimize thresholds
threshold = {}
runtime = {}
for i in np.arange(0.2, 0.85, 0.05): 
    
    start = time.time()
    result = td_coracle(ASV, cbass84_ED50, cbass84_tax, threshold=i)
    duration = time.time()- start
    print(i)
    print(duration)
    threshold[i] = result
    runtime[i] = duration
    """