# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 07:54:11 2022

@author: JohnDoe
"""

import pandas as pd
import time
from tdc import tdc
directory = "/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/"

#%%
### CBASS84 dataset
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

cbass84_ASV = cbass84_ASV.iloc[:, :-5].transpose()

values = cbass84_tax.index.to_series()
for column in cbass84_tax.columns:
    cbass84_tax[column] = cbass84_tax[column].fillna(values)
    
cbass84_tax['ASV'] = values

#%%
from tdc import tdc
start = time.time()
result_cbass84_tdc_wo_unicor = tdc(cbass84_ASV, cbass84_ED50, cbass84_tax, uc = False)
print(time.time()- start)
#%%
start = time.time()
result_cbass84_tdc = tdc(cbass84_ASV, cbass84_ED50, cbass84_tax)
print(time.time()- start)

#%%
cbass84_ASV = cbass84_ASV.transpose().merge(cbass84_tax, left_index=True, right_index=True)
cbass84_ASV_from_family = cbass84_ASV[cbass84_ASV["Family"].isin(("Microbacteriaceae", "Alteromonadaceae", "Marinobacteraceae", "Parvularculaceae", "Moraxellaceae", "ASV0100", "Burkholderiaceae", "Vibrionaceae", "Hyphomonadaceae", "ASV0096"))]
cbass84_ASV_from_family = cbass84_ASV_from_family.transpose()[0:28].astype('int32')
cbass84_from_family = cbass84_ED50.merge(cbass84_ASV_from_family, left_index=True, right_index=True)

#%%
import coracle
start = time.time()
cbass84_ASV_family = cbass84_ASV.groupby(["Family"]).sum()
cbass84_family = cbass84_ED50.merge(cbass84_ASV_family.transpose(), left_index=True, right_index=True)

x = cbass84_family.iloc[:,1:]
y = cbass84_family[['ED50']]


result_f84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)
#%%
start = time.time()
x = cbass84_from_family.iloc[:,1:]
y = cbass84_from_family[['ED50']]


result_ff84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)



#%%
result_f84.to_csv("C:/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/result_tdc_comparison_cbass84/cf84.csv", float_format="%.2g")
result_ff84.to_csv("C:/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/result_tdc_comparison_cbass84/cff84.csv", float_format="%.2g")
result_cbass84_tdc.to_csv("C:/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/result_tdc_comparison_cbass84/tdc84.csv", float_format="%.2g")
result_cbass84_tdc_wo_unicor.to_csv("C:/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/result_tdc_comparison_cbass84/tdc84_wo_unicor.csv", float_format="%.2g")



#%%
"""
#%%
result_cbass84_tdc = {}
for i in [0.2, 0.3, 0.4]:
    for j in [0.2, 0.3, 0.4]:
        result_cbass84_tdc[str([i,j])] = tdc(cbass84_ASV, cbass84_ED50, cbass84_tax, threshold=i, uc=True, uc_threshold=j)

#%%
a = result_cbass84_tdc[str([i,j])].iloc[0, 1:9].mean()

#%%
for i in [0.2, 0.3, 0.4]:
    for j in [0.2, 0.3, 0.4]:
        result_cbass84_tdc[str([i,j])] = result_cbass84_tdc[str([i,j])].iloc[0, 1:9].mean()
#%%
for i in [0.2, 0.3]:
    for j in [0.15, 0.25, 0.35, 0.45]:
        result_cbass84_tdc[str([i,j])] = tdc(cbass84_ASV, cbass84_ED50, cbass84_tax, threshold=i, uc=True, uc_threshold=j).iloc[0, 1:9].mean()
        
#%%
result_cbass84_tdc_1 = {}
for i in [0.16, 0.18, 0.2, 0.22, 0.24]:
    for j in [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]:
        result_cbass84_tdc_1[str([i,j])] = tdc(cbass84_ASV, cbass84_ED50, cbass84_tax, threshold=i, uc=True, uc_threshold=j).iloc[0, 1:9].mean()

"""


#%%
######### CBASS dataset




### CBASS dataset
cbass_meta = pd.read_csv(directory + "CBASS_gradient_metafile.txt", sep= "\s+")
cbass_ASV = pd.read_csv(directory + "ASV_table_QCfiltered.txt", sep= "\s+")
cbass_ED50 = pd.read_csv(directory + "ED50_per_colony.txt", sep= "\s+", index_col=(0))
cbass_tax = cbass_ASV.iloc[:, -7:-1] #get taxonomic information
cbass_ASV = cbass_ASV.iloc[:, :-8]



#cbass_ASV = cbass_ASV.transpose()
cbass_ASV = cbass_ASV.loc[:, cbass_ASV.columns.str.contains('30')]
cbass_ASV.columns = cbass_ASV.columns.to_series().astype(str).str.replace(r'_30.+','')


#change meta
cbass_meta.set_index(["Sample"], inplace = True)
cbass_meta = cbass_meta.loc[cbass_meta.index.str.contains('30'), :]
cbass_meta.index = cbass_meta.index.to_series().astype(str).str.replace(r'_30.+','')
cbass_meta.drop(["Temperature", "Colony"], axis = 1,  inplace = True)


#merge
cbass_ASV = cbass_ASV.transpose()
cbass = cbass_ED50.merge(cbass_meta, left_index=True, right_index=True)
cbass = cbass.merge(cbass_ASV, left_index=True, right_index=True)

#split by species
acropora = cbass.loc[cbass.Species=='Acropora']
pocillopora = cbass.loc[cbass.Species=='Pocillopora']
porites = cbass.loc[cbass.Species=='Porites']
stylophora = cbass.loc[cbass.Species=='Stylophora']

acropora.drop(["Site", "Species"], axis = 1,  inplace = True)
pocillopora.drop(["Site", "Species"], axis = 1,  inplace = True)
porites.drop(["Site", "Species"], axis = 1,  inplace = True)
stylophora.drop(["Site", "Species"], axis = 1,  inplace = True)

acropora_ED50 = pd.DataFrame(acropora["ED50"])
acropora_ASV = acropora.loc[:, (acropora != 0).any(axis=0)]
acropora_ASV = acropora.iloc[:, 1:]


pocillopora_ED50 = pd.DataFrame(pocillopora["ED50"])
pocillopora_ASV = pocillopora.loc[:, (pocillopora != 0).any(axis=0)]
pocillopora_ASV = pocillopora.iloc[:, 1:]

porites_ED50 = pd.DataFrame(porites["ED50"])
porites_ASV = porites.loc[:, (porites != 0).any(axis=0)]
porites_ASV = porites.iloc[:, 1:]

stylophora_ED50 = pd.DataFrame(stylophora["ED50"])
stylophora_ASV = stylophora.loc[:, (stylophora != 0).any(axis=0)]
stylophora_ASV = stylophora.iloc[:, 1:]


#%%
#prepare tax:
values = cbass_tax.index.to_series()
for column in cbass_tax.columns:
    cbass_tax[column] = cbass_tax[column].fillna(values)

"""
cbass_tax.to_csv("C:/Users/JohnDoe/Desktop/data/webserver test/tax.csv", float_format="%.4g")
acropora_ED50.to_csv("C:/Users/JohnDoe/Desktop/data/webserver test/target.csv", float_format="%.4g")
acropora_ASV.to_csv("C:/Users/JohnDoe/Desktop/data/webserver test/ASV.csv", float_format="%.4g")
### test on cbass data with small threshold
"""





#%%
start = time.time()
result_acropora = tdc(acropora_ASV, acropora_ED50, cbass_tax)
print(time.time()- start)
#%%
start = time.time()
result_pocillopora = tdc(pocillopora_ASV, pocillopora_ED50, cbass_tax)
print(time.time()- start)
#%%
start = time.time()
result_porites = tdc(porites_ASV, porites_ED50, cbass_tax)
print(time.time()- start)
#%%
start = time.time()
result_stylophora = tdc(stylophora_ASV, stylophora_ED50, cbass_tax)
print(time.time()- start)




#%%
from td_coracle import td_coracle_wo_uc

### test on cbass data with small threshold
#%%
start = time.time()
result_acropora_wo_uc = tdc(acropora_ASV, acropora_ED50, cbass_tax, uc=False)
print(time.time()- start)
#%%
start = time.time()
result_pocillopora_wo_uc = tdc(pocillopora_ASV, pocillopora_ED50, cbass_tax, uc=False)
print(time.time()- start)
#%%
start = time.time()
result_porites_wo_uc = tdc(porites_ASV, porites_ED50, cbass_tax, uc=False)
print(time.time()- start)
#%%
start = time.time()
result_stylophora_wo_uc = tdc(stylophora_ASV, stylophora_ED50, cbass_tax, uc=False)
print(time.time()- start)
#%%
"""
from td_coracle import td_coracle_wo_uc

### test on cbass data with small threshold
#%%
start = time.time()
result_acropora_wo_uc = td_coracle_wo_uc(acropora_ASV, acropora_ED50, cbass_tax, threshold=0.2)
print(time.time()- start)
#%%
start = time.time()
result_pocillopora_wo_uc = td_coracle_wo_uc(pocillopora_ASV, pocillopora_ED50, cbass_tax, threshold=0.2)
print(time.time()- start)
#%%
start = time.time()
result_porites_wo_uc = td_coracle_wo_uc(porites_ASV, porites_ED50, cbass_tax, threshold=0.2)
print(time.time()- start)
#%%
start = time.time()
result_stylophora_wo_uc = td_coracle_wo_uc(stylophora_ASV, stylophora_ED50, cbass_tax, threshold=0.2)
print(time.time()- start)
"""

#%%
"""
result_acropora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_acropora.csv", float_format="%.4g")
result_pocillopora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_pocillopora.csv", float_format="%.4g")
result_porites.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_porites.csv", float_format="%.4g")
result_stylophora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_stylophora.csv", float_format="%.4g")
result_acropora_wo_uc.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_acropora_wo_uc.csv", float_format="%.4g")
result_pocillopora_wo_uc.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_pocillopora_wo_uc.csv", float_format="%.4g")
result_porites_wo_uc.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_porites_wo_uc.csv", float_format="%.4g")
result_stylophora_wo_uc.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_stylophora_wo_uc.csv", float_format="%.4g")
"""