# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 09:23:45 2024

@author: JohnDoe2Go
"""

import pandas as pd
from unicoracle import unicoracle


directory = "C:/Users/JohnDoe2Go/Downloads/" #use your path to the variables

### CBASS84 dataset
ASV = pd.read_csv(directory + "cbass.csv", index_col=0) #read in cbass ASVs
tax = pd.read_csv(directory + "cbass_tax.csv", index_col=0) #read in cbass taxonomic information
tax['ASV'] = tax.index
### 2. Split Target Variable
y = ASV["ED50"].to_frame()
### 3. Combine ASV with tax
x = ASV.iloc[:,3:]


#%% Fill empty fields with the corresponding asv number
tax_asvfill = tax.copy()

# Fill missing values with the corresponding `asv` value
for col in tax_asvfill.columns[:]:
    tax_asvfill[col] = tax_asvfill[col].replace('', pd.NA)  # Replace empty strings with NaN
    tax_asvfill[col] = tax_asvfill[col].fillna(tax_asvfill.index.to_series())
#%%
import unittest
from unicoracle import unicoracle

#unit testing
class TestWrongInput(unittest.TestCase):
    def test_wrong_input_asv(self): #tests if TypeError is thrown if asv file is not a pandas dataframe
        self.assertRaises(TypeError, unicoracle, 1, y, tax)
    def test_wrong_input_target(self): #tests if TypeError is thrown if target file is not a pandas dataframe, or series
        self.assertRaises(TypeError, unicoracle, x, 1, tax)
    def test_wrong_input_tax(self): #tests if TypeError is thrown if tax file is not a pandas dataframe
        self.assertRaises(TypeError, unicoracle, x, y, 1)
    def test_wrong_input_nfeatures(self): #tests if TypeError is thrown for non-float/integer input for threshold
        self.assertRaises(ValueError, unicoracle, x, y, tax, "1")
    def test_wrong_input_uc(self): #tests if TypeError is thrown for non boolean input for uc
        self.assertRaises(TypeError, unicoracle, x, y, tax, uc=1)
    def test_wrong_input_uc_threshold(self): #tests if TypeError is thrown for non-float/integer input for uc_threshold
        self.assertRaises(TypeError, unicoracle, x, y, tax, uc_threshold="1")
    def test_wrong_nfeatures1(self): #tests ValueError is thrown if threshold not between 0 and 1
        self.assertRaises(ValueError, unicoracle, x, y, tax, n_features=1)
    def test_wrong_nfeatures2(self): #tests ValueError is thrown if threshold not between 0 and 1
        self.assertRaises(ValueError, unicoracle, x, y, tax, n_features=-2.2)
    def test_wrong_uc_threshold1(self): #tests ValueError is thrown if threshold not between 0 and 1
        self.assertRaises(ValueError, unicoracle, x, y, tax, uc_threshold=1.5)
    def test_wrong_uc_threshold2(self): #tests ValueError is thrown if threshold not between 0 and 1
        self.assertRaises(ValueError, unicoracle, x, y, tax, uc_threshold=-1.1)
    def test_wrong_dimensions1(self): #tests ValueError is thrown if asv is not two dimensional
        self.assertRaises(ValueError, unicoracle, y, y, tax)
    def test_wrong_dimensions2(self): #tests ValueError is thrown if tax is not two dimesnional
        self.assertRaises(ValueError, unicoracle, x, y, y)
    def test_wrong_dimensions3(self): #tests ValueError is thrown if asv and tax don't match in their dimensions (number of features)
        self.assertRaises(ValueError, unicoracle, x, y, x)
    def test_wrong_dimensions4(self): #tests ValueError is thrown if asv and target don't match in their dimensions (number of samples)
        self.assertRaises(ValueError, unicoracle, x, tax, tax)
    def test_wrong_input_uc_top_k_type_str(self):  # non-int uc_top_k
        self.assertRaises(TypeError, unicoracle, x, y, tax, uc_top_k="10")
    def test_wrong_input_uc_top_k_type_float(self):  # non-int uc_top_k
        self.assertRaises(TypeError, unicoracle, x, y, tax, uc_top_k=10.5)
    def test_wrong_input_uc_top_k_value_zero(self):  # uc_top_k must be >= 1
        self.assertRaises(ValueError, unicoracle, x, y, tax, uc_top_k=0)
    def test_wrong_input_uc_top_k_value_negative(self):  # uc_top_k must be >= 1
        self.assertRaises(ValueError, unicoracle, x, y, tax, uc_top_k=-5)
    def test_wrong_input_uc_method_value(self):  # only 'pearson' or 'spearman'
        self.assertRaises(ValueError, unicoracle, x, y, tax, uc_method='kendall')
    def test_wrong_input_uc_transformation_value(self):  # only 'raw', 'relative_abundance', 'clr'
        self.assertRaises(ValueError, unicoracle, x, y, tax, uc_transformation='log')
    def test_wrong_input_seed_type_str(self):  # seed must be int or None
        self.assertRaises(TypeError, unicoracle, x, y, tax, seed='42')
    def test_wrong_input_seed_type_float(self):  # seed must be int or None
        self.assertRaises(TypeError, unicoracle, x, y, tax, seed=3.14)
    def test_wrong_input_coracle_kwargs_type(self):  # coracle_kwargs must be dict or None
        self.assertRaises(TypeError, unicoracle, x, y, tax, coracle_kwargs=['rf', 'ridge'])
    def test_xor_uc_threshold_and_top_k_both_set(self):  # XOR: cannot set both
        self.assertRaises(ValueError, unicoracle, x, y, tax, uc_threshold=0.2, uc_top_k=50)
    
#test runner
if __name__ == "__main__":
    unittest.main()

#%%
result_cbass84_unicoracle_wo_unicor = unicoracle(x, y, tax_asvfill, n_features = 100, uc = False)
#%%
import time
from unicoracle import unicoracle 

start_time = time.time()  # Start timer
result_cbass84_unicoracle = unicoracle(x, y, tax_asvfill, n_features = 10)  # Run your function
total_time = time.time()  - start_time# End timer



