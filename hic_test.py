# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 09:23:45 2024

@author: JohnDoe2Go
"""

#%%
import pandas as pd
from hicoracle import hic


directory = "C:/Users/JohnDoe2Go/Downloads/" #use your path to the variables

### CBASS84 dataset
ASV = pd.read_csv(directory + "cbass.csv", index_col=0) #read in cbass ASVs
tax = pd.read_csv(directory + "cbass_tax.csv", index_col=0) #read in cbass taxonomic information
### 2. Split Target Variable
y = ASV["ED50"].to_frame()
### 3. Combine ASV with tax
x = ASV.iloc[:,3:]
#%%
import unittest


#unit testing
class TestWrongInput(unittest.TestCase):
    def test_wrong_input_asv(self): #tests if TypeError is thrown if asv file is not a pandas dataframe
        self.assertRaises(TypeError, hic, 1, y, tax)
    def test_wrong_input_target(self): #tests if TypeError is thrown if target file is not a pandas dataframe, or series
        self.assertRaises(TypeError, hic, x, 1, tax)
    def test_wrong_input_tax(self): #tests if TypeError is thrown if tax file is not a pandas dataframe
        self.assertRaises(TypeError, hic, x, y, 1)
    def test_wrong_input_nfeatures(self): #tests if TypeError is thrown for non-float/integer input for threshold
        self.assertRaises(ValueError, hic, x, y, tax, "1")
    def test_wrong_input_uc(self): #tests if TypeError is thrown for non boolean input for uc
        self.assertRaises(TypeError, hic, x, y, tax, uc=1)
    def test_wrong_input_uc_threshold(self): #tests if TypeError is thrown for non-float/integer input for uc_threshold
        self.assertRaises(TypeError, hic, x, y, tax, uc_threshold="1")
    def test_wrong_nfeatures1(self): #tests ValueError is thrown if threshold not between 0 and 1
        self.assertRaises(ValueError, hic, x, y, tax, n_features=1)
    def test_wrong_nfeatures2(self): #tests ValueError is thrown if threshold not between 0 and 1
        self.assertRaises(ValueError, hic, x, y, tax, n_features=-2.2)
    def test_wrong_uc_threshold1(self): #tests ValueError is thrown if threshold not between 0 and 1
        self.assertRaises(ValueError, hic, x, y, tax, uc_threshold=1.5)
    def test_wrong_uc_threshold2(self): #tests ValueError is thrown if threshold not between 0 and 1
        self.assertRaises(ValueError, hic, x, y, tax, uc_threshold=-1.1)
    def test_wrong_dimensions1(self): #tests ValueError is thrown if asv is not two dimensional
        self.assertRaises(ValueError, hic, y, y, tax)
    def test_wrong_dimensions2(self): #tests ValueError is thrown if tax is not two dimesnional
        self.assertRaises(ValueError, hic, x, y, y)
    def test_wrong_dimensions3(self): #tests ValueError is thrown if asv and tax don't match in their dimensions (number of features)
        self.assertRaises(ValueError, hic, x, y, x)
    def test_wrong_dimensions4(self): #tests ValueError is thrown if asv and target don't match in their dimensions (number of samples)
        self.assertRaises(ValueError, hic, x, tax, tax)

#test runner
if __name__ == "__main__":
    unittest.main()

#%%
result_cbass84_hic_wo_unicor = hic(x, y, tax, n_features = 2, uc = False)
#%%
result_cbass84_hic = hic(x, y, tax, n_features = 100)