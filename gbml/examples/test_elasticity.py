#!/usr/bin/env python

# Test script for gbml elasticity (bulk and shear moduli) predictions

from gbml import elasticity

elasticity.API_KEY = "YOUR_API_KEY_HERE"

mpID = "mp-10003"
(k_value, g_value, caveat_str) = elasticity.predict_k_g(mpID)
print (k_value, g_value, caveat_str)

mpID_list = ["mp-10003","mp-10010","mp-10015","mp-10021","mp-26","mp-10018","mp-19306"]
(matid_list, k_list, g_list, caveat_list) = elasticity.predict_k_g_list(mpID_list)
print (matid_list, k_list, g_list, caveat_list)
