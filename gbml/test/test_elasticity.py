#!/usr/bin/env python

# Test script for gbml elasticity (bulk and shear moduli) predictions

from gbml import elasticity

# Use a Mock query engine to return the data
class MockQE:
    _materials = {
        "mp-26": : {},
        "mp-10003": {},
        "mp-10010": {},
        "mp-10015": {},
        "mp-10021": {},
        "mp-10018": {},
        "mp-19306": {}
    }
    
    def query(mid_list):
        return [ self._materials.get(mid, None) for mid in mid_list ]




mpID = "mp-10003"
(k_value, g_value, caveat_str) = elasticity.predict_k_g(mpID, query_engine=MockQE())
assert (k_value, g_value, caveat_str) == (expected_k_value, expected_g_value, expected_caveat_str)

mpID_list = ["mp-10003","mp-10010","mp-10015","mp-10021","mp-26","mp-10018","mp-19306"]
(matid_list, k_list, g_list, caveat_list) = elasticity.predict_k_g_list(mpID_list, query_engine=MockQE())
assert (matid_list, k_list, g_list, caveat_list) == (expected_matid_list, expected_k_list, expected_g_list, expected_caveat_list)
