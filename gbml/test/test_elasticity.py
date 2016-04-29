#!/usr/bin/env python

# Test script for gbml elasticity (bulk and shear moduli) predictions

from gbml import elasticity

# Use a Mock query engine to return the data
class MockQE:
    _materials = {
        "mp-26": {},
        "mp-10003": {
            u'energy_per_atom': -9.174497691666668,
            u'is_hubbard': False,
            u'material_id': u'mp-10003',
            u'nsites': 12,
            u'pretty_formula': u'Nb4CoSi',
            u'volume': 194.5128160886403
        },
        "mp-10010": {},
        "mp-10015": {},
        "mp-10021": {},
        "mp-10018": {},
        "mp-19306": {}
    }
    
    def query(criteria, properties):
        mid_list = criteria["task_id"]["$in"]
        return [ self._materials.get(mid, None) for mid in mid_list ]


(expected_k_value, expected_g_value, expected_caveat_str) = (175.30512291338607, 84.49987188140813, '')

mpID = "mp-10003"
(k_value, g_value, caveat_str) = elasticity.predict_k_g(mpID, query_engine=MockQE())
assert (k_value, g_value, caveat_str) == (expected_k_value, expected_g_value, expected_caveat_str)

# mpID_list = ["mp-10003","mp-10010","mp-10015","mp-10021","mp-26","mp-10018","mp-19306"]
# (matid_list, k_list, g_list, caveat_list) = elasticity.predict_k_g_list(mpID_list, query_engine=MockQE())
# assert (matid_list, k_list, g_list, caveat_list) == (expected_matid_list, expected_k_list, expected_g_list, expected_caveat_list)
