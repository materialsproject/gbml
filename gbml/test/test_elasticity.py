#!/usr/bin/env python

# Test script for gbml elasticity (bulk and shear moduli) predictions

from gbml import elasticity
import unittest 

from configparser import SafeConfigParser

import os

try:
    parser = ConfigParser.SafeConfigParser()
    cur_dir = os.path.dirname(__file__)
    parser.read('{}/test.ini'.format(cur_dir))
    api_key = parser.get('mprester', 'api_key')
except Exception as e:
    api_key = None

# Use a Mock query engine to return the data
class MockQE(object):
    def __init__(self):
        pass

    _materials = {
        'mp-10003': {u'energy_per_atom': -9.174497691666668,
                        u'is_hubbard': False,
                        u'material_id': u'mp-10003',
                        u'nsites': 12,
                        u'pretty_formula': u'Nb4CoSi',
                        u'volume': 194.5128160886403},
        'mp-10010': {u'energy_per_atom': -6.30060916,
                        u'is_hubbard': False,
                        u'material_id': u'mp-10010',
                        u'nsites': 5,
                        u'pretty_formula': u'Al(CoSi)2',
                        u'volume': 61.957194678711375},
        'mp-10015': {u'energy_per_atom': -8.66025992,
                        u'is_hubbard': False,
                        u'material_id': u'mp-10015',
                        u'nsites': 2,
                        u'pretty_formula': u'SiOs',
                        u'volume': 25.9156062823109},
        'mp-10018': {u'energy_per_atom': -4.0931096,
                        u'is_hubbard': False,
                        u'material_id': u'mp-10018',
                        u'nsites': 1,
                        u'pretty_formula': u'Ac',
                        u'volume': 45.384619900972496},
        'mp-10021': {u'energy_per_atom': -3.026048165,
                        u'is_hubbard': False,
                        u'material_id': u'mp-10021',
                        u'nsites': 2,
                        u'pretty_formula': u'Ga',
                        u'volume': 38.00766563190904},
        'mp-19306': {u'energy_per_atom': -6.709619538571429,
                        u'is_hubbard': True,
                        u'material_id': u'mp-19306',
                        u'nsites': 14,
                        u'pretty_formula': u'Fe3O4',
                        u'volume': 155.34118212181002},
        'mp-26': {u'energy_per_atom': -4.9257722625,
                        u'is_hubbard': False,
                        u'material_id': u'mp-26',
                        u'nsites': 4,
                        u'pretty_formula': u'La',
                        u'volume': 148.5978601715663}

    }

    def query(self, criteria, properties):
        mid_list = criteria["task_id"]["$in"]
        return [ self._materials.get(mid, None) for mid in mid_list ]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


def test_predict_k_g():
    (expected_k_value, expected_g_value, expected_caveat_str) = (175.30512291338607, 84.49987188140813, '')

    mpID = "mp-10003"
    (k_value, g_value, caveat_str) = elasticity.predict_k_g(mpID, query_engine=MockQE())
    assert (k_value, g_value, caveat_str) == (expected_k_value, expected_g_value, expected_caveat_str)

def test_predict_k_g_list():

    (expected_matid_list, expected_k_list, expected_g_list, expected_caveat_list) = (
     ['mp-10003', 'mp-10010', 'mp-10015', 'mp-10018', 'mp-10021', 'mp-19306', 'mp-26'],
     [175.30512291338607, 168.01218642160669, 265.96469661453744,
      45.15072359694464, 68.43138936905679, 136.86585554248228, 55.511505777303256],
     [84.49987188140813, 92.92207342120894, 118.409731828977, 19.816609506500367, 30.473676331990507,
      49.63871682171615, 24.379918816217213],
     ['', '', '',
      'Predictions are likely less reliable for materials containing F-block elements.',
      '',
        'Predictions may be less reliable for materials with non-GGA runs.',
        'Predictions are likely less reliable for materials containing F-block elements.'])

    mpID_list = ['mp-10003', 'mp-10010', 'mp-10015', 'mp-10018', 'mp-10021', 'mp-19306', 'mp-26']
    (matid_list, k_list, g_list, caveat_list) = elasticity.predict_k_g_list(mpID_list, query_engine=MockQE())
    assert (matid_list, k_list, g_list, caveat_list) == (expected_matid_list, expected_k_list, expected_g_list, expected_caveat_list)

@unittest.skipIf(api_key == None, reason="API key not defined")
def test_predict_k_g_remote():

    parser = ConfigParser.SafeConfigParser()
    cur_dir = os.path.dirname(__file__)
    parser.read('{}/test.ini'.format(cur_dir))
    api_key = parser.get('mprester', 'api_key')

    (expected_k_value, expected_g_value, expected_caveat_str) = (175.30512291338607, 84.49987188140813, '')

    mpID = "mp-10003"
    (k_value, g_value, caveat_str) = elasticity.predict_k_g(mpID, api_key)
    assert (k_value, g_value, caveat_str) == (expected_k_value, expected_g_value, expected_caveat_str)

@unittest.skipIf(api_key == None, reason="API key not defined")
def test_predict_k_g_list_remote():

    (expected_matid_list, expected_k_list, expected_g_list, expected_caveat_list) = (
     ['mp-10003', 'mp-10010', 'mp-10015', 'mp-10018', 'mp-10021', 'mp-19306', 'mp-26'],
     [175.30512291338607, 168.01218642160669, 265.96469661453744,
      45.15072359694464, 68.43138936905679, 136.86585554248228, 55.511505777303256],
     [84.49987188140813, 92.92207342120894, 118.409731828977, 19.816609506500367, 30.473676331990507,
      49.63871682171615, 24.379918816217213],
     ['', '', '',
      'Predictions are likely less reliable for materials containing F-block elements.',
      '',
        'Predictions may be less reliable for materials with non-GGA runs.',
        'Predictions are likely less reliable for materials containing F-block elements.'])

    mpID_list = ['mp-10003', 'mp-10010', 'mp-10015', 'mp-10018', 'mp-10021', 'mp-19306', 'mp-26']
    (matid_list, k_list, g_list, caveat_list) = elasticity.predict_k_g_list(mpID_list, api_key)
    assert (matid_list, k_list, g_list, caveat_list) == (expected_matid_list, expected_k_list, expected_g_list, expected_caveat_list)




if __name__ == '__main__':
    unittest.main()
