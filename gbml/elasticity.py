"""
Predict bulk (K) and shear (G) moduli using gbml (GBM-Locfit).
Queries MP db for specified material(s), computes descriptors, and calls gbml.core.predict().
"""

from __future__ import print_function

from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
from pymatgen.ext.matproj import MPRester
import gbml.core # GBM-Locfit core module
import json
import math
import numpy as np
import os
import six
import sys

__author__ = 'Randy Notestine'
__copyright__ = 'Copyright 2016, The Materials Project'
__version__ = '1.00'
__maintainer__ = 'Randy Notestine'
__email__ = 'RNotestine@ucsd.edu'
__date__ = 'March 14, 2016'

API_KEY = MPRester().api_key
CAVEAT_AIAB = 'Unable to estimate cohesive energy for material.'
CAVEAT_F_BLOCK = 'Predictions are likely less reliable for materials containing F-block elements.'
CAVEAT_HUBBARD = 'Predictions may be less reliable for materials with non-GGA runs.'
DATAFILE_AIAB = 'data/element_aiab_energy.json'
DATAFILE_K = 'data/gbml-K-v1.00.data'
DATAFILE_G = 'data/gbml-G-v1.00.data'
VERY_SMALL = 1E-12

def get_element_aiab_energy(element):
    """
    Lookup atom-in-a-box (aiab) energy for specified element.
    Used to estimate cohesive energy of a compound from the compound's VASP energy.
    The elemental atom-in-a-box energies were provided by Wei Chen and Maarten de Jong.
    Returns atom-in-a-box energy for specified element.
    :param element:
    :return: aiab_energy
    """

    element_aiab_energy_dict = None

    try:
        with open(os.path.join(os.path.dirname(__file__),DATAFILE_AIAB),'r') as json_file:
            element_aiab_energy_dict = json.load(json_file)

    finally:
        if element_aiab_energy_dict is None:
            return None

        object = element_aiab_energy_dict.get(element)
        if object is not None:
            return object[0]

        return None


def holder_mean(values, power, weights=None, weights_norm=None):
  """
  Calculate (possibly weighted) Holder (or power) mean
  :param: values: list or array of (typically elemental property) values
  :param: power: the power to which the Holder mean is to be calculated
  :param: weights: option list or array of weights
  :param: weights_norm: optional normalization method to be applied to weights
  :return: holder_mean
  """

  values = np.array(values, dtype=float)
  power = float(power)

  # Make sure weights match length and are normalized
  if weights is None:
    alpha = 1 / len(values)
  else:
    weights = np.array(weights, dtype=float)
    if len(values) != len(weights):
      #warn('Holder.mean returned zero when passed length mis-matched values and weights', UserWarning)
      sys.stderr.write('  Warning: Holder.mean returned zero when passed length mis-matched values and weights\n')
      return 0.0
    if weights_norm is not None:
      if weights_norm == "max" and max(weights) != 1.0:
        weights = weights / max(weights)
      elif weights_norm == "sum" and sum(weights) != 1.0:
        weights = weights / sum(weights)
      else:
        #warn('Holder.mean returned zero when passed unknown weights_norm method', UserWarning)
        sys.stderr.write('  Warning: Holder.mean returned zero when passed unknown weights_norm method\n')
        return 0.0
    alpha = 1 / sum(weights)

  if power == 0.0:  # geometric mean
    if any(abs(value) < VERY_SMALL for value in values):
      return 0.0
    elif any(value < 0 for value in values):
      #warn('Holder.mean returned zero when passed negative value with zero power', UserWarning)
      sys.stderr.write('  Warning: Holder.mean returned zero when passed negative value with zero power\n')
      sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
      return 0.0

    if weights is None:
      return math.exp(alpha * sum(np.log(values)))
    else:
      return math.exp(alpha * sum(weights * np.log(values)))

  elif power == 1.0:  # arithmetic mean
    return np.average(values, weights=weights)

  if any(value < 0 for value in values):
    if power % 2 != 0.0:
      #warn('Holder.mean returned zero when passed negative value with non-even power', UserWarning)
      sys.stderr.write('  Warning: Holder.mean returned zero when passed negative value with non-even power\n')
      sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
      return 0.0

  if weights is None:
    return pow(alpha * sum(np.power(values, power)), 1/power)
  else:
    return pow(alpha * sum(weights * np.power(values, power)), 1/power)


def _get_mp_query(api_key=None, query_engine=None):
    """
    Returns object that can query the MP DB. This is either a local query_engine
    or an MPRester object. We can do this because both MPRester and QueryEngine
    expose the same query interface
    """
    if query_engine:
        return query_engine
    elif api_key:
        return MPRester(api_key)
    else:
        raise Exception("missing API KEY or query engine")



def predict_k_g_list(material_id_list, api_key=API_KEY, query_engine=None):
    """
    Predict bulk (K) and shear (G) moduli for a list of materials.
    :param material_id_list: list of material-ID strings
    :param api_key: The API key used by pymatgen.matproj.rest.MPRester to connect to Materials Project
    :param query_engine: (Optional) QueryEngine object used to query materials instead of MPRester

    :return: (matid_list, predicted_k_list, predicted_g_list, caveats_list)
    Note that len(matid_list) may be less than len(material_id_list),
    if any requested material-IDs are not found.
    """

    if len(material_id_list) == 0 or not isinstance(material_id_list, list):
        return (None, None, None, None)  # material_id_list not properly specified

    lvpa_list = []
    cepa_list = []
    rowH1A_list = []
    rowHn3A_list = []
    xH4A_list = []
    xHn4A_list = []
    matid_list = []
    k_list = []
    g_list = []
    caveats_list = []
    aiab_problem_list = []

    # TODO: figure out if closing the query engine (using 'with' ctx mgr) is an issue
    # If it is a problem then try manually doing a session.close() for MPRester, but ignore for qe

    mpr = _get_mp_query(api_key, query_engine)
    for entry in mpr.query(criteria={"task_id": {"$in": material_id_list}}, properties=
        ["material_id", "pretty_formula", "nsites", "volume", "energy_per_atom", "is_hubbard"]):

        caveats_str = ''
        aiab_flag = False
        f_block_flag = False
        weight_list = []
        energy_list = []
        row_list = []
        x_list = []

        # Construct per-element lists for this material
        composition = Composition(str(entry["pretty_formula"]))
        for element_key, amount in composition.get_el_amt_dict().items():
            element = Element(element_key)
            weight_list.append(composition.get_atomic_fraction(element))
            aiab_energy = get_element_aiab_energy(element_key)  # aiab = atom-in-a-box
            if aiab_energy is None:
                aiab_flag = True
                break
            energy_list.append(aiab_energy)
            if element.block == 'f':
              f_block_flag = True
            row_list.append(element.row)
            x_list.append(element.X)

        # On error, add material to aiab_problem_list and continue with next material
        if aiab_flag:
            aiab_problem_list.append(str(entry["material_id"]))
            continue

        # Check caveats
        if bool(entry["is_hubbard"]):
            if len(caveats_str) > 0: caveats_str += " "
            caveats_str += CAVEAT_HUBBARD
        if f_block_flag:
            if len(caveats_str) > 0: caveats_str += " "
            caveats_str += CAVEAT_F_BLOCK

        # Calculate intermediate weighted averages (WA) for this material
        ewa = np.average(energy_list, weights=weight_list)      # atom-in-a-box energy WA

        print(str(entry["material_id"]))

        # Append descriptors for this material to descriptor lists
        lvpa_list.append(math.log10(float(entry["volume"]) / float(entry["nsites"])))
        cepa_list.append(float(entry["energy_per_atom"]) - ewa)
        rowH1A_list.append(holder_mean(row_list, 1.0, weights=weight_list))
        rowHn3A_list.append(holder_mean(row_list, -3.0, weights=weight_list))
        xH4A_list.append(holder_mean(x_list, 4.0, weights=weight_list))
        xHn4A_list.append(holder_mean(x_list, -4.0, weights=weight_list))
        matid_list.append(str(entry["material_id"]))
        caveats_list.append(caveats_str)

    if isinstance(mpr, MPRester):
        mpr.session.close()

    # Check that at least one valid material was provided
    num_predictions = len(matid_list)
    if num_predictions > 0:
        # Construct descriptor arrays
        if (len(lvpa_list) != num_predictions or len(cepa_list) != num_predictions or
            len(rowH1A_list) != num_predictions or len(rowHn3A_list) != num_predictions or
            len(xH4A_list) != num_predictions or len(xHn4A_list) != num_predictions):
                return (None, None, None, None)
        k_descriptors = np.ascontiguousarray([lvpa_list, rowH1A_list, cepa_list, xHn4A_list],
            dtype=float)
        g_descriptors = np.ascontiguousarray([cepa_list, lvpa_list, rowHn3A_list, xH4A_list],
            dtype=float)

        # Allocate prediction arrays
        k_predictions = np.empty(num_predictions)
        g_predictions = np.empty(num_predictions)

        # Make predictions
        k_filename = os.path.join(os.path.dirname(__file__),DATAFILE_K)
        g_filename = os.path.join(os.path.dirname(__file__),DATAFILE_G)
        gbml.core.predict(k_filename, num_predictions, k_descriptors, k_predictions)
        gbml.core.predict(g_filename, num_predictions, g_descriptors, g_predictions)

        k_list = np.power(10.0, k_predictions).tolist()
        g_list = np.power(10.0, g_predictions).tolist()

    # Append aiab problem cases
    for entry in aiab_problem_list:
        matid_list.append(entry)
        k_list.append(None)
        g_list.append(None)
        caveats_list.append(CAVEAT_AIAB)

    if len(matid_list) == 0:
        return (None, None, None, None)
    else:
        return (matid_list, k_list, g_list, caveats_list)


def predict_k_g(material_id, api_key=API_KEY, query_engine=None):
    """
    Predict bulk (K) and shear (G) moduli for one material.
    :param material_id: material-ID string
    :param api_key: The API key used by pymatgen.matproj.rest.MPRester to connect to Materials Project
    :param query_engine: (Optional) QueryEngine object used to query materials instead of MPRester

    :return: (predicted_k, predicted_g, caveats)
    Note that None may be returned for predicted_k and predicted_g when caveats is not None.
    """

    if len(material_id) == 0 or not isinstance(material_id, six.string_types):
        return (None, None, None)  # material_id not properly specified

    (material_id_list, k_list, g_list, caveats_list) = predict_k_g_list([material_id], api_key, query_engine)

    if material_id_list is None:
        return (None, None, None)  # material_id not found in MP db

    return (k_list[0], g_list[0], caveats_list[0])
