//
//   Copyright 2015, The Materials Project
//

#include <Python.h>
#include <numpy/arrayobject.h>
#include "predict.h"

// package:     gbml
// module:      core
// function:    predict
// python call: gbml.core.predict()

static char core_docstring[] =
    "This module provides an interface to the core GBM-Locfit C functions.";
static char predict_docstring[] =
    "This module makes GBM-Locfit predictions.";

// gbml_predict function definition
static PyObject *gbml_predict(PyObject *self, PyObject *args)
{
  char *filename;
  int nPredictions;
  PyArrayObject *descriptorsArray;
  PyArrayObject *predictionsArray;

  // Parse input tuple
  if (!PyArg_ParseTuple(args, "siO!O!", &filename, &nPredictions, &PyArray_Type, &descriptorsArray,
      &PyArray_Type, &predictionsArray)) return NULL;

  // Get C pointer to descriptors
  double *descriptors = (double *) PyArray_DATA(descriptorsArray);
  double *predictions = (double *) PyArray_DATA(predictionsArray);

  // Call external C function
  int flag = predict(filename, nPredictions, descriptors, predictions);

  if (flag)
    return NULL;
  else
    Py_RETURN_NONE;
};

// method definition
static PyMethodDef core_methods[] = {
    {"predict", gbml_predict, METH_VARARGS, predict_docstring},
    {NULL, NULL, 0, NULL} 
};

static struct PyModuleDef coredef = {
    PyModuleDef_HEAD_INIT,
    "core",
    core_docstring,
    -1,
    core_methods,
};

// init function definition
PyMODINIT_FUNC
PyInit_core(void)
{
  /*m = PyModule_Create(coredef);*/
  import_array();
  return PyModule_Create(&coredef);
};


