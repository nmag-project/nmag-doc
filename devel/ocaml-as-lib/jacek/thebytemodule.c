#include "Python.h"
#include "caml/callback.h"

static PyObject *
pycsquare(PyObject *self, PyObject *args)
{
  // The Python function's argument will be unpacked into arg
  int arg;
  // A handle for the caml function we are going to be calling
  static value *caml_square = NULL;

  // The Python function was passed a Python int. Extract a C int from
  // it, and store it in arg.
  if (!PyArg_ParseTuple(args, "i", &arg))
    return NULL; // If unpacking fails, fail. 
  
  // Get a handle on the Caml function we want to use.
  caml_square = caml_named_value("registered square");
  // Fail if necessary.
  if (caml_square == NULL) return NULL;
    
  // Turn C value into Caml value, 
  // apply caml_square to it,
  // turn the resulting Caml value into a C value,
  // turn the C value into a Python value
  // return it.
  return Py_BuildValue("i", Int_val(callback(*caml_square, Val_int(arg))));
}

static PyMethodDef TheModuleMethods[] = {
  {"pysquare", pycsquare, METH_VARARGS, "test function calling caml code"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initthebytemodule(void)
{
  char *dummy_argv[] = {0};
  caml_startup(dummy_argv);
  (void) Py_InitModule("thebytemodule", TheModuleMethods);
}
