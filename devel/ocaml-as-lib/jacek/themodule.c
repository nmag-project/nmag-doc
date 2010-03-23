#include "Python.h"
#include "caml/callback.h"

static PyObject *
pycsquare(PyObject *self, PyObject *args)
{
  // The Python function's argument will be unpacked into arg
  int arg;
  // A handle for the caml function we are going to be calling
  static value *caml_square = NULL;

  // Get a handle on the Caml function we want to use.
  caml_square = caml_named_value("registered square");
  // Fail if necessary.
  if (caml_square == NULL) return NULL;


  // Here's where all the inter-language communication happens: A
  // Python integer is extractend into C, shoved into Caml, then the
  // returned Caml integer is extracted int C and inserted into
  // Python.
    
  // Extract a C int from the Python function's argument list, and store it in arg.
  if (!PyArg_ParseTuple(args, "i", &arg))
    return NULL; // If unpacking fails, fail. 
  // Turn C int into a Caml int                 (Val_int), 
  // apply caml_square to it                    (callback),
  // turn the resulting Caml int into a C int   (Int_val),
  // turn the C value into a Python value       (Py_BuildValue),
  // return it.
  return Py_BuildValue("i", Int_val(callback(*caml_square, Val_int(arg))));
}

static PyMethodDef TheModuleMethods[] = {
  {"pysquare", pycsquare, METH_VARARGS, "test function calling caml code"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initthemodule(void)
{
  char *dummy_argv[] = {0};
  caml_startup(dummy_argv);
  (void) Py_InitModule("themodule", TheModuleMethods);
}
