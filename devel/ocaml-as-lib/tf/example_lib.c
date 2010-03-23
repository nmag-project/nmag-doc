
#include <caml/callback.h>
#include <caml/mlvalues.h>

#include <stdio.h>

/* This function will go into a C shared object library and will, when
   called, initialize and utilize compiled OCaml code. */

int lib_main(char **argv) {
  /* Initialize the OCaml runtime system (GC, etc.) -- this will
     register some caml functions for being used with the callback
     mechanism:
  */  

  fprintf(stderr,"DDD lib_main.1\n");fflush(stderr);
  caml_main(argv);
  fprintf(stderr,"DDD lib_main.2\n");fflush(stderr);

  fprintf(stderr,"DDD lib_main demo is at: %p\n",*(caml_named_value("caml_callback_demo")));fflush(stderr);

  
  /* Call the ocaml demo() function: */
  caml_callback(*(caml_named_value("caml_callback_demo")),Val_unit);

  fprintf(stderr,"DDD lib_main.3\n");fflush(stderr);

  return 0;
}
