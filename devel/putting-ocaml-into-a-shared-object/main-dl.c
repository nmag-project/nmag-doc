/* (C) 2005 Dr. Thomas Fischbacher */

#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>

extern int tf_init_ocaml(char *);

int main(void)
{
  void *handle;
  int (*f)(char *);
  char *error;

  if(!(handle=dlopen("./caml-in-lib.so", RTLD_LAZY)))
    {
      fprintf(stderr, "AIEE: %s\n", dlerror());
      exit(1);
    }

  dlerror(); /* Clear existing errors */
  *(void **)(&f) = dlsym(handle, "tf_init_ocaml");
  if(0!=(error=dlerror()))
    {
      fprintf (stderr, "AIEE: %s\n", error);
      exit(1);
    }

  printf("\nBefore the initialization of OCaml.\n");

  (*f)("mycaml foo bar");

  printf("After the initialization of OCaml.\n\n");

  return 0;
}

