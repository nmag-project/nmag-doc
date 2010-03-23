
#include <stdio.h>

extern int lib_main(char **);

static char *ml_sys_argv[]={"ocaml",0};

int main(int argc, char **argv) {
  fprintf(stderr,"DDD main.1\n");fflush(stderr);

  lib_main(ml_sys_argv);

  fprintf(stderr,"DDD main.2\n");fflush(stderr);
  return 0;
}
