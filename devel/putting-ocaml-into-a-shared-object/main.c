/* (C) 2005 Dr. Thomas Fischbacher */

#include <stdio.h>

extern int tf_init_ocaml(char *);

int main(void)
{
  printf("\nBefore the initialization of OCaml.\n");

  tf_init_ocaml("mycaml foo bar");

  printf("After the initialization of OCaml.\n\n");

  return 0;
}
