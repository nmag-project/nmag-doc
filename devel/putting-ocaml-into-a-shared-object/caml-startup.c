/* 
   (C) 2005 Dr. Thomas Fischbacher
 */

#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/custom.h>
#include <caml/alloc.h>

#include <sys/types.h>
#include <regex.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>

/* === Arg string decomposition ===
   (This is unrelated to Python)
*/

/* Our regexp:

"^[ \t]*                  # Whitespace
  (([^\'\" \t]|\\.)+|     # Non-Stringdelimiter/whitespace, or backslash and any char
  ('([^'\\]|\\.)*')|      # 'stuff'
  (\"([^\"\\]|\\.)*\"))"  # "stuff"
*/

static char *rx_arg="^[ \t]*((([^\'\" \t]|\\\\.)+|('([^'\\\\]|\\\\.)*')|(\"([^\"\\\\]|\\\\.)*\")))";

static int fetch_next_arg(regex_t *rxc, char *s, char **start_ptr, char **end_ptr)
{
  int status;
  regmatch_t next_arg[2];

  status=regexec(rxc,s,2,next_arg,0);

  if(status==REG_NOMATCH)
    {
      return 0;
    }

  *start_ptr=&s[next_arg[1].rm_so];
  *end_ptr=&s[next_arg[1].rm_eo];

  return next_arg[0].rm_eo;
  /* This is the same as next_arg[1].rm_eo, as the only difference is
     leading whitespace, but this does not matter.
   */
}

static int count_args(regex_t * rxc, char *s)
{
  int nr_args=0, pos=0, match_len;
  char *dummy1,*dummy2;

  for(;;)
    {
      match_len=fetch_next_arg(rxc,&s[pos],&dummy1,&dummy2);
      if(match_len==0)
	{
	  return nr_args;
	}
      else
	{
	  pos+=match_len;
	  nr_args++;
	}
    }
}

static void munch_arg(char *arg, int len)
{
  int i,k;

  /* Remove enclosing quotation marks */
  if(arg[0]=='"' || arg[0]=='\'')
    {
      for(i=1;i<len-1;i++)
	{
	  arg[i-1]=arg[i];
	}
      arg[len-2]=0;
      len-=2;
    }

  /* Eat Backslashes */
  i=0;
  while(i<len)
    {
      if((arg[i]=='\\') && (i<len-1)) /* Backslash followed by another char */
	{
	  for(k=i+1;k<len;k++)
	    {
	      arg[k-1]=arg[k];
	    }
	  arg[len-1]=0;
	  len--;
	}
      i++;
    }
}


static char **parse_args(char *s, int *ret_nr_args)
{
  regex_t compiled_rx;
  int nr_args,i,pos,skip_len,arg_len;
  char **ret;
  char *start_ptr, *end_ptr, *z;

  /* We may think about doing this only once and keeping the regex
     static; however, I do not want to make assumptions about the
     thread-safeness of compiled regexps.
   */

  if(0!=regcomp(&compiled_rx,rx_arg,REG_EXTENDED))
    {
      return 0;
    }
  
  nr_args=count_args(&compiled_rx,s);

  if(0==(ret=malloc((nr_args+1)*sizeof(char *))))
    {
      fprintf(stderr,"malloc() failed!\n");
      exit(1);
    }
  
  pos=0;

  for(i=0;i<nr_args;i++)
    {
      if(0==(skip_len=fetch_next_arg(&compiled_rx,&s[pos],&start_ptr,&end_ptr)))
	{
	  fprintf(stderr,"AIEE! Impossible situation!\n");exit(1);
	}

      arg_len=end_ptr-start_ptr;
      pos+=skip_len;

      if(0==(z=malloc((1+arg_len)*sizeof(char))))
	{
	  fprintf(stderr,"malloc() failure!\n");
	  exit(1);
	}

      strncpy(z,start_ptr,arg_len); z[arg_len]=0;
      munch_arg(z,arg_len); /* remove "", \, and such. */
      ret[i]=z;
    }
  
  ret[nr_args]=0;
  *ret_nr_args=nr_args;
  return ret;
}

/* === End arg string stuff === */

int tf_init_ocaml(char *arg_string)
{
    char **argv;
    int nr_args,i;

    argv=parse_args(arg_string,&nr_args);
    caml_startup(argv);

    for(i=0;i<nr_args;i++)
      {
	free(argv[i]);	    
      }
    free(argv);

    return 1;
}

