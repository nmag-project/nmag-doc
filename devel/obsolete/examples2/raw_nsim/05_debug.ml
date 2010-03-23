
(* ===

nonsplit:

SX=  2 det=-25.000000 coords=[|[|0.; -2.5|]; [|0.; 2.5|]; [|5.; -2.5|]|]: coeff=-2.50000000 femfun #<FEM fun #49
 +   1.00000  (dL_0/dx_2)^2
> integral=3.12500000 contrib=-7.812500
DDD mx add   6   6: -7.81250000

split:

SX=  2 det=-25.000000 coords=[|[|1.; -2.5|]; [|1.; 2.5|]; [|6.; -2.5|]|]: coeff=-2.50000000 femfun #<FEM fun #49
 +   1.00000  (dL_0/dx_2)^2
> integral=6.12500000 contrib=-15.312500
DDD periodic row=[|6; 12|] col=12
DDD mx add   6  12: -15.31250000
DDD mx add  12  12: -15.31250000

=== *)

#use "topfind";;
#require "snippets";;

open Snippets;;

let di3=det_and_inv 3;;

di3 [|[|0.; -2.5;1.0|]; [|0.; 2.5;1.0|]; [|5.; -2.5;1.0|]|];;

di3 [|[|1.; -2.5;1.0|]; [|1.; 2.5;1.0|]; [|6.; -2.5;1.0|]|];;

(* =======================

        Objective Caml version 3.09.2

# Unknown directive `require'.
# Characters 0-13:
  open Snippets;;
  ^^^^^^^^^^^^^
Unbound module Snippets
# Wrong type of argument for directive `use'.
# - : unit = ()
Findlib has been successfully loaded. Additional directives:
  #require "package";;      to load a package
  #list;;                   to list the available packages
  #camlp4o;;                to load camlp4 (standard syntax)
  #camlp4r;;                to load camlp4 (revised syntax)
  #predicates "p,q,...";;   to set these predicates
  Topfind.reset();;         to force that packages will be reloaded
  #thread;;                 to enable threads

- : unit = ()
# /usr/lib/ocaml/3.09.2/unix.cma: loaded
/usr/lib/ocaml/3.09.2/str.cma: loaded
/usr/lib/ocaml/3.09.2/nums.cma: loaded
/usr/lib/ocaml/3.09.2/num-top: added to search path
/usr/lib/ocaml/3.09.2/num-top/num_top.cma: loaded
/usr/lib/ocaml/3.09.2/cryptokit: added to search path
/usr/lib/ocaml/3.09.2/cryptokit/cryptokit.cma: loaded
/usr/local/lib/ocaml/snippets: added to search path
/usr/local/lib/ocaml/snippets/snippets.cma: loaded
# # val di3 : float array array -> float * float array array = <fun>
# 1
  ;;
- : int = 1
# di3 [|[|0.; -2.5;1.0|]; [|0.; 2.5;1.0|]; [|5.; -2.5;1.0|]|];;

- : float * float array array =
(-25., [|[|-0.2; 0.; 0.2|]; [|-0.2; 0.2; 0.|]; [|0.5; 0.5; 0.|]|])
#   di3 [|[|1.; -2.5;1.0|]; [|1.; 2.5;1.0|]; [|6.; -2.5;1.0|]|];;

- : float * float array array =
(-25.0000000000000036,
 [|[|-0.199999999999999956; 2.77555756156289135e-17; 0.199999999999999983|];
   [|-0.199999999999999983; 0.2; 0.|];
   [|0.699999999999999845; 0.499999999999999889; -0.199999999999999956|]|])
# di3 [|[|1.; -2.5;1.0|]; [|1.; 2.5;1.0|]; [|6.; -2.5;1.0|]|];;
- : float * float array array =
(-25.0000000000000036,
 [|[|-0.199999999999999956; 2.77555756156289135e-17; 0.199999999999999983|];
   [|-0.199999999999999983; 0.2; 0.|];
   [|0.699999999999999845; 0.499999999999999889; -0.199999999999999956|]|])
# 
  
  
  di3 [|[|0.; -2.5;1.0|]; [|0.; 2.5;1.0|]; [|5.; -2.5;1.0|]|];;
- : float * float array array =
(-25., [|[|-0.2; 0.; 0.2|]; [|-0.2; 0.2; 0.|]; [|0.5; 0.5; 0.|]|])
# 0.25/.0.49;;
- : float = 0.510204081632653073
# 1;;
- : int = 1
# 
  
  1;;
- : int = 1
# (-15.31250000)/.(-7.81250000);;
- : float = 1.96
# 0.49/.0.25;;
- : float = 1.96
# di3  [|[|-0.199999999999999956; 2.77555756156289135e-17; 0.199999999999999983|];
   [|-0.199999999999999983; 0.2; 0.|];
   [|0.699999999999999845; 0.499999999999999889; -0.199999999999999956|]|];;
    - : float * float array array =
(-0.039999999999999987,
 [|[|1.; -2.5; 1.00000000000000022|]; [|0.999999999999999778; 2.5; 1.|];
   [|6.; -2.5; 1.|]|])
# di3 [|[|-0.2; 0.; 0.2|]; [|-0.2; 0.2; 0.|]; [|0.5; 0.5; 0.|]|];;
- : float * float array array =
(-0.0400000000000000078,
 [|[|0.; -2.5; 1.|]; [|0.; 2.5; 1.|]; [|5.; -2.5; 1.|]|])
# 

========================== *)


