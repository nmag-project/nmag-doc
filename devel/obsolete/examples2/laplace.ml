#use "topfind";;
#require "pyfem2";; (* Just to make sure we have all the stuff... *)

#use "../mesh/mesh_toplevel_printers.ml";;

open Fem;;
open Snippets;;
open Mesh;;

let the_mesh = Fem.square_mesh 10 10 1.0 1.0;;

let elem_phi = make_element 2 ("phi",[||]) 1;;
let elem_rho = make_element 2 ("rho",[||]) 1;;

let mwe_phi = make_mwe "MWE-phi" (fun _ -> elem_phi) the_mesh;;
let mwe_rho = make_mwe "MWE-rho" (fun _ -> elem_rho) the_mesh;;
(* Note that we did not use mwe_sibling here! *)

let diffop_laplace = diffop_from_string "<d/dxj rho || d/dxj phi>, j:2";;

Printf.printf "=== PMX LAPLACE ===\n%!";;

let pmx_laplace = diffop_prematrix diffop_laplace mwe_rho mwe_phi;;

Printf.printf "=== PMX LAPLACE END ===\n%!";;

let app_laplace = prematrix_to_applicator pmx_laplace;;

let ml_laplace = prematrix_to_ml pmx_laplace;;

let ml_laplace_purified = Array.map (Array.map (fun x -> 0.5*.(float_of_int (int_of_float (2.01*.x))))) ml_laplace;;

(* XXX This should give us the zero vector!
# mx_x_vec ml_laplace (Array.make (Array.length ml_laplace) 1.0);;
- : float array = [|0.; 1.5; 1.5; 1.5; 4.; 1.5; 1.5; 1.5; 3.|]
*)


let (bdofs,k_nb,k_nn,mesh_laplace_solver) =
  laplace_solver
    ~dirichlet_bcs:[(Mesh.Body_Nr (-1),(Mesh.Body_Nr 1),
		     fun pos -> pos.(0) = 0.0 || pos.(0) = 9.0
		    )]
    pmx_laplace
;;

let ml_k_nb = Mpi_petsc.matrix_extract k_nb;;
let ml_k_nn = Mpi_petsc.matrix_extract k_nn;;

let ml_neg_dbc_values =
  mx_x_vec ml_k_nb (Array.map (fun x -> if x then 1.0 else 0.0) bdofs);;

let cofield_rho_in = app_laplace (sample_field mwe_phi (fun dof_name dof -> let x = dof.dof_pos.(0) and y = dof.dof_pos.(1) in if x >= 4.0 && x <= 6.0 && y >= 5.0 && y <= 7.0 then 1.0 else 0.0));;


let solved_phi =
  mesh_laplace_solver
    ~bc_fun:(fun dof_name dof -> if dof.dof_pos.(0) = 0.0 then 0.0 else 1.0) (* that's for the DBC DOFs *)
    cofield_rho_in
;;

mesh2d_plot_scalar_field
  ~scale:(fun [|x;y|] -> [|100.0+.x*.50.0;600.0-.y*.50.0|])
  ~plot_order:4
  ~plot_edges:true
  ("phi",[||])
  solved_phi 
  [|(-1.0,[|0.0;0.0;0.5|]);(1.2,[|0.8;0.8;1.0|]);|]
  "/tmp/plot-phi.ps"
;;

(*
Array.iter
  (fun p ->
     mesh2d_plot_scalar_field
       ~scale:(fun [|x;y|] -> [|100.0+.x*.50.0;600.0-.y*.50.0|])
       ~plot_order:4
       ~plot_edges:true
       ("phi",[||])
       (sample_field mwe_phi (fun dof_name dof -> if dof.dof_pos = p.mp_coords then 1.0 else 0.0))
       [|(-1.0,[|0.5;0.2;0.2|]);(1.2,[|1.0;0.6;0.6|]);|]
       (Printf.sprintf "/tmp/plot-part-phi-%f-%f.ps"
	  p.mp_coords.(0) p.mp_coords.(1)))
  the_mesh.mm_points
;;
*)
