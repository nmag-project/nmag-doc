(* 
   (C) 2006 Dr. Thomas Fischbacher

   Implementing code to handle hierarchical octree subdivisions
   of three-dimensional space.

   We might consider using a different tree data structure where every
   box node has two children only, and subdivides in x-, y-, or
   z-direction. This seems feasible, but extra checks would be
   necessary to ensure the boxes stay roughly cubic. Maybe not such a
   good idea.

   But for now, let's stick with a simple octree...

*)

#use "topfind";;
#require "qhull";;
#require "snippets";;

open Snippets;;

type 'a octree =
    {
      oct_center: float array;
      oct_edge_length: float;
      mutable oct_parent: 'a octree option;
      mutable oct_childs: 'a octree option array;
      (* ^ Either [||] or an 8-vector *)
      mutable oct_contents: 'a option;
    }
;;

let print_octree fmt oct =
  Format.fprintf fmt "#<octree center=%s edge_length=%f>"
    (float_array_to_string oct.oct_center) oct.oct_edge_length
;;

#install_printer print_octree;;


let octree_cell_is_fundamental c = (0=Array.length c.oct_childs);;

(* Note that the Seberino NNNPNN method actually does too much work.
   It is true that we have to start at the boxes which are
   childof-childof-largest, but then we iterate over all boxes down to
   the childof-childof-largest level, and introduce the CN whenever we
   encounter a box where every point is farther away than every point
   inside the observer box.

   XXX Need a better explanation!

   Note that both in the CF and the CN, every multipole is being
   accounted for multiple times, once at every level. This evidently
   is a necessity, because this way, we ensure that by taking the
   contribution from a large block, we already accounted for
   everything that also is inside it.

   * CFs are computed small-to-large. At the lowest level, we receive
   contributions from the fundamental multipoles, while everywhere
   else we just combine contributions from lower levels.

   We combine the contributions from small cells
   to ever larger cells (easy: when recursing down through
   cells, just ensure that we first compute CF for the smaller
   cells before we combine them to obtain CF at the outer level.)

   - The CF contribution of a box is easy to understand: it is the
   accumulated multipole moment of everything that is inside, taken
   relative to the box's center.

   * CNs are computed large-to-small.
   
   - What does the CN contribution of a box express? If we take any
   point inside the box and evaluate the box's CN for that point,
   plus include all its contents and the contents of its immediate
   neighbours via simple multipole fields, we get the contribution
   at that particular point.

   The idea is that we can use this property to construct the recursion
   down to lower levels: when we compute the CN of a child box,

   * We include the contribution from its parent box and translate that.
   Furthermore:

   * We have to account for the parent's nearest neighbours.
   So, we iterate over them and use CFCN where appropriate.
   We subdivide them one level (same-size-as-we-are boxes),
   and we skip our own nearest neighbours.

   That's it, basically. How to determine our nearest neighbours?
   Easy: They are so close that they may touch our own circumcircle.

*)

type octree_nav = int array;;
(* An octree navigation description is an array of integers, where
   [n] (in the range 0-7) means: enter the n'th child, and (-1)
   means: step up one level and enter the parent.
*)

let locate_in_octree cell nav =
  let nr_steps = Array.length nav in
  let rec walk cell_now pos =
    if pos=nr_steps then cell_now
    else
      let dir = nav.(pos) in
	if dir = (-1)
	then
	  match cell_now.oct_parent with
	    | None -> failwith "locate_in_octree: tried to enter non-existing parent!"
	    | Some p ->
		walk p (1+pos)
	else if octree_cell_is_fundamental cell_now
	then
	  failwith "locate_in_octree: tried to enter child of a fundamental cell!"
	else
	  match cell_now.oct_childs.(dir) with
	    | None ->
		failwith "locate_in_octree: tried to enter non-existing child!"
	    | Some c ->
		walk c (1+pos)
  in walk cell 0
;;


type ('t_translator,'t_dipole) fmm_cell =
    {
      fmm_dipoles: 't_dipole array; 
      (* Typically, 't_dipole is (int array),
	 using dof_indices for M_x, M_y, M_z dofs.
	 We then can also derive locations from these!
      *)
      (* These will be set when we feed a new field into the octree: *)
      mutable fmm_CF_coeffs: float array;
      mutable fmm_CN_coeffs: float array;
      (* ------ *)
      (* These belong to the initial set-up *)
      mutable fmm_CF_add_local_dipole_contribs: (unit -> unit);
      (* ^ walks all the dipoles, modifies the data in fmm_CF_coeffs *)
      mutable fmm_CF_translator_to_parent: 't_translator;
      mutable fmm_CN_translator_from_parent: 't_translator;
      (* --- here, we have to include some slightly hairy topological stuff --- *)
      mutable fmm_other_boxes_wciocs: ('t_translator,'t_dipole) fmm_cell array;
      (* wciocs = "with contents inside our circumsphere" *)
      mutable fmm_other_boxes_wcoiopcs_and_CFCN_translators:
	(('t_translator,'t_dipole) fmm_cell * 't_translator) array;
      (* wcoiopcs = "with contents only inside
	 our parents' circumsphere (and not in our own)"
      *)
    }
;;

(* Conceptually, we definitely should support non-uniform octrees,
   which means that we subdivide whenever occupation density gets too
   high.
*)

let make_octree ?(max_sites_per_cell=20) ~fun_content ~fun_coords sites =
  let no_childs = [||] in
  let sites_coords = Array.map fun_coords sites in
  let center_of_gravity = float_arrays_avg sites_coords in
    (* XXX Here, we do face a problem if the center of gravity hits the position
       of a dipole, as the far-field translation matrix for distance zero
       is ill-defined!
    *)
  let largest_distance =
    Array.fold_left
      (fun sf x ->
	 let d = sqrt(euclidean_distance_sq center_of_gravity x) in
	   max d sf)
      0.0 sites_coords
  in
  let largest_edge_length=ceil (2.01*.largest_distance) in
    (* XXX Here, we rather should make sure that largest_edge_length
       is a "nice" exactly float-representable number, such as 0.25,
       or 0.06125 or 8.0 or such. Just using "ceiling" is a bad hack
       which will work for most practical purposes, but we REALLY 
       should take the original number and round it up to, say,
       8 binary digits. base-2 logarithmizing will help.

       We actually should have a round_to_n_valid_digits ?(base=10)
       function in snippets.ml!
    *)
  let () = Printf.printf "DDD please remove the edge length hack in make_octree!\n%!" in
  let rec make_box opt_parent nr_child center edge_length contents =
    (* let () = Printf.printf "DDD make_box center=%s\n" (float_array_to_string center) in *)
    let b = {oct_center = center;
	     oct_edge_length = edge_length;
	     oct_parent = opt_parent;
	     oct_childs = no_childs;
	     oct_contents = None;
	    }
    in
    let () =
      match opt_parent with
	| None -> ()
	| Some c -> c.oct_childs.(nr_child) <- Some b
    in
    let () =
      if Array.length contents <= max_sites_per_cell
	(* Atomic box *)
      then
	b.oct_contents <- Some (fun_content contents)
      else
	let () = b.oct_contents <- Some (fun_content [||]) in 
	  (* Everything will go into our childs! *)
      	let () = b.oct_childs <- Array.make 8 None in
	let d_center = edge_length /. 4.0 in
	let subdivided_contents = Array.make 8 [] in
	let () = Array.iter
	  (fun c ->
	     let pos = fun_coords c in
	     let octant =
	         (if pos.(0) >= center.(0) then 1 else 0)
	       + (if pos.(1) >= center.(1) then 2 else 0)
	       + (if pos.(2) >= center.(2) then 4 else 0)
	     in
	       subdivided_contents.(octant) <-
		 c :: subdivided_contents.(octant)
	  )
	  contents
	in
	  Array.iteri
	    (fun nr_octant sub_contents ->
	       (* let () = Printf.printf "DDD OCT making child nr %d - nr contents=%d\n" nr_octant (List.length sub_contents) in *)
	       let v_sub_contents = Array.of_list (List.rev sub_contents) in
		 (* The List.rev is actually not necessary,
		    but helps with debugging. We try to be tidy.
		 *)
		 if Array.length v_sub_contents = 0
		 then () (* We leave None in this place of the child option array. *)
		 else
		   let child_center = 
		     Array.init 3
		       (fun n ->
			  center.(n)
			  +.(if nr_octant land (1 lsl n) <> 0 then 1.0 else -1.0)
			  *.d_center)
		   in
		   let child =
		     make_box (Some b)
		       nr_octant child_center (edge_length/.2.0)
		       v_sub_contents
		   in b.oct_childs.(nr_octant) <- Some child
	    )
	    subdivided_contents
    in
      b
  in make_box None (-1) center_of_gravity largest_edge_length sites
;;

(* We may actually set the FMM octree translator datastructures all in one go,
   but for the sake of modularity, it is much clearer to do one thing
   after another, and start with setting the CF combiners...
*)

let walk_octree ?(top_down=true) f octree =
  let rec walk_td cell =
    let () = f cell in
      for i=0 to Array.length cell.oct_childs-1 do
	match cell.oct_childs.(i) with
	  | None -> ()
	  | Some c -> walk_td c
      done
  in
  let rec walk_bu cell =
    begin
      for i=0 to Array.length cell.oct_childs-1 do
	match cell.oct_childs.(i) with
	  | None -> ()
	  | Some c -> walk_bu c
      done;
      f cell
    end
  in if top_down then walk_td octree else walk_bu octree
;;

(* This is a once-to-be-called setup function *)
let fmm_setup_octree_CF_computation
    ?(max_multipole_order=10)
    ~fun_dipole_position (* dipole -> float array *)
    ~fun_process_dipole_coeffs
    (* Pass on the dipole's coefficients to a continuation.
       Signature: fun dipole cont -> (), where cont is (fun x y z -> ())
    *)
    octree =
  let ht_translator_by_dist = Hashtbl.create 20 in
  let ht_add_ff_coeffs_shifted_dipole_by_dist = Hashtbl.create 20 in
  let translator_by_dist dist =
    memoized ht_translator_by_dist (mp_translator_ff max_multipole_order) dist
  in
  let fun_add_ff_coeffs_shifted_dipole dist =
    let fun_work =
      memoized ht_add_ff_coeffs_shifted_dipole_by_dist
	(mp_add_ff_coeffs_of_shifted_dipole max_multipole_order)
	dist
    in
      fun_work
  in
    walk_octree
      ~top_down:false
      (fun cell ->
	 match cell.oct_contents with
	   | None -> impossible()
	   | Some contents ->
	       begin
		 (* 1.: Setting the translator that will add our coefficients to our parent *)
		 (match cell.oct_parent with
		    | None -> ()
		    | Some parent ->

			let dist = array_pointwise (-.) parent.oct_center cell.oct_center in (* XXX sign right? *)
			let translator = translator_by_dist dist in
			  contents.fmm_CF_translator_to_parent <- translator);
		 (* 2.: Setting the function that adds the CF contributions from the dipoles
		    in this cell *)
		 if Array.length contents.fmm_dipoles = 0
		 then () (* Do nothing - the default dummy function is good enough for us! *)
		 else
		   let dipole_cf_contrib_adders =
		     Array.map
		       (fun dipole ->
			  let dp_pos = fun_dipole_position dipole in
			  let dp_dist = array_pointwise (-.) cell.oct_center dp_pos in (* XXX Note that this sign matches the sign above in let dist = ... *)
			  let adder =
			    fun_add_ff_coeffs_shifted_dipole dp_dist
			      contents.fmm_CF_coeffs (* result_array *)
			  in
			    fun () -> fun_process_dipole_coeffs dipole adder
		       )
		       contents.fmm_dipoles
		   in
		   let fun_CF_add_local_dipole_contribs () =
		     for i=0 to (Array.length dipole_cf_contrib_adders)-1 do
		       let () = dipole_cf_contrib_adders.(i) () in ()
		     done
		   in
		     contents.fmm_CF_add_local_dipole_contribs <- fun_CF_add_local_dipole_contribs
	       end
      )
      octree
;;

(* This is a dynamic run-time computation function *)
let fmm_octree_set_CF_coeffs ?(max_multipole_order=10) octree =
  begin
    (* We start out by clearing the CF data vectors everywhere *)
    walk_octree
      (fun cell ->
	 match cell.oct_contents with
	   | None -> impossible()
	   | Some contents ->
	       for i=0 to Array.length contents.fmm_CF_coeffs-1 do
		 contents.fmm_CF_coeffs.(i) <- 0.0
	       done
      ) octree;
    walk_octree
      ~top_down:false
      (fun cell ->
	 match cell.oct_contents with
	   | None -> impossible()
	   | Some contents ->
	       let () =
		 (if octree_cell_is_fundamental cell
		  then contents.fmm_CF_add_local_dipole_contribs ()
		  else ())
	       in
		 match cell.oct_parent with
		   | None -> ()
		   | Some parent ->
		       (match parent.oct_contents with
			  | None -> impossible ()
			  | Some p_contents ->
			      apply_translator_add_contribs
				p_contents.fmm_CF_coeffs 
				contents.fmm_CF_translator_to_parent
				contents.fmm_CF_coeffs
		       ))
      octree
  end
;;

(* We now are able to map charge distributions to their far-field multipole expansions. *)

(* Next: computing the near-field (taylor) expansion coefficients.

   First problem: mapping far field to near field. We have to keep in mind:

   - Every cell will receive a CN/CN-translated
     near field contribution from its parent.

   - Furthermore, every cell will receive CN/CF contributions from all the cells
     which touch the parent's circumsphere, but not this cell's.

   The second condition is conceptually quite different from the way
   how things are handled in the Seberino/Bertram paper, where they use
   a not-nearest-neighbour-parent-nearest-neighbour scheme. 

   The advantage of our scheme is that

   - It is conceptually simpler
   - It works with non-uniform octree decompositions
   - It may save some work, by including every CF contribution
     at the highest possible level (hence we do not have to consider
     it multiple times at lower levels). (NOTE: think once more about this - 
     this argument might be flawed - right now I think it should hold...)
*)

(* We start out with some harmless setup issues: making all the matrices for
   CN-parent-to-CN-child translation *)

let fmm_setup_octree_CNCN_computation
    ?(max_multipole_order=10)
    octree =
  let ht_translator_by_dist = Hashtbl.create 20 in
  let translator_by_dist dist =
    memoized ht_translator_by_dist (mp_translator_nn max_multipole_order) dist
  in
    walk_octree	(* Actually does not matter if we do this top-down or not... *)
      (fun cell ->
	 match cell.oct_contents with
	   | None -> impossible()
	   | Some contents ->
	       match cell.oct_parent with
		 | None -> ()
		 | Some parent ->
		     
		     let dist = array_pointwise (-.) cell.oct_center parent.oct_center in (* XXX sign right? *)
		     let translator = translator_by_dist dist in
		       contents.fmm_CN_translator_from_parent <- translator
      )
      octree
;;

(* For the FMM, we have to know when one box touches the circumsphere of another box.
   Unfortunately, that is slightly ugly to find out. 

   BUT: there is a better way of doing this - we actually can further improve
   the algorithm here!

   The idea is: boxes contain dipoles. We only have to make sure that
   no dipole from box #2 is contained in the circumsphere of box
   #1. If this holds, we can treat box #2's multipole contents as
   sufficiently far away to be converted into near-field coefficients.

   Hence, we need a helper:
*)

let octree_radius oct = sqrt(3.0)*.0.5*.oct.oct_edge_length;;

let octree_box2_contains_dipole_inside_circumsphere_of_box1
    ~fun_dipole_position (* dipole -> float array *)
    box1 box2 =
  let box1_radius = octree_radius box1
  and box1_center = box1.oct_center
  in
  let rec walk sub_box =
    let r_sub = octree_radius sub_box
    and c_sub = sub_box.oct_center
    in
      if sqrt(euclidean_distance_sq c_sub box1_center) > 1.001*.(r_sub+.box1_radius)
	(* These boxes are so far away from one another that the sub-box just cannot
	   have contents that fall into box1's circumsphere!
	*)
      then false
      else
	(* We have to inspect our contents as well as our children.
	   We do our contents first, as this will often be an array of length 0...
	*)
	match sub_box.oct_contents with
	  | None -> impossible()
	  | Some c ->
	      let dipoles = c.fmm_dipoles in
	      let nr_dipoles = Array.length dipoles in
	      let rec walk_dipoles n =
		if n=nr_dipoles then false
		else 
		  let pos = fun_dipole_position dipoles.(n) in
		    if sqrt(euclidean_distance_sq pos box1_center) <= box1_radius
		    then true (* found a dipole that is inside the circumsphere *)
		    else walk_dipoles (1+n)
	      in
	      let some_dipole_inside_circumsphere = walk_dipoles 0 in
		if some_dipole_inside_circumsphere then true
		else
		  let childs = sub_box.oct_childs in
		  let nr_childs = Array.length childs in
		  let rec walk_childs n =
		    if n=nr_childs then false
		    else
		      match childs.(n) with
			| None -> walk_childs (1+n)
			| Some child ->
			    if walk child then true else walk_childs (1+n)
		  in walk_childs 0
  in walk box2
;;

(* Now for the tricky part: for every box B,
   we have to find all the other boxes which

   - have dipoles inside B's circumsphere

   as well as all the other boxes which
   
   - do not have dipoles inside B's circumsphere, but did have dipoles
     inside B's parents' circumsphere.

   Note that if we proceed this way, a large box will want to include
   contributions from many small "fringe" boxes which approximate the
   outer of its circumsphere as close as possible.

   Is that reasonable?

   After all, do we really need that second step in the FMM?  One
   might argue that, in order to determine the potential at any point
   inside the dipole distribution, we may proceed by recursively
   adding the far-field multipole contributions from all far-away
   boxes, and treat the contents of the too-close boxes specially.

   This certainly would work, but the essence is: by going to the near
   field expansion, we can further reduce the effort, since we do not
   have to walk all boxes, but only consider the one the observer
   point is located in plus all closeby companions.

   (Nevertheless, the "naive" method is very valuable, in particular
   when we want to map the potential close to but outside the largest
   box - and also for consistency checks!)

   So, when we do the near-field method, when do we convert a
   far-field contribution to a near-field one? Eventually, we would
   like to do this at the highest possible level, because (1) wherever
   we go, we eventually do have to include that contribution anyway,
   and (2) the higher we do this the bigger the chance we do not have
   to do it multiple times when splitting.

   There is, however, one objection to including contributions always
   at the highest possible level, namely that this may mean moving
   multipoles very far. So, when we recurse down a box, the
   coefficients from the multipole expansion of a closeby box may
   enter through first having gone through the multipole expansion of
   a much larger box whose center is far away. That certainly is not
   good for numerical accuracy.

   So, we do need an extra element in the algorithm which ensures that
   this does not happen. For this, we adopt the rule that

   the near-field contributions to a box with children must not come
   from other boxes with children which are smaller than the present
   box. (I.e. in such a case, we will leave it to the children to
   include these contributions).

   But actually, this extra rule is an additional screw with which we
   can parametrize and tune the algorithm. Maybe it would be
   interesting to see how bad it actually is to instead follow the
   "recurse down maximally" rule.
*)

(* Setup function: ensure we have got all the topology information
   for the CNCF computation
*)

let fmm_setup_topology octree =
  let ht_translator_by_dist = Hashtbl.create 20 in
  let translator_by_dist dist =
    memoized ht_translator_by_dist (mp_translator_nf max_multipole_order) dist
  in
    walk_octree
      ~top_down:true
      (fun cell ->
	 failwith "XXX WRITE ME!"
      )
      octree
;;



(* === Test Code === *)

type ddd_dipole =
    {
      dp_pos: float array;
      dp_xyz: float array;
    }
;;

let test_dipoles_random =
  Array.init 100 (* 1000 *)
    (fun n ->
       let pos = Array.init 3 (fun _ -> -5.0+.(Random.float 10.0)) in
	 {dp_pos=pos;
	  dp_xyz=Array.make 3 0.0;
	 }
    )
;;

(* This just is some (more or less artistic) spatial distribution of a
   considerable number of dipoles which we use to compare the
   multipole-approximated far field to the far field from
   a direct computation.
*)

let test_dipoles_sculpture =
  let nr=100 in
    Array.init nr
      (fun n ->
	 let s = (float_of_int n)/.(float_of_int nr) in
	   {
	     dp_pos=[|3.0*.cos(8.0*.s)+.5.0*.s*.s;
		      2.0*.sin(3.5*.(s+.0.7))-.3.0*.cos(2.4*.s);
		      -1.0+.2.0*.s+.s*.s*.s
		    |];
	     dp_xyz=[|cos(2.2*.(s-.0.3))+.s+.0.1;
		      cos(2.7*.(s-.0.6))-.4.8*.s+.2.2*.s*.s;
		      sin(s*.s+.1.2)+.0.3;
		    |];
	   })
;;


let test_dipoles_cube =
  [|
    {dp_pos=[|-1.0;-1.0;-1.0|];
     dp_xyz=[|-1.0;-1.0;-1.0|];
    };
    {dp_pos=[|-1.0;-1.0; 1.0|];
     dp_xyz=[|-1.0;-1.0; 1.0|];
    };
    {dp_pos=[|-1.0; 1.0;-1.0|];
     dp_xyz=[|-1.0; 1.0;-1.0|];
    };
    {dp_pos=[|-1.0; 1.0; 1.0|];
     dp_xyz=[|-1.0; 1.0; 1.0|];
    };
    {dp_pos=[| 1.0;-1.0;-1.0|];
     dp_xyz=[| 1.0;-1.0;-1.0|];
    };
    {dp_pos=[| 1.0;-1.0; 1.0|];
     dp_xyz=[| 1.0;-1.0; 1.0|];
    };
    {dp_pos=[| 1.0; 1.0;-1.0|];
     dp_xyz=[| 1.0; 1.0;-1.0|];
    };
    {dp_pos=[| 1.0; 1.0; 1.0|];
     dp_xyz=[| 1.0; 1.0; 1.0|];
    };
  |]
;;

let dipoles_octree ?(max_sites_per_cell=20) ?(ell_max=10) dipoles =
  let nr_coeffs = 2*(ell_max+1)*(ell_max+1) in
  let dp_pos dp = dp.dp_pos in
  let dp_fmm dipoles =
    (* let () = Printf.printf "DDD dp_fmm nr_dipoles=%d\n" (Array.length dipoles) in *)
    {
      fmm_dipoles=dipoles;
      fmm_CF_coeffs=Array.make nr_coeffs 0.0;
      fmm_CN_coeffs=Array.make nr_coeffs 0.0;
      fmm_CF_add_local_dipole_contribs=(fun () -> ());
      fmm_CF_translator_to_parent=[||];
      fmm_CN_translator_from_parent=[||];
      fmm_other_boxes_wciocs=[||];
      fmm_other_boxes_wcoiopcs_and_CFCN_translators=[||];
    }
  in
  let octree =
    make_octree
      ~max_sites_per_cell
      ~fun_content:dp_fmm
      ~fun_coords:dp_pos
      dipoles
  in
  let () =
    fmm_setup_octree_CF_computation
      ~max_multipole_order:ell_max
      ~fun_dipole_position:dp_pos
      ~fun_process_dipole_coeffs:
      (fun dp cont ->
	 let coeffs = dp.dp_xyz in
	   cont coeffs.(0) coeffs.(1) coeffs.(2))
      octree
  in
    octree
;;

let octree_cube = dipoles_octree test_dipoles_cube;;
fmm_octree_set_CF_coeffs ~max_multipole_order:10 octree_cube;;
ddd_mp (match octree_cube.oct_contents with | None -> impossible() | Some c -> c.fmm_CF_coeffs);;

(*
walk_octree
  (fun oct ->
     match oct.oct_contents with
       | None -> ()
       | Some c ->
	   Printf.printf "OCT pos=%s edge=%f *** CN coeffs ***\n%s\n"
	     (float_array_to_string oct.oct_center)
	     oct.oct_edge_length
	     (float_array_to_string (c.fmm_CF_coeffs)))
  octree_cube;;
*)

(*
let octree_random = dipoles_octree test_dipoles_random;;

(* On alpha (1 GHz Intel), interpreted OCaml, 1000 sites, ell_max=10:

   timing dipoles_octree test_dipoles_random;;

   Time passed: 72.426755 sec

   100 sites, ell_max=10: 7.72 sec

   So we roughly scale linearly, and should expect this to take about
   10-15 minutes for 10_000 sites on alpha, in the interpreter.

   Not great, but also not too unreasonable.
*)

(*
timing (fmm_octree_set_CF_coeffs ~max_multipole_order:10) octree_random;;

Time passed: 1.691432 sec

That seems not too unreasonable...
*)

let set_dipoles_random () =
  Array.iter
    (fun dp ->
       begin
	 dp.dp_xyz.(0) <- Random.float 1.0;
	 dp.dp_xyz.(1) <- Random.float 1.0;
	 dp.dp_xyz.(2) <- Random.float 1.0;
       end
    )
    test_dipoles_random
;;

set_dipoles_random();;			       
timing (fmm_octree_set_CF_coeffs ~max_multipole_order:10) octree_random;;

ddd_mp (match octree_random.oct_contents with | None -> impossible() | Some c -> c.fmm_CF_coeffs);;
*)

(*
walk_octree
  (fun oct ->
     match oct.oct_contents with
       | None -> ()
       | Some c ->
	   Printf.printf "OCT pos=%s edge=%f dipoles=%d *** CN coeffs ***\n%s\n"
	     (float_array_to_string oct.oct_center)
	     oct.oct_edge_length (Array.length c.fmm_dipoles)
	     (float_array_to_string (Array.sub (c.fmm_CF_coeffs) 0 10)))
  octree_random;;
*)

(* === Testing the "dipole sculpture" === *)

let potential_of_single_dipole dp pos =
  let dist = array_pointwise (-.) pos dp.dp_pos in
  let dir = dp.dp_xyz in
  let r = sqrt(euclidean_len_sq dist) in
  let sprod = scalar_product dist dir in
    sprod /. (r*.r*.r)
;;

let potential_of_dipoles v_dp pos =
  Array.fold_left (fun sf x -> sf +. potential_of_single_dipole x pos) 0.0 v_dp
;;

let octree_sculpture = dipoles_octree test_dipoles_sculpture;;
fmm_octree_set_CF_coeffs ~max_multipole_order:10 octree_sculpture;;

let vivified_sculpture = (* played too much "nethack", eh? *)
  let (_,vvf_far) = mp_vivificators 10 (match octree_sculpture.oct_contents with | None -> impossible () | Some c -> c.fmm_CF_coeffs) in
    (fun pos ->
       let rel_pos = array_pointwise (-.) pos octree_sculpture.oct_center in
	 vvf_far rel_pos)
;;

let direct_potential_sculpture = potential_of_dipoles test_dipoles_sculpture;;

let test_sculpture pos = (direct_potential_sculpture pos,vivified_sculpture pos);;

(*

# test_sculpture [|5.0;7.0;10.0|];;

- : float * Complex.t =
(0.218917646083779144,
 {Complex.re = 0.218918424883427248;
  Complex.im = 1.30325548055730718e-17})

*)

