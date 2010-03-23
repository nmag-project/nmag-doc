(* (C) 2006 Dr. Thomas Fischbacher

   An implementation of the 3d fast multipole method

 *)

#use "topfind";;
#require "snippets";;

open Snippets;;

(* XXX These are to go into snippets.ml! XXX *)

let int_power x n =
  let rec walk have power_x rest_n =
    if rest_n=0 then have
    else
      let next_power_x=power_x*.power_x in
	walk
	  (if rest_n land 1 <>0 then power_x*.have else have)
	  next_power_x
	  (rest_n lsr 1)
  in
    if n>=0 then walk 1.0 x n
    else 1.0/.(walk 1.0 x (-n))
;;

let binomial n m = (* Note: we are using float here! *)
  let rec walk n_todo m_todo now =
    if m_todo>m then now
    else walk (n_todo-.1.0) (m_todo+.1.0) (now*.n_todo/.m_todo)
  in 
    walk n 1.0 1.0;;

let eval_poly coeffs x = (* Horner Scheme *)
  let nr_coeffs = Array.length coeffs in
  let rec walk pos now =
    let next = now +. coeffs.(pos) in
      if pos = nr_coeffs-1 then next
      else walk (1+pos) (x*.next)
  in
    if nr_coeffs=0 then 0.0 else walk 0 0.0
;;

let derivative_of_poly coeffs =
  let highest_power=(Array.length coeffs)-1 in
    Array.init highest_power
      (fun j -> 
	 let pre_factor = highest_power-j in
	   (float_of_int pre_factor)*.coeffs.(j))
;;

(* XXX this actually is not needed right now! *)
let array_memoized f =
  let r_arr = ref [||] in
    (fun n ->
       let nr_known = Array.length !r_arr in
	 if n < nr_known
	 then !r_arr.(n)
	 else
	   let arr_new =
	     Array.init (1+n) (fun k -> if k<nr_known then !r_arr.(k) else f k)
	   in
	   let () = r_arr:=arr_new in
	     arr_new.(n))
;;

(* XXX unused as well *)
let add_polys p1_coeffs p2_coeffs =
  let len1 = Array.length p1_coeffs
  and len2 = Array.length p2_coeffs
  in
  let len = max len1 len2 in
    Array.init len
      (fun n -> 
	 let power = len-n in
	 let ix1 = len1-power 
	 and ix2 = len2-power in
	   (if ix1<0 then 0.0 else p1_coeffs.(ix1))
	   +.(if ix2<0 then 0.0 else p2_coeffs.(ix2)))
;;


(* === *)

(* One of the problems we are having with legendre_P_lm is that for reasonable L>0, m<0,
   we cannot really evaluate this function on the z-axis, as this would mean taking
   (1-1^2)**(something negative). Unfortunately, the result there is not just zero - 
   the particular limit really does matter!
*)

let legendre_P_lm =
  let rec legendre_P_lm_0 (ell,m) =
    if m<0
    then
      let neg_legendre = legendre_P_lm_0 (ell,-m) in
	fun z -> (neg_legendre z) *.
	  (float_factorial (ell+m) /. float_factorial (ell-m))
	  *.(if m land 1 = 0 then 1.0 else -1.0)
    else
      let f_ell = float_of_int ell in
      let coeff = 1.0/.((int_power 2.0 ell)*.(float_factorial ell)) in
      let x2_minus_1_pow_L =
	Array.init (2*ell+1)
	  (fun n ->
	     let power = 2*ell-n in
	       if power land 1 <>0 then 0.0
	       else 
		 let sign = if (power land 2 = 0) = (ell land 1 = 0)
		   then 1.0 else -1.0
		 in
		   sign*.(binomial f_ell (float_of_int (power/2)))
	  )
      in (* Actually validated that x2_minus_1_pow_L really is correct. *)
      let x2_poly = iterate (ell+m) derivative_of_poly x2_minus_1_pow_L in
      let scaled_x2_poly = Array.map (fun x -> x *.coeff) x2_poly in
	fun x ->
	  (int_power (sqrt(1.0-.x*.x)) m)*.(eval_poly scaled_x2_poly x)
  in
    memoized (Hashtbl.create 10) legendre_P_lm_0
;;

(* Check - Cf. http://en.wikipedia.org/wiki/Associated_Legendre_polynomials

# let x=0.3 in (7.0*.x*.x-.1.0)*.(1.0-.x*.x)*.7.5;;
- : float = -2.52524999999999977
# legendre_P_lm (4,2) 0.3;;
- : float = -2.52525

# legendre_P_lm (4,1) 0.3;;
- : float = -1.69562693051862068
# let x=0.3 in (7.0*.x*.x*.x-.3.0*.x)*.(sqrt(1.0-.x*.x))*.(-2.5);;
- : float = 1.69562693051862068

Considering our differing sign definitions for P_lm, this is right.

*)



(* In the following, when we use spherical harmonics,
   we parametrize them by z=cos(theta) as well as phi.

   Why use phi and not, say, cos(phi)? Answer:
   if we used cos(phi), we would be numerically bad close to
   the zeroes of sin(phi) and vice versa. We could pass
   both at the same time, which would be better, but instead,
   we just settle for using atan2 and phi angles, even if this
   means that we have to go through some possibly unnecessary
   trigonometry...
*)

(* While we actually do not need the unit-normalized Y_lm, but only
   the F_lm and N_lm far-field and near-field tensor functions, it is
   nevertheless nice and useful to also have the Y_lm available...
*)

let y_lm ell m = (* XXX Not tested! *)
  let p_lm = legendre_P_lm (ell,m) in
  let f_l = float_of_int ell and f_m = float_of_int m in
  let coeff =
    (if m land 1 = 0 then 1.0 else -1.0)
    *. sqrt(((2.0*.f_l+.1.0)/.(4.0*.pi))*.
	      (float_factorial (ell-m))/.(float_factorial (ell+m)))
  in 
    fun z phi ->
      Complex.mul
	{Complex.re=coeff*.(p_lm z);Complex.im=0.0}
	{Complex.re=cos(f_m*.phi);Complex.im=sin(f_m*.phi)}
;;

(*   Note: these satisfy F(l,-m) = (-1)^M F(l,m)*   *)

let _rzp0 r z phi = Complex.zero;;

let mp_far_field ell m =
  if ell < 0 then _rzp0 else if abs(m)>ell then _rzp0 else
      let p_lm = legendre_P_lm (ell,m) in
      let f_m = float_of_int m in
	fun r z phi ->
	  let re_factors = (float_factorial (ell-m))*.(p_lm z)*.(r**(float_of_int (-1-ell))) in
	    {Complex.re = re_factors*.cos(f_m*.phi); Complex.im = re_factors*.sin(f_m*.phi)}
;;

let mp_near_field ell m =
  if ell < 0 then _rzp0 else if abs(m)>ell then _rzp0 else
      let p_lm = legendre_P_lm (ell,m) in
      let f_m = float_of_int m in
	fun r z phi ->
	  let re_factors = (p_lm z)*.(r**(float_of_int ell))/.(float_factorial (ell+m)) in
	  let result =
	    {Complex.re = re_factors*.cos(f_m*.phi);
	     Complex.im = -.re_factors*.sin(f_m*.phi)
	    } in
	    result
;;

(*
mp_near_field 1 (-1) 1.0 (-1.0) 3.1415926535897932384626433832;;
*)

(* We need translation of array indices to lexicographical pairs (L,m):

   0         1         2         3         4
   (L=0,m=0) (L=1,m=0) (L=1,m=1) (L=2,m=0) (L=2,m=1) ...

   We do this by lazily building a lookup array:
 *)

let _array_Lm = ref [||];;

let _ensure_array_Lm_is_long_enough max_ix =
  if max_ix < Array.length (!_array_Lm) then ()
  else
    let rec walk_ix_Lm have just_finish_this_ell ix ell m =
      if m > ell then
	if just_finish_this_ell
	then Array.of_list (List.rev have)
	else
	  walk_ix_Lm have (ix>=max_ix) ix (1+ell) (-1-ell)
      else
	walk_ix_Lm ((ell,m)::have) (ix>=max_ix) (1+ix) ell (m+1)
    in _array_Lm := walk_ix_Lm [] (0>=max_ix) 0 0 0
;;
  

let mp_linear_index_to_Lm ix =
  let () = _ensure_array_Lm_is_long_enough ix in
  (!_array_Lm).(ix)
;;

let mp_Lm_to_linear_index ell m =	
  ell*(ell+1)+m;;
  

(* XXX speedup possible through separation of real and imaginary parts!
 *)
let mp_vivificators ell_max =
  let rec walk_ix_lm have ix ell m =
    if ell > ell_max then List.rev have
    else if m > ell then walk_ix_lm have ix (1+ell) (-1-ell)
    else
      let ff = mp_far_field ell m
      and nf = mp_near_field ell m
      in
	walk_ix_lm
	  ((ff,nf)::have) (1+ix) ell (m+1)
  in
  let vvf = Array.of_list (walk_ix_lm [] 0 0 0)
  in
  let nr_vvf_entries = Array.length vvf in
  let vivify want_near v_mp_coeffs r =
    let len_r = sqrt(euclidean_len_sq r) in
    let z = r.(2)/.len_r in
    let phi = atan2 r.(1) r.(0) in
    let rec walk now pos =
      if pos = nr_vvf_entries then now
      else
	let coeff_re = v_mp_coeffs.(pos+pos)
	and coeff_im = v_mp_coeffs.(pos+pos+1)
	in
	let vvf_fun = 
	  if want_near
	  then let (_,x) = vvf.(pos) in x
	  else let (x,_) = vvf.(pos) in x
	in
	walk (Complex.add now
		(Complex.mul
		   {Complex.re=coeff_re;Complex.im=coeff_im}
		   (vvf_fun len_r z phi)))
	  (1+pos)
    in walk {Complex.re=0.0;Complex.im=0.0} 0
  in 
  fun mp_coeffs -> (vivify true mp_coeffs, vivify false mp_coeffs)
;;

(* The multipole coefficients of an elementary dipole, at the origin *)

let mp_coeffs_dipole ell_max v_dipole =
  let nr_coeffs = 2*(ell_max+1)*(ell_max+1) (* real plus imaginary entries *)
  in
  let a = Array.make nr_coeffs 0.0 in
  begin
    a.(2) <- -0.5*.v_dipole.(0); (* cos(m*phi) part of (l,m)=(1,-1) *)
    a.(3) <- -0.5*.v_dipole.(1); (* sin(m*phi) part of (l,m)=(1,-1) *)
    a.(4) <- v_dipole.(2); (* cos(m*phi) part of (l,m)=(1,0) *)
    a.(5) <- 0.0; (* sin(m*phi) part of (l,m)=(1,0) - irrelevant as m=0! *)
    a.(6) <- 0.5*.v_dipole.(0); (* cos(m*phi) part of (l,m)=(1,1) *)
    a.(7) <- -0.5*.v_dipole.(1); (* sin(m*phi) part of (l,m)=(1,1) *)
    a
  end
;;

(* The translator matrices are sparse - we represent them as a vector of sparse rows,
   where every row consists of a vector of (column,coeff_re,coeff_im).

   What may our signs depend on? When do we flip a sign?

   always
   src_m odd
   result_m odd
   src_m > result_m

   So, we get 4 independent sources for sign flips, and these have
   to be used at 8 places where we have signs.

 *)

let mp_make_translator
    ?(debug=false)
    ?(sign_config=Array.make 8 (Array.make 4 true))
    ell_max =
  let nr_coeffs = mp_Lm_to_linear_index (1+ell_max) (-1-ell_max) in
    (fun f_L1m1_L2m2 dist ->
       let r = sqrt(euclidean_len_sq dist) in
       let z = if r=0.0 then 0.0 else dist.(2)/.r in
	 (* Note that we have to prevent division by zero in the r=0 case! *)
       let phi = atan2 dist.(1) dist.(0) in
	 Array.init nr_coeffs
	   (fun ix_lm_result ->
	      let (result_ell,result_m) as result_Lm = mp_linear_index_to_Lm ix_lm_result in
	      let rec make_row have ix_Lm_src =
		let (src_ell,src_m) as src_Lm = mp_linear_index_to_Lm ix_Lm_src in
		  if src_ell > ell_max
		  then Array.of_list (List.rev have)
		  else
		    let coeff_L_m=f_L1m1_L2m2 result_Lm src_Lm r z phi in
		      if Complex.norm coeff_L_m = 0.0 
		      then make_row have (1+ix_Lm_src) (* no entry*)
		      else
			make_row ((ix_Lm_src,coeff_L_m)::have) (1+ix_Lm_src)
	      in
		make_row [] 0))
;;

(* XXX OBSOLETE XXX
let apply_translator translator v_mp =
  let nr_coeffs = Array.length v_mp in
  Array.init nr_coeffs
    (fun n ->
      let is_re = (n land 1 = 0)
      and ix_row = n lsr 1
      in
      let translator_row = translator.(ix_row) in
      let row_len = Array.length translator_row in
      let rec walk_row have pos =
	if pos=row_len then have
	else
	  let (nr_Lm,coeff) = translator_row.(pos) in
	    walk_row
	      (
		if is_re
		then (have+.(coeff.Complex.re)*.v_mp.(nr_Lm+nr_Lm) -.(coeff.Complex.im)*.v_mp.(nr_Lm+nr_Lm+1))
		else (have+. (coeff.Complex.im)*.v_mp.(nr_Lm+nr_Lm) +. (coeff.Complex.re)*.v_mp.(nr_Lm+nr_Lm+1))
	      )
	    (1+pos)
      in walk_row 0.0 0)
;;
*)

let apply_translator_add_contribs v_result translator v_mp =
  let nr_coeffs = Array.length v_mp in
    for n=0 to nr_coeffs-1 do
      let contrib =
	let is_re = (n land 1 = 0)
	and ix_row = n lsr 1
	in
	let translator_row = translator.(ix_row) in
	let row_len = Array.length translator_row in
	let rec walk_row have pos =
	  if pos=row_len then have
	  else
	    let (nr_Lm,coeff) = translator_row.(pos) in
	      walk_row
		(
		  if is_re
		  then (have+.(coeff.Complex.re)*.v_mp.(nr_Lm+nr_Lm) -.(coeff.Complex.im)*.v_mp.(nr_Lm+nr_Lm+1))
		  else (have+. (coeff.Complex.im)*.v_mp.(nr_Lm+nr_Lm) +. (coeff.Complex.re)*.v_mp.(nr_Lm+nr_Lm+1))
		  )
		(1+pos)
	in walk_row 0.0 0
      in
	v_result.(n) <- v_result.(n) +. contrib
    done
;;

let apply_translator translator v_mp =
  let v_result = Array.make (Array.length v_mp) 0.0 in
  let () = apply_translator_add_contribs v_result translator v_mp in
    v_result
;;


(* Again, we are facing the problem that we unnecessarily do duplicate some work.
   (Computing the (L,m) -> index hashes in make_translator).

   But hey - this is our first attempt at the fast multipole method anyway!
   There will be occasion to improve later on!
 *)

let mp_translator_ff ?(debug=false) ell_max dist =
  mp_make_translator ~debug ell_max
    (fun (ell_dst,m_dst) (ell_src,m_src) ->
       mp_near_field (ell_dst-ell_src) (m_dst-m_src))
    (Array.map (fun x -> -.x) dist)
;;

let mp_translator_nn ?(debug=false) ell_max dist =
  mp_make_translator ~debug ell_max
    (fun (ell_dst,m_dst) (ell_src,m_src) ->
       mp_near_field (ell_src-ell_dst) (m_src-m_dst))
    dist
;;


let mp_translator_nf ?(debug=false) ell_max dist =
  mp_make_translator ~debug ell_max
    (fun (ell_dst,m_dst) (ell_src,m_src) ->
       fun r z phi ->
	 let coeff = mp_far_field (ell_dst+ell_src) (m_dst+m_src) r z phi in
	   if ell_dst land 1 <> 0
	   then {Complex.re= -.coeff.Complex.re;Complex.im= -.coeff.Complex.im}
	   else coeff)
    dist
;;

(* We also need multipole coefficients of a shifted dipole, as this is the fundamental
   contribution we will be encountering.
*)

let mp_add_ff_coeffs_of_shifted_dipole ell_max dist =
  let reduced_xlator =
    Array.map
      (fun x ->
	 array_filter
	   (fun (ix_src,_) ->
	      let (ell,_) = mp_linear_index_to_Lm ix_src in ell=1)
	   x)
      (mp_translator_ff ell_max dist)
  in
    fun v_result coeff_x coeff_y coeff_z ->
      for i=0 to Array.length reduced_xlator-1 do
	let contribs = reduced_xlator.(i) in
	  for k=0 to Array.length contribs-1 do
	    let (ix_src,coeff) = contribs.(k) in
	    let (ell_src,m_src) = mp_linear_index_to_Lm ix_src in
	      (* ell_src = 1 as we are using reduced_xlator! *)
	      if m_src=0 then (* using coeff_z *)
		begin
		  v_result.(i+i) <- v_result.(i+i) +. coeff.Complex.re*.coeff_z;
		  v_result.(i+i+1) <- v_result.(i+i+1) +. coeff.Complex.im*.coeff_z;
		end
	      else if m_src=1 then (* using coeff_x - i*coeff_y   XXX check sign *)
		begin
		  v_result.(i+i) <-
		    v_result.(i+i)
		  +. 0.5*.(coeff.Complex.re*.coeff_x +. coeff.Complex.im*.coeff_y);
		  v_result.(i+i+1) <-
		    v_result.(i+i+1)
		  +. 0.5*.(coeff.Complex.im*.coeff_x -. coeff.Complex.re*.coeff_y);
		end
	      else
		begin
		  v_result.(i+i) <-
		    v_result.(i+i)
		  +.0.5*.(-. coeff.Complex.re*.coeff_x +. coeff.Complex.im*.coeff_y);
		  v_result.(i+i+1) <-
		    v_result.(i+i+1)
		  +.0.5*.(-. coeff.Complex.im*.coeff_x -. coeff.Complex.re*.coeff_y);
		end
	  done
	  done
;;





let ddd_mp v_coeffs =
  let nr_coeffs = Array.length v_coeffs in
  let rec walk ell m pos =
    if m > ell then walk (1+ell) (-1-ell) pos
    else
      if pos+pos=nr_coeffs then ()
      else
	let re=v_coeffs.(pos+pos)
	and im=v_coeffs.(pos+pos+1)
	in
	begin
	  Printf.printf "(l,m)=(%d,%d) -- Coeff-Re: %f Coeff-Im: %f\n" ell m re im;
	  walk ell (1+m) (pos+1)
	end
  in walk 0 0 0
;;

let ddd_xlator xlator =
  Array.iteri
    (fun ix_Lm row ->
       let (ell_result,m_result)= mp_linear_index_to_Lm ix_Lm in
       let () = Printf.printf "\n*** (L,m) = (%2d,%2d) ***\n" ell_result m_result in
	 Array.iter
	   (fun (ix_src,c) ->
	      let (ell_src,m_src)= mp_linear_index_to_Lm ix_src in
		Printf.printf "  (%2d,%2d):  #c(%f %f)\n" ell_src m_src c.Complex.re c.Complex.im) 
	 row)
    xlator
;;

(* === Testing Code === *)

(* For now, we only go up to the octupole (L=3) moment *)

let test_ell_max = 3;;

(* Note that our values right now DO NOT converge with increasing max multipole order! *)

let check_mp ?(ell_max=3) v_mp dist0 dist1 dist2 =
  let dist12 = array_pointwise (+.) dist1 dist2 in
  let mp_coeffs =
    apply_translator
      (mp_translator_ff ell_max dist0) (mp_coeffs_dipole ell_max v_mp)
  in
  let mp_xlator_1 = mp_translator_ff ell_max dist1
  and mp_xlator_2 = mp_translator_ff ell_max dist2
  and mp_xlator_12 = mp_translator_ff ell_max dist12
  in
  let xlated_1 = apply_translator mp_xlator_1 mp_coeffs in
  let xlated_2_xlated_1 = apply_translator mp_xlator_2 xlated_1 in
  let xlated_12 = apply_translator mp_xlator_12 mp_coeffs in
    (sqrt(euclidean_distance_sq xlated_2_xlated_1 xlated_12)
     (* ,xlated_2_xlated_1,xlated_12 *) 
    )
;;

check_mp ~ell_max:10 [|2.0;5.0;1.0|] [|0.3;2.1;0.7|] [|3.0;4.0;-2.0|] [|1.2;2.1;4.7|];;
(* This gives:
- : float = 7.65661525105756841e-12

   I take this as a strong indication that this implementation of the algebra does work!
*)

let test_dipole = mp_coeffs_dipole 10 [|0.0;0.0;0.0|];;
mp_add_ff_coeffs_of_shifted_dipole 10 [|1.0;1.3;1.7|] test_dipole 0.5 0.6 0.9;;

sqrt(euclidean_distance_sq
       (apply_translator (mp_translator_ff 10 [|1.0;1.3;1.7|]) (mp_coeffs_dipole 10 [|0.5;0.6;0.9|]))
       test_dipole);;

(* This gives zero, hence this is also fine! *)

