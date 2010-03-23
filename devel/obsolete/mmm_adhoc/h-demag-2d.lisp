;;;; (C) 2006 Dr. Thomas Fischbacher
;;;;
;;;; Quick ad-hoc lisp implementation of the computation of the
;;;; dense M -> H matrix.
;;;;
;;;; We assume a first-order 2d mesh that describes a thin film.
;;;;
;;;; Input:
;;;;
;;;; - array of point coordinates
;;;; - vector of triangles (index-wise)
;;;; - subdivision level
;;;;
;;;; Output: matrix - 3*nr-sites x 3*nr-sites
;;;; (for every site, we map M to H).
;;;;
;;;; It seems as if we do have a problem here: the task of computing
;;;; H_demag for an infinitely thin film with given surface
;;;; magnetization density seems to be somewhat ill-defined.
;;;;
;;;; Can we see this? I think we can. Let us just try to integrate the effect of
;;;; a constant density of z-directed multipoles at a given point in the plane.
;;;;
;;;; We get Integral dr*2pi*r* m_z/r^3 ~ like Integral dr/r^2, from 0 to infinity.
;;;;
;;;; Evidently, this diverges at the lower bound r=0.
;;;;
;;;; So, how do we deal with this? I think this goes away if we
;;;; include the third dimension (film thickness) into our
;;;; considerations: whenever the distance between sub-triangles
;;;; gets too small wrt to their height, we also subdivide the z-direction.


(eval-when (compile load eval)
  (progn
    (clc:clc-require '#:tf-spellbook)
    (use-package '#:tf-spellbook)
    ))



(defun make-h-demag-matrix (points triangles &key
				   (subdivision-level 3)
				   (film-thickness 0.1d0)
				   ; (max-skew 0.03d0) ; 3% error in r makes 10% in 1/r^3
				   (max-skew 0.01d0) ; 3% error in r makes 10% in 1/r^3
				   )
  (declare (optimize (safety 0) (speed 3)))
  ;; (declare (optimize (safety 3) (speed 0)))
  (let* ((film-thickness (coerce film-thickness 'double-float))
	 (max-skew (coerce max-skew 'double-float))
	 (min-r^2-xy
	  (let ((r (/ max-skew film-thickness)))
	    (* r r)))
	 (points (map '(simple-array (simple-array double-float (2)) (*))
		      #'(lambda (p) (map '(simple-array double-float (2))
					 #'(lambda (x) (coerce x 'double-float)) p))
		      points))
	 (triangles
	  (map '(simple-array (simple-array fixnum (3)) (*))
	       #'(lambda (tri)
		   (map '(simple-array fixnum (3)) #'identity tri))
	       triangles))
	 (nr-sites (array-dimension points 0))
	 (nr-dofs (i* 3 nr-sites))
	 (result (make-array `(,nr-dofs ,nr-dofs)
			     :element-type 'double-float
			     :initial-element 0.0d0))
	 (scratch-res (make-array 3
				  :element-type 'double-float
				  :initial-element 0.0d0))
	 (scratch-src (make-array 3
				  :element-type 'double-float
				  :initial-element 0.0d0))
	 (nr-triangles (array-dimension triangles 0))
	 )
    (declare (type (simple-array double-float (*)) scratch-res scratch-src)
	     (type (simple-array double-float (* *)) result)
	     (fixnum nr-triangles nr-sites nr-dofs)
	     (double-float film-thickness max-skew min-r^2-xy)
	     )
    (labels
	(
	 (triangle-area (p0 p1 p2)
	   (declare (type (simple-array double-float (2)) p0 p1 p2))
	   (let ((d1x (f- (aref p1 0) (aref p0 0)))
		 (d1y (f- (aref p1 1) (aref p0 1)))
		 (d2x (f- (aref p2 0) (aref p0 0)))
		 (d2y (f- (aref p2 1) (aref p0 1))))
	     (f* 0.5d0 (f- (f* d1x d2y) (f* d1y d2x)))))
	 ;;
	 (center2 (L0 L1)
	  (declare (type (simple-array double-float (3)) L0 L1))
	  (let ((a (make-array 3
			       :element-type 'double-float
			       :initial-element 0.0d0)))
	    (setf (aref a 0) (f* 0.5d0 (f+ (aref L0 0) (aref L1 0))))
	    (setf (aref a 1) (f* 0.5d0 (f+ (aref L0 1) (aref L1 1))))
	    (setf (aref a 2) (f* 0.5d0 (f+ (aref L0 2) (aref L1 2))))
	    a))
	 ;;
	 (center3 (L0 L1 L2)
	  (declare (type (simple-array double-float (3)) L0 L1 L2))
	  (let ((a (make-array 3
			       :element-type 'double-float
			       :initial-element 0.0d0)))
	    (setf (aref a 0) (f/ (f+ (aref L0 0) (aref L1 0) (aref L2 0)) 3.0d0))
	    (setf (aref a 1) (f/ (f+ (aref L0 1) (aref L1 1) (aref L2 1)) 3.0d0))
	    (setf (aref a 2) (f/ (f+ (aref L0 2) (aref L1 2) (aref L2 2)) 3.0d0))
	    a))
	 ;;
	 (instantiate (target p0 p1 p2 L)
	   (declare (type (simple-array double-float (*)) target p0 p1 p2 L))
	   (setf (aref target 0)
		 (f+ (f* (aref L 0) (aref p0 0))
		     (f* (aref L 1) (aref p1 0))
		     (f* (aref L 2) (aref p2 0))))
	   (setf (aref target 1)
		 (f+ (f* (aref L 0) (aref p0 1))
		     (f* (aref L 1) (aref p1 1))
		     (f* (aref L 2) (aref p2 1))))
	   nil)
	 ;;
	 (distance^2 (p q)
	   (let ((dx (- (aref q 0) (aref p 0)))
		 (dy (- (aref q 1) (aref p 1)))
		 (dz (- (aref q 2) (aref p 2)))
		 )
	     (f+ (f* dx dx) (f* dy dy)  (f* dz dz))))
	 ;;
	 (do-triangle (L0 L1 L2 subdivide f)
	   (if (= 0 subdivide)
	       (funcall f L0 L1 L2)
	     (let ((L01 (center2 L0 L1))
		   (L12 (center2 L1 L2))
		   (L20 (center2 L2 L0))
		   (new-sd (i1- subdivide)))
	       (do-triangle L0  L01 L20 new-sd f)
	       (do-triangle L01 L1  L12 new-sd f)
	       (do-triangle L20 L12 L2  new-sd f)
	       (do-triangle L01 L12 L20 new-sd f))))
	 )
      (let* ((float-sd-level (coerce subdivision-level 'double-float))
	     (element-areas
	      (let ((part (expt 0.25d0 float-sd-level)))
		(map '(simple-array double-float (*))
		     #'(lambda (tri)
			 (f* part
			     (abs
			      (triangle-area (aref points (aref tri 0))
					     (aref points (aref tri 1))
					     (aref points (aref tri 2))))))
		     triangles)))
	     (mini-triangles
	      (let ((tri nil)
		    (L0 (coerce '(1.0d0 0.0d0 0.0d0) '(simple-array double-float (3))))
		    (L1 (coerce '(0.0d0 1.0d0 0.0d0) '(simple-array double-float (3))))
		    (L2 (coerce '(0.0d0 0.0d0 1.0d0) '(simple-array double-float (3)))))
		(do-triangle L0 L1 L2 subdivision-level
			     #'(lambda (small-L0 small-L1 small-L2)
				 (ppush tri
					(coerce (list small-L0 small-L1 small-L2
						      (center3 small-L0 small-L1 small-L2))
						'(simple-array * (*))))))
		(coerce tri '(simple-array * (*)))))
	     (nr-mini-triangles (array-dimension mini-triangles 0)))
	(declare (type (simple-array * (*)) mini-triangles))
      (i-dotimes (nr-tri-res nr-triangles)
        (let* ((tri-res (aref triangles nr-tri-res))
	       (ix0-res (aref tri-res 0))
	       (ix1-res (aref tri-res 1))
	       (ix2-res (aref tri-res 2))
	       (p0-res (aref points ix0-res))
	       (p1-res (aref points ix1-res))
	       (p2-res (aref points ix2-res)))
	  (declare (type (simple-array fixnum (3)) tri-res)
		   (type (simple-array double-float (2)) p0-res p1-res p2-res))
	  ;; this way, we avoid a closure!
	  (i-dotimes (nr-tri-src nr-triangles)
	     (let* ((tri-src (aref triangles nr-tri-src))
		    (ix0-src (aref tri-src 0))
		    (ix1-src (aref tri-src 1))
		    (ix2-src (aref tri-src 2))
		    (p0-src (aref points ix0-src))
		    (p1-src (aref points ix1-src))
		    (p2-src (aref points ix2-src)))
	       (declare (type (simple-array fixnum (3)) tri-src)
			(type (simple-array double-float (2)) p0-src p1-src p2-src))
	       (i-dotimes (nr-mini-tri-res nr-mini-triangles)
		  (let* ((mini-tri-res (aref mini-triangles nr-mini-tri-res))
			 (L-mid-res (aref mini-tri-res 3)))
		    (declare (type (simple-array double-float (3)) L-mid-res)
			     (type (simple-array * (*)) mini-tri-res))
		    (instantiate scratch-res p0-res p1-res p2-res L-mid-res)
		    (i-dotimes (nr-mini-tri-src nr-mini-triangles)
		       (let* ((mini-tri-src (aref mini-triangles nr-mini-tri-src))
			      (L-mid-src (aref mini-tri-src 3)))
			 (declare (type (simple-array double-float (3)) L-mid-src)
				  (type (simple-array * (*)) mini-tri-src))
			 (instantiate scratch-src p0-src p1-src p2-src (aref mini-tri-src 3))
			 (let* ((r^2-xy (distance^2 scratch-res scratch-src)))
			   (if (> r^2-xy 0.0d0)
			       (let* ((nr-z-levels 
				       (if nil
					   (if (f> r^2-xy min-r^2-xy) 1 ; short-cut
					     (ceiling (sqrt (f/ min-r^2-xy r^2-xy))))
					 (expt 2 subdivision-level) ; XXX TEST!
					 ))
				      ;; (nr-z-levels (i* nr-z-levels nr-z-levels)) ; can this help to get divergent behaviour under better control?
				      (z-step (f/ film-thickness (coerce nr-z-levels 'double-float)))
				      )
				 (do ((nr-z-res 0 (i1+ nr-z-res))
				      (z-res 0.0d0 (f+ z-res z-step)))
				     ((i= nr-z-res nr-z-levels))
				   (declare (fixnum nr-z-res) (double-float z-res))
				 (do ((nr-z-src 0 (i1+ nr-z-src))
				      (z-src 0.0d0 (f+ z-src z-step)))
				     ((i= nr-z-src nr-z-levels))
				   (declare (fixnum nr-z-src) (double-float z-src))
				   (let* ((delta-z (f- z-src z-res))
					  (r^2 (f+ r^2-xy (f* delta-z delta-z)))
					  (1/r (f/ (sqrt r^2)))
					  (1/r^2 (f* 1/r 1/r))
					  (1/r^3 (f* 1/r^2 1/r))
					  (1/r^5 (f* 1/r^3 1/r^2)))
				     ;; the dipole field is:
				     ;; v = p/r^3 - 3(p*r)*r/r^5
				     ;; This gives us 3*3 matrix contributions.
				     ;; Note that these will be further weighted with source
				     ;; and result shape functions (just the L-factors).
				     ;; So... we do this in two steps. First, determine the local
				     ;; matrix coefficient...
				     (i-dotimes (dir-h 3)
					(i-dotimes (dir-m 3)
					   (let ((contrib-h-m
						  (f* (aref element-areas nr-tri-res) z-step
						      (aref element-areas nr-tri-src) z-step
						      (f+ (if (i= dir-m dir-h)
							      1/r^3 0.0d0)
							  (f* -3.0d0
							      (f+ (if (i= 2 dir-m) delta-z 0.0d0)
								  (f- (aref scratch-src dir-m)
								      (aref scratch-res dir-m)))
							      (f+ (if (i= 2 dir-h) delta-z 0.0d0)
								  (f- (aref scratch-src dir-h)
								      (aref scratch-res dir-h)))
							      1/r^5)))))
					     (i-dotimes (ix-L-res 3)
						(i-dotimes (ix-L-src 3)
						   (if (i/= (aref tri-res ix-L-res)
							    (aref tri-src ix-L-src))
						       ;; exclude the diagonal (self-interaction)
						       ;; M/H blocks of the matrix.
						       (incf (aref result 
								   (i+ dir-h (i* 3 (aref tri-res ix-L-res)))
								   (i+ dir-m (i* 3 (aref tri-src ix-L-src))))
							     (f* contrib-h-m
								 (aref L-mid-res ix-L-res)
								 (aref L-mid-src ix-L-src)))))))))
				     ))))))))))))))
      result))))

(eval-when (eval)
  (compile 'make-h-demag-matrix))

;; Evidently, we also have to test all this...
;; We just use a simple nmesh mesh for that...

(defparameter demo-points
  ((lambda (z)
     (v-map
      #'(lambda (p)
	  (map '(simple-array double-float (*))
	       #'(lambda (c) (coerce c 'double-float))
	       p))
      z))
   '((0 0) (0 1) (1 0) (1 1))))

(defparameter demo-triangles
  ((lambda (z)
     (v-map #'(lambda (tri) (coerce tri '(simple-array fixnum (3)))) z))
   '((0 1 2) (1 2 3))))

(defparameter m-lvls (v-map #'(lambda (n) (make-h-demag-matrix demo-points demo-triangles :subdivision-level n :film-thickness 0.4d0)) (v-int-range 4)))

#|

(v-map #'(lambda (m) (aref m 0 3)) M-LVLS)

==> gives:

#(-0.001325825214724775d0 -0.0024254812357840496d0 -0.0064556481588382295d0
  -0.010925781010582029d0)

...does not seem to converge. :-(

...with 2^level z-subdivisions:

#(-0.001325825214724775d0 -0.002176877393217904d0 -0.004677672087769297d0
  -0.008438019229360768d0)

Diverges slower, but still diverges.

|#
