;; (C) 2006 Dr. Thomas Fischbacher
;;
;; Some code to get the multipole translation algebra right...
;;
;; Design-wise, this is slow-and-stupid code: do not try to be clever,
;; but try to capture the idea in such a way that every step is
;; checkable and reproducible!

;; First, the associated Legendre polynomials:

;; (eval-when (compile load eval)
;;   (progn
;;     (clc:clc-require :tf-spellbook)
;;     (use-package :tf-spellbook)))

(defun fakt (n)
  (labels
      ((walk (have nn)
	     (if (= nn n) have
	       (walk (* have (1+ nn)) (1+ nn)))))
    (walk 1 0)))

(defmacro hv (ht key &rest args)
  `(gethash ,key ,ht ,@args))

(defmacro mv-bind (&rest stuff)
  `(multiple-value-bind ,@stuff))

(defun hv-inc* (ht key &optional (increment 1))
  (mv-bind (old success) (gethash key ht)
    (let ((new (+ increment (if success old 0))))
      (setf (gethash key ht) new))))
	  

(defun v-map (f &rest args)
  (apply #'map '(simple-array * (*)) f args))

(defun cav (av x &rest args)
  (cdr (apply #'assoc x av args)))

(defun binomial (n k) ; utterly primitive
  (/ (fakt n) (fakt k) (fakt (- n k))))

(defun horner-eval (v-poly x)
  (let ((nr-coeffs (length v-poly)))
    (labels
	((walk (pos now)
           (let ((next (+ now (aref v-poly pos))))
	     (if (= pos (1- nr-coeffs))
		 next
	       (walk (1+ pos) (* x next))))))
      (if (= nr-coeffs 0) 0 (walk 0 0)))))

(defun derivative-of-poly (v-coeffs)
  (let* ((highest-power (1- (length v-coeffs)))
	 (result (make-array highest-power :initial-element 0)))
    (dotimes (j highest-power)
      (setf (aref result j)
	    (* (- highest-power j) (aref v-coeffs j))))
    result))

;; This could use memoization - for now, we do not use that!

(defun church (n f x)
  (labels
      ((walk (nn now)
	 (if (= n nn) now
	   (walk (1+ nn) (funcall f now)))))
    (walk 0 x)))


(defun legendre-P-lm (L m)
  (let* ((coeff (/ (* (expt 2 L) (fakt L))))
	 (x^2-1**L
	  (let ((a (make-array (1+ (* 2 L)))))
	    (dotimes (n (length a))
	      (setf (aref a n)
		    (let ((power (- (+ L L) n)))
		      (if (oddp power) 0
			(let ((sign (if (equal (evenp (/ power 2)) (evenp L)) 1 -1)))
			  (* sign (binomial L (/ power 2))))))))
	    a))
	 (x^2-poly (v-map #'(lambda (x) (* x coeff)) (church (+ L m) #'derivative-of-poly x^2-1**L))))
    #'(lambda (x)
	(* (expt (- 1 (* x x)) (/ m 2.0d0))
	   (horner-eval x^2-poly x)))))

(defun Y-lm (L m theta phi)
  (* (sqrt (* (/ (* 4 pi) (+ L L 1)) (/ (fakt (- L m)) (fakt (+ L m)))))
     (funcall (legendre-P-lm L m) (cos theta))
     (cis (* m phi))))

;; For a few test (L,m) and x, this legendre-P-lm very nicely matches
;; our ML implementation. So, we do the same thing up to here...

(defun mp-far-field (L m r z phi)
  (* (fakt (- L m)) (expt r (- -1 L))
     (funcall (legendre-P-lm L m) z) ; XXX determination of the P-lm should be memoized or moved out!
     (cis (* m phi))))

(defun mp-near-field (L m r z phi)
  (* (/ (fakt (+ L m))) (expt r L)
     (funcall (legendre-P-lm L m) z) ; XXX determination of the P-lm should be memoized or moved out!
     (cis (* -1 m phi))))


(defun mp-translator (L-max f-LmLm)
  (labels
      ((walk-Lm-result (have nr-result L-result m-result)
          (cond
	   ((> L-result L-max) (reverse have))
	   ((> m-result L-result)
	    (walk-Lm-result have nr-result (1+ L-result) (- (1+ L-result))))
	   (t
	    (let ((next-row
		   `((nLm-result #|,nr-result|# ,L-result ,m-result) .
		     ,(labels
			  ((walk-Lm-src (have nr-src L-src m-src)
			     (cond
			      ((> L-src L-max) (reverse have))
			      ((> m-src L-src) (walk-Lm-src have nr-src (1+ L-src) (- (1+ L-src))))
			      (t
			       (let ((coeff (funcall f-LmLm L-result m-result L-src m-src)))
				 (if (= coeff 0)
				     (walk-Lm-src have (1+ nr-src) L-src (1+ m-src))
				   (walk-Lm-src `(((#|,nr-src|# ,L-src ,m-src) . ,coeff) . ,have)
						(1+ nr-src)
						L-src (1+ m-src))))))))
			(walk-Lm-src nil 0 0 0)))))
	      (walk-Lm-result (cons next-row have) (1+ nr-result) L-result (1+ m-result)))))))
    (walk-Lm-result nil 0 0 0)))


(defun mp-translator-ff (L-max dist)
  (let* ((dist (v-map #'- dist))
	 (r (sqrt (reduce #'(lambda (sf x) (+ sf (* x x))) dist :initial-value 0)))
	 (z (/ (aref dist 2) r))
	 (phi (atan (aref dist 1) (aref dist 0)))
	 (f-LmLm
	  #'(lambda (L1 m1 L2 m2)
	      (let ((L1-L2 (- L1 L2))
		    (m1-m2 (- m1 m2)))
		(if (or (< L1-L2 0) (> (abs m1-m2) L1-L2))
		    0
		  (mp-near-field (- L1 L2) (- m1 m2) r z phi))))))
    (mp-translator L-max f-LmLm)))

(defun mp-dipole (v-dipole)
  `(
    ((1 -1) . ,(complex (- (aref v-dipole 0)) (aref v-dipole 1))) ; -x+iy
    ((1  0) . ,(aref v-dipole 2)) ; z component
    ((1 +1) . ,(complex (aref v-dipole 0) (aref v-dipole 1))) ; x+iy
    ))

;; Note: there may be issues with factors of 2 in the (1,1) and (1,-1)
;; components!

(defun apply-translator (xlator mp)
  (mapcan
   #'(lambda (xlator-row)
       (let ((result-Lm (cdar xlator-row))
	     (contribs (cdr xlator-row)))
	 (let ((result-coeff
		(reduce
		 #'(lambda (sf spec)
		     (let ((mp-coeff (cav mp (car spec) :test #'equal)))
		       (if mp-coeff (+ sf (* mp-coeff (cdr spec))) sf)))
		 xlator-row :initial-value 0)))
	   (if (< (abs result-coeff) 1d-10)
	       nil
	     `((,result-Lm . ,result-coeff))))))
   xlator))

(defun mp-dist (mp1 mp2)
  (let ((h (make-hash-table :test 'equal))
	(d 0d0))
    (dolist (c mp1) (setf (hv h (car c)) (cdr c)))
    (dolist (c mp2) (hv-inc* h (car c) (- (cdr c))))
    (maphash
     #'(lambda (k v)
	 (setf d (+ d (abs (* v v)))))
     h)
    (sqrt d)))

(defun mp-print-translator (translator)
  (dolist (row translator)
    (format t "~%*** (L,m)=(~3d,~3d) ***~%" (cadar row) (caddar row))
    (dolist (entry (cdr row))
      (let* ((coeff (cdr entry))
	     (s-coeff
	      (if (= coeff (realpart coeff))
		  (format nil "~10,8F" (realpart coeff))
		(format nil "~10,8F + ~10,8F i"
			(realpart coeff)
			(imagpart coeff)))))
	(format t "  (~3d,~3d) => ~10A~%" (caar entry) (cadar entry) s-coeff)))))


			

(defparameter dipole-1 (mp-dipole #(-3.2d0 1.9d0 0.7d0)))
(defparameter xlate-1 #(4.2d0 2.8d0 -5.1d0))
(defparameter xlate-2 #(1.3d0 1.1d0 2.4d0))
(defparameter xlator-1 (mp-translator-ff 10 xlate-1))
(defparameter xlator-2 (mp-translator-ff 10 xlate-2))
(defparameter xlator-3 (mp-translator-ff 10 (v-map #'+ xlate-1 xlate-2)))

(mp-dist (apply-translator xlator-3 dipole-1)
	 (apply-translator xlator-2 (apply-translator xlator-1 dipole-1)))

;; This gives 2.0160075478027023d-11
;;
;; As our example was totally generic,
;; we can presumably take this as a strong indication that our multipole
;; translation algebra is implemented correctly...

(defparameter dipole-x (mp-dipole #(1 0 0)))
(defparameter xlator-x (mp-translator-ff 3 #(10 0 0)))
(defparameter xlator-2x (mp-translator-ff 3 #(20 0 0)))
;; (apply-translator xlator-x dipole-x)

(defparameter dipole-y (mp-dipole #(0 1 0)))
(defparameter xlator-y (mp-translator-ff 3 #(0 10 0)))

(defparameter xlated-y (apply-translator xlator-y dipole-y))

;;; ================
