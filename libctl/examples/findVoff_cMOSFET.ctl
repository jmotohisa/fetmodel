(set! interactive? false)

(define pi (* (atan 1) 4))
(set! temperature 300) ; temperature
(define-param rr 10e-9)

;; (define p-cMOSFET-Si
;;   (make params-cMOSFET
;; 	(sizes (make NW-radial (radius 6.25e-9) ; nanowire diameter, m
;; 							))
;; 	(ni 1.45e16) ; intrinsic carrier concentration, m^-3	
;;     (Lg 1e-6) ; gate length, m
;; 	(eps-s 11.6) ; dielectric constant
;; 	(dphi 0) ; work function difference, eV
;; 	(tox 1.5e-9) ; gate oxide thickness, m
;; 	(eps-ox 3.9) ; dielectric constant of gate oxide
;; 	(mobility 0.04) ;mobility, m^2/Vs
;; 	))

(define p-cMOSFET-InAs (lambda (r0)
  (make params-cMOSFET
	(sizes (make NW-radial (radius r0) ; nanowire diameter, m
							))
	(ni 1e21) ; intrinsic carrier concentration, m^-3	
    (Lg 1e-6) ; gate length, m
	(eps-s 15.16) ; dielectric constant
	(dphi 0) ; work function difference, eV
	(tox 3e-9) ; gate oxide thickness, m
	(eps-ox 20) ; dielectric constant of gate oxide
	(mobility 0.1) ;mobility, m^2/Vs
	)))

;; (set! FET-params (p-cMOSFET-InAs 5e-10))

(define Ids-density
  (lambda (Vds Vgs p-cMOSFET)
	(let ((r0 (get-radius p-cMOSFET)))
	  (begin
		(set! FET-params p-cMOSFET)
		(/ (Ids-cMOSFET Vds Vgs) (* 2 pi r0))))))

(define Ids2-density
  (lambda (Vds Vgs p-cMOSFET)
	(let ((r0 (get-radius p-cMOSFET)))
	  (begin
		(set! FET-params p-cMOSFET)
		(/ (Ids2-cMOSFET Vds Vgs) (* 2 pi r0))))))

;; (print (Ids-density 0.5 -0.07817 (p-cMOSFET-InAs 10e-9)) "\t"
;; 	   (Ids-density 0.5 (+ -0.0781 0.5) (p-cMOSFET-InAs 10e-9)) "\n")

(define (func-for-find-Voff-cMOSFET  Vds Ioff p-cMOSFET)
  (lambda (Vgs)
	  (- (Ids-density Vds Vgs p-cMOSFET) Ioff)))
(define (func-for-find-Voff2-cMOSFET  Vds Ioff p-cMOSFET)
  (lambda (Vgs)
	  (- (Ids2-density Vds Vgs p-cMOSFET) Ioff)))

;; (print
;; ((func-for-find-Voff-cMOSFET 0.5 1e-4 (p-cMOSFET-InAs 10e-9)) -1) "\t"
;; ((func-for-find-Voff-cMOSFET 0.5 1e-4 (p-cMOSFET-InAs 10e-9))  1) "\n")

(define find-Voff-cMOSFET (lambda (Vds Ioff p-cMOSFET left right)
  (find-root (func-for-find-Voff-cMOSFET Vds Ioff p-cMOSFET)  1e-7 left right)
  ))
(define find-Voff2-cMOSFET (lambda (Vds Ioff p-cMOSFET left right)
  (find-root (func-for-find-Voff2-cMOSFET Vds Ioff p-cMOSFET)  1e-7 left right)
  ))

;; (define Idsmax1-radius
;;   (lambda (r0 Vdsmax Ioff Vov)
;; 	(let ((p-cMOSFET (p-cMOSFET-InAs r0)))
;; 	  (let ((Voff (find-Voff-cMOSFET Vdsmax Ioff p-cMOSFET -1 1)))
;; 		(let  ((Ioff0 (Ids-density Vdsmax Voff p-cMOSFET))
;; 			   (Ion   (Ids-density Vdsmax (+ Voff Vov) p-cMOSFET)))
;; 		  (print r0 "\t" Ioff0 "\t" Voff "\t" Ion "\n"))))))

;; (define Idsmax2-radius
;;   (lambda (r0 Vdsmax Ioff Vov)
;; 	(let ((p-cMOSFET (p-cMOSFET-InAs r0)))
;; 	  (let ((Voff (find-Voff2-cMOSFET Vdsmax Ioff p-cMOSFET -1 1)))
;; 		(let  ((Ioff0 (Ids2-density Vdsmax Voff p-cMOSFET))
;; 			   (Ion   (Ids2-density Vdsmax (+ Voff Vov) p-cMOSFET)))
;; 		  (print r0 "\t" Ioff0 "\t" Voff "\t" Ion "\n"))))))

(define radius-list (interpolate 89 '(10e-9 100e-9)))
;; (for-each (lambda (r0) (Idsmax1-radius r0 0.5 1e-4 0.5)) radius-list)
;; (for-each (lambda (r0) (Idsmax2-radius r0 0.5 1e-4 0.5)) radius-list)

;; (define Idsmax-radius
;;   (lambda (r0 Vdsmax Ioff Vov)
;; 	(let ((p-cMOSFET (p-cMOSFET-InAs r0)))
;; 	  (let ((Voff1 (find-Voff-cMOSFET Vdsmax Ioff p-cMOSFET -1 1))
;; 		(Voff2 (find-Voff2-cMOSFET Vdsmax Ioff p-cMOSFET -1 1)))
;; 		(let  ((Ioff1 (Ids-density Vdsmax Voff1 p-cMOSFET))
;; 			   (Ion1   (Ids-density Vdsmax (+ Voff1 Vov) p-cMOSFET))
;; 			   (Ioff2 (Ids2-density Vdsmax Voff2 p-cMOSFET))
;; 			   (Ion2   (Ids2-density Vdsmax (+ Voff2 Vov) p-cMOSFET)))
;; 		  (print r0 "\t" Ioff1 "\t" Voff1 "\t" Ion1 "\t"
;; 				 Ioff2 "\t" Voff2 "\t" Ion2 "\n"))))))

;; (for-each (lambda (r0) (Idsmax-radius r0 0.5 1e-4 0.5)) radius-list)

;; (define Vgs1 (interpolate 99 (list 0 0.5)))

(define Vgs1 '(0.3))
(define Ids-Vgs 
  (lambda (r0 Vds Ioff)
	(let ((p-cMOSFET (p-cMOSFET-InAs r0)))
	  (let ((Voff1 (find-Voff-cMOSFET Vds Ioff p-cMOSFET -1 1))
			(Voff2 (find-Voff2-cMOSFET Vds Ioff p-cMOSFET -1 1)))
		(for-each (lambda (Vgs)
					(print "r0:\t" Vgs "\t"
						   Voff1 "\t" (Ids-density Vds (+ Voff1 Vgs) p-cMOSFET) "\t"
						   Voff2 "\t" (Ids2-density Vds (+ Voff2 Vgs) p-cMOSFET) "\n"))
				  Vgs1)))))

;; (Ids-Vgs rr 0.5 1e-4)

(print (Ids-density 0.5 0.16 (p-cMOSFET-InAs 100e-9)) "\n")
