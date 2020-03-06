(set! interactive? false)

(define pi (* (atan 1) 4))

(define (func-for-find-Voff-cMOSFET  Vds Ioff p-cMOSFET)
  (lambda (Vgs)
	(let ((radius (get-radius p-cMOSFET)))
		   (begin (set! FET-params p-cMOSFET)
				  (- (/ (func-Ids0-cMOSFET Vds Vgs) (* 2 pi radius)) Ioff)
				  )
		   )))

(define find-Voff-cMOSFET (lambda (Vds Ioff p-cMOSFET left right)
  (find-root (func-for-find-Voff-cMOSFET Vds Ioff p-cMOSFET)  1e-7 left right)
  ))

(set! temperature 300) ; temperature

(define p-cMOSFET-Si
  (make params-cMOSFET
	(sizes (make NW-radial (radius 6.25e-9) ; nanowire diameter, m
							))
	(ni 1.45e16) ; intrinsic carrier concentration, m^-3	
    (Lg 1e-6) ; gate length, m
	(eps-s 11.6) ; dielectric constant
	(dphi 0) ; work function difference, eV
	(tox 1.5e-9) ; gate oxide thickness, m
	(eps-ox 3.9) ; dielectric constant of gate oxide
	(mobility 0.04) ;mobility, m^2/Vs
	))

(define p-cMOSFET-InAs
  (make params-cMOSFET
	(sizes (make NW-radial (radius 10e-9) ; nanowire diameter, m
							))
	(ni 1e21) ; intrinsic carrier concentration, m^-3	
    (Lg 1e-6) ; gate length, m
	(eps-s 15.16) ; dielectric constant
	(dphi 0) ; work function difference, eV
	(tox 3e-9) ; gate oxide thickness, m
	(eps-ox 20) ; dielectric constant of gate oxide
	(mobility 0.1) ;mobility, m^2/Vs
	))

(set! FET-params p-cMOSFET-InAs)

;; (print (find-Voff-cMOSFET 0.5 1e-4 -1 1))

(define Idsmax-radius
  (lambda (r0 Vdsmax Ioff Vov)
	(let ((p-cMOSFET (make params-cMOSFET
					 (sizes (make NW-radial (radius r0) ; nanowire diameter, m
								  ))
					 (ni 1e21) ; intrinsic carrier concentration, m^-3	
					 (Lg 1e-6) ; gate length, m
					 (eps-s 15.15) ; dielectric constant
					 (dphi 0) ; work function difference, eV
					 (tox 3e-9) ; gate oxide thickness, m
					 (eps-ox 20) ; dielectric constant of gate oxide
					 (mobility 0.1) ;mobility, m^2/Vs
					 )))
	  (let ((Voff (find-Voff-cMOSFET Vdsmax Ioff p-cMOSFET -1 1)))
		(let  ((Ioff0 (/ (func-Ids0-cMOSFET Vdsmax Voff) (* 2 pi r0)))
			   (Ion (/ (func-Ids0-cMOSFET Vdsmax (+ Voff Vov)) (* 2 pi r0))))
		  (print r0 "\t" Ioff0 "\t" Voff "\t" Ion "\n"))))))

(define radius-list (interpolate 99 '(10e-9 100e-9)))
(for-each (lambda (r0) (Idsmax-radius r0 0.5 1e-4 0.5)) radius-list)



	
