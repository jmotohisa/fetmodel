(set! temperature 300) ; temperature
(define pi (* 4 (atan 1)))

(define-param r1 10e-9)
(define-param r2 100e-9)

(set! FET-params
  (make params-cMOSFET
	(sizes (make NW-radial (radius r1) ; nanowire diameter, m
							))
	(ni 1e21) ; intrinsic carrier concentration, m^-3	
    (Lg 1e-6) ; gate length, m
	(eps-s 15.15) ; dielectric constant
	(dphi 0) ; work function difference, eV
	(tox 3e-9) ; gate oxide thickness, m
	(eps-ox 20) ; dielectric constant of gate oxide
	(mobility 0.1) ;mobility, m^2/Vs
	))

;; parasitic
(set! params-parasitic
  (make parasitic-component 
	(Rs 50) ; source resistance
    (Rd 50) ; drain resistance
	))

; input variable for MOSFET

;;(set! Cox (Cox-radial epx-ox tox radius))

(define Vgs1 (interpolate 99 (list 0 2)))
(define Vds1 (interpolate 99 (list 0 2)))

(set! interactive? false)

(for-each (lambda (Vgs)
			(print "res1:\t" Vgs "\t"
				   (/ (Ids2-cMOSFET 0.5 Vgs) (* 2 pi (get-radius FET-params)))"\t"
				   (/ (Ids-cMOSFET 0.5 Vgs)  (* 2 pi (get-radius FET-params))) "\n")) Vgs1)

(set! FET-params
  (make params-cMOSFET
	(sizes (make NW-radial (radius r2) ; nanowire diameter, m
							))
	(ni 1e21) ; intrinsic carrier concentration, m^-3	
    (Lg 1e-6) ; gate length, m
	(eps-s 15.15) ; dielectric constant
	(dphi 0) ; work function difference, eV
	(tox 3e-9) ; gate oxide thickness, m
	(eps-ox 20) ; dielectric constant of gate oxide
	(mobility 0.1) ;mobility, m^2/Vs
	))

(for-each (lambda (Vgs)
			(print "res2:\t" Vgs "\t"
				   (/ (Ids2-cMOSFET 0.5 Vgs) (* 2 pi (get-radius FET-params)))"\t"
				   (/ (Ids-cMOSFET 0.5 Vgs)  (* 2 pi (get-radius FET-params))) "\n")) Vgs1)
