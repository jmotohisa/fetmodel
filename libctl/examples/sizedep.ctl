(set! temperature 300) ; temperature

(set! FET-params
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

(for-each (lambda (Vgs) (print "res1:\t" Vgs "\t"
							   (Qcharge2-cMOSFET Vgs) "\t"
							   (Qcharge-cMOSFET Vgs) "\n")) Vgs1)

(set! FET-params
  (make params-cMOSFET
	(sizes (make NW-radial (radius 100e-9) ; nanowire diameter, m
							))
	(ni 1.45e16) ; intrinsic carrier concentration, m^-3	
    (Lg 1e-6) ; gate length, m
	(eps-s 11.6) ; dielectric constant
	(dphi 0) ; work function difference, eV
	(tox 1.5e-9) ; gate oxide thickness, m
	(eps-ox 3.9) ; dielectric constant of gate oxide
	(mobility 0.04) ;mobility, m^2/Vs
	))

(for-each (lambda (Vgs) (print "res2:\t" Vgs "\t"
							   (Qcharge2-cMOSFET Vgs) "\t"
							   (Qcharge-cMOSFET Vgs) "\n")) Vgs1)
