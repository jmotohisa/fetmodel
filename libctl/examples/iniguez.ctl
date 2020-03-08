;; B. Iniguez et al., IEEE Trans. Elec. Dev. 52, No.8, Auguust 2015 (p.1868)

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

;; Fig.2
(for-each (lambda (Vgs) (print "Fig2:\t" Vgs "\t"
							   (Qcharge2-cMOSFET Vgs) "\t"
							   (Qcharge-cMOSFET Vgs) "\n")) Vgs1)

(print "\n")
;; Fig.3; Vds = 0.1 V
(for-each (lambda (Vgs) (print "Fig3:\t" Vgs "\t"
							   (Ids2-cMOSFET 0.1 Vgs) "\t"
							   (Ids-cMOSFET 0.1 Vgs) "\n")) Vgs1)

(print "\n")
;; Fig.4: Vds = 1V;
(for-each (lambda (Vgs) (print "Fig4:\t" Vgs "\t"
							   (Ids2-cMOSFET 1 Vgs) "\t"
							   (Ids-cMOSFET  1 Vgs) "\n")) Vgs1)
