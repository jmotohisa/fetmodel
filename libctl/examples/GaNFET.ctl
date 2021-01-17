;; GaN cMOSFET

(set! temperature 300) ; temperature
;; Parameters for GaN @300K
(define Eg_GaN 3.2)
(define epsS_GaN 8.9)
(define ems_GaN 0.2)
;; (define alpha_GaN = fetmodel.alpha_NP(Eg_GaN, ems_GaN)
(define Nc (* 4.3e20 (sqrt (* 300 300 300))))
(define Nv (* 8.9e21 (sqrt (* 300 300 300))))
(define ni_GaN (* (sqrt (* Nc Nv)) (exp (/ Eg_GaN (* -2 0.026)))))
(print Nc "\t" Nv "\t" ni_GaN "\n")

;; common parameters
(define epsOX 20.)
(define tOX 3e-9)
(define temperature 300)
(define W1 10e-9)
(define W2 10e-9)
(define Cox (Cox-rect epsOX tOX W1 W2))
(define Cc (Cc-rect epsS_GaN W1 W2))
;; alpha_D = 0.
;; alpha_G = 1.
(print Cox "\t" Cc "\n")

(set! FET-params
  (make params-cMOSFET
	(sizes (make NW-radial (radius (/ W1 2)))) ; nanowire diameter, m
	(ni ni_GaN) ; intrinsic carrier concentration, m^-3	
	(Lg 1e-6) ; gate length, m
	(eps-s epsS_GaN) ; dielectric constant
	(dphi 0) ; work function difference, eV
	(tox tOX) ; gate oxide thickness, m
	(eps-ox epsOX) ; dielectric constant of gate oxide
	(mobility 0.1) ;mobility, m^2/Vs
	))

;; parasitic
(set! params-parasitic
  (make parasitic-component 
	(Rs 0) ; source resistance
    (Rd 0) ; drain resistance
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

;; (print "\n")
;; ;; Fig.3; Vds = 0.5 V
(for-each (lambda (Vgs) (print "Fig3:\t" Vgs "\t"
							   (Ids2-cMOSFET 0.5 Vgs) "\t"
							   (Ids-cMOSFET 0.5 Vgs) "\n")) Vgs1)
