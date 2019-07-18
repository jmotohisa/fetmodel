;; B. Iniguez et al., IEEE Trans. Elec. Dev. 52, No.8, Auguust 2015 (p.1868)

(set! radius 6.25e-9) ; nanowire diameter, m
(set! Lg 1e-6) ; gate length, m
(set! eps-semi 11.6) ; dielectric constant
(set! Rs 50) ; source resistance
(set! Rd 50) ; drain resistance
(define-param Vds 1)
(define-param Vgs 0.5) 


; input variable for MOSFET

(set! temp 300) ; temperature
(set! ni 1.45e16) ; intrinsic carrier concentration, m^-3
(set! dphi 0) ; work function difference, eV
(set! tox 1.5e-9) ; gate oxide thickness, m
(set! eps-ox 3.9) ; dielectric constant of gate oxide
(set! mue 0.04) ;mobility, m^2/Vs

(set! Cox (Cox-radial radius tox eps-ox))

;; (define Vgs1 (interpolate 99 (list 0 2)))
;; (define Vds1 (interpolate 99 (list 0 2)))

(set! interactive? false)

(print (func-Ids-cMOSFET Vds Vgs) "\n")
