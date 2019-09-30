;; test of ballisticFET1d (1)

(set! interactive? false)

(define Eg 0.36)
(define epsOX 8.5)
(define epsS 8.9)
(define tOX 20e-9)
(set-param! temperature 300)
(define ems 0.067)
(define W1-NW 10e-9)
(define W2-NW 8e-9)
(define alpha (alpha-NP Eg ems))
(define Cox (Cox-rect epsOX tOX W1-NW W2-NW))
(define Cc  (Cc-rect epsS W1-NW W2-NW))
(define alpha-D0 0)
(define alpha-G1 1)

(define gamma_nm (gamma-nm-rect1dNP alpha ems W1-NW W2-NW 1 1))
(define Enm      (E-nm-NP alpha gamma_nm))
(define alpha_nm (alpha-nm-NP alpha gamma_nm))
(define ems_nm  (ems-nm-NP ems gamma_nm))
(define ene 0)

;; (print alpha "\t" gamma_nm "\t" Enm "\t" alpha_nm "\t" ems_nm "\n")

(set! FET-params
  (make params-cMOSFET
	(sizes (make NW-rect (W1 W1-NW) (W2 W2-NW)))
;;	(ni 1.45e16) ; intrinsic carrier concentration, m^-3	
;;    (Lg 1e-6) ; gate length, m
	(eps-s epsS) ; dielectric constant
;;	(dphi 0) ; work function difference, eV
	(tox tOX) ; gate oxide thickness, m
	(eps-ox epsOX) ; dielectric constant of gate oxide
;;	(mobility 0.04) ;mobility, m^2/Vs
	(alpha-nonparabolicity alpha)
	(nonparabolic? true)
	(effective-mass ems)
	))

(set! ballistic-params
  (make params-ballisticFET
	(Fermi-Energy 0)
	(alpha-D alpha-D0)
	(alpha-G alpha-G1)
	(C-eff (/ (* Cox Cc) (+ Cox Cc)))))

;; (funcval-E00 0 param)

(define E0 (interpolate 99 '(-0.1 0.1)))
(define VGS 0)
(define VDS 0)

(for-each (lambda (e00)
			(print "func:\t" e00 "\t"
				   (func-e0-find e00 VGS VDS) "\n")) E0)
(print "\n")
;; (print (E0-rect1d 0 0) "\n")

;; (define ene (interpolate 99 '(-0.1 0.1)))
;; (for-each (lambda (ene0) ()) ene)

(define VDS 0)
(define VGS (interpolate 99 '(0 1)))
(for-each (lambda (VGS0) (print "E0:\t" VGS0 "\t"
								(E0-rect1d VDS VGS0) "\n")) VGS)

(define VDS (interpolate 99 '(0 1)))
;; (define Ids-func (lambda (VGS) (map (lambda (VDS0) (Ids-ballistic-rect1d VDS0 VGS 0) VDS))))

(for-each (lambda (VDS0)
			(print "Ids:\t" VDS0 "\t"
				   (Ids-ballistic-rect1d VDS0 -0.1 0) "\t"
				   (Ids-ballistic-rect1d VDS0 0 0) "\n")) VDS)
				   
