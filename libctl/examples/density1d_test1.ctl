;; test of density1d

(set! interactive? false)
(define Eg 0.36)
(define ems 0.0671)
(define alpha (alpha-NP Eg ems))
(set-param! temperature 300)

;; (defineene = np.arange(-0.1, 0.1, 0.005)
;; dens = np.empty_like(ene)

;; for i, e in enumerate(ene):
;;     dens[i] = ballistic.density1d_parabollic00(e, alpha, 0.067, temperature)

(define W1-NW 10e-9)
(define W2-NW 8e-9)

(define gamma_nm (gamma-nm-rect1dNP alpha ems W1-NW W2-NW 1 1))
(define Enm      (E-nm-NP alpha gamma_nm))
(define alpha_nm (alpha-nm-NP alpha gamma_nm))
(define ems_nm  (ems-nm-NP ems gamma_nm))

(define ene (interpolate 99 '(-0.1 0.2)))

(print alpha "\t" gamma_nm "\t" Enm "\t" alpha_nm "\t" ems_nm "\n")

;; (print
;;  (density1d0 ene 0 ems temperature)
;;  "\t"
;;  (density1d-NP0 ene 0 0 ems temperature)
;;  "\n")

;; (print 
;;  (density1d0 ene Enm ems_nm temperature)
;;  "\t"
;;  (density1d-NP0 ene Enm 0 ems_nm temperature)
;;  "\n")

;; (print
;;  (density1d-NP0 ene Enm alpha_nm ems_nm temperature)
;;  "\n")
;;

(define dens 
  (map (lambda (x) (density1d0 x Enm ems_nm temperature)) ene))

(define dens2 
  (map (lambda (x) (density1d-NP0 x Enm alpha_nm ems_nm temperature)) ene))

(define dens3
  (map (lambda (x) (density1d-rect1dNP-all0 x alpha ems temperature
												 W1-NW W2-NW 2 2)) ene))
  
;; (for-each (lambda (ene0 dens00 dens20 dens30)
;; 			(print ene0 "\t" dens00 "\t" dens20 "\t" dens30 "\n"))
;; 		  ene dens dens2 dens3)

