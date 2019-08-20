
(set! interactive? false)
(define Eg 3.4)
(define ems 0.2)
(define alpha (alpha-NP Eg ems))

(set-param! temperature 300)

;; (defineene = np.arange(-0.1, 0.1, 0.005)
;; dens = np.empty_like(ene)

;; for i, e in enumerate(ene):
;;     dens[i] = ballistic.density1d_parabollic00(e, alpha, 0.067, temperature)

(define W1 10e-9)
(define W2 8e-9)

(define gamma_nm (gamma-nm-rect1dNP alpha ems W1 W2 1 1))
(define Enm      (E-nm-NP alpha gamma_nm))
(define alpha_nm (alpha-nm-NP alpha gamma_nm))
(define ems_nm  (ems-nm-NP ems gamma_nm))
(define ene 0)

;; (print alpha "\t" gamma_nm "\t" Enm "\t" alpha_nm "\t" ems_nm "\n")

(print
 (density1d0 ene 0 ems temperature)
 "\t"
 (density1d-NP0 ene 0 0 ems temperature)
 "\n")

(print 
 (density1d0 ene Enm ems_nm temperature)
 "\t"
 (density1d-NP0 ene Enm 0 ems_nm temperature)
 "\n")

(print
 (density1d-NP0 ene Enm alpha_nm ems_nm temperature)
 "\n")
