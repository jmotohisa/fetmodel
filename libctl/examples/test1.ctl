
(set! interactive? false)
(define Eg 3.4)
(define ems 0.2)
(define alpha (alphaNP00 Eg ems))
(define temperature 300)

;; (defineene = np.arange(-0.1, 0.1, 0.005)
;; dens = np.empty_like(ene)

;; for i, e in enumerate(ene):
;;     dens[i] = ballistic.density1d_parabollic00(e, alpha, 0.067, temperature)

(define W1 10e-9)
(define W2 8e-9)

(define gamma_nm (gamma-nm00 alpha ems W1 W2 1 1))
(define Enm      (E-nm0 alpha gamma_nm))
(define alpha_nm (alpha-nm0 alpha gamma_nm))
(define ems_nm  (ems-nm0 ems gamma_nm))
(define ene 0)

;; (print alpha "\t" gamma_nm "\t" Enm "\t" alpha_nm "\t" ems_nm "\n")

(print
 (density1d-parabollic00 ene 0 ems temperature)
 "\t"
 (density1d-nonpara00 ene 0 0 ems temperature)
 "\n")

(print 
 (density1d-parabollic00 ene Enm ems_nm temperature)
 "\t"
 (density1d-nonpara00 ene Enm 0 ems_nm temperature)
 "\n")

(print
 (density1d-nonpara00 ene Enm alpha_nm ems_nm temperature)
 "\n")

