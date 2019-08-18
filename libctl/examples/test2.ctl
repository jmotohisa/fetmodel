
(set! interactive? true)
(define Eg 3.4)
(define ems 0.2)
(define alpha (alphaNP00 Eg ems))

(set-param! temperature 300)

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

(define epsOX 8.9)
(define epsS 8.5)
(define tOX 10e-9)

(define Cox (Cox-rect epsOX tOX W1 W2))
(define Cc  (Cc-rect epsS W1 W2))

(define param (make params-ballisticFET0-type
				(Fermi-Energy 0)
				(C-eff (/ (* Cox Cc) (+ Cox Cc)))
				(size-W1 W1)
				(size-W2 W2)
				(effective-mass ems)
				(Vds 0)
				(Vgs 0)
				(alpha-NP alpha)
				))

(funcval-E00 0 param)
;;(E0-func param)
