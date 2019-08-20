
(set! interactive? false)

(define Eg 3.4)
(define ems 0.2)
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
(define ene 0)

(define epsOX 8.9)
(define epsS 8.5)
(define tOX 10e-9)

(define Cox (Cox-rect epsOX tOX W1-NW W2-NW))
(define Cc  (Cc-rect epsS W1-NW W2-NW))

(set! FET-params
  (make params-cMOSFET
	(sizes (make NW-rect (W1 W1-NW) (W2 W1-NW)))
	(ni 1.45e16) ; intrinsic carrier concentration, m^-3	
    (Lg 1e-6) ; gate length, m
	(eps-s epsS) ; dielectric constant
	(dphi 0) ; work function difference, eV
	(tox tOX) ; gate oxide thickness, m
	(eps-ox epsOX) ; dielectric constant of gate oxide
	(mobility 0.04) ;mobility, m^2/Vs
	(alpha-nonparabolicity alpha)
	(effective-mass ems)
	))

(set! ballistic-params
  (make params-ballisticFET
	(Fermi-Energy 0)
	(alpha-D 0)
	(alpha-G 1)
	(C-eff (/ (* Cox Cc) (+ Cox Cc)))))

;; (funcval-E00 0 param)

(print (E0-rect1d 0 0) "\n")
