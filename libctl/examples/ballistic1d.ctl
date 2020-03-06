;; test of ballistic FET (with lower level interface)

;;
;; (define (pw-amp k x0) (lambda (x)
;;   (exp (* 0+1i (vector3-dot k (vector3+ x x0))))))
(define (func-for-findroot-E0-rect1dNP  Vds Vgs p-NWFET p-bFET)
  (lambda (ene0)
	(let ((EFermi (object-property-value p-bFET 'EFermi))
		   (alpha (object-property-value p-NWFET 'alpha))
		   (ems (object-property-value p-NWFET 'ems))
		   (W1 (get-radius p-NWFET))
		   (W2 (get-radius2 p-NWFET))
		   (nmax (object-proeprty-value p-NWFET 'n-max))
		   (mmax (object-proeprty-value p-NWFET 'm-max))
		   (Ceff (object-proeprty-value p-bFET 'C-eff))
		   (alpha-D (object-proeprty-value p-bFET 'alpha-D))
		   (alpha-S (object-proeprty-value p-bFET 'alpha-S)))
	 (let ((nd1_S (density1d-rec1dNP-all0 (- EFermi ene0) alpha ems temperature W1 W2 nmax mmax))
			(nd1_D (density1d-rec1dNP-all0 EFermi alpha ems temperature W1 W2 nmax mmax)))
	  (let ((q0 (* (+ nd1_S nd1_D) elementary-charge (/ (* 2 Ceff)))))
	   (+ ene0 (* alpha-D Vds) (* alpha-S Vgs) (- q0)))))))
    ;; """
    ;; Python implementation of the function to find top of the barrier
    ;; based on fetmodel.density1d_rect1dNP_all0
    ;; Nonparabolic band
    ;; """
    ;; n1d_S = fetmodel.density1d_rect1dNP_all0(
    ;;     p.EFermi - ene0, p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax)
    ;; n1d_D = fetmodel.density1d_rect1dNP_all0(
    ;;     p.EFermi - ene0 - Vds, p.alpha, p.ems, p.temp, p.W1, p.W2, p.nmax, p.mmax)
    ;; q0 = const.elementary_charge * (n1d_S + n1d_D) / (2 * p.Ceff)
    ;; return ene0 + (p.alpha_D * Vds + p.alpha_G * Vgs - q0)


;; """
;; Get top of the barrier (Python implementation)
;; Nonparabolic band
;; """
;; (find-root function tolerance arg-min arg-max)
(define E0-rect1dNP-root (lambda Vds Vgs p-NWFET p-bFET left right)
  (find-root (func-for-findroot-E0-rect1dNP  Vds Vgs p-NWFET p-bFET)  1e-7 left right)
  )
    ;; left0 = -(p.alpha_D*Vds+p.alpha_G*Vgs) - 0.2
    ;; left2=min([left,left0])
    ;; right0=-(p.alpha_D*Vds+p.alpha_G*Vgs)
    ;; right2=max([right,right0])
    ;; e0 = optimize.root_scalar(func_for_findroot_E0_rect1dNP,
    ;;                           args=(Vds, Vgs, p), x0=left2, x1=right2)
    ;; if e0.converged==True:
    ;;     return e0.root
    ;; else:
    ;;     print("EFs convergence error !")
    ;;     print(e0)
    ;;     return 0

(define Ids-ballistic1d-rect1dNP
  (lambda (Vds Vgs p-NWFET p-bFET left right) 
   ((let ((E0 (E0-rect1dNP-root Vds Vgs p-NWFET p-bFET left right))
		  
	  

    ;; e0=E0_rect1dNP_root(Vds,Vgs,p,left,right)
    ;; nlist = np.arange(1, p.nmax+1, dtype=np.int64)
    ;; mlist = np.arange(1, p.mmax+1, dtype=np.int64)
    ;; cur = 0
    ;; for n in nlist:
    ;;     for m in mlist:
    ;;         Enmp = fetmodel.Ep_nm_rect1d(p.ems, p.W1, p.W2, int(n), int(m))
    ;;         gamma_nm = fetmodel.gamma_nm_NP(Enmp, p.alpha)
    ;;         Enm = fetmodel.E_nm_NP(p.alpha, gamma_nm)
    ;;         print('parabolic, nonparabolic',Enmp,Enm)
    ;;         cur1 = func_FD0(EFs-Enm-e0, p.temp)
    ;;         cur2 = func_FD0(EFs-Enm-e0-Vds, p.temp)
    ;;         cur += cur1-cur2

    ;; return (cur*2*const.elementary_charge/const.h*p.temp*const.Boltzmann)
