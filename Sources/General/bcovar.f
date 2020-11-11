!> \file bcovar.f

      SUBROUTINE bcovar_par(lu, lv, tpxc, ier_flag)
      USE vmec_main, fpsi => bvco, p5 => cp5
      USE vmec_params, ONLY: ns4, signgs, pdamp, lamscale, ntmax,
     &                       bsub_bad_js1_flag, arz_bad_value_flag,
     &                       norm_term_flag
      USE realspace, ONLY: pextra1, pextra2, pextra3, pextra4,
     &                     pguu, pguv, pgvv, pru, pzu,
     &                     pr1, prv, pzv, pshalf, pwint, pz1,
     &                     pru0, pzu0, psqrts
      USE vforces, r12 => parmn_o, ru12 => pazmn_e, gsqrt => pazmn_o,
     &             rs => pbzmn_e, zs => pbrmn_e, zu12 => parmn_e,
     &             bsubu_e => pclmn_e, bsubv_e => pblmn_e,
     &             bsubu_o => pclmn_o, bsubv_o => pblmn_o,
     &             bsq => pbzmn_o, phipog => pbrmn_o
      USE xstuff, ONLY: pxc
      USE precon2d, ONLY: ictrl_prec2d, lHess_exact,
     &                    ctor_prec2d
      USE fbal
      USE vmec_input, ONLY: nzeta
      USE vmec_dim, ONLY: ntheta3
      USE parallel_include_module
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(nznt,ns,0:1), INTENT(INOUT) :: lu, lv
      REAL(dp), DIMENSION((1+ntor)*(1+mpol1),1:ns,1:2*ntmax),
     &   INTENT(IN) :: tpxc
      INTEGER, INTENT(inout) :: ier_flag
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!     GENERALLY, IF TEMPORAL CONVERGENCE IS POOR, TRY TO INCREASE PDAMP (< 1)
!     (STORED IN VMEC_PARAMS)
      REAL(dp), PARAMETER :: c1p5 = (one + p5)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: l, js, ndim
      REAL(dp) :: r2, volume, curpol_temp

      REAL(dp) :: arnorm, aznorm, tcon_mul

      REAL(dp) :: bcastton, bcasttoff
      REAL(dp), POINTER, DIMENSION(:,:) :: luu, luv, lvv, tau
      REAL(dp), DIMENSION(:,:), POINTER :: bsupu, bsubuh,
     &                                     bsupv, bsubvh, r12sq
      LOGICAL :: lctor
      INTEGER :: i, j, k, nsmin, nsmax, istat
      REAL(dp) :: wblocal(ns), wbtotal
      REAL(dp) :: wplocal(ns), wptotal
      REAL(dp) :: vptotal
      REAL(dp) :: fnlocal(ns), fntotal
      REAL(dp) :: fn1local(ns), fn1total
      REAL(dp) :: fnLlocal(ns), fnLtotal
!-----------------------------------------------
      IF (irst.EQ.2 .AND. iequi.EQ.0) RETURN

!
!     POINTER ALIAS ASSIGNMENTS

      tau => pextra1(:,:,1)
      luu => pextra2(:,:,1)
      luv => pextra3(:,:,1)
      lvv => pextra4(:,:,1)

      bsupu => luu
      bsubuh => bsubu_o
      bsupv => luv
      bsubvh => bsubv_o
      r12sq => bsq

!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING (MORE THAN ONE) POINTER!
!
      nsmin=t1lglob; nsmax=t1rglob
      pguu(:,nsmin:nsmax) = 0
      pguv(:,nsmin:nsmax) = 0
      pgvv(:,nsmin:nsmax) = 0

!
!     COMPUTE METRIC ELEMENTS GIJ ON HALF MESH
!     FIRST, GIJ = EVEN PART (ON FULL MESH), LIJ = ODD PART (ON FULL MESH)
!     THEN, GIJ(HALF) = < GIJ(even)> + SHALF < GIJ(odd) >

      DO l = nsmin, nsmax
         r12sq(:,l) = psqrts(:,l)*psqrts(:,l)
         pguu(:,l) = pru(:,l,0)*pru(:,l,0) + pzu(:,l,0)*pzu(:,l,0)
     &             + r12sq(:,l)*(pru(:,l,1)*pru(:,l,1)
     &             + pzu(:,l,1)*pzu(:,l,1))
         luu(:,l) = (pru(:,l,0)*pru(:,l,1) +  pzu(:,l,0)*pzu(:,l,1))*2
         phipog(:,l) = 2*pr1(:,l,0)*pr1(:,l,1)
      END DO

      IF (lthreed) THEN
         DO l = nsmin, nsmax
            pguv(:,l) = pru(:,l,0) * prv(:,l,0) + pzu(:,l,0)*pzv(:,l,0)
     &                + r12sq(:,l) * (pru(:,l,1)*prv(:,l,1)
     &                + pzu(:,l,1)*pzv(:,l,1))
            luv(:,l) = pru(:,l,0) * prv(:,l,1) + pru(:,l,1)*prv(:,l,0)
     &               + pzu(:,l,0)*pzv(:,l,1) + pzu(:,l,1)*pzv(:,l,0)
            pgvv(:,l) = prv(:,l,0) * prv(:,l,0) + pzv(:,l,0)*pzv(:,l,0)
     &                + r12sq(:,l) * (prv(:,l,1)*prv(:,l,1)
     &                + pzv(:,l,1)*pzv(:,l,1) )
            lvv(:,l) = (prv(:,l,0) * prv(:,l,1) +
     &                  pzv(:,l,0)*pzv(:,l,1))*2
         END DO
      END IF

      r12sq(:,nsmin:nsmax) = pr1(:,nsmin:nsmax,0)*pr1(:,nsmin:nsmax,0)
     &                     + r12sq(:,nsmin:nsmax)*pr1(:,nsmin:nsmax,1)*
     &                       pr1(:,nsmin:nsmax,1)

      DO l = t1rglob, MAX(t1lglob,2), -1
         pguu(:,l) = p5*(pguu(:,l) + pguu(:,l-1) +
     &             + pshalf(:,l)*(luu(:,l) + luu(:,l-1)))
         r12sq(:,l) = p5*(r12sq(:,l)+r12sq(:,l-1)+pshalf(:,l)*  !Comment: r12sq = r12**2
     &                    (phipog(:,l) + phipog(:,l-1)))
      END DO

      IF (lthreed) THEN
         DO l = t1rglob, MAX(t1lglob,2), -1
            pguv(:,l) = p5*(pguv(:,l) + pguv(:,l-1) +
     &                      pshalf(:,l)*(luv(:,l) + luv(:,l-1)))
            pgvv(:,l) = p5*(pgvv(:,l) + pgvv(:,l-1) +
     &                      pshalf(:,l)*(lvv(:,l) + lvv(:,l-1)))
         END DO
      END IF

      pguv(:,1)=0
      pgvv(:,1)=0

      nsmin = tlglob; nsmax = t1rglob
      DO l = nsmin, nsmax
         tau(:,l) = gsqrt(:,l)
         gsqrt(:,l) = r12(:,l)*tau(:,l)
      END DO
      gsqrt(:,1) = gsqrt(:,2)

      nsmin = MAX(2,tlglob); nsmax = t1rglob
      pgvv(:,nsmin:nsmax) = pgvv(:,nsmin:nsmax)
     &                    + r12sq(:,nsmin:nsmax)
      pgvv(:,1) = 0

!CATCH THIS AFTER WHERE LINE BELOW phipog = 0
      nsmin = MAX(2,tlglob); nsmax = t1rglob
      WHERE (gsqrt(:,nsmin:nsmax) .ne. zero)
         phipog(:,nsmin:nsmax) = one/gsqrt(:,nsmin:nsmax)
      END WHERE
      phipog(:,1) = 0

      vp(1) = 0
      vp(ns+1) = 0
      DO js = nsmin, nsmax
         vp(js) = signgs*SUM(gsqrt(:,js)*pwint(:,js))
      END DO

!
!     COMPUTE CONTRA-VARIANT COMPONENTS OF B (Bsupu,v) ON RADIAL HALF-MESH
!     TO ACCOMODATE LRFP=T CASES, THE OVERALL PHIP FACTOR (PRIOR TO v8.46)
!     HAS BEEN REMOVED FROM PHIPOG, SO NOW PHIPOG == 1/GSQRT!
!
!     NOTE: LU = LAMU == d(LAM)/du, LV = -LAMV == -d(LAM)/dv COMING INTO THIS ROUTINE
!     WILL ADD CHIP IN CALL TO ADD_FLUXES. THE NET BSUPU, BSUPV ARE (PHIPOG=1/GSQRT AS NOTED ABOVE):
!
!          BSUPU = PHIPOG*(chip + LAMV*LAMSCALE),
!          BSUPV = PHIPOG*(phip + LAMU*LAMSCALE)
!

      nsmin = t1lglob
      nsmax = t1rglob

      DO l = nsmin, nsmax
         lu(:,l,:) = lu(:,l,:)*lamscale
         lv(:,l,:) = lv(:,l,:)*lamscale
         lu(:,l,0) = lu(:,l,0) + phipf(l)
      END DO

      nsmin = MAX(2,t1lglob)
      nsmax = t1rglob
      DO l = nsmin, nsmax
         bsupu(:,l) = p5*phipog(:,l) * (lv(:,l,0) + lv(:,l-1,0)
     &              + pshalf(:,l)*(lv(:,l,1) + lv(:,l-1,1)))
         bsupv(:,l) = p5*phipog(:,l) * (lu(:,l,0) + lu(:,l-1,0)
     &              + pshalf(:,l)*(lu(:,l,1) + lu(:,l-1,1)))
      END DO

!v8.49: add ndim points
      IF (rank .EQ. 0) THEN
         bsupu(:,1) =0; ! bsupu(ndim) = 0
         bsupv(:,1) =0; ! bsupv(ndim) = 0
      END IF

!
!     UPDATE IOTA EITHER OF TWO WAYS:
!     1)  FOR ictrl_prec2d = 0, SOLVE THE LINEAR ALGEBRAIC EQUATION <Bsubu> = icurv
!         FOR iotas (after testing, this is preferred way)
!     2)  FOR ictrl_prec2d > 0, EVOLVE IOTAS IN TIME, USING Force-iota  = <Bsubu> - icurv.
!         IOTAS IS "STORED" AT LOCATION LAMBDA-SC(0,0) IN XC-ARRAY
!

!     COMPUTE (IF NEEDED) AND ADD CHIP TO BSUPU
      CALL add_fluxes_par(phipog, bsupu, bsupv, .TRUE.)

!
!     COMPUTE LAMBDA FORCE KERNELS (COVARIANT B COMPONENT bsubu,v) ON RADIAL HALF-MESH
!
      nsmin = t1lglob
      nsmax = t1rglob
      DO l = nsmin, nsmax
         bsubuh(:,l) = pguu(:,l)*bsupu(:,l) + pguv(:,l)*bsupv(:,l)
         bsubvh(:,l) = pguv(:,l)*bsupu(:,l) + pgvv(:,l)*bsupv(:,l)
      END DO

!v8.49
!
!     COMPUTE MAGNETIC AND KINETIC PRESSURE ON RADIAL HALF-MESH
!
      nsmin = t1lglob
      nsmax = t1rglob
      DO l = nsmin, nsmax
         bsq(:,l) = p5*(bsupu(:,l)*bsubuh(:,l) + bsupv(:,l)*bsubvh(:,l))
      END DO

      nsmin = MAX(2,tlglob)
      nsmax = MIN(ns,t1rglob)
      pres(nsmin:nsmax) = mass(nsmin:nsmax)/vp(nsmin:nsmax)**adiabatic
      pres(1)=0

      IF (ictrl_prec2d .LE. 1) THEN
         DO l = tlglob, trglob
            wblocal(l) = SUM(pwint(:,l)*gsqrt(:,l) * bsq(:,l))
            wplocal(l) = vp(l)*pres(l)
         END DO

         CALL Gather1XArray(wblocal)
         wbtotal = SUM(wblocal(2:ns))
         CALL Gather1XArray(wplocal)
         wptotal = SUM(wplocal(2:ns))
         wb = hs*ABS(wbtotal)
         wp = hs*wptotal
      END IF

!     ADD KINETIC PRESSURE TO MAGNETIC PRESSURE
      nsmin = tlglob
      nsmax = t1rglob
      DO l=nsmin, nsmax
         bsq(:,l) = bsq(:,l) + pres(l)
         lvv(:,l) = phipog(:,l)*pgvv(:,l)
      END DO

!SPH122407-MOVED HERE: COMPUTE LAMBDA FULL MESH FORCES
!     NOTE: bsubu_e is used here ONLY as a temporary array

      nsmin = tlglob
      nsmax = MIN(ns - 1, trglob)
      bsubv_e(:,nsmin:nsmax) = p5*(lvv(:,nsmin:nsmax)+
     &                             lvv(:,nsmin+1:nsmax+1))
     &                       * lu(:,nsmin:nsmax,0)
      bsubv_e(:,ns) = p5*lvv(:,ns)*lu(:,ns,0)

      nsmin = tlglob
      nsmax = t1rglob
      DO l = nsmin, nsmax
         lvv(:,l) = lvv(:,l)*pshalf(:,l)
         bsubu_e(:,l) = pguv(:,l)*bsupu(:,l) !*sigma_an(:nrzt) !sigma_an=1 isotropic
      END DO

      nsmin = tlglob
      nsmax = MIN(ns - 1, trglob)
      DO l = nsmin, nsmax
         bsubv_e(:,l) = bsubv_e(:,l)
     &                + p5*((lvv(:,l)+lvv(:,l+1))*lu(:,l,1) +
     &                      bsubu_e(:,l) + bsubu_e(:,l+1))
      END DO
      bsubv_e(:,ns) = bsubv_e(:,ns)
     &              + p5*(lvv(:,ns)*lu(:,ns,1) + bsubu_e(:,ns))

!
!     COMPUTE AVERAGE FORCE BALANCE AND TOROIDAL/POLOIDAL CURRENTS
!
!WAC: UPDATE buco, bvco AFTER pressure called (Gather buco, bvco in calc_fbal_par
      CALL calc_fbal_par(bsubuh, bsubvh)
      rbtor0= c1p5*fpsi(2)  - p5*fpsi(3)
      rbtor = c1p5*fpsi(ns) - p5*fpsi(ns-1)
!
!     (SPH:08/19/04)
!     MUST AVOID BREAKING TRI-DIAGONAL RADIAL COUPLING AT EDGE WHEN USING PRECONDITIONER
!     CTOR IS PASSED TO VACUUM TO COMPUTE EDGE BSQVAC, SO IT CAN ONLY DEPEND ON NS, NS-1
!     THUS, CTOR ~ buco(ns) WORKS, WITH REMAINDER A FIXED CONSTANT.
!
!     ALSO, IF USING FAST SWEEP IN COMPUTE_BLOCKS, MUST MAKE CTOR CONSTANT
!     TO AVOID BREAKING SYMMETRY OF A+(ns-1) AND B-(ns) HESSIAN ELEMENTS
!
!     TO GET CORRECT HESSIAN, USE THE CTOR=ctor_prec2d +... ASSIGNMENT
!     FOR ictrl_prec2d.ne.0 (replace ictrl_prec2d.gt.1 with ictrl_prec2d.ne.0 in IF test below)
!
!

!     NEXT COMPUTE COVARIANT BSUBV COMPONENT ~ lvv ON FULL RADIAL MESH BY AVERAGING HALF-MESH METRICS
!     NOTE: EDGE VALUES AT JS=NS DOWN BY 1/2
!     THIS IS NEEDED FOR NUMERICAL STABILITY

      IF (lHess_exact) THEN
         lctor = lfreeb .AND. ictrl_prec2d.NE.0      !Yields correct hessian near edge
      ELSE
         lctor = lfreeb .AND. ictrl_prec2d.GT.1      !Yields better accuracy in solution
      END IF

      IF (lctor) THEN
         IF (ictrl_prec2d .EQ. 2) THEN
            ctor_prec2d = p5*(buco(ns) - buco(ns1))
         END IF
         ctor = signgs*twopi*(buco(ns)+ctor_prec2d)
      ELSE
         ctor = signgs*twopi*(c1p5*buco(ns) - p5*buco(ns1))
      END IF

!
!     AVERAGE LAMBDA FORCES ONTO FULL RADIAL MESH
!     USE BLENDING FOR bsubv_e FOR NUMERICAL STABILITY NEAR AXIS
!
      nsmin = tlglob
      nsmax = t1rglob
      DO l = nsmin, nsmax
         lvv(:,l) = bdamp(l)
      END DO

      IF (rank.EQ.0) THEN
         IF (ANY(bsubvh(:,1) .ne. zero)) ier_flag = bsub_bad_js1_flag
         IF (ANY(bsubuh(:,1) .ne. zero)) ier_flag = bsub_bad_js1_flag
      END IF

      nsmin = tlglob
      nsmax = MIN(trglob,ns - 1)
      bsubu_e(:,nsmin:nsmax) = p5*(bsubuh(:,nsmin:nsmax) +
     &                             bsubuh(:,nsmin+1:nsmax+1))
      IF (trglob .EQ. ns) bsubu_e(:,ns) = p5*bsubuh(:,ns)

      nsmin = tlglob
      nsmax = MIN(ns - 1,trglob)
      bsubv_e(:,nsmin:nsmax) =
     &   bsubv_e(:,nsmin:nsmax)*lvv(:,nsmin:nsmax) +
     &   p5*(1 - lvv(:,nsmin:nsmax))*(bsubvh(:,nsmin:nsmax) +
     &                                bsubvh(:,nsmin+1:nsmax+1))
      IF (trglob .EQ. ns) THEN
         bsubv_e(:,ns) = bsubv_e(:,ns)*lvv(:,ns)
     &                 + p5*(1-lvv(:,ns))*bsubvh(:,ns)
      END IF

!
!     COMPUTE R,Z AND LAMBDA PRE-CONDITIONING MATRIX
!     ELEMENTS AND FORCE NORMS: NOTE THAT lu=>czmn, lv=>crmn externally
!     SO THIS STORES bsupv in czmn_e, bsupu in crmn_e
!
      nsmin = tlglob
      nsmax = t1rglob
      IF (iequi .EQ. 1) THEN
         lu(:,nsmin:nsmax,0) = bsupv(:,nsmin:nsmax)
         lv(:,nsmin:nsmax,0) = bsupu(:,nsmin:nsmax)
      END IF

!
!     COMPUTE PRECONDITIONING (1D) AND SCALING PARAMETERS
!     NO NEED TO RECOMPUTE WHEN 2D-PRECONDITIONER ON
!

      IF ((MOD(iter2-iter1,ns4).EQ.0 .AND. iequi.EQ.0) .AND.
     &    ictrl_prec2d.EQ.0) THEN
         nsmin = tlglob
         nsmax = t1rglob
         phipog(:,nsmin:nsmax) = phipog(:,nsmin:nsmax)
     &                         * pwint (:,nsmin:nsmax)

         CALL lamcal_par(phipog, pguu, pguv, pgvv)

         CALL precondn_par(bsupv, bsq, gsqrt, r12, zs, zu12,
     &                     pzu(:,:,0), pzu(:,:,1), pz1(:,:,1), arm,
     &                     ard, brm, brd, crd, rzu_fac, cos01)

         CALL precondn_par(bsupv, bsq, gsqrt, r12, rs, ru12,
     &                     pru(:,:,0), pru(:,:,1), pr1(:,:,1), azm,
     &                     azd, bzm, bzd, crd, rru_fac, sin01)

         nsmin = MAX(2,tlglob)
         nsmax = MIN(trglob,ns - 1)
         rzu_fac(nsmin:nsmax) = psqrts(1,nsmin:nsmax)
     &                        * rzu_fac(nsmin:nsmax)
         rru_fac(nsmin:nsmax) = psqrts(1,nsmin:nsmax)
     &                        * rru_fac(nsmin:nsmax)
         frcc_fac(nsmin:nsmax) = one/rzu_fac(nsmin:nsmax)
         rzu_fac(nsmin:nsmax) = rzu_fac(nsmin:nsmax)/2
         fzsc_fac(nsmin:nsmax) =-one/rru_fac(nsmin:nsmax)
         rru_fac(nsmin:nsmax) = rru_fac(nsmin:nsmax)/2

         nsmin = tlglob
         nsmax = t1rglob
         pguu(:,nsmin:nsmax) = pguu(:,nsmin:nsmax)
     &                       * r12(:,nsmin:nsmax)**2

         DO l = MAX(2,tlglob), trglob
            fnlocal(l)  = SUM(pguu(:,l)*pwint(:,l))
            fn1local(l) = SUM(tpxc(2:,l,1:ntmax)**2)
     &                  + SUM(tpxc(1:,l,ntmax+1:2*ntmax)**2)
            fnLlocal(l) = SUM((bsubuh(:,l)**2 + bsubvh(:,l)**2)
     &                  * pwint(:,l))*lamscale**2
         END DO

         CALL Gather1XArray(vp);       vptotal = SUM(vp(2:ns))
         CALL Gather1XArray(fnlocal);  fntotal = SUM(fnlocal(2:ns))
         CALL Gather1XArray(fn1local); fn1total= SUM(fn1local(2:ns))
         CALL Gather1XArray(fnLlocal); fnLtotal= SUM(fnLlocal(2:ns))

         volume = hs*vptotal
         r2 = MAX(wb,wp)/volume
         fnorm = one/(fntotal*(r2*r2))
         fnorm1=one/fn1total
         fnormL = one/fnLtotal
!
!        COMPUTE CONSTRAINT FORCE SCALING FACTOR (TCON)
!        OVERRIDE USER INPUT VALUE HERE
!
         r2 = ns
         tcon0 = MIN(ABS(tcon0), one)                              !!ignore large tcon0 from old-style files
         tcon_mul = tcon0*(1 + r2*(one/60 + r2/(200*120)))

         tcon_mul = tcon_mul/((4*r0scale**2)**2)                   !!Scaling of ard, azd (2*r0scale**2);
                                                                   !!Scaling of cos**2 in alias (4*r0scale**2)
         tcon = tcon0
         DO js = MAX(2,tlglob), MIN(ns-1,trglob)
            arnorm = SUM(pwint(:,js)*pru0(:,js)**2)
            aznorm = SUM(pwint(:,js)*pzu0(:,js)**2)
!            IF (arnorm .eq. zero .or. aznorm .eq. zero) THEN
!               STOP 'arnorm or aznorm=0 in bcovar'
!            END IF
            IF (arnorm .eq. zero .or. aznorm .eq. zero) THEN !SAL 070719
               ier_flag = arz_bad_value_flag
            END IF

            tcon(js) = MIN(ABS(ard(js,1)/arnorm),
     1                     ABS(azd(js,1)/aznorm))*tcon_mul*(32*hs)**2
         END DO
         tcon(ns) = p5*tcon(ns - 1)
         IF (lasym) THEN
            tcon = p5*tcon
         END IF
      ENDIF

      CALL MPI_ALLREDUCE(MPI_IN_PLACE,ier_flag,1,MPI_INTEGER,
     1                MPI_MAX,NS_COMM,MPI_ERR)
      IF (ier_flag .ne. norm_term_flag) RETURN
!
!     COMPUTE COVARIANT BSUBU,V (EVEN, ODD) ON HALF RADIAL MESH
!     FOR FORCE BALANCE AND RETURN (IEQUI=1)
!
!
!     COMPUTE COVARIANT BSUBU,V (EVEN, ODD) ON HALF RADIAL MESH
!     FOR FORCE BALANCE AND RETURN (IEQUI=1)
!

      IF (iequi .EQ. 1) THEN
         nsmin = MAX(tlglob,2)
         nsmax = MIN(trglob,ns - 1)
         DO js = nsmax, nsmin, -1
            bsubvh(:,js) = 2*bsubv_e(:,js) - bsubvh(:,js+1)
         END DO


!     ADJUST <bsubvh> AFTER MESH-BLENDING
         nsmin=MAX(tlglob,2); nsmax=MIN(trglob,ns)
         DO js = nsmin, nsmax
            curpol_temp = fpsi(js) - SUM(bsubvh(:,js)*pwint(:,js))
            bsubvh(:,js) = bsubvh(:,js) + curpol_temp
         END DO

         bsubu_e(:,nsmin:nsmax) = bsubuh(:,nsmin:nsmax)
         bsubv_e(:,nsmin:nsmax) = bsubvh(:,nsmin:nsmax)

         bsubu_o(:,nsmin:nsmax) = pshalf(:,nsmin:nsmax)
     &                          * bsubu_e(:,nsmin:nsmax)
         bsubv_o(:,nsmin:nsmax) = pshalf(:,nsmin:nsmax)
     &                          * bsubv_e(:,nsmin:nsmax)
         RETURN
      END IF

!     MINUS SIGN => HESSIAN DIAGONALS ARE POSITIVE
      nsmin = MAX(tlglob,2)
      nsmax = trglob
      bsubu_e(:,nsmin:nsmax) = -lamscale*bsubu_e(:,nsmin:nsmax)
      bsubv_e(:,nsmin:nsmax) = -lamscale*bsubv_e(:,nsmin:nsmax)
      bsubu_o(:,nsmin:nsmax) = psqrts(:,nsmin:nsmax)
     &                       * bsubu_e(:,nsmin:nsmax)
      bsubv_o(:,nsmin:nsmax) = psqrts(:,nsmin:nsmax)
     &                       * bsubv_e(:,nsmin:nsmax)

!
!     STORE LU * LV COMBINATIONS USED IN FORCES
!

      nsmin = MAX(tlglob,2)
      nsmax = t1rglob
      DO l = nsmin, nsmax
         lvv(:,l) = gsqrt(:,l) !*sigma_an(:,l)
         pguu(:,l) = bsupu(:,l)*bsupu(:,l)*lvv(:,l)
         pguv(:,l) = bsupu(:,l)*bsupv(:,l)*lvv(:,l)
         pgvv(:,l) = bsupv(:,l)*bsupv(:,l)*lvv(:,l)
         lv(:,l,0) = bsq(:,l)*tau(:,l)
         lu(:,l,0) = bsq(:,l)*r12(:,l)
      END DO

      END SUBROUTINE bcovar_par
