!> \file funct3d.f

      SUBROUTINE funct3d_par (lscreen, ier_flag)
      USE vmec_main
      USE vacmod, ONLY: bsqvac, bsqvac0, raxis_nestor, zaxis_nestor, 
     &                  nuv, nuv3
      USE vmec_params, ONLY: ntmax, norm_term_flag
      USE realspace
      USE vforces
      USE xstuff
      USE timer_sub
      USE precon2d, ONLY: ictrl_prec2d, lHess_exact, l_edge
      USE vparams, ONLY: twopi
      USE totzsp_mod
      USE tomnsp_mod
      USE timer_sub
      USE parallel_include_module
      USE parallel_vmec_module, ONLY: SAXLASTNTYPE, ZEROLASTNTYPE,
     &                                SAXPBYLASTNTYPE
      USE blocktridiagonalsolver, ONLY: L_COLSCALE

      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(inout) :: ier_flag
      LOGICAL, INTENT(in) :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l0pi, l, ivacskip
      INTEGER :: nvskip0 = 0
      REAL(dp), DIMENSION(mnmax) ::
     &   rmnc, zmns, lmns, rmns, zmnc, lmnc
      REAL(dp), DIMENSION(:,:,:), POINTER :: lu, lv
      REAL(dp) :: presf_ns, delr_mse, delt0
      REAL(dp) :: tbroadon, tbroadoff
      REAL(dp), EXTERNAL :: pmass
      INTEGER :: i, j, k, nsmin, nsmax, m
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: bcastbuf
      INTEGER, DIMENSION(4) :: bbuf
C-----------------------------------------------
      CALL second0 (tfunon)
!
!     POINTER ALIASES
!

      nfunct3d = nfunct3d + 1
      lu => pczmn;  lv => pcrmn


!     CONVERT ODD M TO 1/SQRT(S) INTERNAL REPRESENTATION
      ACTIVE1: IF (lactive) THEN 
         IF (ictrl_prec2d .EQ. 3) THEN
            CALL SAXPBYLASTNTYPE(one, pxc, one, pxcdot, pgc)
            CALL SAXLASTNTYPE(pgc, pscalxc, pgc)
         ELSE IF (ictrl_prec2d.EQ.1 .AND. l_colscale) THEN
            pgc = (pxc-pxsave)*pcol_scale + pxsave
            CALL SAXLASTNTYPE(pgc, pscalxc, pgc)
         ELSE
            CALL SAXLASTNTYPE(pxc, pscalxc, pgc)
         END IF

!     RIGID BODY SHIFT OF RMNCC(JS.GT.1,0,0) BY DELR_MSE= R00-RAXMSE
!

!     INVERSE FOURIER TRANSFORM TO S,THETA,ZETA SPACE
!     R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
!     FIRST, DO SYMMETRIC [ F(u,v) = F(-u,-v) ] PIECES
!     ON THE RANGE u = 0,pi  and v = 0,2*pi
!

         CALL totzsps_par (pgc, pr1, pru, prv, pz1, pzu, pzv, lu, lv,
     &                     prcon, pzcon, ier_flag)

!
!     ANTI-SYMMETRIC CONTRIBUTIONS TO INVERSE TRANSFORMS
!
         IF (lasym) THEN
            CALL totzspa_par (pgc, parmn, pbrmn, pextra3, pazmn, pbzmn,
     &                        pextra4, pblmn, pclmn, pextra1, pextra2)

!        SUM SYMMETRIC, ANTISYMMETRIC PIECES APPROPRIATELY
!        TO GET R, Z, L, (AND RCON, ZCON) ON FULL RANGE OF u (0 to 2*pi)
            CALL symrzl_par (pr1, pru, prv, pz1, pzu, pzv, lu, lv,
     &                       prcon, pzcon, parmn, pbrmn, pextra3, pazmn,
     &                       pbzmn, pextra4, pblmn, pclmn, pextra1,
     &                       pextra2)
         END IF

!      l0pi = ns*(1 + nzeta*(ntheta2 - 1))        !u = pi, v = 0, js = ns
!      router = r1(ns,0) + r1(ns,1)
!      rinner = r1(l0pi,0) + r1(l0pi,1)
         r00 = pr1(1,1,0)
         z00 = pz1(1,1,0)


!
!     COMPUTE CONSTRAINT RCON, ZCON
!
         nsmin=tlglob; nsmax=trglob
         DO l = nsmin, nsmax
            prcon(:,l,0) = prcon(:,l,0) + prcon(:,l,1)*psqrts(:,l)
            pzcon(:,l,0) = pzcon(:,l,0) + pzcon(:,l,1)*psqrts(:,l)
            pru0(:,l) = pru(:,l,0) + pru(:,l,1)*psqrts(:,l)
            pzu0(:,l) = pzu(:,l,0) +  pzu(:,l,1)*psqrts(:,l)
         END DO

!     COMPUTE RCON0, ZCON0 FOR FIXED BOUNDARY BY SCALING EDGE VALUES
!     SCALE BY POWER OF SQRTS, RATHER THAN USE rcon0 = rcon, etc. THIS
!     PREVENTS A DISCONTINUITY WHEN RESTARTING FIXED BOUNDARY WITH NEW RCON0....
!
!     NOTE: IN ORDER TO MAKE INITIAL CONSTRAINT FORCES SAME FOR FREE/FIXED
!     BOUNDARY, WE SET RCON0,ZCON0 THE SAME INITIALLY, BUT TURN THEM OFF
!     SLOWLY IN FREE-BOUNDARY VACUUM LOOP (BELOW)
!
         IF (ictrl_prec2d .EQ. 2) THEN
            DO l = nsmin, nsmax
               prcon0(:,l) = prcon(:,l,0)
               pzcon0(:,l) = pzcon(:,l,0)
            END DO
         ELSE IF (iter2        .EQ. iter1 .AND.
     &            ivac         .LE. 0     .AND.
     &            ictrl_prec2d .EQ. 0) THEN
#if defined(MPI_OPT)
            ALLOCATE(bcastbuf(2*nznt))
            bcastbuf(1:nznt)=prcon(:,ns,0)
            bcastbuf(nznt+1:2*nznt)=pzcon(:,ns,0)
            CALL second0(tbroadon)
            CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,nranks-1,
     &                     NS_COMM,MPI_ERR)
            CALL second0(tbroadoff)
            broadcast_time = broadcast_time + (tbroadoff-tbroadon)
            prcon(:,ns,0)=bcastbuf(1:nznt)
            pzcon(:,ns,0)=bcastbuf(nznt+1:2*nznt)
            DEALLOCATE(bcastbuf)
#endif
            DO l = nsmin, nsmax
               prcon0(:,l) = prcon(:,ns,0)*psqrts(:,l)**2
               pzcon0(:,l) = pzcon(:,ns,0)*psqrts(:,l)**2
            END DO
         END IF

!
!     COMPUTE S AND THETA DERIVATIVE OF R AND Z AND JACOBIAN ON HALF-GRID
!
         CALL jacobian_par

!     COMPUTE COVARIANT COMPONENTS OF B, MAGNETIC AND KINETIC
!     PRESSURE, AND METRIC ELEMENTS ON HALF-GRID

         CALL second0(tbcovon)
         CALL bcovar_par(lu, lv, pxc, ier_flag)
         CALL second0(tbcovoff)
         bcovar_time=bcovar_time+(tbcovoff - tbcovon)

      END IF ACTIVE1
      
      CALL MPI_BCast( ier_flag, 1, MPI_INTEGER, 0, 
     &                RUNVMEC_COMM_WORLD, MPI_ERR) !SAL 070719

      bbuf(1)=irst; bbuf(2)=iequi; bbuf(3)=ivac; bbuf(4)=iter2
      CALL MPI_BCast(bbuf,4,MPI_INTEGER,0,RUNVMEC_COMM_WORLD,MPI_ERR)
      irst=bbuf(1); iequi=bbuf(2); ivac=bbuf(3); iter2=bbuf(4)
      CALL MPI_BCast(lfreeb,1,MPI_LOGICAL,0,RUNVMEC_COMM_WORLD,MPI_ERR)
      IF (ier_flag .ne. norm_term_flag) RETURN !SAL 070719

      IF (irst.EQ.2 .AND. iequi.EQ.0) THEN
         CALL ZEROLASTNTYPE(pgc)
         GOTO 100
      END IF

      timer(tbcov) = timer(tbcov) + (tbcovoff - tbcovon)

!     COMPUTE VACUUM MAGNETIC PRESSURE AT PLASMA EDGE
!     NOTE: FOR FREE BOUNDARY RUNS, THE VALUE OF RBTOR=R*BTOR
!     AT THE PLASMA EDGE SHOULD BE ADJUSTED TO APPROXIMATELY
!     EQUAL THE VACUUM VALUE. THIS CAN BE DONE BY CHANGING
!     EITHER PHIEDGE OR THE INITIAL CROSS SECTION ACCORDING
!     TO THE SCALING LAW  R*BTOR .EQ. PHIEDGE/(R1 * Z1).

      IF (lfreeb       .AND.
     &    iter2 .GT. 1 .AND.
     &    iequi .EQ. 0) THEN

         IF (ictrl_prec2d.LE.1 .AND. (fsqr + fsqz).LE.1.e-3_dp) 
     &      ivac = ivac+1   !decreased from e-1 to e-3 - sph12/04

         IF (nvskip0 .EQ. 0) nvskip0 = MAX(1, nvacskip)

         IVAC0: IF (ivac .GE. 0) THEN
!SPH OFF: 6.20.17
!           IF INITIALLY ON, TURN OFF rcon0, zcon0 SLOWLY
            IF (lactive) THEN
               IF (ictrl_prec2d .EQ. 2) THEN
                  prcon0(:,nsmin:nsmax) = 0;  pzcon0(:,nsmin:nsmax) = 0
               ELSE IF (ictrl_prec2d .EQ. 0) THEN
                  prcon0(:,nsmin:nsmax) = 0.9_dp*prcon0(:,nsmin:nsmax)
                  pzcon0(:,nsmin:nsmax) = 0.9_dp*pzcon0(:,nsmin:nsmax)
               END IF
            ENDIF
            CALL second0 (tvacon)
            ivacskip = MOD(iter2-iter1,nvacskip)
            IF (ivac .LE. 2) ivacskip = 0

!           EXTEND NVACSKIP AS EQUILIBRIUM CONVERGES
            IF (ivacskip .EQ. 0) THEN
               nvacskip = one/MAX(1.e-1_dp, 1.e11_dp*(fsqr+fsqz))
               nvacskip = MAX(nvacskip, nvskip0)
            END IF

!
!           NORMALLY, WHEN COMPUTING THE HESSIAN, IT IS SUFFICIENT TO
!           COMPUTE THE VARIATIONS IN THE "EXACT" SOLUTION, NOT THE ENTIRE
!           FIELD PERIOD SUM. THUS, FOR ictrl_prec2d >= 2, SET ivacskip = 1
!           FOR ictrl_prec2d = 1 (RUN WITH PRECONDITIONER APPLIED), MUST
!           COMPUTE EXACT VACUUM RESPONSE NOW.
!
!           THE EXCEPTION TO THIS IS IF WE ARE TESTING THE HESSIAN (lHess_exact=T), 
!           THEN MUST USE FULL VACUUM CALCULATION TO COMPUTE IT (ivacskip=0)
!
!           lHess_exact = .FALSE.

            IF (ictrl_prec2d .NE. 0) THEN
               IF (lHess_exact .OR. ictrl_prec2d.EQ.2) THEN        !Accurate Hessian
                  ivacskip = 0
               ELSE
                  ivacskip = 1                                     !Fast vacuum calculation used to compute Hessian
               ENDIF
            ENDIF

!          NOTE: pgc contains correct edge values of r,z,l arrays
!                convert_sym, convert_asym have been applied to m=1 modes
          
            CALL convert_par(rmnc,zmns,lmns,rmns,zmnc,lmnc,pgc)

!          DO NOT UPDATE THIS WHEN USING PRECONDITIONER: BREAKS TRI-DIAGONAL STRUCTURE
            IF (ictrl_prec2d.EQ.0 .OR. ictrl_prec2d.EQ.2) THEN
               raxis_nestor(1:nzeta) = pr1(1:nzeta,1,0)
               zaxis_nestor(1:nzeta) = pz1(1:nzeta,1,0)

               ALLOCATE (bcastbuf(2*nzeta))
               bcastbuf(1:nzeta) = raxis_nestor(1:nzeta)
               bcastbuf(nzeta+1:2*nzeta) = zaxis_nestor(1:nzeta)
               CALL second0(tbroadon)
               CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,0,
     &                        RUNVMEC_COMM_WORLD,MPI_ERR)
               CALL second0(tbroadoff)
               broadcast_time = broadcast_time + (tbroadoff - tbroadon)
               raxis_nestor(1:nzeta) = bcastbuf(1:nzeta)
               zaxis_nestor(1:nzeta) = bcastbuf(nzeta+1:2*nzeta)
               DEALLOCATE (bcastbuf)
            END IF

            ALLOCATE (bcastbuf(2))
            bcastbuf(1)=rbtor
            bcastbuf(2)=ctor
            IF (lactive) THEN
               CALL second0(tbroadon)
               CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,
     &                        nranks-1,NS_COMM,MPI_ERR)
               CALL second0(tbroadoff)
               broadcast_time = broadcast_time + (tbroadoff -tbroadon)
            END IF

            CALL second0(tbroadon)
            IF (vlactive) THEN
               CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,0,
     &                        VAC_COMM,MPI_ERR)
            END IF
            CALL second0(tbroadoff)
            broadcast_time = broadcast_time + (tbroadoff -tbroadon)
            rbtor=bcastbuf(1)
            ctor=bcastbuf(2)
            DEALLOCATE (bcastbuf)

            IF (vlactive) THEN
               IF (ictrl_prec2d .NE. 3 .OR.
     &             l_edge) THEN
                  CALL vacuum_par (rmnc, rmns, zmns, zmnc, xm, xn,
     &                             ctor, rbtor, pwint_ns, ns, ivacskip,
     &                             ivac, mnmax, ier_flag, lscreen)
                  IF (ictrl_prec2d .EQ. 2) bsqvac0 = bsqvac
               ELSE
                  bsqvac = bsqvac0; ier_flag = 0
               END IF
            END IF

            IF (vnranks .LT. nranks) THEN
               CALL MPI_Bcast(bsqvac,SIZE(bsqvac),MPI_REAL8,0,
     &                        NS_COMM,MPI_ERR)
            END IF

            IF (ier_flag .NE. 0) THEN
               RETURN
            END IF
!
!          RESET FIRST TIME FOR SOFT START
!
            IF (ivac .EQ. 1) THEN
               irst = 2;  delt0 = delt
               CALL restart_iter(delt0)
               irst = 1
            END IF

!
!          IN CASE PRESSURE IS NOT ZERO AT EXTRAPOLATED EDGE...
!          UNCOMMENT ALL "RPRES" COMMENTS HERE AND IN BCOVAR, FORCES ROUTINES
!          IF NON-VARIATIONAL FORCES ARE DESIRED
!
!          presf_ns = 1.5_dp*pres(ns) - 0.5_dp*pres(ns1)  
!          MUST NOT BREAK TRI-DIAGONAL RADIAL COUPLING: OFFENDS PRECONDITIONER!
            presf_ns = pmass(hs*(ns-1.5_dp))
            IF (presf_ns .NE. zero) THEN
               presf_ns = (pmass(1._dp)/presf_ns) * pres(ns)
            END IF

            DO l = 1, nznt
               bsqsav(l,3) = 1.5_dp*pbzmn_o(l,ns)
     &                     - 0.5_dp*pbzmn_o(l,ns-1)
               pgcon(l,ns) = bsqvac(l) + presf_ns
               rbsq(l) = pgcon(l,ns)*(pr1(l,ns,0) + pr1(l,ns,1))*ohs
               dbsq(l) = ABS(pgcon(l,ns)-bsqsav(l,3))
            END DO

            IF (ivac .EQ. 1) THEN
               IF (vlactive) THEN
                  bsqsav(:nznt,1) = pbzmn_o(:,ns)
                  bsqsav(:nznt,2) = bsqvac(:nznt)
                  CALL MPI_Bcast(bsqsav(:,1),nznt,MPI_REAL8,
     &                           nranks-1,NS_COMM,MPI_ERR)
               END IF
            ELSE IF (ictrl_prec2d .NE. 3) THEN
               CALL MPI_Bcast(bsqsav(:,1),nznt,MPI_REAL8,
     &                        0,NS_COMM,MPI_ERR)
            END IF

            CALL second0 (tvacoff)
            timer(tvac) = timer(tvac) + (tvacoff - tvacon)
            IF (ictrl_prec2d .GE. 2) THEN
               timer(tvac_2d) = timer(tvac_2d)+ (tvacoff - tvacon)
            END IF
         END IF IVAC0
      END IF
!
!     COMPUTE CONSTRAINT FORCE
!
      ACTIVE2: IF (lactive) THEN
         IF (iequi .NE. 1) THEN
            DO l = nsmin, nsmax
               pextra1(:,l,0) = (prcon(:,l,0) - prcon0(:,l))*pru0(:,l)
     &                        + (pzcon(:,l,0) - pzcon0(:,l))*pzu0(:,l)
            END DO
            CALL alias_par (pgcon, pextra1(:,:,0), pgc, pgc(1+mns),
     &                      pgc(1+2*mns), pextra1(:,:,1))
         ELSE
            IF (lrecon) THEN
               pxc(:ns) = pxc(:ns) + delr_mse
            END IF
            GOTO 100
         END IF

!
!     COMPUTE MHD FORCES ON INTEGER-MESH
!
         CALL forces_par

!     SYMMETRIZE FORCES (in u-v space)
!
         IF (lasym) THEN
            CALL symforce_par (parmn, pbrmn, pcrmn, pazmn, pbzmn,
     &                         pczmn, pblmn, pclmn, prcon, pzcon, pr1,
     &                         pru, prv, pz1, pzu, pzv, pextra3,
     &                         pextra4, pextra1, pextra2)
         END IF
     
!
!     FOURIER-TRANSFORM MHD FORCES TO (M,N)-SPACE
!
         CALL tomnsps_par (pgc, parmn, pbrmn, pcrmn, pazmn, pbzmn,
     &                     pczmn, pblmn, pclmn, prcon, pzcon)

         IF (lasym) THEN
            CALL tomnspa_par (pgc, pr1, pru, prv, pz1, pzu, pzv,
     &                        pextra3, pextra4, pextra1, pextra2)
         END IF

!================================================================
!
!     COMPUTE FORCE RESIDUALS (RAW AND PRECONDITIONED)
!
!================================================================
         CALL second0 (treson)

         CALL SAXLASTNTYPE(pgc, pscalxc, pgc)

         CALL residue_par(pgc, pgc(1+irzloff), pgc(1+2*irzloff))

      END IF ACTIVE2

!NEED THIS ON ALL PROCESSORS IN GROUP (NOT JUST ACTIVE ONES) FOR STOPPING CRITERION IN EVOLVE
#if defined(MPI_OPT)
      IF (gnranks .GT. nranks) THEN
         ALLOCATE(bcastbuf(6))
         bcastbuf(1) = fsqr; bcastbuf(2) = fsqr1
         bcastbuf(3) = fsqz; bcastbuf(4) = fsqz1
         bcastbuf(5) = fsql; bcastbuf(6) = fsql1
         CALL second0(tbroadon)
         CALL MPI_Bcast(bcastbuf,SIZE(bcastbuf),MPI_REAL8,0,
     &                  RUNVMEC_COMM_WORLD,MPI_ERR)
         CALL second0(tbroadoff)
         broadcast_time = broadcast_time + (tbroadoff -tbroadon)
         fsqr = bcastbuf(1); fsqr1 = bcastbuf(2)
         fsqz = bcastbuf(3); fsqz1 = bcastbuf(4)
         fsql = bcastbuf(5); fsql1 = bcastbuf(6)
         DEALLOCATE(bcastbuf)
      END IF
#endif

!     Force new initial axis guess IF ALLOWED (l_moveaxis=T)
      IF (lmove_axis                  .and.
     &    iter2               .eq  .1 .and.
     &    (fsqr + fsqz +fsql) .gt. 1.E2_dp) THEN
         irst = 4
      END IF

      CALL second0 (tresoff)
      timer(tres) = timer(tres) + (tresoff - treson)

 100  CONTINUE

      CALL second0 (tfunoff)
      timer(tfun) = timer(tfun) + (tfunoff - tfunon)
      IF (ictrl_prec2d .GE. 2) THEN
         timer(tfun_2d) = timer(tfun_2d) + (tfunoff - tfunon)
      END IF

      END SUBROUTINE funct3d_par
