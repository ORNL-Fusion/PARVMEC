      SUBROUTINE precal (wint)
      USE vparams, ONLY: zero, one
      USE vacmod
      USE vmec_main, ONLY: mnmax
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), INTENT(in) :: wint(nuv3)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!>  Index of a poloidal grid point.
      INTEGER                                  :: ku
!>  Index of the poloidal grid point at the mirrored angle.
      INTEGER                                  :: kuminus
!>  Index of a toroidal grid point.
      INTEGER                                  :: kv
!>  Index of the toroidal grid point at the mirrored angle.
      INTEGER                                  :: kvminus
!>  Index of a grid point of the surface or of a Chebyshev node, and status
!>  of an allocation.
      INTEGER                                  :: i
!>  Poloidal mode number.
      INTEGER                                  :: m
!>  Toroidal mode number.
      INTEGER                                  :: n
!>  Index of a mode (m,n).
      INTEGER                                  :: mn
!>  Absolute value of the toroidal mode number.
      INTEGER                                  :: n1
!>  Offset of the first grid point of this process.
      INTEGER                                  :: imn
!>  Absolute value of m - n.
      INTEGER                                  :: kmn
!>  Order of a Chebyshev polynomial.
      INTEGER                                  :: l
!>  Status of an allocation.
      INTEGER                                  :: istat1
!>  Smaller of m and n.
      INTEGER                                  :: smn
!>  Degree in the recurrence of the Jacobi polynomial.
      INTEGER                                  :: q
!>  Values at the Chebyshev nodes of the polynomial of each (m,n).
      REAL(dp), DIMENSION(0:mf + nf,0:mf,0:nf) :: cmn
!>  Chebyshev polynomials at the nodes.
      REAL(dp), DIMENSION(0:mf + nf,0:mf + nf) :: cosk
!>  Values of one polynomial at the nodes.
      REAL(dp), DIMENSION(0:mf + nf)           :: pnode
!>  Angle of a Chebyshev node.
      REAL(dp)                                 :: argu
!>  Toroidal mode number times the grid step of v.
      REAL(dp)                                 :: dn1
!>  Chebyshev node.
      REAL(dp)                                 :: tnode
!>  Argument of the Jacobi polynomial at a node.
      REAL(dp)                                 :: xnode
!>  Jacobi polynomial of the degree below the current one.
      REAL(dp)                                 :: pjm1
!>  Jacobi polynomial of the current degree.
      REAL(dp)                                 :: pj
!>  Jacobi polynomial of the degree above the current one.
      REAL(dp)                                 :: pjp1
!>  Twice the degree plus kmn in the recurrence of the Jacobi polynomial.
      REAL(dp)                                 :: tjk
!>  Factor of the next polynomial in the recurrence.
      REAL(dp)                                 :: cnext
!>  Factor of the current polynomial in the recurrence.
      REAL(dp)                                 :: cthis
!>  Factor of the previous polynomial in the recurrence.
      REAL(dp)                                 :: cprev
!>  Time at the start of the routine.
      REAL(dp)                                 :: tprecon
!>  Time at the end of the routine.
      REAL(dp)                                 :: tprecoff
C-----------------------------------------------
!
!     THIS ROUTINE COMPUTES INITIAL CONSTANTS AND ARRAYS
!     NOTE: alu*alv = (2*pi)**2 * wint
!
      CALL second0(tprecon)

      pi2 = 8*ATAN(one)
      pi3 = p5*pi2**3
      pi4 = 2*pi2
      onp = one/nfper
      onp2 = onp*onp
      alu = pi2/nu
      alv = pi2/nv
      alp = pi2*onp
      alvp = onp*alv

!
!     ALLOCATE PERSISTENT ARRAYS. DEALLOCATED IN FILEOUT ROUTINE
!
      IF (nv == 1) THEN         !(AXISYMMETRIC CASE: DO FP SUM TO INTEGRATE IN V)
         CALL image_tables(64)
      ELSE
         CALL image_tables(nfper)
      END IF

      ALLOCATE (sinuv(nuv), cosuv(nuv),
     &          sinu(0:mf,nu), cosu(0:mf,nu), sinv(-nf:nf,nv),
     &          cosv(-nf:nf,nv), sinui(0:mf,nu2), cosui(0:mf,nu2),
     &          cmns(0:(mf+nf),0:mf,0:nf), csign(-nf:nf),
     &          sinu1(nuv3,0:mf), cosu1(nuv3,0:mf),
     &          sinv1(nuv3,0:nf), cosv1(nuv3,0:nf), imirr(nuv),
     &          xmpot(mnpd), xnpot(mnpd), stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error in precal'


!
!     IMIRR(I) GIVES THE INDEX OF THE POINT TWOPI-THETA(I),TWOPI-ZETA(I)
!
      DO ku = 1, nu
         kuminus = MOD(nu + 1 - ku,nu) + 1
         DO kv = 1, nv
            kvminus = MOD(nv + 1 - kv,nv) + 1
            i = kv + nv*(ku - 1)
            imirr(i) = kvminus + nv*(kuminus - 1)
            cosuv(i) = COS(alvp*(kv - 1))
            sinuv(i) = SIN(alvp*(kv - 1))
         END DO
      END DO

      DO m = 0, mf
         l40: DO ku = 1, nu
            cosu(m,ku) = COS(alu*(m*(ku - 1)))
            sinu(m,ku) = SIN(alu*(m*(ku - 1)))
            DO kv = 1, nv
               i = kv + nv*(ku - 1)
               IF (i > nuv3) CYCLE  l40
               cosu1(i,m) = cosu(m,ku)
               sinu1(i,m) = sinu(m,ku)
            END DO
         END DO l40
         DO ku = 1, nu2
            cosui(m,ku) = cosu(m,ku)*alu*alv*2   
            sinui(m,ku) = sinu(m,ku)*alu*alv*2
            IF (ku.eq.1 .or. ku.eq.nu2) cosui(m,ku) = p5*cosui(m,ku)
         END DO
      END DO

      DO n = -nf, nf
         dn1 = alvp*(n*nfper)
         csign(n) = SIGN(one,dn1)
         l50: DO ku = 1, nu
            DO kv = 1, nv
               i = kv + nv*(ku - 1)
               cosv(n,kv) = COS(dn1*(kv - 1))
               sinv(n,kv) = SIN(dn1*(kv - 1))
               IF (i.gt.nuv3 .or. n.lt.0) CYCLE  l50
               cosv1(i,n) = cosv(n,kv)
               sinv1(i,n) = sinv(n,kv)
            END DO
         END DO l50
      END DO

      mn = 0
      imn = nuv3min-1
      numjs_vac=nuv3max-nuv3min+1

      ALLOCATE(sinmni(numjs_vac,mnpd), cosmni(numjs_vac,mnpd),stat=i)
      IF (i .NE. 0) STOP 'Allocation error in scalpot'
      DO n = -nf, nf
         n1 = ABS(n)
         DO m = 0, mf
            mn = mn + 1
            xmpot(mn) = m
            xnpot(mn) = n*nfper
            DO i = nuv3min, nuv3max
               sinmni(i-imn,mn) = wint(i)*(sinu1(i,m)*cosv1(i,n1)
     1                  - csign(n)*cosu1(i,m)*sinv1(i,n1))*(pi2*pi2)
               cosmni(i-imn,mn) = wint(i)*(cosu1(i,m)*cosv1(i,n1)
     1                  + csign(n)*sinu1(i,m)*sinv1(i,n1))*(pi2*pi2)
            END DO
         END DO
      END DO
!
!     COMPUTE CMNS, THE COEFFICIENTS OF THE POLYNOMIALS OF EQ (A13) IN THE J.COMP.PHYS PAPER (PKM),
!     IN THE CHEBYSHEV BASIS: SUM_L C(L;M,N) t**L = SUM_K CMNS(K,M,N) T_K(t).
!
!     THE C(L;M,N) OF EQ (A14), WITH L = 2L' + Kmn (Kmn = |m-n|) AND Kmn <= L <= m+n, ARE THE
!     COEFFICIENTS OF THE RADIAL ZERNIKE POLYNOMIAL
!
!        (-1)**MAX(0,n-m) * t**Kmn * P(1 - 2 t**2),
!
!     P BEING THE JACOBI POLYNOMIAL OF DEGREE MIN(m,n) WITH PARAMETERS (Kmn, 0). IT IS EVALUATED
!     FROM THE THREE-TERM RECURRENCE OF P IN THE DEGREE AT THE mf+nf+1 CHEBYSHEV-GAUSS NODES, WHICH
!     DETERMINE A POLYNOMIAL OF DEGREE m+n <= mf+nf EXACTLY: CMN(I,M,N) HOLDS ITS VALUE AT NODE I,
!     AND COSK(K,I) = COS(K THETA_I).
!
      DO i = 0, mf + nf
         argu = (i + p5)*(p5*pi2)/(mf + nf + 1)
         DO l = 0, mf + nf
            cosk(l,i) = COS(l*argu)
         END DO
         tnode = cosk(1,i)
         xnode = one - 2*tnode*tnode
         DO m = 0, mf
            DO n = 0, nf
               kmn = ABS(m - n)
               smn = MIN(m, n)
               pjm1 = 0
               pj = 1
               IF (smn .ge. 1) THEN
                  pjm1 = 1
                  pj = p5*((kmn + 2)*xnode + kmn)
                  DO q = 1, smn - 1
                     tjk = 2*q + kmn
                     cnext = 2*(q + one)*(q + kmn + one)*tjk
                     cthis = (tjk + 1)*((tjk + 2)*tjk*xnode + kmn*kmn)
                     cprev = 2*(q + kmn)*(tjk + 2)*q
                     pjp1 = (cthis*pj - cprev*pjm1)/cnext
                     pjm1 = pj
                     pj = pjp1
                  END DO
               END IF
               cmn(i,m,n) = pj*tnode**kmn
               IF (MOD(MAX(0,n - m),2) .eq. 1) THEN
                  cmn(i,m,n) = -cmn(i,m,n)
               END IF
            END DO
         END DO
      END DO
!
!     Now combine these into a single polynomial (cmns), Eq. A13), still as values at the nodes.
!     NOTE:  The ALP=2*pi/nfper factor is needed to normalize integral over field periods
!
      DO m = 1,mf
         DO n = 1,nf
            cmns(0:mf+nf,m,n) = p5*alp*(cmn(0:mf+nf,m,n) +
     1      cmn(0:mf+nf,m-1,n) + cmn(0:mf+nf,m,n-1) +
     2      cmn(0:mf+nf,m-1,n-1))
         END DO
      END DO
      cmns(0:mf+nf,1:mf,0) = (p5*alp)*(cmn(0:mf+nf,1:mf,0)
     1                     +           cmn(0:mf+nf,:mf-1,0))
      cmns(0:mf+nf,0,1:nf) = (p5*alp)*(cmn(0:mf+nf,0,1:nf)
     1                     +           cmn(0:mf+nf,0,:nf-1))
      cmns(0:mf+nf,0,0)    = (p5*alp)*(cmn(0:mf+nf,0,0)
     1                     +           cmn(0:mf+nf,0,0))
!
!     CHEBYSHEV COEFFICIENTS FROM THE VALUES AT THE NODES (DISCRETE COSINE TRANSFORM). THE
!     POLYNOMIAL OF (m,n) HAS DEGREE m+n, SO THE COEFFICIENTS ABOVE IT VANISH.
!
      DO m = 0, mf
         DO n = 0, nf
            pnode = cmns(:,m,n)
            cmns(:,m,n) = 0
            DO l = 0, m + n
               cmns(l,m,n) = 2*SUM(pnode*cosk(l,:))/(mf + nf + 1)
            END DO
            cmns(0,m,n) = p5*cmns(0,m,n)
         END DO
      END DO

      numjs_vac=nuv3max-nuv3min+1
!      blksize_scp=mnpd2  

      ALLOCATE (counts_vac(vnranks),disps_vac(vnranks), stat=i)
      IF (i .NE. 0) STOP 'Allocation error in precal'
      DO i = 1, vnranks
         counts_vac(i) = nuv3max_arr(i) - nuv3min_arr(i) + 1
      END DO
      disps_vac(1)=0
      DO i = 2, vnranks
         disps_vac(i) = disps_vac(i - 1) + counts_vac(i - 1)
      END DO

      CALL second0(tprecoff)
      precal_time = precal_time + (tprecoff - tprecon) 

      END SUBROUTINE precal

!-------------------------------------------------------------------------------
!>  @brief Set the number of toroidal images and the tables that depend on it.
!>
!>  nvper is the number of field periods that greenf sums over or, for nv = 1,
!>  the number of images of the evaluation point around the torus over which
!>  it takes the integral in v. tanu and tanv hold 2 TAN(u/2) and 2 TAN(v/2)
!>  of the analytic approximation for the differences of the grid angles, and
!>  sinper and cosper rotate into each period or image.
!>
!>  @param[in] nimages Number of field periods or of toroidal images.
!-------------------------------------------------------------------------------
      SUBROUTINE image_tables(nimages)
      USE vparams, ONLY: epstan
      USE vacmod

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(IN) :: nimages

!  Local Variables
!>  Index of a field period or of an image.
      INTEGER             :: kp
!>  Index of a half step of u.
      INTEGER             :: ku
!>  Index of a half step of v.
      INTEGER             :: kv
!>  Index into the tables of tangents.
      INTEGER             :: i
!>  Size of the tables of tangents.
      INTEGER             :: nuv_tan
!>  Status of the allocation.
      INTEGER             :: istat1
!>  Half the difference of two poloidal grid angles.
      REAL(dp)            :: argu
!>  Half the difference of two toroidal grid angles.
      REAL(dp)            :: argv
!>  Half the toroidal angle of an image.
      REAL(dp)            :: argp
!>  Toroidal angle between two field periods or images.
      REAL(dp)            :: alp_per

!  Local Parameters
!>  One quarter.
      REAL(dp), PARAMETER :: p25 = p5*p5
!>  Stands for the tangent of a quarter turn.
      REAL(dp), PARAMETER :: bigno = 1.e50_dp

!  Start of executable code
      nvper = nimages
      nvp = nv*nvper
      alp_per = pi2/nvper
      IF (nv .EQ. 1) THEN
         nuv_tan = 2*nu*nvper
      ELSE
         nuv_tan = 2*nuv
      END IF

      IF (ALLOCATED(tanu)) THEN
         DEALLOCATE (tanu, tanv, sinper, cosper)
      END IF
      ALLOCATE (tanu(nuv_tan), tanv(nuv_tan), sinper(nvper),
     &          cosper(nvper), stat=istat1)
      IF (istat1 .NE. 0) THEN
         STOP 'allocation error in image_tables'
      END IF

      DO kp = 1, nvper
         cosper(kp) = COS(alp_per*(kp - 1))
         sinper(kp) = SIN(alp_per*(kp - 1))
      END DO

!  Note: the angle difference is pi2*{[nuv + (kup - 1)] - (ku - 1)}. This
!  difference is accounted for by the offset iuoff in the greenf routine.

!  The kp sum below is used only for nv == 1. It performs the v-integral in
!  an axisymmetric plasma.
      i = 0
      DO kp = 1, nvper
         IF (kp .GT. 1 .AND. nv .NE. 1) THEN
            EXIT
         END IF
         argp = p5*alp_per*(kp - 1)
         DO ku = 1, 2*nu
            argu = p5*alu*(ku - 1)
            DO kv = 1, nv
               i = i + 1
               argv = p5*alv*(kv - 1) + argp
               IF (ABS(argu - p25*pi2) .LT. epstan .OR.
     &             ABS(argu - 0.75_dp*pi2) .LT. epstan) THEN
                  tanu(i) = bigno
               ELSE
                  tanu(i) = 2*TAN(argu)
               END IF
               IF (ABS(argv - p25*pi2) .LT. epstan) THEN
                  tanv(i) = bigno
               ELSE
                  tanv(i) = 2*TAN(argv)
               END IF
            END DO
         END DO
      END DO

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Resolve the cross section with the toroidal images of nv = 1.
!>
!>  For nv = 1, greenf takes the integral in v as a sum over nvper images
!>  around the torus. Where the cross section has the radius of curvature rho,
!>  the difference of the Green's function and its analytic approximation
!>  varies over the toroidal angle rho/R, and the sum resolves it with
!>  images_per_radius images for each unit of R/rho. The count is that of the
!>  largest R/rho on the boundary, in steps of 64, and it only grows.
!>
!>  In the quantities of surface, R/rho is 2 ABS(auu)/guu_b**1.5.
!-------------------------------------------------------------------------------
      SUBROUTINE image_count
      USE vacmod
      USE parallel_include_module

      IMPLICIT NONE

!  Local Variables
!>  Number of images that resolves the cross section.
      INTEGER             :: nimages
!>  Largest ratio of R to the radius of curvature on the boundary.
      REAL(dp)            :: curv

!  Local Parameters
!>  Images for each unit of that ratio.
      REAL(dp), PARAMETER :: images_per_radius = 40

!  Start of executable code
      curv = MAXVAL(2*ABS(auu(nuv3min:nuv3max)) /
     &              SQRT(guu_b(nuv3min:nuv3max))**3)
      IF (vlactive) THEN
         CALL MPI_Allreduce(MPI_IN_PLACE, curv, 1, MPI_REAL8, MPI_MAX,
     &                      VAC_COMM, MPI_ERR)
      END IF

      nimages = 64*CEILING(images_per_radius*curv/64)
      IF (nimages .GT. nvper) THEN
         CALL image_tables(nimages)
      END IF

      END SUBROUTINE
