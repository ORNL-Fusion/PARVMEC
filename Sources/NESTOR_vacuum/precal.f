      SUBROUTINE precal (wint)
      USE vparams, ONLY: zero, one, epstan
      USE vacmod
      USE vmec_main, ONLY: mnmax
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), INTENT(in) :: wint(nuv3)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(dp), PARAMETER :: p25 = p5*p5, bigno = 1.e50_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: kp, ku, kuminus, kv, kvminus, i, m, n, mn, n1,
     1   imn, jmn, kmn, l, istat1, smn, nuv_tan, ndim, q, qq
      REAL(dp), DIMENSION(0:mf + nf,0:mf,0:nf) :: cmn
      REAL(dp), DIMENSION(0:mf + nf,0:mf + nf) :: cosk
      REAL(dp), DIMENSION(0:mf + nf) :: pnode
      REAL(dp) :: argu, argv, argp, dn1, alp_per, tnode, xnode,
     1            pjm1, pj, pjp1, tjk, cnext, cthis, cprev,
     2            tprecon, tprecoff
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
         nvper = 64 
         nuv_tan = 2*nu*nvper
      ELSE
         nvper = nfper
         nuv_tan = 2*nuv
      END IF
      
      alp_per = pi2/nvper
      nvp = nv*nvper

      ALLOCATE (tanu(nuv_tan), tanv(nuv_tan), 
     1     sinper(nvper), cosper(nvper), sinuv(nuv), cosuv(nuv),
     2     sinu(0:mf,nu), cosu(0:mf,nu), sinv(-nf:nf,nv),
     3     cosv(-nf:nf,nv), sinui(0:mf,nu2), cosui(0:mf,nu2),
     4     cmns(0:(mf+nf),0:mf,0:nf), csign(-nf:nf),
     5     sinu1(nuv3,0:mf), cosu1(nuv3,0:mf),
     6     sinv1(nuv3,0:nf), cosv1(nuv3,0:nf), imirr(nuv),
     7     xmpot(mnpd), xnpot(mnpd), stat=istat1)
      IF (istat1.ne.0) STOP 'allocation error in precal'


!
!     IMIRR(I) GIVES THE INDEX OF THE POINT TWOPI-THETA(I),TWOPI-ZETA(I)
!
      DO kp = 1, nvper
         cosper(kp) = COS(alp_per*(kp - 1))
         sinper(kp) = SIN(alp_per*(kp - 1))
      END DO

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

!
!     NOTE: ANGLE DIFFERENCE IS PI2*{[NUV + (KUP-1)] - (KU-1)}
!           THIS DIFFERENCE IS ACCOUNTED FOR BY THE OFFSET IUOFF IN GREENF ROUTINE
!
!           THE KP SUM BELOW IS USED ONLY FOR NV == 1. IT PERFORMS THE V-INTEGRAL
!           IN AN AXISYMMETRIC PLASMA
!
      i = 0
      DO kp = 1, nvper
         IF (kp.gt.1 .and. nv.ne.1) EXIT
         argp = p5*alp_per*(kp-1)
         DO ku = 1, 2*nu
            argu = p5*alu*(ku - 1)
            DO kv = 1, nv
               i = i + 1
               argv = p5*alv*(kv - 1) + argp
               IF (ABS(argu - p25*pi2)<epstan .or.
     1         ABS(argu - 0.75_dp*pi2) < epstan) THEN
                  tanu(i) = bigno
               ELSE
                  tanu(i) = 2*TAN(argu)
               ENDIF
               IF (ABS(argv - p25*pi2) < epstan) THEN
                  tanv(i) = bigno
               ELSE
                  tanv(i) = 2*TAN(argv)
               ENDIF
            END DO
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
               IF (MOD(MAX(0,n - m),2) .eq. 1) cmn(i,m,n) = -cmn(i,m,n)
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
