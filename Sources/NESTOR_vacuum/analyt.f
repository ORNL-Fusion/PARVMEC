!*******************************************************************************
!>  @file analyt.f
!>  @brief Contains module @ref analytic.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines routines for compute analytic integration that accounds for the
!>  singularity.
!*******************************************************************************
      MODULE analytic
      USE stel_kinds, ONLY: dp

      IMPLICIT NONE

!  Factors of the five-term recurrence of the Chebyshev moments, indexed by
!  the order k. See moment_factors.
      REAL(dp), DIMENSION(:), ALLOCATABLE, PRIVATE :: fsub2
      REAL(dp), DIMENSION(:), ALLOCATABLE, PRIVATE :: fsub1
      REAL(dp), DIMENSION(:), ALLOCATABLE, PRIVATE :: fdiag
      REAL(dp), DIMENSION(:), ALLOCATABLE, PRIVATE :: fsup1
      REAL(dp), DIMENSION(:), ALLOCATABLE, PRIVATE :: fsup2
      REAL(dp), DIMENSION(:), ALLOCATABLE, PRIVATE :: frhs

      CONTAINS
!-------------------------------------------------------------------------------
!>  @brief Main routine
!>
!>  @param[out] grpmn
!>  @param[out] grpmn
!>  @param[in]  ivacskip
!>  @param[in]  ivacskip
!-------------------------------------------------------------------------------
      SUBROUTINE analyt(grpmn, bvec, ivacskip, ndim)
      USE vacmod
      USE parallel_include_module
      USE timer_sub

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), INTENT(OUT) :: grpmn(mnpd2,nuv3)
      REAL(dp), INTENT(OUT) :: bvec(mnpd,ndim)
      INTEGER, INTENT(IN)   :: ivacskip
      INTEGER, INTENT(IN)   :: ndim

!  local variables
      INTEGER :: l, n, m, k
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: tlp, tlm, slp, slm
      REAL(dp), DIMENSION(:), ALLOCATABLE ::
     &   r0p, r1p, r0m, r1m, sqrtc, sqrta, adp, adm, cma, ra1p, ra1m,
     &   tlps, tlms, slps, slms, tlpm, slpm, ulp1, ulp2, ulm1, ulm2
      REAL(dp) :: sign1, tanalon, tanaloff, ulp, ulm, wl
      REAL(dp) :: sqad1u
      REAL(dp) :: sqad2u
      REAL(dp) :: delt1u
      REAL(dp) :: azp1u
      REAL(dp) :: azm1u
      REAL(dp) :: cma11u
      REAL(dp) :: tlp2
      REAL(dp) :: tlm2

!  Start of executable code
      CALL second0(tanalon)

      ALLOCATE (r0p(nuv3min:nuv3max), r1p(nuv3min:nuv3max),
     &          r0m(nuv3min:nuv3max), r1m(nuv3min:nuv3max),
     &          sqrtc(nuv3min:nuv3max), sqrta(nuv3min:nuv3max),
     &          adp(nuv3min:nuv3max), adm(nuv3min:nuv3max),
     &          cma(nuv3min:nuv3max), ra1p(nuv3min:nuv3max),
     &          ra1m(nuv3min:nuv3max), slpm(nuv3min:nuv3max),
     &          tlpm(nuv3min:nuv3max),
     &          tlp(nuv3min:nuv3max,0:nf + mf),
     &          tlm(nuv3min:nuv3max,0:nf + mf),
     &          slp(nuv3min:nuv3max,0:nf + mf),
     &          slm(nuv3min:nuv3max,0:nf + mf),
     &          tlps(nuv3min:nuv3max), tlms(nuv3min:nuv3max),
     &          slps(nuv3min:nuv3max), slms(nuv3min:nuv3max),
     &          ulp1(nuv3min:nuv3max), ulp2(nuv3min:nuv3max),
     &          ulm1(nuv3min:nuv3max), ulm2(nuv3min:nuv3max), stat = l)
      IF (l .ne. 0) THEN
         STOP 'Allocation error in SUBROUTINE analyt'
      ENDIF

!
!     ALL EQUATIONS REFER TO THE PAPER BY P. MERKEL (PKM)
!     IN J. COMPUT. PHYSICS 66, p83 (1986)
!
!     IN GOING BETWEEN THE COMPLEX NOTATION OF (PKM) AND OUR REAL FORM,
!     NOTE THAT THE INTEGRALS (APPENDIX, PKM) Imn AND Kmn ARE BOTH REAL.
!     THUS, THE SIN(mu-nv) INTEGRALS OF THE SINGULAR PIECE (ANALYTIC CONTRIBUTION)
!     VANISHES.
!
!     THE REQUIRED SOURCE-TERM INTEGRALS ARE (Eq.2.16-2.17):
!
!     BVECS(m,n) = Int< SIN(mu' - nv') han(u',v') >
!     BVECC(m,n) = Int< COS(mu' - nv') han(u',v') >
!
!     Where Int<...> means integration over u (theta) and v (zeta) and
!     summation over field periods. These can be written in terms of PKM integrals
!     Imn(a,b,c), where a(u,v) = guu (g theta-theta), etc.:
!
!     BVECS(m,n) = ALP * Int<SIN(mu' - nv') * F * Im,-n(a,b,c)>
!     BVECC(m,n) = ALP * Int<COS(mu' - nv') * F * Im,-n(a,b,c)>
!
!     Here, F = - BNORM(u',v') is defined in Eq.(2.13), and ALP = (2*pi/nfp).
!
!     Similarly, the analytic part of the matrix A(m,n;m',n') can be written:
!
!     A(m,n;m',n') = (2*pi/nfp) * Int<SIN(mu' - nv')*SIN(m'u' - n'v')
!                              [Km,-n](a',b',c';A',B',C')>
!
!     On EXIT, GRPMN(ip,m,n) = ALP * SIN(ip,m,n) * K[m,-n](ip)
!
!
!     COMPUTE ALL QUANTITIES INDEPENDENT OF THE MODE INDICES L,M,N
!     NOTE: 2b = guv_b HAS FACTOR OF 2 BUILT IN (see SUBROUTINE SURFACE)
!
!     ADP(M): a +(-)2b + c
!     CMA:    c - a
!     DELTA:  4*(ac - b**2)
!     AZP(M): A +(-)2*B + C
!     CMA1:   C - A
!     R1P(M): Coefficient of l*Tl+(-) in eq (A17)
!     R0P(M): Coefficient of l*T(l-1)+(-) in eq (A17)
!     RA1P(M):Coefficient of Tl+(-) in eq (A17)
!
      DO k = nuv3min, nuv3max
         adp(k) = guu_b(k) + guv_b(k) + gvv_b(k) 
         adm(k) = guu_b(k) - guv_b(k) + gvv_b(k) 
         cma(k) = gvv_b(k) - guu_b(k) 
         sqrtc(k)   = two*SQRT(gvv_b(k))
         sqrta(k)   = two*SQRT(guu_b(k))
      END DO

      IF (ivacskip .EQ. 0) THEN

         grpmn(:,nuv3min:nuv3max) = 0

         DO k = nuv3min, nuv3max
            delt1u = adp(k)*adm(k) - cma(k)*cma(k)
            azp1u = auu(k) + auv(k) + avv(k)
            azm1u = auu(k) - auv(k) + avv(k)
            cma11u = avv(k) - auu(k)
            r1p(k) = (azp1u*(delt1u - cma(k)*cma(k))/adp(k)
     &             -  azm1u*adp(k) + two*cma11u*cma(k))/delt1u
            r1m(k) = (azm1u*(delt1u - cma(k)*cma(k))/adm(k)
     &             -  azp1u*adm(k) + two*cma11u*cma(k))/delt1u
            r0p(k) = (-azp1u*adm(k)*cma(k)/adp(k) - azm1u*cma(k)
     &             +  two*cma11u*adm(k))/delt1u
            r0m(k)  = (-azm1u*adp(k)*cma(k)/adm(k) - azp1u*cma(k)
     &              + two*cma11u*adp(k))/delt1u
            ra1p(k) = azp1u/adp(k)
            ra1m(k) = azm1u/adm(k)
         END DO
      ENDIF

!
!     INITIALIZE VECTORS
!
      bvec = 0
!
!     THE POLYNOMIALS OF EQ (A13-A14) ARE EXPANDED IN CHEBYSHEV POLYNOMIALS
!     T_L (SEE PRECAL), SO TLP(M) HOLD THE CHEBYSHEV MOMENTS OF THE KERNEL,
!     Int[-1,1] T_L(t)/SQRT(adp(m) t**2 + 2 cma t + adm(p)) dt,
!     IN PLACE OF THE MONOMIAL MOMENTS TL+(-). T- EXCHANGES ADP AND ADM,
!     WHICH IS guv_b -> -guv_b.
!
      CALL chebyshev_moments(guu_b(nuv3min:nuv3max),
     &                       guv_b(nuv3min:nuv3max),
     &                       gvv_b(nuv3min:nuv3max), tlp)
      CALL chebyshev_moments(guu_b(nuv3min:nuv3max),
     &                       -guv_b(nuv3min:nuv3max),
     &                       gvv_b(nuv3min:nuv3max), tlm)
!
!     COMPUTE SL+ and SL- , Eq (A17) APPLIED TO T_L
!     SLP(M): SL+(-)
!
!     Eq (A17) MAPS THE POLYNOMIAL p TO
!        R1 Int[t p'/SQRT(Q)] + RA1 Int[p/SQRT(Q)] + R0 Int[p'/SQRT(Q)]
!        - (R0 + R1) p(1)/SQRTC + (R0 - R1) p(-1)/SQRTA.
!     FOR p = T_L: T_L' = L U_(L-1) AND t U_(L-1) = (U_L + U_(L-2))/2. ULP(M)
!     ACCUMULATE THE MOMENTS OF U_L = 2 T_L + U_(L-2), WITH ULP(M)1 AND
!     ULP(M)2 THOSE OF U_(L-1) AND U_(L-2).
!
      IF (ivacskip .eq. 0) THEN
         ulp1 = 0
         ulp2 = 0
         ulm1 = 0
         ulm2 = 0
         sign1 = 1
         wl = 1
         DO l = 0, mf + nf
            DO k = nuv3min, nuv3max
               ulp = wl*tlp(k,l) + ulp2(k)
               ulm = wl*tlm(k,l) + ulm2(k)
               slp(k,l) = r1p(k)*l*p5*(ulp + ulp2(k))
     &                  + ra1p(k)*tlp(k,l)
     &                  + r0p(k)*l*ulp1(k)
     &                  - (r1p(k) + r0p(k))/sqrtc(k)
     &                  + sign1*(r0p(k) - r1p(k))/sqrta(k)
               slm(k,l) = r1m(k)*l*p5*(ulm + ulm2(k))
     &                  + ra1m(k)*tlm(k,l)
     &                  + r0m(k)*l*ulm1(k)
     &                  - (r1m(k) + r0m(k))/sqrtc(k)
     &                  + sign1*(r0m(k) - r1m(k))/sqrta(k)
               ulp2(k) = ulp1(k)
               ulp1(k) = ulp
               ulm2(k) = ulm1(k)
               ulm1(k) = ulm
            END DO
            sign1 = -sign1
            wl = 2
         END DO
      ENDIF
!
!     BEGIN MODE NUMBER (m,n) LOOP. THE L-SUM OF EQ (A14) TO COMPUTE THE Imn
!     (and Kmn) INTEGRALS RUNS OVER THE ORDER OF THE CHEBYSHEV POLYNOMIALS AND
!     IS TAKEN FIRST; CMNS(L,M,N) VANISHES FOR L > M+N.
!     TLPS(M): SUM_L CMNS(L,M,N) TL+(-), SLPS(M): SUM_L CMNS(L,M,N) SL+(-)
!
      DO n = 0, nf
         DO m = 0, mf
            tlps = 0
            tlms = 0
            DO l = 0, m + n
               tlps = tlps + cmns(l,m,n)*tlp(:,l)
               tlms = tlms + cmns(l,m,n)*tlm(:,l)
            END DO
            slps = 0
            slms = 0
            IF (ivacskip .eq. 0) THEN
               DO l = 0, m + n
                  slps = slps + cmns(l,m,n)*slp(:,l)
                  slms = slms + cmns(l,m,n)*slm(:,l)
               END DO
            ENDIF

            IF (n.eq.0 .or. m.eq.0) THEN
!
!       1. n = 0 and  m >= 0  OR n > 0 and m = 0
!
               tlpm = tlps + tlms
               slpm = slps + slms
               CALL analysum(grpmn, bvec, slpm, tlpm, m, n,
     &                       ivacskip, ndim)

            ELSE
!
!       2. n>=1  and  m>=1
!
               CALL analysum2(grpmn, bvec, slms, tlms, slps, tlps,
     &                        m, n, ivacskip, ndim)

            ENDIF
         END DO
      END DO

      DEALLOCATE (r0p, r1p, r0m, r1m, sqrtc, sqrta, tlp, tlm, adp,
     &            adm, cma, ra1p, ra1m, slm, slp, tlpm, slpm, tlps,
     &            tlms, slps, slms, ulp1, ulp2, ulm1, ulm2, stat = l)

      CALL second0(tanaloff)
      timer_vac(tanal) = timer_vac(tanal) + (tanaloff-tanalon)
      analyt_time = timer_vac(tanal)

      END SUBROUTINE analyt

!-------------------------------------------------------------------------------
!>  @brief Chebyshev moments of the tangent-plane kernel.
!>
!>  M_k = Int[-1,1] T_k(t)/SQRT(Q(t)) dt for k = 0, ..., mf + nf, where
!>  Q = A t^2 + 2 d t + B with A = a + b2 + c, B = a - b2 + c and d = c - a.
!>
!>  Integrating d/dt (F_k SQRT(Q)) with F_k' = T_k,
!>  F_k = T_(k+1)/(2 (k + 1)) - T_(k-1)/(2 (k - 1)), gives for k >= 2 the
!>  five-term recurrence whose factors @ref moment_factors tabulates. Its
!>  characteristic roots are those of Q((z + 1/z)/2) = 0: a complex pair of
!>  modulus rho > 1 and their reciprocals, so two homogeneous solutions grow
!>  and two decay in either direction and neither a forward nor a backward
!>  pass is stable. The moments decay like k^-2 only, and they are the solution
!>  of the boundary-value problem whose lower boundary values are M_0 and M_1
!>  and whose upper boundary values are the asymptotic
!>  M_k = -(1/SQRT(Q(1)) + (-1)^k/SQRT(Q(-1)))/(k^2 - 1) + O(k^-4),
!>  which a pentadiagonal elimination solves. The rows of its matrix tend to
!>  those of the Toeplitz matrix with the symbol Q(COS(theta)) > 0, and the
!>  elimination runs without pivoting. It is sequential in k and independent
!>  across grid points, so nbatch points are eliminated together, over the
!>  largest of their extents.
!>
!>  @param[in]  a  guu on the local grid points.
!>  @param[in]  b2 guv on the local grid points, with its factor of two. Its
!>                 negative gives the moments of the minus parity.
!>  @param[in]  c  gvv on the local grid points.
!>  @param[out] tl Chebyshev moments by grid point and order.
!-------------------------------------------------------------------------------
      SUBROUTINE chebyshev_moments(a, b2, c, tl)
      USE parallel_include_module
      USE vacmod0, ONLY: mf, nf

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(nuv3min:nuv3max), INTENT(in)         :: a
      REAL(dp), DIMENSION(nuv3min:nuv3max), INTENT(in)         :: b2
      REAL(dp), DIMENSION(nuv3min:nuv3max), INTENT(in)         :: c
      REAL(dp), DIMENSION(nuv3min:nuv3max,0:nf + mf), INTENT(out)
     &   :: tl

!  local parameters
!  Number of grid points eliminated together.
      INTEGER, PARAMETER  :: nbatch = 8
!  A boundary value contaminates the moments below it by rho^-(distance). The
!  problem extends ntail orders above mf + nf, so that this falls below 1e-17
!  at mf + nf. ntail stays within [kMinTail,kMaxTail]; the upper bound binds
!  for rho < 1.01, where the contamination is rho^-kMaxTail of the error of
!  the upper boundary values, which is O(k^-4) + O(rho^-k) at k > kMaxTail.
      REAL(dp), PARAMETER :: kMinBoundaryLogDecay = 39.14394658089878_dp
      INTEGER, PARAMETER  :: kMinTail = 8
      INTEGER, PARAMETER  :: kMaxTail = 4096

!  local variables
!  Rows of the eliminated systems of one batch, (point, k): the two bands
!  above the unit diagonal, and the right-hand side, which the back
!  substitution turns into the moments.
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: up1
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: up2
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: tk
      REAL(dp), DIMENSION(nbatch)           :: aa
      REAL(dp), DIMENSION(nbatch)           :: dd
      REAL(dp), DIMENSION(nbatch)           :: hab
      REAL(dp), DIMENSION(nbatch)           :: sqp
      REAL(dp), DIMENSION(nbatch)           :: sqm
      REAL(dp)                              :: sub2
      REAL(dp)                              :: sub1
      REAL(dp)                              :: rpiv
      REAL(dp)                              :: sign1
      REAL(dp)                              :: semi_major
      REAL(dp)                              :: semi_minor
      REAL(dp)                              :: log_rho
      INTEGER                               :: i0
      INTEGER                               :: i
      INTEGER                               :: j
      INTEGER                               :: k
      INTEGER                               :: kl
      INTEGER                               :: ktop
      INTEGER                               :: ntail
      INTEGER                               :: nb
      INTEGER                               :: istat

!  Start of executable code
      kl = mf + nf
      CALL moment_factors(kl + kMaxTail + 2)
      ALLOCATE (up1(nbatch,0:kl + kMaxTail + 2),
     &          up2(nbatch,0:kl + kMaxTail + 2),
     &          tk(nbatch,0:kl + kMaxTail + 2), stat = istat)
      IF (istat .ne. 0) THEN
         STOP 'Allocation error in SUBROUTINE chebyshev_moments'
      ENDIF

!  M_0 and M_1 enter as rows of the identity.
      up1(:,0:1) = 0
      up2(:,0:1) = 0

      DO i0 = nuv3min, nuv3max, nbatch
         nb = MIN(nbatch, nuv3max - i0 + 1)
         ntail = kMinTail
         DO j = 1, nbatch
!  Points past the last one repeat it.
            i = MIN(i0 + j - 1, nuv3max)
            aa(j) = a(i) + b2(i) + c(i)
            dd(j) = c(i) - a(i)
            hab(j) = 0.5_dp*aa(j) + (a(i) - b2(i) + c(i))
            sqp(j) = 2.0_dp*SQRT(c(i))
            sqm(j) = 2.0_dp*SQRT(a(i))

            tk(j,0) = t0_integral(a(i), b2(i), c(i))
!  From Int (A t + d)/SQRT(Q) dt = SQRT(Q(1)) - SQRT(Q(-1)).
            tk(j,1) = (sqp(j) - sqm(j) - dd(j)*tk(j,0))/aa(j)

!  The roots of Q lie on the ellipse with foci -1 and 1 whose semi-axes are
!  (rho + 1/rho)/2 = (SQRT(Q(1)) + SQRT(Q(-1)))/(2 SQRT(A)) and
!  (rho - 1/rho)/2 = SQRT((2 SQRT(a c) - b2)/A).
            semi_major = 0.5_dp*(sqp(j) + sqm(j))/SQRT(aa(j))
            semi_minor = SQRT(MAX(2.0_dp*SQRT(a(i)*c(i)) - b2(i),
     &                            0.0_dp)/aa(j))
            log_rho = LOG(semi_major + semi_minor)
            IF (log_rho .gt. 0.0_dp) THEN
               ntail = MAX(ntail,
     &                     CEILING(MIN(REAL(kMaxTail,dp),
     &                                 kMinBoundaryLogDecay/log_rho)))
            ELSE
               ntail = kMaxTail
            END IF
         END DO

!  Row k of the eliminated system reads
!  M_k + up1(k) M_(k+1) + up2(k) M_(k+2) = tk(k).
         ktop = kl + ntail
         sign1 = 1
         DO k = 2, ktop
            DO j = 1, nbatch
               sub2 = aa(j)*fsub2(k)
               sub1 = dd(j)*fsub1(k) - up1(j,k - 2)*sub2
               rpiv = 1.0_dp/(hab(j) - aa(j)*fdiag(k)
     &              -         up2(j,k - 2)*sub2 - up1(j,k - 1)*sub1)
               up1(j,k) = (dd(j)*fsup1(k) - up2(j,k - 1)*sub1)*rpiv
               up2(j,k) = aa(j)*fsup2(k)*rpiv
               tk(j,k) = (-(sqp(j) + sign1*sqm(j))*frhs(k)
     &                 -  tk(j,k - 2)*sub2 - tk(j,k - 1)*sub1)*rpiv
            END DO
            sign1 = -sign1
         END DO
         DO k = ktop + 1, ktop + 2
            tk(:,k) = -(1.0_dp/sqp(:) + sign1/sqm(:))*frhs(k)
            sign1 = -sign1
         END DO
         DO k = ktop, 2, -1
            tk(:,k) = tk(:,k) - up1(:,k)*tk(:,k + 1)
     &              -           up2(:,k)*tk(:,k + 2)
         END DO

         tl(i0:i0 + nb - 1,:) = tk(1:nb,0:kl)
      END DO

      DEALLOCATE (up1, up2, tk, stat = istat)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Zeroth moment of the tangent-plane kernel.
!>
!>  T_0 = Int[-1,1] dt/SQRT(A t^2 + 2 d t + B) for A = a + b2 + c,
!>  B = a - b2 + c and d = c - a:
!>  SQRT(A) T_0 = LOG((2 SQRT(c A) + 2 c + b2)/(2 SQRT(a A) - 2 a - b2)).
!>  A term of the quotient that is a difference of nearly equal numbers, as the
!>  denominator is for gvv << guu, is taken in its conjugate form, whose
!>  numerator 4 a c - b2^2 is the determinant of the metric.
!>
!>  @param[in] a  guu
!>  @param[in] b2 guv with its factor of two.
!>  @param[in] c  gvv
!>  @returns T_0
!-------------------------------------------------------------------------------
      PURE FUNCTION t0_integral(a, b2, c)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp)             :: t0_integral
      REAL(dp), INTENT(in) :: a
      REAL(dp), INTENT(in) :: b2
      REAL(dp), INTENT(in) :: c

!  local variables
      REAL(dp)             :: aa
      REAL(dp)             :: root
      REAL(dp)             :: det
      REAL(dp)             :: hi
      REAL(dp)             :: lo
      REAL(dp)             :: sqca
      REAL(dp)             :: sqaa
      REAL(dp)             :: num
      REAL(dp)             :: den

!  Start of executable code
      aa = a + b2 + c
      root = 2.0_dp*SQRT(a*c)
      det = (root - b2)*(root + b2)
      hi = 2.0_dp*c + b2
      lo = 2.0_dp*a + b2
      sqca = 2.0_dp*SQRT(c*aa)
      sqaa = 2.0_dp*SQRT(a*aa)
      IF (hi .ge. 0.0_dp) THEN
         num = sqca + hi
      ELSE
         num = det/(sqca - hi)
      END IF
      IF (lo .ge. 0.0_dp) THEN
         den = det/(sqaa + lo)
      ELSE
         den = sqaa - lo
      END IF
      t0_integral = LOG(num/den)/SQRT(aa)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Tabulate the factors of the five-term moment recurrence.
!>
!>  At order k >= 2 the recurrence reads
!>  A fsup2 M_(k+2) + d fsup1 M_(k+1) + (A/2 + B - A fdiag) M_k
!>  + d fsub1 M_(k-1) + A fsub2 M_(k-2)
!>  = -frhs (SQRT(Q(1)) + (-1)^k SQRT(Q(-1))).
!>  The tables are kept between calls and rebuilt when a larger order is
!>  requested.
!>
!>  @param[in] kmax Largest order needed.
!-------------------------------------------------------------------------------
      SUBROUTINE moment_factors(kmax)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in) :: kmax

!  local variables
      INTEGER             :: k
      INTEGER             :: istat

!  Start of executable code
      IF (ALLOCATED(fsub2)) THEN
         IF (UBOUND(fsub2,1) .ge. kmax) THEN
            RETURN
         END IF
         DEALLOCATE (fsub2, fsub1, fdiag, fsup1, fsup2, frhs)
      END IF

      ALLOCATE (fsub2(2:kmax), fsub1(2:kmax), fdiag(2:kmax),
     &          fsup1(2:kmax), fsup2(2:kmax), frhs(2:kmax),
     &          stat = istat)
      IF (istat .ne. 0) THEN
         STOP 'Allocation error in SUBROUTINE moment_factors'
      ENDIF

      DO k = 2, kmax
         fsub2(k) = (k - 2.0_dp)/(4.0_dp*(k - 1.0_dp))
         fsub1(k) = (2.0_dp*k - 3.0_dp)/(2.0_dp*(k - 1.0_dp))
         fdiag(k) = 1.0_dp/(2.0_dp*(k*k - 1.0_dp))
         fsup1(k) = (2.0_dp*k + 3.0_dp)/(2.0_dp*(k + 1.0_dp))
         fsup2(k) = (k + 2.0_dp)/(4.0_dp*(k + 1.0_dp))
         frhs(k) = 1.0_dp/(k*k - 1.0_dp)
      END DO

      END SUBROUTINE

      END MODULE
