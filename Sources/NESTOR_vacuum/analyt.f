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
      
      IMPLICIT NONE

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
      INTEGER :: l, n, m, k, mn
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: tlp, tlm
      REAL(dp), DIMENSION(:), ALLOCATABLE ::
     &   r0p, r1p, r0m, r1m, sqrtc, sqrta, adp, adm, cma, ra1p, ra1m,
     &   slm, slp, tlpm, slpm
      REAL(dp) :: sign1, tanalon, tanaloff
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
     &          tlp(-1:nf + mf,nuv3min:nuv3max),
     &          tlm(-1:nf + mf,nuv3min:nuv3max),
     &          slm(nuv3min:nuv3max), slp(nuv3min:nuv3max), stat = l)
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
      CALL initialize(adp, adm, cma, sqrtc, sqrta, tlp, tlm)
!
!     BEGIN L-SUM IN EQ (A14) TO COMPUTE Imn (and Kmn) INTEGRALS
!     NOTE THAT IN THE LOOP OVER L BELOW: L == |m - n| + 2L_A14
!     THUS, L BELOW IS THE INDEX OF THE T+- (S+-)
!
      sign1 = 1

      LLOOP: DO l = 0, mf + nf
!
!     COMPUTE SL+ and SL- , Eq (A17)
!     SLP(M): SL+(-)
!
         IF (ivacskip .eq. 0) THEN
            DO k = nuv3min, nuv3max
               slp(k) = (r1p(k)*l + ra1p(k))*tlp(l,k)
     &                + r0p(k)*l*tlp(l - 1,k)
     &                - (r1p(k) + r0p(k))/sqrtc(k)
     &                + sign1*(r0p(k) - r1p(k))/sqrta(k)
               slm(k) = (r1m(k)*l + ra1m(k))*tlm(l,k)
     &                + r0m(k)*l*tlm(l - 1,k)
     &                - (r1m(k) + r0m(k))/sqrtc(k)
     &                + sign1*(r0m(k) - r1m(k))/sqrta(k)
               slpm(k) = slp(k) + slm(k)
            END DO
         ENDIF
         tlpm = tlp(l,:) + tlm(l,:)
!
!     BEGIN MODE NUMBER (m,n) LOOP
!
         DO n = 0, nf
            DO m = 0, mf

               IF (l .EQ. 0) THEN
                  mn = m + mf1*(n+nf) + 1
                  bvec(mn,:) = 0
                  mn = m + mf1*(nf-n) + 1
                  bvec(mn,:) = 0
               END IF

               IF (cmns(l,m,n) .eq. zero) CYCLE
               
               IF (n.eq.0 .or. m.eq.0) THEN
!
!       1. n = 0 and  m >= 0  OR n > 0 and m = 0
!
                  CALL analysum(grpmn, bvec, slpm, tlpm, m, n, l,
     &                          ivacskip, ndim)

               ELSE
!
!       2. n>=1  and  m>=1
!
                  CALL analysum2(grpmn, bvec, slm, tlm(l,:), slp,
     &                           tlp(l,:), m, n, l, ivacskip, ndim)

               ENDIF
            END DO
         END DO
         sign1 = -sign1
      END DO LLOOP

      DEALLOCATE (r0p, r1p, r0m, r1m, sqrtc, sqrta, tlp, tlm, adp,
     &            adm, cma, ra1p, ra1m, slm, slp, tlpm, slpm, stat = l)

      CALL second0(tanaloff)
      timer_vac(tanal) = timer_vac(tanal) + (tanaloff-tanalon)
      analyt_time = timer_vac(tanal)

      END SUBROUTINE analyt

!-------------------------------------------------------------------------------
!>  @brief Initalize original
!>
!>  INITIALIZE T0+ and T0-
!>
!>  TLP(M): TL+(-)
!>  TLP(M)1:T(L-1)+(-)
!>  TLP(M)2:T(L-2)+(-)
!>
!>  @param[in]  adp
!>  @param[in]  adm
!>  @param[in]  cma
!>  @param[in]  sqrtc
!>  @param[in]  sqrta
!>  @param[out] tlp
!>  @param[out] tlm
!-------------------------------------------------------------------------------
      PURE SUBROUTINE initialize(adp, adm, cma, sqrtc, sqrta, tlp, tlm)
      USE parallel_include_module
      USE vacmod0, ONLY: mf, nf

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(nuv3min:nuv3max), INTENT(in)  :: adp
      REAL(dp), DIMENSION(nuv3min:nuv3max), INTENT(in)  :: adm
      REAL(dp), DIMENSION(nuv3min:nuv3max), INTENT(in)  :: cma
      REAL(dp), DIMENSION(nuv3min:nuv3max), INTENT(in)  :: sqrtc
      REAL(dp), DIMENSION(nuv3min:nuv3max), INTENT(in)  :: sqrta
      REAL(dp), DIMENSION(-1:nf + mf,nuv3min:nuv3max), INTENT(out)
     &   :: tlp
      REAL(dp), DIMENSION(-1:nf + mf,nuv3min:nuv3max), INTENT(out)
     &   :: tlm

!  local variables
      REAL(dp)                                          :: sqad1u
      REAL(dp)                                          :: sqad2u
      INTEGER                                           :: k
      REAL(dp)                                          :: t0p
      REAL(dp)                                          :: t0m
      REAL(dp)                                          :: high
      REAL(dp)                                          :: current
      REAL(dp)                                          :: low
      REAL(dp)                                          :: scale
      INTEGER                                           :: sign1

!  Start of executable code
      DO k = nuv3min, nuv3max
         sqad1u = SQRT(adp(k))
         sqad2u = SQRT(adm(k))
         tlp(-1,k) = 0
         tlm(-1,k) = 0
         tlp(0,k)  = log((sqad1u*sqrtc(k) + adp(k) + cma(k)) /
     &                   (sqad1u*sqrta(k) - adp(k) + cma(k)))/sqad1u
         tlm(0,k)  = log((sqad2u*sqrtc(k) + adm(k) + cma(k)) /
     &                   (sqad2u*sqrta(k) - adm(k) + cma(k)))/sqad2u

         CALL recurrence_sum(adm(k), adp(k), sqrtc(k), sqrta(k),
     &                       cma(k), tlp(0,k), tlp(:,k))
         CALL recurrence_sum(adp(k), adm(k), sqrtc(k), sqrta(k),
     &                       cma(k), tlm(0,k), tlm(:,k))
      END DO

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Querey if to use forward or backwards propagation.
!>
!>  Only switch to backward when the forward spurious-mode growth
!>  (|r1 r2| = B/A) would actually exceed double precision over kL steps.
!>  Threshold: forward is considered stable as long as (B/A)^kL < 1e10,
!>  i.e. spurious amplitude stays within ~1e10 of the particular solution.
!>  Near-degenerate kl (|r1|~|r2|~1) fall in the forward branch, where
!>  zero-seed Miller is known to misconverge (spurious modes never damp).
!>  Formula: kL * ln(B/A) < ln(1e10) -> B/A < exp(ln(1e10)/kL).
!>
!>  @param[in] a
!>  @param[in] b
!>  @returns true if log(a/b) > ln(1E10)
!-------------------------------------------------------------------------------
      PURE FUNCTION useBackward(a, b)
      USE stel_kinds
      USE vacmod0, ONLY: mf, nf

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL              :: useBackward
      REAL(dp), INTENT(in) :: a
      REAL(dp), INTENT(in) :: b

!  local parameters
      REAL(dp), PARAMETER  :: kLogGrowthThreshold = 23.0258509299 ! ln(1e10)

!  Start of executable code
      useBackward = a .gt. b   .and.
     &              b .gt. 0.0 .and.
     &              (mf + nf)*log(a/b) .gt. kLogGrowthThreshold

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Recurrence method.
!>
!>  @param[in] a       ad for m or p depending on what parity is computed
!>  @param[in] b       Opposite ad parity to a.
!>  @param[in] sqrtc   Sqrt(c) constant.
!>  @param[in] sqrta   Sqrt(a) constant.
!>  @param[in] sign1   Sign of the parity.
!>  @param[in] fl      Recurrence index.
!>  @param[in] cma
!>  @param[in] current Current value.
!>  @param[in] next    Next value.
!-------------------------------------------------------------------------------
      PURE FUNCTION recurrence(a, b, sqrtc, sqrta, sign1, fl, fl1, fl2,
     &                         cma, current, next)
      USE stel_kinds

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp)             :: recurrence
      REAL(dp), INTENT(in) :: a
      REAL(dp), INTENT(in) :: b
      REAL(dp), INTENT(in) :: sqrtc
      REAL(dp), INTENT(in) :: sqrta
      REAL(dp), INTENT(in) :: sign1
      INTEGER, INTENT(in)  :: fl
      INTEGER, INTENT(in)  :: fl1
      INTEGER, INTENT(in)  :: fl2
      REAL(dp), INTENT(in) :: cma
      REAL(dp), INTENT(in) :: current
      REAL(dp), INTENT(in) :: next

!  Start of executable code
      recurrence = (sqrtc + sign1*sqrta - fl2*cma*next - fl*a*current)
     &           / (b*fl1)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Built the recurrence for a parity.
!>
!>  @param[in]    a     ad for m or p depending on what parity is computed
!>  @param[in]    b     Opposite ad parity to a.
!>  @param[in]    sqrtc Sqrt(c) constant.
!>  @param[in]    sqrta Sqrt(a) constant.
!>  @param[in]    sign1 Sign of the parity.
!>  @param[in]    cma
!>  @param[in]    t0    Inital
!>  @param[inout] tl    Final recurrence.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE recurrence_sum(a, b, sqrtc, sqrta, cma, t0, tl)
      USE stel_kinds
      USE vacmod0, ONLY: mf, nf

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), INTENT(in)                           :: a
      REAL(dp), INTENT(in)                           :: b
      REAL(dp), INTENT(in)                           :: sqrtc
      REAL(dp), INTENT(in)                           :: sqrta
      REAL(dp), INTENT(in)                           :: cma
      REAL(dp), INTENT(in)                           :: t0
      REAL(dp), DIMENSION(-1:nf + mf), INTENT(inout) :: tl

!  local variables
      REAL(dp)                                       :: high
      REAL(dp)                                       :: current
      REAL(dp)                                       :: low
      REAL(dp)                                       :: scale
      REAL(dp)                                       :: sign1
      INTEGER                                        :: l

!  local parameters
      INTEGER, PARAMETER :: kTailExtra = 50

!  Start of executable code
      IF (useBackward(a,b)) THEN
         high = 0.0
         current = 1.0E-300_dp
         IF (MOD(mf + nf + kTailExtra, 2) .eq. 0) THEN
            sign1 = -1
         ELSE
            sign1 = 1
         ENDIF
         DO l = mf + nf + kTailExtra, mf + nf + 2, -1
            low = recurrence(b, a, sqrtc, sqrta, sign1,
     &                       l + 1, l, 2*l + 1, cma, high, current)
            high = current
            current = low
            sign1 = -sign1
         END DO
         DO l = mf + nf + 1, 1, -1
            tl(l - 1) = recurrence(b, a, sqrtc, sqrta, sign1,
     &                       l + 1, l, 2*l + 1, cma, high, current)
            high = current
            current = tl(l - 1)
            sign1 = -sign1
         END DO
         scale = t0/tl(0)
         DO l = 0, mf + nf
            tl(l) = tl(l)*scale
         END DO
      ELSE
         sign1 = 1
         DO l = 0, mf + nf - 1
            sign1 = -sign1
            tl(l + 1) = recurrence(a, b, sqrtc, sqrta, sign1,
     &                             l, l + 1, 2*l + 1, cma,
     &                             tl(l - 1), tl(l))
         END DO
      END IF
      END SUBROUTINE

      END MODULE
