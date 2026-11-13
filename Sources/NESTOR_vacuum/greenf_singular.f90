!*******************************************************************************
!>  @file greenf_singular.f90
!>  @brief Contains module @ref greenf_singular.
!
!  Note separating the Doxygen comment block here so detailed description is
!  found in the Module not the file.
!
!>  Defines the correction of the punctured trapezoidal sums of greenf at the
!>  source point.
!>
!>  About the source point the remainders G - ga1 and dG/dn - ga1*ga2 expand as
!>  f0 + f1 + ... in the steps (s, t) of u and v, with f0 odd and homogeneous
!>  of degree 0 and f1 even and homogeneous of degree 1. The punctured
!>  trapezoidal sum of a remainder times a smooth function b differs from the
!>  integral by L[f1] b + L[f0 s] db/du + L[f0 t] db/dv at third order in the
!>  grid step, where L is the regularized lattice sum. In coordinates in which
!>  the metric is the identity the six functions are |x| times even harmonics
!>  of the polar angle up to order 8, so L follows from Epstein zeta functions
!>  with harmonic weights, which are summed by Ewald's method.
!>
!>  The defects are the leading terms of an expansion in the lattice step over
!>  the radius of curvature of the boundary, so they are faded out where the
!>  grid does not resolve that radius: @ref defect_fade sets the fraction of
!>  them that is applied, one fraction for every source point of the surface.
!*******************************************************************************
      MODULE greenf_singular
      USE stel_kinds, ONLY: dp
      USE stel_constants, ONLY: pi

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: singular_defects
      PUBLIC :: defect_fade

!>  Number of polar angles sampled over a half turn. The six functions are
!>  trigonometric polynomials of degree 8 with period pi in the polar angle.
      INTEGER, PARAMETER  :: nang = 10
!>  Lattice step over the radius of curvature up to which the defects are
!>  applied in full.
      REAL(dp), PARAMETER :: span_full = 0.5_dp
!>  Lattice step over the radius of curvature from which they are dropped.
      REAL(dp), PARAMETER :: span_none = 2.0_dp
!>  Fraction of the defects applied on the current surface.
      REAL(dp)            :: fade = 1.0_dp

      CONTAINS

!-------------------------------------------------------------------------------
!>  @brief Product of two polynomials in (s, t), truncated at degree 4.
!>
!>  Coefficient (i,j) of a polynomial multiplies s**i*t**j.
!>
!>  @param[in] a First polynomial.
!>  @param[in] b Second polynomial.
!>  @returns The product without its terms above degree 4.
!-------------------------------------------------------------------------------
      PURE FUNCTION tmul(a, b)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(0:4,0:4)             :: tmul
      REAL(dp), DIMENSION(0:4,0:4), INTENT(IN) :: a
      REAL(dp), DIMENSION(0:4,0:4), INTENT(IN) :: b

!  Local Variables
      REAL(dp)                                 :: temp
      INTEGER                                  :: i
      INTEGER                                  :: j
      INTEGER                                  :: k

!  Start of executable code
!  Coefficient (k,j - k) of b times column i - j of a goes to column i - k of
!  the product, k rows down.
      tmul = 0.0_dp
      DO k = 0, 4
         DO j = k, 4
            temp = b(k,j - k)
            DO i = j, 4
               tmul(k:4 - (i - k),i - k) = tmul(k:4 - (i - k),i - k)           &
                                         + a(:4 - i,i - j)*temp
            END DO
         END DO
      END DO

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Homogeneous part of a polynomial in (s, t) at a point.
!>
!>  @param[in] a  Polynomial, whose coefficient (i,j) multiplies s**i*t**j.
!>  @param[in] d  Degree of the part.
!>  @param[in] sp Powers 0 to 4 of s at the point.
!>  @param[in] tq Powers 0 to 4 of t at the point.
!>  @returns The sum of the terms of degree d.
!-------------------------------------------------------------------------------
      PURE FUNCTION tpart(a, d, sp, tq)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp)                                 :: tpart
      REAL(dp), DIMENSION(0:4,0:4), INTENT(IN) :: a
      INTEGER, INTENT(IN)                      :: d
      REAL(dp), DIMENSION(0:4), INTENT(IN)     :: sp
      REAL(dp), DIMENSION(0:4), INTENT(IN)     :: tq

!  Local Variables
      INTEGER                                  :: i

!  Start of executable code
      tpart = 0.0_dp
      DO i = 0, d
         tpart = tpart + a(i,d - i)*sp(i)*tq(d - i)
      END DO

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Lagrange-Gauss reduction of a lattice basis.
!>
!>  @param[inout] f1 On exit the shortest vector of the lattice.
!>  @param[inout] f2 On exit the shortest vector that completes the basis.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE reduce_basis(f1, f2)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(2), INTENT(INOUT) :: f1
      REAL(dp), DIMENSION(2), INTENT(INOUT) :: f2

!  Local Variables
      REAL(dp), DIMENSION(2)                :: w
      REAL(dp)                              :: mu
      INTEGER                               :: it
      REAL(dp)                              :: fdotf

!  Start of executable code
      DO it = 1, 100
         fdotf = DOT_PRODUCT(f1, f1)
         IF (fdotf .GT. DOT_PRODUCT(f2, f2)) THEN
            w = f1
            f1 = f2
            f2 = w
            CYCLE
         END IF
         mu = ANINT(DOT_PRODUCT(f1, f2)/fdotf)
         IF (mu .EQ. 0) THEN
            EXIT
         END IF
         f2 = f2 - mu*f1
      END DO

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Coordinates in which the metric of the surface is the identity.
!>
!>  For the steps (s, t) of u and v the coordinates are (ca*s + cb*t, cc*t).
!>
!>  @param[in]  g  Coefficients of the u, mixed and v terms in the quadratic
!>                 form of ga1.
!>  @param[out] ca Coefficient of s in the first coordinate.
!>  @param[out] cb Coefficient of t in the first coordinate.
!>  @param[out] cc Coefficient of t in the second coordinate.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE metric_frame(g, ca, cb, cc)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(3), INTENT(IN) :: g
      REAL(dp), INTENT(OUT)              :: ca
      REAL(dp), INTENT(OUT)              :: cb
      REAL(dp), INTENT(OUT)              :: cc

!  Start of executable code
      ca = SQRT(g(1))
      cb = g(2)/(2*ca)
      cc = SQRT(g(3) - cb*cb)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Upper incomplete gamma functions of half-integer order.
!>
!>  GAMMA(1/2, x) is SQRT(pi)*ERFC(SQRT(x)), and the other orders follow from
!>  GAMMA(a + 1, x) = a*GAMMA(a, x) + x**a*EXP(-x).
!>
!>  @param[in] x Argument, which is positive.
!>  @returns GAMMA(m - 1/2, x) for m = 0, ..., 6.
!-------------------------------------------------------------------------------
      PURE FUNCTION incomplete_gammas(x)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(0:6) :: incomplete_gammas
      REAL(dp), INTENT(IN)     :: x

!  Local Variables
      REAL(dp)                 :: ex
      REAL(dp)                 :: sx
      REAL(dp)                 :: a
      REAL(dp)                 :: xa
      INTEGER                  :: m

!  Start of executable code
      ex = EXP(-x)
      sx = SQRT(x)
      incomplete_gammas(1) = SQRT(pi)*ERFC(sx)
      incomplete_gammas(0) = 2*(ex/sx - incomplete_gammas(1))
      a = 0.5_dp
      xa = sx
      DO m = 2, 6
         incomplete_gammas(m) = a*incomplete_gammas(m - 1) + xa*ex
         a = a + 1
         xa = xa*x
      END DO

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Ewald sums of the harmonic kernels over half of a lattice.
!>
!>  The term of a point x of the lattice is
!>  (x1 + i*x2)**k*X**(1/2 - q - shift)*GAMMA(q + shift - 1/2, X), with
!>  X = pi*|x|**2 and k = 2*q. The points are those of i*b1 + j*b2 with i > 0,
!>  or with i = 0 and j > 0, within the reach of the terms, which fall off as
!>  EXP(-X); x and -x contribute alike.
!>
!>  @param[in] b1    First basis vector of the lattice.
!>  @param[in] b2    Second basis vector of the lattice.
!>  @param[in] shift 0 for the lattice and 2 for its dual.
!>  @returns The sums for q = 0, ..., 4.
!-------------------------------------------------------------------------------
      PURE FUNCTION half_lattice_sums(b1, b2, shift)

      IMPLICIT NONE

!  Declare Arguments
      COMPLEX(dp), DIMENSION(0:4)        :: half_lattice_sums
      REAL(dp), DIMENSION(2), INTENT(IN) :: b1
      REAL(dp), DIMENSION(2), INTENT(IN) :: b2
      INTEGER, INTENT(IN)                :: shift

!  Local Variables
      REAL(dp), DIMENSION(2)             :: w
      REAL(dp), DIMENSION(0:6)           :: g
      REAL(dp)                           :: g11
      REAL(dp)                           :: g22
      REAL(dp)                           :: g12
      REAL(dp)                           :: gdet
      REAL(dp)                           :: r2
      REAL(dp)                           :: x
      REAL(dp)                           :: xinv
      REAL(dp)                           :: xp
      COMPLEX(dp)                        :: z2
      COMPLEX(dp)                        :: zk
      INTEGER                            :: i
      INTEGER                            :: j
      INTEGER                            :: ni
      INTEGER                            :: nj
      INTEGER                            :: q

!  Local Parameters
!  Square of the distance from the origin past which the terms are dropped.
      REAL(dp), PARAMETER                :: reach2 = 7.0_dp

!  Start of executable code
      g11 = DOT_PRODUCT(b1, b1)
      g22 = DOT_PRODUCT(b2, b2)
      g12 = DOT_PRODUCT(b1, b2)
      gdet = g11*g22 - g12*g12
      ni = CEILING(SQRT(reach2*g22/gdet))
      nj = CEILING(SQRT(reach2*g11/gdet))

      half_lattice_sums = 0.0_dp
      DO i = 0, ni
         DO j = -nj, nj
            IF (i .EQ. 0 .AND. j .LE. 0) THEN
               CYCLE
            END IF
            w = i*b1 + j*b2
            r2 = DOT_PRODUCT(w, w)
            IF (r2 .GT. reach2) THEN
               CYCLE
            END IF
            x = pi*r2
            xinv = 1/x
            g = incomplete_gammas(x)
            z2 = CMPLX(w(1), w(2), dp)
            z2 = z2*z2
            zk = 1.0_dp
            xp = SQRT(x)
            DO q = 1, shift
               xp = xp*xinv
            END DO
            DO q = 0, 4
               half_lattice_sums(q) = half_lattice_sums(q)                     &
                                    + zk*xp*g(q + shift)
               zk = zk*z2
               xp = xp*xinv
            END DO
         END DO
      END DO

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Regularized lattice sums of the harmonic kernels.
!>
!>  The sum for k = 2*q is A*Z_k: the regularized sum of
!>  A*|x|**(1-k)*(x1 + i*x2)**k over the lattice spanned by e1 and e2, of cell
!>  area A, less its integral. Z_k is the Epstein zeta function with the
!>  harmonic weight at s = (k - 1)/2, split by Ewald's method between the
!>  lattice of unit cell area and its dual.
!>
!>  @param[in] e1 First basis vector of the lattice.
!>  @param[in] e2 Second basis vector of the lattice.
!>  @returns The sums for k = 0, 2, 4, 6 and 8.
!-------------------------------------------------------------------------------
      PURE FUNCTION harmonic_sums(e1, e2)

      IMPLICIT NONE

!  Declare Arguments
      COMPLEX(dp), DIMENSION(0:4)        :: harmonic_sums
      REAL(dp), DIMENSION(2), INTENT(IN) :: e1
      REAL(dp), DIMENSION(2), INTENT(IN) :: e2

!  Local Variables
      REAL(dp), DIMENSION(2)             :: f1
      REAL(dp), DIMENSION(2)             :: f2
      REAL(dp), DIMENSION(2)             :: d1
      REAL(dp), DIMENSION(2)             :: d2
      REAL(dp), DIMENSION(0:4)           :: pgam
      REAL(dp)                           :: scale
      REAL(dp)                           :: scale3
      REAL(dp)                           :: det
      REAL(dp)                           :: area
      COMPLEX(dp), DIMENSION(0:4)        :: direct
      COMPLEX(dp), DIMENSION(0:4)        :: dual
      COMPLEX(dp)                        :: total
      INTEGER                            :: q

!  Local Parameters
      REAL(dp), PARAMETER                :: inv2pi = 1.0_dp/(2.0_dp*pi)

!  Start of executable code
!  Unit cell area, then the Lagrange-Gauss reduced basis and its dual
      scale = SQRT(ABS(e1(1)*e2(2) - e1(2)*e2(1)))
      f1 = e1/scale
      f2 = e2/scale
      CALL reduce_basis(f1, f2)
      det = f1(1)*f2(2) - f1(2)*f2(1)
      d1 = (/ f2(2), -f2(1) /)/det
      d2 = (/ -f1(2), f1(1) /)/det
      area = ABS(det)

      direct = half_lattice_sums(f1, f2, 0)
      dual = half_lattice_sums(d1, d2, 2)

!  pi**s/GAMMA(s) at s = -1/2, 1/2, ..., 7/2
      pgam(0) = -inv2pi
      DO q = 1, 4
         pgam(q) = pi/(q - 1.5_dp)*pgam(q - 1)
      END DO

!  (-i)**k = (-1)**q; the sums ran over half of each lattice
      scale3 = scale*scale*scale
      total = 2*(direct(0) + dual(0)/area) - 2/(3*area) + 2
      harmonic_sums(0) = area*pgam(0)*total*scale3
      DO q = 1, 4
         total = 2*(direct(q) + (-1)**q*dual(q)/area)
         harmonic_sums(q) = area*pgam(q)*total*scale3
      END DO

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Expansions about the source point of the squared distance from it
!>  and of the displacement from it along its normal.
!>
!>  @param[in]  rt   Taylor coefficients of R in (s, t).
!>  @param[in]  zt   Taylor coefficients of Z in (s, t).
!>  @param[in]  x0   Cartesian position of the source.
!>  @param[in]  phi0 Toroidal angle of the source.
!>  @param[in]  dphi Rate of the toroidal angle in v.
!>  @param[in]  xn   Cartesian normal of the source.
!>  @param[out] r2   Squared distance from the source as a polynomial in (s, t).
!>  @param[out] nr   Displacement from the source along its normal likewise.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE source_series(rt, zt, x0, phi0, dphi, xn, r2, nr)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(0:4,0:4), INTENT(IN)  :: rt
      REAL(dp), DIMENSION(0:4,0:4), INTENT(IN)  :: zt
      REAL(dp), DIMENSION(3), INTENT(IN)        :: x0
      REAL(dp), INTENT(IN)                      :: phi0
      REAL(dp), INTENT(IN)                      :: dphi
      REAL(dp), DIMENSION(3), INTENT(IN)        :: xn
      REAL(dp), DIMENSION(0:4,0:4), INTENT(OUT) :: r2
      REAL(dp), DIMENSION(0:4,0:4), INTENT(OUT) :: nr

!  Local Variables
      REAL(dp), DIMENSION(0:4,0:4)              :: cphi
      REAL(dp), DIMENSION(0:4,0:4)              :: sphi
      REAL(dp), DIMENSION(0:4,0:4)              :: xc
      REAL(dp), DIMENSION(0:4,0:4)              :: yc
      REAL(dp), DIMENSION(0:4,0:4)              :: wc
      REAL(dp), DIMENSION(0:3)                  :: cq
      REAL(dp)                                  :: rate
      INTEGER                                   :: j

!  Local Parameters
      REAL(dp), DIMENSION(0:4), PARAMETER       :: fact = (/ 1, 1, 2, 6, 24 /)

!  Start of executable code
!  COS and SIN of the toroidal angle as polynomials in t, from COS and SIN of
!  phi0 + j*pi/2
      cq = (/ COS(phi0), -SIN(phi0), -COS(phi0), SIN(phi0) /)
      cphi = 0.0_dp
      sphi = 0.0_dp
      rate = 1.0_dp
      DO j = 0, 4
         cphi(0,j) = rate/fact(j)*cq(MOD(j, 4))
         sphi(0,j) = rate/fact(j)*cq(MOD(j + 3, 4))
         rate = rate*dphi
      END DO

!  Cartesian position relative to the source as polynomials in (s, t)
      xc = tmul(rt, cphi)
      yc = tmul(rt, sphi)
      wc = zt
      xc(0,0) = xc(0,0) - x0(1)
      yc(0,0) = yc(0,0) - x0(2)
      wc(0,0) = wc(0,0) - x0(3)

      r2 = tmul(xc, xc) + tmul(yc, yc) + tmul(wc, wc)
      nr = xn(1)*xc + xn(2)*yc + xn(3)*wc

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Leading terms of the remainders in a direction from the source.
!>
!>  @param[in] r2 Squared distance from the source as a polynomial in (s, t).
!>  @param[in] nr Displacement from the source along its normal likewise.
!>  @param[in] g  Coefficients of the u, mixed and v terms in the quadratic
!>                form of ga1.
!>  @param[in] a  Coefficients of the u, mixed and v terms in the quadratic
!>                form of ga2.
!>  @param[in] s  Step of u in the direction.
!>  @param[in] t  Step of v in the direction.
!>  @returns f1, f0*s and f0*t of G - ga1, then those of dG/dn - ga1*ga2.
!-------------------------------------------------------------------------------
      PURE FUNCTION remainder_terms(r2, nr, g, a, s, t)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(6)                   :: remainder_terms
      REAL(dp), DIMENSION(0:4,0:4), INTENT(IN) :: r2
      REAL(dp), DIMENSION(0:4,0:4), INTENT(IN) :: nr
      REAL(dp), DIMENSION(3), INTENT(IN)       :: g
      REAL(dp), DIMENSION(3), INTENT(IN)       :: a
      REAL(dp), INTENT(IN)                     :: s
      REAL(dp), INTENT(IN)                     :: t

!  Local Variables
      REAL(dp), DIMENSION(0:4)                 :: sp
      REAL(dp), DIMENSION(0:4)                 :: tq
      REAL(dp)                                 :: ri
      REAL(dp)                                 :: c3
      REAL(dp)                                 :: d4
      REAL(dp)                                 :: a2
      REAL(dp)                                 :: a3
      REAL(dp)                                 :: a4
      REAL(dp)                                 :: s3t
      REAL(dp)                                 :: t4
      REAL(dp)                                 :: s4
      REAL(dp)                                 :: i32
      REAL(dp)                                 :: i52
      REAL(dp)                                 :: i72
      REAL(dp)                                 :: f0g
      REAL(dp)                                 :: f1g
      REAL(dp)                                 :: f0k
      REAL(dp)                                 :: f1k
      INTEGER                                  :: m

!  Start of executable code
      sp(0) = 1.0_dp
      tq(0) = 1.0_dp
      DO m = 1, 4
         sp(m) = sp(m - 1)*s
         tq(m) = tq(m - 1)*t
      END DO

!  Parts of degree 2, 3 and 4 of the squared distance and of the displacement
!  along the normal
      ri = tpart(r2, 2, sp, tq)
      c3 = tpart(r2, 3, sp, tq)
      d4 = tpart(r2, 4, sp, tq)
      a2 = tpart(nr, 2, sp, tq)
      a3 = tpart(nr, 3, sp, tq)
      a4 = tpart(nr, 4, sp, tq)

!  Quartic parts of the quadratic forms of ga1 and ga2 in 2*TAN(s/2) and
!  2*TAN(t/2)
      s3t = sp(3)*t + s*tq(3)
      t4 = (g(1)*sp(4) + g(3)*tq(4))/6 + g(2)*s3t/12
      s4 = (a(1)*sp(4) + a(3)*tq(4))/6 + a(2)*s3t/12

      i32 = 1/(ri*SQRT(ri))
      i52 = i32/ri
      i72 = i52/ri
      f0g = -c3*i32/2
      f1g = -(d4 - t4)*i32/2 + 3*c3*c3*i52/8
      f0k = a3*i32 - 3*a2*c3*i52/2
      f1k = (a4 - s4)*i32 - 3*(a3*c3 + a2*(d4 - t4))*i52/2                     &
          + 15*a2*c3*c3*i72/8
      remainder_terms = (/ f1g, f0g*s, f0g*t, f1k, f0k*s, f0k*t /)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Harmonics of the leading terms of the remainders on the unit circle
!>  of the coordinates of @ref metric_frame.
!>
!>  @param[in] r2 Squared distance from the source as a polynomial in (s, t).
!>  @param[in] nr Displacement from the source along its normal likewise.
!>  @param[in] g  Coefficients of the u, mixed and v terms in the quadratic
!>                form of ga1.
!>  @param[in] a  Coefficients of the u, mixed and v terms in the quadratic
!>                form of ga2.
!>  @returns The harmonics of order 2*k, k = 0, ..., 4, of the six functions of
!>           @ref remainder_terms.
!-------------------------------------------------------------------------------
      PURE FUNCTION remainder_harmonics(r2, nr, g, a)

      IMPLICIT NONE

!  Declare Arguments
      COMPLEX(dp), DIMENSION(6,0:4)            :: remainder_harmonics
      REAL(dp), DIMENSION(0:4,0:4), INTENT(IN) :: r2
      REAL(dp), DIMENSION(0:4,0:4), INTENT(IN) :: nr
      REAL(dp), DIMENSION(3), INTENT(IN)       :: g
      REAL(dp), DIMENSION(3), INTENT(IN)       :: a

!  Local Variables
      REAL(dp), DIMENSION(6)                   :: phi
      REAL(dp)                                 :: ca
      REAL(dp)                                 :: cb
      REAL(dp)                                 :: cc
      REAL(dp)                                 :: th
      REAL(dp)                                 :: s
      REAL(dp)                                 :: t
      COMPLEX(dp)                              :: e2th
      COMPLEX(dp)                              :: ek
      INTEGER                                  :: n
      INTEGER                                  :: k

!  Start of executable code
      CALL metric_frame(g, ca, cb, cc)

      remainder_harmonics = 0.0_dp
      DO n = 0, nang - 1
!  Steps of u and v at the polar angle th of those coordinates
         th = pi*n/nang
         t = SIN(th)/cc
         s = (COS(th) - cb*t)/ca
         phi = remainder_terms(r2, nr, g, a, s, t)/nang
         e2th = CMPLX(COS(2*th), -SIN(2*th), dp)
         ek = 1.0_dp
         DO k = 0, 4
            remainder_harmonics(:,k) = remainder_harmonics(:,k) + phi*ek
            ek = ek*e2th
         END DO
      END DO

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Defects of the punctured trapezoidal sums at the source point.
!>
!>  The defects are minus the regularized lattice sums, in the units of the
!>  integral over (s, t), the steps of u and v from the source point, times
!>  the fraction set by @ref defect_fade.
!>
!>  @param[in]  rt   Taylor coefficients of R in (s, t).
!>  @param[in]  zt   Taylor coefficients of Z in (s, t).
!>  @param[in]  x0   Cartesian position of the source.
!>  @param[in]  phi0 Toroidal angle of the source.
!>  @param[in]  dphi Rate of the toroidal angle in v.
!>  @param[in]  xn   Cartesian normal of the source.
!>  @param[in]  g    Coefficients of the u, mixed and v terms in the quadratic
!>                   form of ga1.
!>  @param[in]  a    Coefficients of the u, mixed and v terms in the quadratic
!>                   form of ga2.
!>  @param[in]  h    Grid steps in u and v.
!>  @param[out] wg   Defects of G - ga1 for the weights 1, s and t.
!>  @param[out] wk   Defects of dG/dn - ga1*ga2 for the weights 1, s and t.
!-------------------------------------------------------------------------------
      SUBROUTINE singular_defects(rt, zt, x0, phi0, dphi, xn, g, a, h, wg, wk)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(0:4,0:4), INTENT(IN) :: rt
      REAL(dp), DIMENSION(0:4,0:4), INTENT(IN) :: zt
      REAL(dp), DIMENSION(3), INTENT(IN)       :: x0
      REAL(dp), INTENT(IN)                     :: phi0
      REAL(dp), INTENT(IN)                     :: dphi
      REAL(dp), DIMENSION(3), INTENT(IN)       :: xn
      REAL(dp), DIMENSION(3), INTENT(IN)       :: g
      REAL(dp), DIMENSION(3), INTENT(IN)       :: a
      REAL(dp), DIMENSION(2), INTENT(IN)       :: h
      REAL(dp), DIMENSION(3), INTENT(OUT)      :: wg
      REAL(dp), DIMENSION(3), INTENT(OUT)      :: wk

!  Local Variables
      REAL(dp), DIMENSION(0:4,0:4)             :: r2
      REAL(dp), DIMENSION(0:4,0:4)             :: nr
      REAL(dp), DIMENSION(6)                   :: w
      REAL(dp)                                 :: ca
      REAL(dp)                                 :: cb
      REAL(dp)                                 :: cc
      REAL(dp)                                 :: el
      COMPLEX(dp), DIMENSION(6,0:4)            :: coef
      COMPLEX(dp), DIMENSION(0:4)              :: esum
      INTEGER                                  :: k
      INTEGER                                  :: m

!  Start of executable code
      CALL source_series(rt, zt, x0, phi0, dphi, xn, r2, nr)
      coef = remainder_harmonics(r2, nr, g, a)

!  Lattice of the grid in the coordinates of metric_frame
      CALL metric_frame(g, ca, cb, cc)
      esum = harmonic_sums((/ ca*h(1), 0.0_dp /), (/ cb*h(2), cc*h(2) /))

      DO m = 1, 6
         el = REAL(coef(m,0)*esum(0), dp)
         DO k = 1, 4
            el = el + 2*REAL(coef(m,k)*esum(k), dp)
         END DO
         w(m) = -el/(ca*cc)
      END DO
      wg = fade*w(1:3)
      wk = fade*w(4:6)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Set the fraction of the defects applied on the current surface.
!>
!>  At each point of the surface the longer vector of the reduced basis of the
!>  grid lattice is measured against the smaller radius of curvature. The
!>  defects are applied in full while the largest ratio on the surface is
!>  below span_full and dropped above span_none, with a smooth step between.
!>  For nv = 1 the v step is the image spacing, as in greenf.
!-------------------------------------------------------------------------------
      SUBROUTINE defect_fade
      USE vacmod
      USE parallel_include_module

      IMPLICIT NONE

!  Local Variables
      REAL(dp), DIMENSION(2) :: f1
      REAL(dp), DIMENSION(2) :: f2
      REAL(dp)               :: hv
      REAL(dp)               :: ca
      REAL(dp)               :: cb
      REAL(dp)               :: cc
      REAL(dp)               :: b11
      REAL(dp)               :: b12
      REAL(dp)               :: b22
      REAL(dp)               :: kmax
      REAL(dp)               :: ratio
      REAL(dp)               :: x
      INTEGER                :: i

!  Start of executable code
      hv = alv
      IF (nv .EQ. 1) THEN
         hv = pi2/nvper
      END IF

      ratio = 0.0_dp
      DO i = nuv3min, nuv3max
         CALL metric_frame((/ guu_b(i), guv_b(i), gvv_b(i) /), ca, cb, cc)
         f1 = (/ ca*alu, 0.0_dp /)
         f2 = (/ cb*hv, cc*hv /)
         CALL reduce_basis(f1, f2)

!  Second fundamental form in the coordinates of metric_frame, times half the
!  normal
         b11 = auu(i)/(ca*ca)
         b12 = (auv(i)/2 - auu(i)*cb/ca)/(ca*cc)
         b22 = (avv(i) - auv(i)*cb/ca + auu(i)*cb*cb/(ca*ca))/(cc*cc)
         kmax = (ABS(b11 + b22) +                                              &
                 SQRT((b11 - b22)*(b11 - b22) + 4*b12*b12))                    &
              / SQRT(snr(i)*snr(i) + snv(i)*snv(i) + snz(i)*snz(i))
         ratio = MAX(ratio, SQRT(DOT_PRODUCT(f2, f2))*kmax)
      END DO
      IF (vlactive) THEN
         CALL MPI_Allreduce(MPI_IN_PLACE, ratio, 1, MPI_REAL8, MPI_MAX,        &
                            VAC_COMM, MPI_ERR)
      END IF

      x = MIN(1.0_dp, MAX(0.0_dp, (span_none - ratio)/                         &
                                  (span_none - span_full)))
      fade = x*x*(3 - 2*x)

      END SUBROUTINE

      END MODULE
