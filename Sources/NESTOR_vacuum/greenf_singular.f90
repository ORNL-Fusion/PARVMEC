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
!*******************************************************************************
      MODULE greenf_singular
      USE stel_kinds, ONLY: dp
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: singular_defects

!>  Pi.
      REAL(dp), PARAMETER :: pi = 3.14159265358979323846264338_dp
!>  Number of polar angles sampled over a half turn. The six functions are
!>  trigonometric polynomials of degree 8 with period pi in the polar angle.
      INTEGER, PARAMETER :: nang = 10

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
      PURE FUNCTION tmul(a, b) RESULT(c)
      REAL(dp), INTENT(IN) :: a(0:4,0:4), b(0:4,0:4)
      REAL(dp) :: c(0:4,0:4)
      INTEGER :: i, j, k, l
      c = 0
      DO i = 0, 4
         DO j = 0, 4 - i
            IF (a(i,j) .EQ. 0) THEN
               CYCLE
            END IF
            DO k = 0, 4 - i - j
               DO l = 0, 4 - i - j - k
                  c(i+k,j+l) = c(i+k,j+l) + a(i,j)*b(k,l)
               END DO
            END DO
         END DO
      END DO
      END FUNCTION tmul

!-------------------------------------------------------------------------------
!>  @brief Homogeneous part of a polynomial in (s, t) at a point.
!>
!>  @param[in] a  Polynomial, whose coefficient (i,j) multiplies s**i*t**j.
!>  @param[in] d  Degree of the part.
!>  @param[in] sp Powers 0 to 4 of s at the point.
!>  @param[in] tq Powers 0 to 4 of t at the point.
!>  @returns The sum of the terms of degree d.
!-------------------------------------------------------------------------------
      PURE REAL(dp) FUNCTION tpart(a, d, sp, tq)
      REAL(dp), INTENT(IN) :: a(0:4,0:4), sp(0:4), tq(0:4)
      INTEGER, INTENT(IN) :: d
      INTEGER :: i
      tpart = 0
      DO i = 0, d
         tpart = tpart + a(i,d-i)*sp(i)*tq(d-i)
      END DO
      END FUNCTION tpart

!-------------------------------------------------------------------------------
!>  @brief Regularized lattice sums of the harmonic kernels.
!>
!>  esum(q) = A*Z_k, k = 2*q, q = 0, ..., 4: the regularized sum of
!>  A*|x|**(1-k)*(x1 + i*x2)**k over the lattice spanned by e1 and e2, of cell
!>  area A, less its integral. Z_k is the Epstein zeta function with the
!>  harmonic weight at s = (k - 1)/2, split by Ewald's method between the
!>  lattice and its dual; the terms fall off as EXP(-pi*|x|**2) in the lattice
!>  of unit cell area, and x and -x contribute alike.
!>
!>  @param[in]  e1   First basis vector of the lattice.
!>  @param[in]  e2   Second basis vector of the lattice.
!>  @param[out] esum Sums for k = 0, 2, 4, 6 and 8.
!-------------------------------------------------------------------------------
      SUBROUTINE harmonic_sums(e1, e2, esum)
      REAL(dp), INTENT(IN) :: e1(2), e2(2)
      COMPLEX(dp), INTENT(OUT) :: esum(0:4)
      REAL(dp), PARAMETER :: reach2 = 7
      REAL(dp) :: scale, f1(2), f2(2), d1(2), d2(2), b1(2), b2(2), w(2)
      REAL(dp) :: det, area, mu, g11, g22, g12, gdet, r2, x, xinv
      REAL(dp) :: ex, sx, a, xa, xp, g(0:6), pgam(0:4)
      COMPLEX(dp) :: direct(0:4), dual(0:4), z2, zk, total
      INTEGER :: it, side, i, j, ni, nj, q, m

!     unit cell area, then the Lagrange-Gauss reduced basis
      scale = SQRT(ABS(e1(1)*e2(2) - e1(2)*e2(1)))
      f1 = e1/scale
      f2 = e2/scale
      DO it = 1, 100
         IF (SUM(f1*f1) .GT. SUM(f2*f2)) THEN
            w = f1
            f1 = f2
            f2 = w
            CYCLE
         END IF
         mu = ANINT(SUM(f1*f2)/SUM(f1*f1))
         IF (mu .EQ. 0) THEN
            EXIT
         END IF
         f2 = f2 - mu*f1
      END DO
      det = f1(1)*f2(2) - f1(2)*f2(1)
      d1 = (/ f2(2), -f2(1) /)/det
      d2 = (/ -f1(2), f1(1) /)/det
      area = ABS(det)

      direct = 0
      dual = 0
      DO side = 1, 2
         IF (side .EQ. 1) THEN
            b1 = f1
            b2 = f2
         ELSE
            b1 = d1
            b2 = d2
         END IF
         g11 = SUM(b1*b1)
         g22 = SUM(b2*b2)
         g12 = SUM(b1*b2)
         gdet = g11*g22 - g12*g12
         ni = CEILING(SQRT(reach2*g22/gdet))
         nj = CEILING(SQRT(reach2*g11/gdet))
         DO i = 0, ni
            DO j = -nj, nj
               IF (i .EQ. 0 .AND. j .LE. 0) THEN
                  CYCLE
               END IF
               w = i*b1 + j*b2
               r2 = SUM(w*w)
               IF (r2 .GT. reach2) THEN
                  CYCLE
               END IF
               x = pi*r2
               xinv = 1/x
!              g(m) = GAMMA(m - 1/2, x) by ERFC and the recurrence
               ex = EXP(-x)
               sx = SQRT(x)
               g(1) = SQRT(pi)*ERFC(sx)
               g(0) = 2*(ex/sx - g(1))
               a = 0.5_dp
               xa = sx
               DO m = 2, 6
                  g(m) = a*g(m-1) + xa*ex
                  a = a + 1
                  xa = xa*x
               END DO
               z2 = CMPLX(w(1), w(2), dp)**2
               zk = 1
               IF (side .EQ. 1) THEN
!                 P_k x**(-s) GAMMA(s, x)
                  xp = sx
                  DO q = 0, 4
                     direct(q) = direct(q) + zk*xp*g(q)
                     zk = zk*z2
                     xp = xp*xinv
                  END DO
               ELSE
!                 P_k x**(-(k + 3)/2) GAMMA((k + 3)/2, x)
                  xp = sx*xinv*xinv
                  DO q = 0, 4
                     dual(q) = dual(q) + zk*xp*g(q+2)
                     zk = zk*z2
                     xp = xp*xinv
                  END DO
               END IF
            END DO
         END DO
      END DO

!     pi**s/GAMMA(s) at s = -1/2, 1/2, ..., 7/2
      pgam(0) = -1/(2*pi)
      DO q = 1, 4
         pgam(q) = pi/(q - 1.5_dp)*pgam(q-1)
      END DO
      DO q = 0, 4
!        (-i)**k = (-1)**q; the sums ran over half of each lattice
         total = 2*(direct(q) + (-1)**q*dual(q)/area)
         IF (q .EQ. 0) THEN
            total = total - 2/(3*area) + 2
         END IF
         esum(q) = area*pgam(q)*total*scale**3
      END DO
      END SUBROUTINE harmonic_sums

!-------------------------------------------------------------------------------
!>  @brief Defects of the punctured trapezoidal sums at the source point.
!>
!>  The defects are minus the regularized lattice sums, in the units of the
!>  integral over (s, t), the steps of u and v from the source point.
!>
!>  @param[in]  rt   Taylor coefficients of R in (s, t).
!>  @param[in]  zt   Taylor coefficients of Z in (s, t).
!>  @param[in]  x0   Cartesian x of the source.
!>  @param[in]  y0   Cartesian y of the source.
!>  @param[in]  z0   Cartesian z of the source.
!>  @param[in]  phi0 Toroidal angle of the source.
!>  @param[in]  dphi Rate of the toroidal angle in v.
!>  @param[in]  xn   Cartesian x component of the normal of the source.
!>  @param[in]  yn   Cartesian y component of the normal of the source.
!>  @param[in]  zn   Cartesian z component of the normal of the source.
!>  @param[in]  guu  Coefficient of the u term in the quadratic form of ga1.
!>  @param[in]  guv  Coefficient of the mixed term in the quadratic form of ga1.
!>  @param[in]  gvv  Coefficient of the v term in the quadratic form of ga1.
!>  @param[in]  auu  Coefficient of the u term in the quadratic form of ga2.
!>  @param[in]  auv  Coefficient of the mixed term in the quadratic form of ga2.
!>  @param[in]  avv  Coefficient of the v term in the quadratic form of ga2.
!>  @param[in]  hu   Grid step in u.
!>  @param[in]  hv   Grid step in v.
!>  @param[out] wg   Defects of G - ga1 for the weights 1, s and t.
!>  @param[out] wk   Defects of dG/dn - ga1*ga2 for the weights 1, s and t.
!-------------------------------------------------------------------------------
      SUBROUTINE singular_defects(rt, zt, x0, y0, z0, phi0, dphi,             &
                                  xn, yn, zn, guu, guv, gvv, auu, auv, avv,   &
                                  hu, hv, wg, wk)
      REAL(dp), INTENT(IN) :: rt(0:4,0:4), zt(0:4,0:4), x0, y0, z0, phi0,    &
                              dphi, xn, yn, zn, guu, guv, gvv, auu, auv,    &
                              avv, hu, hv
      REAL(dp), INTENT(OUT) :: wg(3), wk(3)
      REAL(dp), PARAMETER :: fact(0:4) = (/ 1, 1, 2, 6, 24 /)
      REAL(dp), DIMENSION(0:4,0:4) :: cphi, sphi, xc, yc, wc, r2, nr
      REAL(dp) :: rate, ca, cb, cc, th, s, t, sp(0:4), tq(0:4), ri, c3, d4
      REAL(dp) :: a2, a3, a4, s3t, t4, s4, i32, i52, i72, f0g, f1g, f0k
      REAL(dp) :: f1k, phi(6), e1(2), e2(2), el, cq(0:3)
      COMPLEX(dp) :: esum(0:4), coef(6,0:4), e2th, ek
      INTEGER :: j, n, k, m

!     Cartesian position relative to the source and the normal component,
!     as polynomials in (s, t)
!     COS and SIN of phi0 + j*pi/2
      cq = (/ COS(phi0), -SIN(phi0), -COS(phi0), SIN(phi0) /)
      cphi = 0
      sphi = 0
      rate = 1
      DO j = 0, 4
         cphi(0,j) = rate/fact(j)*cq(MOD(j, 4))
         sphi(0,j) = rate/fact(j)*cq(MOD(j + 3, 4))
         rate = rate*dphi
      END DO
      xc = tmul(rt, cphi)
      yc = tmul(rt, sphi)
      wc = zt
      xc(0,0) = xc(0,0) - x0
      yc(0,0) = yc(0,0) - y0
      wc(0,0) = wc(0,0) - z0
      r2 = tmul(xc, xc) + tmul(yc, yc) + tmul(wc, wc)
      nr = xn*xc + yn*yc + zn*wc

!     (ca*s + cb*t, cc*t) are coordinates in which the metric is the identity
      ca = SQRT(guu)
      cb = guv/(2*ca)
      cc = SQRT(gvv - cb*cb)

!     harmonics of the six functions on the unit circle of those coordinates
      coef = 0
      DO n = 0, nang - 1
         th = pi*n/nang
         t = SIN(th)/cc
         s = (COS(th) - cb*t)/ca
         sp(0) = 1
         tq(0) = 1
         DO m = 1, 4
            sp(m) = sp(m-1)*s
            tq(m) = tq(m-1)*t
         END DO
         ri = tpart(r2, 2, sp, tq)
         c3 = tpart(r2, 3, sp, tq)
         d4 = tpart(r2, 4, sp, tq)
         a2 = tpart(nr, 2, sp, tq)
         a3 = tpart(nr, 3, sp, tq)
         a4 = tpart(nr, 4, sp, tq)
!        quartic parts of the quadratic forms of ga1 and ga2 in 2*TAN(s/2)
!        and 2*TAN(t/2)
         s3t = sp(3)*t + s*tq(3)
         t4 = (guu*sp(4) + gvv*tq(4))/6 + guv*s3t/12
         s4 = (auu*sp(4) + avv*tq(4))/6 + auv*s3t/12
         i32 = 1/(ri*SQRT(ri))
         i52 = i32/ri
         i72 = i52/ri
         f0g = -c3*i32/2
         f1g = -(d4 - t4)*i32/2 + 3*c3*c3*i52/8
         f0k = a3*i32 - 3*a2*c3*i52/2
         f1k = (a4 - s4)*i32 - 3*(a3*c3 + a2*(d4 - t4))*i52/2                 &
             + 15*a2*c3*c3*i72/8
         phi = (/ f1g, f0g*s, f0g*t, f1k, f0k*s, f0k*t /)/nang
         e2th = CMPLX(COS(2*th), -SIN(2*th), dp)
         ek = 1
         DO k = 0, 4
            coef(:,k) = coef(:,k) + phi*ek
            ek = ek*e2th
         END DO
      END DO

      e1 = (/ ca*hu, 0.0_dp /)
      e2 = (/ cb*hv, cc*hv /)
      CALL harmonic_sums(e1, e2, esum)
      DO m = 1, 6
         el = REAL(coef(m,0)*esum(0), dp)
         DO k = 1, 4
            el = el + 2*REAL(coef(m,k)*esum(k), dp)
         END DO
         el = -el/(ca*cc)
         IF (m .LE. 3) THEN
            wg(m) = el
         ELSE
            wk(m-3) = el
         END IF
      END DO
      END SUBROUTINE singular_defects

      END MODULE greenf_singular
