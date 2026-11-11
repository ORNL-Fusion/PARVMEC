!*******************************************************************************
!>  @file greenf_kink.f90
!>  @brief Contains module @ref greenf_kink.
!
!  Note separating the Doxygen comment block here so detailed description is
!  found in the Module not the file.
!
!>  Defines the values that the analytic approximations of greenf take on the
!>  nodes about their kinks.
!>
!>  In S = 2*TAN(s/2) and T = 2*TAN(t/2), s and t the steps of u and v from the
!>  source point, the approximations are g = Q**(-1/2) and k = P*Q**(-3/2),
!>  Q the quadratic form of ga1 and P that of ga2. They vanish with a kink on
!>  the lines s = pi and t = pi, where the Green's function is smooth, so the
!>  trapezoidal sums of greenf are in error there.
!>
!>  On a line of the grid that crosses a kink, the sums of g and k times 1, COS
!>  and SIN of the angle along the line are set against the integrals, which
!>  are elementary, and the differences are made up by values on the node of
!>  the kink and its two neighbours on the line, or on the two nodes on either
!>  side of the kink where it lies midway between nodes. Near the source the
!>  sums of a line also hold the error of the peak of the approximations,
!>  which that of the Green's function cancels; there the line takes the
!>  values of a kink between straight slopes, those of the line through the
!>  source.
!>
!>  With every line exact the sum over the grid still misses the integral
!>  about the point where the kinks cross. In the steps x and y of u and v
!>  from that point g tends to |x*y|/(4*SQRT(a*y**2 + b*x*y + c*x**2)). By
!>  Poisson summation the value missing at the crossing node is the sum of the
!>  transform of that function over the dual lattice off its axes, the axes
!>  being the lines of the grid, and the transform is a mixed second
!>  derivative in the indices of the lattice. k follows by differentiating in
!>  a, b and c.
!*******************************************************************************
      MODULE greenf_kink
      USE stel_kinds, ONLY: dp
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: kink_line, kink_corner

!>  Pi.
      REAL(dp), PARAMETER :: pi = 3.14159265358979323846264338_dp
!>  Resolution of the peak at the source on a line, as the logarithm of the
!>  error of the kink over the error of the peak in the sums of the line, up
!>  to which the line takes the values of the line through the source.
      REAL(dp), PARAMETER :: peak_source = 10
!>  Resolution from which the line takes the values of its own sums.
      REAL(dp), PARAMETER :: peak_sums = 14
!>  Further resolution required of a line before its values take the ratio of
!>  its own sums for the weights 1 and EXP(-i*t) in place of that of a kink
!>  between straight slopes. In that ratio the error of the peak counts twice
!>  the square of the harmonic the values are used for.
      REAL(dp), PARAMETER :: peak_spread = 6
!>  Number of points of the dual lattice summed in each index at the crossing
!>  of the kinks; the rest is integrated by the midpoint rule.
      INTEGER, PARAMETER :: ndual = 8

      CONTAINS

!-------------------------------------------------------------------------------
!>  @brief Integrals of the approximations along a line of the grid.
!>
!>  On the line of the tangent s of the other angle, Q = c*T**2 + b*s*T
!>  + a*s**2 and P = cc*T**2 + bb*s*T + aa*s**2. With t0 = 2*i,
!>  z = 4*c*Q(t0)/Q'(t0)**2 and G(z) = ATANH(SQRT(z))/SQRT(z), the integral of
!>  dT/((T - t0)*SQRT(Q)) over the real line is J = -4*SQRT(c)*G/Q'(t0), and
!>  dJ/dt0 = 2*SQRT(c)*(G - 1)/Q(t0). The integral of g over a period of t is
!>  2*AIMAG(J) and that of g*EXP(-i*t) is -4*dJ/dt0. For k,
!>  P/((T - t0)*Q**(3/2)) is reduced to 1/((T - t0)*SQRT(Q)) and a derivative.
!>  Where z is small the differences that vanish with Q(t0) are summed as
!>  series.
!>
!>  @param[in]  s  Tangent of the other angle on the line.
!>  @param[in]  a  Coefficient of s**2 in Q.
!>  @param[in]  b  Coefficient of s*T in Q.
!>  @param[in]  c  Coefficient of T**2 in Q.
!>  @param[in]  aa Coefficient of s**2 in P.
!>  @param[in]  bb Coefficient of s*T in P.
!>  @param[in]  cc Coefficient of T**2 in P.
!>  @param[out] g0 Integral of g over a period of t.
!>  @param[out] g1 Integral of g*EXP(-i*t).
!>  @param[out] k0 Integral of k.
!>  @param[out] k1 Integral of k*EXP(-i*t).
!-------------------------------------------------------------------------------
      PURE SUBROUTINE kink_moments(s, a, b, c, aa, bb, cc, g0, g1, k0, k1)
      REAL(dp), INTENT(IN) :: s, a, b, c, aa, bb, cc
      REAL(dp), INTENT(OUT) :: g0, k0
      COMPLEX(dp), INTENT(OUT) :: g1, k1
      COMPLEX(dp), PARAMETER :: t0 = (0.0_dp, 2.0_dp)
      COMPLEX(dp) :: q0, q1, z, w, gm, cj, cjp, rq, bq, p0, p1, alpha
      REAL(dp) :: rc, delta

      rc = SQRT(c)
      q0 = c*t0*t0 + b*s*t0 + a*s*s
      q1 = 2*c*t0 + b*s
      delta = s*s*(4*a*c - b*b)/4
      z = 4*c*q0/(q1*q1)
!     gm = G - 1, cjp = dJ/dt0, rq = (J - SQRT(c)*Q'/delta)/Q and
!     bq = (dJ/dt0 - 2*c*SQRT(c)/delta - Q'*rq)/Q at t0
      IF (ABS(z) .LT. 0.01_dp) THEN
         gm = z*(1/3.0_dp + z*(1/5.0_dp + z*(1/7.0_dp + z*(1/9.0_dp           &
                                                         + z/11.0_dp))))
         cjp = 8*c*rc/(q1*q1)*(1/3.0_dp + z*(1/5.0_dp + z*(1/7.0_dp         &
                                       + z*(1/9.0_dp + z/11.0_dp))))
         rq = 16*c*rc/q1**3*(2/3.0_dp + z*(4/5.0_dp + z*(6/7.0_dp           &
                                      + z*(8/9.0_dp + z*10/11.0_dp))))
         bq = -32*c*c*rc/q1**4*(2/5.0_dp + z*(4/7.0_dp + z*(6/9.0_dp        &
                                         + z*(8/11.0_dp + z*10/13.0_dp))))
      ELSE
         w = SQRT(z)
         gm = LOG((1 + w)/(1 - w))/(2*w) - 1
         cjp = 2*rc*gm/q0
         rq = (-4*rc*(gm + 1)/q1 - rc*q1/delta)/q0
         bq = (cjp - 2*c*rc/delta - q1*rq)/q0
      END IF
      cj = -4*rc*(gm + 1)/q1
      p0 = cc*t0*t0 + bb*s*t0 + aa*s*s
      p1 = 2*cc*t0 + bb*s
      alpha = (c*(cc*t0 + bb*s) - b*s*cc/2)/delta
      g0 = 2*AIMAG(cj)
      g1 = -4*cjp
      k0 = 2*AIMAG(p0*rq + 2*alpha/rc)
      k1 = -4*(p1*rq + p0*bq + 2*c*cc/(delta*rc))
      END SUBROUTINE kink_moments

!-------------------------------------------------------------------------------
!>  @brief Weights of the sums of a line against the values of the line through
!>  the source.
!>
!>  Along the line the approximations are singular where Q vanishes, at
!>  T = s*(-b +- i*SQRT(4*a*c - b**2))/(2*c). For |T| > 2 that point lies over
!>  the kink and the sums of the line hold no other error. Otherwise it lies
!>  over the source at a height eta in t, and the sums hold the error of the
!>  peak there, which falls as EXP(-n*eta) with the number n of nodes while
!>  that of the kink falls as h**2.
!>
!>  @param[in]  s  Tangent of the other angle on the line.
!>  @param[in]  a  Coefficient of s**2 in Q.
!>  @param[in]  b  Coefficient of s*T in Q.
!>  @param[in]  c  Coefficient of T**2 in Q.
!>  @param[in]  n  Number of nodes on the line.
!>  @param[in]  h  Step of the angle along the line.
!>  @param[out] w0 Weight of the sum of the line for the weight 1: 1 where the
!>                 line takes the values of its own sums, 0 where it takes
!>                 those of the line through the source.
!>  @param[out] w1 Weight of the sum for the weight EXP(-i*t) likewise.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE kink_resolved(s, a, b, c, n, h, w0, w1)
      REAL(dp), INTENT(IN) :: s, a, b, c, h
      INTEGER, INTENT(IN) :: n
      REAL(dp), INTENT(OUT) :: w0, w1
      REAL(dp) :: x, y, res, t

      x = -b*s/(2*c)
      y = ABS(s)*SQRT(4*a*c - b*b)/(2*c)
      IF (x*x + y*y .GE. 4) THEN
         w0 = 1
         w1 = 1
      ELSE
         res = n*LOG(((2 + y)**2 + x*x)/((2 - y)**2 + x*x))/2 + 2*LOG(h)
         t = MIN(1.0_dp, MAX(0.0_dp, (res - peak_source)                    &
                                     /(peak_sums - peak_source)))
         w0 = t*t*(3 - 2*t)
         t = MIN(1.0_dp, MAX(0.0_dp, (res - peak_source - peak_spread)      &
                                     /(peak_sums - peak_source)))
         w1 = t*t*(3 - 2*t)
      END IF
      END SUBROUTINE kink_resolved

!-------------------------------------------------------------------------------
!>  @brief Values on the nodes about a kink that make up the errors of a line.
!>
!>  The values add h*SUM(v*EXP(-i*n*t)) to the sum of the line for the weight
!>  EXP(-i*n*t), and cancel its error for n = 0 and n = 1, hence for n = -1.
!>  A line of fewer than five nodes takes the values for n = 0 alone.
!>
!>  @param[in]  e0  Sum less integral for the weight 1.
!>  @param[in]  e1  Sum less integral for the weight EXP(-i*t).
!>  @param[in]  h   Step of the angle along the line.
!>  @param[in]  n   Number of nodes on the line.
!>  @param[in]  mid True where the kink lies midway between two nodes.
!>  @param[out] v   Values for the nodes at -h, 0 and h from the kink, or at
!>                  -3*h/2, -h/2, h/2 and 3*h/2 where it lies between nodes.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE kink_spread(e0, e1, h, n, mid, v)
      REAL(dp), INTENT(IN) :: e0, h
      COMPLEX(dp), INTENT(IN) :: e1
      INTEGER, INTENT(IN) :: n
      LOGICAL, INTENT(IN) :: mid
      REAL(dp), INTENT(OUT) :: v(4)
      REAL(dp) :: even, odd, inner

      even = 0
      odd = 0
      IF (mid) THEN
         IF (n .GE. 5) THEN
            even = (REAL(e1, dp) + e0*COS(h/2))                              &
                 / (2*h*(COS(3*h/2) - COS(h/2)))
            odd = -AIMAG(e1)/(2*h*SIN(h/2))
         END IF
         inner = -e0/(2*h) - even
         v = (/ even, inner - odd, inner + odd, even /)
      ELSE
         IF (n .GE. 5) THEN
            even = -(e0 + REAL(e1, dp))/(2*h*(1 - COS(h)))
            odd = -AIMAG(e1)/(2*h*SIN(h))
         END IF
         v = (/ even - odd, -e0/h - 2*even, even + odd, 0.0_dp /)
      END IF
      END SUBROUTINE kink_spread

!-------------------------------------------------------------------------------
!>  @brief Values of the approximations on the nodes about the kink of a line.
!>
!>  For a line of v, s is the tangent of the step of u and a, b, c and aa, bb,
!>  cc are the coefficients of ga1 and ga2 in the order (uu, uv, vv). For a
!>  line of u, s is the tangent of the step of v and the order is (vv, uv, uu).
!>  The sums run over the nodes of the line but the source.
!>
!>  A kink of slopes -+S on a node leaves the sum short by h**2*S/6 for the
!>  weight 1 and over by h**2*S/6*(1 + h**2/20) for EXP(-i*t); midway between
!>  nodes it leaves it over by h**2*S/12 and short by h**2*S/12*(1 + 7*h**2/80).
!>  On the line through the source S is 1/(4*SQRT(c)) for g and cc/c times
!>  that for k.
!>
!>  @param[in]  s   Tangent of the other angle on the line.
!>  @param[in]  a   Coefficient of s**2 in the quadratic form of ga1.
!>  @param[in]  b   Coefficient of s*T in the quadratic form of ga1.
!>  @param[in]  c   Coefficient of T**2 in the quadratic form of ga1.
!>  @param[in]  aa  Coefficient of s**2 in the quadratic form of ga2.
!>  @param[in]  bb  Coefficient of s*T in the quadratic form of ga2.
!>  @param[in]  cc  Coefficient of T**2 in the quadratic form of ga2.
!>  @param[in]  n   Number of nodes on the line.
!>  @param[in]  h   Step of the angle along the line.
!>  @param[in]  mid True where the kink lies midway between two nodes.
!>  @param[in]  s0  Sums of ga1 and ga1*ga2 over the nodes of the line.
!>  @param[in]  s1  Sums of ga1 and ga1*ga2 times EXP(-i*t).
!>  @param[out] vg  Values of ga1 for the nodes about the kink, in the order
!>                  of the angle.
!>  @param[out] vk  Values of ga1*ga2 for the same nodes.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE kink_line(s, a, b, c, aa, bb, cc, n, h, mid, s0, s1,    &
                                vg, vk)
      REAL(dp), INTENT(IN) :: s, a, b, c, aa, bb, cc, h, s0(2)
      COMPLEX(dp), INTENT(IN) :: s1(2)
      INTEGER, INTENT(IN) :: n
      LOGICAL, INTENT(IN) :: mid
      REAL(dp), INTENT(OUT) :: vg(4), vk(4)
      REAL(dp) :: w0, w1, e0(2), ratio, g0, k0
      COMPLEX(dp) :: e1(2), g1, k1

      e0(1) = h*h/(24*SQRT(c))
      e0(2) = e0(1)*cc/c
      IF (mid) THEN
         e0 = e0/2
         ratio = -(1 + 7*h*h/80)
      ELSE
         e0 = -e0
         ratio = -(1 + h*h/20)
      END IF
      CALL kink_resolved(s, a, b, c, n, h, w0, w1)
      IF (w0 .GT. 0) THEN
         CALL kink_moments(s, a, b, c, aa, bb, cc, g0, g1, k0, k1)
         e0 = w0*(h*s0 - (/ g0, k0 /)) + (1 - w0)*e0
         e1 = w1*(h*s1 - (/ g1, k1 /)) + (1 - w1)*ratio*e0
      ELSE
         e1 = ratio*e0
      END IF
      CALL kink_spread(e0(1), e1(1), h, n, mid, vg)
      CALL kink_spread(e0(2), e1(2), h, n, mid, vk)
      END SUBROUTINE kink_line

!-------------------------------------------------------------------------------
!>  @brief Values of the approximations at the crossing of the kinks.
!>
!>  The value of g is hu*hv/(8*pi**2) times the sum over n, m >= 1 of the mixed
!>  derivative of U = Q**(-1/2), Q = a*hv**2*n**2 -+ b*hu*hv*n*m
!>  + c*hu**2*m**2, for both signs. The value of k is that of U = P*Q**(-3/2)
!>  with P from aa, bb and cc. The terms to ndual in each index are summed and
!>  the rest integrated by the midpoint rule, which gives derivatives of U at
!>  ndual + 1/2. Where the kink of v lies midway between two nodes the terms
!>  alternate in m and the rest in m is half the next term.
!>
!>  @param[in]  a   Coefficient of the u term in the quadratic form of ga1.
!>  @param[in]  b   Coefficient of the mixed term in the quadratic form of ga1.
!>  @param[in]  c   Coefficient of the v term in the quadratic form of ga1.
!>  @param[in]  aa  Coefficient of the u term in the quadratic form of ga2.
!>  @param[in]  bb  Coefficient of the mixed term in the quadratic form of ga2.
!>  @param[in]  cc  Coefficient of the v term in the quadratic form of ga2.
!>  @param[in]  hu  Grid step in u.
!>  @param[in]  hv  Grid step in v.
!>  @param[in]  mid True where the kink of v lies midway between two nodes.
!>  @param[out] cg  Value of ga1 for the crossing node, or for the two nodes on
!>                  either side of the crossing together.
!>  @param[out] ck  Value of ga1*ga2 likewise.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE kink_corner(a, b, c, aa, bb, cc, hu, hv, mid, cg, ck)
      REAL(dp), INTENT(IN) :: a, b, c, aa, bb, cc, hu, hv
      LOGICAL, INTENT(IN) :: mid
      REAL(dp), INTENT(OUT) :: cg, ck
      REAL(dp) :: qa, qb, qc, pa, pb, pc, edge, tot(2), u(2), un(2), um(2),  &
                  unm(2)
      INTEGER :: n, m, side

      edge = ndual + 0.5_dp
      tot = 0
      DO side = 1, 2
         qa = a*hv*hv
         qb = (3 - 2*side)*b*hu*hv
         qc = c*hu*hu
         pa = aa*hv*hv
         pb = (3 - 2*side)*bb*hu*hv
         pc = cc*hu*hu
         DO n = 1, ndual
            DO m = 1, ndual
               CALL dual(REAL(n, dp), REAL(m, dp), u, un, um, unm)
               IF (mid .AND. MOD(m, 2) .EQ. 1) THEN
                  tot = tot - unm
               ELSE
                  tot = tot + unm
               END IF
            END DO
            CALL dual(REAL(n, dp), edge, u, un, um, unm)
            IF (mid) THEN
               tot = tot + (-1)**(ndual + 1)*unm/2
            ELSE
               tot = tot - un
            END IF
         END DO
         DO m = 1, ndual
            CALL dual(edge, REAL(m, dp), u, un, um, unm)
            IF (mid .AND. MOD(m, 2) .EQ. 1) THEN
               tot = tot + um
            ELSE
               tot = tot - um
            END IF
         END DO
         CALL dual(edge, edge, u, un, um, unm)
         IF (mid) THEN
            tot = tot - (-1)**(ndual + 1)*um/2
         ELSE
            tot = tot + u
         END IF
      END DO
      cg = hu*hv*tot(1)/(8*pi*pi)
      ck = hu*hv*tot(2)/(8*pi*pi)

      CONTAINS

!-------------------------------------------------------------------------------
!>  @brief U for g and k and its derivatives at a point of the dual lattice.
!>
!>  @param[in]  xn  Index dual to u.
!>  @param[in]  xm  Index dual to v.
!>  @param[out] u   U for g and k.
!>  @param[out] un  Derivative of U in n.
!>  @param[out] um  Derivative of U in m.
!>  @param[out] unm Mixed derivative of U.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE dual(xn, xm, u, un, um, unm)
      REAL(dp), INTENT(IN) :: xn, xm
      REAL(dp), INTENT(OUT) :: u(2), un(2), um(2), unm(2)
      REAL(dp) :: q, qn, qm, p, pn, pm, r, r3, r5

      q = qa*xn*xn - qb*xn*xm + qc*xm*xm
      qn = 2*qa*xn - qb*xm
      qm = 2*qc*xm - qb*xn
      p = pa*xn*xn - pb*xn*xm + pc*xm*xm
      pn = 2*pa*xn - pb*xm
      pm = 2*pc*xm - pb*xn
      r = 1/SQRT(q)
      r3 = r/q
      r5 = r3/q
      u(1) = r
      un(1) = -qn*r3/2
      um(1) = -qm*r3/2
      unm(1) = 3*qn*qm*r5/4 + qb*r3/2
      u(2) = p*r3
      un(2) = pn*r3 - 3*p*qn*r5/2
      um(2) = pm*r3 - 3*p*qm*r5/2
      unm(2) = -pb*r3 - 3*(pm*qn + pn*qm - p*qb)*r5/2                        &
             + 15*p*qn*qm*r5/(4*q)
      END SUBROUTINE dual

      END SUBROUTINE kink_corner

      END MODULE greenf_kink
