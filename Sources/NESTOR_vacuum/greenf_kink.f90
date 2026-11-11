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
!>  Q = a*S**2 + b*S*T + c*T**2 the quadratic form of ga1 and P that of ga2.
!>  They vanish with a kink on the lines s = pi and t = pi, where the Green's
!>  function is smooth, so the trapezoidal sums of greenf are in error there.
!>
!>  On each line of the grid across a kink, values on the node of the kink and
!>  its two neighbours on the line, or on the two nodes on either side of the
!>  kink where it lies midway between nodes, make up the errors of the sums of
!>  g and k along the line for the weights 1, COS and SIN of the angle. Where
!>  the grid step resolves the width of the kink on the line, the errors are
!>  those of the expansion of the approximations about the kink. Where it does
!>  not, about the point where the kinks cross, the sums of the line are set
!>  against the integrals, which are elementary. Near the source the sums also
!>  hold the error of the peak of the approximations, which that of the
!>  Green's function cancels, and they are not used.
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
      PUBLIC :: kink_weights, kink_steps, kink_moments, kink_line,          &
                kink_corner

!>  Pi.
      REAL(dp), PARAMETER :: pi = 3.14159265358979323846264338_dp
!>  Product h*p of the grid step and the modulus p of @ref kink_weights up to
!>  which a line takes the values of the expansion of its kink.
      REAL(dp), PARAMETER :: width_series = 0.5_dp
!>  Product h*p from which the line takes the values of its own sums.
      REAL(dp), PARAMETER :: width_sums = 1
!>  Resolution of the peak at the source on a line, as the logarithm of the
!>  error of the kink over the error of the peak in the sums of the line, up
!>  to which the sums are not used.
      REAL(dp), PARAMETER :: peak_source = 10
!>  Resolution from which the sums of the line are used in full.
      REAL(dp), PARAMETER :: peak_sums = 14
!>  Further resolution required of a line before its values take the ratio of
!>  its own sums for the weights 1 and EXP(-i*t) in place of that of the
!>  expansion. In that ratio the error of the peak counts twice the square of
!>  the harmonic the values are used for.
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
      COMPLEX(dp) :: q0, q1, iq0, iq1, z, w, gm, cj, cjp, rq, bq, p0, p1,   &
                     alpha, ser, sr, sb
      REAL(dp) :: rc, delta, idel
      INTEGER :: j

!     Q, P and their derivatives at t0 = 2*i
      rc = SQRT(c)
      q0 = CMPLX(a*s*s - 4*c, 2*b*s, dp)
      q1 = CMPLX(b*s, 4*c, dp)
      p0 = CMPLX(aa*s*s - 4*cc, 2*bb*s, dp)
      p1 = CMPLX(bb*s, 4*cc, dp)
      delta = s*s*(4*a*c - b*b)/4
      idel = 1/delta
      iq1 = CONJG(q1)/(REAL(q1, dp)**2 + AIMAG(q1)**2)
      z = 4*c*q0*iq1*iq1
!     gm = G - 1, cjp = dJ/dt0, rq = (J - SQRT(c)*Q'/delta)/Q and
!     bq = (dJ/dt0 - 2*c*SQRT(c)/delta - Q'*rq)/Q at t0
      IF (REAL(z, dp)**2 + AIMAG(z)**2 .LT. 1.0E-4_dp) THEN
         ser = 0
         sr = 0
         sb = 0
         DO j = 6, 0, -1
            ser = ser*z + 1/(2*j + 3.0_dp)
            sr = sr*z + (2*j + 2)/(2*j + 3.0_dp)
            sb = sb*z + (2*j + 2)/(2*j + 5.0_dp)
         END DO
         gm = z*ser
         cjp = 8*c*rc*iq1*iq1*ser
         rq = 16*c*rc*iq1*iq1*iq1*sr
         bq = -32*c*c*rc*iq1*iq1*iq1*iq1*sb
      ELSE
         iq0 = CONJG(q0)/(REAL(q0, dp)**2 + AIMAG(q0)**2)
         w = SQRT(z)
         gm = 1 - w
         gm = (1 + w)*CONJG(gm)/(REAL(gm, dp)**2 + AIMAG(gm)**2)
         gm = LOG(gm)*CONJG(w)/(2*(REAL(w, dp)**2 + AIMAG(w)**2)) - 1
         cjp = 2*rc*gm*iq0
         rq = (-4*rc*(gm + 1)*iq1 - rc*idel*q1)*iq0
         bq = (cjp - 2*c*rc*idel - q1*rq)*iq0
      END IF
      cj = -4*rc*(gm + 1)*iq1
      alpha = (c*CMPLX(bb*s, 2*cc, dp) - b*s*cc/2)*idel
      g0 = 2*AIMAG(cj)
      g1 = -4*cjp
      k0 = 2*AIMAG(p0*rq + 2*alpha/rc)
      k1 = -4*(p1*rq + p0*bq + 2*c*cc*idel/rc)
      END SUBROUTINE kink_moments

!-------------------------------------------------------------------------------
!>  @brief Weights of the sums of a line against the expansion of its kink.
!>
!>  Along the line the approximations are singular where Q vanishes, at
!>  T = s*(-b +- i*SQRT(4*a*c - b**2))/(2*c), of modulus p. The approximations
!>  leave the straight slopes of the kink at a distance 4/p from it, so the
!>  expansion of the kink in the grid step h holds while h*p is small, and the
!>  sums of the line are used where it is not. For p < 2 the singular point
!>  lies over the source at a height eta in t, and the sums hold the error of
!>  the peak there as well, which falls as EXP(-n*eta) with the number n of
!>  nodes while that of the kink falls as h**2.
!>
!>  @param[in]  s  Tangent of the other angle on the line.
!>  @param[in]  a  Coefficient of s**2 in Q.
!>  @param[in]  b  Coefficient of s*T in Q.
!>  @param[in]  c  Coefficient of T**2 in Q.
!>  @param[in]  n  Number of nodes on the line.
!>  @param[in]  h  Step of the angle along the line.
!>  @param[out] w0 Weight of the sum of the line for the weight 1: 1 where the
!>                 line takes the values of its own sums, 0 where it takes
!>                 those of the expansion.
!>  @param[out] w1 Weight of the sum for the weight EXP(-i*t) likewise.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE kink_weights(s, a, b, c, n, h, w0, w1)
      REAL(dp), INTENT(IN) :: s, a, b, c, h
      INTEGER, INTENT(IN) :: n
      REAL(dp), INTENT(OUT) :: w0, w1
      REAL(dp) :: x, y, res, t, wk

      x = -b*s/(2*c)
      y = ABS(s)*SQRT(4*a*c - b*b)/(2*c)
      t = MIN(1.0_dp, MAX(0.0_dp, (h*SQRT(x*x + y*y) - width_series)        &
                                  /(width_sums - width_series)))
      wk = t*t*(3 - 2*t)
      w0 = wk
      w1 = wk
      IF (wk .GT. 0 .AND. x*x + y*y .LT. 4) THEN
         res = n*LOG(((2 + y)**2 + x*x)/((2 - y)**2 + x*x))/2 + 2*LOG(h)
         t = MIN(1.0_dp, MAX(0.0_dp, (res - peak_source)                    &
                                     /(peak_sums - peak_source)))
         w0 = wk*t*t*(3 - 2*t)
         t = MIN(1.0_dp, MAX(0.0_dp, (res - peak_source - peak_spread)      &
                                     /(peak_sums - peak_source)))
         w1 = wk*t*t*(3 - 2*t)
      END IF
      END SUBROUTINE kink_weights

!-------------------------------------------------------------------------------
!>  @brief Constants of the lines of one direction for a source point.
!>
!>  Across a kink on a node the trapezoidal sum of |x|, |x|**3 and |x|**5 less
!>  the integral is -h**2/6, h**4/60 and -h**6/126; across a kink midway
!>  between two nodes it is h**2/12, -7*h**4/480 and 31*h**6/4032.
!>
!>  @param[in]  h   Step of the angle along the lines.
!>  @param[in]  n   Number of nodes on a line. A line of fewer than five nodes
!>                  takes the values for the weight 1 alone.
!>  @param[in]  mid True where the kink lies midway between two nodes.
!>  @param[in]  c   Coefficient of T**2 in Q on the lines.
!>  @param[out] st  In 1 to 3, COS(h/2), 1/(2*h*(COS(3*h/2) - COS(h/2))) and
!>                  1/(2*h*SIN(h/2)) where the kink lies between nodes; 1,
!>                  1/(2*h*(1 - COS(h))) and 1/(2*h*SIN(h)) where it lies on a
!>                  node; the last two zero for the weight 1 alone. In 4 to 6,
!>                  the sums less the integrals of |x|, |x|**3 and |x|**5 over
!>                  4*SQRT(c). In 7, the ratio of the sums less the integrals
!>                  of |x|*EXP(-i*t) and |x|.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE kink_steps(h, n, mid, c, st)
      REAL(dp), INTENT(IN) :: h, c
      INTEGER, INTENT(IN) :: n
      LOGICAL, INTENT(IN) :: mid
      REAL(dp), INTENT(OUT) :: st(7)

      st(1:3) = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
      IF (mid) THEN
         st(1) = COS(h/2)
         IF (n .GE. 5) THEN
            st(2) = 1/(2*h*(COS(3*h/2) - COS(h/2)))
            st(3) = 1/(2*h*SIN(h/2))
         END IF
         st(4:6) = (/ h**2/12, -7*h**4/480, 31*h**6/4032 /)
      ELSE
         IF (n .GE. 5) THEN
            st(2) = 1/(2*h*(1 - COS(h)))
            st(3) = 1/(2*h*SIN(h))
         END IF
         st(4:6) = (/ -h**2/6, h**4/60, -h**6/126 /)
      END IF
      st(7) = st(5)/(2*st(4)) - 1
      st(4:6) = st(4:6)/(4*SQRT(c))
      END SUBROUTINE kink_steps

!-------------------------------------------------------------------------------
!>  @brief Values on the nodes about a kink that make up the errors of a line.
!>
!>  The values add h*SUM(v*EXP(-i*n*t)) to the sum of the line for the weight
!>  EXP(-i*n*t), and cancel its error for n = 0 and n = 1, hence for n = -1.
!>
!>  @param[in]  e0  Sum less integral for the weight 1.
!>  @param[in]  e1  Sum less integral for the weight EXP(-i*t).
!>  @param[in]  h   Step of the angle along the line.
!>  @param[in]  mid True where the kink lies midway between two nodes.
!>  @param[in]  st  Constants of the lines of this direction, from
!>                  @ref kink_steps.
!>  @param[out] v   Values for the nodes at -h, 0 and h from the kink, or at
!>                  -3*h/2, -h/2, h/2 and 3*h/2 where it lies between nodes.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE kink_spread(e0, e1, h, mid, st, v)
      REAL(dp), INTENT(IN) :: e0, h, st(7)
      COMPLEX(dp), INTENT(IN) :: e1
      LOGICAL, INTENT(IN) :: mid
      REAL(dp), INTENT(OUT) :: v(4)
      REAL(dp) :: even, odd, inner

      odd = -AIMAG(e1)*st(3)
      IF (mid) THEN
         even = (REAL(e1, dp) + e0*st(1))*st(2)
         inner = -e0/(2*h) - even
         v = (/ even, inner - odd, inner + odd, even /)
      ELSE
         even = -(e0 + REAL(e1, dp))*st(2)
         v = (/ even - odd, -e0/h - 2*even, even + odd, 0.0_dp /)
      END IF
      END SUBROUTINE kink_spread

!-------------------------------------------------------------------------------
!>  @brief Values of the approximations on the nodes about the kink of a line.
!>
!>  For a line of v, s is the tangent of the step of u and a, b, c and aa, bb,
!>  cc are the coefficients of ga1 and ga2 in the order (uu, uv, vv). For a
!>  line of u, s is the tangent of the step of v and the order is (vv, uv, uu).
!>  The sums run over the nodes of the line, and they and the integrals are
!>  read only where w0 is not zero.
!>
!>  In w = 1/T the approximations are |w|/SQRT(c) times a power series: that
!>  of (1 + be*w + ga*w**2)**(-1/2) for g and that of
!>  (cc/c + bp*w + ap*w**2)*(1 + be*w + ga*w**2)**(-3/2) for k, with
!>  be = b*s/c, ga = a*s**2/c, bp = bb*s/c and ap = aa*s**2/c. In the distance
!>  x along the line from the kink w = -TAN(x/2)/2, which turns them into
!>  |x|/(4*SQRT(c)) times a series in x, taken to x**4. The sum less the
!>  integral for the weight 1 follows term by term from those of |x|, |x|**3
!>  and |x|**5, and that for EXP(-i*t) = -COS(x) + i*SIN(x) likewise.
!>
!>  @param[in]  s   Tangent of the other angle on the line.
!>  @param[in]  a   Coefficient of s**2 in the quadratic form of ga1.
!>  @param[in]  b   Coefficient of s*T in the quadratic form of ga1.
!>  @param[in]  c   Coefficient of T**2 in the quadratic form of ga1.
!>  @param[in]  aa  Coefficient of s**2 in the quadratic form of ga2.
!>  @param[in]  bb  Coefficient of s*T in the quadratic form of ga2.
!>  @param[in]  cc  Coefficient of T**2 in the quadratic form of ga2.
!>  @param[in]  h   Step of the angle along the line.
!>  @param[in]  mid True where the kink lies midway between two nodes.
!>  @param[in]  st  Constants of the lines of this direction, from
!>                  @ref kink_steps.
!>  @param[in]  w0  Weight of the sums of the line for the weight 1, from
!>                  @ref kink_weights.
!>  @param[in]  w1  Weight of the sums of the line for the weight EXP(-i*t).
!>  @param[in]  s0  Sums of ga1 and ga1*ga2 over the nodes of the line.
!>  @param[in]  s1  Sums of ga1 and ga1*ga2 times EXP(-i*t).
!>  @param[in]  m0  Integrals of ga1 and ga1*ga2 along the line, from
!>                  @ref kink_moments.
!>  @param[in]  m1  Integrals of ga1 and ga1*ga2 times EXP(-i*t).
!>  @param[out] vg  Values of ga1 for the nodes about the kink, in the order
!>                  of the angle.
!>  @param[out] vk  Values of ga1*ga2 for the same nodes.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE kink_line(s, a, b, c, aa, bb, cc, h, mid, st, w0, w1,   &
                                s0, s1, m0, m1, vg, vk)
      REAL(dp), INTENT(IN) :: s, a, b, c, aa, bb, cc, h, st(7), w0, w1,     &
                              s0(2), m0(2)
      COMPLEX(dp), INTENT(IN) :: s1(2), m1(2)
      LOGICAL, INTENT(IN) :: mid
      REAL(dp), INTENT(OUT) :: vg(4), vk(4)
      REAL(dp) :: be, ga, bp, ap, d(4), sw(2,0:4), sx(2,0:4), e0(2), x0(2)
      COMPLEX(dp) :: e1(2)

!     series in w of g and k over |w|/SQRT(c); d is that of the power -3/2
      be = b*s/c
      ga = a*s*s/c
      bp = bb*s/c
      ap = aa*s*s/c
      sw(1,0) = 1
      sw(1,1) = -be/2
      sw(1,2) = -ga/2 + 3*be*be/8
      sw(1,3) = be*(3*ga/4 - 5*be*be/16)
      sw(1,4) = 3*ga*ga/8 + be*be*(35*be*be/128 - 15*ga/16)
      d(1) = -3*be/2
      d(2) = -3*ga/2 + 15*be*be/8
      d(3) = be*(15*ga/4 - 35*be*be/16)
      d(4) = 15*ga*ga/8 + be*be*(315*be*be/128 - 105*ga/16)
      sw(2,0) = cc/c
      sw(2,1) = bp + sw(2,0)*d(1)
      sw(2,2) = ap + bp*d(1) + sw(2,0)*d(2)
      sw(2,3) = ap*d(1) + bp*d(2) + sw(2,0)*d(3)
      sw(2,4) = ap*d(2) + bp*d(3) + sw(2,0)*d(4)
!     series in x of g and k over |x|/(4*SQRT(c))
      sx(:,0) = sw(:,0)
      sx(:,1) = -sw(:,1)/4
      sx(:,2) = sw(:,0)/12 + sw(:,2)/16
      sx(:,3) = -sw(:,1)/24 - sw(:,3)/64
      sx(:,4) = sw(:,0)/120 + sw(:,2)/64 + sw(:,4)/256
!     sums less integrals for the weights 1 and -COS(x) + i*SIN(x)
      e0 = st(4)*sx(:,0) + st(5)*sx(:,2) + st(6)*sx(:,4)
      e1 = CMPLX(-st(4)*sx(:,0) - st(5)*(sx(:,2) - sx(:,0)/2)               &
                 - st(6)*(sx(:,4) - sx(:,2)/2 + sx(:,0)/24),                &
                 st(5)*sx(:,1) + st(6)*(sx(:,3) - sx(:,1)/6), dp)
      IF (w0 .GT. 0) THEN
         x0 = h*s0 - m0
         e1 = w1*(h*s1 - m1) + (1 - w1)*(e1 + st(7)*w0*(x0 - e0))
         e0 = w0*x0 + (1 - w0)*e0
      END IF
      CALL kink_spread(e0(1), e1(1), h, mid, st, vg)
      CALL kink_spread(e0(2), e1(2), h, mid, st, vk)
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
      REAL(dp), DIMENSION(ndual + 1) :: xm, wm, q, qn, qm, p, pn, pm, r, r2,  &
                                        r3, r5, u1, un1, um1, unm1, u2,    &
                                        un2, um2, unm2
      REAL(dp) :: qa, qb, qc, pa, pb, pc, xn, half, tot(2)
      INTEGER :: n, m, side

!     the last element lies at ndual + 1/2 in m
      DO m = 1, ndual
         xm(m) = m
         wm(m) = 1
         IF (mid .AND. MOD(m, 2) .EQ. 1) THEN
            wm(m) = -1
         END IF
      END DO
      xm(ndual + 1) = ndual + 0.5_dp
      wm(ndual + 1) = 0
      half = (-1)**(ndual + 1)*0.5_dp
      qa = a*hv*hv
      qc = c*hu*hu
      pa = aa*hv*hv
      pc = cc*hu*hu
      tot = 0
      DO side = 1, 2
         qb = (3 - 2*side)*b*hu*hv
         pb = (3 - 2*side)*bb*hu*hv
         DO n = 1, ndual + 1
            xn = n
            IF (n .GT. ndual) THEN
               xn = ndual + 0.5_dp
            END IF
!           U for g and k and its derivatives at the points (xn, xm)
            q = qa*xn*xn - qb*xn*xm + qc*xm*xm
            qn = 2*qa*xn - qb*xm
            qm = 2*qc*xm - qb*xn
            p = pa*xn*xn - pb*xn*xm + pc*xm*xm
            pn = 2*pa*xn - pb*xm
            pm = 2*pc*xm - pb*xn
            r = 1/SQRT(q)
            r2 = r*r
            r3 = r*r2
            r5 = r3*r2
            u1 = r
            un1 = -qn*r3/2
            um1 = -qm*r3/2
            unm1 = 3*qn*qm*r5/4 + qb*r3/2
            u2 = p*r3
            un2 = pn*r3 - 3*p*qn*r5/2
            um2 = pm*r3 - 3*p*qm*r5/2
            unm2 = -pb*r3 - 3*(pm*qn + pn*qm - p*qb)*r5/2                    &
                 + 15*p*qn*qm*r5*r2/4
            m = ndual + 1
            IF (n .LE. ndual) THEN
               tot(1) = tot(1) + SUM(wm*unm1)
               tot(2) = tot(2) + SUM(wm*unm2)
               IF (mid) THEN
                  tot(1) = tot(1) + half*unm1(m)
                  tot(2) = tot(2) + half*unm2(m)
               ELSE
                  tot(1) = tot(1) - un1(m)
                  tot(2) = tot(2) - un2(m)
               END IF
            ELSE
               tot(1) = tot(1) - SUM(wm*um1)
               tot(2) = tot(2) - SUM(wm*um2)
               IF (mid) THEN
                  tot(1) = tot(1) - half*um1(m)
                  tot(2) = tot(2) - half*um2(m)
               ELSE
                  tot(1) = tot(1) + u1(m)
                  tot(2) = tot(2) + u2(m)
               END IF
            END IF
         END DO
      END DO
      cg = hu*hv*tot(1)/(8*pi*pi)
      ck = hu*hv*tot(2)/(8*pi*pi)
      END SUBROUTINE kink_corner

      END MODULE greenf_kink
