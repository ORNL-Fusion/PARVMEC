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
!>
!>  @ref kink_setup, @ref kink_grid, @ref kink_image and @ref kink_image_sums
!>  place the values on the grid of greenf for one source point.
!*******************************************************************************
      MODULE greenf_kink
      USE stel_kinds, ONLY: dp
      USE stel_constants, ONLY: pi

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: kink_context
      PUBLIC :: kink_setup
      PUBLIC :: kink_grid
      PUBLIC :: kink_image
      PUBLIC :: kink_image_sums

!>  Product h*p of the grid step and the modulus p of @ref kink_weights up to
!>  which a line takes the values of the expansion of its kink.
      REAL(dp), PARAMETER :: width_series = 0.5_dp
!>  Product h*p from which the line takes the values of its own sums.
      REAL(dp), PARAMETER :: width_sums = 1.0_dp
!>  Resolution of the peak at the source on a line, as the logarithm of the
!>  error of the kink over the error of the peak in the sums of the line, up
!>  to which the sums are not used.
      REAL(dp), PARAMETER :: peak_source = 10.0_dp
!>  Resolution from which the sums of the line are used in full.
      REAL(dp), PARAMETER :: peak_sums = 14.0_dp
!>  Further resolution required of a line before its values take the ratio of
!>  its own sums for the weights 1 and EXP(-i*t) in place of that of the
!>  expansion. In that ratio the error of the peak counts twice the square of
!>  the harmonic the values are used for.
      REAL(dp), PARAMETER :: peak_spread = 6.0_dp
!>  Number of points of the dual lattice summed in each index at the crossing
!>  of the kinks; the rest is integrated by the midpoint rule.
      INTEGER, PARAMETER  :: ndual = 8

!-------------------------------------------------------------------------------
!>  Quantities of the values about the kinks for one source point. Those with
!>  a last index of two hold the lines of v first and the lines of u second.
!-------------------------------------------------------------------------------
      TYPE :: kink_context
!>  Index in u of the source point, from zero.
         INTEGER                                    :: ku0
!>  Index in v of the source point, from zero.
         INTEGER                                    :: kv0
!>  Offset in tanu of the steps of u from the source point.
         INTEGER                                    :: iuoff
!>  Offset in tanv of the steps of v from the source point.
         INTEGER                                    :: ivoff
!>  Offset of the nodes on the kink of u.
         INTEGER                                    :: kua
!>  Index in u, before it wraps, of the first node about the kink of a line
!>  of u.
         INTEGER                                    :: kfu
!>  Index in v, before it wraps, of the first node about the kink of a line
!>  of v.
         INTEGER                                    :: kfv
!>  True where the kink of a line lies midway between two nodes.
         LOGICAL, DIMENSION(2)                      :: mid
!>  Step of the angle along a line. For nv = 1 that of v is the spacing of
!>  the images.
         REAL(dp), DIMENSION(2)                     :: h
!>  Coefficients of the quadratic form of ga1 in the order of @ref kink_line.
         REAL(dp), DIMENSION(3,2)                   :: g
!>  Coefficients of the quadratic form of ga2 in the same order.
         REAL(dp), DIMENSION(3,2)                   :: p
!>  Constants of the lines, from @ref kink_steps.
         REAL(dp), DIMENSION(7,2)                   :: st
!>  EXP(-i*t) at the steps t of the angle along a line from the source.
         COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE   :: w
!>  For nv = 1, sums of ga1 and ga1*ga2 over the images by node of u.
         REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: rows
!>  For nv = 1, weights of those sums, from @ref kink_weights.
         REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: wrow
!>  Tangent of the line whose integrals are kept, by its step from the
!>  source. Zero where none are kept.
         REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: skept
!>  Integrals kept for the weight 1.
         REAL(dp), DIMENSION(:,:,:), ALLOCATABLE    :: m0
!>  Integrals kept for the weight EXP(-i*t).
         COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: m1
      END TYPE

      CONTAINS

!-------------------------------------------------------------------------------
!>  @brief Smooth step from zero to one.
!>
!>  @param[in] x Argument.
!>  @returns Zero for x below zero, one for x above one and x**2*(3 - 2*x)
!>           between.
!-------------------------------------------------------------------------------
      PURE FUNCTION smooth_step(x)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp)             :: smooth_step
      REAL(dp), INTENT(IN) :: x

!  Local Variables
      REAL(dp)             :: t

!  Start of executable code
      t = MIN(1.0_dp, MAX(0.0_dp, x))
      smooth_step = t*t*(3 - 2*t)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Factors of the integrals of the approximations along a line.
!>
!>  With t0 = 2*i, z = 4*c*Q(t0)/Q'(t0)**2 and G(z) = ATANH(SQRT(z))/SQRT(z),
!>  the integral of dT/((T - t0)*SQRT(Q)) over the real line is
!>  J = -4*SQRT(c)*G/Q'(t0), and dJ/dt0 = 2*SQRT(c)*(G - 1)/Q(t0). Where z is
!>  small the differences that vanish with Q(t0) are summed as series.
!>
!>  @param[in] q0    Q at t0.
!>  @param[in] q1    Derivative of Q at t0.
!>  @param[in] c     Coefficient of T**2 in Q.
!>  @param[in] delta Discriminant of Q over four, with its sign reversed.
!>  @returns J, dJ/dt0, (J - SQRT(c)*Q'/delta)/Q and
!>           (dJ/dt0 - 2*c*SQRT(c)/delta - Q'*(J - SQRT(c)*Q'/delta)/Q)/Q at
!>           t0.
!-------------------------------------------------------------------------------
      PURE FUNCTION moment_factors(q0, q1, c, delta)

      IMPLICIT NONE

!  Declare Arguments
      COMPLEX(dp), DIMENSION(4) :: moment_factors
      COMPLEX(dp), INTENT(IN)   :: q0
      COMPLEX(dp), INTENT(IN)   :: q1
      REAL(dp), INTENT(IN)      :: c
      REAL(dp), INTENT(IN)      :: delta

!  Local Variables
      COMPLEX(dp)               :: iq0
      COMPLEX(dp)               :: iq1
      COMPLEX(dp)               :: z
      COMPLEX(dp)               :: w
      COMPLEX(dp)               :: gm
      COMPLEX(dp)               :: cjp
      COMPLEX(dp)               :: rq
      COMPLEX(dp)               :: bq
      COMPLEX(dp)               :: ser
      COMPLEX(dp)               :: sr
      COMPLEX(dp)               :: sb
      REAL(dp)                  :: rc
      REAL(dp)                  :: idel
      INTEGER                   :: j

!  Start of executable code
      rc = SQRT(c)
      idel = 1/delta
      iq1 = CONJG(q1)/(REAL(q1, dp)*REAL(q1, dp) + AIMAG(q1)*AIMAG(q1))
      z = 4*c*q0*iq1*iq1

!  gm is G - 1, cjp the second factor, rq the third and bq the fourth
      IF (REAL(z, dp)*REAL(z, dp) + AIMAG(z)*AIMAG(z) .LT. 1.0E-4_dp) THEN
         ser = 0.0_dp
         sr = 0.0_dp
         sb = 0.0_dp
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
         iq0 = CONJG(q0)/(REAL(q0, dp)*REAL(q0, dp) + AIMAG(q0)*AIMAG(q0))
         w = SQRT(z)
         gm = 1 - w
         gm = (1 + w)*CONJG(gm)/(REAL(gm, dp)*REAL(gm, dp) +                   &
                                 AIMAG(gm)*AIMAG(gm))
         gm = LOG(gm)*CONJG(w)/(2*(REAL(w, dp)*REAL(w, dp) +                   &
                                   AIMAG(w)*AIMAG(w))) - 1
         cjp = 2*rc*gm*iq0
         rq = (-4*rc*(gm + 1)*iq1 - rc*idel*q1)*iq0
         bq = (cjp - 2*c*rc*idel - q1*rq)*iq0
      END IF

      moment_factors(1) = -4*rc*(gm + 1)*iq1
      moment_factors(2) = cjp
      moment_factors(3) = rq
      moment_factors(4) = bq

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Integrals of the approximations along a line of the grid.
!>
!>  On the line of the tangent s of the other angle, Q = c*T**2 + b*s*T
!>  + a*s**2 and P = cc*T**2 + bb*s*T + aa*s**2. With the factors of
!>  @ref moment_factors at t0 = 2*i, the integral of g over a period of t is
!>  2*AIMAG(J) and that of g*EXP(-i*t) is -4*dJ/dt0. For k,
!>  P/((T - t0)*Q**(3/2)) is reduced to 1/((T - t0)*SQRT(Q)) and a derivative.
!>
!>  @param[in]  s  Tangent of the other angle on the line.
!>  @param[in]  g  Coefficients a, b and c of Q.
!>  @param[in]  p  Coefficients aa, bb and cc of P.
!>  @param[out] m0 Integrals of g and of k over a period of t.
!>  @param[out] m1 Integrals of g*EXP(-i*t) and of k*EXP(-i*t).
!-------------------------------------------------------------------------------
      PURE SUBROUTINE kink_moments(s, g, p, m0, m1)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), INTENT(IN)                   :: s
      REAL(dp), DIMENSION(3), INTENT(IN)     :: g
      REAL(dp), DIMENSION(3), INTENT(IN)     :: p
      REAL(dp), DIMENSION(2), INTENT(OUT)    :: m0
      COMPLEX(dp), DIMENSION(2), INTENT(OUT) :: m1

!  Local Variables
      COMPLEX(dp), DIMENSION(4)              :: f
      COMPLEX(dp)                            :: q0
      COMPLEX(dp)                            :: q1
      COMPLEX(dp)                            :: p0
      COMPLEX(dp)                            :: p1
      COMPLEX(dp)                            :: alpha
      REAL(dp)                               :: rc
      REAL(dp)                               :: delta

!  Start of executable code
!  Q, P and their derivatives at t0 = 2*i
      q0 = CMPLX(g(1)*s*s - 4*g(3), 2*g(2)*s, dp)
      q1 = CMPLX(g(2)*s, 4*g(3), dp)
      p0 = CMPLX(p(1)*s*s - 4*p(3), 2*p(2)*s, dp)
      p1 = CMPLX(p(2)*s, 4*p(3), dp)
      rc = SQRT(g(3))
      delta = s*s*(4*g(1)*g(3) - g(2)*g(2))/4

      f = moment_factors(q0, q1, g(3), delta)
      alpha = (g(3)*CMPLX(p(2)*s, 2*p(3), dp) - g(2)*s*p(3)/2)/delta
      m0(1) = 2*AIMAG(f(1))
      m1(1) = -4*f(2)
      m0(2) = 2*AIMAG(p0*f(3) + 2*alpha/rc)
      m1(2) = -4*(p1*f(3) + p0*f(4) + 2*g(3)*p(3)/(delta*rc))

      END SUBROUTINE

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
!>  @param[in] s Tangent of the other angle on the line.
!>  @param[in] g Coefficients a, b and c of Q.
!>  @param[in] n Number of nodes on the line.
!>  @param[in] h Step of the angle along the line.
!>  @returns The weights of the sums of the line for the weights 1 and
!>           EXP(-i*t): 1 where the line takes the values of its own sums, 0
!>           where it takes those of the expansion.
!-------------------------------------------------------------------------------
      PURE FUNCTION kink_weights(s, g, n, h)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(2)             :: kink_weights
      REAL(dp), INTENT(IN)               :: s
      REAL(dp), DIMENSION(3), INTENT(IN) :: g
      INTEGER, INTENT(IN)                :: n
      REAL(dp), INTENT(IN)               :: h

!  Local Variables
      REAL(dp)                           :: x
      REAL(dp)                           :: y
      REAL(dp)                           :: res
      REAL(dp)                           :: wk

!  Start of executable code
      x = -g(2)*s/(2*g(3))
      y = ABS(s)*SQRT(4*g(1)*g(3) - g(2)*g(2))/(2*g(3))
      wk = smooth_step((h*SQRT(x*x + y*y) - width_series)/                     &
                       (width_sums - width_series))
      kink_weights = wk
      IF (wk .GT. 0 .AND. x*x + y*y .LT. 4) THEN
         res = n*LOG(((2 + y)*(2 + y) + x*x)/((2 - y)*(2 - y) + x*x))/2 +      &
               2*LOG(h)
         kink_weights(1) = wk*smooth_step((res - peak_source)/                 &
                                          (peak_sums - peak_source))
         kink_weights(2) = wk*smooth_step((res - peak_source - peak_spread)/   &
                                          (peak_sums - peak_source))
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Constants of the lines of one direction for a source point.
!>
!>  Across a kink on a node the trapezoidal sum of |x|, |x|**3 and |x|**5 less
!>  the integral is -h**2/6, h**4/60 and -h**6/126; across a kink midway
!>  between two nodes it is h**2/12, -7*h**4/480 and 31*h**6/4032.
!>
!>  @param[in] h   Step of the angle along the lines.
!>  @param[in] n   Number of nodes on a line. A line of fewer than five nodes
!>                 takes the values for the weight 1 alone.
!>  @param[in] mid True where the kink lies midway between two nodes.
!>  @param[in] c   Coefficient of T**2 in Q on the lines.
!>  @returns In 1 to 3, COS(h/2), 1/(2*h*(COS(3*h/2) - COS(h/2))) and
!>           1/(2*h*SIN(h/2)) where the kink lies between nodes; 1,
!>           1/(2*h*(1 - COS(h))) and 1/(2*h*SIN(h)) where it lies on a node;
!>           the last two zero for the weight 1 alone. In 4 to 6, the sums
!>           less the integrals of |x|, |x|**3 and |x|**5 over 4*SQRT(c). In
!>           7, the ratio of the sums less the integrals of |x|*EXP(-i*t) and
!>           |x|.
!-------------------------------------------------------------------------------
      PURE FUNCTION kink_steps(h, n, mid, c)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(7) :: kink_steps
      REAL(dp), INTENT(IN)   :: h
      INTEGER, INTENT(IN)    :: n
      LOGICAL, INTENT(IN)    :: mid
      REAL(dp), INTENT(IN)   :: c

!  Local Variables
      REAL(dp)               :: h2
      REAL(dp)               :: h4
      REAL(dp)               :: h6

!  Start of executable code
      h2 = h*h
      h4 = h2*h2
      h6 = h4*h2
      kink_steps(1:3) = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
      IF (mid) THEN
         kink_steps(1) = COS(h/2)
         IF (n .GE. 5) THEN
            kink_steps(2) = 1/(2*h*(COS(3*h/2) - COS(h/2)))
            kink_steps(3) = 1/(2*h*SIN(h/2))
         END IF
         kink_steps(4:6) = (/ h2/12, -7*h4/480, 31*h6/4032 /)
      ELSE
         IF (n .GE. 5) THEN
            kink_steps(2) = 1/(2*h*(1 - COS(h)))
            kink_steps(3) = 1/(2*h*SIN(h))
         END IF
         kink_steps(4:6) = (/ -h2/6, h4/60, -h6/126 /)
      END IF
      kink_steps(7) = kink_steps(5)/(2*kink_steps(4)) - 1
      kink_steps(4:6) = kink_steps(4:6)/(4*SQRT(c))

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Values on the nodes about a kink that make up the errors of a line.
!>
!>  The values add h*SUM(v*EXP(-i*n*t)) to the sum of the line for the weight
!>  EXP(-i*n*t), and cancel its error for n = 0 and n = 1, hence for n = -1.
!>
!>  @param[in] e0  Sum less integral for the weight 1.
!>  @param[in] e1  Sum less integral for the weight EXP(-i*t).
!>  @param[in] h   Step of the angle along the line.
!>  @param[in] mid True where the kink lies midway between two nodes.
!>  @param[in] st  Constants of the lines of this direction, from
!>                 @ref kink_steps.
!>  @returns The values for the nodes at -h, 0 and h from the kink, or at
!>           -3*h/2, -h/2, h/2 and 3*h/2 where it lies between nodes.
!-------------------------------------------------------------------------------
      PURE FUNCTION kink_spread(e0, e1, h, mid, st)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(4)             :: kink_spread
      REAL(dp), INTENT(IN)               :: e0
      COMPLEX(dp), INTENT(IN)            :: e1
      REAL(dp), INTENT(IN)               :: h
      LOGICAL, INTENT(IN)                :: mid
      REAL(dp), DIMENSION(7), INTENT(IN) :: st

!  Local Variables
      REAL(dp)                           :: even
      REAL(dp)                           :: odd
      REAL(dp)                           :: inner

!  Start of executable code
      odd = -AIMAG(e1)*st(3)
      IF (mid) THEN
         even = (REAL(e1, dp) + e0*st(1))*st(2)
         inner = -e0/(2*h) - even
         kink_spread = (/ even, inner - odd, inner + odd, even /)
      ELSE
         even = -(e0 + REAL(e1, dp))*st(2)
         kink_spread = (/ even - odd, -e0/h - 2*even, even + odd, 0.0_dp /)
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Expansion of the approximations about the kink of a line.
!>
!>  In w = 1/T the approximations are |w|/SQRT(c) times a power series: that
!>  of (1 + be*w + ga*w**2)**(-1/2) for g and that of
!>  (cc/c + bp*w + ap*w**2)*(1 + be*w + ga*w**2)**(-3/2) for k, with
!>  be = b*s/c, ga = a*s**2/c, bp = bb*s/c and ap = aa*s**2/c. In the distance
!>  x along the line from the kink w = -TAN(x/2)/2, which turns them into
!>  |x|/(4*SQRT(c)) times a series in x.
!>
!>  @param[in] s Tangent of the other angle on the line.
!>  @param[in] g Coefficients a, b and c of Q.
!>  @param[in] p Coefficients aa, bb and cc of P.
!>  @returns The coefficients of x**0 to x**4 of the series of g and of k.
!-------------------------------------------------------------------------------
      PURE FUNCTION kink_expansion(s, g, p)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(2,0:4)         :: kink_expansion
      REAL(dp), INTENT(IN)               :: s
      REAL(dp), DIMENSION(3), INTENT(IN) :: g
      REAL(dp), DIMENSION(3), INTENT(IN) :: p

!  Local Variables
      REAL(dp), DIMENSION(2,0:4)         :: sw
      REAL(dp), DIMENSION(4)             :: d
      REAL(dp)                           :: be
      REAL(dp)                           :: ga
      REAL(dp)                           :: bp
      REAL(dp)                           :: ap

!  Start of executable code
!  Series in w of g and k over |w|/SQRT(c); d is that of the power -3/2
      be = g(2)*s/g(3)
      ga = g(1)*s*s/g(3)
      bp = p(2)*s/g(3)
      ap = p(1)*s*s/g(3)
      sw(1,0) = 1.0_dp
      sw(1,1) = -be/2
      sw(1,2) = -ga/2 + 3*be*be/8
      sw(1,3) = be*(3*ga/4 - 5*be*be/16)
      sw(1,4) = 3*ga*ga/8 + be*be*(35*be*be/128 - 15*ga/16)
      d(1) = -3*be/2
      d(2) = -3*ga/2 + 15*be*be/8
      d(3) = be*(15*ga/4 - 35*be*be/16)
      d(4) = 15*ga*ga/8 + be*be*(315*be*be/128 - 105*ga/16)
      sw(2,0) = p(3)/g(3)
      sw(2,1) = bp + sw(2,0)*d(1)
      sw(2,2) = ap + bp*d(1) + sw(2,0)*d(2)
      sw(2,3) = ap*d(1) + bp*d(2) + sw(2,0)*d(3)
      sw(2,4) = ap*d(2) + bp*d(3) + sw(2,0)*d(4)

!  Series in x of g and k over |x|/(4*SQRT(c))
      kink_expansion(:,0) = sw(:,0)
      kink_expansion(:,1) = -sw(:,1)/4
      kink_expansion(:,2) = sw(:,0)/12 + sw(:,2)/16
      kink_expansion(:,3) = -sw(:,1)/24 - sw(:,3)/64
      kink_expansion(:,4) = sw(:,0)/120 + sw(:,2)/64 + sw(:,4)/256

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Values of the approximations on the nodes about the kink of a line.
!>
!>  For a line of v, s is the tangent of the step of u and g and p hold the
!>  coefficients of ga1 and ga2 in the order (uu, uv, vv). For a line of u, s
!>  is the tangent of the step of v and the order is (vv, uv, uu).
!>
!>  The sum less the integral of the expansion of @ref kink_expansion follows
!>  term by term from those of |x|, |x|**3 and |x|**5 for the weight 1, and
!>  likewise for EXP(-i*t) = -COS(x) + i*SIN(x). The sums of the line less its
!>  integrals take their place by the weights of @ref kink_weights.
!>
!>  @param[in] s   Tangent of the other angle on the line.
!>  @param[in] g   Coefficients a, b and c of Q.
!>  @param[in] p   Coefficients aa, bb and cc of P.
!>  @param[in] h   Step of the angle along the line.
!>  @param[in] mid True where the kink lies midway between two nodes.
!>  @param[in] st  Constants of the lines of this direction, from
!>                 @ref kink_steps.
!>  @param[in] w   Weights of the sums of the line, from @ref kink_weights.
!>  @param[in] x0  Sums of ga1 and of ga1*ga2 over the nodes of the line, times
!>                 h, less the integrals. Read where w is not zero.
!>  @param[in] x1  The same for the sums and integrals times EXP(-i*t).
!>  @returns The values of ga1 and of ga1*ga2 for the nodes about the kink, in
!>           the order of the angle.
!-------------------------------------------------------------------------------
      PURE FUNCTION kink_line(s, g, p, h, mid, st, w, x0, x1)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(4,2)              :: kink_line
      REAL(dp), INTENT(IN)                  :: s
      REAL(dp), DIMENSION(3), INTENT(IN)    :: g
      REAL(dp), DIMENSION(3), INTENT(IN)    :: p
      REAL(dp), INTENT(IN)                  :: h
      LOGICAL, INTENT(IN)                   :: mid
      REAL(dp), DIMENSION(7), INTENT(IN)    :: st
      REAL(dp), DIMENSION(2), INTENT(IN)    :: w
      REAL(dp), DIMENSION(2), INTENT(IN)    :: x0
      COMPLEX(dp), DIMENSION(2), INTENT(IN) :: x1

!  Local Variables
      REAL(dp), DIMENSION(2,0:4)            :: sx
      REAL(dp), DIMENSION(2)                :: e0
      COMPLEX(dp), DIMENSION(2)             :: e1
      INTEGER                               :: k

!  Start of executable code
!  Sums less integrals of the expansion for the weights 1 and
!  -COS(x) + i*SIN(x)
      sx = kink_expansion(s, g, p)
      e0 = st(4)*sx(:,0) + st(5)*sx(:,2) + st(6)*sx(:,4)
      e1 = CMPLX(-st(4)*sx(:,0) - st(5)*(sx(:,2) - sx(:,0)/2) -                &
                 st(6)*(sx(:,4) - sx(:,2)/2 + sx(:,0)/24),                     &
                 st(5)*sx(:,1) + st(6)*(sx(:,3) - sx(:,1)/6), dp)

      IF (w(1) .GT. 0) THEN
         e1 = w(2)*x1 + (1 - w(2))*(e1 + st(7)*w(1)*(x0 - e0))
         e0 = w(1)*x0 + (1 - w(1))*e0
      END IF

      DO k = 1, 2
         kink_line(:,k) = kink_spread(e0(k), e1(k), h, mid, st)
      END DO

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Summands of the values at the crossing of the kinks.
!>
!>  U is Q**(-1/2) for g and P*Q**(-3/2) for k, with
!>  Q = q(1)*n**2 - q(2)*n*m + q(3)*m**2 and P likewise from p.
!>
!>  @param[in] q  Coefficients of Q.
!>  @param[in] p  Coefficients of P.
!>  @param[in] xn Index n of the dual lattice.
!>  @returns U, dU/dn, dU/dm and d2U/(dn dm) of g and of k at m = 1, ...,
!>           ndual and at m = ndual + 1/2.
!-------------------------------------------------------------------------------
      PURE FUNCTION dual_summands(q, p, xn)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(ndual + 1,4,2) :: dual_summands
      REAL(dp), DIMENSION(3), INTENT(IN) :: q
      REAL(dp), DIMENSION(3), INTENT(IN) :: p
      REAL(dp), INTENT(IN)               :: xn

!  Local Variables
      REAL(dp), DIMENSION(ndual + 1)     :: xm
      REAL(dp), DIMENSION(ndual + 1)     :: qv
      REAL(dp), DIMENSION(ndual + 1)     :: qn
      REAL(dp), DIMENSION(ndual + 1)     :: qm
      REAL(dp), DIMENSION(ndual + 1)     :: pv
      REAL(dp), DIMENSION(ndual + 1)     :: pn
      REAL(dp), DIMENSION(ndual + 1)     :: pm
      REAL(dp), DIMENSION(ndual + 1)     :: r
      REAL(dp), DIMENSION(ndual + 1)     :: r2
      REAL(dp), DIMENSION(ndual + 1)     :: r3
      REAL(dp), DIMENSION(ndual + 1)     :: r5
      INTEGER                            :: m

!  Start of executable code
      DO m = 1, ndual
         xm(m) = m
      END DO
      xm(ndual + 1) = ndual + 0.5_dp

      qv = q(1)*xn*xn - q(2)*xn*xm + q(3)*xm*xm
      qn = 2*q(1)*xn - q(2)*xm
      qm = 2*q(3)*xm - q(2)*xn
      pv = p(1)*xn*xn - p(2)*xn*xm + p(3)*xm*xm
      pn = 2*p(1)*xn - p(2)*xm
      pm = 2*p(3)*xm - p(2)*xn
      r = 1/SQRT(qv)
      r2 = r*r
      r3 = r*r2
      r5 = r3*r2

      dual_summands(:,1,1) = r
      dual_summands(:,2,1) = -qn*r3/2
      dual_summands(:,3,1) = -qm*r3/2
      dual_summands(:,4,1) = 3*qn*qm*r5/4 + q(2)*r3/2
      dual_summands(:,1,2) = pv*r3
      dual_summands(:,2,2) = pn*r3 - 3*pv*qn*r5/2
      dual_summands(:,3,2) = pm*r3 - 3*pv*qm*r5/2
      dual_summands(:,4,2) = -p(2)*r3 - 3*(pm*qn + pn*qm - pv*q(2))*r5/2 +     &
                             15*pv*qn*qm*r5*r2/4

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Sum over a quadrant of the dual lattice at the crossing of the kinks.
!>
!>  The mixed derivative of U is summed over n, m = 1, ..., ndual and the rest
!>  is integrated by the midpoint rule, which gives derivatives of U at
!>  ndual + 1/2. Where the kink of v lies midway between two nodes the terms
!>  alternate in m and the rest in m is half the next term.
!>
!>  @param[in] q   Coefficients of Q in @ref dual_summands.
!>  @param[in] p   Coefficients of P in @ref dual_summands.
!>  @param[in] mid True where the kink of v lies midway between two nodes.
!>  @returns The sums of g and of k.
!-------------------------------------------------------------------------------
      PURE FUNCTION dual_sums(q, p, mid)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(2)             :: dual_sums
      REAL(dp), DIMENSION(3), INTENT(IN) :: q
      REAL(dp), DIMENSION(3), INTENT(IN) :: p
      LOGICAL, INTENT(IN)                :: mid

!  Local Variables
      REAL(dp), DIMENSION(ndual + 1,4,2) :: u
      REAL(dp), DIMENSION(ndual)         :: wm
      INTEGER                            :: n
      INTEGER                            :: m
      INTEGER                            :: k

!  Local Parameters
!  Factor of the next term in the rest of an alternating sum.
      REAL(dp), PARAMETER                :: half = 0.5_dp - MOD(ndual + 1, 2)

!  Start of executable code
      wm = 1.0_dp
      IF (mid) THEN
         wm(1:ndual:2) = -1.0_dp
      END IF
      m = ndual + 1

      dual_sums = 0.0_dp
      DO n = 1, ndual
         u = dual_summands(q, p, REAL(n, dp))
         DO k = 1, 2
            dual_sums(k) = dual_sums(k) + SUM(wm*u(1:ndual,4,k))
            IF (mid) THEN
               dual_sums(k) = dual_sums(k) + half*u(m,4,k)
            ELSE
               dual_sums(k) = dual_sums(k) - u(m,2,k)
            END IF
         END DO
      END DO

!  Rest of the sum in n
      u = dual_summands(q, p, ndual + 0.5_dp)
      DO k = 1, 2
         dual_sums(k) = dual_sums(k) - SUM(wm*u(1:ndual,3,k))
         IF (mid) THEN
            dual_sums(k) = dual_sums(k) - half*u(m,3,k)
         ELSE
            dual_sums(k) = dual_sums(k) + u(m,1,k)
         END IF
      END DO

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Values of the approximations at the crossing of the kinks.
!>
!>  The value of g is hu*hv/(8*pi**2) times the sum over n, m >= 1 of the mixed
!>  derivative of U = Q**(-1/2), Q = a*hv**2*n**2 -+ b*hu*hv*n*m
!>  + c*hu**2*m**2, for both signs. The value of k is that of U = P*Q**(-3/2)
!>  with P from aa, bb and cc.
!>
!>  @param[in] g   Coefficients of the u, mixed and v terms in the quadratic
!>                 form of ga1.
!>  @param[in] p   Coefficients of the u, mixed and v terms in the quadratic
!>                 form of ga2.
!>  @param[in] hu  Grid step in u.
!>  @param[in] hv  Grid step in v.
!>  @param[in] mid True where the kink of v lies midway between two nodes.
!>  @returns The values of ga1 and of ga1*ga2 for the crossing node, or for
!>           the two nodes on either side of the crossing together.
!-------------------------------------------------------------------------------
      PURE FUNCTION kink_corner(g, p, hu, hv, mid)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(2)             :: kink_corner
      REAL(dp), DIMENSION(3), INTENT(IN) :: g
      REAL(dp), DIMENSION(3), INTENT(IN) :: p
      REAL(dp), INTENT(IN)               :: hu
      REAL(dp), INTENT(IN)               :: hv
      LOGICAL, INTENT(IN)                :: mid

!  Local Variables
      REAL(dp), DIMENSION(3)             :: q
      REAL(dp), DIMENSION(3)             :: pq

!  Start of executable code
      q = (/ g(1)*hv*hv, g(2)*hu*hv, g(3)*hu*hu /)
      pq = (/ p(1)*hv*hv, p(2)*hu*hv, p(3)*hu*hu /)
      kink_corner = dual_sums(q, pq, mid)
      q(2) = -q(2)
      pq(2) = -pq(2)
      kink_corner = kink_corner + dual_sums(q, pq, mid)
      kink_corner = hu*hv*kink_corner/(8*pi*pi)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief EXP(-i*t) at the steps t of an angle from the source point.
!>
!>  @param[in] n  Number of nodes of the angle.
!>  @param[in] h  Grid step of the angle.
!>  @param[in] k0 Index of the source point, from zero.
!>  @returns EXP(-i*h*(k - k0)) for the nodes k = 0, ..., n - 1.
!-------------------------------------------------------------------------------
      PURE FUNCTION phase_table(n, h, k0)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(IN)             :: n
      COMPLEX(dp), DIMENSION(0:n - 1) :: phase_table
      REAL(dp), INTENT(IN)            :: h
      INTEGER, INTENT(IN)             :: k0

!  Local Variables
      COMPLEX(dp)                     :: wstep
      INTEGER                         :: k

!  Start of executable code
      wstep = CMPLX(COS(h), -SIN(h), dp)
      phase_table(k0) = 1.0_dp
      DO k = 1, n - 1
         phase_table(MOD(k0 + k, n)) = phase_table(MOD(k0 + k - 1, n))*wstep
      END DO

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Set the quantities of the values about the kinks for a source point.
!>
!>  @param[in]    ip Index of the source point.
!>  @param[inout] kc Quantities of the values about the kinks.
!-------------------------------------------------------------------------------
      SUBROUTINE kink_setup(ip, kc)
      USE vacmod, ONLY: nu, nv, nuv, nvper, alu, alv, pi2, tanu, guu_b,        &
                        guv_b, gvv_b, auu, auv, avv

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(IN)                :: ip
      TYPE (kink_context), INTENT(INOUT) :: kc

!  Local Variables
      INTEGER                            :: nmax
      INTEGER                            :: i

!  Start of executable code
      kc%ku0 = (ip - 1)/nv
      kc%kv0 = MOD(ip - 1, nv)
      kc%iuoff = nuv - nv*kc%ku0
      kc%ivoff = nuv + 1 - ip
      kc%kua = MOD(kc%ku0 + nu/2, nu)*nv
      kc%kfu = kc%ku0 + nu/2 - 1
      kc%kfv = kc%kv0 + (nv - 2)/2

!  A line of v has the coefficients in the order (uu, uv, vv) and a line of u
!  in the order (vv, uv, uu)
      kc%g(:,1) = (/ guu_b(ip), guv_b(ip), gvv_b(ip) /)
      kc%p(:,1) = (/ auu(ip), auv(ip), avv(ip) /)
      kc%g(:,2) = kc%g(3:1:-1,1)
      kc%p(:,2) = kc%p(3:1:-1,1)

      nmax = MAX(nu, nv, nvper)
      ALLOCATE (kc%w(0:nmax - 1,2), kc%rows(nu,2), kc%wrow(nu,2),              &
                kc%skept(0:nmax/2,2), kc%m0(2,0:nmax/2,2),                     &
                kc%m1(2,0:nmax/2,2))
      kc%skept = 0.0_dp

      kc%h(2) = alu
      kc%mid(2) = .FALSE.
      kc%st(:,2) = kink_steps(alu, nu, .FALSE., guu_b(ip))
      kc%w(0:nu - 1,2) = phase_table(nu, alu, kc%ku0)
      IF (nv .GT. 1) THEN
         kc%h(1) = alv
         kc%mid(1) = MOD(nv, 2) .EQ. 1
         kc%st(:,1) = kink_steps(alv, nv, kc%mid(1), gvv_b(ip))
         kc%w(0:nv - 1,1) = phase_table(nv, alv, kc%kv0)
      ELSE
!  The lines of v are the sums over the images, which take the weight 1 alone
         kc%h(1) = pi2/nvper
         kc%mid(1) = .FALSE.
         kc%st(:,1) = kink_steps(kc%h(1), nvper, .FALSE., gvv_b(ip))
         kc%rows = 0.0_dp
         DO i = 1, nu
            kc%wrow(i,:) = kink_weights(tanu(i + kc%iuoff), kc%g(:,1), nvper,  &
                                        kc%h(1))
         END DO
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Integrals of the approximations along a line of the grid, kept for
!>  the line at the opposite step from the source.
!>
!>  The integrals of the two lines are equal for the weight 1 and conjugate
!>  for the weight EXP(-i*t).
!>
!>  @param[inout] kc  Quantities of the values about the kinks.
!>  @param[in]    fam 1 for a line of v and 2 for a line of u.
!>  @param[in]    d   Number of steps of the line from the source.
!>  @param[in]    s   Tangent of the other angle on the line.
!>  @param[out]   m0  Integrals of ga1 and of ga1*ga2 along the line.
!>  @param[out]   m1  Integrals of ga1 and of ga1*ga2 times EXP(-i*t).
!-------------------------------------------------------------------------------
      PURE SUBROUTINE kink_pair(kc, fam, d, s, m0, m1)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (kink_context), INTENT(INOUT)     :: kc
      INTEGER, INTENT(IN)                    :: fam
      INTEGER, INTENT(IN)                    :: d
      REAL(dp), INTENT(IN)                   :: s
      REAL(dp), DIMENSION(2), INTENT(OUT)    :: m0
      COMPLEX(dp), DIMENSION(2), INTENT(OUT) :: m1

!  Start of executable code
      IF (kc%skept(d,fam) .EQ. 0) THEN
         CALL kink_moments(s, kc%g(:,fam), kc%p(:,fam), kc%m0(:,d,fam),        &
                           kc%m1(:,d,fam))
         kc%skept(d,fam) = s
      END IF

      m0 = kc%m0(:,d,fam)
      IF (s*kc%skept(d,fam) .GT. 0) THEN
         m1 = kc%m1(:,d,fam)
      ELSE
         m1 = CONJG(kc%m1(:,d,fam))
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Values of the approximations on the nodes about the kink of a line
!>  of the grid.
!>
!>  @param[inout] kc  Quantities of the values about the kinks.
!>  @param[in]    fam 1 for a line of v and 2 for a line of u.
!>  @param[in]    d   Number of steps of the line from the source.
!>  @param[in]    s   Tangent of the other angle on the line.
!>  @param[in]    f1  ga1 on the nodes of the line.
!>  @param[in]    f2  Ratio of ga1*ga2 to ga1 on the nodes of the line.
!>  @param[out]   v   Values of ga1 and of ga1*ga2 for the nodes about the
!>                    kink, in the order of the angle.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE line_values(kc, fam, d, s, f1, f2, v)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (kink_context), INTENT(INOUT)     :: kc
      INTEGER, INTENT(IN)                    :: fam
      INTEGER, INTENT(IN)                    :: d
      REAL(dp), INTENT(IN)                   :: s
      REAL(dp), DIMENSION(:), INTENT(IN)     :: f1
      REAL(dp), DIMENSION(:), INTENT(IN)     :: f2
      REAL(dp), DIMENSION(4,2), INTENT(OUT)  :: v

!  Local Variables
      REAL(dp), DIMENSION(2)                 :: w
      REAL(dp), DIMENSION(2)                 :: s0
      REAL(dp), DIMENSION(2)                 :: m0
      REAL(dp), DIMENSION(2)                 :: x0
      COMPLEX(dp), DIMENSION(2)              :: s1
      COMPLEX(dp), DIMENSION(2)              :: m1
      COMPLEX(dp), DIMENSION(2)              :: x1
      REAL(dp)                               :: fk
      INTEGER                                :: n
      INTEGER                                :: k

!  Start of executable code
      n = SIZE(f1)
      w = kink_weights(s, kc%g(:,fam), n, kc%h(fam))
      x0 = 0.0_dp
      x1 = 0.0_dp
      IF (w(1) .GT. 0) THEN
         s0 = 0.0_dp
         s1 = 0.0_dp
         DO k = 1, n
            fk = f1(k)*f2(k)
            s0(1) = s0(1) + f1(k)
            s0(2) = s0(2) + fk
            s1(1) = s1(1) + f1(k)*kc%w(k - 1,fam)
            s1(2) = s1(2) + fk*kc%w(k - 1,fam)
         END DO
         CALL kink_pair(kc, fam, d, s, m0, m1)
         x0 = kc%h(fam)*s0 - m0
         x1 = kc%h(fam)*s1 - m1
      END IF

      v = kink_line(s, kc%g(:,fam), kc%p(:,fam), kc%h(fam), kc%mid(fam),       &
                    kc%st(:,fam), w, x0, x1)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Values about the kinks on the grid of one field period.
!>
!>  The nodes about the kink of each line of v and of each line of u take the
!>  values of @ref kink_line, and the node where the kinks cross takes that of
!>  @ref kink_corner.
!>
!>  @param[in]    ga1    First analytic approximation on the grid.
!>  @param[in]    ga2    Ratio of the second analytic approximation to the
!>                       first.
!>  @param[inout] kc     Quantities of the values about the kinks.
!>  @param[inout] delgr  Difference of the Green's function and ga1.
!>  @param[inout] delgrp Difference of its normal derivative and ga1*ga2.
!-------------------------------------------------------------------------------
      SUBROUTINE kink_grid(ga1, ga2, kc, delgr, delgrp)
      USE vacmod, ONLY: nu, nv, nuv, tanu, tanv

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(:), INTENT(IN)    :: ga1
      REAL(dp), DIMENSION(:), INTENT(IN)    :: ga2
      TYPE (kink_context), INTENT(INOUT)    :: kc
      REAL(dp), DIMENSION(:), INTENT(INOUT) :: delgr
      REAL(dp), DIMENSION(:), INTENT(INOUT) :: delgrp

!  Local Variables
      REAL(dp), DIMENSION(4,2)              :: v
      REAL(dp), DIMENSION(2)                :: c
      INTEGER                               :: ku
      INTEGER                               :: kv
      INTEGER                               :: i0
      INTEGER                               :: i
      INTEGER                               :: k
      INTEGER                               :: d

!  Start of executable code
!  Lines of v, but that on the kink of u
      DO ku = 0, nu - 1
         i0 = ku*nv
         IF (i0 .EQ. kc%kua) THEN
            CYCLE
         END IF
         d = MOD(ku - kc%ku0 + nu, nu)
         CALL line_values(kc, 1, MIN(d, nu - d), tanu(i0 + 1 + kc%iuoff),      &
                          ga1(i0 + 1:i0 + nv), ga2(i0 + 1:i0 + nv), v)
         DO k = 1, 4
            i = i0 + 1 + MOD(kc%kfv + k - 1, nv)
            delgr(i) = delgr(i) - v(k,1)
            delgrp(i) = delgrp(i) - v(k,2)
         END DO
      END DO

!  Lines of u, but that on the kink of v
      DO kv = 1, nv
         IF (.NOT.kc%mid(1) .AND. kv .EQ. 1 + MOD(kc%kfv + 1, nv)) THEN
            CYCLE
         END IF
         d = MOD(kv - 1 - kc%kv0 + nv, nv)
         CALL line_values(kc, 2, MIN(d, nv - d), tanv(kv + kc%ivoff),          &
                          ga1(kv:nuv:nv), ga2(kv:nuv:nv), v)
         DO k = 1, 3
            i = MOD(kc%kfu + k - 1, nu)*nv + kv
            delgr(i) = delgr(i) - v(k,1)
            delgrp(i) = delgrp(i) - v(k,2)
         END DO
      END DO

!  Crossing of the kinks
      c = kink_corner(kc%g(:,1), kc%p(:,1), kc%h(2), kc%h(1), kc%mid(1))
      IF (kc%mid(1)) THEN
         DO k = 1, 2
            i = kc%kua + 1 + MOD(kc%kfv + k, nv)
            delgr(i) = delgr(i) - c(1)/2
            delgrp(i) = delgrp(i) - c(2)/2
         END DO
      ELSE
         i = kc%kua + 1 + MOD(kc%kfv + 1, nv)
         delgr(i) = delgr(i) - c(1)
         delgrp(i) = delgrp(i) - c(2)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Values about the kink of u on one image, for nv = 1.
!>
!>  The line of u of the image takes the values of @ref kink_line, unless the
!>  image lies on the kink of v, and the sums over the images of the lines of
!>  v that use them take the terms of the image.
!>
!>  @param[in]    kp     Index of the image.
!>  @param[in]    ga1    First analytic approximation on the image.
!>  @param[in]    ga2    Ratio of the second analytic approximation to the
!>                       first.
!>  @param[inout] kc     Quantities of the values about the kinks.
!>  @param[inout] delgr  Difference of the Green's function and ga1.
!>  @param[inout] delgrp Difference of its normal derivative and ga1*ga2.
!-------------------------------------------------------------------------------
      SUBROUTINE kink_image(kp, ga1, ga2, kc, delgr, delgrp)
      USE vacmod, ONLY: nu, nvper, tanv

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(IN)                   :: kp
      REAL(dp), DIMENSION(:), INTENT(IN)    :: ga1
      REAL(dp), DIMENSION(:), INTENT(IN)    :: ga2
      TYPE (kink_context), INTENT(INOUT)    :: kc
      REAL(dp), DIMENSION(:), INTENT(INOUT) :: delgr
      REAL(dp), DIMENSION(:), INTENT(INOUT) :: delgrp

!  Local Variables
      REAL(dp), DIMENSION(4,2)              :: v
      INTEGER                               :: i
      INTEGER                               :: k

!  Start of executable code
      IF (2*(kp - 1) .NE. nvper) THEN
         DO i = 1, nu
            IF (kc%wrow(i,1) .GT. 0) THEN
               kc%rows(i,1) = kc%rows(i,1) + ga1(i)
               kc%rows(i,2) = kc%rows(i,2) + ga1(i)*ga2(i)
            END IF
         END DO

         CALL line_values(kc, 2, MIN(kp - 1, nvper - kp + 1),                  &
                          tanv(1 + kc%ivoff + 2*nu*(kp - 1)), ga1, ga2, v)
         DO k = 1, 3
            i = 1 + MOD(kc%kfu + k - 1, nu)
            delgr(i) = delgr(i) - v(k,1)
            delgrp(i) = delgrp(i) - v(k,2)
         END DO
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Values about the kink of v and at the crossing, for nv = 1.
!>
!>  The lines of v are the sums over the images, for which the modes do not
!>  vary, so each node of u takes the sum of the values of @ref kink_line, and
!>  the node on the kink of u takes the value of @ref kink_corner.
!>
!>  @param[inout] kc     Quantities of the values about the kinks.
!>  @param[inout] delgr  Difference of the Green's function and ga1.
!>  @param[inout] delgrp Difference of its normal derivative and ga1*ga2.
!-------------------------------------------------------------------------------
      SUBROUTINE kink_image_sums(kc, delgr, delgrp)
      USE vacmod, ONLY: nu, tanu

      IMPLICIT NONE

!  Declare Arguments
      TYPE (kink_context), INTENT(INOUT)    :: kc
      REAL(dp), DIMENSION(:), INTENT(INOUT) :: delgr
      REAL(dp), DIMENSION(:), INTENT(INOUT) :: delgrp

!  Local Variables
      REAL(dp), DIMENSION(4,2)              :: v
      REAL(dp), DIMENSION(2)                :: c
      REAL(dp), DIMENSION(2)                :: m0
      REAL(dp), DIMENSION(2)                :: x0
      COMPLEX(dp), DIMENSION(2)             :: m1
      COMPLEX(dp), DIMENSION(2)             :: x1
      REAL(dp)                              :: s
      INTEGER                               :: i
      INTEGER                               :: d

!  Start of executable code
      x1 = 0.0_dp
      DO i = 1, nu
         IF (i .EQ. kc%kua + 1) THEN
            CYCLE
         END IF
         s = tanu(i + kc%iuoff)
         x0 = 0.0_dp
         IF (kc%wrow(i,1) .GT. 0) THEN
            d = MOD(i - 1 - kc%ku0 + nu, nu)
            CALL kink_pair(kc, 1, MIN(d, nu - d), s, m0, m1)
            x0 = kc%h(1)*kc%rows(i,:) - m0
         END IF
         v = kink_line(s, kc%g(:,1), kc%p(:,1), kc%h(1), .FALSE., kc%st(:,1),  &
                       kc%wrow(i,:), x0, x1)
         delgr(i) = delgr(i) - SUM(v(:,1))
         delgrp(i) = delgrp(i) - SUM(v(:,2))
      END DO

      c = kink_corner(kc%g(:,1), kc%p(:,1), kc%h(2), kc%h(1), .FALSE.)
      delgr(kc%kua + 1) = delgr(kc%kua + 1) - c(1)
      delgrp(kc%kua + 1) = delgrp(kc%kua + 1) - c(2)

      END SUBROUTINE

      END MODULE
