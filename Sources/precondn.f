!> \file precondn.f

      ! SKS-RANGE: All OUT arrays computed correctly between
      ! [tlglob, trglob] 
      SUBROUTINE precondn_par(lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo,
     &                    xodd, axm, axd, bxm, bxd,
     &                    cx, eqfactor, trigmult)
      USE vmec_main
      USE vmec_params, ONLY: signgs
      USE realspace, ONLY: pshalf, pwint
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nznt,ns), INTENT(in) ::
     &   lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo, xodd
      REAL(rprec), DIMENSION(ns+1,2), INTENT(out) ::
     &   axm, axd, bxm, bxd
      REAL(rprec), DIMENSION(ns+1), INTENT(out) :: cx
      REAL(rprec), DIMENSION(ns), INTENT(out) :: eqfactor
      REAL(rprec), DIMENSION(nznt), INTENT(in) :: trigmult
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, l, lk
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ax, bx
      !REAL(rprec) :: temp(ns+1)
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: temp
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: ptau, ptau2
      REAL(rprec) :: t1, t2, t3, pfactor

      INTEGER :: nsmin, nsmax,  i, j, k, numjs
      INTEGER, DIMENSION(:), ALLOCATABLE :: ldisps, lcounts
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: sbuf

      nsmin=tlglob; nsmax=t1rglob
      ! Correct incoming lu1, bsq, gsqrt,   , xue, xuo, xodd, trigmult 
      ! between [t1lglob, t1rglob]

      ALLOCATE (ax(ns+1,4), bx(ns+1,4),
     &          ptau(nznt), ptau2(nznt), temp(ns+1))
      ax = 0
      bx = 0
      cx = 0
      temp = 0
      pfactor = -4*r0scale**2        !restored in v8.51

      nsmin = MAX(2,tlglob)
      nsmax = t1rglob
      DO js = nsmin, nsmax
!
!     COMPUTE DOMINANT (1/DELTA-S)**2 PRECONDITIONING
!     MATRIX ELEMENTS
!
         lk = 0
         DO k = 1, ntheta3
            DO j = 1, nzeta
               lk = lk + 1
               t1 = pfactor*r12(lk,js)*bsq(lk,js)
               ptau2(lk) = r12(lk,js)*t1/gsqrt(lk,js)
               t1 = t1*pwint(lk,js)
               temp(js) = temp(js) + t1*trigmult(lk)*xu12(lk,js)
               ptau(lk) = r12(lk,js)*t1/gsqrt(lk,js)
               t1 = xu12(lk,js)*ohs
               t2 = cp25*(xue(lk,js)/pshalf(lk,js)+xuo(lk,js))
     &            / pshalf(lk,js)
               t3 = cp25*(xue(lk,js-1)/pshalf(lk,js) +
     &                    xuo(lk,js-1))/pshalf(lk,js)
               ax(js,1) = ax(js,1) + ptau(lk)*t1*t1
               ax(js,2) = ax(js,2) + ptau(lk)*(-t1+t3)*(t1+t2)
               ax(js,3) = ax(js,3) + ptau(lk)*(t1+t2)*(t1+t2)
               ax(js,4) = ax(js,4) + ptau(lk)*(-t1+t3)*(-t1+t3)
            END DO
         END DO
!
!       COMPUTE PRECONDITIONING MATRIX ELEMENTS FOR M**2, N**2 TERMS
!
         lk = 0
         DO k = 1, ntheta3
            DO j = 1, nzeta
               lk = lk+1
               t1 = cp5*(xs(lk,js) + cp5*xodd(lk,js)/pshalf(lk,js))
               t2 = cp5*(xs(lk,js) + cp5*xodd(lk,js-1)/pshalf(lk,js))
               bx(js,1) = bx(js,1) + ptau(lk)*t1*t2
               bx(js,2) = bx(js,2) + ptau(lk)*t1*t1
               bx(js,3) = bx(js,3) + ptau(lk)*t2*t2
               cx(js) = cx(js)
     &                + cp25*pfactor*lu1(lk,js)**2 *
     &                  gsqrt(lk,js)*pwint(lk,js)
            END DO
         END DO

      END DO

      nsmin = MAX(2,tlglob)
      nsmax = t1rglob
      temp(1) = 0
      temp(nsmin:nsmax) = temp(nsmin:nsmax)/vp(nsmin:nsmax);
      temp(ns+1) = 0

      nsmin = t1lglob
      nsmax = t1rglob
      DO js = nsmin, nsmax
         axm(js,1) = -ax(js,1)
         axd(js,1) =  ax(js,1) + ax(js+1,1)
         axm(js,2) =  ax(js,2) * sm(js) * sp(js-1)
         axd(js,2) =  ax(js,3)*sm(js)**2 + ax(js+1,4)*sp(js)**2
         bxm(js,1) =  bx(js,1)
         bxm(js,2) =  bx(js,1) * sm(js) * sp(js-1)
         bxd(js,1) =  bx(js,2) + bx(js+1,3)
         bxd(js,2) =  bx(js,2)*sm(js)**2 + bx(js+1,3)*sp(js)**2
         cx(js)    =  cx(js) + cx(js+1)
         temp(js)  =  signgs*(temp(js) + temp(js+1))
      END DO

      nsmin = MAX(2,tlglob)
      nsmax = MIN(ns-1,trglob)
      eqfactor(nsmin:nsmax) = axd(nsmin:nsmax,2)*hs*hs/
     &                        temp(nsmin:nsmax)
      eqfactor(1) = 0
      eqfactor(ns) = 0
      axm(ns+1,:) = 0
      axd(ns+1,:) = 0
      bxm(ns+1,:) = 0
      bxd(ns+1,:) = 0
      DEALLOCATE (ax, bx, ptau, ptau2, temp)

      END SUBROUTINE precondn_par
