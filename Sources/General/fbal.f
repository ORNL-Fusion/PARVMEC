!> \file fbal.f

      MODULE fbal
      USE stel_kinds, ONLY: dp
      REAL(dp), DIMENSION(:), ALLOCATABLE :: rzu_fac, rru_fac,
     1  frcc_fac, fzsc_fac

      CONTAINS

      SUBROUTINE calc_fbal_par(bsubu, bsubv)
      USE vmec_main, ONLY: buco, bvco, equif, iequi,
     1                     jcurv, jcuru, chipf, vp, pres, 
     2                     phipf, vpphi, presgrad, ohs
      USE vmec_params, ONLY: signgs
      USE vmec_dim, ONLY: ns, nrzt, nznt, ns1
      USE realspace, ONLY: pwint, phip
      USE vmec_input, ONLY: nzeta
      USE vmec_dim, ONLY: ntheta3
      USE parallel_include_module 
      IMPLICIT NONE
!-----------------------------------------------
      REAL(dp), INTENT(in) :: bsubu(nznt,ns),
     1                        bsubv(nznt,ns)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: js, lk
      INTEGER :: nsmin, nsmax
!-----------------------------------------------

      nsmin=t1lglob; nsmax=t1rglob
      DO js = nsmin, nsmax
        buco(js) = SUM(bsubu(:,js)*pwint(:,js))
        bvco(js) = SUM(bsubv(:,js)*pwint(:,js))
      END DO

      CALL Gather1XArray(bvco)
      CALL Gather1XArray(buco)

!     FROM AMPERE'S LAW, JcurX are angle averages of jac*JsupX, so
!                        JcurX = (dV/ds)/twopi**2 <JsupX> where <...> is flux surface average
      !nsmin=MAX(2,t1lglob); nsmax=MIN(t1rglob,ns-1)
      nsmin=MAX(2,tlglob); nsmax=MIN(trglob,ns-1)
      DO js = nsmin, nsmax
         jcurv(js) = (signgs*ohs)*(buco(js+1) - buco(js))
         jcuru(js) =-(signgs*ohs)*(bvco(js+1) - bvco(js))
         vpphi(js) = (vp(js+1) + vp(js))/2
         presgrad(js) = (pres(js+1) - pres(js))*ohs
         equif(js) = (-phipf(js)*jcuru(js) + chipf(js)*jcurv(js))
     1                /vpphi(js) + presgrad(js)
      END DO
      equif(1) = 0
      equif(ns) = 0

      !SKS-RANGE: All LHS's computed correctly in [t1lglob, trglob]

      END SUBROUTINE calc_fbal_par

      END MODULE fbal
