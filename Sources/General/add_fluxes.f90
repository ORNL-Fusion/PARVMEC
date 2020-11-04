!> \file add_fluxes.f90

      SUBROUTINE add_fluxes_par(overg, bsupu, bsupv, lcurrent)
      USE vmec_main
      USE realspace, ONLY: pwint, pguu, pguv, pchip, pphip
      USE vmec_input, ONLY: nzeta
      USE vmec_dim, ONLY: ntheta3
      USE parallel_include_module
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(nznt,ns), INTENT(in)    :: overg
      REAL(rprec), DIMENSION(nznt,ns), INTENT(inout) :: bsupu, bsupv
      LOGICAL, INTENT(in) :: lcurrent
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: p5=0.5_dp, c1p5=1.5_dp
      REAL(rprec), PARAMETER :: iotaped = 0.10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: js, l
      REAL(rprec) :: top, bot

      INTEGER :: i, j, k, nsmin, nsmax, lnsnum, istat
!-----------------------------------------------
!
!     ADD MAGNETIC FLUX (CHIP, PHIP) TERMS TO BSUPU=-OVERG*LAM_V, BSUPV=OVERG*LAM_U
!     COMPUTE FLUX FROM ITOR = <B_u>, ITOR(s) = integrated toroidal current (icurv)
!     IF ncurr == 1, AND ictrl_prec2d != 0, COMPUTE FORCE IN TOMNSP TO UPDATE chips
!

      IF (.NOT.lcurrent .OR. ncurr.EQ.0) GOTO 100

      nsmin=MAX(2,t1lglob); nsmax=t1rglob
      DO js = nsmin, nsmax
        top = icurv(js)
        bot = 0
        DO j=1,nznt
          top = top - pwint(j,js)*(pguu(j,js)*bsupu(j,js) &
            + pguv(j,js)*bsupv(j,js))
          bot = bot + pwint(j,js)*overg(j,js)*pguu(j,js)
        END DO
        IF (bot.ne.zero) chips(js) = top/bot
        IF (phips(js).ne.zero)  iotas(js) = chips(js)/phips(js)
      END DO

 100  CONTINUE

      nsmin=MAX(2,t1lglob); nsmax=t1rglob
!     CHANGE THIS FOR lRFP = T  (solve for phips?)
      IF (ncurr .EQ. 0) THEN
         chips(nsmin:nsmax) = iotas(nsmin:nsmax)*phips(nsmin:nsmax)
      ELSE IF (.NOT.lcurrent) THEN
         WHERE (phips(nsmin:nsmax) .NE. zero) &
           iotas(nsmin:nsmax) = chips(nsmin:nsmax)/phips(nsmin:nsmax)
      END IF

      DO js = nsmin, nsmax
         pchip(:,js) = chips(js)
      END DO

      nsmin=MAX(2,t1lglob); nsmax=MIN(ns-1,trglob)
      IF (t1lglob .eq. 1 .and. trglob .gt. 2) THEN
         chipf(1) = c1p5*chips(2) - p5*chips(3)
      ELSE IF (t1lglob .eq. 1) THEN
         chipf(1) = chips(2)
      END IF
      chipf(nsmin:nsmax) = (chips(nsmin:nsmax) + chips(nsmin+1:nsmax+1))/2
      IF (nsmax.EQ.ns) chipf(ns)    = c1p5*chips(ns)- p5*chips(ns-1)

!     Do not compute iota too near origin
      IF(trglob_arr(1).LE.2) THEN
#if defined(MPI_OPT)
        CALL MPI_Bcast(iotas(3),1,MPI_REAL8,1,NS_COMM,MPI_ERR)
#endif
      END IF
      IF (lrfp) THEN
         IF (nsmin.EQ.2) iotaf(1)  = one/(c1p5/iotas(2) - p5/iotas(3))
         IF (nsmax.EQ.ns-1) iotaf(ns)=one/(c1p5/iotas(ns)-p5/iotas(ns-1))
         DO js = MAX(2,t1lglob), MIN(ns-1,t1rglob)
            iotaf(js) = 2.0_dp/(one/iotas(js) + one/iotas(js+1))
         END DO
      ELSE
        IF (nsmin.EQ.2) iotaf(1)  = c1p5*iotas(2) - p5*iotas(3)
        IF (nsmax.EQ.ns-1) iotaf(ns)=c1p5*iotas(ns) - p5*iotas(ns-1)
        DO js = MAX(2,t1lglob), MIN(ns-1,trglob)
          iotaf(js) = p5*(iotas(js) + iotas(js+1))
        END DO
      END IF

      nsmin=MAX(1,t1lglob); nsmax=MIN(ns,t1rglob)
      bsupu(:,nsmin:nsmax) = bsupu(:,nsmin:nsmax)+pchip(:,nsmin:nsmax)*overg(:,nsmin:nsmax)

      END SUBROUTINE add_fluxes_par
