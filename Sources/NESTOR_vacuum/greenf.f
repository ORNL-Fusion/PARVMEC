      SUBROUTINE greenf (delgr, delgrp, ip)
      USE vacmod
      USE vparams, ONLY: one
      USE parallel_include_module
      USE timer_sub
      USE greenf_singular, ONLY: singular_defects
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ip
      REAL(dp), DIMENSION(nuv), INTENT(OUT) :: delgr, delgrp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(2) :: ilow, ihigh
      INTEGER :: ivoff, iskip, iuoff, i, kp, nloop
      INTEGER :: kua, kva, kv0, kvm, iup, ium, ivp, ivm
      REAL(dp), DIMENSION(:), ALLOCATABLE ::
     1    ftemp, gsave, htemp, ga1, ga2, dsave
      REAL(dp):: xip, yip, xper, yper,
     1    sxsave, sysave, tgreenon, tgreenoff
      REAL(dp) :: hv, su, du, sv, dv
      REAL(dp) :: wg(3), wk(3), dphi, area
C-----------------------------------------------
!
!     ON ENTRANCE, IP IS THE INDEX OF THE PRIMED MESH POINT (lies in 1st field period)
!
!     ON EXIT, DELGR IS THE DIFFERENCE OF "GREEN'S FUNCTION"
!     AND ANALYTIC APPROXIMATION, SUMMED OVER ALL FIELD PERIODS
!     DELGRP IS DIFFERENCE OF DERIVATIVE OF "GREEN'S FUNCTION"
!     AND ANALYTIC APPROXIMATION.
!
!     BOTH THESE QUANTITIES ARE COMPUTED FOR ALL UNPRIMED U,V POINTS IN ONE FIELD PERIOD,
!     FOR THIS FIXED PRIMED POINT (IP).
!
      CALL second0(tgreenon)

      ALLOCATE (ftemp(nuv), gsave(nuv), htemp(nuv), ga1(nuv), ga2(nuv),
     1          dsave(nuv), stat=i)
      IF (i .NE. 0) STOP 'allocation error in greenf'

!
!     COMPUTE OFFSETS FOR U,V ANGLE DIFFERENCES AND CONSTANTS
!
      ilow(1) = 1
      ilow(2) = ip + 1
      ihigh(1) = ip - 1
      ihigh(2) = nuv
      ivoff = nuv + 1 - ip
      iskip = (ip - 1)/nv
      iuoff = nuv - nv*iskip
      xip = rcosuv(ip)             !x == r*COS(ip), in 1st field period
      yip = rsinuv(ip)             !y == r*SIN(ip), in 1st field period
      delgr  = 0
      delgrp = 0

!
!     COMPUTE FIELD-PERIOD INVARIANT VECTORS
!
!     NOTE: |x - x'|**2 = gsave - 2*[x*x' + y*y']
!
      DO i = 1, nuv
        gsave(i) = rzb2(ip) + rzb2(i) - 2*z1b(ip)*z1b(i)
        dsave(i) = drv(ip) + z1b(i)*snz(ip)
      END DO

!
!     SUM OVER FIELD-PERIODS (NVPER=NFPER) OR INTEGRATE OVER NV (NVPER=64) IF NV == 1
!
!     NOTE THE SURFACE NORMAL SNORM == Xu cross Xv = NP*[SNR, SNV, SNZ]
!     IS PERIODIC ON EACH FIELD PERIOD
!
      DO kp = 1, nvper
         xper = xip*cosper(kp) - yip*sinper(kp)              !x(ip) in field period kp
         yper = yip*cosper(kp) + xip*sinper(kp)              !y(ip) in field period kp
         sxsave = (snr(ip)*xper - snv(ip)*yper)/r1b(ip)
         sysave = (snr(ip)*yper + snv(ip)*xper)/r1b(ip)

         IF (kp.EQ.1 .OR. nv.EQ.1) THEN

!        INITIALIZE ANALYTIC APPROXIMATIONS GA1, GA2
            DO i = 1, nuv
               ga1(i) = tanu(i+iuoff)*(guu_b(ip)*tanu(i+iuoff) 
     1                               + guv_b(ip)*tanv(i+ivoff))
     2                + gvv_b(ip)*tanv(i+ivoff)*tanv(i+ivoff)
               ga2(i) = tanu(i+iuoff)*(auu(ip)*tanu(i+iuoff)
     1                + auv(ip)*tanv(i+ivoff))
     2                + avv(ip)*tanv(i+ivoff)*tanv(i+ivoff)
            END DO

            DO nloop = 1, 2
               IF (kp.GT.1 .AND. nloop.EQ.2) CYCLE 
               DO i = ilow(nloop), ihigh(nloop)
                 ga2(i) = ga2(i)/ga1(i)
                 ga1(i) = one/SQRT(ga1(i))
                 ftemp(i) = one/(gsave(i) 
     1                    - 2*(xper*rcosuv(i) + yper*rsinuv(i)))
                 htemp(i) = SQRT(ftemp(i))
                 delgrp(i) = delgrp(i) - ga2(i)*ga1(i) 
     1                     + ftemp(i)*htemp(i)*
     2              (rcosuv(i)*sxsave + rsinuv(i)*sysave + dsave(i))
                 delgr(i) = delgr(i) + htemp(i) - ga1(i)
               END DO
            END DO

!
!           GA1 AND GA1*GA2 VANISH WITH A KINK WHERE U - U' OR V - V' IS PI,
!           WHERE THE GREEN'S FUNCTION IS SMOOTH. THE TRAPEZOIDAL SUM ACROSS
!           A KINK OF SLOPE S IS OFF BY H**2*S/6, SO THE NODES ON THE KINK
!           TAKE H*S/6 IN PLACE OF ZERO, AND THE NODE WHERE THE U AND V KINKS
!           CROSS TAKES THE V VALUE. FOR ODD NV THE V KINK LIES MIDWAY
!           BETWEEN TWO NODES, AND EACH OF THEM TAKES H*S/24 LESS.
!
            hv = alv
            IF (nv .eq. 1) THEN
               hv = pi2/nvper
            END IF
            su = alu/(24*SQRT(guu_b(ip)))
            du = su*auu(ip)/guu_b(ip)
            sv = hv/(24*SQRT(gvv_b(ip)))
            dv = sv*avv(ip)/gvv_b(ip)
            kua = MOD(iskip + nu/2, nu)*nv
            IF (nv .eq. 1) THEN
               IF (2*(kp - 1) .eq. nvper) THEN
                  delgr = delgr - sv
                  delgrp = delgrp - dv
               ELSE
                  delgr(kua + 1) = delgr(kua + 1) - su
                  delgrp(kua + 1) = delgrp(kua + 1) - du
               END IF
            ELSE
               delgr(kua + 1:kua + nv) = delgr(kua + 1:kua + nv) - su
               delgrp(kua + 1:kua + nv) = delgrp(kua + 1:kua + nv) - du
               kv0 = MOD(ip - 1, nv)
               IF (MOD(nv, 2) .eq. 0) THEN
                  kva = MOD(kv0 + nv/2, nv) + 1
                  delgr(kua + kva) = delgr(kua + kva) + su
                  delgrp(kua + kva) = delgrp(kua + kva) + du
                  delgr(kva:nuv:nv) = delgr(kva:nuv:nv) - sv
                  delgrp(kva:nuv:nv) = delgrp(kva:nuv:nv) - dv
               ELSE
                  DO kvm = kv0 + (nv - 1)/2, kv0 + (nv + 1)/2
                     kva = MOD(kvm, nv) + 1
                     delgr(kva:nuv:nv) = delgr(kva:nuv:nv) + sv/4
                     delgrp(kva:nuv:nv) = delgrp(kva:nuv:nv) + dv/4
                     delgr(kua + kva) = delgr(kua + kva) - sv/4
                     delgrp(kua + kva) = delgrp(kua + kva) - dv/4
                  END DO
               END IF
            END IF

!
!           THE SUM SKIPS THE SOURCE POINT, WHERE THE REMAINDERS TEND TO
!           A DIRECTION-DEPENDENT LIMIT AND RISE LINEARLY; ITS DEFECTS FOR
!           THE WEIGHTS 1, U - U' AND V - V' GO ONTO THE SOURCE NODE AND
!           ITS NEIGHBOURS. FOR NV = 1 THE V STEP IS THE IMAGE SPACING IN
!           THE TOROIDAL ANGLE, ALONG WHICH THE MODES DO NOT VARY.
!
            IF (kp .EQ. 1) THEN
               kv0 = MOD(ip - 1, nv)
               dphi = onp
               IF (nv .EQ. 1) THEN
                  dphi = one
               END IF
               CALL singular_defects(rtay(:,:,ip), ztay(:,:,ip), xip,
     &                               yip, z1b(ip), alvp*kv0, dphi,
     &                               sxsave, sysave, snz(ip),
     &                               guu_b(ip), guv_b(ip), gvv_b(ip),
     &                               auu(ip), auv(ip), avv(ip), alu,
     &                               hv, wg, wk)
               area = alu*hv
               iup = ip + nv
               IF (iup .GT. nuv) THEN
                  iup = iup - nuv
               END IF
               ium = ip - nv
               IF (ium .LT. 1) THEN
                  ium = ium + nuv
               END IF
               delgr(ip) = delgr(ip) + wg(1)/area
               delgrp(ip) = delgrp(ip) + wk(1)/area
               delgr(iup) = delgr(iup) + wg(2)/(2*alu*area)
               delgr(ium) = delgr(ium) - wg(2)/(2*alu*area)
               delgrp(iup) = delgrp(iup) + wk(2)/(2*alu*area)
               delgrp(ium) = delgrp(ium) - wk(2)/(2*alu*area)
               IF (nv .GT. 1) THEN
                  ivp = ip + 1
                  IF (kv0 .EQ. nv - 1) THEN
                     ivp = ivp - nv
                  END IF
                  ivm = ip - 1
                  IF (kv0 .EQ. 0) THEN
                     ivm = ivm + nv
                  END IF
                  delgr(ivp) = delgr(ivp) + wg(3)/(2*hv*area)
                  delgr(ivm) = delgr(ivm) - wg(3)/(2*hv*area)
                  delgrp(ivp) = delgrp(ivp) + wk(3)/(2*hv*area)
                  delgrp(ivm) = delgrp(ivm) - wk(3)/(2*hv*area)
               END IF
            END IF

            IF (kp.EQ.nvper .AND. nv.EQ.1) THEN
                delgrp = delgrp/nvper
                delgr  = delgr /nvper
            END IF

            ivoff = ivoff + 2*nu
            ihigh(1) = nuv

         ELSE
            DO i = 1, nuv
              ftemp(i) = one/(gsave(i)  
     1                 -  2*(xper*rcosuv(i) + yper*rsinuv(i)))
              htemp(i) = SQRT(ftemp(i))
              delgrp(i) = delgrp(i) + ftemp(i)*htemp(i)*
     1           (rcosuv(i)*sxsave + rsinuv(i)*sysave + dsave(i))
              delgr(i) = delgr(i) + htemp(i)
           END DO
         ENDIF
      END DO

      DEALLOCATE (ftemp, gsave, htemp, ga1, ga2, dsave, stat=i)

      CALL second0(tgreenoff)
      timer_vac(tgreenf) = timer_vac(tgreenf) + (tgreenoff-tgreenon)
      greenf_time = timer_vac(tgreenf)

      END SUBROUTINE greenf

