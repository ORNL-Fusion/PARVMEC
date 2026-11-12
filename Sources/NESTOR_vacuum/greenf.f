      SUBROUTINE greenf (delgr, delgrp, ip)
      USE vacmod
      USE vparams, ONLY: one
      USE parallel_include_module
      USE timer_sub
      USE greenf_kink, ONLY: kink_weights, kink_steps, kink_moments,
     &                       kink_line, kink_corner
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
      INTEGER :: kua, kv0, kfu, kfv, kr, kc, i0, k, iup, ium, ivp, ivm
      INTEGER :: nc
      REAL(dp), DIMENSION(:), ALLOCATABLE ::
     1    ftemp, gsave, htemp, ga1, ga2, dsave
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rows, cs
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: c0
      COMPLEX(dp), DIMENSION(:), ALLOCATABLE :: wu, wv
      COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: c1
      REAL(dp):: xip, yip, xper, yper,
     1    sxsave, sysave, tgreenon, tgreenoff
      REAL(dp) :: hv, gk, w0, w1, s0(2), m0(2), vg(4), vk(4), cg, ck,
     &            stu(7), stv(7)
      COMPLEX(dp) :: s1(2), m1(2), wstep
      LOGICAL :: lmid
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

      nc = MAX(nu, nv, nvper)/2
      ALLOCATE (ftemp(nuv), gsave(nuv), htemp(nuv), ga1(nuv), ga2(nuv),
     1          dsave(nuv), rows(nuv,4), wu(0:nu-1), wv(nv),
     2          cs(0:nc,2), c0(2,0:nc,2), c1(2,0:nc,2), stat=i)
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
!           GA1 AND GA1*GA2 VANISH WITH A KINK WHERE U - U' OR V - V' IS
!           PI, WHERE THE GREEN'S FUNCTION IS SMOOTH, SO THE TRAPEZOIDAL
!           SUMS ACROSS THE KINKS ARE IN ERROR. THE NODES ABOUT THE KINKS
!           TAKE THE VALUES OF GREENF_KINK IN PLACE OF ZERO: ON EACH LINE
!           OF THE GRID THEY MAKE THE SUMS OF GA1 AND GA1*GA2 TIMES 1, COS
!           AND SIN OF THE ANGLE EQUAL TO THE INTEGRALS, AND THE NODE
!           WHERE THE KINKS CROSS TAKES WHAT THE SUM OVER THE GRID STILL
!           MISSES. FOR NV = 1 THE LINES OF V ARE THE SUMS OVER THE IMAGES.
!
            hv = alv
            IF (nv .EQ. 1) THEN
               hv = pi2/nvper
            END IF
            kua = MOD(iskip + nu/2, nu)*nv
            IF (nv .EQ. 1) THEN
               CALL kink_images
            ELSE
               CALL kink_grid
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
               CALL singular_defects(rtay(:,:,ip), ztay(:,:,ip),
     &                               (/ xip, yip, z1b(ip) /),
     &                               alvp*kv0, dphi,
     &                               (/ sxsave, sysave, snz(ip) /),
     &                               (/ guu_b(ip), guv_b(ip),
     &                                  gvv_b(ip) /),
     &                               (/ auu(ip), auv(ip), avv(ip) /),
     &                               (/ alu, hv /), wg, wk)
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

      DEALLOCATE (ftemp, gsave, htemp, ga1, ga2, dsave, rows, wu, wv,
     &            cs, c0, c1, stat=i)

      CALL second0(tgreenoff)
      timer_vac(tgreenf) = timer_vac(tgreenf) + (tgreenoff-tgreenon)
      greenf_time = timer_vac(tgreenf)

      CONTAINS

!-----------------------------------------------------------------------
!     VALUES ABOUT THE KINKS ON THE GRID OF ONE FIELD PERIOD
!-----------------------------------------------------------------------
      SUBROUTINE kink_grid
      lmid = MOD(nv, 2) .EQ. 1
      kv0 = MOD(ip - 1, nv)
!     FIRST OF THE NODES ABOUT THE KINK ON A LINE OF V AND ON A LINE OF U
      kfv = kv0 + (nv - 2)/2
      kfu = iskip + nu/2 - 1
!     EXP(-i*t) AT THE STEPS OF V AND OF U FROM THE SOURCE
      wstep = CMPLX(COS(alv), -SIN(alv), dp)
      wv(1 + kv0) = 1
      DO k = 1, nv - 1
         wv(1 + MOD(kv0 + k, nv)) = wv(1 + MOD(kv0 + k - 1, nv))*wstep
      END DO
      wstep = CMPLX(COS(alu), -SIN(alu), dp)
      wu(iskip) = 1
      DO k = 1, nu - 1
         wu(MOD(iskip + k, nu)) = wu(MOD(iskip + k - 1, nu))*wstep
      END DO
      CALL kink_steps(hv, nv, lmid, gvv_b(ip), stv)
      CALL kink_steps(alu, nu, .FALSE., guu_b(ip), stu)
!     NO LINE HAS ITS INTEGRALS IN KINK_PAIR YET
      cs = 0

!     LINES OF V, BUT THAT ON THE KINK OF U
      DO kr = 0, nu - 1
         i0 = kr*nv
         IF (i0 .EQ. kua) THEN
            CYCLE
         END IF
         CALL kink_weights(tanu(i0 + 1 + iuoff), guu_b(ip), guv_b(ip),
     &                     gvv_b(ip), nv, hv, w0, w1)
         s0 = 0
         s1 = 0
         IF (w0 .GT. 0) THEN
            DO kc = 1, nv
               i = i0 + kc
               gk = ga1(i)*ga2(i)
               s0(1) = s0(1) + ga1(i)
               s0(2) = s0(2) + gk
               s1(1) = s1(1) + ga1(i)*wv(kc)
               s1(2) = s1(2) + gk*wv(kc)
            END DO
            k = MOD(kr - iskip + nu, nu)
            CALL kink_pair(1, MIN(k, nu - k), tanu(i0 + 1 + iuoff),
     &                     guu_b(ip), guv_b(ip), gvv_b(ip), auu(ip),
     &                     auv(ip), avv(ip))
         END IF
         CALL kink_line(tanu(i0 + 1 + iuoff), guu_b(ip), guv_b(ip),
     &                  gvv_b(ip), auu(ip), auv(ip), avv(ip), hv, lmid,
     &                  stv, w0, w1, s0, s1, m0, m1, vg, vk)
         DO k = 1, 4
            i = i0 + 1 + MOD(kfv + k - 1, nv)
            delgr(i) = delgr(i) - vg(k)
            delgrp(i) = delgrp(i) - vk(k)
         END DO
      END DO

!     LINES OF U, BUT THAT ON THE KINK OF V
      DO kc = 1, nv
         IF (.NOT.lmid .AND. kc .EQ. 1 + MOD(kfv + 1, nv)) THEN
            CYCLE
         END IF
         CALL kink_weights(tanv(kc + ivoff), gvv_b(ip), guv_b(ip),
     &                     guu_b(ip), nu, alu, w0, w1)
         s0 = 0
         s1 = 0
         IF (w0 .GT. 0) THEN
            DO kr = 0, nu - 1
               i = kr*nv + kc
               gk = ga1(i)*ga2(i)
               s0(1) = s0(1) + ga1(i)
               s0(2) = s0(2) + gk
               s1(1) = s1(1) + ga1(i)*wu(kr)
               s1(2) = s1(2) + gk*wu(kr)
            END DO
            k = MOD(kc - 1 - kv0 + nv, nv)
            CALL kink_pair(2, MIN(k, nv - k), tanv(kc + ivoff),
     &                     gvv_b(ip), guv_b(ip), guu_b(ip), avv(ip),
     &                     auv(ip), auu(ip))
         END IF
         CALL kink_line(tanv(kc + ivoff), gvv_b(ip), guv_b(ip),
     &                  guu_b(ip), avv(ip), auv(ip), auu(ip), alu,
     &                  .FALSE., stu, w0, w1, s0, s1, m0, m1, vg, vk)
         DO k = 1, 3
            i = MOD(kfu + k - 1, nu)*nv + kc
            delgr(i) = delgr(i) - vg(k)
            delgrp(i) = delgrp(i) - vk(k)
         END DO
      END DO

!     CROSSING OF THE KINKS
      CALL kink_corner(guu_b(ip), guv_b(ip), gvv_b(ip), auu(ip),
     &                 auv(ip), avv(ip), alu, hv, lmid, cg, ck)
      IF (lmid) THEN
         DO k = 1, 2
            i = kua + 1 + MOD(kfv + k, nv)
            delgr(i) = delgr(i) - cg/2
            delgrp(i) = delgrp(i) - ck/2
         END DO
      ELSE
         i = kua + 1 + MOD(kfv + 1, nv)
         delgr(i) = delgr(i) - cg
         delgrp(i) = delgrp(i) - ck
      END IF
      END SUBROUTINE kink_grid

!-----------------------------------------------------------------------
!     VALUES ABOUT THE KINKS FOR NV = 1, IMAGE BY IMAGE
!-----------------------------------------------------------------------
      SUBROUTINE kink_images
      kfu = iskip + nu/2 - 1
      IF (kp .EQ. 1) THEN
         rows = 0
         wstep = CMPLX(COS(alu), -SIN(alu), dp)
         wu(iskip) = 1
         DO k = 1, nu - 1
            wu(MOD(iskip + k, nu)) = wu(MOD(iskip + k - 1, nu))*wstep
         END DO
         CALL kink_steps(hv, nvper, .FALSE., gvv_b(ip), stv)
         CALL kink_steps(alu, nu, .FALSE., guu_b(ip), stu)
         cs = 0
!        WEIGHTS OF THE SUMS OVER THE IMAGES, KEPT BESIDE THE SUMS
         DO kr = 0, nu - 1
            CALL kink_weights(tanu(kr + 1 + iuoff), guu_b(ip),
     &                        guv_b(ip), gvv_b(ip), nvper, hv,
     &                        rows(kr + 1,3), rows(kr + 1,4))
         END DO
      END IF

!     LINE OF U OF THIS IMAGE, BUT THAT ON THE KINK OF V
      IF (2*(kp - 1) .NE. nvper) THEN
         CALL kink_weights(tanv(1 + ivoff), gvv_b(ip), guv_b(ip),
     &                     guu_b(ip), nu, alu, w0, w1)
         s0 = 0
         s1 = 0
         DO i = 1, nuv
            IF (rows(i,3) .GT. 0) THEN
               rows(i,1) = rows(i,1) + ga1(i)
               rows(i,2) = rows(i,2) + ga1(i)*ga2(i)
            END IF
            IF (w0 .GT. 0) THEN
               gk = ga1(i)*ga2(i)
               s0(1) = s0(1) + ga1(i)
               s0(2) = s0(2) + gk
               s1(1) = s1(1) + ga1(i)*wu(i - 1)
               s1(2) = s1(2) + gk*wu(i - 1)
            END IF
         END DO
         IF (w0 .GT. 0) THEN
            CALL kink_pair(2, MIN(kp - 1, nvper - kp + 1),
     &                     tanv(1 + ivoff), gvv_b(ip), guv_b(ip),
     &                     guu_b(ip), avv(ip), auv(ip), auu(ip))
         END IF
         CALL kink_line(tanv(1 + ivoff), gvv_b(ip), guv_b(ip),
     &                  guu_b(ip), avv(ip), auv(ip), auu(ip), alu,
     &                  .FALSE., stu, w0, w1, s0, s1, m0, m1, vg, vk)
         DO k = 1, 3
            i = 1 + MOD(kfu + k - 1, nu)
            delgr(i) = delgr(i) - vg(k)
            delgrp(i) = delgrp(i) - vk(k)
         END DO
      END IF

!     LINES OF V, WHICH ARE THE SUMS OVER THE IMAGES, AND THE CROSSING
      IF (kp .EQ. nvper) THEN
         s1 = 0
         DO i = 1, nuv
            IF (i .EQ. kua + 1) THEN
               CYCLE
            END IF
            s0 = rows(i,1:2)
            IF (rows(i,3) .GT. 0) THEN
               k = MOD(i - ip + nu, nu)
               CALL kink_pair(1, MIN(k, nu - k), tanu(i + iuoff),
     &                        guu_b(ip), guv_b(ip), gvv_b(ip), auu(ip),
     &                        auv(ip), avv(ip))
            END IF
            CALL kink_line(tanu(i + iuoff), guu_b(ip), guv_b(ip),
     &                     gvv_b(ip), auu(ip), auv(ip), avv(ip), hv,
     &                     .FALSE., stv, rows(i,3), rows(i,4), s0, s1,
     &                     m0, m1, vg, vk)
            delgr(i) = delgr(i) - SUM(vg)
            delgrp(i) = delgrp(i) - SUM(vk)
         END DO
         CALL kink_corner(guu_b(ip), guv_b(ip), gvv_b(ip), auu(ip),
     &                    auv(ip), avv(ip), alu, hv, .FALSE., cg, ck)
         delgr(kua + 1) = delgr(kua + 1) - cg
         delgrp(kua + 1) = delgrp(kua + 1) - ck
      END IF
      END SUBROUTINE kink_images

!-----------------------------------------------------------------------
!     INTEGRALS OF THE APPROXIMATIONS ALONG THE LINE AT THE STEP D FROM
!     THE SOURCE IN THE FAMILY FAM OF LINES. THOSE OF THE LINE AT THE
!     OPPOSITE STEP ARE THEIR CONJUGATES.
!-----------------------------------------------------------------------
      SUBROUTINE kink_pair(fam, d, s, a, b, c, aa, bb, cc)
      INTEGER, INTENT(in) :: fam, d
      REAL(dp), INTENT(in) :: s, a, b, c, aa, bb, cc
      IF (cs(d,fam) .EQ. 0) THEN
         CALL kink_moments(s, a, b, c, aa, bb, cc, c0(1,d,fam),
     &                     c1(1,d,fam), c0(2,d,fam), c1(2,d,fam))
         cs(d,fam) = s
      END IF
      m0 = c0(:,d,fam)
      IF (s*cs(d,fam) .GT. 0) THEN
         m1 = c1(:,d,fam)
      ELSE
         m1 = CONJG(c1(:,d,fam))
      END IF
      END SUBROUTINE kink_pair

      END SUBROUTINE greenf

