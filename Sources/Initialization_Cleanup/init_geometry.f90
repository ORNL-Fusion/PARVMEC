      MODULE INIT_GEOMETRY

      LOGICAL     :: lflip

      CONTAINS

      SUBROUTINE flip_theta(rmn, zmn, lmn)
      USE vmec_main
      USE vmec_params, ONLY: ntmax, rcc, rss, zsc, zcs,                 &
                                    zcc, zss, rsc, rcs
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(0:ntor,0:mpol1,ntmax),                     &
         INTENT(inout) :: rmn, zmn
      REAL(rprec), DIMENSION(0:ntor,0:mpol1,ntmax),                     &
        INTENT(inout), OPTIONAL :: lmn
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: n, m
      REAL(rprec) :: mul1
      LOGICAL :: l_lmn
!-----------------------------------------------
!
!     FLIP THETA -> PI - THETA (INITIALLY, TO MAKE JACOBIAN < 0)
!
      mul1=-1
      l_lmn = PRESENT(lmn)
      DO m=1,mpol1
         DO n=0,ntor
            rmn(n,m,rcc) = mul1*rmn(n,m,rcc)
            zmn(n,m,zsc) =-mul1*zmn(n,m,zsc)
            IF (l_lmn) lmn(n,m,zsc) =-mul1*lmn(n,m,zsc)
            IF (lthreed) THEN
               rmn(n,m,rss) =-mul1*rmn(n,m,rss)
               zmn(n,m,zcs) = mul1*zmn(n,m,zcs)
               IF (l_lmn) lmn(n,m,zcs) = mul1*lmn(n,m,zcs)
            END IF
            IF (lasym) THEN
               rmn(n,m,rsc) =-mul1*rmn(n,m,rsc)
               zmn(n,m,zcc) = mul1*zmn(n,m,zcc)
               IF (l_lmn) lmn(n,m,zcc) = mul1*lmn(n,m,zcc)
               IF (lthreed) THEN
                  rmn(n,m,rcs) = mul1*rmn(n,m,rcs)
                  zmn(n,m,zss) =-mul1*zmn(n,m,zss)
                  IF (l_lmn) lmn(n,m,zss) =-mul1*lmn(n,m,zss)
               END IF
            END IF
         END DO

         mul1 = -mul1
 
      END DO

      END SUBROUTINE flip_theta

      SUBROUTINE reset_boundary
      USE vmec_main
      USE vmec_params, ONLY: rcc, rss, zsc, zcs, zcc, zss, rsc, rcs,    &
                             signgs
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: m, n, mj, ni, isgn, ioff, joff
      REAL(rprec) :: delta, orient, trc, tzc
      REAL(rprec), DIMENSION(:,:), POINTER ::                           &
         rbcc, rbss, rbcs, rbsc, zbcs, zbsc, zbcc, zbss
      REAL(rprec), ALLOCATABLE :: temp(:)
!-----------------------------------------------
!
!     SETS rmn_bdy, zmn_bdy, lflip AND signgs FROM THE BOUNDARY
!     COEFFICIENTS rbc, zbs, rbs, zbc OF THE INPUT. THE lasym ROTATION
!     CHANGES THOSE FOUR IN PLACE; A SECOND CALL FINDS delta = 0
!
!     CONVERT TO REPRESENTATION WITH RBS(m=1) = ZBC(m=1). THE m=1
!     DETERMINANT, WHICH NO ROTATION OF THETA CHANGES, IS NEGATIVE FOR A
!     BOUNDARY THAT FLIP_THETA REVERSES BELOW; THAT ONE IS ROTATED SO
!     THAT THE REPRESENTATION HOLDS AFTER THE FLIP
!
      IF (lasym) THEN
         orient = SUM(rbc(-ntor:ntor,1))*SUM(zbs(-ntor:ntor,1))         &
                - SUM(rbs(-ntor:ntor,1))*SUM(zbc(-ntor:ntor,1))
         IF (orient .lt. zero) THEN
            delta = -ATAN2(rbs(0,1) + zbc(0,1), zbs(0,1) - rbc(0,1))
         ELSE IF (rbc(0,1) + zbs(0,1) .gt. zero) THEN
            delta = ATAN((rbs(0,1) - zbc(0,1))/                         &
                         (rbc(0,1) + zbs(0,1)))
         ELSE
            delta = ATAN2(rbs(0,1) - zbc(0,1), rbc(0,1) + zbs(0,1))
         END IF
         IF (delta .ne. zero) THEN
            DO m = 0,mpol1
               DO n = -ntor,ntor
                  trc = rbc(n,m)*COS(m*delta) + rbs(n,m)*SIN(m*delta)
                  rbs(n,m) = rbs(n,m)*COS(m*delta)                      &
                           - rbc(n,m)*SIN(m*delta)
                  rbc(n,m) = trc
                  tzc = zbc(n,m)*COS(m*delta) + zbs(n,m)*SIN(m*delta)
                  zbs(n,m) = zbs(n,m)*COS(m*delta)                      &
                           - zbc(n,m)*SIN(m*delta)
                  zbc(n,m) = tzc
               END DO
            END DO
         END IF
      END IF

!
!     CONVERT TO INTERNAL REPRESENTATION OF MODES
!
!     R = RBCC*COS(M*U)*COS(N*V) + RBSS*SIN(M*U)*SIN(N*V)
!       + RBCS*COS(M*U)*SIN(N*V) + RBSC*SIN(M*U)*COS(N*V)
!     Z = ZBCS*COS(M*U)*SIN(N*V) + ZBSC*SIN(M*U)*COS(N*V)
!       + ZBCC*COS(M*U)*COS(N*V) + ZBSS*SIN(M*U)*SIN(N*V)
!
!
!     POINTER ASSIGNMENTS (NOTE: INDICES START AT 1, NOT 0, FOR POINTERS, EVEN THOUGH
!                          THEY START AT ZERO FOR RMN_BDY)
!     ARRAY STACKING ORDER DETERMINED HERE
!
      rbcc => rmn_bdy(:,:,rcc)
      zbsc => zmn_bdy(:,:,zsc)
      IF (lthreed) THEN
         rbss => rmn_bdy(:,:,rss)
         zbcs => zmn_bdy(:,:,zcs)
      END IF

      IF (lasym) THEN
         rbsc => rmn_bdy(:,:,rsc)
         zbcc => zmn_bdy(:,:,zcc)
         IF (lthreed) THEN
            rbcs => rmn_bdy(:,:,rcs)
            zbss => zmn_bdy(:,:,zss)
         END IF
      END IF

      rmn_bdy = 0
      zmn_bdy = 0

      ioff = LBOUND(rbcc,1)
      joff = LBOUND(rbcc,2)

      DO m = 0, mpol1
         mj = m + joff
         IF (lfreeb .and.                                               &
             (mfilter_fbdy.gt.1 .and. m.gt.mfilter_fbdy)) THEN
            CYCLE
         END IF
         DO n = -ntor, ntor
            IF (lfreeb .and.                                            &
                (nfilter_fbdy.gt.0 .and. ABS(n).gt.nfilter_fbdy)) THEN
               CYCLE
            END IF
            ni = ABS(n) + ioff
            IF (n .eq. 0) THEN
               isgn = 0
            ELSE IF (n .gt. 0) THEN
               isgn = 1
            ELSE
               isgn = -1
            END IF
            rbcc(ni,mj) = rbcc(ni,mj) + rbc(n,m)
            IF (m .gt. 0) THEN
               zbsc(ni,mj) = zbsc(ni,mj) + zbs(n,m)
            END IF

            IF (lthreed) THEN
               IF (m .gt. 0) THEN
                  rbss(ni,mj) = rbss(ni,mj) + isgn*rbc(n,m)
               END IF
               zbcs(ni,mj) = zbcs(ni,mj) - isgn*zbs(n,m)
            END IF

            IF (lasym) THEN
               IF (m .gt. 0) THEN
                  rbsc(ni,mj) = rbsc(ni,mj) + rbs(n,m)
               END IF
               zbcc(ni,mj) = zbcc(ni,mj) + zbc(n,m)
               IF (lthreed) THEN
                  rbcs(ni,mj) = rbcs(ni,mj) - isgn*rbs(n,m)
                  IF (m .gt. 0) THEN
                     zbss(ni,mj) = zbss(ni,mj) + isgn*zbc(n,m)
                  END IF
               END IF
            END IF
         END DO
      END DO

!
!     CHECK SIGN OF JACOBIAN (SHOULD BE SAME AS SIGNGS)
!     FOR lasym, ORIENT IS THE m=1 DETERMINANT TAKEN ABOVE
!
      mj = 1 + joff
      IF (.not.lasym) THEN
         orient = SUM(rbcc(1:ntor1,mj))*SUM(zbsc(1:ntor1,mj))
      END IF
      lflip = (orient .lt. zero)
      signgs = -1
      IF (lflip) THEN
         CALL flip_theta(rmn_bdy, zmn_bdy)
      END IF

!
!     CONVERT TO INTERNAL FORM FOR (CONSTRAINED) m=1 MODES
!     INTERNALLY, FOR m=1: XC(rss) = .5(RSS+ZCS), XC(zcs) = .5(RSS-ZCS)
!     WITH XC(zcs) -> 0 FOR POLAR CONSTRAINT
!     (see convert_sym, convert_asym in totzsp_mod file)
!
      IF (lconm1 .and. (lthreed .or. lasym)) THEN
         ALLOCATE (temp(SIZE(rbcc,1)))
         IF (lthreed) THEN
            temp = rbss(:,mj)
            rbss(:,mj) = cp5*(temp(:) + zbcs(:,mj))
            zbcs(:,mj) = cp5*(temp(:) - zbcs(:,mj))
         END IF
         IF (lasym) THEN
            temp = rbsc(:,mj)
            rbsc(:,mj) = cp5*(temp(:) + zbcc(:,mj))
            zbcc(:,mj) = cp5*(temp(:) - zbcc(:,mj))
         END IF
         DEALLOCATE (temp)
      END IF

      END SUBROUTINE reset_boundary

      END MODULE INIT_GEOMETRY
