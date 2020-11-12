!> \file totzsp_mod.f

      MODULE totzsp_mod
      USE vmec_main
      USE timer_sub
      IMPLICIT NONE

      INTEGER, PARAMETER, PRIVATE :: m0=0, m1=1, n0=0
      REAL(dp), ALLOCATABLE, PRIVATE :: work1(:,:,:), work2(:,:)
      REAL(dp), PRIVATE :: cosmux, sinmux

      CONTAINS

!-------------------------------------------------------------------------------
!>  @brief Convert symmetric quantities from Fourier space to real space.
!>
!>  Forier transforms between Fourier space and real space. Computes quantities
!>  for R, dR/du, dR/dv, Z, dZ/du, dZ/dv, dlambda/du and dlambda/dv. Non
!>  derivative quantities are trans formed via
!>
!>    A_real = A_mnc*cos(mu - nv) + A_mns*sin(mu - nv)                       (1)
!>
!>  Derivatives with respect to u are transformed as
!>
!>    dA_real/du = -m*A_mnc*sin(mu - nv) + m*A_mns*cos(mu - nv)              (2)
!>
!>  Derivatives with respect to v are transformed as
!>
!>    dA_real/dv = n*A_mnc*sin(mu - nv) - m*A_mns*cos(mu - nv)               (3)
!>
!>  @param[inout] rzl_array Fourier amplitudes for Rmnc, Zmns and Lmns for
!>                          lasym false. When lasym is true, this also contains
!>                          Rmns, Zmnc, Lmnc.
!>  @paran[out]   r11       Real space R.
!>  @param[out]   ru1       Real space dR/du.
!>  @param[out]   rv1       Real space dR/dz.
!>  @param[out]   z11       Real space Z.
!>  @param[out]   zu1       Real space dZ/du.
!>  @param[out]   zv1       Real space dZ/dv.
!>  @param[out]   lu1       Real space dlambda/du.
!>  @param[out]   lv1       Real space dlambda/dv.
!>  @param[out]   rcn1      Unknown R quantity.
!>  @param[out]   zcn1      Unknown Z quantity.
!>  @param[out]   ier_flag  Status of the transform. Takes the value of
!>                          @ref r01_bad_value_flag if rmnc(0,1) is zero.
!>
!>  @note FIXME Figure out what rcn1 and zcn1 are.
!-------------------------------------------------------------------------------
      SUBROUTINE totzsps_par(rzl_array, r11, ru1, rv1, z11, zu1, zv1,
     &                       lu1, lv1, rcn1, zcn1, ier_flag)
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcc, rss, zsc, zcs,
     &                       r01_bad_value_flag
      USE precon2d, ONLY: ictrl_prec2d
      USE parallel_include_module
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,3*ntmax),
     &   TARGET, INTENT(INOUT) :: rzl_array
      REAL(dp), DIMENSION(nzeta,ntheta3,ns,0:1),
     &   INTENT(out) :: r11, ru1,
     &   rv1, z11, zu1,  zv1, lu1, lv1, rcn1, zcn1
      INTEGER, INTENT(inout) :: ier_flag
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: n, m, mparity, k, i, l
      INTEGER :: ioff, joff, mj, ni, nsz
      INTEGER :: nsmin, nsmax, js
      REAL(dp), DIMENSION(:,:,:), POINTER ::
     &   rmncc, rmnss, zmncs, zmnsc, lmncs, lmnsc
      REAL(dp) :: tbroadon, tbroadoff
!-----------------------------------------------
      CALL second0(tffton)

      nsmin = t1lglob
      nsmax = t1rglob

      rmncc=>rzl_array(:,:,:,rcc)                  !COS(mu) COS(nv)
      zmnsc=>rzl_array(:,:,:,zsc+ntmax)            !SIN(mu) COS(nv)
      lmnsc=>rzl_array(:,:,:,zsc+2*ntmax)          !SIN(mu) COS(nv)
      IF (lthreed) THEN
         rmnss=>rzl_array(:,:,:,rss)               !SIN(mu) SIN(nv)
         zmncs=>rzl_array(:,:,:,zcs+ntmax)         !COS(mu) SIN(nv)
         lmncs=>rzl_array(:,:,:,zcs+2*ntmax)       !COS(mu) SIN(nv)
      END IF
      rzl_array(:,m1,1,:) = rzl_array(:,m1,2,:)

      ioff = LBOUND(rmncc,1)
      joff = LBOUND(rmncc,2)
      IF (lthreed) THEN
         CALL convert_sym_par(rmnss(:,m1+joff,:), zmncs(:,m1+joff,:),
     &                        nsmin, nsmax)
      END IF

!
!     ORIGIN EXTRAPOLATION OF M=0 MODES FOR LAMBDA 
!
      IF (lthreed .AND. jlam(m0) .GT. 1) THEN
         lmncs(:,m0+joff,1) = lmncs(:,m0+joff,2)
      END IF

      ALLOCATE (work1(nzeta,12,nsmin:nsmax), stat=i)
      IF (i .ne. 0) THEN
         STOP 'Allocation error in VMEC2000 totzsps'
      END IF

      DO js = nsmin, nsmax
         r11(:,:,js,:) = 0
         ru1(:,:,js,:) = 0
         rv1(:,:,js,:) = 0
         rcn1(:,:,js,:) = 0
         zcn1(:,:,js,:) = 0
         z11(:,:,js,:) = 0
         zu1(:,:,js,:) = 0
         zv1(:,:,js,:) = 0
         lu1(:,:,js,:) = 0
         lv1(:,:,js,:) = 0
         DO m = 0, mpol1
            mparity = MOD(m,2)
            mj = m + joff
            work1(:,:,js) = 0
!
!        INVERSE TRANSFORM IN N-ZETA, FOR FIXED M
!
            DO n = 0, ntor
               ni = n + ioff
               DO k = 1, nzeta
                  work1(k,1,js) = work1(k,1,js)
     &                          + rmncc(ni,mj,js)*cosnv(k,n)
                  work1(k,6,js) = work1(k,6,js)
     &                          + zmnsc(ni,mj,js)*cosnv(k,n)
                  work1(k,10,js) = work1(k,10,js)
     &                           + lmnsc(ni,mj,js)*cosnv(k,n)

                  IF (.NOT.lthreed) CYCLE

                  work1(k,4,js) = work1(k,4,js)
     &                          + rmnss(ni,mj,js)*cosnvn(k,n)
                  work1(k,7,js) = work1(k,7,js)
     &                          + zmncs(ni,mj,js)*cosnvn(k,n)
                  work1(k,11,js) = work1(k,11,js)
     &                           + lmncs(ni,mj,js)*cosnvn(k,n)

                  work1(k,2,js) = work1(k,2,js)
     &                          + rmnss(ni,mj,js)*sinnv(k,n)
                  work1(k,5,js) = work1(k,5,js)
     &                          + zmncs(ni,mj,js)*sinnv(k,n)
                  work1(k,9,js) = work1(k,9,js)
     &                          + lmncs(ni,mj,js)*sinnv(k,n)

                  work1(k,3,js) = work1(k,3,js)
     &                          + rmncc(ni,mj,js)*sinnvn(k,n)
                  work1(k,8,js) = work1(k,8,js)
     &                          + zmnsc(ni,mj,js)*sinnvn(k,n)
                  work1(k,12,js) = work1(k,12,js)
     &                           + lmnsc(ni,mj,js)*sinnvn(k,n)
               END DO
            END DO

!
!        INVERSE TRANSFORM IN M-THETA, FOR ALL RADIAL, ZETA VALUES
!
            l = 0
            DO i = 1, ntheta2
               cosmux = xmpq(m,1)*cosmu(i,m)
               sinmux = xmpq(m,1)*sinmu(i,m)

               r11(:,i,js,mparity) = r11(:,i,js,mparity)
     &                             + work1(:,1,js)*cosmu(i,m)
               ru1(:,i,js,mparity) = ru1(:,i,js,mparity)
     &                             + work1(:,1,js)*sinmum(i,m)
               rcn1(:,i,js,mparity) = rcn1(:,i,js,mparity)
     &                              + work1(:,1,js)*cosmux

               z11(:,i,js,mparity) = z11(:,i,js,mparity)
     &                             + work1(:,6,js)*sinmu(i,m)
               zu1(:,i,js,mparity) = zu1(:,i,js,mparity)
     &                             + work1(:,6,js)*cosmum(i,m)
               zcn1(:,i,js,mparity) = zcn1(:,i,js,mparity)
     &                              + work1(:,6,js)*sinmux

               lu1(:,i,js,mparity) = lu1(:,i,js,mparity)
     &                             + work1(:,10,js)*cosmum(i,m)

               IF (.not.lthreed) CYCLE

               r11(:,i,js,mparity) = r11(:,i,js,mparity)
     &                             + work1(:,2,js)*sinmu(i,m)
               ru1(:,i,js,mparity) = ru1(:,i,js,mparity)
     &                             + work1(:,2,js)*cosmum(i,m)
               rcn1(:,i,js,mparity) = rcn1(:,i,js,mparity)
     &                              + work1(:,2,js)*sinmux

               rv1(:,i,js,mparity) = rv1(:,i,js,mparity)
     &                             + work1(:,3,js)*cosmu(i,m)
     &                             + work1(:,4,js)*sinmu(i,m)
               z11(:,i,js,mparity) = z11(:,i,js,mparity)
     &                             + work1(:,5,js)*cosmu(i,m)

               zu1(:,i,js,mparity) = zu1(:,i,js,mparity)
     &                             + work1(:,5,js)*sinmum(i,m)
               zcn1(:,i,js,mparity) = zcn1(:,i,js,mparity)
     &                              + work1(:,5,js)*cosmux
               zv1(:,i,js,mparity) = zv1(:,i,js,mparity)
     &                             + work1(:,7,js)*cosmu(i,m)
     &                             + work1(:,8,js)*sinmu(i,m)

               lu1(:,i,js,mparity) = lu1(:,i,js,mparity)
     &                             + work1(:,9,js)*sinmum(i,m)
               lv1(:,i,js,mparity) = lv1(:,i,js,mparity)
     &                             - (work1(:,11,js)*cosmu(i,m)
     &                             + work1(:,12,js)*sinmu(i,m))
            END DO
         END DO
      END DO

      DEALLOCATE (work1)

      z01(nsmin:nsmax) = zmnsc(n0+ioff,m1+joff,nsmin:nsmax)
      r01(nsmin:nsmax) = rmncc(n0+ioff,m1+joff,nsmin:nsmax)
      IF (lactive) THEN
         IF (rank.EQ.0 .AND. r01(1).EQ.zero) THEN
            ier_flag = r01_bad_value_flag
         ELSE IF (rank.EQ.0 .AND. r01(1).NE.zero) THEN
            dkappa = z01(1)/r01(1)
         END IF
         CALL second0(tbroadon)
         CALL MPI_Bcast(dkappa,1, MPI_REAL8,0,NS_COMM,MPI_ERR)
         CALL second0(tbroadoff)
         broadcast_time = broadcast_time + (tbroadoff - tbroadon)
      END IF

      CALL second0(tfftoff)
      totzsps_time = totzsps_time + (tfftoff - tffton) 
      timer(tfft) = timer(tfft) + (tfftoff - tffton)

      END SUBROUTINE totzsps_par

!-------------------------------------------------------------------------------
!>  @brief Convert asymmetric quantities from Fourier space to real space.
!>
!>  Forier transforms between Fourier space and real space. Computes quantities
!>  for R, dR/du, dR/dv, Z, dZ/du, dZ/dv, dlambda/du and dlambda/dv. Non
!>  derivative quantities are trans formed via
!>
!>    A_real = A_mnc*cos(mu - nv) + A_mns*sin(mu - nv)                       (1)
!>
!>  Derivatives with respect to u are transformed as
!>
!>    dA_real/du = -m*A_mnc*sin(mu - nv) + m*A_mns*cos(mu - nv)              (2)
!>
!>  Derivatives with respect to v are transformed as
!>
!>    dA_real/dv = n*A_mnc*sin(mu - nv) - m*A_mns*cos(mu - nv)               (3)
!>
!>  @param[inout] rzl_array Fourier amplitudes for Rmnc, Zmns and Lmns for
!>                          lasym false. When lasym is true, this also contains
!>                          Rmns, Zmnc, Lmnc.
!>  @paran[out]   r11       Real space R.
!>  @param[out]   ru1       Real space dR/du.
!>  @param[out]   rv1       Real space dR/dz.
!>  @param[out]   z11       Real space Z.
!>  @param[out]   zu1       Real space dZ/du.
!>  @param[out]   zv1       Real space dZ/dv.
!>  @param[out]   lu1       Real space dlambda/du.
!>  @param[out]   lv1       Real space dlambda/dv.
!>  @param[out]   rcn1      Unknown R quantity.
!>  @param[out]   zcn1      Unknown Z quantity.
!>  @param[out]   ier_flag  Status of the transform. Takes the value of
!>                          @ref r01_bad_value_flag if rmnc(0,1) is zero.
!>
!>  @note FIXME Figure out what rcn1 and zcn1 are.
!-------------------------------------------------------------------------------
      SUBROUTINE totzspa_par(rzl_array, r11, ru1, rv1, z11, zu1, zv1,
     1                       lu1, lv1, rcn1, zcn1)
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcs, rsc, zcc, zss
      USE precon2d, ONLY: ictrl_prec2d
      USE parallel_include_module
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(dp), DIMENSION(nzeta,ntheta3,ns,0:1), INTENT(out) ::
     1   r11, ru1, rv1, z11, zu1, zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, n, mparity, k, i, l, j1
      INTEGER :: ioff, joff, mj, ni
      INTEGER :: nsmin, nsmax, js
      REAL(dp), DIMENSION(:,:,:), POINTER ::
     1           rmncs, rmnsc, zmncc, zmnss, lmncc, lmnss
C-----------------------------------------------
      CALL second0(tffton)
      nsmin = t1lglob
      nsmax = t1rglob

      rmnsc => rzl_array(:,:,:,rsc)               !!SIN(mu) COS(nv)
      zmncc => rzl_array(:,:,:,zcc+ntmax)         !!COS(mu) COS(nv)
      lmncc => rzl_array(:,:,:,zcc+2*ntmax)       !!COS(mu) COS(nv)
      IF (lthreed) THEN
         rmncs => rzl_array(:,:,:,rcs)               !!COS(mu) SIN(nv)
         zmnss => rzl_array(:,:,:,zss+ntmax)         !!SIN(mu) SIN(nv)
         lmnss => rzl_array(:,:,:,zss+2*ntmax)       !!SIN(mu) SIN(nv)
      END IF

!
!     CONVERT FROM INTERNAL XC REPRESENTATION FOR m=1 MODES, R+(at rsc) = .5(rsc + zcc),
!     R-(at zcc) = .5(rsc - zcc), TO REQUIRED rsc, zcc FORMS
!
      ioff = LBOUND(rmnsc,1)
      joff = LBOUND(rmnsc,2)
      CALL convert_asym_par(rmnsc(:,m1+joff,:), zmncc(:,m1+joff,:),
     &                      nsmin, nsmax)

      z00b = zmncc(ioff,joff,ns)

      ALLOCATE (work1(nzeta,12,nsmin:nsmax), stat=i)
      IF (i .NE. 0) THEN
         STOP 'Allocation error in VMEC totzspa'
      END IF

!
!     INITIALIZATION BLOCK
!
      IF (jlam(m0) .gt. 1) THEN
         lmncc(:,m0+joff,1) = lmncc(:,m0+joff,2)
      END IF

      DO js = nsmin, nsmax
         r11(:,:,js,:) = 0
         ru1(:,:,js,:) = 0
         rv1(:,:,js,:) = 0
         rcn1(:,:,js,:) = 0
         zcn1(:,:,js,:) = 0
         z11(:,:,js,:) = 0
         zu1(:,:,js,:) = 0
         zv1(:,:,js,:) = 0
         lu1(:,:,js,:) = 0
         lv1(:,:,js,:) = 0
         DO m = 0, mpol1
            mparity = MOD(m,2)
            mj = m+joff
            work1(:,:,js) = 0
            j1 = jmin1(m)

            DO n = 0, ntor
               ni = n+ioff
               DO k = 1, nzeta
                  work1(k,1,js) = work1(k,1,js)
     &                          + rmnsc(ni,mj,js)*cosnv(k,n)
                  work1(k,6,js) = work1(k,6,js)
     &                          + zmncc(ni,mj,js)*cosnv(k,n)
                  work1(k,10,js) = work1(k,10,js)
     &                           + lmncc(ni,mj,js)*cosnv(k,n)

                  IF (.NOT.lthreed) CYCLE

                  work1(k,2,js) = work1(k,2,js)
     &                          + rmncs(ni,mj,js)*sinnv(k,n)
                  work1(k,3,js) = work1(k,3,js)
     &                          + rmnsc(ni,mj,js)*sinnvn(k,n)
                  work1(k,4,js) = work1(k,4,js)
     &                          + rmncs(ni,mj,js)*cosnvn(k,n)
                  work1(k,5,js) = work1(k,5,js)
     &                          + zmnss(ni,mj,js)*sinnv(k,n)
                  work1(k,7,js) = work1(k,7,js)
     &                          + zmnss(ni,mj,js)*cosnvn(k,n)
                  work1(k,8,js) = work1(k,8,js)
     &                          + zmncc(ni,mj,js)*sinnvn(k,n)
                  work1(k,9,js) = work1(k,9,js)
     &                          + lmnss(ni,mj,js)*sinnv(k,n)
                  work1(k,11,js) = work1(k,11,js)
     &                           + lmnss(ni,mj,js)*cosnvn(k,n)
                  work1(k,12,js) = work1(k,12,js)
     &                           + lmncc(ni,mj,js)*sinnvn(k,n)
               END DO
            END DO

!
!        INVERSE TRANSFORM IN M-THETA
!
            DO i = 1, ntheta2
               cosmux = xmpq(m,1)*cosmu(i,m)
               sinmux = xmpq(m,1)*sinmu(i,m)

               r11(:,i,js,mparity) = r11(:,i,js,mparity)
     &                             + work1(:,1,js)*sinmu(i,m)
               ru1(:,i,js,mparity) = ru1(:,i,js,mparity)
     &                             + work1(:,1,js)*cosmum(i,m)
               z11(:,i,js,mparity) = z11(:,i,js,mparity)
     &                             + work1(:,6,js)*cosmu(i,m)
               zu1(:,i,js,mparity) = zu1(:,i,js,mparity)
     &                             + work1(:,6,js)*sinmum(i,m)
               lu1(:,i,js,mparity) = lu1(:,i,js,mparity)
     &                             + work1(:,10,js)*sinmum(i,m)
               rcn1(:,i,js,mparity) = rcn1(:,i,js,mparity)
     &                              + work1(:,1,js)*sinmux
               zcn1(:,i,js,mparity) = zcn1(:,i,js,mparity)
     &                              + work1(:,6,js)*cosmux

               IF (.not.lthreed) CYCLE

               r11(:,i,js,mparity) = r11(:,i,js,mparity)
     &                             + work1(:,2,js)*cosmu(i,m)
               ru1(:,i,js,mparity) = ru1(:,i,js,mparity)
     &                             + work1(:,2,js)*sinmum(i,m)
               z11(:,i,js,mparity) = z11(:,i,js,mparity)
     &                             + work1(:,5,js)*sinmu(i,m)
               zu1(:,i,js,mparity) = zu1(:,i,js,mparity)
     &                             + work1(:,5,js)*cosmum(i,m)
               lu1(:,i,js,mparity) = lu1(:,i,js,mparity)
     &                             + work1(:,9,js)*cosmum(i,m)
               rcn1(:,i,js,mparity) = rcn1(:,i,js,mparity)
     &                              + work1(:,2,js)*cosmux
               zcn1(:,i,js,mparity) = zcn1(:,i,js,mparity)
     &                              + work1(:,5,js)*sinmux
               rv1(:,i,js,mparity) = rv1(:,i,js,mparity)
     &                             + work1(:,3,js)*sinmu(i,m)
     &                             + work1(:,4,js)*cosmu(i,m)
               zv1(:,i,js,mparity) = zv1(:,i,js,mparity)
     &                             + work1(:,7,js)*sinmu(i,m)
     &                             + work1(:,8,js)*cosmu(i,m)
               lv1(:,i,js,mparity) = lv1(:,i,js,mparity)
     &                             - work1(:,11,js)*sinmu(i,m)
     &                             - work1(:,12,js)*cosmu(i,m)
            END DO
         END DO
      END DO

      DEALLOCATE (work1)

      CALL second0(tfftoff)
      totzspa_time = totzspa_time + (tfftoff - tffton)
      timer(tfft) = timer(tfft) + (tfftoff - tffton)

      END SUBROUTINE totzspa_par

      SUBROUTINE convert_sym_par(rmnss, zmncs, nsmin, nsmax)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: nsmin, nsmax
      REAL(dp), DIMENSION(0:ntor,ns), INTENT(INOUT) :: rmnss, zmncs
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,nsmin:nsmax) :: temp
C-----------------------------------------------
!
!     CONVERT FROM INTERNAL REPRESENTATION TO "PHYSICAL" RMNSS, ZMNCS FOURIER FORM
!     (for lconm1, rss = zmncs)
!
      IF (lconm1) THEN
         temp(:,nsmin:nsmax) = rmnss(:,nsmin:nsmax)
         rmnss(:,nsmin:nsmax) = temp(:,nsmin:nsmax)
     &                        + zmncs(:,nsmin:nsmax)
         zmncs(:,nsmin:nsmax) = temp(:,nsmin:nsmax)
     &                        - zmncs(:,nsmin:nsmax)
      END IF

      END SUBROUTINE convert_sym_par

      SUBROUTINE convert_asym_par(rmnsc, zmncc, nsmin, nsmax)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: nsmin, nsmax
      REAL(dp), DIMENSION(0:ntor,ns), INTENT(INOUT) :: rmnsc, zmncc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,nsmin:nsmax) :: temp
C-----------------------------------------------
!
!     CONVERT FROM INTERNAL REPRESENTATION TO RMNSC, ZMNCC FOURIER FORM
!
      IF (lconm1) THEN
         temp(:,nsmin:nsmax) = rmnsc(:,nsmin:nsmax)
         rmnsc(:,nsmin:nsmax) = temp(:,nsmin:nsmax)
     &                        + zmncc(:,nsmin:nsmax)
         zmncc(:,nsmin:nsmax) = temp(:,nsmin:nsmax)
     &                        - zmncc(:,nsmin:nsmax)
      END IF

      END SUBROUTINE convert_asym_par

      END MODULE totzsp_mod
