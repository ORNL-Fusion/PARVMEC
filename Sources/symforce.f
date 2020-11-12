!> \file symforce.f

      SUBROUTINE symforce_par(ars, brs, crs, azs, bzs, czs, bls, cls,
     &                        rcs, zcs, ara, bra, cra, aza, bza, cza,
     &                        bla, cla, rca, zca)
      USE vmec_main, p5 => cp5
      USE realspace, ONLY: ireflect_par
      USE parallel_include_module
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(nzeta,ntheta3,ns,0:1),
     &   INTENT(inout) :: ars, brs, crs, azs, bzs, czs,
     &   bls, cls, rcs, zcs
      REAL(dp), DIMENSION(nzeta,ntheta3,ns,0:1), INTENT(out) ::
     &   ara, bra, cra, aza, bza, cza, bla, cla, rca, zca
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mpar, ir, i, jk, jka
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ars_0, brs_0, azs_0,
     &              bzs_0, bls_0, rcs_0, zcs_0, crs_0, czs_0, cls_0
      INTEGER :: nsmin, nsmax, j, k
C-----------------------------------------------
      CALL second0(tforon)
      nsmin=t1lglob
      nsmax=t1rglob

      ALLOCATE (ars_0(nzeta,ns), brs_0(nzeta,ns), azs_0(nzeta,ns),
     &          bzs_0(nzeta,ns), bls_0(nzeta,ns), rcs_0(nzeta,ns),
     &          zcs_0(nzeta,ns), crs_0(nzeta,ns), czs_0(nzeta,ns),
     &          cls_0(nzeta,ns), stat=ir)

!
!       SYMMETRIZE FORCES ON RESTRICTED THETA INTERVAL (0 <= u <= pi)
!       SO COS,SIN INTEGRALS CAN BE PERFORMED. FOR EXAMPLE,
!
!       ARS(v,u) = .5*( ARS(v,u) + ARS(-v,-u) )     ! * COS(mu - nv)
!       ARA(v,u) = .5*( ARS(v,u) - ARS(-v,-u) )     ! * SIN(mu - nv)
!
      DO k = nsmin, nsmax
         DO mpar = 0, 1
            DO i = 1, ntheta2
               ir = ntheta1 + 2 - i                 !-theta
               IF (i .eq. 1) THEN
                  ir = 1
               END IF
               DO j = 1, nzeta
                  jka = ireflect_par(j)                !-zeta
                  ara(j,i,k,mpar) = p5*(ars(j,i,k,mpar) -
     &                                  ars(jka,ir,k,mpar))
                  ars_0(j,k)      = p5*(ars(j,i,k,mpar) +
     &                                  ars(jka,ir,k,mpar))
                  bra(j,i,k,mpar) = p5*(brs(j,i,k,mpar) +
     &                                  brs(jka,ir,k,mpar))
                  brs_0(j,k)      = p5*(brs(j,i,k,mpar) -
     &                                  brs(jka,ir,k,mpar))
                  aza(j,i,k,mpar) = p5*(azs(j,i,k,mpar) +
     &                                  azs(jka,ir,k,mpar))
                  azs_0(j,k)      = p5*(azs(j,i,k,mpar) -
     &                                  azs(jka,ir,k,mpar))
                  bza(j,i,k,mpar) = p5*(bzs(j,i,k,mpar) -
     &                                  bzs(jka,ir,k,mpar))
                  bzs_0(j,k)      = p5*(bzs(j,i,k,mpar) +
     &                                  bzs(jka,ir,k,mpar))
                  bla(j,i,k,mpar) = p5*(bls(j,i,k,mpar) -
     &                                  bls(jka,ir,k,mpar))
                  bls_0(j,k)      = p5*(bls(j,i,k,mpar) +
     &                                  bls(jka,ir,k,mpar))
                  rca(j,i,k,mpar) = p5*(rcs(j,i,k,mpar) -
     &                                  rcs(jka,ir,k,mpar))
                  rcs_0(j,k)      = p5*(rcs(j,i,k,mpar) +
     &                                  rcs(jka,ir,k,mpar))
                  zca(j,i,k,mpar) = p5*(zcs(j,i,k,mpar) +
     &                                  zcs(jka,ir,k,mpar))
                  zcs_0(j,k)      = p5*(zcs(j,i,k,mpar) -
     &                                  zcs(jka,ir,k,mpar))
               END DO

               ars(:,i,k,mpar) = ars_0(:,k)
               brs(:,i,k,mpar) = brs_0(:,k)
               azs(:,i,k,mpar) = azs_0(:,k)
               bzs(:,i,k,mpar) = bzs_0(:,k)
               bls(:,i,k,mpar) = bls_0(:,k)
               rcs(:,i,k,mpar) = rcs_0(:,k)
               zcs(:,i,k,mpar) = zcs_0(:,k)

               IF (lthreed) THEN
                  DO j = 1, nzeta
                     jka = ireflect_par(j)
                     cra(j,i,k,mpar)= p5*(crs(j,i,k,mpar) +
     &                                    crs(jka,ir,k,mpar))
                     crs_0(j,k)     = p5*(crs(j,i,k,mpar) -
     &                                    crs(jka,ir,k,mpar))
                     cza(j,i,k,mpar)= p5*(czs(j,i,k,mpar) -
     &                                    czs(jka,ir,k,mpar))
                     czs_0(j,k)     = p5*(czs(j,i,k,mpar) +
     &                                    czs(jka,ir,k,mpar))
                     cla(j,i,k,mpar)= p5*(cls(j,i,k,mpar) -
     &                                    cls(jka,ir,k,mpar))
                     cls_0(j,k)     = p5*(cls(j,i,k,mpar) +
     &                                    cls(jka,ir,k,mpar))
                  END DO

                  crs(:,i,k,mpar) = crs_0(:,k)
                  czs(:,i,k,mpar) = czs_0(:,k)
                  cls(:,i,k,mpar) = cls_0(:,k)
               END IF

            END DO
         END DO
      END DO

      DEALLOCATE (ars_0, brs_0, azs_0, bzs_0, bls_0,
     1            rcs_0, zcs_0, crs_0, czs_0, cls_0, stat=ir)

      CALL second0(tforoff)
      symforces_time = symforces_time + (tforoff - tforon)
      timer(tfor) = timer(tfor) + (tforoff - tforon)

      END SUBROUTINE symforce_par

      SUBROUTINE symoutput(bsq, gsqrt , bsubu , bsubv ,bsupu,
     &                     bsupv, bsubs,
     &                     bsqa, gsqrta, bsubua, bsubva, bsupua,
     &                     bsupva, bsubsa)

      USE vmec_main, p5 => cp5
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(ns*nzeta,ntheta3), INTENT(inout) ::
     &   bsq, gsqrt, bsubu, bsubv, bsupu, bsupv, bsubs
      REAL(dp), DIMENSION(ns*nzeta,ntheta3), INTENT(out)   ::
     &   bsqa,gsqrta,bsubua,bsubva,bsupua,bsupva,bsubsa
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ir, i, jk, jka
      REAL(dp), DIMENSION(ns*nzeta) :: bsq_0, gsqrt_0, bsubu_0,
     &    bsubv_0, bsupu_0, bsupv_0, bsubs_0
C-----------------------------------------------

!
!       SYMMETRIZE FORCES ON RESTRICTED THETA INTERVAL (0 <= u <= pi)
!       SO COS,SIN INTEGRALS CAN BE PERFORMED. FOR EXAMPLE,
!
!       BSQ-S(v,u) = .5*( BSQ(v,u) + BSQ(-v,-u) )     ! * COS(mu - nv)
!       BSQ-A(v,u) = .5*( BSQ(v,u) - BSQ(-v,-u) )     ! * SIN(mu - nv)
!
!       FOR BSUBS, THIS IS REVERSED, S-PIECE ~ SIN, A-PIECE ~ COS
!
      DO i = 1, ntheta2
         ir = ntheta1 + 2 - i                 !-theta
         IF (i == 1) THEN
            ir = 1
         END IF
         DO jk = 1, ns*nzeta
            jka = ireflect(jk)                !-zeta
            bsqa(jk,i)      = p5*(bsq(jk,i)      - bsq(jka,ir))
            bsq_0(jk)       = p5*(bsq(jk,i)      + bsq(jka,ir))
            gsqrta(jk,i)    = p5*(gsqrt(jk,i)    - gsqrt(jka,ir))
            gsqrt_0(jk)     = p5*(gsqrt(jk,i)    + gsqrt(jka,ir))
            bsubua(jk,i)    = p5*(bsubu(jk,i)    - bsubu(jka,ir))
            bsubu_0(jk)     = p5*(bsubu(jk,i)    + bsubu(jka,ir))
            bsubva(jk,i)    = p5*(bsubv(jk,i)    - bsubv(jka,ir))
            bsubv_0(jk)     = p5*(bsubv(jk,i)    + bsubv(jka,ir))
            bsupua(jk,i)    = p5*(bsupu(jk,i)    - bsupu(jka,ir))
            bsupu_0(jk)     = p5*(bsupu(jk,i)    + bsupu(jka,ir))
            bsupva(jk,i)    = p5*(bsupv(jk,i)    - bsupv(jka,ir))
            bsupv_0(jk)     = p5*(bsupv(jk,i)    + bsupv(jka,ir))
! Dominant symmetry reversed
            bsubsa(jk,i)    = p5*(bsubs(jk,i)    + bsubs(jka,ir))
            bsubs_0(jk)     = p5*(bsubs(jk,i)    - bsubs(jka,ir))
         END DO

         bsq(:,i)      = bsq_0(:)
         gsqrt(:,i)    = gsqrt_0(:)
         bsubu(:,i)    = bsubu_0(:)
         bsubv(:,i)    = bsubv_0(:)
         bsupu(:,i)    = bsupu_0(:)
         bsupv(:,i)    = bsupv_0(:)
         bsubs(:,i)    = bsubs_0(:)
      END DO

      END SUBROUTINE symoutput

!  Put the surface routines in a separate subroutine since these the quantites
!  these work on only exist on free boundary runs.
      SUBROUTINE symoutput_sur(bsubu, bsubv, bsupu, bsupv,                     &
     &                         bsubua, bsubva, bsupua, bsupva)
      USE vmec_main, p5 => cp5

      IMPLICIT NONE

C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(nzeta,ntheta3), INTENT(inout) ::
     1   bsubu, bsubv, bsupu, bsupv
      REAL(dp), DIMENSION(nzeta,ntheta2), INTENT(out)   ::
     1   bsubua, bsubva, bsupua, bsupva

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ir, i, jk, jka
      REAL(dp), DIMENSION(nzeta) :: bsubu_0, bsubv_0,
     1                                 bsupu_0, bsupv_0
C-----------------------------------------------

!
!       SYMMETRIZE FORCES ON RESTRICTED THETA INTERVAL (0 <= u <= pi)
!       SO COS,SIN INTEGRALS CAN BE PERFORMED. FOR EXAMPLE,
!
!       BSQ-S(v,u) = .5*( BSQ(v,u) + BSQ(-v,-u) )     ! * COS(mu - nv)
!       BSQ-A(v,u) = .5*( BSQ(v,u) - BSQ(-v,-u) )     ! * SIN(mu - nv)
!
!       FOR BSUBS, THIS IS REVERSED, S-PIECE ~ SIN, A-PIECE ~ COS
!
      ir = 1 !-theta
      DO i = 1, ntheta2
         jka = 1 !-zeta
         DO jk = 1, nzeta
            bsubua(jk,i) = p5*(bsubu(jk,i) - bsubu(jka,ir))
            bsubu_0(jk)  = p5*(bsubu(jk,i) + bsubu(jka,ir))
            bsubva(jk,i) = p5*(bsubv(jk,i) - bsubv(jka,ir))
            bsubv_0(jk)  = p5*(bsubv(jk,i) + bsubv(jka,ir))
            bsupua(jk,i) = p5*(bsupu(jk,i) - bsupu(jka,ir))
            bsupu_0(jk)  = p5*(bsupu(jk,i) + bsupu(jka,ir))
            bsupva(jk,i) = p5*(bsupv(jk,i) - bsupv(jka,ir))
            bsupv_0(jk)  = p5*(bsupv(jk,i) + bsupv(jka,ir))
            jka = nzeta - jk + 1

         END DO
         ir = ntheta3 - i + 1

         bsubu(:,i) = bsubu_0(:)
         bsubv(:,i) = bsubv_0(:)
         bsupu(:,i) = bsupu_0(:)
         bsupv(:,i) = bsupv_0(:)

      END DO

      END SUBROUTINE symoutput_sur
