!> \file symrzl.f

      SUBROUTINE symrzl_par(r1s, rus, rvs, z1s, zus, zvs, lus, lvs,
     1    rcons, zcons, r1a, rua, rva, z1a, zua, zva, lua, lva, rcona,
     2    zcona)
      USE vmec_main
      USE realspace, ONLY: ireflect_par
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(nzeta,ntheta3,ns,0:1), INTENT(inout) ::
     1   r1s, rus, rvs, z1s, zus, zvs, lus, lvs, rcons, zcons
      REAL(dp), DIMENSION(nzeta,ntheta3,ns,0:1), INTENT(in) ::
     1   r1a, rua, rva, z1a, zua, zva, lua, lva, rcona, zcona
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER  :: mpar, ir, i, jk, jka, n2
      INTEGER  :: j, k, nsmin, nsmax
      REAL(dp) :: tsymon, tsymoff
C-----------------------------------------------
!
!     FIRST SUM SYMMETRIC, ANTISYMMETRIC PIECES ON EXTENDED INTERVAL, THETA = [PI,2*PI]
!
      CALl second0(tsymon)
      nsmin = t1lglob
      nsmax = t1rglob

      DO k = nsmin, nsmax
         DO mpar = 0, 1
            DO i = 1 + ntheta2, ntheta1
               ir = ntheta1 + 2 - i                 !-theta
               DO j = 1, nzeta
                  jka = ireflect_par(j)                !-zeta
                  r1s(j,i,k,mpar) = r1s(jka,ir,k,mpar)
     &                            - r1a(jka,ir,k,mpar)
                  rus(j,i,k,mpar) = rua(jka,ir,k,mpar)
     &                            - rus(jka,ir,k,mpar)
                  z1s(j,i,k,mpar) = z1a(jka,ir,k,mpar)
     &                            - z1s(jka,ir,k,mpar)
                  zus(j,i,k,mpar) = zus(jka,ir,k,mpar)
     &                            - zua(jka,ir,k,mpar)
                  lus(j,i,k,mpar) = lus(jka,ir,k,mpar)
     &                            - lua(jka,ir,k,mpar)
                  rcons(j,i,k,mpar) = rcons(jka,ir,k,mpar)
     &                              - rcona(jka,ir,k,mpar)
                  zcons(j,i,k,mpar) = zcona(jka,ir,k,mpar)
     &                              - zcons(jka,ir,k,mpar)
               END DO
               IF (lthreed) THEN
                  DO j = 1, nzeta
                     jka = ireflect_par(j)                !-zeta
                     rvs(j,i,k,mpar) = rva(jka,ir,k,mpar)
     &                               - rvs(jka,ir,k,mpar)
                     zvs(j,i,k,mpar) = zvs(jka,ir,k,mpar)
     &                               - zva(jka,ir,k,mpar)
                     lvs(j,i,k,mpar) = lvs(jka,ir,k,mpar)
     &                               - lva(jka,ir,k,mpar)
                  END DO
               END IF
            END DO

!
!        NOW SUM SYMMETRIC, ANTISYMMETRIC PIECES FOR THETA = [0,PI]
!
            n2 = ntheta2
            r1s(:,:n2,k,mpar) = r1s(:,:n2,k,mpar) + r1a(:,:n2,k,mpar)
            rus(:,:n2,k,mpar) = rus(:,:n2,k,mpar) + rua(:,:n2,k,mpar)
            z1s(:,:n2,k,mpar) = z1s(:,:n2,k,mpar) + z1a(:,:n2,k,mpar)
            zus(:,:n2,k,mpar) = zus(:,:n2,k,mpar) + zua(:,:n2,k,mpar)
            lus(:,:n2,k,mpar) = lus(:,:n2,k,mpar) + lua(:,:n2,k,mpar)
            rcons(:,:n2,k,mpar) = rcons(:,:n2,k,mpar)
     &                          + rcona(:,:n2,k,mpar)
            zcons(:,:n2,k,mpar) = zcons(:,:n2,k,mpar)
     &                          + zcona(:,:n2,k,mpar)
            IF (lthreed) THEN
               rvs(:,:n2,k,mpar) = rvs(:,:n2,k,mpar)
     &                           + rva(:,:n2,k,mpar)
               zvs(:,:n2,k,mpar) = zvs(:,:n2,k,mpar)
     &                           + zva(:,:n2,k,mpar)
               lvs(:,:n2,k,mpar) = lvs(:,:n2,k,mpar)
     &                           + lva(:,:n2,k,mpar)
            END IF

         END DO
      END DO

      CALl second0(tsymoff)
      symrzl_time = symrzl_time + (tsymoff - tsymon)

      END SUBROUTINE symrzl_par
