!> \file alias.f

      SUBROUTINE alias_par(gcons, ztemp, gcs, gsc, gcc, gss)
      USE vmec_main
      USE realspace, ONLY:ireflect_par, psqrts
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), PARAMETER :: p5 = 0.5_dp
      REAL(dp), DIMENSION(nzeta,ntheta3,ns), INTENT(out) :: gcons
      REAL(dp), DIMENSION(nzeta,ntheta3,ns), INTENT(in)  :: ztemp
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns) :: gcs, gsc, gcc, gss
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, i, ir, jk, jka, n, k, js, l, j
      INTEGER :: nsmin, nsmax
      INTEGER :: jcount, kk
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: work
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: gcona
      REAL(dp) :: talon, taloff
C-----------------------------------------------
      CALL second0(talon)

      nsmin=tlglob; nsmax=t1rglob

      ALLOCATE (work(4,nzeta,ns), gcona(nzeta,ntheta3,ns))

      gcons(:,:,nsmin:nsmax) = 0
      gcona(:,:,nsmin:nsmax) = 0

      gcs(:ntor,:mpol1,nsmin:nsmax) = 0
      gsc(:ntor,:mpol1,nsmin:nsmax) = 0
      gcc(:ntor,:mpol1,nsmin:nsmax) = 0
      gss(:ntor,:mpol1,nsmin:nsmax) = 0

      !BEGIN M-LOOP
      DO js = nsmin, nsmax
        DO m = 1, mpol1 - 1

          work(:,:,js) = 0
          DO i = 1, ntheta2
            DO k = 1, nzeta
              work(1,k,js) = work(1,k,js) + ztemp(k,i,js)*cosmui(i,m)
              work(2,k,js) = work(2,k,js) + ztemp(k,i,js)*sinmui(i,m)
            END DO

            IF (.not.lasym) CYCLE
            ir = ntheta1 + 2 - i
            IF (i .eq. 1) ir = 1
            DO k = 1, nzeta
              kk=ireflect_par(k)
              work(3,k,js) = work(3,k,js) + 
     1                         ztemp(kk,ir,js)*cosmui(i,m)
              work(4,k,js) = work(4,k,js) + 
     1                         ztemp(kk,ir,js)*sinmui(i,m)
            END DO
          END DO

          IF(js.GT.1) THEN
            DO n = 0, ntor
              DO k = 1, nzeta
                IF (.not.lasym) THEN
                  gcs(n,m,js) = gcs(n,m, js) + tcon(js)*work(1,k,js)*
     1               sinnv(k,n)
                  gsc(n,m,js) = gsc(n,m,js) + tcon(js)*work(2,k,js)*
     1               cosnv(k,n)
                ELSE
                  gcs(n,m,js) = gcs(n,m,js) + p5*tcon(js)*sinnv(k,n)*
     1               (work(1,k,js)-work(3,k,js))
                  gsc(n,m,js) = gsc(n,m,js) + p5*tcon(js)*cosnv(k,n)*
     1               (work(2,k,js)-work(4,k,js))
                  gss(n,m,js) = gss(n,m,js) + p5*tcon(js)*sinnv(k,n)*
     1               (work(2,k,js)+work(4,k,js))
                  gcc(n,m,js) = gcc(n,m,js) + p5*tcon(js)*cosnv(k,n)*
     1               (work(1,k,js)+work(3,k,js))
                END IF
              END DO
            END DO
          END IF
!
!        INVERSE FOURIER TRANSFORM DE-ALIASED GCON
!
          work(:,:,js) = 0

          IF(js.GT.1) THEN
            DO n = 0, ntor
              DO k = 1, nzeta
                work(3,k,js) = work(3,k,js) + gcs(n,m,js)*sinnv(k,n)
                work(4,k,js) = work(4,k,js) + gsc(n,m,js)*cosnv(k,n)
                IF (.not.lasym) CYCLE
                work(1,k,js) = work(1,k,js) + gcc(n,m,js)*cosnv(k,n)
                work(2,k,js) = work(2,k,js) + gss(n,m,js)*sinnv(k,n)
              END DO
            END DO
          END IF

          nsmin=tlglob; nsmax=t1rglob
          DO i = 1, ntheta2
            DO k = 1, nzeta
              gcons(k,i,js) = gcons(k,i,js) + (work(3,k,js)*cosmu(i,m)
     1                     + work(4,k,js)*sinmu(i,m))*faccon(m)
            END DO
            IF (.not.lasym) CYCLE
            DO k = 1, nzeta
              gcona(k,i,js) = gcona(k,i,js) + (work(1,k,js)*cosmu(i,m)
     1                     + work(2,k,js)*sinmu(i,m))*faccon(m)
            END DO
          END DO

        END DO
      END DO
      !END M-LOOP

      IF (lasym) THEN

        !EXTEND GCON INTO THETA = PI,2*PI DOMAIN
        DO js = nsmin, nsmax
          DO i = 1 + ntheta2, ntheta1
            ir = ntheta1 + 2 - i
            DO k = 1, nzeta
              kk=ireflect_par(k)
              gcons(k,i,js) = -gcons(kk,ir,js) + gcona(kk,ir,js)
            END DO
          END DO
        END DO

        !ADD SYMMETRIC, ANTI-SYMMETRIC PIECES IN THETA = 0,PI DOMAIN
        gcons(:,:ntheta2,nsmin:nsmax) = gcons(:,:ntheta2,nsmin:nsmax)
     1                                 + gcona(:,:ntheta2,nsmin:nsmax)

      END IF

      DEALLOCATE (work, gcona)

      CALL second0(taloff)
      alias_time = alias_time + (taloff-talon)

      END SUBROUTINE alias_par
