!> \file spectrum.f

      SUBROUTINE spectrum_par(rmn, zmn)
      USE parallel_include_module
      USE vmec_main
      USE vmec_params, ONLY: mscale, nscale, ntmax, rss, zcs, rsc, zcc
      USE totzsp_mod, ONLY:  convert_sym_par, convert_asym_par

      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,ntmax), 
     1             INTENT(inout) :: rmn, zmn
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: m1 = 1
      INTEGER :: js, ntype, n, m, nsmin, nsmax
      REAL(dp), DIMENSION(ns) :: t1, dnumer, denom
      REAL(dp) :: scale
C-----------------------------------------------
!
!     CONVERT FROM INTERNAL XC REPRESENTATION FOR m=1 MODES, R+(at rsc) = .5(rsc + zcc),
!     R-(at zcc) = .5(rsc - zcc), TO REQUIRED rsc, zcc FORMS
!
      nsmin=MAX(2,tlglob); nsmax=MIN(t1rglob,ns)

      IF (lthreed) THEN
         CALL convert_sym_par(rmn(:,m1,:,rss), zmn(:,m1,:,zcs),
     &                        nsmin, nsmax)
      END IF
      IF (lasym) THEN
         CALL convert_asym_par(rmn(:,m1,:,rsc), zmn(:,m1,:,zcc),
     &                         nsmin, nsmax)
      END IF

      dnumer(nsmin:nsmax) = zero
      denom(nsmin:nsmax) = zero
      DO ntype = 1,ntmax
         DO n = 0,ntor
            DO m = 1,mpol1
               scale = (mscale(m)*nscale(n))**2
               DO js = nsmin,nsmax
                  t1(js) = (rmn(n,m,js,ntype)**2 +
     &                      zmn(n,m,js,ntype)**2)*scale
               END DO
               dnumer(nsmin:nsmax) = dnumer(nsmin:nsmax)
     &                             + t1(nsmin:nsmax)*xmpq(m,3)
               denom(nsmin:nsmax) = denom (nsmin:nsmax)
     &                            + t1(nsmin:nsmax)*xmpq(m,2)
            END DO
         END DO
      END DO

      specw(nsmin:nsmax) = dnumer(nsmin:nsmax)/denom(nsmin:nsmax)

      END SUBROUTINE spectrum_par
