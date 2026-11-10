      SUBROUTINE tolicu (torcur)
      USE vparams, ONLY: mu0
      USE vacmod
      USE vmec_main, ONLY: ivac
      USE biotsavart
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), INTENT(IN) :: torcur
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(dp), PARAMETER :: vertex_density = 16
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, SAVE :: nsub = 0
      INTEGER :: i, j, kper, kv, nvsub
      REAL(dp) :: current(1), ttolion, ttolioff, dmin, argv, weight
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xpts
      REAL(dp), DIMENSION(:), ALLOCATABLE :: rax, zax, cosax, sinax
C-----------------------------------------------
!
!     COMPUTE WIRE SEGMENTS (DO NOT CLOSE LOOP, CLOSURE DONE IN biotsavart ROUTINES)
!
!     THE FILAMENT IS A POLYGON THROUGH THE MAGNETIC AXIS. THE BOUNDARY POINTS
!     LIE IN THE TOROIDAL PLANES OF THE GRID, WHERE A POLYGON WITH ITS VERTICES
!     IN THOSE PLANES ALONE HAS THE FIELD OF A CORNER, NOT OF AN ARC, UNLESS ITS
!     SIDES ARE SHORT AGAINST THE DISTANCE FROM THE AXIS TO THE BOUNDARY. EACH
!     SIDE IS DIVIDED INTO NSUB, THE AXIS BETWEEN THE PLANES BEING ITS
!     TRIGONOMETRIC INTERPOLANT OVER THE FIELD PERIOD.
!
!     NSUB GIVES THE POLYGON vertex_density VERTICES FOR EACH UNIT OF THE RATIO
!     OF THE AXIS RADIUS TO THAT DISTANCE. IT IS SET AT THE FIRST VACUUM CALL
!     OF A RUN (IVAC = 0, OR NSUB UNSET IN A RUN RESTARTED AT IVAC > 0) AND
!     KEPT, SO THAT THE FIELD OF THE FILAMENT IS A CONTINUOUS FUNCTION OF THE
!     PLASMA SHAPE DURING THE ITERATIONS.
!
      CALL second0(ttolion)

      IF (ivac .LE. 0 .OR. nsub .EQ. 0) THEN
         dmin = HUGE(dmin)
         DO i = 1, nuv
            kv = 1 + MOD(i - 1,nv)
            dmin = MIN(dmin, (r1b(i) - raxis_nestor(kv))**2
     &                     + (z1b(i) - zaxis_nestor(kv))**2)
         END DO
         nsub = MAX(1, CEILING(vertex_density*MAXVAL(raxis_nestor)
     &                         /(SQRT(dmin)*nvp)))
      END IF
      nvsub = nv*nsub

      ALLOCATE (xpts(3,nvsub*nvper), rax(nvsub), zax(nvsub),
     &          cosax(nvsub), sinax(nvsub), stat=i)
      IF (i .ne. 0) STOP ' allocation error in tolicu'

      current = torcur/mu0

!
!     AXIS AT THE VERTICES OF ONE FIELD PERIOD. WEIGHT IS THE PERIODIC
!     INTERPOLATION KERNEL OF THE NV PLANES.
!
      DO j = 1, nvsub
         IF (MOD(j - 1,nsub) .EQ. 0) THEN
            kv = 1 + (j - 1)/nsub
            rax(j) = raxis_nestor(kv)
            zax(j) = zaxis_nestor(kv)
            cosax(j) = cosuv(kv)
            sinax(j) = sinuv(kv)
         ELSE
            rax(j) = 0
            zax(j) = 0
            DO kv = 1, nv
               argv = p5*pi2*((j - 1)/REAL(nvsub,dp)
     &                      - (kv - 1)/REAL(nv,dp))
               IF (MOD(nv,2) .EQ. 0) THEN
                  weight = SIN(nv*argv)/(nv*TAN(argv))
               ELSE
                  weight = SIN(nv*argv)/(nv*SIN(argv))
               END IF
               rax(j) = rax(j) + weight*raxis_nestor(kv)
               zax(j) = zax(j) + weight*zaxis_nestor(kv)
            END DO
            argv = (pi2/nvper)*(j - 1)/REAL(nvsub,dp)
            cosax(j) = COS(argv)
            sinax(j) = SIN(argv)
         END IF
      END DO

      i = 0
      DO kper = 1, nvper
         DO j = 1, nvsub
            i = i + 1
            xpts(1,i) = rax(j)*(cosper(kper)*cosax(j)
     &                -         sinper(kper)*sinax(j))
            xpts(2,i) = rax(j)*(sinper(kper)*cosax(j)
     &                +         cosper(kper)*sinax(j))
            xpts(3,i) = zax(j)
         END DO
      END DO
 
!
!     INITIALIZE COIL-RELATED QUANTITIES
!     
      CALL initialize_biotsavart (current, xpt=xpts)

      CALL second0(ttolioff)
      s_tolicu_time = s_tolicu_time + (ttolioff - ttolion)

      END SUBROUTINE tolicu
