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
!>  Vertices of the filament for each unit of the ratio of the axis radius
!>  to the least distance from the axis to the boundary.
      REAL(dp), PARAMETER :: vertex_density = 16
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!>  Index of a boundary point or of a vertex, and status of the allocation.
      INTEGER                               :: i
!>  Index of a vertex in one field period.
      INTEGER                               :: j
!>  Index of a field period.
      INTEGER                               :: kper
!>  Index of a toroidal plane.
      INTEGER                               :: kv
!>  Number of vertices in one field period.
      INTEGER                               :: nvsub
!>  Current of the filament.
      REAL(dp), DIMENSION(1)                :: current
!>  Time at the start of the routine.
      REAL(dp)                              :: ttolion
!>  Time at the end of the routine.
      REAL(dp)                              :: ttolioff
!>  Least squared distance from the axis to the boundary in a plane.
      REAL(dp)                              :: dmin
!>  Half the angle from a plane to a vertex, in the toroidal angle of one
!>  field period, then the toroidal angle of the vertex.
      REAL(dp)                              :: argv
!>  Weight of a plane in the interpolant of the axis at a vertex.
      REAL(dp)                              :: weight
!>  Cartesian positions of the vertices.
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xpts
!>  R of the axis at the vertices of one field period.
      REAL(dp), DIMENSION(:), ALLOCATABLE   :: rax
!>  Z of the axis at the vertices of one field period.
      REAL(dp), DIMENSION(:), ALLOCATABLE   :: zax
!>  Cosine of the toroidal angle of the vertices of one field period.
      REAL(dp), DIMENSION(:), ALLOCATABLE   :: cosax
!>  Sine of the toroidal angle of the vertices of one field period.
      REAL(dp), DIMENSION(:), ALLOCATABLE   :: sinax
C-----------------------------------------------

!  Compute wire segments (do not close the loop, closure is done in the
!  biotsavart routines).

!  The filament is a polygon through the magnetic axis. The boundary points
!  lie in the toroidal planes of the grid, where a polygon with its vertices
!  in those planes alone has the field of a corner, not of an arc, unless
!  its sides are short against the distance from the axis to the boundary.
!  Each side is divided into nsub_axis, the axis between the planes being
!  its trigonometric interpolant over the field period.

!  nsub_axis gives the polygon vertex_density vertices for each unit of the
!  ratio of the axis radius to that distance. It is set at the first vacuum
!  call of a run (ivac = 0, or nsub_axis unset in a run restarted at
!  ivac > 0) and kept, so that the field of the filament is a continuous
!  function of the plasma shape during the iterations.
      CALL second0(ttolion)

      IF (ivac .LE. 0 .OR. nsub_axis .EQ. 0) THEN
         dmin = HUGE(dmin)
         DO i = 1, nuv
            kv = 1 + MOD(i - 1,nv)
            dmin = MIN(dmin, (r1b(i) - raxis_nestor(kv))**2 +
     &                       (z1b(i) - zaxis_nestor(kv))**2)
         END DO
         nsub_axis = MAX(1,
     &                   CEILING(vertex_density*MAXVAL(raxis_nestor) /
     &                           (SQRT(dmin)*nvp)))
      END IF
      nvsub = nv*nsub_axis

      ALLOCATE (xpts(3,nvsub*nvper), rax(nvsub), zax(nvsub),
     &          cosax(nvsub), sinax(nvsub), stat=i)
      IF (i .ne. 0) STOP ' allocation error in tolicu'

      current = torcur/mu0

!  Axis at the vertices of one field period. weight is the periodic
!  interpolation kernel of the nv planes.
      DO j = 1, nvsub
         IF (MOD(j - 1,nsub_axis) .EQ. 0) THEN
            kv = 1 + (j - 1)/nsub_axis
            rax(j) = raxis_nestor(kv)
            zax(j) = zaxis_nestor(kv)
            cosax(j) = cosuv(kv)
            sinax(j) = sinuv(kv)
         ELSE
            rax(j) = 0
            zax(j) = 0
            IF (MOD(nv,2) .EQ. 0) THEN
               DO kv = 1, nv
                  argv = p5*pi2*((j - 1)/REAL(nvsub,dp) -
     &                           (kv - 1)/REAL(nv,dp))
                  weight = SIN(nv*argv)/(nv*TAN(argv))
                  rax(j) = rax(j) + weight*raxis_nestor(kv)
                  zax(j) = zax(j) + weight*zaxis_nestor(kv)
               END DO
            ELSE
               DO kv = 1, nv
                  argv = p5*pi2*((j - 1)/REAL(nvsub,dp) -
     &                           (kv - 1)/REAL(nv,dp))
                  weight = SIN(nv*argv)/(nv*SIN(argv))
                  rax(j) = rax(j) + weight*raxis_nestor(kv)
                  zax(j) = zax(j) + weight*zaxis_nestor(kv)
               END DO
            END IF
            argv = (pi2/nvper)*(j - 1)/REAL(nvsub,dp)
            cosax(j) = COS(argv)
            sinax(j) = SIN(argv)
         END IF
      END DO

      i = 0
      DO kper = 1, nvper
         DO j = 1, nvsub
            i = i + 1
            xpts(1,i) = rax(j)*(cosper(kper)*cosax(j) -
     &                          sinper(kper)*sinax(j))
            xpts(2,i) = rax(j)*(sinper(kper)*cosax(j) +
     &                          cosper(kper)*sinax(j))
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
