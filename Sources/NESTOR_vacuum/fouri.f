!*******************************************************************************
!>  @file extcurz.f
!>  Contains module @ref extcurz.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Fourier transform of the unprimed mesh for the source.
!*******************************************************************************
      MODULE fouri_mod
      USE vacmod
      USE parallel_include_module
      USE timer_sub

      IMPLICIT NONE

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief General Fourier transform of the unprimed mesh for the source.
!>
!>  amatrix(,1) = A(Sin)(Sin')
!>  amatrix(,2) = A(Sin)(Cos')
!>  amatrix(,3) = A(Cos)(Sin')
!>  amatrix(,4) = A(Cos)(Cos')
!>
!>  Argument of the trig function: m'u' - n'v' corresponding to the primed mesh
!>  in PKM eq.(2.14), and mu  - nv (unprimed mesh) in eq.(2.16)
!>
!>  On entry, grpmn(mn,...) has been Fourier-transformed over the unprimed
!>  coordinates (in fourp), so the first index of grpmn corresponds to the first
!>  index of Amatrix. Thus, the Fourier transforms or grpmn here are over the
!>  primed mesh.
!>
!>  In contrast, the integral of the source term over the primed mesh was
!>  already done in scalpot, so here the Fourier transforms are over the
!>  unprimed mesh for the source.
!>
!>  @param[in]    grpmn   Fourier transform of the prime mesh.
!>  @param[in]    gsource Source term.
!>  @param[out]   amatrix A Matrix.
!>  @param[out]   amatrix Square matrix.
!>  @param[inout] bvec    Magnetic field vector.
!>  @param[inout] bvecNS  Magnetic field neumann source.
!>  @param[in]    ndim    Number of dimensions.
!-------------------------------------------------------------------------------
      SUBROUTINE fouri(grpmn, gsource, amatrix, amatsq, bvec, bvecNS,
     &                 ndim)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(mnpd2,nuv3), INTENT(IN)         :: grpmn
      REAL(dp), DIMENSION(nuv), INTENT(IN)                :: gsource
      REAL(dp), DIMENSION(mnpd,mnpd,ndim**2), INTENT(OUT) :: amatrix
      REAL(dp), DIMENSION(mnpd2,mnpd2), INTENT(OUT)       :: amatsq
      REAL(dp), DIMENSION(mnpd,ndim), INTENT(INOUT)       :: bvec
      REAL(dp), DIMENSION(mnpd,ndim), INTENT(INOUT)       :: bvecNS
      INTEGER, INTENT(IN)                                 :: ndim

!  Local Variables
      REAL(dp)                                            :: tfourion
      REAL(dp)                                            :: tfourioff

!  Local Parameters
!  interior  (int_ext=-1), exterior  (int_ext=+1)  neumann problem
      REAL(dp), PARAMETER :: int_ext = 1

!  Start of executable code.
      CALL second0(tfourion)

      IF (.NOT.lasym) THEN
         CALL fouri_sym(grpmn, gsource, amatrix, amatsq, bvecNS, ndim)
      ELSE
         CALL fouri_asym(grpmn, gsource, amatrix, amatsq, bvecNS, ndim)
      END IF

!  Add (still not reduced) analytic and non-singular parts.
      bvec = bvec + bvecNS

      CALL second0(tfourioff)
      timer_vac(tfouri) = timer_vac(tfouri) + (tfourioff-tfourion)
      fouri_time = timer_vac(tfouri)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Fourier transform for symmetric equilibria.
!>
!>  @param[in]    grpmn   Fourier transform of the prime mesh.
!>  @param[in]    gsource Source term.
!>  @param[out]   amatrix A Matrix.
!>  @param[out]   amatrix Square matrix.
!>  @param[inout] bvecNS  Magnetic field neumann source.
!>  @param[in]    ndim    Number of dimensions.
!-------------------------------------------------------------------------------
      SUBROUTINE fouri_sym(grpmn, gsource, amatrix, amatsq, bvecNS,
     &                     ndim)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(mnpd2,nuv3), INTENT(IN)         :: grpmn
      REAL(dp), DIMENSION(nuv), INTENT(IN)                :: gsource
      REAL(dp), DIMENSION(mnpd,mnpd,ndim**2), INTENT(OUT) :: amatrix
      REAL(dp), DIMENSION(mnpd2,mnpd2), INTENT(OUT)       :: amatsq
      REAL(dp), DIMENSION(mnpd,ndim), INTENT(INOUT)       :: bvecNS
      INTEGER, INTENT(IN)                                 :: ndim

!  Local Variables
      INTEGER                                             :: i
      INTEGER                                             :: j
      INTEGER                                             :: k
      INTEGER                                             :: m
      INTEGER                                             :: mn
      INTEGER                                             :: n
      REAL(dp), ALLOCATABLE, DIMENSION(:)                 :: source
      REAL(dp)                                            :: ton
      REAL(dp)                                            :: toff

!  Local Parameters
!  interior  (int_ext=-1), exterior  (int_ext=+1)  neumann problem
      REAL(dp), PARAMETER :: int_ext = 1

!  Start of executable code.
      ALLOCATE (source(nuv3), stat=i)
      IF (i .NE. 0) STOP 'Allocation error in fouri_sym'

!  Stellarator-symmetrize the source terms with respect to u,v, and -u,-v.
!  Index (1) is anti-symmetric and index (2) is symmetric.
!  GSOURCE = - (2pi)**2 B * dS (h - hsing) * WINT
!  WINT: needed to normalize integral over angles to unity
!  Initialize running-sum arrays.
      bvecNS = 0
      amatrix = 0

!  Perform the m,n transforms.
      source(nuv3min:nuv3max) = p5*onp*(gsource(nuv3min:nuv3max)
     &                        - gsource(imirr(nuv3min:nuv3max)))

      j = 0
      mn = nf
      DO i = nuv3min, nuv3max
         mn = 0
         j = j + 1
         DO n = -nf, -1
            mn = mn + 1
            DO m = 1, mf
               mn = mn + 1
               bvecNS(mn,1) = bvecNS(mn,1) + sinmni(j,mn)*source(i)

               amatrix(:,mn,1) = amatrix(:,mn,1)
     &                         + sinmni(j,mn)*grpmn(:mnpd,i)
            END DO
         END DO
         DO n = 0, nf
            mn = mn + 1
            bvecNS(mn,1) = bvecNS(mn,1) + sinmni(j,mn)*source(i)

            amatrix(:,mn,1) = amatrix(:,mn,1)
     &                      + sinmni(j,mn)*grpmn(:mnpd,i)
            DO m = 1, mf
               mn = mn + 1
               bvecNS(mn,1) = bvecNS(mn,1) + sinmni(j,mn)*source(i)

               amatrix(:,mn,1) = amatrix(:,mn,1)
     &                         + sinmni(j,mn)*grpmn(:mnpd,i)
            END DO
         END DO
      END DO

      IF (vlactive) THEN
         CALL second0(ton)
         CALL MPI_Allreduce(MPI_IN_PLACE, amatrix, SIZE(amatrix),
     &                      MPI_REAL8, MPI_SUM, VAC_COMM, MPI_ERR)
         CALL second0(toff)
         allreduce_time = allreduce_time + (toff - ton)
      END IF

      DEALLOCATE (source, stat=i)

!  Sanity checks
!  Index of m,n=(0,0)
!      mn0 = 1 + mf1*nf
!      IF (ANY(bvec(1:mn0-mf1:mf1,:) .NE. 0._dp)) STOP 'BVEC != 0'
!      IF (ANY(amatrix(:,1:mn0-mf1:mf1,:) .ne. 0._dp)) STOP 'AMAT1 != 0'
!      IF (ANY(amatrix(1:mn0-mf1:mf1,:,:) .ne. 0._dp)) STOP 'AMAT2 != 0'

!  Add diagonal terms to the A matrix [The first term in eq(3.2) of PKM]
      DO mn = 1, mnpd
         amatrix(mn,mn,1) = amatrix(mn,mn,1) + pi3*int_ext
      END DO

!  Put elements into a square matrix.
      amatsq(:mnpd,:mnpd) = amatrix(:,:,1) !Sin-Sin'

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Fourier transform for asymmetric equilibria.
!>
!>  @param[in]    grpmn   Fourier transform of the prime mesh.
!>  @param[in]    gsource Source term.
!>  @param[out]   amatrix A Matrix.
!>  @param[out]   amatrix Square matrix.
!>  @param[inout] bvecNS  Magnetic field neumann source.
!>  @param[in]    ndim    Number of dimensions.
!-------------------------------------------------------------------------------
      SUBROUTINE fouri_asym(grpmn, gsource, amatrix, amatsq, bvecNS,
     &                      ndim)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(mnpd2,nuv3), INTENT(IN)         :: grpmn
      REAL(dp), DIMENSION(nuv), INTENT(IN)                :: gsource
      REAL(dp), DIMENSION(mnpd,mnpd,ndim**2), INTENT(OUT) :: amatrix
      REAL(dp), DIMENSION(mnpd2,mnpd2), INTENT(OUT)       :: amatsq
      REAL(dp), DIMENSION(mnpd,ndim), INTENT(INOUT)       :: bvecNS
      INTEGER, INTENT(IN)                                 :: ndim

!  Local Variables
      INTEGER                                             :: i
      INTEGER                                             :: j
      INTEGER                                             :: k
      INTEGER                                             :: m
      INTEGER                                             :: mn
      INTEGER                                             :: mn0
      INTEGER                                             :: n
      REAL(dp), ALLOCATABLE, DIMENSION(:)                 :: source
      REAL(dp)                                            :: ton
      REAL(dp)                                            :: toff

!  Local Parameters
!  interior  (int_ext=-1), exterior  (int_ext=+1)  neumann problem
      REAL(dp), PARAMETER :: int_ext = 1

!  Start of executable code.
      ALLOCATE (source(nuv3), stat=i)
      IF (i .NE. 0) STOP 'Allocation error in fouri_asym'

!  Stellarator-symmetrize the source terms with respect to u,v, and -u,-v.
!  Index (1) is anti-symmetric and index (2) is symmetric.
!  GSOURCE = - (2pi)**2 B * dS (h - hsing) * WINT
!  WINT: needed to normalize integral over angles to unity
!  Initialize running-sum arrays.
      bvecNS = 0
      amatrix = 0

!  Perform the m,n transforms.
      source(nuv3min:nuv3max) = onp*gsource(nuv3min:nuv3max)

      j = 0
      mn = nf
      DO i = nuv3min, nuv3max
         mn = 0
         j = j + 1
         DO n = -nf, -1
            mn = mn + 1
            DO m = 1, mf
               mn = mn + 1
               bvecNS(mn,1) = bvecNS(mn,1) + sinmni(j,mn)*source(i)
               bvecNS(mn,2) = bvecNS(mn,2) + cosmni(j,mn)*source(i)

               amatrix(:,mn,1) = amatrix(:,mn,1)
     &                         + sinmni(j,mn)*grpmn(:mnpd,i)
               amatrix(:,mn,2) = amatrix(:,mn,2)
     &                            + cosmni(j,mn)*grpmn(:mnpd,i)
               amatrix(:,mn,3) = amatrix(:,mn,3)
     &                            + sinmni(j,mn)*grpmn(mnpd+1:,i)
               amatrix(:,mn,4) = amatrix(:,mn,4)
     &                            + cosmni(j,mn)*grpmn(mnpd+1:,i)
            END DO
         END DO
         DO n = 0, nf
            mn = mn + 1
            bvecNS(mn,1) = bvecNS(mn,1) + sinmni(j,mn)*source(i)
            bvecNS(mn,2) = bvecNS(mn,2) + cosmni(j,mn)*source(i)

            amatrix(:,mn,1) = amatrix(:,mn,1)
     &                      + sinmni(j,mn)*grpmn(:mnpd,i)
            amatrix(:,mn,2) = amatrix(:,mn,2)
     &                      + cosmni(j,mn)*grpmn(:mnpd,i)
            amatrix(:,mn,3) = amatrix(:,mn,3)
     &                      + sinmni(j,mn)*grpmn(mnpd+1:,i)
            amatrix(:,mn,4) = amatrix(:,mn,4)
     &                      + cosmni(j,mn)*grpmn(mnpd+1:,i)
            DO m = 1, mf
               mn = mn + 1
               bvecNS(mn,1) = bvecNS(mn,1) + sinmni(j,mn)*source(i)
               bvecNS(mn,2) = bvecNS(mn,2) + cosmni(j,mn)*source(i)

               amatrix(:,mn,1) = amatrix(:,mn,1)
     &                         + sinmni(j,mn)*grpmn(:mnpd,i)
               amatrix(:,mn,2) = amatrix(:,mn,2)
     &                         + cosmni(j,mn)*grpmn(:mnpd,i)
               amatrix(:,mn,3) = amatrix(:,mn,3)
     &                         + sinmni(j,mn)*grpmn(mnpd+1:,i)
               amatrix(:,mn,4) = amatrix(:,mn,4)
     &                         + cosmni(j,mn)*grpmn(mnpd+1:,i)
            END DO
         END DO
      END DO

      IF (vlactive) THEN
         CALL second0(ton)
         CALL MPI_Allreduce(MPI_IN_PLACE, amatrix, SIZE(amatrix),
     &                      MPI_REAL8, MPI_SUM, VAC_COMM, MPI_ERR)
         CALL second0(toff)
         allreduce_time = allreduce_time + (toff - ton)
      END IF

      DEALLOCATE (source, stat=i)

!  Index of m,n=(0,0)
      mn0 = 1 + mf1*nf

!  SANITY CHECKS
!      IF (ANY(bvec(1:mn0-mf1:mf1,:) .NE. 0._dp)) STOP 'BVEC != 0'
!      IF (ANY(amatrix(:,1:mn0-mf1:mf1,:) .ne. 0._dp)) STOP 'AMAT1 != 0'
!      IF (ANY(amatrix(1:mn0-mf1:mf1,:,:) .ne. 0._dp)) STOP 'AMAT2 != 0'

!  Add diagonal terms to the A matrix [The first term in eq(3.2) of PKM]
      DO mn = 1, mnpd
         amatrix(mn,mn,1) = amatrix(mn,mn,1) + pi3*int_ext
         amatrix(mn,mn,4) = amatrix(mn,mn,4) + pi3*int_ext
      END DO
      amatrix(mn0,mn0,4) = amatrix(mn0,mn0,4) + pi3*int_ext ! COS(0-0) mode *2

!  Put elements into a square matrix.
      amatsq(:mnpd,:mnpd) = amatrix(:,:,1)               !Sin-Sin'
      amatsq(:mnpd,1+mnpd:mnpd2) = amatrix(:,:,2)        !Sin-Cos'
      amatsq(1+mnpd:mnpd2,:mnpd) = amatrix(:,:,3)        !Cos-Sin'
      amatsq(1+mnpd:mnpd2,1+mnpd:mnpd2) = amatrix(:,:,4) !Cos-Cos'

      END SUBROUTINE

      END MODULE
