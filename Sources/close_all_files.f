!> \file close_all_files.f

      SUBROUTINE close_all_files
      USE vparams, ONLY nthreed
      IMPLICIT NONE
C-----------------------------------------------

      IF (nthreed .gt. 0) CLOSE (nthreed)

      END SUBROUTINE close_all_files
