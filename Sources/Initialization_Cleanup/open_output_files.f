!> \file open_output_files.f

      SUBROUTINE open_output_files (extension, lscreen,
     1           lfirst, lwrite)
      USE safe_open_mod
      USE vparams, ONLY: nmac, nthreed, nmac0, nthreed0
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      CHARACTER(LEN=*) :: extension
      LOGICAL :: lscreen, lfirst, lwrite
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iread, inthreed=0, imac0=0
      CHARACTER(LEN=120) :: threed1_file
C-----------------------------------------------
!
!     OPEN FILES FOR READING, WRITING
!
      threed1_file = 'threed1.'//extension

!      PRINT *,'lwrite: ', lwrite,' Threed1 file: ', TRIM(threed1_file)

      IF (lwrite) THEN
         INQUIRE(FILE=threed1_file, OPENED=lfirst)
!         PRINT *,' lfirst: ', lfirst
         lfirst = .not.lfirst
         IF (.NOT.lfirst) RETURN

         IF (lscreen) WRITE (*, '(33('' -''))')
         nthreed = nthreed0
         CALL safe_open(nthreed, iread, threed1_file, 'new',
     1                 'formatted')
            IF (iread .ne. 0) THEN
               IF (iseq .eq. 0 .and. lscreen) PRINT *,
     1         ' VMEC OUTPUT FILES ALREADY EXIST: OVERWRITING THEM ...'
               CALL safe_open(nthreed, inthreed, threed1_file,
     1         'replace', 'formatted')
         ENDIF

         IF (inthreed.ne.0 .or. imac0.ne.0) THEN
            PRINT *,' nthreed = ',      nthreed,
     1              ' istat_threed = ', inthreed
            PRINT *, 'Error opening output file in VMEC ',
     1               'open_output_files'
            STOP 10
         ENDIF
      ENDIF

      END SUBROUTINE open_output_files
