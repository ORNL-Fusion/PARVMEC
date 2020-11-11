!> \file vmec.f

      PROGRAM vmec

      USE vmec_input
      USE safe_open_mod
      USE vparams, ONLY: nlog, nlog0, nthreed
      USE vmec_params, ONLY: bad_jacobian_flag,
     &    restart_flag, readin_flag, timestep_flag,
     &    output_flag, cleanup_flag,
     &    norm_term_flag, successful_term_flag
      USE parallel_include_module
      USE parallel_vmec_module, ONLY: InitializeParallel,
     &                                FinalizeParallel
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      CHARACTER(LEN=*), PARAMETER ::
     &    increase_niter = "Try increasing NITER",
     &    bad_jacobian = "The jacobian was non-definite!"
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: numargs, ierr_vmec, index_end,
     &   iopen, isnml, iread,
     &   index_dat, iunit, ncount, nsteps, i
      INTEGER :: ictrl(4)
      CHARACTER(LEN=120) :: input_file, reset_file_name, arg
      CHARACTER(LEN=120) :: log_file
      CHARACTER(LEN=120), DIMENSION(10) :: command_arg
      LOGICAL :: lscreen
      REAL(dp) :: ton, toff
      REAL(dp) :: totalton, totaltoff

C-----------------------------------------------
!***
!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a beta version of the PROGRAM VMEC, which is currently
!       under development by S. P. Hirshman at the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a beta version, this program is subject to change
!       and improvement without notice.
!
!       1. CODE SYNOPSIS
!
!       THIS PROGRAM - VMEC (Variational Moments Equilibrium Code)  -
!       SOLVES THREE-DIMENSIONAL MHD EQUILIBRIUM EQUATIONS USING
!       FOURIER SPECTRAL (MOMENTS) METHODS. A CYLINDRICAL COORDINATE
!       REPRESENTATION IS USED (R-Z COORDINATES). THE POLOIDAL
!       ANGLE VARIABLE IS RENORMALIZED THROUGH THE STREAM FUNCTION
!       LAMBDA, WHICH IS SELF-CONSISTENTLY DETERMINED AND DIFFERENCED
!       VARIATIONALLY ON THE HALF-RADIAL MESH. THE POLOIDAL ANGLE IS
!       DETERMINED BY MINIMIZING <M> = m**2 S(m) , WHERE S(m) =
!       Rm**2 + Zm**2 . AN EVEN-ODD DECOMPOSITION IN THE POLOIDAL MODE
!       NO. OF R,Z, AND LAMDA IS USED TO IMPROVE RADIAL RESOLUTION.
!       A FREE-BOUNDARY OPTION IS AVAILABLE (FOR lfreeb=T), WITH A
!       USER-SUPPLIED DATA-FILE "MGRID" NEEDED TO COMPUTE THE PLASMA
!       VACUUM FIELD COMPONENTS BR, BPHI, BZ (see SUBROUTINE BECOIL)
!
!       THE MAGNETIC FIELD IS REPRESENTED INTERNALLY AS FOLLOWS:
!
!       B(s,u,v) = grad(phiT) X ( grad(u) + grad(lambda) ) +
!
!                  iota(s) * grad(v) X grad(phiT)
!
!       WHERE phiT is the toroidal flux (called phi in code) and
!       u,v are the poloidal, toroidal angles, respectively.
!
!       2. ADDITIONAL CODES REQUIRED
!       For the fixed boundary calculation, the user must provide the Fourier
!       coefficients for the plasma boundary (the last surface outside of which
!       the pressure gradient vanishes). For ALL but the simplest geometry, the
!       SCRUNCH code (available from R. Wieland), based on the DESCUR curve-fitting
!       code, can be used to produce the optimized VMEC Fourier representation for
!       an arbritrary closed boundary (it need not be a 'star-like' DOmain, nor
!       need it possess vertical, or 'stellarator', symmetry).
!
!       For the free boundary calculation, the MAKEGRID code (available upon
!       request) is needed to create a binary Green''s FUNCTION table for the
!       vacuum magnetic field(s) and, IF data analysis is to be done, flux and
!       field loops as well. The user provides a SUBROUTINE (BFIELD) which can be
!       called at an arbitrary spatial location and which should RETURN the three
!       cylindrical components of the vacuum field at that point. (Similary,
!       locations of diagnostic flux loops, Rogowski coils, etc. are required IF
!       equilibrium reconstruction is to be done.)
!
!       Plotting is handled by a stand-alone package, PROUT.NCARG (written by
!       R. M. Wieland). It uses NCAR-graphics calls and reads the primary VMEC output
!       file, WOUT.EXT, WHERE 'EXT' is the command-line extension of the INPUT file.
!
!
!       3. UNIX SCRIPT SETUP PARAMETERS
!       The VMEC source code (vmec.lsqh) is actually a UNIX script file which uses
!       the C-precompiler to produce both the machine-specific Fortran source and a
!       make-file specific to ANY one of the following platforms:
!
!       IBM-RISC6000, CRAY, ALPHA (DEC-STATION), HP-UX WORKSTATION,
!       WINDOWS-NT, DEC-VMS
!
!       Additional platforms are easy to add to the existing script as required.
!
!
!       4. FORTRAN PARAMETER STATEMENTS set by user
!       In the Fortran-90 version of VMEC these PARAMETER statements have
!       been replaced by dynamic memory allocation. So the user should set the
!       run-time parameters ns (through ns_array), mpol, ntor in the NAMELIST INDATA.
!
!
!       Added features since last edition (see vmec_params for revision history list)
!       1. Implemented preconditioning algorithm for R,Z
!       2. The physical (unpreconditioned) residuals are used
!          to determine the level of convergence
!       3. The original (MOMCON) scaling of lambda is used, i.e.,
!          Bsupu = phip*(iota - lamda[sub]v)/SQRT(g). This is needed to
!          maintain consistency with the time-stepper for arbitrary PHIP.
!
!       WRITTEN BY S. P. HIRSHMAN (8/28/85 - REVISED 3/1/86) BASED ON
!       1. S. P. Hirshman and J. C. Whitson, Phys. Fluids 26, 3553 (1983).
!       2. S. P. Hirshman and H. K. Meier, Phys. Fluids 28, 1387 (1985).
!       3. S. P. Hirshman and D. K. Lee, Comp. Phys. Comm. 39, 161 (1986).
!

!     Local variables
!
!     ictrl:   array(4) of control variables for running "runvmec" routine
!              see "runvmec" for a description
!

!
!     Read in command-line arguments to get input file or sequence file,
!     screen display information, and restart information
!
      INTERFACE
         SUBROUTINE runvmec(ictrl_array, input_file0,
     &                      lscreen, reset_file_name)
         IMPLICIT NONE
         INTEGER, INTENT(inout), TARGET :: ictrl_array(5)
         LOGICAL, INTENT(in) :: lscreen
         CHARACTER(LEN=*), INTENT(in) :: input_file0
         CHARACTER(LEN=*), OPTIONAL :: reset_file_name
         END SUBROUTINE runvmec
      END INTERFACE

      CALL InitializeParallel

      CALL second0(totalton)
      ton = totalton

!     command line parsing
      CALL getcarg(1, command_arg(1), numargs)
      DO iseq = 2, numargs
         CALL getcarg(iseq, command_arg(iseq), numargs)
      END DO

      CALL second0(toff)
      get_args_time = get_args_time + (toff -ton)

      lscreen = .false.
      IF(grank.EQ.0) lscreen = .true.
      reset_file_name = " "

      IF (numargs .lt. 1) THEN
         STOP 'Invalid command line'
      ELSE IF (command_arg(1).eq.'-h' .or. command_arg(1).eq.'/h') THEN
         PRINT *,
     &   ' ENTER INPUT FILE NAME OR INPUT-FILE SUFFIX ON COMMAND LINE'
         PRINT *
         PRINT *,' For example: '
         PRINT *,'    xvmec input.tftr OR xvmec tftr ',
     &           'OR xvmec ../input.tftr'
         PRINT *
         PRINT *,' Additional (optional) command arguments are',
     &           ' allowed:'
         PRINT *
         PRINT *,'  xvmec <filename> [noscreen] [reset=reset_wout_file]'
         PRINT *
         PRINT *,' noscreen: supresses all output to screen ',
     &           ' (default, or "screen", displays output)'
         PRINT *,' reset: name of reset wout file ',
     &           ' (defaults to none)'
         STOP
      ELSE
         DO iseq = 2, MIN(numargs,10)
            arg = command_arg(iseq)
            IF (TRIM(arg) .eq. 'noscreen' THEN
               lscreen = .false.
            END IF

            index_end = INDEX(arg, "reset=")
            IF (index_end .gt. 0) then
               reset_file_name = arg(index_end+6:)
            end if
         END DO
      END IF

!
!     Determine input file name.
!
      input_file = TRIM(command_arg(1))

!
!     CALL EQUILIBRIUM SOLVER
!
      ictrl = 0
      ictrl(1) = restart_flag + readin_flag  + timestep_flag
     &         + output_flag  + cleanup_flag                ! Sets all flags

      CALL runvmec(ictrl, input_file, lscreen,
     &             reset_file_name)

      ierr_vmec = ictrl(2)
      SELECT CASE (ierr_vmec)
         CASE (bad_jacobian_flag)                               !Bad jacobian even after axis reset and ns->3
            IF (grank .EQ. 0) THEN
               IF (lscreen) WRITE (6, '(/,1x,a)') bad_jacobian
               WRITE (nthreed, '(/,1x,a)') bad_jacobian
               END IF
         CASE DEFAULT
      END SELECT

      CALL second0(totaltoff)
      total_time = total_time + (totaltoff - totalton)

      toff = totaltoff

      CALL FinalizeParallel

      END PROGRAM vmec
