!*******************************************************************************
!>  @file restart.f
!>  @brief Contains module @ref restart.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module for reading and writting a VMEC restart file.
!*******************************************************************************
      MODULE restart
      USE stel_kinds
      USE xstuff
      USE netcdf_inc
      USE profiler
      USE v3_utilities

      IMPLICIT NONE

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) m grid base class.
!
!*******************************************************************************
      TYPE :: restart_class
!>  File id
         INTEGER   :: file_id
!>  xc dimid.
         INTEGER   :: dimid
!>  xc variable id.
         INTEGER   :: varid
!>  Flag to indicate if VMEC was restarted.
         LOGICAL   :: restarted
      CONTAINS
         PROCEDURE :: read => restart_read
         PROCEDURE :: write => restart_write
         FINAL     :: restart_destruct
      END TYPE

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref restart_class object.
!>
!>  @param[in] file_name Name of the restart file.
!>  @returns A @ref restart_class object.
!-------------------------------------------------------------------------------
      FUNCTION restart_construct_new(file_name)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (restart_class), POINTER :: restart_construct_new
      CHARACTER (len=*), INTENT(in)  :: file_name

!  local variables
      REAL (rprec)                   :: start_time
      INTEGER                        :: status

!  Start of executable code.
      start_time = profiler_get_start_time()

      ALLOCATE(restart_construct_new)

      restart_construct_new%restarted = .false.

      status = nf_create(TRIM(file_name), NF_CLOBBER,                          &
     &                   restart_construct_new%file_id)
      CALL assert_eq(status, 0, nf_strerror(status))

      status = nf_def_dim(restart_construct_new%file_id, 'xc_dim',             &
     &                    SIZE(xc), restart_construct_new%dimid)
      CALL assert_eq(status, 0, nf_strerror(status))

      status = nf_def_var(restart_construct_new%file_id, 'xc',                 &
     &                    NF_DOUBLE, 1,                                        &
     &                    (/ restart_construct_new%dimid /),                   &
     &                    restart_construct_new%varid)
      CALL assert_eq(status, 0, nf_strerror(status))

      status = nf_enddef(restart_construct_new%file_id)
      CALL assert_eq(status, 0, nf_strerror(status))

      CALL profiler_set_stop_time('restart_construct_new', start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct a @ref restart_class object.
!>
!>  @param[in] file_name Name of the restart file.
!>  @returns A @ref restart_class object.
!-------------------------------------------------------------------------------
      FUNCTION restart_construct_open(file_name)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (restart_class), POINTER :: restart_construct_open
      CHARACTER (len=*), INTENT(in)  :: file_name

!  local variables
      REAL (rprec)                   :: start_time
      INTEGER                        :: status
      INTEGER                        :: xc_len

!  Start of executable code.
      start_time = profiler_get_start_time()

      ALLOCATE(restart_construct_open)

      restart_construct_open%restarted = .true.

      status = nf_open(TRIM(file_name), NF_WRITE,                              &
     &                 restart_construct_open%file_id)
      CALL assert_eq(status, 0, nf_strerror(status))

      status = nf_inq_dimid(restart_construct_open%file_id, 'xc_dim',          &
     &                      restart_construct_open%dimid)
      CALL assert_eq(status, 0, nf_strerror(status))

      status = nf_inq_dimlen(restart_construct_open%file_id,                   &
     &                       restart_construct_open%dimid,                     &
     &                       xc_len)
      CALL assert_eq(status, 0, nf_strerror(status))

      CALL assert_eq(xc_len, SIZE(xc),                                         &
     &               'Cannot restart xc array size mismatch.')

      status = nf_inq_varid(restart_construct_open%file_id, 'xc',              &
     &                      restart_construct_open%varid)
      CALL assert_eq(status, 0, nf_strerror(status))

      CALL profiler_set_stop_time('restart_construct_open', start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Destruct a @ref restart_class object.
!>
!>  @param[inout] this A @ref restart_class object instance.
!-------------------------------------------------------------------------------
      SUBROUTINE restart_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (restart_class), INTENT(inout) :: this

!  local variables
      INTEGER                             :: status

!  Start of executable code
      status = nf_close(this%file_id)
      CALL assert_eq(status, 0, nf_strerror(status))

      END SUBROUTINE

!*******************************************************************************
!  UTILITIES SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Write the retart file.
!>
!>  @param[in] this A @ref restart_class object instance.
!-------------------------------------------------------------------------------
      SUBROUTINE restart_write(this)
 
      IMPLICIT NONE

!  Declare Arguments
      CLASS (restart_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                      :: start_time
      INTEGER                           :: status

!  Start of executable code.
      start_time = profiler_get_start_time()

      status = NF_PUT_VAR_DOUBLE(this%file_id, this%varid, xc)
      CALL assert_eq(status, 0, nf_strerror(status))

      status = nf_sync(this%file_id)
      CALL assert_eq(status, 0, nf_strerror(status))

      CALL profiler_set_stop_time('restart_write', start_time)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write the retart file.
!>
!>  @param[in] this A @ref restart_class object instance.
!-------------------------------------------------------------------------------
      SUBROUTINE restart_read(this)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (restart_class), INTENT(in) :: this

!  local variables
      REAL (rprec)                      :: start_time
      INTEGER                           :: status

!  Start of executable code.
      start_time = profiler_get_start_time()

      status = NF_GET_VAR_DOUBLE(this%file_id, this%varid, xc)
      CALL assert_eq(status, 0, nf_strerror(status))

      CALL profiler_set_stop_time('restart_read', start_time)

      END SUBROUTINE

      END MODULE
