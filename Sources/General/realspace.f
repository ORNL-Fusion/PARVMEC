      MODULE realspace
      USE stel_kinds, ONLY: dp
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), DIMENSION(:,:), ALLOCATABLE ::
     1   r1, ru, rv, zu, zv, rcon, zcon
      REAL(dp), DIMENSION(:,:), ALLOCATABLE, TARGET :: z1
#ifdef _ANIMEC
      REAL(dp), DIMENSION(:), ALLOCATABLE :: 
     1   pperp, ppar, onembc, pp1, pp2, pp3
#endif
      REAL(dp), DIMENSION(:), ALLOCATABLE :: guu, guv, gvv, sigma_an,
     1   ru0, zu0, gcon, rcon0, zcon0, phip, chip, shalf, sqrts, wint
      REAL(dp), DIMENSION(:,:), ALLOCATABLE, TARGET ::
     1   extra1, extra2, extra3, extra4
      INTEGER, DIMENSION(:), ALLOCATABLE :: ireflect_par
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: pru, pz1, pzu, pr1
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: prv, pzv
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: prcon, pzcon
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pgcon
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pshalf
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::
     1   pextra1, pextra2, pextra3, pextra4
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pguu, pguv, pgvv
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pwint, pchip, pphip
      REAL(dp), DIMENSION(:), ALLOCATABLE :: pwint_ns
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pru0, pzu0
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: prcon0, pzcon0
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: psqrts

!  Add variables for odd prefactor (oddpf_)
! _h       On the half-grid
! _f       On the full grid
! _ds_h    Derivative wrt s, on the half grid
      REAL(dp), DIMENSION(:), ALLOCATABLE :: oddpf_h
      REAL(dp), DIMENSION(:), ALLOCATABLE :: oddpf_f
      REAL(dp), DIMENSION(:), ALLOCATABLE :: oddpf_ds_h
      
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: poddpf_h
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: poddpf_f
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: poddpf_ds_h

      END MODULE realspace
