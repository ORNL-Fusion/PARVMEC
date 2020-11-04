!> \file precon2d.f

      MODULE precon2d
      USE stel_kinds, ONLY: dp
      USE vmec_dim
      USE vmec_params
      USE vparams, ONLY: nthreed, one, zero
      USE vmec_input, ONLY: ntor, nzeta, lfreeb, lasym
      USE timer_sub
      USE safe_open_mod
      USE directaccess
      USE parallel_include_module

      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PRIVATE, PARAMETER :: sp = dp
      INTEGER, PRIVATE :: ntyptot, m_2d, n_2d, ntype_2d
      INTEGER, PRIVATE, ALLOCATABLE :: ipiv_blk(:,:)
      INTEGER, PRIVATE :: mblk_size
      INTEGER, PRIVATE :: mystart(3), myend(3) 
      LOGICAL, PRIVATE :: FIRSTPASS=.TRUE.
      REAL(sp),PRIVATE, ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: 
     1    block_diag, block_plus, block_mins,
     2    block_dsave, block_msave, block_psave
      REAL(sp),PRIVATE, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) ::
     1    block_diag_sw, block_plus_sw, block_mins_sw
   
      REAL(dp), PRIVATE, DIMENSION(:,:,:,:), ALLOCATABLE :: gc_save
      REAL(dp) :: ctor_prec2d
      INTEGER :: ictrl_prec2d
      LOGICAL :: lHess_exact = .TRUE.,  !FALSE -> FASTER, LESS ACCURATE VACUUM CALCULATION OF HESSIAN
     1           l_backslv = .FALSE.,
     2           l_comp_prec2D = .TRUE., 
     3           l_edge  = .FALSE.,                !=T IF EDGE PERTURBED
     4           edge_mesh(3)

      LOGICAL, PARAMETER, PRIVATE :: lscreen = .FALSE.
      INTEGER, PARAMETER, PRIVATE :: jstart(3) = (/1,2,3/)

      PRIVATE :: swap_forces, reswap_forces 

!
!     Direct-Access (swap to disk) stuff
!      CHARACTER(LEN=3)   :: FlashDrive ="F:\"
      CHARACTER(LEN=3)   :: FlashDrive =""
      CHARACTER(LEN=128) :: ScratchFile=""
      INTEGER, PARAMETER :: blmin=1, bldia=2, blpls=3
      REAL(sp), ALLOCATABLE, DIMENSION(:,:,:) :: DataItem
      INTEGER            :: iunit_dacess=10
      LOGICAL :: lswap2disk = .FALSE.                                 !Set internally if blocks do not fit in memory
      INTEGER, PARAMETER :: LOWER=3,DIAG=2,UPPER=1

C-----------------------------------------------
!
!     SP:        forces single precision for blocks (smaller size)
!     ICTRL_PREC2D: controls initialization and application of 2d block preconditioner
!                   = 0, no preconditioner applied
!                   = 1, apply preconditioner
!                   = 2, initial call of funct3d to set up residue vector, store saved vectors,
!                        and (for .not.lasym case), call LAMBLKS routine
!                   = 3, radial jog vector is being computed to calculate hessian elements
!     L_BACKSLV: if true, test that Hessian is inverted correctly by back-solving
!     LHESS_EXACT : if true, edge value of ctor (in bcovar) is computed as a constant to make
!                   Hessian symmetric. Also, sets ivacskip=0 in call to vacuum in computation of
!                   Hessian. However, the ivacskip=0 option is (very) slow and found not to be necessary
!                   in practice. Set this true primarily for debugging purposes (check Ap ~ -p in MatVec 
!                   routine, for example, in GMRes module)
!

      CONTAINS

      SUBROUTINE swap_forces(gc, temp, mblk, nblocks)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: mblk, nblocks
      REAL(dp), DIMENSION(nblocks,mblk), INTENT(in)  :: gc
      REAL(dp), DIMENSION(mblk,nblocks), INTENT(out) :: temp
C-----------------------------------------------
!
!     reorders forces (gc) array prior to applying 
!     block-tridiagonal pre-conditioner. on exit, temp is the reordered array
!     flip sign so eigenvalue is negative (corresponding to damping)
!
      temp = -TRANSPOSE(gc)
    
      END SUBROUTINE swap_forces

      SUBROUTINE reswap_forces(temp, gc, mblk, nblocks)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: mblk, nblocks
      REAL(dp), DIMENSION(nblocks,mblk), INTENT(inout)  :: gc
      REAL(dp), DIMENSION(mblk,nblocks), INTENT(in) :: temp
C-----------------------------------------------
!
!     Following application of block pre-conditioner, restores original
!     order of forces (gc) array previously ordered by call to "swap_forces"
!
      gc = TRANSPOSE(temp)

      END SUBROUTINE reswap_forces

      SUBROUTINE block_precond_par(gc)
      USE blocktridiagonalsolver, ONLY: SetMatrixRHS
      USE blocktridiagonalsolver, ONLY: BackwardSolve 
      USE blocktridiagonalsolver, ONLY: GetSolutionVector
      USE parallel_include_module
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,ntyptot) :: gc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: mblk, istat, globrow
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: tmp
      REAL(dp) :: ton, toff
      REAL(dp), DIMENSION (:,:), ALLOCATABLE :: solvec
C-----------------------------------------------
      IF (.NOT.lactive) THEN
         RETURN
      END IF
!
!     Applies 2D block-preconditioner to forces vector (gc)
!
      IF (ntyptot .LE. 0) THEN
         STOP 'ntyptot must be > 0'
      END IF

      mblk = ntyptot*mnsize

!     Apply preconditioner to temp, using LU factors stored in block_... matrices

      ALLOCATE (tmp(mblk,ns), stat=istat)
      CALL tolastns(gc, tmp)
      tmp(:,tlglob:trglob) = -tmp(:,tlglob:trglob)
      
      DO globrow=tlglob, trglob
         CALL SetMatrixRHS(globrow,tmp(:,globrow))
      END DO
      DEALLOCATE (tmp, stat=istat)

      CALL second0(ton)
      CALL BackwardSolve
      CALL second0(toff)
      bcyclic_backwardsolve_time=bcyclic_backwardsolve_time+(toff-ton)
       
      ALLOCATE (solvec(mblk,ns), stat=istat)
      IF (istat .NE. 0) THEN
         STOP 'Allocation error in block_precond_par before gather'
      END IF

      DO globrow=tlglob, trglob
         CALL GetSolutionVector (globrow,solvec(:,globrow))
      END DO

      CALL tolastntype(solvec, gc)

      CALL Gather4XArray(gc)
      
      DEALLOCATE (solvec)
   
      END SUBROUTINE block_precond_par

      SUBROUTINE compute_blocks_par (xc, xcdot, gc)
      USE blocktridiagonalsolver, ONLY: ForwardSolve
      USE parallel_include_module
      USE vmec_main, ONLY: iter2
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp),DIMENSION(0:ntor,0:mpol1,ns,3*ntmax) :: xc, gc, xcdot
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(dp), PARAMETER :: p5 = 0.5_dp
      INTEGER :: m, n, i, ntype, istat, mblk, ibsize, iunit
      REAL(dp) :: time_on, time_off, bsize, tprec2don, tprec2doff
      REAL(dp) :: ton, toff
      CHARACTER(LEN=100):: label
      LOGICAL, PARAMETER :: lscreen = .false.
      INTEGER :: j, k, l
C-----------------------------------------------
!
!     COMPUTES THE JACOBIAN BLOCKS block_mins, block_diag, block_plus
!     USING EITHER A SLOW - BUT RELIABLE - "JOG" TECHNIQUE, OR
!     USING PARTIAL ANALYTIC FORMULAE.
!
!     THE SUBROUTINE lam_blks IS CALLED FROM BCOVAR TO COMPUTE 
!     THE ANALYTIC BLOCK ELEMENTS
!
      CALL second0(tprec2don)
      
      IF (l_backslv .and. sp.ne.dp) THEN
         STOP 'Should set sp = dp!'
      END IF

      ntyptot = SIZE(gc,4)
      IF (ntyptot .NE. 3*ntmax) THEN
         STOP ' NTYPTOT != 3*ntmax'
      END IF
      mblk = ntyptot*mnsize

      bsize = REAL(mblk*mblk, dp)*3*KIND(block_diag)
      IF (bsize .gt. HUGE(mblk)) THEN
         WRITE (6, *) ' bsize: ', bsize, ' exceeds HUGE(int): ',
     &                HUGE(mblk)
!        WRITE (6, *) ' Blocks will be written to disk.'
!        lswap2disk = .TRUE.
      ELSE
         lswap2disk = .FALSE.
      END IF

      bsize = bsize*ns
      IF (bsize .lt. 1.E6_dp) THEN
         ibsize = bsize/1.E1_dp
         label = " Kb"
      ELSE IF (bsize .lt. 1.E9_dp) THEN
         ibsize = bsize/1.E4_dp
         label = " Mb"
      ELSE
         ibsize = bsize/1.E7_dp
         label = " Gb"
      END IF

      DO i = 1,2
         IF (i .eq. 1) THEN
            iunit = 6
         END IF
         IF (i .eq. 2) THEN
            iunit = nthreed
         END IF
         IF (grank.EQ.0) THEN
            WRITE (iunit, '(/,2x,a,i5,a,/,2x,a,i5,a)')
     &         'Initializing 2D block preconditioner at ', iter2,
     &         ' iterations',
     &         'Estimated time to compute Hessian = ',
     &         3*ntyptot*mnsize,' VMEC time steps'
            WRITE (iunit, '(2x,a,i4,a,f12.2,a)') 'Block dim: ', mblk,
     &         '^2  Preconditioner size: ', REAL(ibsize)/100,
     &         TRIM(label)
         END IF
      END DO
!
!     COMPUTE AND STORE BLOCKS (MN X MN) FOR PRECONDITIONER
!
      CALL second0(time_on)

      ALLOCATE (gc_save(0:ntor,0:mpol1,ns,ntyptot), stat=istat)
      IF (istat .NE. 0) THEN
         STOP 'Allocation error: gc_save in compute_blocks'
      END IF

      IF (ALLOCATED(block_diag)) THEN
         DEALLOCATE (block_diag, block_plus, block_mins, stat=istat)
         IF (istat .ne. 0) THEN
            STOP 'Deallocation error in compute blocks'
         END IF
      END IF

!
!     GENERAL (SLOWER BY 2/3 THAN SYMMETRIC VERSION) METHOD: ASSUMES NO SYMMETRIES OF R, Z COEFFICIENTS
!
      CALL sweep3_blocks_par (xc, xcdot, gc)
      IF (lactive) THEN
         CALL compute_col_scaling_par
      END IF

      ictrl_prec2d = 1                 !Signals funct3d (residue) to use block preconditioner

      CALL second0(time_off)
      IF (grank .EQ. 0) THEN
         WRITE (6,1000) time_off - time_on
         WRITE (nthreed,1000) time_off - time_on
      END IF

!
!     FACTORIZE HESSIAN 
!
      CALL second0(time_on)
      IF (ALLOCATED(ipiv_blk)) THEN
         DEALLOCATE(ipiv_blk, stat=ntype)
      END IF
      ALLOCATE (ipiv_blk(mblk,ns), stat=ntype)     
      IF (ntype .ne. 0) THEN
         STOP 'Allocation error2 in block_precond_par'
      END IF

      CALL second0(ton)
      IF (lactive) THEN
         CALL ForwardSolve
      END IF

      CALL second0(time_off)
      toff = time_off
      bcyclic_forwardsolve_time = bcyclic_forwardsolve_time
     &                          + (toff - ton)

      IF (grank.EQ.0) THEN
         WRITE(6,1001) time_off - time_on
         WRITE(nthreed,1001) time_off - time_on
      END IF

      IF (.NOT.l_backslv) THEN
         DEALLOCATE (gc_save)
      END IF

      CALL second0(tprec2doff)

      timer(tprec2d) = timer(tprec2d) + (tprec2doff - tprec2don)
      compute_blocks_time = compute_blocks_time
     &                    + (tprec2doff - tprec2don)

1000  FORMAT(1x,' Time to compute blocks: ',f10.2,' s')
1001  FORMAT(1x,' Time to factor blocks:  ',f10.2,' s')

      END SUBROUTINE compute_blocks_par

      SUBROUTINE sweep3_blocks_par(xc, xcdot, gc)
      USE vmec_main, ONLY: ncurr, r01, z01, lthreed, chips, delt0r
      USE blocktridiagonalsolver, ONLY: Initialize, SetBlockRowCol
      USE blocktridiagonalsolver, ONLY: WriteBlocks
      USE parallel_vmec_module, ONLY: MPI_STAT
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,ntyptot) :: xc, xcdot, gc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: js, js1, istat, mesh, lamtype, rztype, icol
      INTEGER :: nsmin, nsmax
      INTEGER :: lastrank, left, right
      REAL(dp) :: eps, hj, hj_scale
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: diag_val
      REAL(dp) :: ton, toff
C-----------------------------------------------
!
!     COMPUTE FORCE "RESPONSE" TO PERTURBATION AT EVERY 3rd RADIAL POINT
!     FOR EACH MESH STARTING AT js=1,2,3, RESPECTIVELY
!
      CALL second0(ton)
      ALLOCATE(diag_val(ns), stat=istat)
      diag_val=zero
      IF (grank.EQ.0) THEN
        WRITE (6, *)
     &   " Using non-symmetric sweep to compute Hessian elements"
      END IF

      eps = SQRT(EPSILON(eps))
      eps = eps/10
      rztype = 2*ntmax
      lamtype = rztype + 1

      n_2d = 0
      m_2d = 0

      ALLOCATE(DataItem(0:ntor,0:mpol1,1:3*ntmax), stat=istat)
      IF (istat .ne. 0) THEN
         STOP 'Allocation error in sweep3_blocks'
      END IF

!
!     CALL FUNCT3D FIRST TIME TO STORE INITIAL UN-PRECONDITIONED FORCES
!     THIS WILL CALL LAMBLKS (SO FAR, ONLY IMPLEMENTED FOR lasym = false)
!
      ictrl_prec2d = 2                 !Signals funct3d that preconditioner is being initialized
      CALL funct3d_par(lscreen, istat)
      IF (istat .NE. 0) THEN
         PRINT *,' ier_flag = ', istat,
     1           ' in SWEEP3_BLOCKS_PAR call to funct3d_par'
         STOP
      ENDIF

      nsmin = t1lglob
      nsmax = t1rglob
      xcdot(:,:,nsmin:nsmax,:) = 0
!      PRINT *,'rank: ', rank,' vrank: ', vrank,' nsmin: ',nsmin,
!     1        ' nsmax: ', nsmax

      IF (FIRSTPASS) THEN
         edge_mesh = .FALSE.
         FIRSTPASS = .FALSE.
         mblk_size = (ntor + 1)*(mpol1 + 1)*3*ntmax
         IF (mblk_size .NE. ntmaxblocksize) THEN
            STOP 'wrong mblk_size in precon2d!'
         END IF
         CALL Initialize(.FALSE.,ns,mblk_size)
         myend = nsmax
!Align starting pt in (nsmin,nsmax)
         DO mesh = 1, 3
            icol = MOD(jstart(mesh) - nsmin, 3)
            IF (icol .LT. 0) THEN
               icol = icol + 3
            END IF
            mystart(mesh) = nsmin + icol
            IF (MOD(jstart(mesh) - ns, 3) .EQ. 0) THEN     ! .AND. nsmax.EQ.ns)
               edge_mesh(mesh) = .TRUE.
            END IF
         END DO
      END IF

!     STORE chips in xc 
#if defined(CHI_FORCE)
      IF (ncurr .EQ. 1) THEN
         xc(0,0,nsmin:nsmax,lamtype) = chips(nsmin:nsmax)
      END IF
#endif
      left = rank - 1
      IF (rank .EQ. 0) THEN
         left = MPI_PROC_NULL
      END IF
      right = rank + 1
      IF (rank .EQ. nranks - 1) THEN
         right = MPI_PROC_NULL
      END IF

      CALL PadSides(xc)
      CALL PadSides(gc)
      CALL PadSides(xcdot)

      CALL restart_iter(delt0r)
      gc_save(:,:,nsmin:nsmax,:) = gc(:,:,nsmin:nsmax,:) 
      ictrl_prec2d = 3                 !Signals funct3d that preconditioner is being computed

      CALL MPI_COMM_SIZE(NS_COMM,lastrank,MPI_ERR)
      lastrank = lastrank - 1
!
!     FIRST DO R00 JOG TO LOAD DIAG_VAL ARRAY (DO NOT RELY ON IT BEING THE FIRST JOG)
!
      m_2d=0
      n_2d=0
      ntype_2d = rcc

!     APPLY JOG
      hj = eps * MAX(ABS(r01(ns)), ABS(z01(ns)))
      IF (nranks .GT. 1) THEN
         CALL MPI_BCAST(hj,1,MPI_REAL8,lastrank,NS_COMM,MPI_ERR)
         CALL MPI_BCAST(edge_mesh,3,MPI_LOGICAL,lastrank,NS_COMM,              &
     &                  MPI_ERR)
      END IF
      DO js = mystart(1), myend(1), 3
         xcdot(n_2d,m_2d,js,ntype_2d) = hj
      END DO

      istat = 0
      CALL funct3d_par (lscreen, istat)

      LACTIVE0: IF (lactive) THEN
         IF (nranks .GT. 1) THEN
            CALL Gather4XArray(gc)
         END IF
         IF (istat .NE. 0) THEN
            STOP 'Error computing Hessian jog!'
         END IF
!     CLEAR JOG AND STORE BLOCKS FOR THIS JOG
         xcdot(:,:,nsmin:nsmax,:) = 0
         DO js = mystart(1), myend(1), 3
            DataItem = (gc(:,:,js,:) - gc_save(:,:,js,:))/hj
            diag_val(js) = DataItem(0,0,ntype_2d)
         END DO

         IF (nranks .GT. 1) THEN
            icol = 0
            IF (trglob_arr(1) .LT. 4) THEN
               icol = 1
            END IF
            CALL MPI_BCAST(diag_val(4), 1, MPI_REAL8, icol, NS_COMM,
     &                     MPI_ERR)
            icol = nranks - 1
            IF (tlglob_arr(nranks) .GT. ns - 3) THEN
               icol = nranks-2
            END IF
            CALL MPI_BCAST(diag_val(ns - 3), 1, MPI_REAL8, icol,
     &                     NS_COMM, MPI_ERR)
         END IF
         IF (diag_val(1) .EQ. zero) THEN
            diag_val(1)  = diag_val(4)
         END IF
         IF (diag_val(ns) .EQ. zero) THEN
            diag_val(ns) = diag_val(ns - 3)
         END IF

         IF (nranks .GT. 1) THEN
            CALL MPI_Sendrecv(diag_val(trglob), 1, MPI_REAL8, right, 1,
     &                        diag_val(t1lglob), 1, MPI_REAL8, left, 1,
     &                        NS_COMM, MPI_STAT, MPI_ERR)
         END IF
         DO js = mystart(2), myend(2), 3
            diag_val(js) = diag_val(js - 1)
         END DO

         hj_scale = MAX(ABS(r01(ns)), ABS(z01(ns)))

         IF (nranks .GT. 1) THEN
            CALL MPI_Sendrecv(diag_val(trglob), 1, MPI_REAL8, right, 1,
     &                        diag_val(t1lglob), 1, MPI_REAL8, left, 1,
     &                        NS_COMM, MPI_STAT, MPI_ERR)
            CALL MPI_BCAST(hj_scale, 1, MPI_REAL8, lastrank, NS_COMM,
     &                   MPI_ERR)
         END IF

         DO js = mystart(3), myend(3), 3
            diag_val(js) = diag_val(js - 1)
         END DO

         IF (ANY(diag_val(tlglob:trglob) .EQ. zero)) THEN
            PRINT *, 'For rank: ', rank, ' some diag_val == 0'
            STOP
         END IF
      END IF LACTIVE0
!
!     PERFORM "JOGS" FOR EACH VARIABLE AT EVERY 3rd RADIAL POINT ACROSS MESH
!     FOR ntyp = (Rcc, Rss, Rsc, Rcs, Zsc, Zcs, Zcc, Zss)
!     AND EVERY n2d (toroidal mode index) and EVERY m2d (poloidal mode index)

      icol=0

      NTYPE2D: DO ntype_2d = 1, ntyptot
         hj = eps
         IF (ntype_2d .LT. lamtype) THEN
            hj = hj*hj_scale
         END IF

         M2D: DO m_2d = 0, mpol1
            
            N2D: DO n_2d = 0, ntor

               icol = icol + 1
            
               MESH_3PT: DO mesh = 1,3
 
!              APPLY JOG TO ACTIVE PROCESSORS
                  IF (lactive) THEN
                     DO js = mystart(mesh), myend(mesh), 3
                        xcdot(n_2d,m_2d,js,ntype_2d) = hj
                     END DO
                     IF (m_2d.GT.0 .AND. mystart(mesh).EQ.1) THEN
                        xcdot(n_2d,m_2d,1,ntype_2d) = 0
                     END IF
                  END IF

                  l_edge = edge_mesh(mesh)
                  CALL funct3d_par (lscreen, istat)
                  IF (istat .NE. 0) STOP 'Error computing Hessian jog!'
              
!
!              COMPUTE PRECONDITIONER (HESSIAN) ELEMENTS. LINEARIZED EQUATIONS
!              OF FORM (FIXED mn FOR SIMPLICITY):
!
!              F(j-1) = a(j-1)x(j-2) + d(j-1)x(j-1) + b(j-1)x(j)
!              F(j)   =                a(j)x(j-1)   + d(j)  x(j) + b(j)  x(j+1)
!              F(j+1) =                               a(j+1)x(j) + d(j+1)x(j+1) + b(j+1)x(j+2)
!
!              HESSIAN IS H(k,j) == dF(k)/dx(j); aj == block_mins; dj == block_diag; bj = block_plus
!
!              THUS, A PERTURBATION (in xc) AT POSITION js PRODUCES THE FOLLOWING RESULTS:
!
!                     d(js)   = dF(js  )/hj(js)
!                     b(js-1) = dF(js-1)/hj(js)
!                     a(js+1) = dF(js+1)/hj(js)
!
!
                  LACTIVE1: IF (lactive) THEN
                     SKIP3_MESH: DO js = mystart(mesh), myend(mesh), 3

!                 CLEAR JOG AND STORE BLOCKS FOR THIS JOG
                        xcdot(n_2d,m_2d,js,ntype_2d) = 0

                  !block_mins(js+1) 
                       js1 = js+1
                        IF (tlglob.LE.js1 .AND. js1.LE.trglob) THEN
                           DataItem =
     &                        (gc(:,:,js1,:)-gc_save(:,:,js1,:))/hj
                           CALL SetBlockRowCol(js1,icol,DataItem,LOWER)
                        END IF

                  !block_diag(js)
                        IF (tlglob.LE.js .AND. js.LE.trglob) THEN
                           DataItem = (gc(:,:,js,:)
     &                              - gc_save(:,:,js,:))/hj

                           IF (rank .EQ. lastrank .AND.
     &                         js   .EQ. ns       .AND.
     &                         .NOT.lfreeb        .AND.
     &                         ANY(DataItem(:,:,
     &                                      1:rztype) .NE. zero)) THEN
                              STOP 'DIAGONAL BLOCK AT EDGE != 0'
                           END IF

!Levenberg-like offset - do NOT apply here if applied in colscaling routine
                           IF (ntype_2d .GE. lamtype) THEN
                              DataItem(n_2d,m_2d,ntype_2d) =
     &                           1.0001_dp*DataItem(n_2d,m_2d,ntype_2d)
                           END IF

                           IF (DataItem(n_2d,m_2d,
     &                                  ntype_2d) .EQ. zero) THEN
                              DataItem(n_2d,m_2d,ntype_2d) =
     &                           diag_val(js)
                           END IF

                           CALL SetBlockRowCol(js,icol,DataItem,DIAG)
                        END IF

                 !block_plus(js-1)
                        js1 = js - 1
                        IF (tlglob .LE. js1 .AND. js1 .LE. trglob) THEN
                           DataItem = (gc(:,:,js1,:) -
     &                                 gc_save(:,:,js1,:))/hj
                           CALL SetBlockRowCol(js1,icol,DataItem,UPPER)
                        END IF
                     END DO SKIP3_MESH
                  END IF LACTIVE1

               END DO MESH_3PT
            END DO N2D
         END DO M2D
      END DO NTYPE2D

      l_edge = .FALSE.

      DEALLOCATE(DataItem, diag_val)
      CALL second0(toff)
      fill_blocks_time=fill_blocks_time + (toff - ton)
 
      END SUBROUTINE sweep3_blocks_par


      SUBROUTINE compute_col_scaling_par
      USE xstuff, ONLY: pcol_scale
      USE blocktridiagonalsolver, ONLY: GetColSum, ParallelScaling
      USE parallel_vmec_module, ONLY: ToLastNtype, CopyLastNtype
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nsmin, nsmax
      REAL(dp), ALLOCATABLE  :: tmp(:)
      REAL(dp), ALLOCATABLE  :: colsum(:,:)
      REAL(dp), PARAMETER    :: levmarq_param = 1.E-6_dp
C-----------------------------------------------

!FOR NO COL SCALING - col-scaling not working well yet (8.1.17)
      pcol_scale = 1
      RETURN

!BE SURE TO TURN OFF LEV_MARQ SCALING IN SUBROUTINE sweep3_blocks_par
      nsmin = tlglob;  nsmax = trglob

      ALLOCATE (colsum(mblk_size,nsmin:nsmax))
      CALL GetColSum(colsum)
      CALL VectorCopyPar (colsum, pcol_scale)
      CALL ParallelScaling(levmarq_param,colsum)

      DEALLOCATE(colsum)

!Convert to internal PARVMEC format
      ALLOCATE (tmp(ntmaxblocksize*ns))
      CALL tolastntype(pcol_scale,tmp)
      CALL copylastntype(tmp,pcol_scale)
      DEALLOCATE(tmp)

      END SUBROUTINE compute_col_scaling_par

      SUBROUTINE VectorCopyPar (colsum, colscale)
      USE blocktridiagonalsolver, ONLY: rank
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN)  :: colsum(mblk_size,tlglob:trglob)
      REAL(dp), INTENT(OUT) :: colscale(mblk_size,ns)

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER               :: js, M1
      INTEGER               :: MPI_STAT(MPI_STATUS_SIZE)
!-----------------------------------------------
 
      DO js = tlglob, trglob 
        colscale(:,js) = colsum(:,js)
      END DO

      M1 = mblk_size

! Get left boundary elements (tlglob-1)
      IF (rank.LT.nranks-1) THEN
        CALL MPI_Send(colsum(:,trglob),M1,MPI_REAL8,    
     1                rank+1,1,NS_COMM,MPI_ERR)
      END IF
      IF (rank.GT.0) THEN
        CALL MPI_Recv(colscale(:,tlglob-1),M1,      
     1                MPI_REAL8,rank-1,1,NS_COMM,MPI_STAT,MPI_ERR)
      END IF

! Get right boundary elements (trglob+1)
      IF (rank.GT.0) THEN
        CALL MPI_Send(colsum(:,tlglob),M1,MPI_REAL8,
     1                rank-1,1,NS_COMM,MPI_ERR)
      END IF
      IF (rank.LT.nranks-1) THEN
        CALL MPI_Recv(colscale(:,trglob+1),M1,MPI_REAL8,
     1                rank+1,1,NS_COMM,MPI_STAT,MPI_ERR)
      END IF

      END SUBROUTINE VectorCopyPar


      SUBROUTINE free_mem_precon
      INTEGER :: istat

      istat=0
      IF (ALLOCATED(block_diag)) 
     1    DEALLOCATE (block_diag, block_plus, block_mins, stat=istat)
      IF (istat .ne. 0) STOP 'Deallocation error-1 in free_mem_precon'
      
      istat=0
      IF (ALLOCATED(block_diag_sw))
     1   DEALLOCATE (block_diag_sw, block_plus_sw, block_mins_sw, 
     2                stat=istat)
      IF (istat .ne. 0) STOP 'Deallocation error-2 in free_mem_precon'

      istat=0
      IF (ALLOCATED(ipiv_blk)) DEALLOCATE (ipiv_blk, stat=istat)     
      IF (istat .ne. 0) STOP 'Deallocation error-3 in free_mem_precon'

      END SUBROUTINE free_mem_precon

      END MODULE precon2d
