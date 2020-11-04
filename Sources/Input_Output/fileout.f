!> \file fileout.f

      SUBROUTINE fileout_par(iseq, ictrl_flag, ier_flag, lscreen)
      USE vmec_main
      USE parallel_include_module
      USE xstuff, ONLY: pxc, pgc, pxsave, pscalxc,
     &                  xc, gc, xsave, scalxc
      USE vmec_params, ONLY: uminus, output_flag, cleanup_flag,
     &                       more_iter_flag, norm_term_flag,
     &                       signgs, successful_term_flag,
     &                       mscale, nscale
      USE realspace
      USE vacmod, ONLY: bsqvac, brv, bphiv, bzv, nv, nuv3,
     &                  bsupu_sur, bsupv_sur, bsubu_sur, bsubv_sur
      USE vac_persistent
      USE vforces
      USE precon2d, ONLY: ictrl_prec2d
      USE timer_sub
#ifdef _HBANGLE
      USE angle_constraints, ONLY: free_multipliers, getrz
#endif
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: iseq, ictrl_flag
      INTEGER, INTENT(inout) :: ier_flag
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      LOGICAL :: loutput !SAL 070719
      INTEGER :: i, j, ij, k, js, jk, lk, lt, lz, jcount
      REAL(dp) :: tfileon, tfileoff
      REAL(dp), ALLOCATABLE :: buffer(:,:), tmp(:,:)
      integer :: irst0, istat
      LOGICAL :: lwrite, lterm
      integer :: loc_ier_flag
      integer:: istat1=0
      REAL(dp), ALLOCATABLE :: br_out(:), bz_out(:)
      CHARACTER(LEN=*), PARAMETER, DIMENSION(0:14) :: werror = (/
     &   'EXECUTION TERMINATED NORMALLY                            ', ! norm_term_flag
     &   'INITIAL JACOBIAN CHANGED SIGN (IMPROVE INITIAL GUESS)    ', ! bad_jacobian_flag
     &   'FORCE RESIDUALS EXCEED FTOL: MORE ITERATIONS REQUIRED    ', ! more_iter_flag
     &   'VMEC INDATA ERROR: NCURR.ne.1 but BLOAT.ne.1.            ', !
     &   'MORE THAN 75 JACOBIAN ITERATIONS (DECREASE DELT)         ', ! jac75_flag
     &   'ERROR READING INPUT FILE OR NAMELIST                     ', ! input_error_flag
     &   'NEW AXIS GUESS STILL FAILED TO GIVE GOOD JACOBIAN        ', !
     &   'PHIEDGE HAS WRONG SIGN IN VACUUM SUBROUTINE              ', ! phiedge_error_flag
     &   'NS ARRAY MUST NOT BE ALL ZEROES                          ', ! ns_error_flag
     &   'ERROR READING MGRID FILE                                 ', ! misc_error_flag
     &   'VAC-VMEC I_TOR MISMATCH : BOUNDARY MAY ENCLOSE EXT. COIL ', !
     &   'SUCCESSFUL VMEC CONVERGENCE                              ', ! successful_term_flag
     &   'BSUBU OR BSUBV JS=1 COMPONENT NON-ZERO                   ', ! bsub_bad_js1_flag
     &   'RMNC N=0, M=1 IS ZERO                                    ', ! r01_bad_value_flag
     &   'ARNORM OR AZNORM EQUAL ZERO IN BCOVAR                    '  ! arz_bad_value_flag
     &   /)
      CHARACTER(LEN=*), PARAMETER ::
     &    Warning = " Error deallocating global memory FILEOUT"
      REAL(dp) :: tmpxc, rmssum
C-----------------------------------------------
      CALL second0(tfileon)

      istat = norm_term_flag ! jons: define initial value for istat
      ! or better ier_flag ? but then bcovar would not successfully execute last branch,
      ! where ier_flag.ne.norm_term_flag is checked...
      ! or modify bcovar to also accept successful_term_flag ?

      INFILEOUT=.TRUE.

      loutput = (IAND(ictrl_flag, output_flag) .ne. 0) !SAL 070719

      IF (lfreeb .AND. vlactive .AND. loutput) THEN !SAL070719
         ALLOCATE(buffer(numjs_vac, 7), tmp(nznt, 7),stat=i)
         IF (i .NE. 0) CALL STOPMPI(440)
         buffer(:,1) = brv(nuv3min:nuv3max)
         buffer(:,2) = bphiv(nuv3min:nuv3max)
         buffer(:,3) = bzv(nuv3min:nuv3max)
         buffer(:,4) = bsupu_sur(nuv3min:nuv3max)
         buffer(:,5) = bsupv_sur(nuv3min:nuv3max)
         buffer(:,6) = bsubu_sur(nuv3min:nuv3max)
         buffer(:,7) = bsubv_sur(nuv3min:nuv3max)

         DO i = 1, 7
            CALL MPI_Allgatherv(buffer(:,i), numjs_vac, MPI_REAL8,
     &                          tmp(:,i), counts_vac, disps_vac,
     &                          MPI_REAL8, VAC_COMM,MPI_ERR)
         END DO
         DEALLOCATE(buffer)

         brv = tmp(:,1);
         bphiv = tmp(:,2);
         bzv = tmp(:,3)
         bsupu_sur = tmp(:,4);
         bsupv_sur = tmp(:,5)
         bsubu_sur = tmp(:,6);
         bsubv_sur = tmp(:,7)
         DEALLOCATE(tmp)
      END IF
!
!     COMPUTE ARRAY FOR REFLECTING v = -v (ONLY needed for lasym)
!
      jcount = 0
      DO k = 1, nzeta
         jk = nzeta + 2 - k
         IF (k .eq. 1) jk = 1
         DO js = 1,ns
            jcount = jcount + 1
            ireflect(jcount) = js+ns*(jk - 1)           !Index for -zeta[k]
         ENDDO
      END DO

!     INDEX FOR u = -u (need for lasym integration in wrout)
      lk = 0
      IF (.NOT.ALLOCATED(uminus)) ALLOCATE (uminus(nznt))
      DO lt = 1, ntheta2
         k = ntheta1-lt+2                  
         IF (lt .eq. 1) k = 1             !u=-0 => u=0
         DO lz = 1, nzeta
            lk = lk + 1
            uminus(lk) = k                !(-u), for u = 0,pi
         END DO
      END DO

!
!     COMPUTE REMAINING COVARIANT COMPONENT OF B (BSUBS),
!     CYLINDRICAL COMPONENTS OF B (BR, BPHI, BZ), AND
!     AVERAGE EQUILIBRIUM PROPERTIES AT END OF RUN
!

      iequi = 1
      lterm = ier_flag .eq. norm_term_flag  .or.
     &        ier_flag .eq. successful_term_flag
      lwrite = lterm .or. ier_flag.eq.more_iter_flag
      loutput = (IAND(ictrl_flag, output_flag) .ne. 0)
      loc_ier_flag = ier_flag
      if (ier_flag .eq. successful_term_flag) THEN
          loc_ier_flag = norm_term_flag
      end if

      IF (lwrite .AND. loutput) THEN
!
!     The sign of the jacobian MUST multiply phi to get the physically
!     correct toroidal flux
!
         phi(1) = zero
         DO js = 2, ns
            phi(js) = phi(js-1) + phip(js)
         END DO
         phi = (signgs*twopi*hs)*phi

!        Must save irst value if in "restart" mode
         irst0 = irst
         CALL funct3d_par (lscreen, istat)
         fo_funct3d_time = timer(tfun)

!        Write out any special files here
!        CALL dump_special

         irst = irst0
      END IF
      
      IF (grank.LT.nranks .AND. loutput) THEN
         CALL Gather1XArray(vp)
         CALL Gather1XArray(iotas)
         CALL Gather1XArray(phips)
         CALL Gather1XArray(phipf)
         CALL Gather1XArray(chips)
         CALL Gather1XArray(mass)
         CALL Gather1XArray(icurv)
         CALL Gather1XArray(specw)
         CALL Gather1XArray(bdamp)
         CALL Gather1XArray(sm)
         CALL Gather1XArray(sp)

         CALL Gather2XArray(pphip)
         CALL Parallel2Serial2X(pphip, phip)
         CALL Gather2XArray(psqrts)
         CALL Parallel2Serial2X(psqrts, sqrts)
         CALL Gather2XArray(pshalf)
         CALL Parallel2Serial2X(pshalf, shalf)
         CALL Gather2XArray(pwint)
         CALL Parallel2Serial2X(pwint, wint)

         CALL Gather4XArray(pxc)
         CALL Parallel2Serial4X(pxc,xc)
         CALL Gather4XArray(pscalxc)
         CALL Parallel2Serial4X(pscalxc,scalxc)

         ! vforces
         CALL Gather2XArray(parmn)
         CALL Parallel2Serial2X(parmn, armn)
         CALL Gather2XArray(pazmn)
         CALL Parallel2Serial2X(pazmn, azmn)
         CALL Gather2XArray(pbrmn)
         CALL Parallel2Serial2X(pbrmn, brmn)
         CALL Gather2XArray(pbzmn)
         CALL Parallel2Serial2X(pbzmn, bzmn)
         CALL Gather2XArray(pcrmn)
         CALL Parallel2Serial2X(pcrmn, crmn)
         CALL Gather2XArray(pczmn)
         CALL Parallel2Serial2X(pczmn, czmn)

         CALL Gather2XArray(pblmn)
         CALL Parallel2Serial2X(pblmn, blmn)
         CALL Gather2XArray(pclmn)
         CALL Parallel2Serial2X(pclmn, clmn)

         ! even and odd separately to fill separate serial arrays
         CALL Gather2XArray(parmn_o)
         CALL Parallel2Serial2X(parmn_o, armn_o)
         CALL Gather2XArray(parmn_e)
         CALL Parallel2Serial2X(parmn_e, armn_e)
         CALL Gather2XArray(pazmn_o)
         CALL Parallel2Serial2X(pazmn_o, azmn_o)
         CALL Gather2XArray(pazmn_e)
         CALL Parallel2Serial2X(pazmn_e, azmn_e)

         CALL Gather2XArray(pbrmn_o)
         CALL Parallel2Serial2X(pbrmn_o, brmn_o)
         CALL Gather2XArray(pbrmn_e)
         CALL Parallel2Serial2X(pbrmn_e, brmn_e)
         CALL Gather2XArray(pbzmn_o)
         CALL Parallel2Serial2X(pbzmn_o, bzmn_o)
         CALL Gather2XArray(pbzmn_e)
         CALL Parallel2Serial2X(pbzmn_e, bzmn_e)

         CALL Gather2XArray(pcrmn_o)
         CALL Parallel2Serial2X(pcrmn_o, crmn_o)
         CALL Gather2XArray(pcrmn_e)
         CALL Parallel2Serial2X(pcrmn_e, crmn_e)
         CALL Gather2XArray(pczmn_o)
         CALL Parallel2Serial2X(pczmn_o, czmn_o)
         CALL Gather2XArray(pczmn_e)
         CALL Parallel2Serial2X(pczmn_e, czmn_e)

         ! probably not needed since blmn and clmn
         ! as a whole contain bsubu and bsubv here
         CALL Gather2XArray(pblmn_o)
         CALL Parallel2Serial2X(pblmn_o, blmn_o)
         CALL Gather2XArray(pblmn_e)
         CALL Parallel2Serial2X(pblmn_e, blmn_e)
         CALL Gather2XArray(pclmn_o)
         CALL Parallel2Serial2X(pclmn_o, clmn_o)
         CALL Gather2XArray(pclmn_e)
         CALL Parallel2Serial2X(pclmn_e, clmn_e)
 
         ! needed implicitly for bss
         CALL Gather2XArray(prv(:,:,0))
         CALL Gather2XArray(prv(:,:,1))
         CALL Parallel2Serial2X(prv(:,:,0), rv(:,0))
         CALL Parallel2Serial2X(prv(:,:,1), rv(:,1))
         CALL Gather2XArray(pzv(:,:,0))
         CALL Gather2XArray(pzv(:,:,1))
         CALL Parallel2Serial2X(pzv(:,:,0), zv(:,0))
         CALL Parallel2Serial2X(pzv(:,:,1), zv(:,1))

         call gather2xarray(pru0)
         call parallel2serial2x(pru0, ru0)
         call gather2xarray(pzu0)
         call parallel2serial2x(pzu0, zu0)

         ! jperp2
         call gather2xarray(pguu)
         call parallel2serial2x(pguu, guu)
         call gather2xarray(pguv)
         call parallel2serial2x(pguv, guv)
         call gather2xarray(pgvv)
         call parallel2serial2x(pgvv, gvv)

         ! for equif
         call Gather1XArray(vpphi)
         call Gather1XArray(jcurv)
         call Gather1XArray(chipf)
         call Gather1XArray(jcuru)
         call Gather1XArray(phipf)
         call Gather1XArray(pres)
         call Gather1XArray(presgrad)

         call Gather1XArray(buco)
         call Gather1XArray(bvco)

         call Gather1XArray(equif)

         CALL Gather1XArray(iotaf)

         CALL second0(tfileoff)
      END IF
      fo_prepare_time = fo_prepare_time + (tfileoff-tfileon)

      IF (grank .EQ. 0) THEN
         IF (lwrite .AND. loutput) THEN
           ALLOCATE(br_out(nrzt), bz_out(nrzt), stat=istat)
           gc = xc
#ifdef _HBANGLE
           CALL getrz(gc)
#endif
           CALL eqfor(br_out, bz_out, clmn, blmn, rcon(1,1),
     &                gc, ier_flag)
         END IF
         
         
         
!         CALL free_mem_precon
!        
!        Call WROUT to write output or error message if lwrite = false
!        
         
         IF (loutput .AND. ASSOCIATED(bzmn_o)) THEN
            CALL wrout(bzmn_o, azmn_o, clmn, blmn, crmn_o, czmn_e,
     &                 crmn_e, xsave, gc, loc_ier_flag, lwrite
#ifdef _ANIMEC
     &                 ,brmn_o, sigma_an, ppar, pperp, onembc, pp1, pp2,
     &                 pp3
#endif    
     &                )
         
            IF (ntor .EQ. 0) THEN
               CALL write_dcon (xc)
            END IF
         
            IF (lscreen .and. ier_flag.ne.more_iter_flag) 
     &      PRINT 120, TRIM(werror(loc_ier_flag))
            IF (lscreen .and. lterm) THEN
               IF (grank.EQ.0) THEN
                  PRINT 10, TRIM(input_extension), ijacob
               END IF
            END IF
         
            IF (nthreed .gt. 0) THEN
               WRITE (nthreed,120) TRIM(werror(loc_ier_flag))
               IF (.not. lterm) GOTO 1000
               WRITE (nthreed, 10) TRIM(input_extension), ijacob
               IF (rank.EQ.0) THEN
                  CALL write_times(nthreed, lscreen, lfreeb, lrecon,
     &                             ictrl_prec2d .ne. 0)
         
               IF (grank.EQ.0) THEN 
                 WRITE(nthreed,*)
                 WRITE(nthreed,'(1x,a,i4)') 'NO. OF PROCS:  ',gnranks
                 WRITE(nthreed,101)         'LPRECOND    :  ',LPRECOND
                 WRITE(nthreed,101)         'LV3FITCALL  :  ',LV3FITCALL
               END IF
 101     FORMAT(1x,a,l4)
               END IF
            END IF
         END IF
         
   10    FORMAT(' FILE : ',a,/,' NUMBER OF JACOBIAN RESETS = ',i4,/)
  120    FORMAT(/1x,a,/)
         
!        
!        TESTING READ_WOUT MODULE WRITING ROUTINES
!        
         IF (ALLOCATED(br_out)) THEN
!           IF (lscreen) CALL TestWout(xc, br_out, bz_out, crmn_e, czmn_e)
            DEALLOCATE (br_out, bz_out)
         END IF
         
!        END TEST
         
 1000    CONTINUE
         
!        
!        DEALLOCATE GLOBAL MEMORY AND CLOSE FILES
!        
         IF (IAND(ictrl_flag, cleanup_flag) .eq. 0 .or.
     &       ier_flag                       .eq. more_iter_flag) THEN
            RETURN
         END IF
         
         IF (ALLOCATED(cosmu))
     &     DEALLOCATE(cosmu, sinmu, cosmum, sinmum, cosmui, cosmumi,
     &                sinmui, sinmumi, cosnv, sinnv, cosnvn, sinnvn,
     &                cosmui3, cosmumi3, cos01, sin01, stat=istat1)
         IF (istat1 .ne. 0) PRINT *, Warning // "#1"
#ifdef _HBANGLE
         CALL free_multipliers
#endif   
         
         IF (ALLOCATED(xm)) DEALLOCATE (xm, xn, ixm, xm_nyq, xn_nyq, 
     &      jmin3, mscale, nscale, uminus, stat=istat1)
         IF (istat1 .ne. 0) PRINT *, Warning // "#2"
         
         IF (ALLOCATED(tanu))
     &      DEALLOCATE(tanu, tanv, sinper, cosper, sinuv, cosuv, sinu,
     &               cosu, sinv, cosv, sinui, cosui, cmns, csign, sinu1,
     &               cosu1, sinv1, cosv1, imirr, xmpot, xnpot,
     &               stat=istat1)
         IF (istat1 .ne. 0) PRINT *, Warning // "#3"
         
         CALL close_all_files
      ENDIF

      CALL free_persistent_mem

      CALL free_mem_funct3d_par
      CALL free_mem_ns_par(.false.)
      CALL free_mem_nunv
      
      !CALL MPI_Barrier(NS_COMM, MPI_ERR) !SAL 070719
      CALL second0(tfileoff)
      fileout_time = fileout_time + (tfileoff-tfileon)
      fo_par_call_time = fileout_time

      END SUBROUTINE fileout_par
