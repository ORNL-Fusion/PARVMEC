!> \file jacobian.f

      SUBROUTINE jacobian_par
      USE vmec_input, ONLY: nzeta
      USE vmec_main, ONLY: ohs, nrzt, irst, nznt, iter2
      USE vmec_params, ONLY: meven, modd
      USE realspace
      USE vmec_dim, ONLY: ns, ntheta3
      USE vforces, pr12 => parmn_o, pzu12 => parmn_e, pru12 => pazmn_e,
     &             prs => pbzmn_e, pzs => pbrmn_e, ptau => pazmn_o
      USE parallel_include_module

      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(dp), PARAMETER :: zero = 0, p5 = 0.5_dp, p25 = p5*p5
      REAL(dp) :: dphids
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, nsmin, nsmax, lnsnum
      REAL(dp) :: ltaumax, ltaumin
      REAL(dp) :: taumax, taumin
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: temp
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: minarr, maxarr
      REAL(dp) :: t1, t2, tjacon, tjacoff
C-----------------------------------------------
      CALL second0(tjacon)

      nsmin=MAX(2,tlglob); nsmax=t1rglob;
      !nsmin=tlglob; nsmax=t1rglob; ! jons; this fixes some issues in threed1, but leads to segfault for solovev2 input at ns=64 multigrid step
      dphids = p25
      irst = 1

      DO i = nsmin, nsmax
         pru12(:,i) = p5*(pru(:,i,meven) + pru(:,i-1,meven) +
     &                    pshalf(:,i)*(pru(:,i,modd) +
     &                                 pru(:,i-1,modd)))
         pzs(:,i)   = ohs*(pz1(:,i,meven) - pz1(:,i-1,meven) +
     &                     pshalf(:,i)*(pz1(:,i,modd) -
     &                                  pz1(:,i-1,modd)))
         ptau(:,i) = pru12(:,i)*pzs(:,i)
     &             + dphids*(pru(:,i,modd)*pz1(:,i,modd) +
     &                       pru(:,i-1,modd)*pz1(:,i-1,modd) +
     &                       (pru(:,i,meven)*pz1(:,i,modd) +
     &                        pru(:,i-1,meven)*pz1(:,i-1,modd)) /
     &                       pshalf(:,i))
      END DO

      DO i = nsmin, nsmax
         pzu12(:,i) = p5*(pzu(:,i,meven) + pzu(:,i-1,meven) +
     &                    pshalf(:,i)*(pzu(:,i,modd) +
     &                                 pzu(:,i-1,modd)))
         prs(:,i)   = ohs*(pr1(:,i,meven) - pr1(:,i-1,meven) +
     &                     pshalf(:,i)*(pr1(:,i,modd) -
     &                                  pr1(:,i-1,modd)))
         pr12(:,i)  = p5*(pr1(:,i,meven) + pr1(:,i-1,meven) +
     &                    pshalf(:,i)*(pr1(:,i,modd) +
     &                                 pr1(:,i-1,modd)))
         ptau(:,i) = ptau(:,i) - prs(:,i)*pzu12(:,i)
     &             - dphids*(pzu(:,i,modd)*pr1(:,i,modd) +
     &                       pzu(:,i-1,modd)*pr1(:,i-1,modd) +
     &                       (pzu(:,i,meven)*pr1(:,i,modd) +
     &                        pzu(:,i-1,meven)*pr1(:,i-1,modd)) /
     &                       pshalf(:,i))
      END DO

      ALLOCATE(temp(1:nznt))
      temp(:)=ptau(:,2)
      ptau(:,1)=temp(:)
      DEALLOCATE(temp)

      ltaumax=MAXVAL(ptau(:,nsmin:nsmax))
      ltaumin=MINVAL(ptau(:,nsmin:nsmax))
!      ltaumax=MAXVAL(ptau(:,tlglob:trglob))
!      ltaumin=MINVAL(ptau(:,tlglob:trglob))

      taumax=ltaumax
      taumin=ltaumin

      IF (nranks.GT.1.AND.grank.LT.nranks) THEN
         CALL second0(t1)
         CALL MPI_Allreduce(ltaumax,taumax,1,MPI_REAL8,
     &                      MPI_MAX,NS_COMM,MPI_ERR)
         CALL MPI_Allreduce(ltaumin,taumin,1,MPI_REAL8,
     &                      MPI_MIN,NS_COMM,MPI_ERR)
         CALL second0(t2)
         allreduce_time = allreduce_time + (t2-t1)
      END IF

      IF (taumax*taumin .lt. zero) THEN
         irst = 2
      END IF

      CALL second0(tjacoff)
      jacobian_time=jacobian_time+(tjacoff-tjacon)

      END SUBROUTINE jacobian_par
