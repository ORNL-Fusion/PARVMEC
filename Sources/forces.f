!> \file forces.f

      SUBROUTINE forces_par
      USE vmec_main, p5 => cp5
      USE realspace
      USE vforces, ru12 => pazmn_e, zu12 => parmn_e,
     &             pazmn_e => pazmn_e, parmn_e => parmn_e,
     &             lv_e => pcrmn_e, lu_e => pczmn_e, lu_o => pczmn_o,
     &             pcrmn_e => pcrmn_e, pczmn_e => pczmn_e,
     &             pczmn_o => pczmn_o
      USE parallel_include_module
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(dp), PARAMETER :: p25 = p5*p5, dshalfds=p25
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l, ndim
      INTEGER :: i, j, k, nsmin, nsmax
      REAL(dp), DIMENSION(:,:), POINTER ::
     &    bsqr, gvvs, guvs, guus
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: bcastbuf
C-----------------------------------------------
      IF (.NOT.lactive .AND. .NOT.lfreeb) RETURN

      CALL second0 (tforon)

      ndim = 1+nrzt

      nsmin=tlglob; nsmax=t1rglob

!     POINTER ALIASES
      bsqr => pextra1(:,:,1);  gvvs => pextra2(:,:,1)
      guvs => pextra3(:,:,1);  guus => pextra4(:,:,1)

      lu_e(:,1) = 0; lv_e(:,1) = 0
      pguu(:,1)  = 0; pguv(:,1)  = 0; pgvv(:,1) = 0

      DO l = nsmin, nsmax
         guus(:,l) = pguu(:,l)*pshalf(:,l)
         guvs(:,l) = pguv(:,l)*pshalf(:,l)
         gvvs(:,l) = pgvv(:,l)* pshalf(:,l)

         parmn_e(:,l) = ohs*zu12(:,l)*lu_e(:,l)
         pazmn_e(:,l) =-ohs*ru12(:,l)*lu_e(:,l)
         pbrmn_e(:,l) = pbrmn_e(:,l)*lu_e(:,l)
         pbzmn_e(:,l) =-pbzmn_e(:,l)*lu_e(:,l)
         bsqr(:,l)    = dshalfds*lu_e(:,l)/pshalf(:,l)

         parmn_o(:,l) = parmn_e(:,l)*pshalf(:,l)
         pazmn_o(:,l) = pazmn_e(:,l)*pshalf(:,l)
         pbrmn_o(:,l) = pbrmn_e(:,l)*pshalf(:,l)
         pbzmn_o(:,l) = pbzmn_e(:,l)*pshalf(:,l)
      END DO

!
!     CONSTRUCT CYLINDRICAL FORCE KERNELS
!     NOTE: presg(ns+1) == 0, AND WILL BE "FILLED IN" AT EDGE
!     FOR FREE-BOUNDARY BY RBSQ
!
!DIR$ IVDEP
      nsmin=tlglob; nsmax=MIN(ns-1,trglob)
      DO l = nsmin, nsmax
         pguu(:,l) = p5*(pguu(:,l) + pguu(:,l+1))
         pgvv(:,l) = p5*(pgvv(:,l) + pgvv(:,l+1))
         bsqr(:,l) = bsqr(:,l) +     bsqr(:,l+1)
         guus(:,l) = p5*(guus(:,l) + guus(:,l+1))
         gvvs(:,l) = p5*(gvvs(:,l) + gvvs(:,l+1))
      END DO

      IF (trglob .ge. ns) THEN
         pguu(:,ns) = p5*pguu(:,ns)
         pgvv(:,ns) = p5*pgvv(:,ns)
         guus(:,ns) = p5*guus(:,ns)
         gvvs(:,ns) = p5*gvvs(:,ns)
      END IF

!DIR$ IVDEP
      nsmin=tlglob; nsmax=MIN(ns-1,trglob)
      DO l = nsmin, nsmax
         parmn_e(:,l) = parmn_e(:,l+1) - parmn_e(:,l)
     &                + p5*(lv_e(:,l) + lv_e(:,l+1))
         pazmn_e(:,l) = pazmn_e(:,l+1) - pazmn_e(:,l)
         pbrmn_e(:,l) = p5*(pbrmn_e(:,l) + pbrmn_e(:,l+1))
         pbzmn_e(:,l) = p5*(pbzmn_e(:,l) + pbzmn_e(:,l+1))
      END DO

      parmn_e(:,ns) = - parmn_e(:,ns) + p5*lv_e(:,ns)
      pazmn_e(:,ns) = - pazmn_e(:,ns)
      pbrmn_e(:,ns) = p5*pbrmn_e(:,ns)
      pbzmn_e(:,ns) = p5*pbzmn_e(:,ns)

      nsmin=tlglob; nsmax=t1rglob
      DO l = nsmin, nsmax
         parmn_e(:,l) = parmn_e(:,l)
     &                - (gvvs(:,l)*pr1(:,l,1) + pgvv(:,l)*pr1(:,l,0))
         pbrmn_e(:,l) = pbrmn_e(:,l) + bsqr(:,l)*pz1(:,l,1)
     &                - (guus(:,l)*pru(:,l,1) + pguu(:,l)*pru(:,l,0))
         pbzmn_e(:,l) = pbzmn_e(:,l) - (bsqr(:,l)*pr1(:,l,1)
     &                +  guus(:,l)*pzu(:,l,1) + pguu(:,l)*pzu(:,l,0))
         lv_e(:,l) = lv_e(:,l)*pshalf(:,l)
         lu_o(:,l) = dshalfds*lu_e(:,l)
      END DO

      nsmin=tlglob; nsmax=MIN(ns-1,trglob)
!DIR$ IVDEP
      DO l = nsmin, nsmax
         parmn_o(:,l) = parmn_o(:,l+1) - parmn_o(:,l)
     &                - pzu(:,l,0)*bsqr(:,l)
     &                + p5*(lv_e(:,l)+lv_e(:,l+1))
         pazmn_o(:,l) = pazmn_o(:,l+1) - pazmn_o(:,l)
     &                + pru(:,l,0)*bsqr(:,l)
         pbrmn_o(:,l) = p5*(pbrmn_o(:,l) +  pbrmn_o(:,l+1))
         pbzmn_o(:,l) = p5*(pbzmn_o(:,l) +  pbzmn_o(:,l+1))
         lu_o(:,l)    = lu_o(:,l) + lu_o(:,l+1)
      END DO

      parmn_o(:,ns) = - parmn_o(:,ns) - pzu(:,ns,0)*bsqr(:,ns)
     &                + p5*lv_e(:,ns)
      pazmn_o(:,ns) = - pazmn_o(:,ns) + pru(:,ns,0)*bsqr(:,ns)
      pbrmn_o(:,ns) = p5*pbrmn_o(:,ns)
      pbzmn_o(:,ns) = p5*pbzmn_o(:,ns)
      lu_o(:,ns)   = lu_o(:,ns)

      nsmin=tlglob; nsmax=trglob
      DO l = nsmin, nsmax
         pguu(:,l) = pguu(:,l) * psqrts(:,l)**2
         bsqr(:,l) = pgvv(:,l) * psqrts(:,l)**2
      END DO

      DO l = nsmin, nsmax
         parmn_o(:,l) = parmn_o(:,l) - (pzu(:,l,1)*lu_o(:,l)
     &                + bsqr(:,l)*pr1(:,l,1) + gvvs(:,l)*pr1(:,l,0))
         pazmn_o(:,l) = pazmn_o(:,l) +  pru(:,l,1)*lu_o(:,l)
         pbrmn_o(:,l) = pbrmn_o(:,l) +  pz1(:,l,1)*lu_o(:,l)
     &                -(pguu(:,l)*pru(:,l,1) + guus(:,l)*pru(:,l,0))
         pbzmn_o(:,l) = pbzmn_o(:,l) - (pr1(:,l,1)*lu_o(:,l)
     &                + pguu(:,l)*pzu(:,l,1) + guus(:,l)*pzu(:,l,0))
      END DO

      IF (lthreed) THEN
!DIR$ IVDEP
         nsmin=tlglob; nsmax=MIN(ns-1,trglob)
         DO l = nsmin, nsmax
            pguv(:,l)  = p5*(pguv(:,l) + pguv(:,l+1))
            guvs(:,l) = p5*(guvs(:,l) + guvs(:,l+1))
         END DO
         pguv(:,ns) = p5*pguv(:,ns)
         guvs(:,ns) = p5*guvs(:,ns)

         nsmin=tlglob; nsmax=trglob
         DO l = nsmin, nsmax
            pbrmn_e(:,l) = pbrmn_e(:,l)
     &                   - (pguv(:,l)*prv(:,l,0) + guvs(:,l)*prv(:,l,1))
            pbzmn_e(:,l) = pbzmn_e(:,l)
     &                   - (pguv(:,l)*pzv(:,l,0) + guvs(:,l)*pzv(:,l,1))
            pcrmn_e(:,l) = pguv(:,l)*pru(:,l,0) + pgvv(:,l)*prv(:,l,0)
     &                   + gvvs(:,l)*prv(:,l,1) + guvs(:,l)*pru(:,l,1)
            pczmn_e(:,l) = pguv(:,l)*pzu(:,l,0) + pgvv(:,l)*pzv(:,l,0)
     &                   + gvvs(:,l)*pzv(:,l,1) + guvs(:,l)*pzu(:,l,1)
            pguv(:,l) = pguv(:,l)*psqrts(:,l)*psqrts(:,l)
            pbrmn_o(:,l) = pbrmn_o(:,l)
     &                   - (guvs(:,l)*prv(:,l,0) + pguv(:,l)*prv(:,l,1))
            pbzmn_o(:,l) = pbzmn_o(:,l)
     &                   - (guvs(:,l)*pzv(:,l,0) + pguv(:,l)*pzv(:,l,1))
            pcrmn_o(:,l) = guvs(:,l)*pru(:,l,0) + gvvs(:,l)*prv(:,l,0)
     &                   + bsqr(:,l)*prv(:,l,1) + pguv(:,l)*pru(:,l,1)
             pczmn_o(:,l) = guvs(:,l)*pzu(:,l,0) + gvvs(:,l)*pzv(:,l,0)
     &                    + bsqr(:,l)*pzv(:,l,1) + pguv(:,l)*pzu(:,l,1)
         END DO
      ENDIF

!
!     ASSIGN EDGE FORCES (JS = NS) FOR FREE BOUNDARY CALCULATION
!
      IF (ivac .GE. 1) THEN
         DO k = 1, ntheta3
            DO j = 1, nzeta
               l = (k-1)*nzeta + j
               parmn_e(l,ns) = parmn_e(l,ns) + pzu0(l,ns)*rbsq(l)
               parmn_o(l,ns) = parmn_o(l,ns) + pzu0(l,ns)*rbsq(l)
               pazmn_e(l,ns) = pazmn_e(l,ns) - pru0(l,ns)*rbsq(l)
               pazmn_o(l,ns) = pazmn_o(l,ns) - pru0(l,ns)*rbsq(l)
            END DO
         END DO
      ENDIF

!
!     COMPUTE CONSTRAINT FORCE KERNELS
!
      DO l = nsmin, nsmax
         prcon(:,l,0) = (prcon(:,l,0)-prcon0(:,l))*pgcon(:,l)
         pzcon(:,l,0) = (pzcon(:,l,0)-pzcon0(:,l))*pgcon(:,l)
         pbrmn_e(:,l) = pbrmn_e(:,l) + prcon(:,l,0)
         pbzmn_e(:,l) = pbzmn_e(:,l) + pzcon(:,l,0)
         pbrmn_o(:,l) = pbrmn_o(:,l)+ prcon(:,l,0)*psqrts(:,l)
         pbzmn_o(:,l) = pbzmn_o(:,l)+ pzcon(:,l,0)*psqrts(:,l)
         prcon(:,l,0) = pru0(:,l) * pgcon(:,l)
         pzcon(:,l,0) = pzu0(:,l) * pgcon(:,l)
         prcon(:,l,1) = prcon(:,l,0) * psqrts(:,l)
         pzcon(:,l,1) = pzcon(:,l,0) * psqrts(:,l)
      END DO

      CALL second0 (tforoff)
      timer(tfor) = timer(tfor) + (tforoff - tforon)
      forces_time = timer(tfor)

      END SUBROUTINE forces_par
