!> \file wrout.f

      SUBROUTINE wrout(bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu,
     1                 rzl_array, gc_array, ier_flag)
! ... from SPH 2009-10-05; changes for modB sine-harmonics included
      USE vmec_input, ONLY: ns_array, ftol_array
      USE vmec_params
      USE vmec_main
      USE vmercier
      USE vmec_persistent
      USE vparams, p5 => cp5, two => c2p0
      USE vac_persistent
      USE xstuff
      USE vmec_io
      USE realspace, ONLY: phip, chip, gsqrta=>z1, z1=>z1
      USE totzsp_mod
      USE vforces, ONLY: bsupua=>brmn_e, bsupva=>czmn_o, bsqa=>bzmn_e,
     1                   bsubsa=>armn_e, bsubua=>azmn_e, bsubva=>armn_o
      USE vacmod, ONLY: potvac, mnpd, xmpot, xnpot

      USE ezcdf
      USE read_wout_mod, ONLY: Compute_Currents,
	1  vn_version, vn_extension, vn_mgrid,
     1  vn_magen, vn_therm, vn_gam, vn_maxr, vn_minr, vn_maxz, vn_fp,
     2  vn_radnod, vn_polmod, vn_tormod, vn_maxmod, vn_maxit, vn_actit,
     3  vn_asym, vn_recon, vn_free, vn_error, vn_aspect, vn_beta,
     4  vn_pbeta, vn_tbeta, vn_abeta, vn_b0, vn_rbt0, vn_maxmod_nyq,
     5  vn_rbt1, vn_sgs, vn_lar, vn_modB, vn_ctor, vn_amin, vn_Rmaj,
     5  vn_potsin, vn_potcos, vn_maxpot, vn_xmpot, vn_xnpot,             !diagno/extender output (SPH071414)
     6  vn_vol, vn_mse, vn_thom, vn_ac, vn_ai, vn_am,
     6  vn_pmass_type, vn_pcurr_type, vn_piota_type,
     6  vn_am_aux_s, vn_am_aux_f, vn_ac_aux_s, vn_ac_aux_f,
     6  vn_ai_aux_s, vn_ai_aux_f,
     6  vn_ftolv, vn_fsqr, vn_fsqz, vn_fsql,
     7  vn_pmod, vn_tmod, vn_pmod_nyq, vn_tmod_nyq,
     7  vn_racc, vn_zacs, vn_racs, vn_zacc, vn_iotaf, vn_qfact,
     8  vn_presf, vn_phi, vn_phipf, vn_jcuru, vn_jcurv, vn_iotah,
     8  vn_chi, vn_chipf,
     9  vn_mass, vn_presh, vn_betah, vn_buco, vn_bvco, vn_vp, vn_specw,
     A  vn_phip, vn_jdotb, vn_bdotb, vn_overr, vn_bgrv, vn_merc,
     B  vn_mshear, vn_mwell, vn_mcurr, vn_mgeo, vn_equif, vn_fsq,
     C  vn_wdot, vn_extcur, vn_curlab, vn_rmnc, vn_zmns, vn_lmns,
     D  vn_gmnc, vn_bmnc, vn_bsubumnc, vn_bsubvmnc, vn_bsubsmns,
     E  vn_bsupumnc, vn_bsupvmnc, vn_rmns, vn_zmnc, vn_lmnc, vn_gmns,
     F  vn_bmns, vn_bsubumns, vn_bsubvmns, vn_bsubsmnc, vn_bsupumns,
     G  vn_bsupvmns, vn_rbc, vn_zbs, vn_rbs, vn_zbc,
     H  ln_version, ln_extension, ln_mgrid,

     &  vn_bsubumnc_sur, vn_bsubvmnc_sur,                      !MRC 10-15-15
     &  vn_bsupumnc_sur, vn_bsupvmnc_sur,
     &  vn_bsubumns_sur, vn_bsubvmns_sur,
     &  vn_bsupumns_sur, vn_bsupvmns_sur,
     &  vn_currumnc, vn_currumns, vn_currvmnc, vn_currvmns,    !MRC 8-12-16

     1  ln_magen, ln_therm, ln_gam, ln_maxr, ln_minr, ln_maxz, ln_fp,
     2  ln_radnod, ln_polmod, ln_tormod, ln_maxmod, ln_maxit, ln_actit,
     2  ln_maxpot, ln_potsin, ln_potcos,
     3  ln_asym, ln_recon, ln_free, ln_error, ln_aspect, ln_beta,
     4  ln_pbeta, ln_tbeta, ln_abeta, ln_b0, ln_rbt0, ln_maxmod_nyq,
     5  ln_rbt1, ln_sgs, ln_lar, ln_modB, ln_ctor, ln_amin, ln_Rmaj,
     6  ln_mse, ln_thom, ln_flp, ln_nobd, ln_nbset, ln_next, ln_nbfld,
     7  ln_pmod, ln_tmod, ln_pmod_nyq, ln_tmod_nyq, ln_racc, ln_zacs,
     7  ln_racs, ln_zacc, ln_iotaf, ln_qfact, ln_am, ln_ac, ln_ai,
     7  ln_pmass_type, ln_pcurr_type, ln_piota_type,
     7  ln_am_aux_s, ln_am_aux_f, ln_ac_aux_s, ln_ac_aux_f,
     7  ln_ai_aux_s, ln_ai_aux_f, ln_chi, ln_chipf,
     8  ln_presf, ln_phi, ln_phipf, ln_jcuru, ln_jcurv, ln_iotah,
     9  ln_mass, ln_presh, ln_betah, ln_buco, ln_bvco, ln_vp, ln_specw,
     A  ln_vol, ln_phip, ln_jdotb, ln_bdotb, ln_bgrv, ln_merc,
     B  ln_mshear, ln_mwell, ln_mcurr, ln_mgeo, ln_equif, ln_fsq,
     C  ln_wdot, ln_extcur, ln_curlab, ln_rmnc, ln_zmns, ln_lmns,
     D  ln_gmnc, ln_bmnc, ln_bsubumnc, ln_bsubvmnc, ln_bsubsmns,
     E  ln_bsupumnc, ln_bsupvmnc, ln_rmns, ln_zmnc, ln_lmnc, ln_gmns,
     F  ln_bmns, ln_bsubumns, ln_bsubvmns, ln_bsubsmnc, ln_bsupumns,
     G  ln_bsupvmns, ln_rbc, ln_zbs, ln_rbs, ln_zbc,

     &  ln_bsubumnc_sur, ln_bsubvmnc_sur,          !MRC 10-15-15
     &  ln_bsupumnc_sur, ln_bsupvmnc_sur,
     &  ln_bsubumns_sur, ln_bsubvmns_sur,
     &  ln_bsupumns_sur, ln_bsupvmns_sur,
     &  ln_currumnc, ln_currumns, ln_currvmnc, ln_currvmns

!------------------DEC$ ELSE !to use safe_open_mod in any case (J.Geiger)

      USE safe_open_mod
      USE mgrid_mod

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: ier_flag
      REAL(dp), DIMENSION(mnmax,ns,3*MAX(ntmax/2,1)),           !reverse ns, mnmax for backwards compatibility
     1   INTENT(inout), TARGET :: rzl_array, gc_array
      REAL(dp), DIMENSION(ns,nznt), INTENT(inout) ::
     1   bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu
      REAL(dp) :: qfact(ns)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(dp), PARAMETER :: c1p5 = 1.5_dp
      CHARACTER(LEN=*), PARAMETER, DIMENSION(1) ::
     1             r1dim = (/'radius'/), mn1dim = (/'mn_mode'/),
     2             mn2dim = (/'mn_mode_nyq'/),
     2             mnpotdim = (/'mn_mode_pot'/),
     3             currg = (/'ext_current'/),
     4             currl = (/'current_label'/)
      CHARACTER(LEN=*), DIMENSION(2), PARAMETER ::
     1             r2dim = (/'mn_mode','radius '/),
     1             r3dim = (/'mn_mode_nyq','radius     '/)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: j, js, jlk, mn, lk, iasym,
     1           m, n, k, iwout0, n1, nwout, istat, i, indx1(1),
     2           mnmax_nyq0, mnyq0, nnyq0, nwout2   ! nwout2 by J.Geiger
     3          ,isgn, js2, nfort      !for diagno 1.5
      REAL(dp) :: dmult, tcosi, tsini, vversion, sgn, tmult,
     1            presfactor, ftolx1, d_bsupumn, d_bsupvmn   ! diagno 1.5
      REAL(dp), POINTER, DIMENSION(:,:) :: rmnc, rmns, zmns,
     1   zmnc, lmns, lmnc
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) ::
     1   gmnc, bmnc, gmns, bmns,
     2   bsubumnc, bsubvmnc, bsubsmns, bsubumns, bsubvmns, bsubsmnc,
     3   currumnc, currvmnc, currumns, currvmns
      REAL(dp), DIMENSION(mnmax) :: rmnc1, zmns1, lmns1,
     1   rmns1, zmnc1, lmnc1, bmodmn, bmodmn1
      REAL(dp), DIMENSION(:), ALLOCATABLE :: gmn, bmn,
     1   bsubumn, bsubvmn, bsubsmn, bsupumn, bsupvmn
      REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubumnc_sur  !MRC 10-15-15
      REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubvmnc_sur
      REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupumnc_sur
      REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupvmnc_sur
      REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubumns_sur
      REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubvmns_sur
      REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupumns_sur
      REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupvmns_sur
      REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubua_sur, bsubva_sur
      REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupua_sur, bsupva_sur

      CHARACTER(LEN=120) :: wout_file, wout2_file         ! wout2_file by J.Geiger
      CHARACTER(LEN=120) :: fort_file   ! fort_file for diagno 1.5
      REAL(dp), DIMENSION(:), ALLOCATABLE :: xfinal
      REAL(dp), DIMENSION(:), POINTER ::   xm_nyq0, xn_nyq0
!     ELIMINATE THESE EVENTUALLY
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) ::
     1   bsupumnc, bsupumns, bsupvmnc, bsupvmns

      LOGICAL :: lcurr
      INTEGER :: nmin0     ! J Geiger:   Added for diagno-file

!-----------------------------------------------
      CALL second0 (twouton)
!
!  Pointer assignments for storage arrays
!
      n1 = MAX(1,ntmax/2)
      rmnc => rzl_array(:,:,1)            !!store COS(mu-nv) components
      zmns => rzl_array(:,:,1+n1)         !!store SIN(mu-nv)
      lmns => rzl_array(:,:,1+2*n1)       !!store SIN(mu-nv)

      IF (lasym) THEN
         rmns => gc_array(:,:,1)            !!store SIN(mu-nv)
         zmnc => gc_array(:,:,1+n1)         !!store COS(mu-nv)
         lmnc => gc_array(:,:,1+2*n1)       !!store COS(mu-nv)
      END IF

!
!     THIS SUBROUTINE CREATES THE FILE WOUT.IT CONTAINS THE CYLINDRICAL COORDINATE SPECTRAL
!     COEFFICIENTS RMN,ZMN (full), LMN (half_mesh - CONVERTED FROM
!     INTERNAL full REPRESENTATION), AS WELL AS COEFFICIENTS (ON NYQ MESH) FOR COMPUTED
!     QUANTITIES:
!
!     BSQ, BSUPU,V, BSUBU,V, GSQRT (HALF); BSUBS (FULL-CONVERTED IN JXBFORCE)
!
      mnmax_nyq0 = mnmax_nyq
      mnyq0 = mnyq
      nnyq0 = nnyq
      xm_nyq0 => xm_nyq; xn_nyq0 => xn_nyq

      ALLOCATE (gmn(mnmax_nyq0), bmn(mnmax_nyq0),
     1   bsubumn(mnmax_nyq0), bsubvmn(mnmax_nyq0), bsubsmn(mnmax_nyq0),
     2   bsupumn(mnmax_nyq0), bsupvmn(mnmax_nyq0),
     6   stat=istat)

      IF (lfreeb) THEN        !MRC 10-15-15
         ALLOCATE (bsubua_sur(nzeta*ntheta2), bsubva_sur(nzeta*ntheta2))
         ALLOCATE (bsupua_sur(nzeta*ntheta2), bsupva_sur(nzeta*ntheta2))

         ALLOCATE (bsubumnc_sur(mnmax_nyq0), bsubvmnc_sur(mnmax_nyq0))
         ALLOCATE (bsupumnc_sur(mnmax_nyq0), bsupvmnc_sur(mnmax_nyq0))
         IF (lasym) THEN
            ALLOCATE (bsubumns_sur(mnmax_nyq0),                                &
     &                bsubvmns_sur(mnmax_nyq0))
            ALLOCATE (bsupumns_sur(mnmax_nyq0),                                &
     &                bsupvmns_sur(mnmax_nyq0))
         END IF
      END IF

      ALLOCATE (gmnc(mnmax_nyq0,ns), bmnc(mnmax_nyq0,ns),
     1          bsubumnc(mnmax_nyq0,ns), bsubvmnc(mnmax_nyq0,ns),
     2          bsubsmns(mnmax_nyq0,ns), bsupumnc(mnmax_nyq0,ns),
     3          bsupvmnc(mnmax_nyq0,ns),
     4          currumnc(mnmax_nyq0,ns), currvmnc(mnmax_nyq0,ns),
     9          stat=istat)
      IF (lasym) THEN
      ALLOCATE (gmns(mnmax_nyq0,ns), bmns(mnmax_nyq0,ns),
     1          bsubumns(mnmax_nyq0,ns), bsubvmns(mnmax_nyq0,ns),
     2          bsubsmnc(mnmax_nyq0,ns), bsupumns(mnmax_nyq0,ns),
     3          bsupvmns(mnmax_nyq0,ns),
     4          currumns(mnmax_nyq0,ns), currvmns(mnmax_nyq0,ns),
     8          stat=istat)
      END IF
      IF (istat .ne. 0) STOP 'Error allocating arrays in VMEC WROUT'

!      IF (nextcur .eq. 0) THEN
!         DO j = SIZE(extcur), 1, -1
!           IF (extcur(j) .ne. zero) THEN
!               nextcur = j
!               EXIT
!            END IF
!         END DO
!      END IF

! ftol info evaluated here!
      indx1=MAXLOC(ns_array)
      ftolx1=ftol_array(indx1(1))

!     NYQUIST FREQUENCY REQUIRES FACTOR OF 1/2
      IF (mnyq .ne. 0) cosmui(:,mnyq) = p5*cosmui(:,mnyq)
      IF (nnyq .ne. 0) cosnv (:,nnyq) = p5*cosnv (:,nnyq)

      wout_file = version_
      READ (wout_file, *) vversion

      wout_file = 'wout_' // TRIM(input_extension) // '.nc'
      CALL cdf_open(nwout,wout_file,'w',iwout0)
      IF (iwout0 .ne. 0) STOP 'Error opening wout.nc file VMEC WROUT'

!================================
! Define Variables
!================================
!  Scalars
      CALL cdf_define(nwout, vn_version, vversion)
      CALL cdf_define(nwout, vn_extension, input_extension)
      CALL cdf_define(nwout, vn_mgrid, mgrid_file)
      CALL cdf_define(nwout, vn_pcurr_type, pcurr_type)
      CALL cdf_define(nwout, vn_pmass_type, pmass_type)
      CALL cdf_define(nwout, vn_piota_type, piota_type)
      CALL cdf_define(nwout, vn_magen, wb)
      CALL cdf_define(nwout, vn_therm, wp)
      CALL cdf_define(nwout, vn_gam, adiabatic)
      CALL cdf_define(nwout, vn_maxr, rmax_surf)
      CALL cdf_define(nwout, vn_minr, rmin_surf)
      CALL cdf_define(nwout, vn_maxz, zmax_surf)
      CALL cdf_define(nwout, vn_fp, nfp)
      CALL cdf_define(nwout, vn_radnod, ns)
      CALL cdf_define(nwout, vn_polmod, mpol)
      CALL cdf_define(nwout, vn_tormod, ntor)
      CALL cdf_define(nwout, vn_maxmod, mnmax)
      CALL cdf_define(nwout, vn_maxmod_nyq, mnmax_nyq0)
      CALL cdf_define(nwout, vn_maxit, iter2)
      CALL cdf_define(nwout, vn_actit, itfsq)
      CALL cdf_define(nwout, vn_asym, lasym)
      CALL cdf_define(nwout, vn_free, lfreeb)
      CALL cdf_define(nwout, vn_error, ier_flag)
      CALL cdf_define(nwout, vn_aspect, aspect)
      CALL cdf_define(nwout, vn_beta, betatot)
      CALL cdf_define(nwout, vn_pbeta, betapol)
      CALL cdf_define(nwout, vn_tbeta, betator)
      CALL cdf_define(nwout, vn_abeta, betaxis)
      CALL cdf_define(nwout, vn_b0, b0)
      CALL cdf_define(nwout, vn_rbt0, rbtor0)
      CALL cdf_define(nwout, vn_rbt1, rbtor)
      CALL cdf_define(nwout, vn_sgs, NINT(signgs))
      CALL cdf_define(nwout, vn_lar, IonLarmor)
      CALL cdf_define(nwout, vn_modB, volAvgB)
      CALL cdf_define(nwout, vn_ctor, ctor)
      CALL cdf_define(nwout, vn_amin, Aminor_p)
      CALL cdf_define(nwout, vn_Rmaj, Rmajor_p)
      CALL cdf_define(nwout, vn_vol, volume_p)
      CALL cdf_define(nwout, vn_ftolv, ftolx1)
      CALL cdf_define(nwout, vn_fsql, fsql)
      CALL cdf_define(nwout, vn_fsqr, fsqr)
      CALL cdf_define(nwout, vn_fsqz, fsqz)

      CALL cdf_define(nwout, vn_nextcur, nextcur)
      CALL cdf_define(nwout, vn_extcur, extcur(1:nextcur),
     1                dimname=currg)
      CALL cdf_define(nwout, vn_mgmode, mgrid_mode)
      IF (lfreeb) THEN
         CALL cdf_define(nwout, vn_maxpot, mnpd)
      END IF

! 1D Arrays

      CALL cdf_define(nwout, vn_pmod, xm, dimname=mn1dim)
      CALL cdf_setatt(nwout, vn_pmod, ln_pmod)
      CALL cdf_define(nwout, vn_tmod, xn, dimname=mn1dim)
      CALL cdf_setatt(nwout, vn_tmod, ln_tmod)
      CALL cdf_define(nwout, vn_pmod_nyq, xm_nyq0, dimname=mn2dim)
      CALL cdf_setatt(nwout, vn_pmod_nyq, ln_pmod_nyq)
      CALL cdf_define(nwout, vn_tmod_nyq, xn_nyq0, dimname=mn2dim)
      CALL cdf_setatt(nwout, vn_tmod_nyq, ln_tmod_nyq)

      CALL cdf_define(nwout, vn_racc, raxis_cc(0:ntor),
     1                dimname=(/'n_tor'/))
      CALL cdf_setatt(nwout, vn_racc, ln_racc)
      CALL cdf_define(nwout, vn_zacs, zaxis_cs(0:ntor),
     1                dimname=(/'n_tor'/))
      CALL cdf_setatt(nwout, vn_zacs, ln_zacs)
      IF (lasym) THEN
         CALL cdf_define(nwout, vn_racs, raxis_cs(0:ntor),
     1                dimname=(/'n_tor'/))
         CALL cdf_setatt(nwout, vn_racs, ln_racs)
         CALL cdf_define(nwout, vn_zacc, zaxis_cc(0:ntor),
     1                dimname=(/'n_tor'/))
         CALL cdf_setatt(nwout, vn_zacc, ln_zacc)
      END IF

      j = SIZE(am)-1
      CALL cdf_define(nwout, vn_am, am(0:j),
     1                dimname=(/'preset'/))
      j = SIZE(ac)-1
      CALL cdf_define(nwout, vn_ac, ac(0:j),
     1                dimname=(/'preset'/))
      j = SIZE(ai)-1
      CALL cdf_define(nwout, vn_ai, ai(0:j),
     1                dimname=(/'preset'/))

      j = SIZE(am_aux_s)
      CALL cdf_define(nwout, vn_am_aux_s, am_aux_s(1:j),
     1                dimname=(/'ndfmax'/))
      j = SIZE(am_aux_f)
      CALL cdf_define(nwout, vn_am_aux_f, am_aux_f(1:j),
     1                dimname=(/'ndfmax'/))
      j = SIZE(ai_aux_s)
      CALL cdf_define(nwout, vn_ai_aux_s, ai_aux_s(1:j),
     1                dimname=(/'ndfmax'/))
      j = SIZE(ai_aux_f)
      CALL cdf_define(nwout, vn_ai_aux_f, ai_aux_f(1:j),
     1                dimname=(/'ndfmax'/))
      j = SIZE(ac_aux_s)
      CALL cdf_define(nwout, vn_ac_aux_s, ac_aux_s(1:j),
     1                dimname=(/'ndfmax'/))
      j = SIZE(ac_aux_f)
      CALL cdf_define(nwout, vn_ac_aux_f, ac_aux_f(1:j),
     1                dimname=(/'ndfmax'/))

      CALL cdf_define(nwout, vn_iotaf, iotaf(1:ns),
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_iotaf, ln_iotaf)

      qfact=HUGE(qfact)
      WHERE (iotaf(1:ns) .NE. zero) qfact=one/iotaf(1:ns)

      CALL cdf_define(nwout, vn_qfact, qfact(1:ns),
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_iotaf, ln_qfact)
      CALL cdf_define(nwout, vn_presf, presf,
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_presf, ln_presf, units='Pa')
      CALL cdf_define(nwout, vn_phi, phi,
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_phi, ln_phi, units='wb')
      CALL cdf_define(nwout, vn_phipf,
     1                phipf, dimname=r1dim)
      CALL cdf_setatt(nwout, vn_phipf, ln_phipf)
      CALL cdf_define(nwout, vn_chi, chi,
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_chi, ln_chi, units='wb')
      CALL cdf_define(nwout, vn_chipf,
     1                phipf, dimname=r1dim)
      CALL cdf_setatt(nwout, vn_chipf, ln_chipf)
      CALL cdf_define(nwout, vn_jcuru,
     1                jcuru, dimname=r1dim)
      CALL cdf_define(nwout, vn_jcurv,
     1                jcurv, dimname=r1dim)

      CALL cdf_define(nwout, vn_iotah, iotas(1:ns),
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_iotah, ln_iotah)
      CALL cdf_define(nwout, vn_mass, mass,
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_mass, ln_mass)
      CALL cdf_define(nwout, vn_presh, pres(1:ns),
     1                dimname=r1dim)
      CALL cdf_setatt(nwout, vn_presh, ln_presh, units='Pa')
      CALL cdf_define(nwout, vn_betah, beta_vol,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_buco, buco,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_bvco, bvco,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_vp, vp(1:ns),
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_specw, specw,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_phip,
     1                phips(1:ns), dimname=r1dim)
      CALL cdf_define(nwout, vn_overr,
     2                overr(1:ns), dimname=r1dim)

      CALL cdf_define(nwout, vn_jdotb, jdotb,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_bdotb, bdotb,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_bgrv, bdotgradv,
     1                dimname=r1dim)

      CALL cdf_define(nwout, vn_merc, Dmerc,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_mshear, Dshear,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_mwell, Dwell,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_mcurr, Dcurr,
     1                dimname=r1dim)
      CALL cdf_define(nwout, vn_mgeo,
     1                Dgeod, dimname=r1dim)
      CALL cdf_define(nwout, vn_equif,
     1                equif, dimname=r1dim)

      CALL cdf_define(nwout, vn_fsq, fsqt(1:nstore_seq),
     1                dimname=(/'time'/))
      CALL cdf_define(nwout, vn_wdot, wdot(1:nstore_seq),
     1                dimname=(/'time'/))

      IF (lfreeb) THEN
         CALL cdf_define(nwout, vn_potsin, potvac(1:mnpd),
     1                   dimname=mnpotdim)
         CALL cdf_setatt(nwout, vn_potsin, ln_potsin)
         CALL cdf_define(nwout, vn_xmpot, xmpot(1:mnpd),
     1                   dimname=mnpotdim)
         CALL cdf_define(nwout, vn_xnpot, xnpot(1:mnpd),
     1                   dimname=mnpotdim)
         IF (lasym) THEN
            CALL cdf_define(nwout, vn_potcos,
     1                      potvac(1+mnpd:2*mnpd), dimname=mnpotdim)
            CALL cdf_setatt(nwout, vn_potcos, ln_potcos)
         END IF

         IF (nextcur.gt.0 .and. ALLOCATED(curlabel)) THEN
         CALL cdf_define(nwout, vn_curlab,
     1        curlabel(1:nextcur), dimname=currl)
         END IF
      ENDIF

! 2D Arrays
      CALL cdf_define(nwout, vn_rmnc, rmnc, dimname=r2dim)
      CALL cdf_setatt(nwout, vn_rmnc, ln_rmnc, units='m')
      CALL cdf_define(nwout, vn_zmns, zmns, dimname=r2dim)
      CALL cdf_setatt(nwout, vn_zmns, ln_zmns, units='m')
      CALL cdf_define(nwout, vn_lmns, lmns, dimname=r2dim)
      CALL cdf_setatt(nwout, vn_lmns, ln_lmns)
      CALL cdf_define(nwout, vn_gmnc, gmnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_gmnc, ln_gmnc)
      CALL cdf_define(nwout, vn_bmnc, bmnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bmnc, ln_bmnc)
      CALL cdf_define(nwout, vn_bsubumnc, bsubumnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bsubumnc, ln_bsubumnc)
      CALL cdf_define(nwout, vn_bsubvmnc, bsubvmnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bsubvmnc, ln_bsubvmnc)
      CALL cdf_define(nwout, vn_bsubsmns, bsubsmns, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_bsubsmns, ln_bsubsmns)

      CALL cdf_define(nwout, vn_currumnc, currumnc, dimname=r3dim)    !MRC 8-12-16
      CALL cdf_setatt(nwout, vn_currumnc, ln_currumnc)
      CALL cdf_define(nwout, vn_currvmnc, currvmnc, dimname=r3dim)
      CALL cdf_setatt(nwout, vn_currvmnc, ln_currvmnc)

      IF (lfreeb) THEN
         CALL cdf_define(nwout, vn_bsubumnc_sur, bsubumnc_sur,                 &
     &                   dimname=mn2dim)
         CALL cdf_setatt(nwout, vn_bsubumnc_sur, ln_bsubumnc_sur)
         CALL cdf_define(nwout, vn_bsubvmnc_sur, bsubvmnc_sur,                 &
     &                   dimname=mn2dim)
         CALL cdf_setatt(nwout, vn_bsubvmnc_sur, ln_bsubvmnc_sur)
         CALL cdf_define(nwout, vn_bsupumnc_sur, bsupumnc_sur,                 &
     &                   dimname=mn2dim)
         CALL cdf_setatt(nwout, vn_bsupumnc_sur, ln_bsupumnc_sur)
         CALL cdf_define(nwout, vn_bsupvmnc_sur, bsupvmnc_sur,                 &
     &                   dimname=mn2dim)
         CALL cdf_setatt(nwout, vn_bsupvmnc_sur, ln_bsupvmnc_sur)
      END IF

!     ELIMINATE THESE EVENTUALLY: DON'T NEED THEM - CAN COMPUTE FROM GSQRT
      CALL cdf_define(nwout, vn_bsupumnc, bsupumnc, dimname=r3dim)
      CALL cdf_define(nwout, vn_bsupvmnc, bsupvmnc, dimname=r3dim)
!     IF (lfreeb) THEN
!         CALL cdf_define(nwout, vn_rbc, rbc,
!    1                dimname=(/'n_mode','m_mode'/))
!         CALL cdf_setatt(nwout, vn_rbc, ln_rbc, units='m')
!         CALL cdf_define(nwout, vn_zbs, zbs,
!    1                dimname=(/'n_mode','m_mode'/))
!         CALL cdf_setatt(nwout, vn_zbs, ln_zbs, units='m')
!        IF (lasym) THEN
!           CALL cdf_define(nwout, vn_rbs, rbs,
!    1                dimname=(/'n_mode','m_mode'/))
!           CALL cdf_define(nwout, vn_zbc, zbc,
!    1                dimname=(/'n_mode','m_mode'/))
!        END IF
!     END IF

      IF (lasym) then
         CALL cdf_define(nwout, vn_rmns, rmns, dimname=r2dim)
         CALL cdf_setatt(nwout, vn_rmns, ln_rmns, units='m')
         CALL cdf_define(nwout, vn_zmnc, zmnc, dimname=r2dim)
         CALL cdf_setatt(nwout, vn_zmnc, ln_zmnc, units='m')
         CALL cdf_define(nwout, vn_lmnc, lmnc, dimname=r2dim)
         CALL cdf_setatt(nwout, vn_lmnc, ln_lmnc)
         CALL cdf_define(nwout, vn_gmns, gmns, dimname=r3dim)
         CALL cdf_setatt(nwout, vn_gmns, ln_gmns)
         CALL cdf_define(nwout, vn_bmns, bmns, dimname=r3dim)
         CALL cdf_setatt(nwout, vn_bmns, ln_bmns)
         CALL cdf_define(nwout, vn_bsubumns, bsubumns, dimname=r3dim)
         CALL cdf_setatt(nwout, vn_bsubumns, ln_bsubumns)
         CALL cdf_define(nwout, vn_bsubvmns, bsubvmns, dimname=r3dim)
         CALL cdf_setatt(nwout, vn_bsubvmns, ln_bsubvmns)
         CALL cdf_define(nwout, vn_bsubsmnc, bsubsmnc, dimname=r3dim)
         CALL cdf_setatt(nwout, vn_bsubsmnc, ln_bsubsmnc)

         CALL cdf_define(nwout, vn_currumns, currumns, dimname=r3dim)
         CALL cdf_setatt(nwout, vn_currumns, ln_currumns)
         CALL cdf_define(nwout, vn_currvmns, currvmns, dimname=r3dim)
         CALL cdf_setatt(nwout, vn_currvmns, ln_currvmns)

         IF (lfreeb) THEN
            CALL cdf_define(nwout, vn_bsubumns_sur, bsubumns_sur,                 &
     &                      dimname=mn2dim)
            CALL cdf_setatt(nwout, vn_bsubumns_sur, ln_bsubumns_sur)
            CALL cdf_define(nwout, vn_bsubvmns_sur, bsubvmns_sur,                 &
     &                      dimname=mn2dim)
            CALL cdf_setatt(nwout, vn_bsubvmns_sur, ln_bsubvmns_sur)
            CALL cdf_define(nwout, vn_bsupumns_sur, bsupumns_sur,                 &
     &                      dimname=mn2dim)
            CALL cdf_setatt(nwout, vn_bsupumns_sur, ln_bsupumns_sur)
            CALL cdf_define(nwout, vn_bsupvmns_sur, bsupvmns_sur,                 &
     &                      dimname=mn2dim)
            CALL cdf_setatt(nwout, vn_bsupvmns_sur, ln_bsupvmns_sur)
         END IF

!        ELIMINATE THESE EVENTUALLY: DON'T NEED THEM
         CALL cdf_define(nwout, vn_bsupumns, bsupumns, dimname=r3dim)
         CALL cdf_define(nwout, vn_bsupvmns, bsupvmns, dimname=r3dim)
      end if ! lasym

!================================
! Write Variables
!================================

! Scalars
      CALL cdf_write(nwout, vn_version, vversion)
      CALL cdf_write(nwout, vn_extension, input_extension)
      CALL cdf_write(nwout, vn_mgrid, mgrid_file)
      CALL cdf_write(nwout, vn_pcurr_type, pcurr_type)
      CALL cdf_write(nwout, vn_piota_type, piota_type)
      CALL cdf_write(nwout, vn_pmass_type, pmass_type)
      CALL cdf_write(nwout, vn_magen, wb)
      CALL cdf_write(nwout, vn_therm, wp)
      CALL cdf_write(nwout, vn_gam, adiabatic)
      CALL cdf_write(nwout, vn_maxr, rmax_surf)
      CALL cdf_write(nwout, vn_minr, rmin_surf)
      CALL cdf_write(nwout, vn_maxz, zmax_surf)
      CALL cdf_write(nwout, vn_fp, nfp)
      CALL cdf_write(nwout, vn_radnod, ns)
      CALL cdf_write(nwout, vn_polmod, mpol)
      CALL cdf_write(nwout, vn_tormod, ntor)
      CALL cdf_write(nwout, vn_maxmod, mnmax)
      CALL cdf_write(nwout, vn_maxmod_nyq, mnmax_nyq0)
      CALL cdf_write(nwout, vn_maxit, iter2)
      CALL cdf_write(nwout, vn_actit, itfsq)
      CALL cdf_write(nwout, vn_asym, lasym)
      CALL cdf_write(nwout, vn_free, lfreeb)
      CALL cdf_write(nwout, vn_error, ier_flag)
!
      CALL cdf_write(nwout, vn_aspect, aspect)
      CALL cdf_write(nwout, vn_beta, betatot)
      CALL cdf_write(nwout, vn_pbeta, betapol)
      CALL cdf_write(nwout, vn_tbeta, betator)
      CALL cdf_write(nwout, vn_abeta, betaxis)
      CALL cdf_write(nwout, vn_b0, b0)
      CALL cdf_write(nwout, vn_rbt0, rbtor0)
      CALL cdf_write(nwout, vn_rbt1, rbtor)
      CALL cdf_write(nwout, vn_sgs, NINT(signgs))
      CALL cdf_write(nwout, vn_lar, IonLarmor)
      CALL cdf_write(nwout, vn_modB, volAvgB)
      CALL cdf_write(nwout, vn_ctor, ctor/mu0)
      CALL cdf_write(nwout, vn_amin, Aminor_p)
      CALL cdf_write(nwout, vn_rmaj, Rmajor_p)
      CALL cdf_write(nwout, vn_vol, volume_p)
      CALL cdf_write(nwout, vn_ftolv, ftolx1)
      CALL cdf_write(nwout, vn_fsql, fsql)
      CALL cdf_write(nwout, vn_fsqr, fsqr)
      CALL cdf_write(nwout, vn_fsqz, fsqz)

      CALL cdf_write(nwout, vn_nextcur, nextcur)
      IF (nextcur .gt. 0) THEN
         CALL cdf_write(nwout, vn_extcur, extcur(1:nextcur))
         CALL cdf_write(nwout, vn_mgmode, mgrid_mode)
      ENDIF
      IF (lfreeb) THEN
         CALL cdf_write(nwout, vn_maxpot, mnpd)
         IF (nextcur.gt.0 .and. ALLOCATED(curlabel))
     1   CALL cdf_write(nwout, vn_curlab, curlabel(1:nextcur))
      END IF

! 1D Arrays
      CALL cdf_write(nwout, vn_pmod, xm)
      CALL cdf_write(nwout, vn_tmod, xn)
      CALL cdf_write(nwout, vn_pmod_nyq, xm_nyq0)
      CALL cdf_write(nwout, vn_tmod_nyq, xn_nyq0)

      IF (lfreeb) THEN
         CALL cdf_write(nwout, vn_potsin, potvac(1:mnpd))
         IF (lasym)                                                            &
     &      CALL cdf_write(nwout, vn_potcos, potvac(1+mnpd:2*mnpd))
         CALL cdf_write(nwout, vn_xmpot, xmpot)
         CALL cdf_write(nwout, vn_xnpot, xnpot)
      END IF

!---------------------DEC$ ENDIF

      ALLOCATE (xfinal(neqs), stat=js)
      IF (js .NE. 0) STOP 'Allocation error for xfinal in WROUT!'
      xfinal = xc
!
!     MUST CONVERT m=1 MODES... FROM INTERNAL TO PHYSICAL FORM
!     Extrapolation of m=0 Lambda (cs) modes, which are not evolved at j=1, done in CONVERT
!
      lk = ns*ntor1
      IF (lthreed) CALL convert_sym_par  (xfinal(1+mns*(rss-1)+lk),
     1                                xfinal(1+irzloff+mns*(zcs-1)+lk),
     2                                1, ns)
      IF (lasym)   CALL convert_asym_par (xfinal(1+mns*(rsc-1)+lk),
     1                                xfinal(1+irzloff+mns*(zcc-1)+lk),
     2                                1, ns)
!
!     CONVERT TO rmnc, zmns, lmns, etc EXTERNAL representation (without internal mscale, nscale)
!     IF B^v ~ phip + lamu, MUST DIVIDE BY phipf(js) below to maintain old-style format
!     THIS COULD BE A PROBLEM FOR RFP WHERE PHIPF->0 INSIDE THE PLASMA!
!
      RADIUS1: DO js = 1, ns

         CALL convert (rmnc1, zmns1, lmns1, rmns1, zmnc1, lmnc1,
     1                         xfinal, js)

         rmnc(:,js) = rmnc1(:)
         zmns(:,js) = zmns1(:)
         lmns(:,js) = (lmns1(:)/phipf(js)) * lamscale
         IF (lasym) THEN
            rmns(:,js) = rmns1(:)
            zmnc(:,js) = zmnc1(:)
            lmnc(:,js) = (lmnc1(:)/phipf(js)) * lamscale
         END IF

      END DO RADIUS1

      DEALLOCATE (xfinal)

!
!     INTERPOLATE LAMBDA ONTO HALF-MESH FOR BACKWARDS CONSISTENCY WITH EARLIER VERSIONS OF VMEC
!     AND SMOOTHS POSSIBLE UNPHYSICAL "WIGGLE" ON RADIAL MESH
!

      WHERE (NINT(xm) .le. 1) lmns(:,1) = lmns(:,2)
      DO js = ns,2,-1
         WHERE (MOD(NINT(xm),2) .eq. 0)
            lmns(:,js) = p5*(lmns(:,js) + lmns(:,js-1))
         ELSEWHERE
            lmns(:,js) = p5*(sm(js)*lmns(:,js) + sp(js-1)*lmns(:,js-1))
         END WHERE
      END DO

      lmns(:,1) = 0
      raxis_cc(0:ntor) = rmnc(1:ntor+1,1)
      zaxis_cs(0:ntor) = zmns(1:ntor+1,1)

      IF (lasym) then
         WHERE (NINT(xm) .le. 1) lmnc(:,1) = lmnc(:,2)
         DO js = ns,2,-1
            WHERE (MOD(NINT(xm),2) .eq. 0)
               lmnc(:,js) = p5*(lmnc(:,js) + lmnc(:,js-1))
            ELSEWHERE
               lmnc(:,js) = p5*(sm(js)*lmnc(:,js)
     &                      + sp(js-1)*lmnc(:,js-1))
            END WHERE
         END DO

         lmnc(:,1) = 0;
         raxis_cs(0:ntor) = rmns(1:ntor+1,1)
         zaxis_cc(0:ntor) = zmnc(1:ntor+1,1)
      end if ! lasym


!SPH100209: COMPUTE |B| = SQRT(|B|**2) and store in bsq, bsqa
      DO js = 2, ns
         bsq(js,:nznt) = SQRT(2*ABS(bsq(js,:nznt)-pres(js)))
      END DO

      tmult = p5/r0scale**2
!SPH: FIXED THIS 03-05-07 TO CALL symmetrization routine
      IF (lasym) THEN
!Changed integration norm in fixaray, SPH012314
         tmult = 2*tmult
         bsubs(1,:) = 0
         CALL symoutput (bsq,   gsqrt,  bsubu,  bsubv,  bsupu,
     1                   bsupv,  bsubs,
     4                   bsqa,  gsqrta, bsubua, bsubva, bsupua,
     5                   bsupva, bsubsa)

         IF (lfreeb) THEN     !MRC  10-15-15
            CALL symoutput_sur(bsubu_sur, bsubv_sur,                           &
     &                         bsupu_sur, bsupv_sur,                           &
     &                         bsubua_sur, bsubva_sur,                         &
     &                         bsupua_sur, bsupva_sur)
         END IF
      END IF

!         DO js = 2, ns
!            WRITE (200, *) 'JS: ', js, 'BSUBU, BSUBV'
!            WRITE (200, '(1p,6e12.4)') bsubu(js,:), bsubv(js,:)
!         END DO

      RADIUS2: DO js = 2, ns
         gmn = 0
         bmn = 0
         bsubumn = 0
         bsubvmn = 0
         bsubsmn = 0
         bsupumn = 0
         bsupvmn = 0

         MN2: DO mn = 1, mnmax_nyq0
            n = NINT(xn_nyq0(mn))/nfp
            m = NINT(xm_nyq0(mn))
            n1 = ABS(n)
            dmult = mscale(m)*nscale(n1)*tmult
            IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
            sgn = SIGN(1, n)
            lk = 0
            DO j = 1, ntheta2
               DO k = 1, nzeta
                  lk = lk + 1
                  tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     1                       sgn*sinmui(j,m)*sinnv(k,n1))          !cos(mu - nv)
                  tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     1                       sgn*cosmui(j,m)*sinnv(k,n1))          !sin(mu - nv)
                  bmn(mn) = bmn(mn) + tcosi*bsq(js,lk)
                  gmn(mn) = gmn(mn) + tcosi*gsqrt(js,lk)
                  bsubumn(mn) = bsubumn(mn) + tcosi*bsubu(js,lk)
                  bsubvmn(mn) = bsubvmn(mn) + tcosi*bsubv(js,lk)
                  bsubsmn(mn) = bsubsmn(mn) + tsini*bsubs(js,lk)
                  bsupumn(mn) = bsupumn(mn) + tcosi*bsupu(js,lk)
                  bsupvmn(mn) = bsupvmn(mn) + tcosi*bsupv(js,lk)
               END DO
            END DO
         END DO MN2

         IF (js .eq. ns/2) bmodmn = bmn(1:mnmax)
         IF (js .eq. ns) bmodmn1 = bmn(1:mnmax)
         gmnc(:,js) = gmn(:)
         bmnc(:,js) = bmn(:)
         bsubumnc(:,js) = bsubumn(:)
         bsubvmnc(:,js) = bsubvmn(:)
         bsubsmns(:,js) = bsubsmn(:)
         bsupumnc(:,js) = bsupumn(:)
         bsupvmnc(:,js) = bsupvmn(:)
      END DO RADIUS2

      IF (lfreeb) THEN    !MRC    10-15-15
         bsubumnc_sur = 0
         bsubvmnc_sur = 0
         bsupumnc_sur = 0
         bsupvmnc_sur = 0
         DO mn = 1, mnmax_nyq0
            n = NINT(xn_nyq0(mn))/nfp
            m = NINT(xm_nyq0(mn))
            n1 = ABS(n)
            dmult = mscale(m)*nscale(n1)*tmult
            IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
            sgn = SIGN(1, n)
            lk = 0
            DO j = 1, ntheta2
               DO k = 1, nzeta
                  lk = lk + 1
                  tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     1                       sgn*sinmui(j,m)*sinnv(k,n1))
                  bsubumnc_sur(mn) = bsubumnc_sur(mn)                          &
     &                             + tcosi*bsubu_sur(lk)
                  bsubvmnc_sur(mn) = bsubvmnc_sur(mn)                          &
     &                             + tcosi*bsubv_sur(lk)
                  bsupumnc_sur(mn) = bsupumnc_sur(mn)                          &
     &                             + tcosi*bsupu_sur(lk)
                  bsupvmnc_sur(mn) = bsupvmnc_sur(mn)                          &
     &                             + tcosi*bsupv_sur(lk)
               END DO
            END DO
         END DO
      END IF

      gmnc(:,1) = 0; bmnc(:,1) = 0;
      bsubumnc(:,1) = 0
      bsubvmnc(:,1) = 0
      bsubsmns(:,1) = 2*bsubsmns(:,2) - bsubsmns(:,3)
      bsupumnc(:,1) = 0;  bsupvmnc(:,1) = 0

      IF (lasym) then
         RADIUS3: DO js = 2, ns
            gmn = 0
            bmn = 0
            bsubumn = 0
            bsubvmn = 0
            bsubsmn = 0
            bsupumn = 0
            bsupvmn = 0

            MN3: DO mn = 1, mnmax_nyq0
               n = NINT(xn_nyq0(mn))/nfp
               m = NINT(xm_nyq0(mn))
               n1 = ABS(n)
               dmult = mscale(m)*nscale(n1)*tmult
               IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
               sgn = SIGN(1, n)
               lk = 0
               jlk = js
               DO j = 1, ntheta2
                  DO k = 1, nzeta
                     lk = lk + 1
                     tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     1                          sgn*sinmui(j,m)*sinnv(k,n1))
                     tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     1                          sgn*cosmui(j,m)*sinnv(k,n1))
                     bmn(mn) = bmn(mn) + tsini*bsqa(jlk)
                     gmn(mn) = gmn(mn) + tsini*gsqrta(jlk,0)
                     bsubumn(mn) = bsubumn(mn) + tsini*bsubua(jlk)
                     bsubvmn(mn) = bsubvmn(mn) + tsini*bsubva(jlk)
                     bsubsmn(mn) = bsubsmn(mn) + tcosi*bsubsa(jlk)
                     bsupumn(mn) = bsupumn(mn) + tsini*bsupua(jlk)
                     bsupvmn(mn) = bsupvmn(mn) + tsini*bsupva(jlk)

                     jlk = jlk+ns
                  END DO
               END DO
            END DO MN3

            gmns(:,js) = gmn(:)
            bmns(:,js) = bmn(:)
            bsubumns(:,js) = bsubumn(:)
            bsubvmns(:,js) = bsubvmn(:)
            bsubsmnc(:,js) = bsubsmn(:)
            bsupumns(:,js) = bsupumn(:)
            bsupvmns(:,js) = bsupvmn(:)
         END DO RADIUS3

         gmns(:,1) = 0; bmns(:,1) = 0
         bsubumns(:,1) = 0
         bsubvmns(:,1) = 0
         bsubsmnc(:,1) = 2*bsubsmnc(:,2) - bsubsmnc(:,3)
         bsupumns(:,1) = 0;  bsupvmns(:,1) = 0

         IF (lfreeb) THEN        !MRC  10-15-15
            bsubumns_sur = 0
            bsubvmns_sur = 0
            bsupumns_sur = 0
            bsupvmns_sur = 0

            DO mn = 1, mnmax_nyq0
               n = NINT(xn_nyq0(mn))/nfp
               m = NINT(xm_nyq0(mn))
               n1 = ABS(n)
               dmult = mscale(m)*nscale(n1)*tmult
               IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
               sgn = SIGN(1, n)
               lk = 0
               DO j = 1, ntheta2
                  DO k = 1, nzeta
                     lk = lk + 1
                     tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     1                          sgn*cosmui(j,m)*sinnv(k,n1))
                     bsubumns_sur(mn) = bsubumns_sur(mn)                          &
     &                                + tsini*bsubua_sur(lk)
                     bsubvmns_sur(mn) = bsubvmns_sur(mn)                          &
     &                                + tsini*bsubva_sur(lk)
                     bsupumns_sur(mn) = bsupumns_sur(mn)                          &
     &                                + tsini*bsupua_sur(lk)
                     bsupvmns_sur(mn) = bsupvmns_sur(mn)                          &
     &                                + tsini*bsupva_sur(lk)
                  END DO
               END DO
            END DO
         END IF
      end if ! lasym

      CALL Compute_Currents(bsubsmnc, bsubsmns, bsubumnc, bsubumns,            &
     &                      bsubvmnc, bsubvmns,                                &
     &                      xm_nyq0, xn_nyq0, mnmax_nyq0, lasym, ns,           &
     &                      currumnc, currvmnc, currumns, currvmns)

#ifdef _DEBUG
      WRITE (333, *) '    JS     M*B_S     GRAD(B_U)    J^V'
      DO mn = 1, mnmax_nyq0
         WRITE (333,'(2(a,i4))') 'm=',  INT(xm_nyq0(mn)),
     1                           ' n=', INT(xn_nyq0(mn))/nfp
         DO js = 2,ns-1
            tmult=-xm_nyq0(mn)*bsubsmns(mn,js) +
     1                    ohs*(bsubumnc(mn,js+1)-bsubumnc(mn,js))
            WRITE (333,'(i6,1p,3e12.4)') js,
     1                  bsubsmns(mn,js)*xm_nyq0(mn),
     2             ohs*(bsubumnc(mn,js+1)-bsubumnc(mn,js)),
     3             tmult
         END DO
      END DO

      WRITE(333,*) version_
      IF (lasym) THEN
         WRITE(333,2002) 'mn', 'rmnc', 'rmns', 'zmnc', 'zmns',                 &
     &                         'lmnc', 'lmns', 'gmnc', 'gmns',                 &
     &                         'bmnc', 'bmns',                                 &
     &                         'bsubumnc', 'bsubumns',                         &
     &                         'bsubvmnc', 'bsubvmns',                         &
     &                         'bsubsmnc', 'bsubsmns',                         &
     &                         'bsupumnc', 'bsupumns',                         &
     &                         'bsupvmnc', 'bsupvmns'
      ELSE
         WRITE(333,2000) 'mn', 'rmnc', 'lmns', 'gmnc', 'bmnc',                 &
     &                         'bsubumnc', 'bsubvmnc',                         &
     &                         'bsubsmns',                                     &
     &                         'bsupumnc', 'bsupvmnc'
      END IF
      DO mn = 1, mnmax
         IF (lasym) THEN
            WRITE(333,2003) mn, rmnc(mn,ns/2), rmns(mn,ns/2),                  &
     &                          zmnc(mn,ns/2), zmns(mn,ns/2),                  &
     &                          lmnc(mn,ns/2), lmns(mn,ns/2),                  &
     &                          gmnc(mn,ns/2), gmns(mn,ns/2),                  &
     &                          bmnc(mn,ns/2), bmns(mn,ns/2),                  &
     &                          bsubumnc(mn,ns/2), bsubumns(mn,ns/2),          &
     &                          bsubvmnc(mn,ns/2), bsubvmns(mn,ns/2),          &
     &                          bsubsmnc(mn,ns/2), bsubsmns(mn,ns/2),          &
     &                          bsupumnc(mn,ns/2), bsupumns(mn,ns/2),          &
     &                          bsupvmnc(mn,ns/2), bsupvmns(mn,ns/2)
         ELSE
            WRITE(333,2001) mn, rmnc(mn,ns/2), lmns(mn,ns/2),                  &
     &                          gmnc(mn,ns/2), bmnc(mn,ns/2),                  &
     &                          bsubumnc(mn,ns/2), bsubvmnc(mn,ns/2),          &
     &                          bsubsmns(mn,ns/2),                             &
     &                          bsupumnc(mn,ns/2), bsupvmnc(mn,ns/2)
         END IF
      END DO
2000  FORMAT(a2,10(2x,a12))
2001  FORMAT(i2,10(2x,e12.5))
2002  FORMAT(a2,20(2x,a12))
2003  FORMAT(i2,20(2x,es12.5))
#endif
!
!     WRITE OUT ARRAYS
!
      CALL cdf_write(nwout, vn_racc, raxis_cc(0:ntor))
      CALL cdf_write(nwout, vn_zacs, zaxis_cs(0:ntor))
      CALL cdf_write(nwout, vn_rmnc, rmnc)
      CALL cdf_write(nwout, vn_zmns, zmns)
      CALL cdf_write(nwout, vn_lmns, lmns)
      CALL cdf_write(nwout, vn_gmnc, gmnc)              !Half mesh
      CALL cdf_write(nwout, vn_bmnc, bmnc)              !Half mesh
      CALL cdf_write(nwout, vn_bsubumnc, bsubumnc)      !Half mesh
      CALL cdf_write(nwout, vn_bsubvmnc, bsubvmnc)      !Half mesh
      CALL cdf_write(nwout, vn_bsubsmns, bsubsmns)      !Full mesh

      CALL cdf_write(nwout, vn_currumnc, currumnc)      !MRK 8-12-16
      CALL cdf_write(nwout, vn_currvmnc, currvmnc)

!     GET RID OF THESE EVENTUALLY: DON'T NEED THEM (can express in terms of lambdas)
      CALL cdf_write(nwout, vn_bsupumnc, bsupumnc)
      CALL cdf_write(nwout, vn_bsupvmnc, bsupvmnc)

      IF (lfreeb) THEN        !MRC    10-15-15
         CALL cdf_write(nwout, vn_bsubumnc_sur, bsubumnc_sur)
         CALL cdf_write(nwout, vn_bsubvmnc_sur, bsubvmnc_sur)
         CALL cdf_write(nwout, vn_bsupumnc_sur, bsupumnc_sur)
         CALL cdf_write(nwout, vn_bsupvmnc_sur, bsupvmnc_sur)
      END IF

!     FULL-MESH quantities
!     NOTE: jdotb is in units_of_A (1/mu0 incorporated in jxbforce...)
!     prior to version 6.00, this was output in internal VMEC units...

      j = SIZE(am)-1
      CALL cdf_write(nwout, vn_am, am(0:j))
      j = SIZE(ac)-1
      CALL cdf_write(nwout, vn_ac, ac(0:j))
      j = SIZE(ai)-1
      CALL cdf_write(nwout, vn_ai, ai(0:j))

      j = SIZE(am_aux_s)
      CALL cdf_write(nwout, vn_am_aux_s, am_aux_s(1:j))
      j = SIZE(am_aux_f)
      CALL cdf_write(nwout, vn_am_aux_f, am_aux_f(1:j))
      j = SIZE(ac_aux_s)
      CALL cdf_write(nwout, vn_ac_aux_s, ac_aux_s(1:j))
      j = SIZE(ac_aux_f)
      CALL cdf_write(nwout, vn_ac_aux_f, ac_aux_f(1:j))
      j = SIZE(ai_aux_s)
      CALL cdf_write(nwout, vn_ai_aux_s, ai_aux_s(1:j))
      j = SIZE(ai_aux_f)
      CALL cdf_write(nwout, vn_ai_aux_f, ai_aux_f(1:j))

      CALL cdf_write(nwout, vn_iotaf, iotaf(1:ns))
      CALL cdf_write(nwout, vn_qfact, qfact(1:ns))
      CALL cdf_write(nwout, vn_presf, presf/mu0)
      CALL cdf_write(nwout, vn_phi, phi)
      CALL cdf_write(nwout, vn_phipf, twopi*signgs*phipf)
      CALL cdf_write(nwout, vn_chi, chi)
      CALL cdf_write(nwout, vn_chipf, twopi*signgs*chipf)
      CALL cdf_write(nwout, vn_jcuru, jcuru/mu0)
      CALL cdf_write(nwout, vn_jcurv, jcurv/mu0)
      CALL cdf_write(nwout, vn_jdotb, jdotb)
      CALL cdf_write(nwout, vn_bdotb, bdotb)
      CALL cdf_write(nwout, vn_bgrv, bdotgradv)

!     HALF-MESH quantities
      iotas(1) = 0; mass(1) = 0; pres(1) = 0; phip(1) = 0;
      buco(1) = 0; bvco(1) = 0; vp(1) = 0; overr(1) = 0;  specw(1) = 1
      beta_vol(1) = 0
      CALL cdf_write(nwout, vn_iotah, iotas(1:ns))
      CALL cdf_write(nwout, vn_mass, mass/mu0)
      CALL cdf_write(nwout, vn_presh, pres(1:ns)/mu0)
      CALL cdf_write(nwout, vn_betah, beta_vol)
      CALL cdf_write(nwout, vn_buco, buco)
      CALL cdf_write(nwout, vn_bvco, bvco)
      CALL cdf_write(nwout, vn_vp, vp(1:ns))
      CALL cdf_write(nwout, vn_specw, specw)
      CALL cdf_write(nwout, vn_phip, phips(1:ns))
      CALL cdf_write(nwout, vn_overr, overr(1:ns))

!     MERCIER_CRITERION
      CALL cdf_write(nwout, vn_merc, Dmerc)
      CALL cdf_write(nwout, vn_mshear, Dshear)
      CALL cdf_write(nwout, vn_mwell, Dwell)
      CALL cdf_write(nwout, vn_mcurr, Dcurr)
      CALL cdf_write(nwout, vn_mgeo, Dgeod)
      CALL cdf_write(nwout, vn_equif, equif)

      CALL cdf_write(nwout, vn_fsq, fsqt(1:nstore_seq))
      CALL cdf_write(nwout, vn_wdot, wdot(1:nstore_seq))

!--------------------DEC$ ENDIF

      IF (lasym) THEN
         CALL cdf_write(nwout, vn_racs, raxis_cs(0:ntor))
         CALL cdf_write(nwout, vn_zacc, zaxis_cc(0:ntor))
         CALL cdf_write(nwout, vn_rmns, rmns)
         CALL cdf_write(nwout, vn_zmnc, zmnc)
         CALL cdf_write(nwout, vn_lmnc, lmnc)
         CALL cdf_write(nwout, vn_gmns, gmns)
         CALL cdf_write(nwout, vn_bmns, bmns)
         CALL cdf_write(nwout, vn_bsubumns, bsubumns)
         CALL cdf_write(nwout, vn_bsubvmns, bsubvmns)
         CALL cdf_write(nwout, vn_bsubsmnc, bsubsmnc)

         CALL cdf_write(nwout, vn_currumns, currumns)     !MRC  8-12-16
         CALL cdf_write(nwout, vn_currvmns, currvmns)

!     GET RID OF THESE EVENTUALLY: DON'T NEED THEM
         CALL cdf_write(nwout, vn_bsupumns, bsupumns)
         CALL cdf_write(nwout, vn_bsupvmns, bsupvmns)

         IF (lfreeb) THEN     !MRC    10-15-15
            CALL cdf_write(nwout, vn_bsubumns_sur, bsubumns_sur)
            CALL cdf_write(nwout, vn_bsubvmns_sur, bsubvmns_sur)
            CALL cdf_write(nwout, vn_bsupumns_sur, bsupumns_sur)
            CALL cdf_write(nwout, vn_bsupvmns_sur, bsupvmns_sur)
         END IF
      END IF

      CALL cdf_close(nwout)

!
!     RESTORE nyq ENDPOINT VALUES
!
      IF (mnyq .ne. 0) cosmui(:,mnyq) = 2*cosmui(:,mnyq)
      IF (nnyq .ne. 0) cosnv (:,nnyq) = 2*cosnv (:,nnyq)

!
! DEALLOCATIONS ! J Geiger: these have been moved downwards.
!
      IF (ALLOCATED(gmnc)) DEALLOCATE(gmnc, bmnc, bsubumnc, bsubvmnc,
     1                                bsubsmns, bsupumnc, bsupvmnc
     3                                )
      IF (ALLOCATED(gmns)) DEALLOCATE(gmns, bmns, bsubumns, bsubvmns,
     1                                bsubsmnc, bsupumns, bsupvmns
     3                                )
! J Geiger: check also for allocation.
      IF (ALLOCATED(gmn)) DEALLOCATE (gmn, bmn, bsubumn, bsubvmn,
     1             bsubsmn, bsupumn, bsupvmn,
     4             stat=istat)

      IF (ALLOCATED(bsubumnc_sur)) THEN
         DEALLOCATE(bsubumnc_sur, bsubvmnc_sur)
         DEALLOCATE(bsupumnc_sur, bsupvmnc_sur)
      END IF
      IF (ALLOCATED(bsubumns_sur)) THEN
         DEALLOCATE(bsubumns_sur, bsubvmns_sur)
         DEALLOCATE(bsupumns_sur, bsupvmns_sur)
      END IF
      IF (ALLOCATED(bsubua_sur)) THEN
         DEALLOCATE(bsubua_sur, bsubva_sur)
         DEALLOCATE(bsupua_sur, bsupva_sur)
      END IF

      rzl_array = 0

      CALL second0 (twoutoff)
      timer(twout) = timer(twout) + twoutoff - twouton
      fo_wrout_time = timer(twout)

      END SUBROUTINE wrout