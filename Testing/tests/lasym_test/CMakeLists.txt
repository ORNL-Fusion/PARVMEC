#  Copy input file from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test_sym.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test_sym.vmec
                COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test_asym.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test_asym.vmec
                COPYONLY)

#  Test Fixed boundary equilibrium.
add_test (NAME    vmec_sym_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test_sym.vmec)
set_tests_properties (vmec_sym_test
                      PROPERTIES
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)
add_test (NAME    vmec_asym_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test_asym.vmec)
set_tests_properties (vmec_asym_test
                      PROPERTIES
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)

#  Check woutfile to ensure that values match.
add_test (NAME    vmec_asym_check_aspect_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=aspect -tol=5.0E-16)
set_tests_properties (vmec_asym_check_aspect_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_b0_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=b0 -tol=2.2E-15)
set_tests_properties (vmec_asym_check_b0_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_bdotb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=bdotb -tol=4.6E-15)
set_tests_properties (vmec_asym_check_bdotb_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_bdotgradv_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=bdotgradv -tol=2.1E-14)
set_tests_properties (vmec_asym_check_bdotgradv_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_beta_vol_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=beta_vol -tol=1.1E-15)
set_tests_properties (vmec_asym_check_beta_vol_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_betapol_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=betapol -tol=1.1E-13)
set_tests_properties (vmec_asym_check_betapol_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_betator_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=betator -tol=2.5E-17)
set_tests_properties (vmec_asym_check_betator_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_betatotal_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=betatotal -tol=2.0E-17)
set_tests_properties (vmec_asym_check_betatotal_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_betaxis_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=betaxis -tol=3.1E-16)
set_tests_properties (vmec_asym_check_betaxis_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_bmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=bmnc -tol=9.0E-15)
set_tests_properties (vmec_asym_check_bmnc_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_bsubsmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=bsubsmns -tol=2.1E-15)
set_tests_properties (vmec_asym_check_bsubsmns_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_bsubumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=bsubumnc -tol=2.2E-15)
set_tests_properties (vmec_asym_check_bsubumnc_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_bsubvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=bsubvmnc -tol=6.5E-15)
set_tests_properties (vmec_asym_check_bsubvmnc_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_bsupumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=bsupumnc -tol=3.7E-14)
set_tests_properties (vmec_asym_check_bsupumnc_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_bsupvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=bsupvmnc -tol=2.0E-14)
set_tests_properties (vmec_asym_check_bsupvmnc_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_buco_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=buco -tol=4.4E-18)
set_tests_properties (vmec_asym_check_buco_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_bvco_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=bvco -tol=4.0E-15)
set_tests_properties (vmec_asym_check_bvco_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_chi_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=chi -tol=9.4E-17)
set_tests_properties (vmec_asym_check_chi_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_chipf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=chipf -tol=1.3E-15)
set_tests_properties (vmec_asym_check_chipf_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_currumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=currumnc -tol=2.9E-7)
set_tests_properties (vmec_asym_check_currumnc_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_currvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=currvmnc -tol=6.9E-8)
set_tests_properties (vmec_asym_check_currvmnc_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_DCurr_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=DCurr -tol=5.1E-11)
set_tests_properties (vmec_asym_check_DCurr_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_DGeod_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=DGeod -tol=3.7E-11)
set_tests_properties (vmec_asym_check_DGeod_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_DMerc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=DMerc -tol=3.7E-11)
set_tests_properties (vmec_asym_check_DMerc_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_DShear_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=DShear -tol=3.0E-11)
set_tests_properties (vmec_asym_check_DShear_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_DWell_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=DWell -tol=1.0E-20)
set_tests_properties (vmec_asym_check_DWell_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_fsql_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=fsql -tol=1.0E-20)
set_tests_properties (vmec_asym_check_fsql_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_fsqr_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=fsqr -tol=1.0E-20)
set_tests_properties (vmec_asym_check_fsqr_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_fsqt_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=fsqt -tol=5.0E-19)
set_tests_properties (vmec_asym_check_fsqt_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_fsqz_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=fsqz -tol=1.0E-20)
set_tests_properties (vmec_asym_check_fsqz_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_gmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=gmnc -tol=4.7E-15)
set_tests_properties (vmec_asym_check_gmnc_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_IonLarmor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=IonLarmor -tol=6.0E-15)
set_tests_properties (vmec_asym_check_IonLarmor_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_iotaf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=iotaf -tol=2.5E-14)
set_tests_properties (vmec_asym_check_iotaf_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_jcuru_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=jcuru -tol=2.0E-7)
set_tests_properties (vmec_asym_check_jcuru_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_jcurv_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=jcurv -tol=7.8E-11)
set_tests_properties (vmec_asym_check_jcurv_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_jdotb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=jdotb -tol=7.0E-8)
set_tests_properties (vmec_asym_check_jdotb_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_lmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=lmns -tol=2.0E-13)
set_tests_properties (vmec_asym_check_lmns_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_mnmax_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=mnmax -tol=1.0E-20)
set_tests_properties (vmec_asym_check_mnmax_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_mnmax_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=mnmax_nyq -tol=1.0E-20)
set_tests_properties (vmec_asym_check_mnmax_nyq_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_phipf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=phipf -tol=1.0E-20)
set_tests_properties (vmec_asym_check_phipf_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_phips_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=phips -tol=1.0E-20)
set_tests_properties (vmec_asym_check_phips_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_pres_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=pres -tol=1.0E-20)
set_tests_properties (vmec_asym_check_pres_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_presf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=presf -tol=1.0E-20)
set_tests_properties (vmec_asym_check_presf_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_q_factor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=q_factor -tol=2.1E-13)
set_tests_properties (vmec_asym_check_q_factor_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_rbtor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=rbtor -tol=3.0E-15)
set_tests_properties (vmec_asym_check_rbtor_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_rbtor0_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=rbtor0 -tol=6.0E-16)
set_tests_properties (vmec_asym_check_rbtor0_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_rmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=rmnc -tol=7.3E-15)
set_tests_properties (vmec_asym_check_rmnc_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_specw_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=specw -tol=4.5E-16)
set_tests_properties (vmec_asym_check_specw_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_volavgB_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=volavgB -tol=2.0E-16)
set_tests_properties (vmec_asym_check_volavgB_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_volume_p_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=volume_p -tol=3.0E-16)
set_tests_properties (vmec_asym_check_volume_p_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_vp_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=vp -tol=1.5E-15)
set_tests_properties (vmec_asym_check_vp_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_wb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=wb -tol=6.6E-19)
set_tests_properties (vmec_asym_check_wb_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_wdot_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=wdot -tol=5.9E-16)
set_tests_properties (vmec_asym_check_wdot_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_wp_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=wp -tol=1.0E-20)
set_tests_properties (vmec_asym_check_wp_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_xm_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=xm -tol=1.0E-20)
set_tests_properties (vmec_asym_check_xm_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_xm_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=xm_nyq -tol=1.0E-20)
set_tests_properties (vmec_asym_check_xm_nyq_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_xn_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=xn -tol=1.0E-20)
set_tests_properties (vmec_asym_check_xn_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_xn_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=xn_nyq -tol=1.0E-20)
set_tests_properties (vmec_asym_check_xn_nyq_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
add_test (NAME    vmec_asym_check_zmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_sym.vmec.nc -wout_file2=wout_test_asym.vmec.nc -quantity=zmns -tol=2.7E-15)
set_tests_properties (vmec_asym_check_zmns_test
                      PROPERTIES
                      DEPENDS "vmec_sym_test;vmec_asym_test")
