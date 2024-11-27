#  Copy input file from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test_init.vmec
                COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test_reset.vmec
                COPYONLY)

#  Test Fixed boundary equilibrium.
add_test (NAME    vmec_reset_test_init
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test_init.vmec)
add_test (NAME    vmec_reset_test_reset
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test_reset.vmec reset=wout_test_init.vmec.nc)
set_tests_properties (vmec_reset_test_reset
                      PROPERTIES
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)

#  Check woutfile to ensure that values match.
add_test (NAME    vmec_reset_check_aspect_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=aspect -tol=1.0E-20)
set_tests_properties (vmec_reset_check_aspect_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_b0_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=b0 -tol=1.0E-20)
set_tests_properties (vmec_reset_check_b0_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_bdotb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=bdotb -tol=1.0E-20)
set_tests_properties (vmec_reset_check_bdotb_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_bdotgradv_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=bdotgradv -tol=1.0E-20)
set_tests_properties (vmec_reset_check_bdotgradv_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_beta_vol_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=beta_vol -tol=1.0E-20)
set_tests_properties (vmec_reset_check_beta_vol_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_betapol_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=betapol -tol=1.0E-20)
set_tests_properties (vmec_reset_check_betapol_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_betator_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=betator -tol=1.0E-20)
set_tests_properties (vmec_reset_check_betator_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_betatotal_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=betatotal -tol=1.0E-20)
set_tests_properties (vmec_reset_check_betatotal_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_betaxis_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=betaxis -tol=1.0E-20)
set_tests_properties (vmec_reset_check_betaxis_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_bmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=bmnc -tol=1.0E-20)
set_tests_properties (vmec_reset_check_bmnc_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_bsubsmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=bsubsmns -tol=1.0E-20)
set_tests_properties (vmec_reset_check_bsubsmns_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_bsubumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=bsubumnc -tol=1.0E-20)
set_tests_properties (vmec_reset_check_bsubumnc_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_bsubvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=bsubvmnc -tol=1.0E-20)
set_tests_properties (vmec_reset_check_bsubvmnc_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_bsupumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=bsupumnc -tol=1.0E-20)
set_tests_properties (vmec_reset_check_bsupumnc_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_bsupvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=bsupvmnc -tol=1.0E-20)
set_tests_properties (vmec_reset_check_bsupvmnc_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_buco_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=buco -tol=1.0E-20)
set_tests_properties (vmec_reset_check_buco_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_bvco_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=bvco -tol=1.0E-20)
set_tests_properties (vmec_reset_check_bvco_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_chi_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=chi -tol=1.0E-20)
set_tests_properties (vmec_reset_check_chi_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_chipf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=chipf -tol=1.0E-20)
set_tests_properties (vmec_reset_check_chipf_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_currumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=currumnc -tol=1.0E-20)
set_tests_properties (vmec_reset_check_currumnc_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_currvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=currvmnc -tol=1.0E-20)
set_tests_properties (vmec_reset_check_currvmnc_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_DCurr_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=DCurr -tol=1.0E-20)
set_tests_properties (vmec_reset_check_DCurr_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_DGeod_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=DGeod -tol=1.0E-20)
set_tests_properties (vmec_reset_check_DGeod_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_DMerc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=DMerc -tol=1.0E-20)
set_tests_properties (vmec_reset_check_DMerc_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_DShear_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=DShear -tol=1.0E-20)
set_tests_properties (vmec_reset_check_DShear_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_DWell_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=DWell -tol=1.0E-20)
set_tests_properties (vmec_reset_check_DWell_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_fsql_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=fsql -tol=1.0E-20)
set_tests_properties (vmec_reset_check_fsql_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_fsqr_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=fsqr -tol=1.0E-20)
set_tests_properties (vmec_reset_check_fsqr_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_fsqt_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=fsqt -tol=1.0E-20)
set_tests_properties (vmec_reset_check_fsqt_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_fsqz_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=fsqz -tol=1.0E-20)
set_tests_properties (vmec_reset_check_fsqz_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_gmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=gmnc -tol=1.0E-20)
set_tests_properties (vmec_reset_check_gmnc_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_IonLarmor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=IonLarmor -tol=1.0E-20)
set_tests_properties (vmec_reset_check_IonLarmor_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_iotaf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=iotaf -tol=1.0E-20)
set_tests_properties (vmec_reset_check_iotaf_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_jcuru_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=jcuru -tol=1.0E-20)
set_tests_properties (vmec_reset_check_jcuru_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_jcurv_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=jcurv -tol=1.0E-20)
set_tests_properties (vmec_reset_check_jcurv_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_jdotb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=jdotb -tol=1.0E-20)
set_tests_properties (vmec_reset_check_jdotb_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_lmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=lmns -tol=1.0E-20)
set_tests_properties (vmec_reset_check_lmns_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_mnmax_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=mnmax -tol=1.0E-20)
set_tests_properties (vmec_reset_check_mnmax_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_mnmax_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=mnmax_nyq -tol=1.0E-20)
set_tests_properties (vmec_reset_check_mnmax_nyq_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_phipf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=phipf -tol=1.0E-20)
set_tests_properties (vmec_reset_check_phipf_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_phips_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=phips -tol=1.0E-20)
set_tests_properties (vmec_reset_check_phips_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_pres_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=pres -tol=1.0E-20)
set_tests_properties (vmec_reset_check_pres_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_presf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=presf -tol=1.0E-20)
set_tests_properties (vmec_reset_check_presf_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_q_factor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=q_factor -tol=1.0E-20)
set_tests_properties (vmec_reset_check_q_factor_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_rbtor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=rbtor -tol=1.0E-20)
set_tests_properties (vmec_reset_check_rbtor_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_rbtor0_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=rbtor0 -tol=1.0E-20)
set_tests_properties (vmec_reset_check_rbtor0_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_rmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=rmnc -tol=1.0E-20)
set_tests_properties (vmec_reset_check_rmnc_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_specw_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=specw -tol=1.0E-20)
set_tests_properties (vmec_reset_check_specw_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_volavgB_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=volavgB -tol=1.0E-20)
set_tests_properties (vmec_reset_check_volavgB_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_volume_p_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=volume_p -tol=1.0E-20)
set_tests_properties (vmec_reset_check_volume_p_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_vp_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=vp -tol=1.0E-20)
set_tests_properties (vmec_reset_check_vp_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_wb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=wb -tol=1.0E-20)
set_tests_properties (vmec_reset_check_wb_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_wdot_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=wdot -tol=1.0E-20)
set_tests_properties (vmec_reset_check_wdot_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_wp_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=wp -tol=1.0E-20)
set_tests_properties (vmec_reset_check_wp_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_xm_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=xm -tol=1.0E-20)
set_tests_properties (vmec_reset_check_xm_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_xm_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=xm_nyq -tol=1.0E-20)
set_tests_properties (vmec_reset_check_xm_nyq_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_xn_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=xn -tol=1.0E-20)
set_tests_properties (vmec_reset_check_xn_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_xn_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=xn_nyq -tol=1.0E-20)
set_tests_properties (vmec_reset_check_xn_nyq_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
add_test (NAME    vmec_reset_check_zmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_init.vmec.nc -wout_file2=wout_test_reset.vmec.nc -quantity=zmns -tol=1.0E-20)
set_tests_properties (vmec_reset_check_zmns_test
                      PROPERTIES
                      DEPENDS "vmec_reset_test_init;vmec_reset_test_reset")
