#  Copy input file from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test_forward.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test_forward.vmec
                COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test_reversed.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test_reversed.vmec
                COPYONLY)

#  Test an asymmetric boundary with theta running either way around it.
add_test (NAME    vmec_forward_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test_forward.vmec)
set_tests_properties (vmec_forward_test
                      PROPERTIES
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)
add_test (NAME    vmec_reversed_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test_reversed.vmec)
set_tests_properties (vmec_reversed_test
                      PROPERTIES
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)

#  Check woutfile to ensure that values match.
add_test (NAME    vmec_reversed_check_aspect_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=aspect -tol=1.4E-13)
set_tests_properties (vmec_reversed_check_aspect_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_b0_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=b0 -tol=4.0E-13)
set_tests_properties (vmec_reversed_check_b0_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bdotb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bdotb -tol=3.9E-12)
set_tests_properties (vmec_reversed_check_bdotb_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bdotgradv_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bdotgradv -tol=2.1E-11)
set_tests_properties (vmec_reversed_check_bdotgradv_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_beta_vol_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=beta_vol -tol=8.5E-13)
set_tests_properties (vmec_reversed_check_beta_vol_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_betapol_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=betapol -tol=4.5E-12)
set_tests_properties (vmec_reversed_check_betapol_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_betator_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=betator -tol=1.1E-15)
set_tests_properties (vmec_reversed_check_betator_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_betatotal_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=betatotal -tol=3.5E-16)
set_tests_properties (vmec_reversed_check_betatotal_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_betaxis_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=betaxis -tol=4.8E-14)
set_tests_properties (vmec_reversed_check_betaxis_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bmnc -tol=1.5E-12)
set_tests_properties (vmec_reversed_check_bmnc_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bmns -tol=9.0E-13)
set_tests_properties (vmec_reversed_check_bmns_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bsubsmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bsubsmnc -tol=2.6E-12)
set_tests_properties (vmec_reversed_check_bsubsmnc_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bsubsmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bsubsmns -tol=2.2E-12)
set_tests_properties (vmec_reversed_check_bsubsmns_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bsubumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bsubumnc -tol=6.0E-13)
set_tests_properties (vmec_reversed_check_bsubumnc_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bsubumns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bsubumns -tol=8.2E-13)
set_tests_properties (vmec_reversed_check_bsubumns_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bsubvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bsubvmnc -tol=7.1E-13)
set_tests_properties (vmec_reversed_check_bsubvmnc_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bsubvmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bsubvmns -tol=5.8E-13)
set_tests_properties (vmec_reversed_check_bsubvmns_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bsupumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bsupumnc -tol=1.7E-11)
set_tests_properties (vmec_reversed_check_bsupumnc_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bsupumns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bsupumns -tol=2.2E-11)
set_tests_properties (vmec_reversed_check_bsupumns_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bsupvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bsupvmnc -tol=2.8E-12)
set_tests_properties (vmec_reversed_check_bsupvmnc_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bsupvmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bsupvmns -tol=2.0E-12)
set_tests_properties (vmec_reversed_check_bsupvmns_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_buco_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=buco -tol=3.5E-16)
set_tests_properties (vmec_reversed_check_buco_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_bvco_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=bvco -tol=1.7E-13)
set_tests_properties (vmec_reversed_check_bvco_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_chi_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=chi -tol=7.5E-14)
set_tests_properties (vmec_reversed_check_chi_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_chipf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=chipf -tol=1.1E-12)
set_tests_properties (vmec_reversed_check_chipf_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_currumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=currumnc -tol=3.3E-5)
set_tests_properties (vmec_reversed_check_currumnc_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_currumns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=currumns -tol=1.7E-5)
set_tests_properties (vmec_reversed_check_currumns_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_currvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=currvmnc -tol=1.1E-5)
set_tests_properties (vmec_reversed_check_currvmnc_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_currvmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=currvmns -tol=1.7E-5)
set_tests_properties (vmec_reversed_check_currvmns_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_DCurr_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=DCurr -tol=3.2E-8)
set_tests_properties (vmec_reversed_check_DCurr_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_DGeod_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=DGeod -tol=2.2E-9)
set_tests_properties (vmec_reversed_check_DGeod_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_DMerc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=DMerc -tol=1.5E-8)
set_tests_properties (vmec_reversed_check_DMerc_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_DShear_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=DShear -tol=1.9E-8)
set_tests_properties (vmec_reversed_check_DShear_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_DWell_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=DWell -tol=1.0E-20)
set_tests_properties (vmec_reversed_check_DWell_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_fsql_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=fsql -tol=7.7E-24)
set_tests_properties (vmec_reversed_check_fsql_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_fsqr_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=fsqr -tol=6.9E-22)
set_tests_properties (vmec_reversed_check_fsqr_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_fsqt_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=fsqt -tol=3.1E-16)
set_tests_properties (vmec_reversed_check_fsqt_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_fsqz_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=fsqz -tol=1.4E-21)
set_tests_properties (vmec_reversed_check_fsqz_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_gmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=gmnc -tol=4.4E-12)
set_tests_properties (vmec_reversed_check_gmnc_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_gmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=gmns -tol=4.2E-12)
set_tests_properties (vmec_reversed_check_gmns_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_IonLarmor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=IonLarmor -tol=3.5E-16)
set_tests_properties (vmec_reversed_check_IonLarmor_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_iotaf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=iotaf -tol=2.1E-11)
set_tests_properties (vmec_reversed_check_iotaf_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_jcuru_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=jcuru -tol=3.3E-6)
set_tests_properties (vmec_reversed_check_jcuru_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_jcurv_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=jcurv -tol=5.8E-9)
set_tests_properties (vmec_reversed_check_jcurv_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_jdotb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=jdotb -tol=1.4E-6)
set_tests_properties (vmec_reversed_check_jdotb_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_lmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=lmnc -tol=2.2E-10)
set_tests_properties (vmec_reversed_check_lmnc_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_lmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=lmns -tol=1.2E-10)
set_tests_properties (vmec_reversed_check_lmns_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_mnmax_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=mnmax -tol=1.0E-20)
set_tests_properties (vmec_reversed_check_mnmax_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_mnmax_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=mnmax_nyq -tol=1.0E-20)
set_tests_properties (vmec_reversed_check_mnmax_nyq_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_phipf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=phipf -tol=1.0E-20)
set_tests_properties (vmec_reversed_check_phipf_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_phips_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=phips -tol=1.0E-20)
set_tests_properties (vmec_reversed_check_phips_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_pres_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=pres -tol=1.0E-20)
set_tests_properties (vmec_reversed_check_pres_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_presf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=presf -tol=1.0E-20)
set_tests_properties (vmec_reversed_check_presf_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_q_factor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=q_factor -tol=2.4E-10)
set_tests_properties (vmec_reversed_check_q_factor_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_rbtor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=rbtor -tol=1.8E-13)
set_tests_properties (vmec_reversed_check_rbtor_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_rbtor0_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=rbtor0 -tol=2.0E-14)
set_tests_properties (vmec_reversed_check_rbtor0_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_rmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=rmnc -tol=2.7E-12)
set_tests_properties (vmec_reversed_check_rmnc_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_rmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=rmns -tol=1.6E-12)
set_tests_properties (vmec_reversed_check_rmns_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_specw_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=specw -tol=1.8E-12)
set_tests_properties (vmec_reversed_check_specw_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_volavgB_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=volavgB -tol=5.6E-15)
set_tests_properties (vmec_reversed_check_volavgB_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_volume_p_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=volume_p -tol=2.3E-14)
set_tests_properties (vmec_reversed_check_volume_p_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_vp_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=vp -tol=1.3E-12)
set_tests_properties (vmec_reversed_check_vp_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_wb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=wb -tol=4.4E-17)
set_tests_properties (vmec_reversed_check_wb_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_wdot_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=wdot -tol=3.0E-14)
set_tests_properties (vmec_reversed_check_wdot_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_wp_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=wp -tol=3.4E-19)
set_tests_properties (vmec_reversed_check_wp_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_xm_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=xm -tol=1.0E-20)
set_tests_properties (vmec_reversed_check_xm_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_xm_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=xm_nyq -tol=1.0E-20)
set_tests_properties (vmec_reversed_check_xm_nyq_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_xn_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=xn -tol=1.0E-20)
set_tests_properties (vmec_reversed_check_xn_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_xn_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=xn_nyq -tol=1.0E-20)
set_tests_properties (vmec_reversed_check_xn_nyq_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_zmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=zmnc -tol=2.8E-12)
set_tests_properties (vmec_reversed_check_zmnc_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
add_test (NAME    vmec_reversed_check_zmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_forward.vmec.nc -wout_file2=wout_test_reversed.vmec.nc -quantity=zmns -tol=2.3E-12)
set_tests_properties (vmec_reversed_check_zmns_test
                      PROPERTIES
                      DEPENDS "vmec_forward_test;vmec_reversed_test")
