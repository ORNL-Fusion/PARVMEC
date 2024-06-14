#  Copy input file from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/coils.test
                ${CMAKE_CURRENT_BINARY_DIR}/coils.test
                COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test_serial.vmec
                COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test_parallel.vmec
                COPYONLY)

#  Test Free boundary equilibrium.
add_test (NAME    vmec_mgrid_test
          COMMAND $<TARGET_PROPERTY:mgrid,BINARY_DIR>/mgrid coils.test)
add_test (NAME    vmec_free_boundary_serial_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test_serial.vmec)
set_tests_properties (vmec_free_boundary_serial_test
                      PROPERTIES
                      DEPENDS vmec_mgrid_test)
add_test (NAME    vmec_free_boundary_parallel_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test_parallel.vmec)
set_tests_properties (vmec_free_boundary_parallel_test
                      PROPERTIES
                      DEPENDS vmec_mgrid_test
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)

#  Check woutfile to ensure that values match.
add_test (NAME    vmec_free_boundary_check_aspect_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=aspect -tol=7.4E-14)
set_tests_properties (vmec_free_boundary_check_aspect_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_b0_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=b0 -tol=7.3E-15)
set_tests_properties (vmec_free_boundary_check_b0_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_bdotb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bdotb -tol=1.1E-13)
set_tests_properties (vmec_free_boundary_check_bdotb_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_bdotgradv_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bdotgradv -tol=2.4E-13)
set_tests_properties (vmec_free_boundary_check_bdotgradv_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_beta_vol_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=beta_vol -tol=2.2E-16)
set_tests_properties (vmec_free_boundary_check_beta_vol_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_betapol_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=betapol -tol=6.6E-15)
set_tests_properties (vmec_free_boundary_check_betapol_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_betator_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=betator -tol=2.0E-17)
set_tests_properties (vmec_free_boundary_check_betator_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_betatotal_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=betatotal -tol=2.1E-17)
set_tests_properties (vmec_free_boundary_check_betatotal_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_betaxis_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=betaxis -tol=7.3E-17)
set_tests_properties (vmec_free_boundary_check_betaxis_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_bmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bmnc -tol=1.7E-13)
set_tests_properties (vmec_free_boundary_check_bmnc_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_bsubsmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bsubsmns -tol=1.6E-13)
set_tests_properties (vmec_free_boundary_check_bsubsmns_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_bsubumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bsubumnc -tol=9.6E-14)
set_tests_properties (vmec_free_boundary_check_bsubumnc_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_bsubvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bsubvmnc -tol=8.4E-14)
set_tests_properties (vmec_free_boundary_check_bsubvmnc_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_bsupumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bsupumnc -tol=8.8E-12)
set_tests_properties (vmec_free_boundary_check_bsupumnc_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_bsupvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bsupvmnc -tol=3.1E-13)
set_tests_properties (vmec_free_boundary_check_bsupvmnc_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_buco_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=buco -tol=5.9E-17)
set_tests_properties (vmec_free_boundary_check_buco_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_bvco_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bvco -tol=2.0E-14)
set_tests_properties (vmec_free_boundary_check_bvco_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_chi_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=chi -tol=2.2E-15)
set_tests_properties (vmec_free_boundary_check_chi_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_chipf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=chipf -tol=3.5E-14)
set_tests_properties (vmec_free_boundary_check_chipf_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_currumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=currumnc -tol=2.4E-6)
set_tests_properties (vmec_free_boundary_check_currumnc_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_currvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=currvmnc -tol=4.8E-7)
set_tests_properties (vmec_free_boundary_check_currvmnc_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_DCurr_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=DCurr -tol=2.9E-9)
set_tests_properties (vmec_free_boundary_check_DCurr_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_DGeod_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=DGeod -tol=2.6E-9)
set_tests_properties (vmec_free_boundary_check_DGeod_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_DMerc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=DMerc -tol=3.7E-9)
set_tests_properties (vmec_free_boundary_check_DMerc_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_DShear_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=DShear -tol=8.5E-10)
set_tests_properties (vmec_free_boundary_check_DShear_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_DWell_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=DWell -tol=3.3E-10)
set_tests_properties (vmec_free_boundary_check_DWell_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_fsql_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=fsql -tol=1.0E-20)
set_tests_properties (vmec_free_boundary_check_fsql_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_fsqr_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=fsqr -tol=1.0E-20)
set_tests_properties (vmec_free_boundary_check_fsqr_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_fsqt_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=fsqt -tol=3.0E-17)
set_tests_properties (vmec_free_boundary_check_fsqt_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_fsqz_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=fsqz -tol=1.0E-20)
set_tests_properties (vmec_free_boundary_check_fsqz_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_gmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=gmnc -tol=7.8E-14)
set_tests_properties (vmec_free_boundary_check_gmnc_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_IonLarmor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=IonLarmor -tol=2.2-16)
set_tests_properties (vmec_free_boundary_check_IonLarmor_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_iotaf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=iotaf -tol=9.9E-13)
set_tests_properties (vmec_free_boundary_check_iotaf_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_jcuru_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=jcuru -tol=5.1E-7)
set_tests_properties (vmec_free_boundary_check_jcuru_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_jcurv_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=jcurv -tol=2.0E-9)
set_tests_properties (vmec_free_boundary_check_jcurv_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_jdotb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=jdotb -tol=1.8E-6)
set_tests_properties (vmec_free_boundary_check_jdotb_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_lmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=lmns -tol=2.5E-12)
set_tests_properties (vmec_free_boundary_check_lmns_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_mnmax_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=mnmax -tol=1.0E-20)
set_tests_properties (vmec_free_boundary_check_mnmax_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_mnmax_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=mnmax_nyq -tol=1.0E-20)
set_tests_properties (vmec_free_boundary_check_mnmax_nyq_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_phipf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=phipf -tol=1.0E-20)
set_tests_properties (vmec_free_boundary_check_phipf_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_phips_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=phips -tol=1.0E-20)
set_tests_properties (vmec_free_boundary_check_phips_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_pres_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=pres -tol=1.0E-20)
set_tests_properties (vmec_free_boundary_check_pres_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_presf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=presf -tol=1.0E-20)
set_tests_properties (vmec_free_boundary_check_presf_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_q_factor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=q_factor -tol=6.2E-13)
set_tests_properties (vmec_free_boundary_check_q_factor_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_rbtor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=rbtor -tol=2.0E-14)
set_tests_properties (vmec_free_boundary_check_rbtor_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_rbtor0_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=rbtor0 -tol=5.4E-15)
set_tests_properties (vmec_free_boundary_check_rbtor0_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_rmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=rmnc -tol=1.5E-13)
set_tests_properties (vmec_free_boundary_check_rmnc_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_specw_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=specw -tol=6.2E-12)
set_tests_properties (vmec_free_boundary_check_specw_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_volavgB_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=volavgB -tol=2.7E-14)
set_tests_properties (vmec_free_boundary_check_volavgB_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_volume_p_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=volume_p -tol=2.9E-14)
set_tests_properties (vmec_free_boundary_check_volume_p_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_vp_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=vp -tol=2.8E-15)
set_tests_properties (vmec_free_boundary_check_vp_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_wb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=wb -tol=1.1E-17)
set_tests_properties (vmec_free_boundary_check_wb_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_wdot_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=wdot -tol=4.2E-15)
set_tests_properties (vmec_free_boundary_check_wdot_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_wp_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=wp -tol=2.1E-20)
set_tests_properties (vmec_free_boundary_check_wp_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_xm_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=xm -tol=1.0E-20)
set_tests_properties (vmec_free_boundary_check_xm_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_xm_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=xm_nyq -tol=1.0E-20)
set_tests_properties (vmec_free_boundary_check_xm_nyq_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_xn_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=xn -tol=1.0E-20)
set_tests_properties (vmec_free_boundary_check_xn_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_xn_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=xn_nyq -tol=1.0E-20)
set_tests_properties (vmec_free_boundary_check_xn_nyq_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
add_test (NAME    vmec_free_boundary_check_zmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=zmns -tol=1.5E-13)
set_tests_properties (vmec_free_boundary_check_zmns_test
                      PROPERTIES
                      DEPENDS "vmec_free_boundary_serial_test;vmec_free_boundary_parallel_test")
