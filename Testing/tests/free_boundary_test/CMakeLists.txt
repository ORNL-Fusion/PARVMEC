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
add_test (NAME    vmec_free_boundary_parallel_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test_parallel.vmec)

#  Check woutfile to ensure that values match.
add_test (NAME    vmec_free_boundary_check_aspect_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=aspect -tol=5.0E-14)
add_test (NAME    vmec_free_boundary_check_b0_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=b0 -tol=7.3E-15)
add_test (NAME    vmec_free_boundary_check_bdotb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bdotb -tol=3.1E-14)
add_test (NAME    vmec_free_boundary_check_bdotgradv_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bdotgradv -tol=7.0E-14)
add_test (NAME    vmec_free_boundary_check_beta_vol_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=beta_vol -tol=6.0E-17)
add_test (NAME    vmec_free_boundary_check_betapol_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=betapol -tol=4.0E-15)
add_test (NAME    vmec_free_boundary_check_betator_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=betator -tol=2.0E-17)
add_test (NAME    vmec_free_boundary_check_betatotal_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=betatotal -tol=6.0E-18)
add_test (NAME    vmec_free_boundary_check_betaxis_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=betaxis -tol=7.0E-17)
add_test (NAME    vmec_free_boundary_check_bmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bmnc -tol=5.0E-14)
add_test (NAME    vmec_free_boundary_check_bsubsmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bsubsmns -tol=5.0E-14)
add_test (NAME    vmec_free_boundary_check_bsubumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bsubumnc -tol=3.0E-14)
add_test (NAME    vmec_free_boundary_check_bsubvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bsubvmnc -tol=3.4E-14)
add_test (NAME    vmec_free_boundary_check_bsupumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bsupumnc -tol=4.0E-12)
add_test (NAME    vmec_free_boundary_check_bsupvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bsupvmnc -tol=1.0E-13)
add_test (NAME    vmec_free_boundary_check_buco_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=buco -tol=3.7E-17)
add_test (NAME    vmec_free_boundary_check_bvco_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=bvco -tol=2.0E-14)
add_test (NAME    vmec_free_boundary_check_chi_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=chi -tol=6.0E-16)
add_test (NAME    vmec_free_boundary_check_chipf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=chipf -tol=2.0E-14)
add_test (NAME    vmec_free_boundary_check_currumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=currumnc -tol=1.2E-6)
add_test (NAME    vmec_free_boundary_check_currvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=currvmnc -tol=2.0E-7)
add_test (NAME    vmec_free_boundary_check_DCurr_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=DCurr -tol=9.0E-10)
add_test (NAME    vmec_free_boundary_check_DGeod_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=DGeod -tol=1.2E-9)
add_test (NAME    vmec_free_boundary_check_DMerc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=DMerc -tol=1.9E-9)
add_test (NAME    vmec_free_boundary_check_DShear_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=DShear -tol=4.0E-10)
add_test (NAME    vmec_free_boundary_check_DWell_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=DWell -tol=2.0E-10)
add_test (NAME    vmec_free_boundary_check_fsql_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=fsql -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_fsqr_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=fsqr -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_fsqt_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=fsqt -tol=3.0E-17)
add_test (NAME    vmec_free_boundary_check_fsqz_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=fsqz -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_gmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=gmnc -tol=3.0E-14)
add_test (NAME    vmec_free_boundary_check_IonLarmor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=IonLarmor -tol=8.0E-17)
add_test (NAME    vmec_free_boundary_check_iotaf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=iotaf -tol=3.0E-13)
add_test (NAME    vmec_free_boundary_check_jcuru_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=jcuru -tol=3.1E-7)
add_test (NAME    vmec_free_boundary_check_jcurv_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=jcurv -tol=1.2E-9)
add_test (NAME    vmec_free_boundary_check_jdotb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=jdotb -tol=7.0E-7)
add_test (NAME    vmec_free_boundary_check_lmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=lmns -tol=9.0E-13)
add_test (NAME    vmec_free_boundary_check_mnmax_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=mnmax -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_mnmax_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=mnmax_nyq -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_phipf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=phipf -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_phips_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=phips -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_pres_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=pres -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_presf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=presf -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_q_factor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=q_factor -tol=2.0E-13)
add_test (NAME    vmec_free_boundary_check_rbtor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=rbtor -tol=2.0E-14)
add_test (NAME    vmec_free_boundary_check_rbtor0_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=rbtor0 -tol=4.8E-15)
add_test (NAME    vmec_free_boundary_check_rmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=rmnc -tol=5.0E-14)
add_test (NAME    vmec_free_boundary_check_specw_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=specw -tol=2.0E-12)
add_test (NAME    vmec_free_boundary_check_volavgB_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=volavgB -tol=8.0E-15)
add_test (NAME    vmec_free_boundary_check_volume_p_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=volume_p -tol=9.0E-15)
add_test (NAME    vmec_free_boundary_check_vp_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=vp -tol=8.3E-16)
add_test (NAME    vmec_free_boundary_check_wb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=wb -tol=2.9E-18)
add_test (NAME    vmec_free_boundary_check_wdot_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=wdot -tol=3.0E-15)
add_test (NAME    vmec_free_boundary_check_wp_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=wp -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_xm_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=xm -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_xm_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=xm_nyq -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_xn_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=xn -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_xn_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=xn_nyq -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_zmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_serial.vmec.nc -wout_file2=wout_test_parallel.vmec.nc -quantity=zmns -tol=5.0E-14)
