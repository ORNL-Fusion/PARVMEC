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
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPI_OVERSUBSCRIBE_FLAG} $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test_parallel.vmec)

#  Check woutfile to ensure that values match.
add_test (NAME    vmec_free_boundary_check_aspect_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=aspect -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_b0_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=b0 -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_bdotb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=bdotb -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_bdotgradv_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=bdotgradv -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_beta_vol_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=beta_vol -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_betapol_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=betapol -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_betator_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=betator -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_betatotal_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=betatotal -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_betaxis_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=betaxis -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_bmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=bmnc -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_bsubsmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=bsubsmns -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_bsubumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=bsubumnc -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_bsubvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=bsubvmnc -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_bsupumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=bsupumnc -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_bsupvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=bsupvmnc -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_buco_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=buco -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_bvco_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=bvco -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_chi_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=chi -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_chipf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=chipf -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_currumnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=currumnc -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_currvmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=currvmnc -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_DCurr_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=DCurr -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_DGeod_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=DGeod -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_DMerc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=DMerc -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_DShear_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=DShear -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_DWell_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=DWell -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_fsql_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=fsql -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_fsqr_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=fsqr -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_fsqt_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=fsqt -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_fsqz_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=fsqz -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_gmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=gmnc -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_IonLarmor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=IonLarmor -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_iotaf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=iotaf -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_jcuru_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=jcuru -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_jcurv_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=jcurv -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_jdotb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=jdotb -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_lmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=lmns -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_mnmax_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=mnmax -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_mnmax_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=mnmax_nyq -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_phipf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=phipf -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_phips_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=phips -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_pres_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=pres -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_presf_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=presf -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_q_factor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=q_factor -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_rbtor_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=rbtor -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_rbtor0_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=rbtor0 -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_rmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=rmnc -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_specw_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=specw -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_volabgB_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=volabgB -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_volume_p_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=volume_p -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_vp_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=vp -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_wb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=wb -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_wdot_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=wdot -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_wp_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=wp -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_xm_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=xm -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_xm_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=xm_nyq -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_xn_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=xn -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_xn_nyq_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=xn_nyq -tol=1.0E-20)
add_test (NAME    vmec_free_boundary_check_zmns_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=zmns -tol=1.0E-20)
