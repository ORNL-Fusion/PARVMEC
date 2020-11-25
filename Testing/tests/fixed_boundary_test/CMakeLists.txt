#  Copy input file from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test_serial.vmec
                COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test_parallel.vmec
                COPYONLY)

#  Test Fixed boundary equilibrium.
add_test (NAME    vmec_fixed_boundary_serial_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test_serial.vmec)
add_test (NAME    vmec_fixed_boundary_parallel_test
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPI_OVERSUBSCRIBE_FLAG} $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test_parallel.vmec)

#  Check woutfile to ensure that values match.
add_test (NAME    vmec_fixed_boundary_check_aspect_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=aspect -tol=1.0E-20)
add_test (NAME    vmec_fixed_boundary_check_b0_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=b0 -tol=1.0E-20)
add_test (NAME    vmec_fixed_boundary_check_bdotb_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=bdotb -tol=1.0E-20)
add_test (NAME    vmec_fixed_boundary_check_bdotgradv_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=bdotgradv -tol=1.0E-20)
add_test (NAME    vmec_fixed_boundary_check_beta_vol_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=beta_vol -tol=1.0E-20)
add_test (NAME    vmec_fixed_boundary_check_betapol_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=betapol -tol=1.0E-20)
add_test (NAME    vmec_fixed_boundary_check_betator_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=betator -tol=1.0E-20)
add_test (NAME    vmec_fixed_boundary_check_betatotal_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=betatotal -tol=1.0E-20)
add_test (NAME    vmec_fixed_boundary_check_betaxis_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=betaxis -tol=1.0E-20)
add_test (NAME    vmec_fixed_boundary_check_bmnc_test
          COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=input.test_serial.vmec -wout_file2=input.test_parallel.vmec -quantity=bmnc -tol=1.0E-20)
