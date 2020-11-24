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
