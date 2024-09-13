#  Run the reconstruction task. This just tests if the reconstruction executes correctly.
#  Subsequent tests will examine the results of the reconstruction for correct
#  results.

#  Copy input file from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.test.vmec
                ${CMAKE_CURRENT_BINARY_DIR}/input.test.vmec
                COPYONLY)

add_test (NAME niter_clean_wout_test
          COMMAND rm -f wout_test.vmec.nc)

add_test (NAME    niter_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test.vmec)
set_tests_properties (niter_test
                      PROPERTIES
                      DEPENDS niter_clean_wout_test
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)

add_test (NAME niter_check_wout_test
          COMMAND ${Python3_EXECUTABLE} ${TEST_SCRIPT_DIR}/check_file.py --file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc)
set_tests_properties (niter_check_wout_test
                      PROPERTIES
                      DEPENDS niter_test)
