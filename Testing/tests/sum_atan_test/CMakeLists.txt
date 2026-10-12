#  Copy input files from the source to the build directory.
foreach (name test_positive test_negative test_positive_q0 test_negative_q0
               test_iota_positive test_iota_negative
               test_iota_positive_q0 test_iota_negative_q0)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.${name}.vmec
                    ${CMAKE_CURRENT_BINARY_DIR}/input.${name}.vmec
                    COPYONLY)
    add_test (NAME    vmec_sum_atan_${name}
              COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.${name}.vmec)
#  Runs in one directory share parvmecinfo.txt and timings.txt.
    set_tests_properties (vmec_sum_atan_${name}
                          PROPERTIES
                          RESOURCE_LOCK sum_atan_test)
endforeach ()

#  The same sum_atan current or iota profile written with the signs of the
#  amplitude and the coefficient both reversed gives the same equilibrium.
foreach (pair test_positive:test_negative test_positive_q0:test_negative_q0
              test_iota_positive:test_iota_negative
              test_iota_positive_q0:test_iota_negative_q0)
    string (REPLACE ":" ";" pair ${pair})
    list (GET pair 0 first)
    list (GET pair 1 second)
    foreach (quantity ctor jcurv iotaf)
        add_test (NAME    vmec_sum_atan_${second}_check_${quantity}_test
                  COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_${first}.vmec.nc -wout_file2=wout_${second}.vmec.nc -quantity=${quantity} -tol=1.0E-20)
        set_tests_properties (vmec_sum_atan_${second}_check_${quantity}_test
                              PROPERTIES
                              DEPENDS "vmec_sum_atan_${first};vmec_sum_atan_${second}")
    endforeach ()
endforeach ()
