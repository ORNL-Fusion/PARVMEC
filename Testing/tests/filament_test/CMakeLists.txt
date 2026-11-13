#  Copy input files from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/coils.test
                ${CMAKE_CURRENT_BINARY_DIR}/coils.test
                COPYONLY)

#  Vacuum field of an axisymmetric coil set on twelve toroidal planes.
add_test (NAME    vmec_filament_mgrid_test
          COMMAND $<TARGET_PROPERTY:mgrid,BINARY_DIR>/mgrid coils.test)

#  A free boundary tokamak on the axisymmetric grid and on the twelve planes.
foreach (name test_axisymmetric test_planes)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/input.${name}.vmec
                    ${CMAKE_CURRENT_BINARY_DIR}/input.${name}.vmec
                    COPYONLY)
    add_test (NAME    vmec_filament_${name}
              COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.${name}.vmec)
#  Runs in one directory share parvmecinfo.txt and timings.txt.
    set_tests_properties (vmec_filament_${name}
                          PROPERTIES
                          DEPENDS vmec_filament_mgrid_test
                          RESOURCE_LOCK filament_test)
endforeach ()

#  The plasma and the coils are axisymmetric, so the two grids give the same
#  equilibrium up to the toroidal resolution of the vacuum solution.
foreach (check rmin_surf:5.0E-3 rmax_surf:5.0E-3 iotaf:2.0E-3)
    string (REPLACE ":" ";" check ${check})
    list (GET check 0 quantity)
    list (GET check 1 tolerance)
    add_test (NAME    vmec_filament_check_${quantity}_test
              COMMAND $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xwout_diff -wout_file1=wout_test_axisymmetric.vmec.nc -wout_file2=wout_test_planes.vmec.nc -quantity=${quantity} -tol=${tolerance})
    set_tests_properties (vmec_filament_check_${quantity}_test
                          PROPERTIES
                          DEPENDS "vmec_filament_test_axisymmetric;vmec_filament_test_planes")
endforeach ()
