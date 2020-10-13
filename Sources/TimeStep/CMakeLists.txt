target_sources(vmec
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/eqsolve.f
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/evolve.f
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gmres_mod.f
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interp.f
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/parallel_include_module.f90
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/parallel_vmec_module.f90
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/restart.f
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/runvmec.f
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vmec_history.f
)

target_sources(xvmec
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vmec.f
)
