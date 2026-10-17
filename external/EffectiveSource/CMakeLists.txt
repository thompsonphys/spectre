set(LIBRARY EffectiveSource)

add_library(${LIBRARY} 
  kerr-circular.c 
  kerr-equatorial-ctx.c 
  kerr-equatorial-coeffs-ctx.c 
  kerr-equatorial-dtcoeffs-ctx.c 
  kerr-equatorial-dttcoeffs-ctx.c
  kerr-equatorial-mmode.c
  kerr-equatorial-numerical.c
  kerr-equatorial-point-ctx.c
)

target_link_libraries(
  ${LIBRARY}
  PRIVATE
  GSL::gsl
  )

# The m-modes evaluate the coefficients of the field point and the analytic
# contraction or numerical azimuthal sum at every (grid point, time) sample,
# and the numerical reference in the unit test evaluates the effective source
# at many azimuths, so these files dominate the run time even in Debug builds.
# They compile quickly with optimizations.
# The coefficient files need a lot of memory to optimize and only run once per
# time sample, so they keep the default flags.
set_source_files_properties(
  kerr-equatorial-ctx.c
  kerr-equatorial-mmode.c
  kerr-equatorial-numerical.c
  kerr-equatorial-point-ctx.c
  PROPERTIES COMPILE_OPTIONS "$<$<CONFIG:Debug>:-O2>"
  )

target_include_directories(${LIBRARY} SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
