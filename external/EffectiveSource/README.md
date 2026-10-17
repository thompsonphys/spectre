This is a modified version of https://github.com/barrywardell/EffectiveSource
(revision 07a31cedfcc895fbdaab47de88825234a03f7cec).
This code is licensed under MIT. The license can be found in LICENSE.md.

Summary of modifications:

- Added CMakeLists.txt.
- Added C++ compatibility header `effsource.hpp`.
- Added the numerical m-modes of the equatorial effective source and singular
  field: `kerr-equatorial-numerical.c`, `kerr-equatorial-point-ctx.c`,
  `effsource_equatorial_numerical.h`, `effsource_equatorial_point.h` and
  `effsource_dd.h`, copied from the `phi_resum` branch of
  https://github.com/thompsonphys/EffectiveSource (revision
  789d634a9dcb4747228499af09b668d50888999a). `kerr-equatorial-point-ctx.c` and
  `effsource_equatorial_point.h` are generated there by
  `codegen/point_codegen.py` from `kerr-equatorial-ctx.c`.
- Added the new headers to the C++ compatibility header
  `effsource_equatorial.hpp`.
- Added the analytic m-modes of the equatorial effective source and singular
  field: `kerr-equatorial-mmode.c` and `effsource_equatorial_mmode.h`, copied
  from the `updated_analytic_m` branch of
  https://github.com/thompsonphys/EffectiveSource (revision
  4d6f2a334bf33a8b93316b384a011ac78c788e07), and added
  `effsource_equatorial_mmode.h` to the C++ compatibility header
  `effsource_equatorial.hpp`.
