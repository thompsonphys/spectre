\cond NEVER
Distributed under the MIT License.
See LICENSE.txt for details.
\endcond
# %Running CCE {#tutorial_cce}

\tableofcontents

The basic instructions for getting up and running with a stand-alone
CCE using external data are:
- Clone spectre and build the CharacteristicExtract target
- At this point (provided the build succeeds) you should have the
  executable `bin/CharacteristicExtract` in the build directory; You can now
  run it using an input file to provide options. The input file
  `tests/InputFiles/Cce/CharacteristicExtract.yaml` from the spectre
  source tree can help you get started with writing your input file.
  There are a few important notes there:
  - For resolution, the example input file has lmax (`Cce.LMax`) of 20, and
    filter lmax (`Filtering.FilterLMax`) of 18; that may run a bit slow for
    basic tests, but this value yields the best precision-to-run-time ratio
    for a typical BBH system. Note that precision doesn't improve above lmax 24,
    filter 22 (be sure to update the filter as you update lmax -- it should
    generally be around 2 lower than the maximum l to reduce possible aliasing).
  - If you want to just run through the end of the provided worldtube data,
    you can just omit the `EndTime` option and the executable will figure it
    out from the worldtube file.
  - The `ScriOutputDensity` adds extra interpolation points to the output,
    which is useful for finite-difference derivatives on the output data, but
    otherwise it'll just unnecessarily inflate the output files, so if you
    don't need the extra points, best just set it to 1.
  - If you're extracting at 100M or less (which isn't currently recommended due
    to the junk radiation being much worse), best to keep the `TargetStepSize`
    as .5 for 100M and perhaps even lower yet for nearer extraction.
  - The `InitializeJ` in the example file uses `ConformalFactor` which has been
    found to perform better than the other schemes implemented so far.
    Other schemes for `InitializeJ` can be found in namespace
    \link Cce::InitializeJ \endlink.  Some of these are:
    - `ZeroNonSmooth`: Make `J` vanish
    - `NoIncomingRadiation`: Make \f$\Psi_0 = 0\f$; this does not actually lead
      to no incoming radiation, since \f$\Psi_0\f$ and \f$\Psi_4\f$ both include
      incoming and outgoing radiation.
    - `InverseCubic`: Ansatz where \f$J = A/r + B/r^3\f$
    - `ConformalFactor`: Try to make initial time coordinate as inertial as
      possible at \f$\mathscr{I}^+\f$.

- An example of an appropriate submission command for slurm systems is:
  ```
  srun -n 1 -c 1 path/to/build/bin/CharacteristicExtract ++ppn 3 \
  --input-file path/to/input.yaml
  ```
  CCE doesn't currently scale to more than 4 cores, so those slurm options are
  best.
- CCE will work faster if the input worldtube hdf5 file is chunked in small
  numbers of complete rows.
  This is relevant because by default, SpEC writes its worldtube files
  chunked along full time-series columns, which is efficient for
  compression, but not for reading in to SpECTRE -- in that case,
  it is recommended to rechunk the input file before running CCE
  for maximum performance. This can be done, for instance, using h5py
  (you will need to fill in filenames appropriate to your case in place
  of "BondiCceR0050.h5" and "RechunkBondiCceR0050.h5"):
  ```py
  import h5py
  input_file = "BondiCceR0050.h5"
  output_file = "RechunkBondiCceR0050.h5"
  with h5py.File(input_file,'r') as input_h5,\
    h5py.File(output_file, 'w') as output_h5:
    for dset in input_h5:
        if("Version" in dset):
            output_h5[dset] = input_h5[dset][()]
            continue
        number_of_columns = input_h5[dset][()].shape[1]
        output_h5.create_dataset(dset, data=input_h5[dset],
                                 maxshape=(None, number_of_columns),
                                 chunks=(4, number_of_columns), dtype='d')
        for attribute in input_h5[dset].attrs.keys():
            output_h5[dset].attrs[attribute] = input_h5[dset].attrs[attribute]
  ```
- Note, these `*.h5` files required as an `input_file` (and also needed
  in the `.yaml` file) will need to be obtained from an SXS member who has
  completed a relevant `SpEC` run.
- The output data will be written as spin-weighted spherical harmonic
  modes, one physical quantity per dataset, and each row will
  have the time value followed by the real and imaginary parts
  of the complex modes in m-varies-fastest order.

Once you have the reduction data output file from a successful CCE run, you can
confirm the integrity of the h5 file and its contents by running
```
h5ls -r CharacteristicExtractReduction.h5
```

For the reduction file produced by a successful run, the output of the `h5ls`
should resemble

```
/SpectreR0100.cce                 Group
/SpectreR0100.cce/EthInertialRetardedTime Dataset {26451/Inf, 163}
/SpectreR0100.cce/News            Dataset {26451/Inf, 163}
/SpectreR0100.cce/Psi0            Dataset {26451/Inf, 163}
/SpectreR0100.cce/Psi1            Dataset {26451/Inf, 163}
/SpectreR0100.cce/Psi2            Dataset {26451/Inf, 163}
/SpectreR0100.cce/Psi3            Dataset {26451/Inf, 163}
/SpectreR0100.cce/Psi4            Dataset {26451/Inf, 163}
/SpectreR0100.cce/Strain          Dataset {26451/Inf, 163}
/src.tar.gz              Dataset {7757329}
```

\note Prior to
[this Pull Request](https://github.com/sxs-collaboration/spectre/pull/5985), the
output of `h5ls` looked like this
```
/                        Group
/Cce                     Group
/Cce/EthInertialRetardedTime.dat Dataset {3995/Inf, 163}
/Cce/News.dat                 Dataset {3995/Inf, 163}
/Cce/Psi0.dat                 Dataset {3995/Inf, 163}
/Cce/Psi1.dat                 Dataset {3995/Inf, 163}
/Cce/Psi2.dat                 Dataset {3995/Inf, 163}
/Cce/Psi3.dat                 Dataset {3995/Inf, 163}
/Cce/Psi4.dat                 Dataset {3995/Inf, 163}
/Cce/Strain.dat               Dataset {3995/Inf, 163}
/src.tar.gz               Dataset {3750199}
```

The `Strain` represents the asymptotic transverse-traceless contribution
to the metric scaled by the Bondi radius (to give the asymptotically leading
part), the `News` represents the first derivative of the strain, and each
of the `Psi...` datasets represent the Weyl scalars, each scaled by the
appropriate factor of the Bondi-Sachs radius to retrieve the asymptotically
leading contribution.

The `EthInertialRetardedTime` is a diagnostic dataset that represents the
angular derivative of the inertial retarded time, which determines the
coordinate transformation that is performed at scri+.

If you'd like to visualize the output of a CCE run, we offer a
[CLI](py/cli.html) that will produce a plot of all of quantities except
`EthInertialRetardedTime`. To see how to use this cli, run

```
spectre plot cce -h
```

If you'd like to do something more complicated than just make a quick plot,
you'll have to load in the output data yourself using `h5py` or our
`spectre.IO.H5` bindings.

\note The CLI can also plot the "old" version of CCE output, described above.
Pass `--cce-group Cce` to the CLI. This option is only for backwards
compatibility with the old CCE output and is not supported for the current
version of output. This options is deprecated and will be removed in the future.

### Input data formats

The worldtube data must be constructed as spheres of constant coordinate
radius, and (for the time being) written to a filename of the format
`...CceRXXXX.h5`, where the `XXXX` is to be replaced by the integer for which
the extraction radius is equal to `XXXX`M. For instance, a 100M extraction
should have filename `...CceR0100.h5`. This scheme of labeling files with the
extraction radius is constructed for compatibility with SpEC worldtube data.

There are two possible formats of the input data, one based on the Cauchy metric
at finite radius, and one based on Bondi data. The metric data format must be
provided as spherical harmonic modes with the following datasets:
- `gxx.dat`, `gxy.dat`, `gxz.dat`, `gyy.dat`, `gyz.dat`, `gzz.dat`
- `Drgxx.dat`, `Drgxy.dat`, `Drgxz.dat`, `Drgyy.dat`, `Drgyz.dat`, `Drgzz.dat`
- `Dtgxx.dat`, `Dtgxy.dat`, `Dtgxz.dat`, `Dtgyy.dat`, `Dtgyz.dat`, `Dtgzz.dat`
- `Shiftx.dat`, `Shifty.dat`, `Shiftz.dat`
- `DrShiftx.dat`, `DrShifty.dat`, `DrShiftz.dat`
- `DtShiftx.dat`, `DtShifty.dat`, `DtShiftz.dat`
- `Lapse.dat`
- `DrLapse.dat`
- `DtLapse.dat`

In this format, each row must start with the time stamp, and the remaining
values are the complex modes in m-varies-fastest format. That is,
```
"time", "Lapse_Re(0,0)", "Lapse_Im(0,0)",
"Lapse_Re(1,1)", "Lapse_Im(1,1)", "Lapse_Re(1,0)", "Lapse_Im(1,0)",
"Lapse_Re(1,-1)", "Lapse_Im(1,-1)",
"Lapse_Re(2,2)", "Lapse_Im(2,2)", "Lapse_Re(2,1)", "Lapse_Im(2,1)",
"Lapse_Re(2,0)", "Lapse_Im(2,0)", "Lapse_Re(2,-1)", "Lapse_Im(2,-1)",
"Lapse_Re(2,-2)", "Lapse_Im(2,-2)"
```
Each dataset in the file must also have an attribute named "Legend" which
is an ASCII-encoded null-terminated variable-length string. That is, the HDF5
type is:
```
DATATYPE  H5T_STRING {
  STRSIZE H5T_VARIABLE;
  STRPAD H5T_STR_NULLTERM;
  CSET H5T_CSET_ASCII;
  CTYPE H5T_C_S1;
}
```
This can be checked for a dataset by running
```
h5dump -a DrLapse.dat/Legend CceR0150.h5
```

The second format is Bondi-Sachs metric component data.
This format is far more space-efficient (by around a factor of 4), and SpECTRE
provides a separate executable for converting to the Bondi-Sachs worldtube
format, `ReduceCceWorldtube`.
The `ReduceCceWorldtube` executable should be run on a Cauchy metric worldtube
file, and will produce a corresponding 'reduced' Bondi-Sachs worldtube file.
The basic command-line arguments for the executable are:
```
ReduceCceWorldtube --input-file CceR0050.h5 --output-file BondiCceR0050.h5\
 --lmax_factor 3
```
The argument `--lmax_factor` determines the factor by which the resolution at
which the boundary computation that is run will exceed the resolution of the
input and output files.
Empirically, we have found that `lmax_factor` of 3 is sufficient to achieve
roundoff precision in all boundary data we have attempted, and an `lmax_factor`
of 2 is usually sufficient to vastly exceed the precision of the simulation that
provided the boundary dataset.

The format is similar to the metric components, except in spin-weighted
spherical harmonic modes, and the real (spin-weight-0) quantities omit the
redundant negative-m modes and imaginary parts of m=0 modes.
The quantities that must be provided by the Bondi-Sachs metric data format are:
- `Beta.dat`
- `DrJ.dat`
- `DuR.dat`
- `H.dat`
- `J.dat`
- `Q.dat`
- `R.dat`
- `U.dat`
- `W.dat`

The Bondi-Sachs data may also be used directly for CCE input data.
To specify that the input type is in 'reduced' Bondi-Sachs form, use:
```
...
Cce:
  H5IsBondiData: True
...
```
Otherwise, for the Cauchy metric data format, use:
```
...
Cce:
  H5IsBondiData: False
...
```
