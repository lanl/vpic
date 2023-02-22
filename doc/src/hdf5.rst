======================
HDF5 Dump Code Usage
======================
1. How to get the HDF5 dump code for VPIC?
   Before the merge to https://github.com/lanl/vpic: 
   > git clone https://github.com/brtnfld/vpic.git
   > git check dump_strategy

   Note: please check status of merge at the Pull request here: 
   https://github.com/lanl/vpic/pull/144
   After the merge, you may get it from https://github.com/lanl/vpic

2. How to compile the HDF5 dump code for VPIC?
   General build step can be found here: https://github.com/lanl/vpic
   Here we use the reference-Debug and as example 

   > mkdir build
   > ../arch/reference-Debug
   > make
   > cp ../sample/harrisHDF5 ./ 
   > ./bin/vpic harrisHDF5

  Trouble shooting for test cases:

  Error: Undefined symbols for architecture x86_64: "_H5Dclose" 
  Solution: you may manually add below flag for HDF5 to bin/vpic file
            this is because the cmake sometimes can not find proper HDF5 file and link libraries 

  -I/Users/dbin/work/soft/hdf5-1.10.7/build/include -L/Users/dbin/work/soft/hdf5-1.10.7/build/lib -lhdf5  -lhdf5_hl -lz


3. How to run VPIC with the HDF5 dump code ?
   > ./harrisHDF5.Darwin 
      *** Initializing
      Booting with 1 threads (pipelines) and 1 (MPI) ranks
      ./harrisHDF5(64)[0]: Defaulting to mass_ratio of 1 and seed of 0
      ./harrisHDF5(65)[0]: For Custom Usage: ./harrisHDF5.Darwin mass_ratio seed
      ./harrisHDF5(75)[0]: Computing simulation parameters
      
      ... ...

      user_particle_injection |   0% 0.0e+00 4.0e+00 0.0e+00 |   0% 4.3e-05 4.8e+02 9.0e-08
      user_current_injection |   0% 0.0e+00 4.0e+00 0.0e+00 |   0% 2.7e-05 4.8e+02 5.7e-08
      user_field_injection |   0% 0.0e+00 4.0e+00 0.0e+00 |   0% 2.0e-05 4.8e+02 4.1e-08
      user_diagnostics |  46% 1.4e-01 5.0e+00 2.7e-02 |  15% 2.9e+00 4.8e+02 5.9e-03

      normal exit

4. How to inspect the results?

Particle data is stored in "particle_hdf5" directory.
Within the "particle_hdf5" directory, it has time steps, from T.0 , T.1, ....
Within each time step, it has two files electron_[time step].h5 and  ion_[time step].h5
> tree particle_hdf5
particle_hdf5
├── T.0
│   ├── electron_0.h5
│   └── ion_0.h5
      ... ....

├── T.480
    ├── electron_480.h5
    └── ion_480.h5

For each particle file, you can use h5dump to inspect its results:

> h5dump -A particle_hdf5/T.0/electron_0.h5
HDF5 "particle_hdf5/T.0/electron_0.h5" {
GROUP "/" {
   GROUP "Timestep_0" {
      DATASET "dX" {
         DATATYPE  H5T_IEEE_F32LE
         DATASPACE  SIMPLE { ( 131072 ) / ( 131072 ) }
      }
      
      ... ...

      DATASET "uz" {
         DATATYPE  H5T_IEEE_F32LE
         DATASPACE  SIMPLE { ( 131072 ) / ( 131072 ) }
}}}}


Hydro data is stored in "hydro_hdf5" directory.
It has the same structure as the "particle_hdf5" directory
By default, it also has two extra files "hydro-electron.xdmf" and "hydro-ion.xdmf"
The "hydro-electron.xdmf"  and "hydro-ion.xdmf" are used to plot results of hydro data
One can also use h5dump to inspect the results of each file.

> tree hydro_hdf5
hydro_hdf5
├── T.0
│   ├── hydro_electron_0.h5
│   └── hydro_ion_0.h5
   
   ... ...

├── T.480
│   ├── hydro_electron_480.h5
│   └── hydro_ion_480.h5
├── hydro-electron.xdmf
└── hydro-ion.xdmf


Hydro data is stored in "field_hdf5" directory.
It has the same structure as the "hydro_hdf5" directory
By default, it also has two extra files "hdf5_field.xdmf".
The "hdf5_field.xdmf" is used to plot results of field data
One can also use h5dump to inspect the results of each file.


> tree field_hdf5
field_hdf5
├── T.0
│   └── fields_0.h5

  ... ...

├── T.480
│   └── fields_480.h5
└── hdf5_field.xdmf


5. Other Notes on the HDF5 dump code (ongoing work)
   * How to choose different output structure?
   * How to select different attributes of fields for dump?
   
   The src/vpic/dump_strategy.h comes with different option to 
   explore different output structures.

    By default, below lines are commented out. 
    You can remove the comment symbols to use compound data structure for field/hydro/particle data.
    But you may need to re-compile the vpic and the deck file 
   // #define HAS_FIELD_COMP 1
   // #define HAS_HYDRO_COMP 1
   // #define HAS_PARTICLE_COMP 1

   One can also try the async I/O in HDF5 to dump data  (not tested recently)
   //#define H5_ASYNC 1
   #ifdef H5_ASYNC
   #include "h5_vol_external_async_native.h"
   #endif

========================================
HDF5 Dump Code's XDMF file for ParaView
========================================
1. XDMF file 
   One can use the ParaView to open "hdf5_field.xdmf" to view the filed data. 
   One can use the ParaView to open "hydro-electron.xdmf"  and "hydro-ion.xdmf" to view the hydro data. 

==========================
HDF5 Dump Code Test Cases
==========================
1. Test cases covered by the HDF5 Dump Code

The HDF5 have below test case for HDF5 dump code
//Filed related test cases
- hdf5_fields
- hdf5_fields_restart
- hdf5_fields_readhdf5
- hdf5_fields_location
- hdf5_fields_location_restart
- hdf5_fields_location_readhdf5
- hdf5_fields_many_steps
- hdf5_fields_many_steps_restart
- hdf5_fields_many_steps_readhdf5

//Hydro related test cases
- hdf5_hydro
- hdf5_hydro_readhdf5

//Particle related test cases
- buffered_particle_dump_readhdf5

You can run all these test cases with below command
> make test

Running tests...
Test project /Users/dbin/work/soft/vpic-doc/vpic/build
      Start  1: rng
 1/37 Test  #1: rng ...............................   Passed    0.03 sec
      Start  2: array_index
 2/37 Test  #2: array_index .......................   Passed    0.07 sec
      Start  3: buffered_particle_dump
 3/37 Test  #3: buffered_particle_dump ............   Passed    0.09 sec
      Start  4: buffered_particle_dump_restart
 4/37 Test  #4: buffered_particle_dump_restart ....   Passed    0.09 sec
      Start  5: buffered_particle_dump_readhdf5
 5/37 Test  #5: buffered_particle_dump_readhdf5 ...   Passed    3.53 sec
      Start  6: hdf5_fields
 6/37 Test  #6: hdf5_fields .......................   Passed    0.11 sec
      Start  7: hdf5_fields_restart
 7/37 Test  #7: hdf5_fields_restart ...............   Passed    0.09 sec
      Start  8: hdf5_fields_readhdf5
 8/37 Test  #8: hdf5_fields_readhdf5 ..............   Passed    0.16 sec
      Start  9: hdf5_fields_location
 9/37 Test  #9: hdf5_fields_location ..............   Passed    0.11 sec
      Start 10: hdf5_fields_location_restart
10/37 Test #10: hdf5_fields_location_restart ......   Passed    0.09 sec
      Start 11: hdf5_fields_location_readhdf5
11/37 Test #11: hdf5_fields_location_readhdf5 .....   Passed    0.16 sec
      Start 12: hdf5_fields_many_steps
12/37 Test #12: hdf5_fields_many_steps ............   Passed    0.12 sec
      Start 13: hdf5_fields_many_steps_restart
13/37 Test #13: hdf5_fields_many_steps_restart ....   Passed    0.11 sec
      Start 14: hdf5_fields_many_steps_readhdf5
14/37 Test #14: hdf5_fields_many_steps_readhdf5 ...   Passed    0.18 sec
      Start 15: hdf5_hydro
15/37 Test #15: hdf5_hydro ........................   Passed    0.09 sec
      Start 16: hdf5_hydro_readhdf5
16/37 Test #16: hdf5_hydro_readhdf5 ...............   Passed    0.16 sec
      Start 17: accel
17/37 Test #17: accel .............................   Passed    0.06 sec
      Start 18: cyclo
18/37 Test #18: cyclo .............................   Passed    0.06 sec
      Start 19: inbndj
19/37 Test #19: inbndj ............................   Passed    0.07 sec
      Start 20: interpe
20/37 Test #20: interpe ...........................   Passed    0.06 sec
      Start 21: outbndj
21/37 Test #21: outbndj ...........................   Passed    1.28 sec
      Start 22: pcomm
22/37 Test #22: pcomm .............................   Passed    0.10 sec
      Start 23: simple
23/37 Test #23: simple ............................   Passed    0.21 sec
      Start 24: dump
24/37 Test #24: dump ..............................   Passed    0.07 sec
      Start 25: reconnection_test
25/37 Test #25: reconnection_test .................   Passed    0.69 sec
      Start 26: parallel
26/37 Test #26: parallel ..........................   Passed    0.22 sec
      Start 27: threaded
27/37 Test #27: threaded ..........................   Passed    0.28 sec
      Start 28: generate_restore
28/37 Test #28: generate_restore ..................   Passed    0.07 sec
      Start 29: perform_restore
29/37 Test #29: perform_restore ...................   Passed    0.06 sec
      Start 30: test_collision
30/37 Test #30: test_collision ....................   Passed    0.38 sec
      Start 31: test_collision_script
31/37 Test #31: test_collision_script .............   Passed    0.13 sec
      Start 32: array_syntax
32/37 Test #32: array_syntax ......................   Passed    0.07 sec
      Start 33: weibel_driver
33/37 Test #33: weibel_driver .....................   Passed    0.75 sec
      Start 34: 3d_test
34/37 Test #34: 3d_test ...........................   Passed    1.15 sec
      Start 35: weibel_driver_non_vaccum
35/37 Test #35: weibel_driver_non_vaccum ..........   Passed    0.74 sec
      Start 36: 3d_test_non_vaccum
36/37 Test #36: 3d_test_non_vaccum ................   Passed    1.11 sec
      Start 37: 3d_test_non_vaccum_threaded
37/37 Test #37: 3d_test_non_vaccum_threaded .......   Passed    0.29 sec

100% tests passed, 0 tests failed out of 37

Total Test time (real) =  13.07 sec



2. Trouble shooting for test cases
Error 1:  ModuleNotFoundError: No module named 'h5py'
Solution: pip install h5py

==========================
Todo
==========================



