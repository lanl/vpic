# Tracer particles

## Background

In VPIC, as in real life different particles of one species are
indistinguishable. But being able to follow a particular particle over time can
be very helpful in physics investigation. Making particles artificially unique
and following them over time, often coupled with high cadence output of that
subset of particles, is commonly known as _tracer particles_. Users have
implemented this feature (multiple times) in their input decks, but VPIC now
has build-in support for this feature

## Create tracers

Usually only a limited number of particles are tracer particles since we can
not examine a trillion of particles anyway. And the extra information on tracer
particles introduces overhead that needs to be considered. So VPIC offers to
create a limited number of tracer particles in a new species to get insights
into the behavior of the parent species. These tracer particles will have
unique IDs to allow re-identification at later time steps. This implies that
VPIC should be compiled with _-DGLOBAL_PARTICLE_ID=ON_. The input deck will
then create a regular species. In the case of electrons that might look like

    species_t * electron = define_species( "electron", -ec, me, 1.5*Ne/nproc(), nmovers, sort_interval, sort_type);

After that all electron will be created by (many) calls to

    inject_particle( electron, x, y, z, ux, uy, uz, we, 0, 0 );

Once all electrons exist we can promote a small number of them to tracer
particles by calling

    species_t * elec_tracer  = make_tracers_by_percentage(electron, 1.0, Tracertype::Copy, "electron_tracer");

This would promote 1% of all electron to tracers. In the newly created
*elec_tracer* species they will have uniques IDs and can be written do disk
independently from all the other electrons. There is a number of these creation
functions, the name in the last argument is optional (without VPIC will
automatically chose a unique but boring name) and there is two options for the
second to last parameter: it can be Tracertype::Copy or Tracertype::Move. If
the user selects _Copy_ the selected fraction of particles will be copied to
the newly created species. Since they are duplicates they will have zero
statistical weight. The downside is that now both the original particle in the
electron species and the tracer in elec_tracer need memory and have to be
pushed separately. The upside is that the electron species unmodified.
Alternatively the user can select to _Move_ the particle in which case the
weight is retained and only one copy of the particle exists. The downside is
that the particles in *elec_tracer* do not automatically contribute to e.g. the
hydro output of all electrons and need to be added back in based on the hydro
output of the elec_tracer species. For a more complete list of these creation
functions and a documentation of their interface have a look at
src/species_advance/species_advance.cc and search for function that start with
tracerspecies_*.

## Annotate tracers

If the code is compiled with _-DVPIC_PARTICLE_ANNOTATION_ set to float or
double, then the use can requests extra buffer space per particle for custom
annotations. A call to

    elec_tracer->allocate_annotation_buffer(17)

will request enough space for 17 float variables per particle that are
initalized to 0. Calls to elec_tracer->get_annotation(i, j) will get the
annotation of particle i (between 0 and elec_tracer->np) in slot j (beween 0
and 17). Calls to elec_tracer->set_annotation(i, j, v) will set that variable
to the value in v. elec_tracer->increment_annotation(i, j, a) will increment
the value by a.

There is also the convenience functions

  void interpolate_fields_annotation(const char * sp_name,
    const interpolator_array_t * ia, int where_ex, int where_ey, int where_ez,
    int where_bx, int where_by, int where_bz);

to interpolate the field values to the particle locations and

  void interpolate_hydro_annotation(const char * sp_name, const hydro_array_t * ha,
                                    int where_jx,  int where_jy,  int where_jz,  int where_rho,
                                    int where_px,  int where_py,  int where_pz,  int where_ke,
                                    int where_txx, int where_tyy, int where_tzz,
                                    int where_tyz, int where_tzx, int where_txy);

to interpolate the hydro quantities to the particle locations and set the
annotations listed in the where_ variabels to those values. Some (but not all)
where_ locations can be set to -1 in which case the associated interpolated
value is not stored in any annotation slot.

## Output tracers

To output tracers (along with their global ID and annotation values) there is a
set of function for buffered output. The output bufferes are meant to be
allocated once, filled frequently over a number of timesteps and then dumped
once and cleared for further use. This limits the number of output operations
and improves the performance since the number of tracer particles is typically
small and IO performance would be dominated by metadata operations. Output is
performed to a single global HDF5 file per buffer dump. The functions to
interact with the output code from the input deck are

  void init_buffered_particle_dump(const char * speciesname, const int N_timesteps, const double safety_factor = 2.0);

to allocate a buffer for the named species with enough memory to handle
buffering for N_timesteps with up to safety_factor times as many particles of
that species as currently allocated localy.

  void accumulate_buffered_particle_dump(const char * speciesname, const int frame);

to copy the current state of the particles of the names species (including
annotation) into the slot enumerated by frame. Note that 0 <= frame <
N_timesteps has to be satisfied. This also implies that frames should be number
consequtively even if buffering is only performed e.g. every second step.

  void write_buffered_particle_dump(const char * speciesname);

will before the actual output to a single, global HDF5 file. (Further calls to
this function in later timesteps will use a new global file).

  void clear_buffered_particle_dump(const char * speciesname);

Clears the buffer for reuse in subsequent timesteps. If all frames are
overwritten before the next output it is not strictly necessary, but the
function is cheap (compared to the actual output) and can prevent confusion if
not all frames are overwritten.
