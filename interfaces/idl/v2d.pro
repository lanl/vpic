   ;==========================================================================
   ; Particle postprocessor.  Works off of particle data file prepared using
   ; pppp.c particle pre-postprocessor.  
   ;
   pro v2d, vdata
     fname = dialog_pickfile(/read, filter = 'particle*')
     openr, unit1, fname, /get_lun, /swap_if_big_endian
     phasedata = fltarr(256,256)
     vdata=assoc(unit1,phasedata)
   end
   ;=========================================================================

