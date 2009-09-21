SECTIONS
{
  OVERLAY :
  {
    /* The stuff is .text_continued is really stuff we didn't want as part of
       the overlays.  So, we throw it the sole segment of its own region */
    .text_continued {
      :root_segment.o (.text.main)
      /opt/cell/toolchain/lib/gcc/spu/4.1.1/libgcc.a:_unpack_df.o (.text)
      /opt/cell/toolchain/lib/gcc/spu/4.1.1/libgcc.a:_df_to_si.o (.text)
      /opt/cell/toolchain/lib/gcc/spu/4.1.1/libgcc.a:_divdi3.o (.text)
      /opt/cell/toolchain/lib/gcc/spu/4.1.1/libgcc.a:_muldi3.o (.text)
      /opt/cell/toolchain/lib/gcc/spu/4.1.1/libgcc.a:fast_divdf3.o (.text)
      /opt/cell/sysroot/usr/spu/lib/libc.a:lib_a-memset.o (.text.memset)
    }
  }

  OVERLAY :
  {
    .ovly_compiler {
      /opt/cell/sysroot/usr/spu/lib/libgloss.a:exit.o (.text)
      /opt/cell/sysroot/usr/spu/lib/libc.a:lib_a-exit.o (.text.exit)
      /opt/cell/sysroot/usr/spu/lib/libc.a:lib_a-__call_atexit.o (.text.__call_exitprocs)
      /opt/cell/sysroot/usr/spu/lib/libc.a:lib_a-__atexit.o (.text.__register_exitproc)
      /opt/cell/sysroot/usr/spu/lib/libc.a:lib_a-atexit.o (.text.atexit)
      :/opt/cell/sysroot/usr/spu/lib/crti.o (.init)
      :/opt/cell/toolchain/lib/gcc/spu/4.1.1/crtbegin.o (.init)
      :/opt/cell/toolchain/lib/gcc/spu/4.1.1/crtend.o (.init)
      :/opt/cell/sysroot/usr/spu/lib/crtn.o (.init)
      :/opt/cell/toolchain/lib/gcc/spu/4.1.1/crtbegin.o (.text)
      :/opt/cell/toolchain/lib/gcc/spu/4.1.1/crtend.o (.text)
      :/opt/cell/sysroot/usr/spu/lib/crti.o (.fini)
      :/opt/cell/toolchain/lib/gcc/spu/4.1.1/crtbegin.o (.fini)
      :/opt/cell/sysroot/usr/spu/lib/crtn.o (.fini)
      /opt/cell/sysroot/usr/spu/lib/libc.a:lib_a-assert.o (.text.__assert_func)
      /opt/cell/sysroot/usr/spu/lib/libc.a:lib_a-stdio.o (.text.__sinit)
      /opt/cell/sysroot/usr/spu/lib/libc.a:lib_a-stdio.o (.text.__check_init)
      /opt/cell/sysroot/usr/spu/lib/libc.a:lib_a-fprintf.o (.text)
      /opt/cell/sysroot/usr/spu/lib/libgloss.a:syscalls.o (.text)
      /opt/cell/sysroot/usr/spu/lib/libc.a:lib_a-stack_reg_va.o (.text)
      /opt/cell/sysroot/usr/spu/lib/libc.a:lib_a-assert.o (.text.__assert)
      /opt/cell/sysroot/usr/spu/lib/libc.a:lib_a-stdio.o (.text.__sfp)
      /opt/cell/sysroot/usr/spu/lib/libc.a:lib_a-fclose.o (.text.fclose)
      /opt/cell/sysroot/usr/spu/lib/libc.a:lib_a-stdio.o (.text.__cleanup)
    }
    .ovly_particle_advance {
      :advance_p_pipeline_spu.o (.text._SPUEAR_advance_p_pipeline_spu)
    }
    .ovly_field_advance {
      :clear_accumulators_pipeline_spu.o (.text._SPUEAR_clear_accumulators_pipeline_spu)
      :reduce_accumulators_pipeline_spu.o (.text._SPUEAR_reduce_accumulators_pipeline_spu)
      :unload_accumulator_pipeline_spu.o (.text._SPUEAR_unload_accumulator_pipeline_spu)
      :load_interpolator_pipeline_spu.o (.text._SPUEAR_load_interpolator_pipeline_spu)
    }
    .ovly_other {
      :sort_p_pipeline_spu.o (.text._SPUEAR_coarse_sort_pipeline_spu)
      :sort_p_pipeline_spu.o (.text._SPUEAR_subsort_pipeline_spu)
      :sort_p_pipeline_spu.o (.text._SPUEAR_coarse_count_pipeline_spu)
    }
    .ovly_collision {
      :langevin_pipeline_spu.o (.text._SPUEAR_langevin_pipeline_spu)
      /opt/cell/sysroot/usr/spu/lib/libm.a:lib_a-wf_exp.o (.text.expf)
      /opt/cell/sysroot/usr/spu/lib/libm.a:lib_a-wf_log.o (.text.logf)
    } 
  }
}
INSERT AFTER .text;
