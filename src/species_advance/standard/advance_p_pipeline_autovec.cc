#include "autovec_helper.h" 

#include "immintrin.h"

void
advance_p_pipeline( advance_p_pipeline_args_t * args,
                    int pipeline_rank,
                    int n_pipeline ) {/*{{{*/
    particle_t           * ALIGNED(128) p0 = args->p0;
    accumulator_t        * ALIGNED(128) a0 = args->a0;
    const interpolator_t * ALIGNED(128) f0 = args->f0;
    const grid_t *                      g  = args->g;

    particle_mover_t     * ALIGNED(16)  pm;
    const interpolator_t * ALIGNED(16)  f;
    float                * ALIGNED(16)  a;

    int n = args->np;

    const float qdt_2mc        = args->qdt_2mc;
    const float cdt_dx         = args->cdt_dx;
    const float cdt_dy         = args->cdt_dy;
    const float cdt_dz         = args->cdt_dz;
    const float qsp            = args->qsp;
    const float one            = 1.;
    const float one_third      = 1./3.;
    const float two_fifteenths = 2./15.;

    int ii;

    int itmp, nm, max_nm;

    float move_holders_x[VLEN];
    float move_holders_y[VLEN];
    float move_holders_z[VLEN];
    int move_holders_i[VLEN];

    const __m512i linear = _mm512_set_epi32( 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    const __m512 zeroes = _mm512_set1_ps( 0.0f );
    const int ind_stride = sizeof(interpolator_t) / sizeof(float);

    DECLARE_ALIGNED_ARRAY( particle_mover_t, 1, local_pm, 1 );

    // Determine which quads of particles quads this pipeline processes
    DISTRIBUTE( args->np, 16, pipeline_rank, n_pipeline, itmp, n );
    int particle_offset = itmp;
    particle_t* first_p = args->p0 + particle_offset;

    // Determine which movers are reserved for this pipeline
    // Movers (16 bytes) should be reserved for pipelines in at least
    // multiples of 8 such that the set of particle movers reserved for
    // a pipeline is 128-byte aligned and a multiple of 128-byte in
    // size.  The host is guaranteed to get enough movers to process its
    // particles with this allocation.

    max_nm = args->max_nm - (args->np&15);
    if( max_nm<0 ) max_nm = 0;
    DISTRIBUTE( max_nm, 8, pipeline_rank, n_pipeline, itmp, max_nm );
    if( pipeline_rank==n_pipeline ) max_nm = args->max_nm - itmp;
    pm   = args->pm + itmp;
    nm   = 0;
    int n_ignored = 0;

    //printf(" Processing particles: n = %d -- itmp = %d\n", n, itmp);

    // Determine which accumulator array to use
    // The host gets the first accumulator array

    if( pipeline_rank!=n_pipeline )
        a0 += (1+pipeline_rank)*
            POW2_CEIL((args->nx+2)*(args->ny+2)*(args->nz+2),2);

    simd_mover_queue_t<VLEN> move_queue;
    simd_mover_queue_t<VLEN> copy_queue;

    // Process particles for this pipeline
    //for (particle_t* p = first_p; p < first_p + n; p+=VLEN)
    for ( int i = 0; i < n; i+=VLEN )
    {
        particle_t * p = first_p + i;
        // Gather the field data.
        float ex[VLEN], dexdy[VLEN], dexdz[VLEN], d2exdydz[VLEN] __attribute__((aligned(ALIGNMENT)));
        float ey[VLEN], deydz[VLEN], deydx[VLEN], d2eydzdx[VLEN] __attribute__((aligned(ALIGNMENT)));
        float ez[VLEN], dezdx[VLEN], dezdy[VLEN], d2ezdxdy[VLEN] __attribute__((aligned(ALIGNMENT)));
        float _cbx[VLEN], dcbxdx[VLEN] __attribute__((aligned(ALIGNMENT)));
        float _cby[VLEN], dcbydy[VLEN] __attribute__((aligned(ALIGNMENT)));
        float _cbz[VLEN], dcbzdz[VLEN] __attribute__((aligned(ALIGNMENT)));

        bool outbnd[VLEN] __attribute__((aligned(ALIGNMENT)));
        float mya[12][VLEN] __attribute__((aligned(ALIGNMENT)));

//#define KNL 0

#if KNL
        // /*
        {
            //Intrinsic Gathers {{{
            // Each interpolator_t is made up of 20 floats.
            // Gather across the different interpolators into the arrays.
            // Need to do a masked gather, to prevent the last bit from gathering where there aren't particles.

            const int to_process = ( (i + VLEN) > n ) ? ( n - i ) : VLEN;
            int gather_indices_arr[VLEN] __attribute__((aligned(ALIGNMENT)));
            for ( int v = 0; v < VLEN; v++ ){
                if ( i + v >= n ) continue;
                gather_indices_arr[v] = (p+v)->i * ind_stride;
            }
            const __m512i gather_indices = _mm512_load_epi32( gather_indices_arr );

            // Build mask:
            const __mmask16 gather_mask = _mm512_cmplt_epi32_mask(linear, _mm512_set1_epi32(to_process));

            { // gathers, instead of load / transpose / store {{{
                // Perform gathers:
                { // ex {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].ex, 4);
                    _mm512_mask_store_ps(&ex[0], gather_mask, gathered_values);
                } // }}}
                { // dexdy {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].dexdy, 4);
                    _mm512_mask_store_ps(&dexdy[0], gather_mask, gathered_values);
                } // }}}
                { // dexdz {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].dexdz, 4);
                    _mm512_mask_store_ps(&dexdz[0], gather_mask, gathered_values);
                } // }}}
                { // d2exdydz {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].d2exdydz, 4);
                    _mm512_mask_store_ps(&d2exdydz[0], gather_mask, gathered_values);
                } // }}}

                { // ey {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].ey, 4);
                    _mm512_mask_store_ps(&ey[0], gather_mask, gathered_values);
                } // }}}
                { // deydz {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].deydz, 4);
                    _mm512_mask_store_ps(&deydz[0], gather_mask, gathered_values);
                } // }}}
                { // deydx {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].deydx, 4);
                    _mm512_mask_store_ps(&deydx[0], gather_mask, gathered_values);
                } // }}}
                { // d2eydzdx {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].d2eydzdx, 4);
                    _mm512_mask_store_ps(&d2eydzdx[0], gather_mask, gathered_values);
                } // }}}

                { // ez {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].ez, 4);
                    _mm512_mask_store_ps(&ez[0], gather_mask, gathered_values);
                } // }}}
                { // dezdx {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].dezdx, 4);
                    _mm512_mask_store_ps(&dezdx[0], gather_mask, gathered_values);
                } // }}}
                { // dezdy {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].dezdy, 4);
                    _mm512_mask_store_ps(&dezdy[0], gather_mask, gathered_values);
                } // }}}
                { // d2ezdxdy {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].d2ezdxdy, 4);
                    _mm512_mask_store_ps(&d2ezdxdy[0], gather_mask, gathered_values);
                } // }}}

                { // cbx {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].cbx, 4);
                    _mm512_mask_store_ps(&_cbx[0], gather_mask, gathered_values);
                } // }}}
                { // dcbxdx {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].dcbxdx, 4);
                    _mm512_mask_store_ps(&dcbxdx[0], gather_mask, gathered_values);
                } // }}}
                { // cby {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].cby, 4);
                    _mm512_mask_store_ps(&_cby[0], gather_mask, gathered_values);
                } // }}}
                { // dcbydy {{{
                    const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].dcbydy, 4);
                    _mm512_mask_store_ps(&dcbydy[0], gather_mask, gathered_values);
                } // }}}
            } // }}}

            // These two need to be gathers.
            { // cbz {{{
                const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].cbz, 4);
                _mm512_mask_store_ps(&_cbz[0], gather_mask, gathered_values);
            } // }}}
            { // dcbzdz {{{
                const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].dcbzdz, 4);
                _mm512_mask_store_ps(&dcbzdz[0], gather_mask, gathered_values);
            } // }}}
        } //  }}}
        // */
#else 
        { // Loops for gathers {{{
        //#pragma omp simd
#pragma novector
        for (int v = 0; v < VLEN; ++v)
        {
            //if ( p + v  >= first_p + n ) continue;
            if ( i + v >= n ) continue;
            int ii = (p+v)->i;
            ex[v] = f0[ii].ex;
            dexdy[v] = f0[ii].dexdy;
            dexdz[v] = f0[ii].dexdz;
            d2exdydz[v] = f0[ii].d2exdydz;
            ey[v] = f0[ii].ey;
            deydz[v] = f0[ii].deydz;
            deydx[v] = f0[ii].deydx;
            d2eydzdx[v] = f0[ii].d2eydzdx;
            ez[v] = f0[ii].ez;
            dezdx[v] = f0[ii].dezdx;
            dezdy[v] = f0[ii].dezdy;
            d2ezdxdy[v] = f0[ii].d2ezdxdy;
            _cbx[v] = f0[ii].cbx;
            dcbxdx[v] = f0[ii].dcbxdx;
            _cby[v] = f0[ii].cby;
            dcbydy[v] = f0[ii].dcbydy;
        }

        #pragma omp simd
        for (int v = 0; v < VLEN; ++v)
        {
            //if ( p + v  >= first_p + n ) continue;
            if ( i + v >= n ) continue;
            int ii = (p+v)->i;
            _cbz[v] = f0[ii].cbz;
            dcbzdz[v] = f0[ii].dcbzdz;
        }
        } // }}}
#endif

        #pragma omp simd
        for (int v = 0; v < VLEN; ++v)
        {
            outbnd[v] = false;

            if ( i + v >= n ) continue;
            int ii = (p+v)->i;

            float dx = (p+v)->dx;
            float dy = (p+v)->dy;
            float dz = (p+v)->dz;


            // TODO: this is a horrible way of writing f0[ii]
            const float hax  = qdt_2mc*(    ( ex[v]    + dy*dexdy[v]    ) +
                    dz*( dexdz[v] + dy*d2exdydz[v] ) );
            const float hay  = qdt_2mc*(    ( ey[v]    + dz*deydz[v]    ) +
                    dx*( deydx[v] + dz*d2eydzdx[v] ) );
            const float haz  = qdt_2mc*(    ( ez[v]    + dx*dezdx[v]    ) +
                    dy*( dezdy[v] + dx*d2ezdxdy[v] ) );

            float cbx = _cbx[v];
            float cby = _cby[v];
            float cbz = _cbz[v];
            float ux = (p+v)->ux;                                     // Load momentum
            float uy = (p+v)->uy;
            float uz = (p+v)->uz;
            float q = (p+v)->w;

            ux  += hax;                                           // Half advance E
            uy  += hay;
            uz  += haz;

            // TODO: There's really no reason to re-use these variable names...
            float v0, v1, v2, v3, v4, v5;
            v0   = qdt_2mc/sqrtf(one + (ux*ux + (uy*uy + uz*uz)));
            /**/                                      // Boris - scalars
            v1   = cbx*cbx + (cby*cby + cbz*cbz);
            v2   = (v0*v0)*v1;
            v3   = v0*(one+v2*(one_third+v2*two_fifteenths));
            v4   = v3/(one+v1*(v3*v3));
            v4  += v4;
            v0   = ux + v3*( uy*cbz - uz*cby );       // Boris - uprime
            v1   = uy + v3*( uz*cbx - ux*cbz );
            v2   = uz + v3*( ux*cby - uy*cbx );

            ux  += v4*( v1*cbz - v2*cby );            // Boris - rotation
            uy  += v4*( v2*cbx - v0*cbz );
            uz  += v4*( v0*cby - v1*cbx );

            ux  += hax;                               // Half advance E
            uy  += hay;
            uz  += haz;

            (p+v)->ux = ux;                               // Store momentum
            (p+v)->uy = uy;
            (p+v)->uz = uz;

            v0   = one/sqrtf(one + (ux*ux+ (uy*uy + uz*uz)));
            /**/                                      // Get norm displacement
            ux  *= cdt_dx;
            uy  *= cdt_dy;
            uz  *= cdt_dz;
            ux  *= v0;
            uy  *= v0;
            uz  *= v0;

            v0   = dx + ux;                           // Streak midpoint (inbnds)
            v1   = dy + uy;
            v2   = dz + uz;

            v3   = v0 + ux;                           // New position
            v4   = v1 + uy;
            v5   = v2 + uz;

            outbnd[v] = (v3 > 1.0f or v3 < -1.0f or v4 > 1.0f or v4 < -1.0f or v5 > 1.0f or v5 < -1.0f);

            float mask_val = (outbnd[v]) ? float(0.0) : float(1.0);

            // Accumulate current of inbnd particles
            q = mask_val * q * qsp;
            (p+v)->dx = mask_val * v3 + (1.0f - mask_val) * dx;
            (p+v)->dy = mask_val * v4 + (1.0f - mask_val) * dy;
            (p+v)->dz = mask_val * v5 + (1.0f - mask_val) * dz;

            dx = v0;                                // Streak midpoint
            dy = v1;
            dz = v2;
            v5 = q*ux*uy*uz*one_third;              // Compute correction

            // TODO: Better hack that doesn't rely on the ACCUMULATE macro.
#   define ACCUMULATE_J(X,Y,Z,offset)                                 \
            v4  = q*u##X;   /* v2 = q ux                            */        \
            v1  = v4*d##Y;  /* v1 = q ux dy                         */        \
            v0  = v4-v1;    /* v0 = q ux (1-dy)                     */        \
            v1 += v4;       /* v1 = q ux (1+dy)                     */        \
            v4  = one+d##Z; /* v4 = 1+dz                            */        \
            v2  = v0*v4;    /* v2 = q ux (1-dy)(1+dz)               */        \
            v3  = v1*v4;    /* v3 = q ux (1+dy)(1+dz)               */        \
            v4  = one-d##Z; /* v4 = 1-dz                            */        \
            v0 *= v4;       /* v0 = q ux (1-dy)(1-dz)               */        \
            v1 *= v4;       /* v1 = q ux (1+dy)(1-dz)               */        \
            v0 += v5;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */        \
            v1 -= v5;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */        \
            v2 -= v5;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */        \
            v3 += v5;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */        \
            mya[offset+0][v] = v0;                     \
            mya[offset+1][v] = v1;                     \
            mya[offset+2][v] = v2;                     \
            mya[offset+3][v] = v3;
            ACCUMULATE_J( x,y,z, 0 ); // accumulate ex, dexdy, dexdz, d2exdydz
            ACCUMULATE_J( y,z,x, 4 ); // accumulate ey, deydz, deydx, d2eydzdx
            ACCUMULATE_J( z,x,y, 8 ); // accumulate ez, dezdx, dezdy, d2ezdxdy
#   undef ACCUMULATE_J

            move_holders_x[v] = ux;
            move_holders_y[v] = uy;
            move_holders_z[v] = uz;
            move_holders_i[v] = particle_offset + i + v;
        }

        move_queue.enqueue(move_holders_x, move_holders_y, move_holders_z, move_holders_i, outbnd);

        // TODO: Look at the mya[12][VLEN] vs mya[VLEN][12]
        for (int v = 0; v < VLEN; ++v)
        {
            if ( i + v >= n ) continue;
            const int ii = (p+v)->i;
            float* a  = (float *)( a0 + ii );       // Get accumulator
            for (int o = 0; o < 16; ++o)
            {
                a[o] += mya[o][v];
            }
        }

        // /* // New move_p loop
        if ( move_queue.size() >= VLEN ) {
            int to_move = move_queue.dequeue(move_holders_x, move_holders_y, move_holders_z, move_holders_i);

            bool mover_returns[VLEN] __attribute__((aligned(ALIGNMENT)));

            move_p_array(p0, &move_holders_x[0], &move_holders_y[0], &move_holders_z[0], &move_holders_i[0],
                         a0, g, qsp, &mover_returns[0]);

            copy_queue.enqueue(move_holders_x, move_holders_y, move_holders_z, move_holders_i, mover_returns);

//#pragma vector always
//            for ( int v = 0; v < VLEN; v++ ) {
//                if ( mover_returns[v] ) {
//                    if ( nm < max_nm ) {
//                        pm[nm].dispx = move_holders_x[v];
//                        pm[nm].dispy = move_holders_y[v];
//                        pm[nm].dispz = move_holders_z[v];
//                        pm[nm].i = move_holders_i[v];
//                        nm++;
//                    } else {
//                        n_ignored++;
//                    }
//                }
//            }
        } // */

        if ( copy_queue.size() >= VLEN ) {
            int to_copy = copy_queue.dequeue(move_holders_x, move_holders_y, move_holders_z, move_holders_i);
            int ignored = 0;

#pragma omp simd simdlen(VLEN) reduction(+:ignored)
            for ( int v = 0; v < VLEN; v++ ){
                if ( nm + v < max_nm ) {
                    pm[nm + v].dispx = move_holders_x[v];
                    pm[nm + v].dispy = move_holders_y[v];
                    pm[nm + v].dispz = move_holders_z[v];
                    pm[nm + v].i = move_holders_i[v];
                } else {
                    ignored++;
                }
            }

            nm += (to_copy - ignored);
            n_ignored += ignored;
        }
    }

    // /* // New move_p loop
    while ( move_queue.size() > 0 ) {
        int to_move = move_queue.dequeue(move_holders_x, move_holders_y, move_holders_z, move_holders_i);

        bool mover_returns[VLEN] __attribute__((aligned(ALIGNMENT)));

        move_p_array(p0, &move_holders_x[0], &move_holders_y[0], &move_holders_z[0], &move_holders_i[0],
                a0, g, qsp, &mover_returns[0]);


        copy_queue.enqueue(move_holders_x, move_holders_y, move_holders_z, move_holders_i, mover_returns);

//#pragma vector always
//        for ( int v = 0; v < VLEN; v++ ) {
//            if ( mover_returns[v] ) {
//                if ( nm < max_nm ) {
//                    pm[nm].dispx = move_holders_x[v];
//                    pm[nm].dispy = move_holders_y[v];
//                    pm[nm].dispz = move_holders_z[v];
//                    pm[nm].i = move_holders_i[v];
//                    nm++;
//                } else {
//                    n_ignored++;
//                }
//            }
//        }
    } // */

    while ( copy_queue.size() > 0 ) {
        int to_copy = copy_queue.dequeue(move_holders_x, move_holders_y, move_holders_z, move_holders_i);
        int ignored = 0;

#pragma omp simd simdlen(VLEN) reduction(+:ignored)
        for ( int v = 0; v < VLEN; v++ ){
            if ( v >= to_copy ) continue;

            if ( nm + v < max_nm ) {
                pm[nm + v].dispx = move_holders_x[v];
                pm[nm + v].dispy = move_holders_y[v];
                pm[nm + v].dispz = move_holders_z[v];
                pm[nm + v].i = move_holders_i[v];
            } else {
                ignored++;
            }
        }

        nm += (to_copy - ignored);
        n_ignored += ignored;
    }

    args->seg[pipeline_rank].pm        = pm;
    args->seg[pipeline_rank].max_nm    = max_nm;
    args->seg[pipeline_rank].nm        = nm;
    args->seg[pipeline_rank].n_ignored = n_ignored;
}/*}}}*/
