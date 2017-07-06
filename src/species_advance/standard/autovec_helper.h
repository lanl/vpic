// FIXME: PARTICLE MOVERS NEED TO BE OVERALLOCATED IN STRUCTORS TO
// ACCOUNT FOR SPLITTING THE MOVER ARRAY BETWEEN HOST AND PIPELINES

#define VLEN 16
#define ALIGNMENT 64

#include "immintrin.h"
#define MOVE_P_ARRAY

template <int VLEN_INT, int ALIGNMENT_INT=64>
class simd_mover_queue_t {/*{{{*/
    public:
        simd_mover_queue_t(){/*{{{*/
            _size = 0;
            _data_x = (float *) _mm_malloc(2 * VLEN_INT * sizeof(float), ALIGNMENT_INT);
            _data_y = (float *) _mm_malloc(2 * VLEN_INT * sizeof(float), ALIGNMENT_INT);
            _data_z = (float *) _mm_malloc(2 * VLEN_INT * sizeof(float), ALIGNMENT_INT);
            _data_i = (int *) _mm_malloc(2 * VLEN_INT * sizeof(int), ALIGNMENT_INT);
        }/*}}}*/

        ~simd_mover_queue_t(){/*{{{*/
            _size = 0;
            if ( _data_x ) _mm_free(_data_x);
            if ( _data_y ) _mm_free(_data_y);
            if ( _data_z ) _mm_free(_data_z);
            if ( _data_i ) _mm_free(_data_i);
        }/*}}}*/

        void enqueue(float v_x[VLEN_INT], float v_y[VLEN_INT], float v_z[VLEN_INT], int v_i[VLEN_INT]){/*{{{*/
            float* restrict d_x = _data_x;
            float* restrict d_y = _data_y;
            float* restrict d_z = _data_z;
            int* restrict d_i = _data_i;
            int j = _size;

#pragma vector always
            for ( int i = 0; i < VLEN_INT; i++ ){
                   d_x[j + i] = v_x[i];
                   d_y[j + i] = v_y[i];
                   d_z[j + i] = v_z[i];
                   d_i[j + i] = v_i[i];
            }

            _size += VLEN_INT;

            //assert(_size <= VLEN_INT);
        }/*}}}*/

        void enqueue(float v_x[VLEN_INT], float v_y[VLEN_INT], float v_z[VLEN_INT], int v_i[VLEN_INT], bool m[VLEN_INT]){/*{{{*/
            float* restrict d_x = _data_x;
            float* restrict d_y = _data_y;
            float* restrict d_z = _data_z;
            int* restrict d_i = _data_i;

//            unsigned short mask_int = 0;
//
//#pragma omp simd simdlen(VLEN_INT)
//            for ( int v = 0; v < VLEN_INT; v++ )
//            {
//                mask_int ^= m[v] << v;
//            }
//
//            const __mmask16 mask = _mm512_int2mask( mask_int );
//
//            const int count = _popcnt32(mask);
//
//            _mm512_mask_compressstoreu_ps( &d_x[_size], mask, _mm512_load_ps( &v_x[0] ) );
//            _mm512_mask_compressstoreu_ps( &d_y[_size], mask, _mm512_load_ps( &v_y[0] ) );
//            _mm512_mask_compressstoreu_ps( &d_z[_size], mask, _mm512_load_ps( &v_z[0] ) );
//            _mm512_mask_compressstoreu_epi32( &d_i[_size], mask, _mm512_load_epi32( &v_i[0] ) );
//
//            _size += count;

            int j = _size;

#pragma vector always
            for ( int v = 0; v < VLEN_INT; v++ )
            {
                if ( m[v] ) {
                    d_x[j] = v_x[v];
                    d_y[j] = v_y[v];
                    d_z[j] = v_z[v];
                    d_i[j] = v_i[v];
                    j++;
                }
            }

            _size = j;

            //assert(_size <= VLEN_INT);
        }/*}}}*/

        int dequeue(float v_x[VLEN_INT], float v_y[VLEN_INT], float v_z[VLEN_INT], int v_i[VLEN_INT]){/*{{{*/
            float* restrict d_x = _data_x;
            float* restrict d_y = _data_y;
            float* restrict d_z = _data_z;
            int* restrict d_i = _data_i;
            int ndequeue = ( _size <= VLEN_INT ) ? _size : VLEN_INT;

#pragma vector always
            for ( int i = ndequeue; i < VLEN_INT; i++ ){
                v_x[i] = -9999;
                v_y[i] = -9999;
                v_z[i] = -9999;
                v_i[i] = -1;
            }

            ///*
#pragma vector always
            for ( int i = 0; i < ndequeue; i++){
                v_x[i] = d_x[i];
                v_y[i] = d_y[i];
                v_z[i] = d_z[i];
                v_i[i] = d_i[i];
            } // */
            //memcpy(&v[0], &d[0], sizeof(struct move_data_t) * ndequeue);

            _size -= ndequeue;

            if ( _size > 0 ) {
                // /*
#pragma vector always
                for (int i = 0; i < VLEN_INT; i++ ){
                    d_x[i] = d_x[i + VLEN_INT];
                    d_y[i] = d_y[i + VLEN_INT];
                    d_z[i] = d_z[i + VLEN_INT];
                    d_i[i] = d_i[i + VLEN_INT];
                }// */
                //memcpy(&d[0], &d[VLEN_INT], sizeof(struct move_data_t) * VLEN_INT);
            }

            return ndequeue;
        }/*}}}*/

        bool empty() const {/*{{{*/
            return ( _size == 0 );
        }/*}}}*/

        int size() const {/*{{{*/
            return _size;
        }/*}}}*/

    protected:
        int _size;
        float * _data_x = NULL;
        float * _data_y = NULL;
        float * _data_z = NULL;
        int* _data_i = NULL;
};/*}}}*/

void move_p_array(
    particle_t* p0,
    float * mover_x, float * mover_y, float * mover_z, int * mover_i,
    accumulator_t* a0,
    const grid_t* g,
    const float qsp, bool * ret_vals)
{
    // Determine which particles need to be moved by looking for i == -1.
    // Assume that the movers won't be marked "inuse".
    int nactive = 0;
    int active[VLEN], inuse[VLEN] __attribute__((aligned(ALIGNMENT)));

    float p_dx[VLEN], p_dy[VLEN], p_dz[VLEN], p_ux[VLEN], p_uy[VLEN], p_uz[VLEN], p_w[VLEN] __attribute__((aligned(ALIGNMENT)));
    int p_i[VLEN] __attribute__((aligned(ALIGNMENT)));
    float p_q[VLEN] __attribute__((aligned(ALIGNMENT)));

    __assume_aligned(mover_x, ALIGNMENT);
    __assume_aligned(mover_y, ALIGNMENT);
    __assume_aligned(mover_z, ALIGNMENT);
    __assume_aligned(mover_i, ALIGNMENT);
    __assume_aligned(ret_vals, ALIGNMENT);

    #pragma omp simd reduction(+:nactive)
    for (int v = 0; v < VLEN; ++v)
    {
        active[v] = (mover_i[v] >= 0);
        nactive += (active[v]) ? 1 : 0;
        inuse[v] = 0;
        ret_vals[v] = 0;
        p_i[v] = -1;
    }

    // Gather particle data.
    #pragma omp simd simdlen(VLEN)
    for ( int v = 0; v < VLEN; ++v)
    {
        if ( not active[v] ) continue;
        int i = mover_i[v];
        p_dx[v] = (p0+i)->dx;
        p_dy[v] = (p0+i)->dy;
        p_dz[v] = (p0+i)->dz;

        p_ux[v] = (p0+i)->ux;
        p_uy[v] = (p0+i)->uy;
        p_uz[v] = (p0+i)->uz;
        p_w[v] = (p0+i)->w;
        p_q[v] = (p0+i)->w * qsp;
        p_i[v] = (p0+i)->i;
    }

    // Keep looping until every particle in the block has been moved to completion.
    //
    // Note that:
    // "continue" from original "for (;;)" is equivalent to "continue" here
    // "break" and "return" from original "for (;;)" are equivalent to "active[v]=false; continue;" here
    while (nactive > 0)
    {
        float mya[VLEN][16] __attribute__((aligned(ALIGNMENT)));
        float s_midx[VLEN], s_midy[VLEN], s_midz[VLEN] __attribute__((aligned(ALIGNMENT)));
        float s_dispx[VLEN], s_dispy[VLEN], s_dispz[VLEN] __attribute__((aligned(ALIGNMENT)));
        float s_dirx[VLEN], s_diry[VLEN], s_dirz[VLEN] __attribute__((aligned(ALIGNMENT)));
        float v0[VLEN], v1[VLEN], v2[VLEN], v3[VLEN], v4[VLEN], v5[VLEN], v6[VLEN] __attribute__((aligned(ALIGNMENT)));
        int axis[VLEN] __attribute__((aligned(ALIGNMENT)));

        for ( int o = 0; o < 16; o++ )
        {
#pragma omp simd
            for ( int v = 0; v < 16; v++ )
            {
                mya[o][v] = 0.0f;
            }
        }

        #pragma omp simd simdlen(VLEN)
        for (int v = 0; v < VLEN; ++v)
        {

            if ( not active[v] ) continue;

            int i = mover_i[v];

            s_midx[v] = p_dx[v];
            s_midy[v] = p_dy[v];
            s_midz[v] = p_dz[v];

            s_dispx[v] = mover_x[v];
            s_dispy[v] = mover_y[v];
            s_dispz[v] = mover_z[v];

            // s_dir unrolled because [axis] would have been a gather
            // expect it will be cheaper to update all three axes
            s_dirx[v] = (s_dispx[v] > 0) ? 1 : -1;
            s_diry[v] = (s_dispy[v] > 0) ? 1 : -1;
            s_dirz[v] = (s_dispz[v] > 0) ? 1 : -1;

            // Compute the twice the fractional distance to each potential
            // streak/cell face intersection.
            v0[v] = (s_dispx[v] == 0) ? 3.4e38f : ( s_dirx[v] - s_midx[v] ) / s_dispx[v];
            v1[v] = (s_dispy[v] == 0) ? 3.4e38f : ( s_diry[v] - s_midy[v] ) / s_dispy[v];
            v2[v] = (s_dispz[v] == 0) ? 3.4e38f : ( s_dirz[v] - s_midz[v] ) / s_dispz[v];

            // Determine the fractional length and axis of current streak. The
            // streak ends on either the first face intersected by the
            // particle track or at the end of the particle track.
            v3[v] = 2.0f;
            axis[v] = 3;
            if (v0[v] < v3[v]) v3[v] = v0[v], axis[v] = 0;
            if (v1[v] < v3[v]) v3[v] = v1[v], axis[v] = 1;
            if (v2[v] < v3[v]) v3[v] = v2[v], axis[v] = 2;
            v3[v] *= 0.5f;

            // Compute the midpoint and the normalized displacement of the streak
            s_dispx[v] *= v3[v];
            s_dispy[v] *= v3[v];
            s_dispz[v] *= v3[v];
            s_midx[v] += s_dispx[v];
            s_midy[v] += s_dispy[v];
            s_midz[v] += s_dispz[v];

            // Accumulate the streak.  Note: accumulator values are 4 times
            // the total physical charge that passed through the appropriate
            // current quadrant in a time-step
            v5[v] = p_q[v] * s_dispx[v] * s_dispy[v] * s_dispz[v] * ( 1.0f / 3.0f );

#   define accumulate_j(X,Y,Z, offset)                                        \
            v4[v]  = p_q[v] * s_disp##X[v];    /* v2 = q ux                            */  \
            v1[v]  = v4[v] * s_mid##Y[v];    /* v1 = q ux dy                         */  \
            v0[v]  = v4[v] - v1[v];          /* v0 = q ux (1-dy)                     */  \
            v1[v] += v4[v];             /* v1 = q ux (1+dy)                     */  \
            v4[v]  = 1+s_mid##Z[v];     /* v4 = 1+dz                            */  \
            v2[v]  = v0[v] * v4[v];          /* v2 = q ux (1-dy)(1+dz)               */  \
            v3[v]  = v1[v] * v4[v];          /* v3 = q ux (1+dy)(1+dz)               */  \
            v4[v]  = 1 - s_mid##Z[v];     /* v4 = 1-dz                            */  \
            v0[v] *= v4[v];             /* v0 = q ux (1-dy)(1-dz)               */  \
            v1[v] *= v4[v];             /* v1 = q ux (1+dy)(1-dz)               */  \
            v0[v] += v5[v];             /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */  \
            v1[v] -= v5[v];             /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */  \
            v2[v] -= v5[v];             /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */  \
            v3[v] += v5[v];             /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */  \
            mya[v][offset+0] += v0[v];                                           \
            mya[v][offset+1] += v1[v];                                           \
            mya[v][offset+2] += v2[v];                                           \
            mya[v][offset+3] += v3[v];
            accumulate_j(x,y,z,0);
            accumulate_j(y,z,x,4);
            accumulate_j(z,x,y,8);
#   undef accumulate_j
        }

        for ( int v = 0; v < VLEN; v++ )
        {
            if ( not active[v] ) continue;

            float * restrict a = ( float * ) ( a0 + p_i[v] );

#pragma vector always
            for ( int o = 0; o < 16; o++ )
            {
                a[o] += mya[v][o];
            }
        }

#pragma omp simd simdlen(VLEN)
        for ( int v = 0; v < VLEN; v++ )
        {
            if ( not active[v] ) continue;
            // Compute the remaining particle displacment
            mover_x[v] -= s_dispx[v];
            mover_y[v] -= s_dispy[v];
            mover_z[v] -= s_dispz[v];

            // Compute the new particle offset
            p_dx[v] += s_dispx[v] + s_dispx[v];
            p_dy[v] += s_dispy[v] + s_dispy[v];
            p_dz[v] += s_dispz[v] + s_dispz[v];

            // If an end streak, return success (should be ~50% of the time)
            if (axis[v] == 3)
            {
                active[v] = false; // this lane is done
                continue;
            }

            // Determine if the particle crossed into a local cell or if it
            // hit a boundary and convert the coordinate system accordingly.
            // Note: Crossing into a local cell should happen ~50% of the
            // time; hitting a boundary is usually a rare event.  Note: the
            // entry / exit coordinate for the particle is guaranteed to be
            // +/-1 _exactly_ for the particle.
            v0[v] = s_dirx[v];
            v1[v] = s_diry[v];
            v2[v] = s_dirz[v];

            // Avoid roundoff fiascos -- put the particle _exactly_ on the boundary.
            p_dx[v] = (axis[v] == 0) ? v0[v] : p_dx[v];
            p_dy[v] = (axis[v] == 1) ? v1[v] : p_dy[v];
            p_dz[v] = (axis[v] == 2) ? v2[v] : p_dz[v];

            int face = axis[v];
            face += (axis[v] == 0 && v0[v] > 0) ? 3 : 0;
            face += (axis[v] == 1 && v1[v] > 0) ? 3 : 0;
            face += (axis[v] == 2 && v2[v] > 0) ? 3 : 0;
            //int neighbor = g->neighbor[ 6*b->i[v] + face ];
            int neighbor = g->neighbor[ 6 * p_i[v] + face ];

            // Hit a reflecting boundary condition.  Reflect the particle
            // momentum and remaining displacement and keep moving the
            // particle.
            p_ux[v] *= (neighbor == reflect_particles && axis[v] == 0) ? -1.0f : 1.0f;
            p_uy[v] *= (neighbor == reflect_particles && axis[v] == 1) ? -1.0f : 1.0f;
            p_uz[v] *= (neighbor == reflect_particles && axis[v] == 2) ? -1.0f : 1.0f;
            mover_x[v] *= (neighbor == reflect_particles && axis[v] == 0) ? -1.0f : 1.0f;
            mover_y[v] *= (neighbor == reflect_particles && axis[v] == 1) ? -1.0f : 1.0f;
            mover_z[v] *= (neighbor == reflect_particles && axis[v] == 2) ? -1.0f : 1.0f;
            if (neighbor == reflect_particles)
            {
                continue;
            }

            // Cannot handle the boundary condition here.  Save the updated
            // particle position, face it hit and update the remaining
            // displacement in the particle mover.
            p_i[v] = (neighbor < g->rangel || neighbor > g->rangeh) ? 8*p_i[v] + face : p_i[v];
            if (neighbor < g->rangel || neighbor > g->rangeh)
            {
                inuse[v] = true;
                active[v] = false;
                continue;
            }

            // Crossed into a normal voxel.  Update the voxel index, convert the
            // particle coordinate system and keep moving the particle.
            p_i[v] = neighbor - g->rangel; // Compute local index of neighbor
            /**/                            // Note: neighbor - g->rangel < 2^31/6
            p_dx[v] = (axis[v] == 0) ? -v0[v] : p_dx[v];
            p_dy[v] = (axis[v] == 1) ? -v1[v] : p_dy[v];
            p_dz[v] = (axis[v] == 2) ? -v2[v] : p_dz[v];
        }

        nactive = 0;
        #pragma omp simd simdlen(VLEN) reduction(+:nactive)
        for (int v = 0; v < VLEN; ++v)
        {
            nactive += (active[v]) ? 1 : 0;
        }
    }

    // Set return values, and return data to particle structures.
#pragma omp simd simdlen(VLEN)
    for ( int v = 0; v < VLEN; v++ )
    {
        int i = mover_i[v];

        ret_vals[v] = 0;

        if ( i < 0 ) continue;

        ret_vals[v] = inuse[v];

        (p0+i)->dx = p_dx[v];
        (p0+i)->dy = p_dy[v];
        (p0+i)->dz = p_dz[v];

        (p0+i)->ux = p_ux[v];
        (p0+i)->uy = p_uy[v];
        (p0+i)->uz = p_uz[v];
        (p0+i)->w = p_w[v];
        (p0+i)->i = p_i[v];
    }
}
