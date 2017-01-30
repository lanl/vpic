// FIXME: PARTICLE MOVERS NEED TO BE OVERALLOCATED IN STRUCTORS TO
// ACCOUNT FOR SPLITTING THE MOVER ARRAY BETWEEN HOST AND PIPELINES

#define VLEN 16
#define ALIGNMENT 64

#include "immintrin.h"
#define MOVE_P_ARRAY

const __m512 missing_ps = _mm512_set1_ps( -9999.9f );
const __m512i missing_i = _mm512_set1_epi32( -1 );

template <int VLEN_INT, int ALIGNMENT_INT=64>
class simd_mover_queue_t {/*{{{*/
    public:
        simd_mover_queue_t(){/*{{{*/
            _size = 0;
            _data_ps = (float *) _mm_malloc(2 * VLEN_INT * ps_fields * sizeof(float), ALIGNMENT_INT);
            _data_i = (int *) _mm_malloc(2 * VLEN_INT * i_fields * sizeof(int), ALIGNMENT_INT);
        }/*}}}*/

        ~simd_mover_queue_t(){/*{{{*/
            _size = 0;
            if ( _data_ps ) _mm_free(_data_ps);
            if ( _data_i ) _mm_free(_data_i);
        }/*}}}*/

        void enqueue(float v_x[VLEN_INT], float v_y[VLEN_INT], float v_z[VLEN_INT], int v_i[VLEN_INT]){/*{{{*/
            float* restrict d_x = &_data_ps[field_stride * x_id];
            float* restrict d_y = &_data_ps[field_stride * y_id];
            float* restrict d_z = &_data_ps[field_stride * z_id];
            int* restrict d_i = &_data_i[field_stride * i_id];
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
            float* restrict d_x = &_data_ps[field_stride * x_id];
            float* restrict d_y = &_data_ps[field_stride * y_id];
            float* restrict d_z = &_data_ps[field_stride * z_id];
            int* restrict d_i = &_data_i[field_stride * i_id];

            unsigned short mask_int = 0;

#pragma omp simd simdlen(VLEN_INT)
            for ( int v = 0; v < VLEN_INT; v++ )
            {
                mask_int ^= m[v] << v;
            }

            const __mmask16 mask = _mm512_int2mask( mask_int );

            const int count = _popcnt32(mask);

            _mm512_mask_compressstoreu_ps( &d_x[_size], mask, _mm512_load_ps( &v_x[0] ) );
            _mm512_mask_compressstoreu_ps( &d_y[_size], mask, _mm512_load_ps( &v_y[0] ) );
            _mm512_mask_compressstoreu_ps( &d_z[_size], mask, _mm512_load_ps( &v_z[0] ) );
            _mm512_mask_compressstoreu_epi32( &d_i[_size], mask, _mm512_load_epi32( &v_i[0] ) );

            _size += count;

            //assert(_size <= VLEN_INT);
        }/*}}}*/

        int dequeue(float v_x[VLEN_INT], float v_y[VLEN_INT], float v_z[VLEN_INT], int v_i[VLEN_INT]){/*{{{*/
            float* restrict d_x = &_data_ps[field_stride * x_id];
            float* restrict d_y = &_data_ps[field_stride * y_id];
            float* restrict d_z = &_data_ps[field_stride * z_id];
            int* restrict d_i = &_data_i[field_stride * i_id];
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
        const int field_stride = 2 * VLEN_INT;
        const int ps_fields = 3;
        const int i_fields = 1;
        enum ps_ids  {x_id, y_id, z_id};
        enum i_ids {i_id};

        float * _data_ps = NULL;
        int* _data_i;
};/*}}}*/

void move_p_array(
    particle_t* p0,
    float * mover_x, float * mover_y, float * mover_z, int * mover_i,
    accumulator_t* a0,
    const grid_t* g,
    const float qsp, bool * ret_vals)
{
    //printf("-- New call --\n");

    // Determine which particles need to be moved by looking for i == -1.
    // Assume that the movers won't be marked "inuse".
    int nactive = 0;
    bool active[VLEN], inuse[VLEN] __attribute__((aligned(ALIGNMENT)));

    float p_dx[VLEN], p_dy[VLEN], p_dz[VLEN], p_ux[VLEN], p_uy[VLEN], p_uz[VLEN], p_w[VLEN] __attribute__((aligned(ALIGNMENT)));
    int p_i[VLEN] __attribute__((aligned(ALIGNMENT)));

    #pragma omp simd reduction(+:nactive)
    for (int v = 0; v < VLEN; ++v)
    {
        active[v] = (mover_i[v] != -1);
        nactive += (active[v]) ? 1 : 0;
        inuse[v] = false;
        ret_vals[v] = false;
    }

    // Gather particle data.
    #pragma omp simd simdlen(VLEN)
    for ( int v = 0; v < VLEN; ++v)
    {
        int i = mover_i[v];
        p_dx[v] = (p0+i)->dx;
        p_dy[v] = (p0+i)->dy;
        p_dz[v] = (p0+i)->dz;

        p_ux[v] = (p0+i)->ux;
        p_uy[v] = (p0+i)->uy;
        p_uz[v] = (p0+i)->uz;
        p_w[v] = (p0+i)->w;
        p_i[v] = (p0+i)->i;
    }

    //printf("Starting iterations, nactive = %d\n", nactive);

    // Keep looping until every particle in the block has been moved to completion.
    //
    // Note that:
    // "continue" from original "for (;;)" is equivalent to "continue" here
    // "break" and "return" from original "for (;;)" are equivalent to "active[v]=false; continue;" here
    while (nactive > 0)
    {
        //printf(" --- New IT ---\n");
        #pragma omp simd simdlen(VLEN)
        for (int v = 0; v < VLEN; ++v)
        {

            if ( not active[v] ) continue;

            int i = mover_i[v];

            //float q = qsp * b->w[v];

            //float s_midx = b->dx[v];
            //float s_midy = b->dy[v];
            //float s_midz = b->dz[v];

            float q = qsp * p_w[v];
            float s_midx = p_dx[v];
            float s_midy = p_dy[v];
            float s_midz = p_dz[v];

            float s_dispx = mover_x[v];
            float s_dispy = mover_y[v];
            float s_dispz = mover_z[v];

            // s_dir unrolled because [axis] would have been a gather
            // expect it will be cheaper to update all three axes
            float s_dirx = (s_dispx > 0) ? 1 : -1;
            float s_diry = (s_dispy > 0) ? 1 : -1;
            float s_dirz = (s_dispz > 0) ? 1 : -1;

            // Compute the twice the fractional distance to each potential
            // streak/cell face intersection.
            float v0 = (s_dispx == 0) ? 3.4e38f : (s_dirx-s_midx)/s_dispx;
            float v1 = (s_dispy == 0) ? 3.4e38f : (s_diry-s_midy)/s_dispy;
            float v2 = (s_dispz == 0) ? 3.4e38f : (s_dirz-s_midz)/s_dispz;

            // Determine the fractional length and axis of current streak. The
            // streak ends on either the first face intersected by the
            // particle track or at the end of the particle track.
            float v3 = 2.0f;
            int axis = 3;
            if (v0 < v3) v3 = v0, axis = 0;
            if (v1 < v3) v3 = v1, axis = 1;
            if (v2 < v3) v3 = v2, axis = 2;
            v3 *= 0.5f;

            // Compute the midpoint and the normalized displacement of the streak
            s_dispx *= v3;
            s_dispy *= v3;
            s_dispz *= v3;
            s_midx += s_dispx;
            s_midy += s_dispy;
            s_midz += s_dispz;

            // Accumulate the streak.  Note: accumulator values are 4 times
            // the total physical charge that passed through the appropriate
            // current quadrant in a time-step
            float v4;
            float v5 = q*s_dispx*s_dispy*s_dispz*(1.0f/3.0f);
            //float* a = (a0 + b->i[v])->data;
            float* a = (float *)(a0 + p_i[v]);

            //v4  = q*s_dispx;    [> v2 = q ux                            <]
            //v1  = v4*s_midy;    [> v1 = q ux dy                         <]
            //v0  = v4-v1;          [> v0 = q ux (1-dy)                     <]
            //v1 += v4;             [> v1 = q ux (1+dy)                     <]
            //v4  = 1+s_midz;     [> v4 = 1+dz                            <]
            //v2  = v0*v4;          [> v2 = q ux (1-dy)(1+dz)               <]
            //v3  = v1*v4;          [> v3 = q ux (1+dy)(1+dz)               <]
            //v4  = 1-s_midz;     [> v4 = 1-dz                            <]
            //v0 *= v4;             [> v0 = q ux (1-dy)(1-dz)               <]
            //v1 *= v4;             [> v1 = q ux (1+dy)(1-dz)               <]
            //v0 += v5;             [> v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] <]
            //v1 -= v5;             [> v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] <]
            //v2 -= v5;             [> v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] <]
            //v3 += v5;             [> v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] <]
            //a[0] += v0;
            //a[1] += v1;
            //a[2] += v2;
            //a[3] += v3;


#   define accumulate_j(X,Y,Z)                                                \
            v4  = q*s_disp##X;    /* v2 = q ux                            */  \
            v1  = v4*s_mid##Y;    /* v1 = q ux dy                         */  \
            v0  = v4-v1;          /* v0 = q ux (1-dy)                     */  \
            v1 += v4;             /* v1 = q ux (1+dy)                     */  \
            v4  = 1+s_mid##Z;     /* v4 = 1+dz                            */  \
            v2  = v0*v4;          /* v2 = q ux (1-dy)(1+dz)               */  \
            v3  = v1*v4;          /* v3 = q ux (1+dy)(1+dz)               */  \
            v4  = 1-s_mid##Z;     /* v4 = 1-dz                            */  \
            v0 *= v4;             /* v0 = q ux (1-dy)(1-dz)               */  \
            v1 *= v4;             /* v1 = q ux (1+dy)(1-dz)               */  \
            v0 += v5;             /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */  \
            v1 -= v5;             /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */  \
            v2 -= v5;             /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */  \
            v3 += v5;             /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */  \
            a[0] += v0;                                                       \
            a[1] += v1;                                                       \
            a[2] += v2;                                                       \
            a[3] += v3
            accumulate_j(x,y,z); a += 4;
//            a = (a0 + b->i[v] + 4 - lowest)->data;
            accumulate_j(y,z,x); a += 4;
//            a = (a0 + b->i[v] + 8 - lowest)->data;
            accumulate_j(z,x,y);
#   undef accumulate_j

            // Compute the remaining particle displacment
            //pm->dispx[v] -= s_dispx;
            //pm->dispy[v] -= s_dispy;
            //pm->dispz[v] -= s_dispz;
            mover_x[v] -= s_dispx;
            mover_y[v] -= s_dispy;
            mover_z[v] -= s_dispz;

            // Compute the new particle offset
            //b->dx[v] += s_dispx + s_dispx;
            //b->dy[v] += s_dispy + s_dispy;
            //b->dz[v] += s_dispz + s_dispz;
            p_dx[v] += s_dispx + s_dispx;
            p_dy[v] += s_dispy + s_dispy;
            p_dz[v] += s_dispz + s_dispz;

            // If an end streak, return success (should be ~50% of the time)
            if (axis == 3)
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
            v0 = s_dirx;
            v1 = s_diry;
            v2 = s_dirz;

            // Avoid roundoff fiascos -- put the particle _exactly_ on the boundary.
            //b->dx[v] = (axis == 0) ? v0 : b->dx[v];
            //b->dy[v] = (axis == 1) ? v1 : b->dy[v];
            //b->dz[v] = (axis == 2) ? v2 : b->dz[v];

            p_dx[v] = (axis == 0) ? v0 : p_dx[v];
            p_dy[v] = (axis == 1) ? v1 : p_dy[v];
            p_dz[v] = (axis == 2) ? v2 : p_dz[v];

            int face = axis;
            face += (axis == 0 && v0 > 0) ? 3 : 0;
            face += (axis == 1 && v1 > 0) ? 3 : 0;
            face += (axis == 2 && v2 > 0) ? 3 : 0;
            //int neighbor = g->neighbor[ 6*b->i[v] + face ];
            int neighbor = g->neighbor[ 6 * p_i[v] + face ];

            // Hit a reflecting boundary condition.  Reflect the particle
            // momentum and remaining displacement and keep moving the
            // particle.
            //b->ux[v] *= (neighbor == reflect_particles && axis == 0) ? -1.0f : 1.0f;
            //b->uy[v] *= (neighbor == reflect_particles && axis == 1) ? -1.0f : 1.0f;
            //b->uz[v] *= (neighbor == reflect_particles && axis == 2) ? -1.0f : 1.0f;
            //pm->dispx[v] *= (neighbor == reflect_particles && axis == 0) ? -1.0f : 1.0f;
            //pm->dispy[v] *= (neighbor == reflect_particles && axis == 1) ? -1.0f : 1.0f;
            //pm->dispz[v] *= (neighbor == reflect_particles && axis == 2) ? -1.0f : 1.0f;

            p_ux[v] *= (neighbor == reflect_particles && axis == 0) ? -1.0f : 1.0f;
            p_uy[v] *= (neighbor == reflect_particles && axis == 1) ? -1.0f : 1.0f;
            p_uz[v] *= (neighbor == reflect_particles && axis == 2) ? -1.0f : 1.0f;
            mover_x[v] *= (neighbor == reflect_particles && axis == 0) ? -1.0f : 1.0f;
            mover_y[v] *= (neighbor == reflect_particles && axis == 1) ? -1.0f : 1.0f;
            mover_z[v] *= (neighbor == reflect_particles && axis == 2) ? -1.0f : 1.0f;
            if (neighbor == reflect_particles)
            {
                continue;
            }

            // Cannot handle the boundary condition here.  Save the updated
            // particle position, face it hit and update the remaining
            // displacement in the particle mover.
            //b->i[v] = (neighbor < g->rangel || neighbor > g->rangeh) ? 8*b->i[v] + face : b->i[v];
            p_i[v] = (neighbor < g->rangel || neighbor > g->rangeh) ? 8*p_i[v] + face : p_i[v];
            if (neighbor < g->rangel || neighbor > g->rangeh)
            {
                inuse[v] = true;
                active[v] = false;
                continue;
            }

            // Crossed into a normal voxel.  Update the voxel index, convert the
            // particle coordinate system and keep moving the particle.
            //b->i[v] = neighbor - g->rangel; // Compute local index of neighbor
            //[><]                            // Note: neighbor - g->rangel < 2^31/6
            //b->dx[v] = (axis == 0) ? -v0 : b->dx[v];
            //b->dy[v] = (axis == 1) ? -v1 : b->dy[v];
            //b->dz[v] = (axis == 2) ? -v2 : b->dz[v];

            p_i[v] = neighbor - g->rangel; // Compute local index of neighbor
            /**/                            // Note: neighbor - g->rangel < 2^31/6
            p_dx[v] = (axis == 0) ? -v0 : p_dx[v];
            p_dy[v] = (axis == 1) ? -v1 : p_dy[v];
            p_dz[v] = (axis == 2) ? -v2 : p_dz[v];

        }

        nactive = 0;
        #pragma omp simd simdlen(VLEN) reduction(+:nactive)
        for (int v = 0; v < VLEN; ++v)
        {
            nactive += (active[v]) ? 1 : 0;
        }

        //printf("Finished iteration, nactive = %d\n", nactive);
    }

    // Set return values, and return data to particle structures.
#pragma omp simd simdlen(VLEN)
    for ( int v = 0; v < VLEN; v++ )
    {
        int i = mover_i[v];

        ret_vals[v] = inuse[v];

        if ( mover_i[v] == -1 ) continue;

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
