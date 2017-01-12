// FIXME: PARTICLE MOVERS NEED TO BE OVERALLOCATED IN STRUCTORS TO
// ACCOUNT FOR SPLITTING THE MOVER ARRAY BETWEEN HOST AND PIPELINES

#define VLEN 16
#define ALIGNMENT 64

template <int VLEN_INT, int ALIGNMENT_INT=64>
class simd_mover_queue_t {/*{{{*/
    public:
        simd_mover_queue_t(){/*{{{*/
            _size = 0;
            _data_x = (float*) _mm_malloc(2 * VLEN_INT * sizeof(float), ALIGNMENT_INT);
            _data_y = (float*) _mm_malloc(2 * VLEN_INT * sizeof(float), ALIGNMENT_INT);
            _data_z = (float*) _mm_malloc(2 * VLEN_INT * sizeof(float), ALIGNMENT_INT);
            _data_i = (int*) _mm_malloc(2 * VLEN_INT * sizeof(int), ALIGNMENT_INT);
        }/*}}}*/

        ~simd_mover_queue_t(){/*{{{*/
            _size = 0;
            if ( _data_x ) _mm_free(_data_x);
            if ( _data_y ) _mm_free(_data_y);
            if ( _data_z ) _mm_free(_data_z);
            if ( _data_i ) _mm_free(_data_i);
        }/*}}}*/

        void __declspec(noinline) enqueue(float v_x[VLEN_INT], float v_y[VLEN_INT], float v_z[VLEN_INT], int v_i[VLEN_INT]){/*{{{*/
            float* d_x = _data_x;
            float* d_y = _data_y;
            float* d_z = _data_z;
            int* d_i = _data_i;
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

        void __declspec(noinline) enqueue(float v_x[VLEN_INT], float v_y[VLEN_INT], float v_z[VLEN_INT], int v_i[VLEN_INT], bool m[VLEN_INT]){/*{{{*/
            float* d_x = _data_x;
            float* d_y = _data_y;
            float* d_z = _data_z;
            int* d_i = _data_i;
            int j = _size;

#pragma vector always
            for ( int i = 0; i < VLEN_INT; i++ ){
                if ( m[i] ) {
                    d_x[j] = v_x[i];
                    d_y[j] = v_y[i];
                    d_z[j] = v_z[i];
                    d_i[j] = v_i[i];
                    j++;
                }
            }

            _size = j;

            //assert(_size <= VLEN_INT);
        }/*}}}*/

        int __declspec(noinline) dequeue(float v_x[VLEN_INT], float v_y[VLEN_INT], float v_z[VLEN_INT], int v_i[VLEN_INT]){/*{{{*/
            float* d_x = _data_x;
            float* d_y = _data_y;
            float* d_z = _data_z;
            int* d_i = _data_i;
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
        float* _data_x = NULL;
        float* _data_y = NULL;
        float* _data_z = NULL;
        int* _data_i = NULL;
};/*}}}*/


