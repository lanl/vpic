#ifndef _v16_portable_h_
#define _v16_portable_h_

#ifndef IN_v16_h
#error "Do not include v16_portable.h directly; use v16.h"
#endif

#define V16_ACCELERATION
#define V16_PORTABLE_ACCELERATION

#include <math.h>

#ifndef ALIGNED
#define ALIGNED( n )
#endif

// This does not work with gcc 5.3.1 and the -fopenmp-simd
// flag.  Does not seem to work with -fopenmp either.  Not
// sure why.  It does work with the Intel compiler.  Need
// to try later versions of gcc.
// #define ALWAYS_VECTORIZE _Pragma( "omp simd" )

// #define ALWAYS_VECTORIZE _Pragma( "simd" )

#define ALWAYS_VECTORIZE _Pragma( "simd" ) _Pragma( "vector aligned" )

#define ALWAYS_INLINE __attribute__( ( always_inline ) )

namespace v16
{
class v16;
class v16int;
class v16float;

////////////////
// v16 base class

class v16
{
    friend class v16int;
    friend class v16float;

    // v16 miscellaneous friends

    friend inline int any( const v16& a ) ALWAYS_INLINE;
    friend inline int all( const v16& a ) ALWAYS_INLINE;

    template <int n>
    friend inline v16 splat( const v16& a ) ALWAYS_INLINE;

    template <int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7,
              int i8, int i9, int i10, int i11, int i12, int i13, int i14,
              int i15>
    friend inline v16 shuffle( const v16& a ) ALWAYS_INLINE;

    friend inline void swap( v16& a, v16& b ) ALWAYS_INLINE;
    friend inline void transpose( v16& a00, v16& a01, v16& a02, v16& a03,
                                  v16& a04, v16& a05, v16& a06, v16& a07,
                                  v16& a08, v16& a09, v16& a10, v16& a11,
                                  v16& a12, v16& a13, v16& a14,
                                  v16& a15 ) ALWAYS_INLINE;

    // v16int miscellaneous friends

    friend inline v16 czero( const v16int& c, const v16& a ) ALWAYS_INLINE;
    friend inline v16 notczero( const v16int& c, const v16& a ) ALWAYS_INLINE;
    friend inline v16 merge( const v16int& c, const v16& a,
                             const v16& b ) ALWAYS_INLINE;

    // v16 memory manipulation friends

    friend inline void load_16x1( const void* ALIGNED( 64 ) p,
                                  v16& a ) ALWAYS_INLINE;
    friend inline void store_16x1( const v16& a,
                                   void* ALIGNED( 64 ) p ) ALWAYS_INLINE;
    friend inline void stream_16x1( const v16& a,
                                    void* ALIGNED( 64 ) p ) ALWAYS_INLINE;
    friend inline void clear_16x1( void* ALIGNED( 64 ) dst ) ALWAYS_INLINE;
    friend inline void copy_16x1( void* ALIGNED( 64 ) dst,
                                  const void* ALIGNED( 64 ) src ) ALWAYS_INLINE;
    friend inline void swap_16x1( void* ALIGNED( 64 ) a,
                                  void* ALIGNED( 64 ) b ) ALWAYS_INLINE;

    // v16 transposed memory manipulation friends
    // Note: Half aligned values are permissible in the 16x2_tr variants.

    friend inline void
    load_16x1_tr( const void* a00, const void* a01, const void* a02,
                  const void* a03, const void* a04, const void* a05,
                  const void* a06, const void* a07, const void* a08,
                  const void* a09, const void* a10, const void* a11,
                  const void* a12, const void* a13, const void* a14,
                  const void* a15, v16& a ) ALWAYS_INLINE;
    friend inline void
    load_16x2_tr( const void* ALIGNED( 8 ) a00, const void* ALIGNED( 8 ) a01,
                  const void* ALIGNED( 8 ) a02, const void* ALIGNED( 8 ) a03,
                  const void* ALIGNED( 8 ) a04, const void* ALIGNED( 8 ) a05,
                  const void* ALIGNED( 8 ) a06, const void* ALIGNED( 8 ) a07,
                  const void* ALIGNED( 8 ) a08, const void* ALIGNED( 8 ) a09,
                  const void* ALIGNED( 8 ) a10, const void* ALIGNED( 8 ) a11,
                  const void* ALIGNED( 8 ) a12, const void* ALIGNED( 8 ) a13,
                  const void* ALIGNED( 8 ) a14, const void* ALIGNED( 8 ) a15,
                  v16& a, v16& b ) ALWAYS_INLINE;
    friend inline void
    load_16x3_tr( const void* ALIGNED( 64 ) a00, const void* ALIGNED( 64 ) a01,
                  const void* ALIGNED( 64 ) a02, const void* ALIGNED( 64 ) a03,
                  const void* ALIGNED( 64 ) a04, const void* ALIGNED( 64 ) a05,
                  const void* ALIGNED( 64 ) a06, const void* ALIGNED( 64 ) a07,
                  const void* ALIGNED( 64 ) a08, const void* ALIGNED( 64 ) a09,
                  const void* ALIGNED( 64 ) a10, const void* ALIGNED( 64 ) a11,
                  const void* ALIGNED( 64 ) a12, const void* ALIGNED( 64 ) a13,
                  const void* ALIGNED( 64 ) a14, const void* ALIGNED( 64 ) a15,
                  v16& a, v16& b, v16& c ) ALWAYS_INLINE;
    friend inline void
    load_16x4_tr( const void* ALIGNED( 64 ) a00, const void* ALIGNED( 64 ) a01,
                  const void* ALIGNED( 64 ) a02, const void* ALIGNED( 64 ) a03,
                  const void* ALIGNED( 64 ) a04, const void* ALIGNED( 64 ) a05,
                  const void* ALIGNED( 64 ) a06, const void* ALIGNED( 64 ) a07,
                  const void* ALIGNED( 64 ) a08, const void* ALIGNED( 64 ) a09,
                  const void* ALIGNED( 64 ) a10, const void* ALIGNED( 64 ) a11,
                  const void* ALIGNED( 64 ) a12, const void* ALIGNED( 64 ) a13,
                  const void* ALIGNED( 64 ) a14, const void* ALIGNED( 64 ) a15,
                  v16& a, v16& b, v16& c, v16& d ) ALWAYS_INLINE;
    friend inline void
    load_16x8_tr( const void* ALIGNED( 64 ) a00, const void* ALIGNED( 64 ) a01,
                  const void* ALIGNED( 64 ) a02, const void* ALIGNED( 64 ) a03,
                  const void* ALIGNED( 64 ) a04, const void* ALIGNED( 64 ) a05,
                  const void* ALIGNED( 64 ) a06, const void* ALIGNED( 64 ) a07,
                  const void* ALIGNED( 64 ) a08, const void* ALIGNED( 64 ) a09,
                  const void* ALIGNED( 64 ) a10, const void* ALIGNED( 64 ) a11,
                  const void* ALIGNED( 64 ) a12, const void* ALIGNED( 64 ) a13,
                  const void* ALIGNED( 64 ) a14, const void* ALIGNED( 64 ) a15,
                  v16& a, v16& b, v16& c, v16& d, v16& e, v16& f, v16& g,
                  v16& h ) ALWAYS_INLINE;
    friend inline void
    load_16x16_tr( const void* ALIGNED( 64 ) a00, const void* ALIGNED( 64 ) a01,
                   const void* ALIGNED( 64 ) a02, const void* ALIGNED( 64 ) a03,
                   const void* ALIGNED( 64 ) a04, const void* ALIGNED( 64 ) a05,
                   const void* ALIGNED( 64 ) a06, const void* ALIGNED( 64 ) a07,
                   const void* ALIGNED( 64 ) a08, const void* ALIGNED( 64 ) a09,
                   const void* ALIGNED( 64 ) a10, const void* ALIGNED( 64 ) a11,
                   const void* ALIGNED( 64 ) a12, const void* ALIGNED( 64 ) a13,
                   const void* ALIGNED( 64 ) a14, const void* ALIGNED( 64 ) a15,
                   v16& b00, v16& b01, v16& b02, v16& b03, v16& b04, v16& b05,
                   v16& b06, v16& b07, v16& b08, v16& b09, v16& b10, v16& b11,
                   v16& b12, v16& b13, v16& b14, v16& b15 ) ALWAYS_INLINE;
    friend inline void load_16x8_tr_p(
        const void* ALIGNED( 64 ) a00, const void* ALIGNED( 64 ) a01,
        const void* ALIGNED( 64 ) a02, const void* ALIGNED( 64 ) a03,
        const void* ALIGNED( 64 ) a04, const void* ALIGNED( 64 ) a05,
        const void* ALIGNED( 64 ) a06, const void* ALIGNED( 64 ) a07, v16& a,
        v16& b, v16& c, v16& d, v16& e, v16& f, v16& g, v16& h ) ALWAYS_INLINE;
    friend inline void load_16x16_tr_p(
        const void* ALIGNED( 64 ) a00, const void* ALIGNED( 64 ) a01,
        const void* ALIGNED( 64 ) a02, const void* ALIGNED( 64 ) a03,
        const void* ALIGNED( 64 ) a04, const void* ALIGNED( 64 ) a05,
        const void* ALIGNED( 64 ) a06, const void* ALIGNED( 64 ) a07,
        const void* ALIGNED( 64 ) a08, const void* ALIGNED( 64 ) a09,
        const void* ALIGNED( 64 ) a10, const void* ALIGNED( 64 ) a11,
        const void* ALIGNED( 64 ) a12, const void* ALIGNED( 64 ) a13,
        const void* ALIGNED( 64 ) a14, const void* ALIGNED( 64 ) a15, v16& b00,
        v16& b01, v16& b02, v16& b03, v16& b04, v16& b05, v16& b06, v16& b07,
        v16& b08, v16& b09, v16& b10, v16& b11, v16& b12, v16& b13, v16& b14,
        v16& b15 ) ALWAYS_INLINE;

    friend inline void store_16x1_tr( const v16& a, void* a00, void* a01,
                                      void* a02, void* a03, void* a04,
                                      void* a05, void* a06, void* a07,
                                      void* a08, void* a09, void* a10,
                                      void* a11, void* a12, void* a13,
                                      void* a14, void* a15 ) ALWAYS_INLINE;
    friend inline void store_16x2_tr(
        const v16& a, const v16& b, void* ALIGNED( 8 ) a00,
        void* ALIGNED( 8 ) a01, void* ALIGNED( 8 ) a02, void* ALIGNED( 8 ) a03,
        void* ALIGNED( 8 ) a04, void* ALIGNED( 8 ) a05, void* ALIGNED( 8 ) a06,
        void* ALIGNED( 8 ) a07, void* ALIGNED( 8 ) a08, void* ALIGNED( 8 ) a09,
        void* ALIGNED( 8 ) a10, void* ALIGNED( 8 ) a11, void* ALIGNED( 8 ) a12,
        void* ALIGNED( 8 ) a13, void* ALIGNED( 8 ) a14,
        void* ALIGNED( 8 ) a15 ) ALWAYS_INLINE;
    friend inline void
    store_16x3_tr( const v16& a, const v16& b, const v16& c,
                   void* ALIGNED( 64 ) a00, void* ALIGNED( 64 ) a01,
                   void* ALIGNED( 64 ) a02, void* ALIGNED( 64 ) a03,
                   void* ALIGNED( 64 ) a04, void* ALIGNED( 64 ) a05,
                   void* ALIGNED( 64 ) a06, void* ALIGNED( 64 ) a07,
                   void* ALIGNED( 64 ) a08, void* ALIGNED( 64 ) a09,
                   void* ALIGNED( 64 ) a10, void* ALIGNED( 64 ) a11,
                   void* ALIGNED( 64 ) a12, void* ALIGNED( 64 ) a13,
                   void* ALIGNED( 64 ) a14,
                   void* ALIGNED( 64 ) a15 ) ALWAYS_INLINE;
    friend inline void
    store_16x4_tr( const v16& a, const v16& b, const v16& c, const v16& d,
                   void* ALIGNED( 64 ) a00, void* ALIGNED( 64 ) a01,
                   void* ALIGNED( 64 ) a02, void* ALIGNED( 64 ) a03,
                   void* ALIGNED( 64 ) a04, void* ALIGNED( 64 ) a05,
                   void* ALIGNED( 64 ) a06, void* ALIGNED( 64 ) a07,
                   void* ALIGNED( 64 ) a08, void* ALIGNED( 64 ) a09,
                   void* ALIGNED( 64 ) a10, void* ALIGNED( 64 ) a11,
                   void* ALIGNED( 64 ) a12, void* ALIGNED( 64 ) a13,
                   void* ALIGNED( 64 ) a14,
                   void* ALIGNED( 64 ) a15 ) ALWAYS_INLINE;
    friend inline void
    store_16x8_tr( const v16& a, const v16& b, const v16& c, const v16& d,
                   const v16& e, const v16& f, const v16& g, const v16& h,
                   void* ALIGNED( 64 ) a00, void* ALIGNED( 64 ) a01,
                   void* ALIGNED( 64 ) a02, void* ALIGNED( 64 ) a03,
                   void* ALIGNED( 64 ) a04, void* ALIGNED( 64 ) a05,
                   void* ALIGNED( 64 ) a06, void* ALIGNED( 64 ) a07,
                   void* ALIGNED( 64 ) a08, void* ALIGNED( 64 ) a09,
                   void* ALIGNED( 64 ) a10, void* ALIGNED( 64 ) a11,
                   void* ALIGNED( 64 ) a12, void* ALIGNED( 64 ) a13,
                   void* ALIGNED( 64 ) a14,
                   void* ALIGNED( 64 ) a15 ) ALWAYS_INLINE;
    friend inline void store_16x16_tr(
        const v16& b00, const v16& b01, const v16& b02, const v16& b03,
        const v16& b04, const v16& b05, const v16& b06, const v16& b07,
        const v16& b08, const v16& b09, const v16& b10, const v16& b11,
        const v16& b12, const v16& b13, const v16& b14, const v16& b15,
        void* ALIGNED( 64 ) a00, void* ALIGNED( 64 ) a01,
        void* ALIGNED( 64 ) a02, void* ALIGNED( 64 ) a03,
        void* ALIGNED( 64 ) a04, void* ALIGNED( 64 ) a05,
        void* ALIGNED( 64 ) a06, void* ALIGNED( 64 ) a07,
        void* ALIGNED( 64 ) a08, void* ALIGNED( 64 ) a09,
        void* ALIGNED( 64 ) a10, void* ALIGNED( 64 ) a11,
        void* ALIGNED( 64 ) a12, void* ALIGNED( 64 ) a13,
        void* ALIGNED( 64 ) a14, void* ALIGNED( 64 ) a15 ) ALWAYS_INLINE;
    friend inline void
    store_16x8_tr_p( const v16& a, const v16& b, const v16& c, const v16& d,
                     const v16& e, const v16& f, const v16& g, const v16& h,
                     void* ALIGNED( 64 ) a00, void* ALIGNED( 64 ) a01,
                     void* ALIGNED( 64 ) a02, void* ALIGNED( 64 ) a03,
                     void* ALIGNED( 64 ) a04, void* ALIGNED( 64 ) a05,
                     void* ALIGNED( 64 ) a06,
                     void* ALIGNED( 64 ) a07 ) ALWAYS_INLINE;
    friend inline void store_16x16_tr_p(
        const v16& b00, const v16& b01, const v16& b02, const v16& b03,
        const v16& b04, const v16& b05, const v16& b06, const v16& b07,
        const v16& b08, const v16& b09, const v16& b10, const v16& b11,
        const v16& b12, const v16& b13, const v16& b14, const v16& b15,
        void* ALIGNED( 64 ) a00, void* ALIGNED( 64 ) a01,
        void* ALIGNED( 64 ) a02, void* ALIGNED( 64 ) a03,
        void* ALIGNED( 64 ) a04, void* ALIGNED( 64 ) a05,
        void* ALIGNED( 64 ) a06, void* ALIGNED( 64 ) a07,
        void* ALIGNED( 64 ) a08, void* ALIGNED( 64 ) a09,
        void* ALIGNED( 64 ) a10, void* ALIGNED( 64 ) a11,
        void* ALIGNED( 64 ) a12, void* ALIGNED( 64 ) a13,
        void* ALIGNED( 64 ) a14, void* ALIGNED( 64 ) a15 ) ALWAYS_INLINE;

  protected:
    union {
        int i[16];
        float f[16];
    };

  public:
    v16() {} // Default constructor

    v16( const v16& a ) // Copy constructor
    {
        ALWAYS_VECTORIZE
        for ( int j = 0; j < 16; j++ )
            i[j] = a.i[j];
    }

    ~v16() {} // Default destructor
};

// v16 miscellaneous functions

inline int any( const v16& a )
{
    return a.i[0] || a.i[1] || a.i[2] || a.i[3] || a.i[4] || a.i[5] || a.i[6] ||
           a.i[7] || a.i[8] || a.i[9] || a.i[10] || a.i[11] || a.i[12] ||
           a.i[13] || a.i[14] || a.i[15];
}

inline int all( const v16& a )
{
    return a.i[0] && a.i[1] && a.i[2] && a.i[3] && a.i[4] && a.i[5] && a.i[6] &&
           a.i[7] && a.i[8] && a.i[9] && a.i[10] && a.i[11] && a.i[12] &&
           a.i[13] && a.i[14] && a.i[15];
}

template <int n>
inline v16 splat( const v16& a )
{
    v16 b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.i[j] = a.i[n];

    return b;
}

template <int i00, int i01, int i02, int i03, int i04, int i05, int i06,
          int i07, int i08, int i09, int i10, int i11, int i12, int i13,
          int i14, int i15>
inline v16 shuffle( const v16& a )
{
    v16 b;

    b.i[0] = a.i[i00];
    b.i[1] = a.i[i01];
    b.i[2] = a.i[i02];
    b.i[3] = a.i[i03];
    b.i[4] = a.i[i04];
    b.i[5] = a.i[i05];
    b.i[6] = a.i[i06];
    b.i[7] = a.i[i07];
    b.i[8] = a.i[i08];
    b.i[9] = a.i[i09];
    b.i[10] = a.i[i10];
    b.i[11] = a.i[i11];
    b.i[12] = a.i[i12];
    b.i[13] = a.i[i13];
    b.i[14] = a.i[i14];
    b.i[15] = a.i[i15];

    return b;
}

#define sw( x, y ) x ^= y, y ^= x, x ^= y

inline void swap( v16& a, v16& b )
{
    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        sw( a.i[j], b.i[j] );
}

inline void transpose( v16& a00, v16& a01, v16& a02, v16& a03, v16& a04,
                       v16& a05, v16& a06, v16& a07, v16& a08, v16& a09,
                       v16& a10, v16& a11, v16& a12, v16& a13, v16& a14,
                       v16& a15 )
{
    sw( a00.i[1], a01.i[0] );
    sw( a00.i[2], a02.i[0] );
    sw( a00.i[3], a03.i[0] );
    sw( a00.i[4], a04.i[0] );
    sw( a00.i[5], a05.i[0] );
    sw( a00.i[6], a06.i[0] );
    sw( a00.i[7], a07.i[0] );
    sw( a00.i[8], a08.i[0] );
    sw( a00.i[9], a09.i[0] );
    sw( a00.i[10], a10.i[0] );
    sw( a00.i[11], a11.i[0] );
    sw( a00.i[12], a12.i[0] );
    sw( a00.i[13], a13.i[0] );
    sw( a00.i[14], a14.i[0] );
    sw( a00.i[15], a15.i[0] );
    sw( a01.i[2], a02.i[1] );
    sw( a01.i[3], a03.i[1] );
    sw( a01.i[4], a04.i[1] );
    sw( a01.i[5], a05.i[1] );
    sw( a01.i[6], a06.i[1] );
    sw( a01.i[7], a07.i[1] );
    sw( a01.i[8], a08.i[1] );
    sw( a01.i[9], a09.i[1] );
    sw( a01.i[10], a10.i[1] );
    sw( a01.i[11], a11.i[1] );
    sw( a01.i[12], a12.i[1] );
    sw( a01.i[13], a13.i[1] );
    sw( a01.i[14], a14.i[1] );
    sw( a01.i[15], a15.i[1] );
    sw( a02.i[3], a03.i[2] );
    sw( a02.i[4], a04.i[2] );
    sw( a02.i[5], a05.i[2] );
    sw( a02.i[6], a06.i[2] );
    sw( a02.i[7], a07.i[2] );
    sw( a02.i[8], a08.i[2] );
    sw( a02.i[9], a09.i[2] );
    sw( a02.i[10], a10.i[2] );
    sw( a02.i[11], a11.i[2] );
    sw( a02.i[12], a12.i[2] );
    sw( a02.i[13], a13.i[2] );
    sw( a02.i[14], a14.i[2] );
    sw( a02.i[15], a15.i[2] );
    sw( a03.i[4], a04.i[3] );
    sw( a03.i[5], a05.i[3] );
    sw( a03.i[6], a06.i[3] );
    sw( a03.i[7], a07.i[3] );
    sw( a03.i[8], a08.i[3] );
    sw( a03.i[9], a09.i[3] );
    sw( a03.i[10], a10.i[3] );
    sw( a03.i[11], a11.i[3] );
    sw( a03.i[12], a12.i[3] );
    sw( a03.i[13], a13.i[3] );
    sw( a03.i[14], a14.i[3] );
    sw( a03.i[15], a15.i[3] );
    sw( a04.i[5], a05.i[4] );
    sw( a04.i[6], a06.i[4] );
    sw( a04.i[7], a07.i[4] );
    sw( a04.i[8], a08.i[4] );
    sw( a04.i[9], a09.i[4] );
    sw( a04.i[10], a10.i[4] );
    sw( a04.i[11], a11.i[4] );
    sw( a04.i[12], a12.i[4] );
    sw( a04.i[13], a13.i[4] );
    sw( a04.i[14], a14.i[4] );
    sw( a04.i[15], a15.i[4] );
    sw( a05.i[6], a06.i[5] );
    sw( a05.i[7], a07.i[5] );
    sw( a05.i[8], a08.i[5] );
    sw( a05.i[9], a09.i[5] );
    sw( a05.i[10], a10.i[5] );
    sw( a05.i[11], a11.i[5] );
    sw( a05.i[12], a12.i[5] );
    sw( a05.i[13], a13.i[5] );
    sw( a05.i[14], a14.i[5] );
    sw( a05.i[15], a15.i[5] );
    sw( a06.i[7], a07.i[6] );
    sw( a06.i[8], a08.i[6] );
    sw( a06.i[9], a09.i[6] );
    sw( a06.i[10], a10.i[6] );
    sw( a06.i[11], a11.i[6] );
    sw( a06.i[12], a12.i[6] );
    sw( a06.i[13], a13.i[6] );
    sw( a06.i[14], a14.i[6] );
    sw( a06.i[15], a15.i[6] );
    sw( a07.i[8], a08.i[7] );
    sw( a07.i[9], a09.i[7] );
    sw( a07.i[10], a10.i[7] );
    sw( a07.i[11], a11.i[7] );
    sw( a07.i[12], a12.i[7] );
    sw( a07.i[13], a13.i[7] );
    sw( a07.i[14], a14.i[7] );
    sw( a07.i[15], a15.i[7] );
    sw( a08.i[9], a09.i[8] );
    sw( a08.i[10], a10.i[8] );
    sw( a08.i[11], a11.i[8] );
    sw( a08.i[12], a12.i[8] );
    sw( a08.i[13], a13.i[8] );
    sw( a08.i[14], a14.i[8] );
    sw( a08.i[15], a15.i[8] );
    sw( a09.i[10], a10.i[9] );
    sw( a09.i[11], a11.i[9] );
    sw( a09.i[12], a12.i[9] );
    sw( a09.i[13], a13.i[9] );
    sw( a09.i[14], a14.i[9] );
    sw( a09.i[15], a15.i[9] );
    sw( a10.i[11], a11.i[10] );
    sw( a10.i[12], a12.i[10] );
    sw( a10.i[13], a13.i[10] );
    sw( a10.i[14], a14.i[10] );
    sw( a10.i[15], a15.i[10] );
    sw( a11.i[12], a12.i[11] );
    sw( a11.i[13], a13.i[11] );
    sw( a11.i[14], a14.i[11] );
    sw( a11.i[15], a15.i[11] );
    sw( a12.i[13], a13.i[12] );
    sw( a12.i[14], a14.i[12] );
    sw( a12.i[15], a15.i[12] );
    sw( a13.i[14], a14.i[13] );
    sw( a13.i[15], a15.i[13] );
    sw( a14.i[15], a15.i[14] );
}

#undef sw

// v16 memory manipulation functions

inline void load_16x1( const void* ALIGNED( 64 ) p, v16& a )
{
    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        a.i[j] = ( (const int* ALIGNED( 64 ))p )[j];
}

inline void store_16x1( const v16& a, void* ALIGNED( 64 ) p )
{
    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        ( (int* ALIGNED( 64 ))p )[j] = a.i[j];
}

inline void stream_16x1( const v16& a, void* ALIGNED( 64 ) p )
{
    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        ( (int* ALIGNED( 64 ))p )[j] = a.i[j];
}

inline void clear_16x1( void* ALIGNED( 64 ) p )
{
    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        ( (int* ALIGNED( 64 ))p )[j] = 0;
}

// FIXME: Ordering semantics
inline void copy_16x1( void* ALIGNED( 64 ) dst, const void* ALIGNED( 64 ) src )
{
    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        ( (int* ALIGNED( 64 ))dst )[j] = ( (const int* ALIGNED( 64 ))src )[j];
}

inline void swap_16x1( void* ALIGNED( 64 ) a, void* ALIGNED( 64 ) b )
{
    int t;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
    {
        t = ( (int* ALIGNED( 64 ))a )[j];
        ( (int* ALIGNED( 64 ))a )[j] = ( (int* ALIGNED( 64 ))b )[j];
        ( (int* ALIGNED( 64 ))b )[j] = t;
    }
}

// v16 transposed memory manipulation functions

inline void load_16x1_tr( const void* a00, const void* a01, const void* a02,
                          const void* a03, const void* a04, const void* a05,
                          const void* a06, const void* a07, const void* a08,
                          const void* a09, const void* a10, const void* a11,
                          const void* a12, const void* a13, const void* a14,
                          const void* a15, v16& a )
{
    a.i[0] = ( (const int*)a00 )[0];
    a.i[1] = ( (const int*)a01 )[0];
    a.i[2] = ( (const int*)a02 )[0];
    a.i[3] = ( (const int*)a03 )[0];
    a.i[4] = ( (const int*)a04 )[0];
    a.i[5] = ( (const int*)a05 )[0];
    a.i[6] = ( (const int*)a06 )[0];
    a.i[7] = ( (const int*)a07 )[0];
    a.i[8] = ( (const int*)a08 )[0];
    a.i[9] = ( (const int*)a09 )[0];
    a.i[10] = ( (const int*)a10 )[0];
    a.i[11] = ( (const int*)a11 )[0];
    a.i[12] = ( (const int*)a12 )[0];
    a.i[13] = ( (const int*)a13 )[0];
    a.i[14] = ( (const int*)a14 )[0];
    a.i[15] = ( (const int*)a15 )[0];
}

inline void
load_16x2_tr( const void* ALIGNED( 8 ) a00, const void* ALIGNED( 8 ) a01,
              const void* ALIGNED( 8 ) a02, const void* ALIGNED( 8 ) a03,
              const void* ALIGNED( 8 ) a04, const void* ALIGNED( 8 ) a05,
              const void* ALIGNED( 8 ) a06, const void* ALIGNED( 8 ) a07,
              const void* ALIGNED( 8 ) a08, const void* ALIGNED( 8 ) a09,
              const void* ALIGNED( 8 ) a10, const void* ALIGNED( 8 ) a11,
              const void* ALIGNED( 8 ) a12, const void* ALIGNED( 8 ) a13,
              const void* ALIGNED( 8 ) a14, const void* ALIGNED( 8 ) a15,
              v16& a, v16& b )
{
    a.i[0] = ( (const int* ALIGNED( 8 ))a00 )[0];
    b.i[0] = ( (const int* ALIGNED( 8 ))a00 )[1];

    a.i[1] = ( (const int* ALIGNED( 8 ))a01 )[0];
    b.i[1] = ( (const int* ALIGNED( 8 ))a01 )[1];

    a.i[2] = ( (const int* ALIGNED( 8 ))a02 )[0];
    b.i[2] = ( (const int* ALIGNED( 8 ))a02 )[1];

    a.i[3] = ( (const int* ALIGNED( 8 ))a03 )[0];
    b.i[3] = ( (const int* ALIGNED( 8 ))a03 )[1];

    a.i[4] = ( (const int* ALIGNED( 8 ))a04 )[0];
    b.i[4] = ( (const int* ALIGNED( 8 ))a04 )[1];

    a.i[5] = ( (const int* ALIGNED( 8 ))a05 )[0];
    b.i[5] = ( (const int* ALIGNED( 8 ))a05 )[1];

    a.i[6] = ( (const int* ALIGNED( 8 ))a06 )[0];
    b.i[6] = ( (const int* ALIGNED( 8 ))a06 )[1];

    a.i[7] = ( (const int* ALIGNED( 8 ))a07 )[0];
    b.i[7] = ( (const int* ALIGNED( 8 ))a07 )[1];

    a.i[8] = ( (const int* ALIGNED( 8 ))a08 )[0];
    b.i[8] = ( (const int* ALIGNED( 8 ))a08 )[1];

    a.i[9] = ( (const int* ALIGNED( 8 ))a09 )[0];
    b.i[9] = ( (const int* ALIGNED( 8 ))a09 )[1];

    a.i[10] = ( (const int* ALIGNED( 8 ))a10 )[0];
    b.i[10] = ( (const int* ALIGNED( 8 ))a10 )[1];

    a.i[11] = ( (const int* ALIGNED( 8 ))a11 )[0];
    b.i[11] = ( (const int* ALIGNED( 8 ))a11 )[1];

    a.i[12] = ( (const int* ALIGNED( 8 ))a12 )[0];
    b.i[12] = ( (const int* ALIGNED( 8 ))a12 )[1];

    a.i[13] = ( (const int* ALIGNED( 8 ))a13 )[0];
    b.i[13] = ( (const int* ALIGNED( 8 ))a13 )[1];

    a.i[14] = ( (const int* ALIGNED( 8 ))a14 )[0];
    b.i[14] = ( (const int* ALIGNED( 8 ))a14 )[1];

    a.i[15] = ( (const int* ALIGNED( 8 ))a15 )[0];
    b.i[15] = ( (const int* ALIGNED( 8 ))a15 )[1];
}

inline void
load_16x3_tr( const void* ALIGNED( 64 ) a00, const void* ALIGNED( 64 ) a01,
              const void* ALIGNED( 64 ) a02, const void* ALIGNED( 64 ) a03,
              const void* ALIGNED( 64 ) a04, const void* ALIGNED( 64 ) a05,
              const void* ALIGNED( 64 ) a06, const void* ALIGNED( 64 ) a07,
              const void* ALIGNED( 64 ) a08, const void* ALIGNED( 64 ) a09,
              const void* ALIGNED( 64 ) a10, const void* ALIGNED( 64 ) a11,
              const void* ALIGNED( 64 ) a12, const void* ALIGNED( 64 ) a13,
              const void* ALIGNED( 64 ) a14, const void* ALIGNED( 64 ) a15,
              v16& a, v16& b, v16& c )
{
    a.i[0] = ( (const int* ALIGNED( 64 ))a00 )[0];
    b.i[0] = ( (const int* ALIGNED( 64 ))a00 )[1];
    c.i[0] = ( (const int* ALIGNED( 64 ))a00 )[2];

    a.i[1] = ( (const int* ALIGNED( 64 ))a01 )[0];
    b.i[1] = ( (const int* ALIGNED( 64 ))a01 )[1];
    c.i[1] = ( (const int* ALIGNED( 64 ))a01 )[2];

    a.i[2] = ( (const int* ALIGNED( 64 ))a02 )[0];
    b.i[2] = ( (const int* ALIGNED( 64 ))a02 )[1];
    c.i[2] = ( (const int* ALIGNED( 64 ))a02 )[2];

    a.i[3] = ( (const int* ALIGNED( 64 ))a03 )[0];
    b.i[3] = ( (const int* ALIGNED( 64 ))a03 )[1];
    c.i[3] = ( (const int* ALIGNED( 64 ))a03 )[2];

    a.i[4] = ( (const int* ALIGNED( 64 ))a04 )[0];
    b.i[4] = ( (const int* ALIGNED( 64 ))a04 )[1];
    c.i[4] = ( (const int* ALIGNED( 64 ))a04 )[2];

    a.i[5] = ( (const int* ALIGNED( 64 ))a05 )[0];
    b.i[5] = ( (const int* ALIGNED( 64 ))a05 )[1];
    c.i[5] = ( (const int* ALIGNED( 64 ))a05 )[2];

    a.i[6] = ( (const int* ALIGNED( 64 ))a06 )[0];
    b.i[6] = ( (const int* ALIGNED( 64 ))a06 )[1];
    c.i[6] = ( (const int* ALIGNED( 64 ))a06 )[2];

    a.i[7] = ( (const int* ALIGNED( 64 ))a07 )[0];
    b.i[7] = ( (const int* ALIGNED( 64 ))a07 )[1];
    c.i[7] = ( (const int* ALIGNED( 64 ))a07 )[2];

    a.i[8] = ( (const int* ALIGNED( 64 ))a08 )[0];
    b.i[8] = ( (const int* ALIGNED( 64 ))a08 )[1];
    c.i[8] = ( (const int* ALIGNED( 64 ))a08 )[2];

    a.i[9] = ( (const int* ALIGNED( 64 ))a09 )[0];
    b.i[9] = ( (const int* ALIGNED( 64 ))a09 )[1];
    c.i[9] = ( (const int* ALIGNED( 64 ))a09 )[2];

    a.i[10] = ( (const int* ALIGNED( 64 ))a10 )[0];
    b.i[10] = ( (const int* ALIGNED( 64 ))a10 )[1];
    c.i[10] = ( (const int* ALIGNED( 64 ))a10 )[2];

    a.i[11] = ( (const int* ALIGNED( 64 ))a11 )[0];
    b.i[11] = ( (const int* ALIGNED( 64 ))a11 )[1];
    c.i[11] = ( (const int* ALIGNED( 64 ))a11 )[2];

    a.i[12] = ( (const int* ALIGNED( 64 ))a12 )[0];
    b.i[12] = ( (const int* ALIGNED( 64 ))a12 )[1];
    c.i[12] = ( (const int* ALIGNED( 64 ))a12 )[2];

    a.i[13] = ( (const int* ALIGNED( 64 ))a13 )[0];
    b.i[13] = ( (const int* ALIGNED( 64 ))a13 )[1];
    c.i[13] = ( (const int* ALIGNED( 64 ))a13 )[2];

    a.i[14] = ( (const int* ALIGNED( 64 ))a14 )[0];
    b.i[14] = ( (const int* ALIGNED( 64 ))a14 )[1];
    c.i[14] = ( (const int* ALIGNED( 64 ))a14 )[2];

    a.i[15] = ( (const int* ALIGNED( 64 ))a15 )[0];
    b.i[15] = ( (const int* ALIGNED( 64 ))a15 )[1];
    c.i[15] = ( (const int* ALIGNED( 64 ))a15 )[2];
}

inline void
load_16x4_tr( const void* ALIGNED( 64 ) a00, const void* ALIGNED( 64 ) a01,
              const void* ALIGNED( 64 ) a02, const void* ALIGNED( 64 ) a03,
              const void* ALIGNED( 64 ) a04, const void* ALIGNED( 64 ) a05,
              const void* ALIGNED( 64 ) a06, const void* ALIGNED( 64 ) a07,
              const void* ALIGNED( 64 ) a08, const void* ALIGNED( 64 ) a09,
              const void* ALIGNED( 64 ) a10, const void* ALIGNED( 64 ) a11,
              const void* ALIGNED( 64 ) a12, const void* ALIGNED( 64 ) a13,
              const void* ALIGNED( 64 ) a14, const void* ALIGNED( 64 ) a15,
              v16& a, v16& b, v16& c, v16& d )
{
    a.i[0] = ( (const int* ALIGNED( 64 ))a00 )[0];
    b.i[0] = ( (const int* ALIGNED( 64 ))a00 )[1];
    c.i[0] = ( (const int* ALIGNED( 64 ))a00 )[2];
    d.i[0] = ( (const int* ALIGNED( 64 ))a00 )[3];

    a.i[1] = ( (const int* ALIGNED( 64 ))a01 )[0];
    b.i[1] = ( (const int* ALIGNED( 64 ))a01 )[1];
    c.i[1] = ( (const int* ALIGNED( 64 ))a01 )[2];
    d.i[1] = ( (const int* ALIGNED( 64 ))a01 )[3];

    a.i[2] = ( (const int* ALIGNED( 64 ))a02 )[0];
    b.i[2] = ( (const int* ALIGNED( 64 ))a02 )[1];
    c.i[2] = ( (const int* ALIGNED( 64 ))a02 )[2];
    d.i[2] = ( (const int* ALIGNED( 64 ))a02 )[3];

    a.i[3] = ( (const int* ALIGNED( 64 ))a03 )[0];
    b.i[3] = ( (const int* ALIGNED( 64 ))a03 )[1];
    c.i[3] = ( (const int* ALIGNED( 64 ))a03 )[2];
    d.i[3] = ( (const int* ALIGNED( 64 ))a03 )[3];

    a.i[4] = ( (const int* ALIGNED( 64 ))a04 )[0];
    b.i[4] = ( (const int* ALIGNED( 64 ))a04 )[1];
    c.i[4] = ( (const int* ALIGNED( 64 ))a04 )[2];
    d.i[4] = ( (const int* ALIGNED( 64 ))a04 )[3];

    a.i[5] = ( (const int* ALIGNED( 64 ))a05 )[0];
    b.i[5] = ( (const int* ALIGNED( 64 ))a05 )[1];
    c.i[5] = ( (const int* ALIGNED( 64 ))a05 )[2];
    d.i[5] = ( (const int* ALIGNED( 64 ))a05 )[3];

    a.i[6] = ( (const int* ALIGNED( 64 ))a06 )[0];
    b.i[6] = ( (const int* ALIGNED( 64 ))a06 )[1];
    c.i[6] = ( (const int* ALIGNED( 64 ))a06 )[2];
    d.i[6] = ( (const int* ALIGNED( 64 ))a06 )[3];

    a.i[7] = ( (const int* ALIGNED( 64 ))a07 )[0];
    b.i[7] = ( (const int* ALIGNED( 64 ))a07 )[1];
    c.i[7] = ( (const int* ALIGNED( 64 ))a07 )[2];
    d.i[7] = ( (const int* ALIGNED( 64 ))a07 )[3];

    a.i[8] = ( (const int* ALIGNED( 64 ))a08 )[0];
    b.i[8] = ( (const int* ALIGNED( 64 ))a08 )[1];
    c.i[8] = ( (const int* ALIGNED( 64 ))a08 )[2];
    d.i[8] = ( (const int* ALIGNED( 64 ))a08 )[3];

    a.i[9] = ( (const int* ALIGNED( 64 ))a09 )[0];
    b.i[9] = ( (const int* ALIGNED( 64 ))a09 )[1];
    c.i[9] = ( (const int* ALIGNED( 64 ))a09 )[2];
    d.i[9] = ( (const int* ALIGNED( 64 ))a09 )[3];

    a.i[10] = ( (const int* ALIGNED( 64 ))a10 )[0];
    b.i[10] = ( (const int* ALIGNED( 64 ))a10 )[1];
    c.i[10] = ( (const int* ALIGNED( 64 ))a10 )[2];
    d.i[10] = ( (const int* ALIGNED( 64 ))a10 )[3];

    a.i[11] = ( (const int* ALIGNED( 64 ))a11 )[0];
    b.i[11] = ( (const int* ALIGNED( 64 ))a11 )[1];
    c.i[11] = ( (const int* ALIGNED( 64 ))a11 )[2];
    d.i[11] = ( (const int* ALIGNED( 64 ))a11 )[3];

    a.i[12] = ( (const int* ALIGNED( 64 ))a12 )[0];
    b.i[12] = ( (const int* ALIGNED( 64 ))a12 )[1];
    c.i[12] = ( (const int* ALIGNED( 64 ))a12 )[2];
    d.i[12] = ( (const int* ALIGNED( 64 ))a12 )[3];

    a.i[13] = ( (const int* ALIGNED( 64 ))a13 )[0];
    b.i[13] = ( (const int* ALIGNED( 64 ))a13 )[1];
    c.i[13] = ( (const int* ALIGNED( 64 ))a13 )[2];
    d.i[13] = ( (const int* ALIGNED( 64 ))a13 )[3];

    a.i[14] = ( (const int* ALIGNED( 64 ))a14 )[0];
    b.i[14] = ( (const int* ALIGNED( 64 ))a14 )[1];
    c.i[14] = ( (const int* ALIGNED( 64 ))a14 )[2];
    d.i[14] = ( (const int* ALIGNED( 64 ))a14 )[3];

    a.i[15] = ( (const int* ALIGNED( 64 ))a15 )[0];
    b.i[15] = ( (const int* ALIGNED( 64 ))a15 )[1];
    c.i[15] = ( (const int* ALIGNED( 64 ))a15 )[2];
    d.i[15] = ( (const int* ALIGNED( 64 ))a15 )[3];
}

inline void
load_16x8_tr( const void* ALIGNED( 64 ) a00, const void* ALIGNED( 64 ) a01,
              const void* ALIGNED( 64 ) a02, const void* ALIGNED( 64 ) a03,
              const void* ALIGNED( 64 ) a04, const void* ALIGNED( 64 ) a05,
              const void* ALIGNED( 64 ) a06, const void* ALIGNED( 64 ) a07,
              const void* ALIGNED( 64 ) a08, const void* ALIGNED( 64 ) a09,
              const void* ALIGNED( 64 ) a10, const void* ALIGNED( 64 ) a11,
              const void* ALIGNED( 64 ) a12, const void* ALIGNED( 64 ) a13,
              const void* ALIGNED( 64 ) a14, const void* ALIGNED( 64 ) a15,
              v16& a, v16& b, v16& c, v16& d, v16& e, v16& f, v16& g, v16& h )
{
    a.i[0] = ( (const int* ALIGNED( 64 ))a00 )[0];
    b.i[0] = ( (const int* ALIGNED( 64 ))a00 )[1];
    c.i[0] = ( (const int* ALIGNED( 64 ))a00 )[2];
    d.i[0] = ( (const int* ALIGNED( 64 ))a00 )[3];
    e.i[0] = ( (const int* ALIGNED( 64 ))a00 )[4];
    f.i[0] = ( (const int* ALIGNED( 64 ))a00 )[5];
    g.i[0] = ( (const int* ALIGNED( 64 ))a00 )[6];
    h.i[0] = ( (const int* ALIGNED( 64 ))a00 )[7];

    a.i[1] = ( (const int* ALIGNED( 64 ))a01 )[0];
    b.i[1] = ( (const int* ALIGNED( 64 ))a01 )[1];
    c.i[1] = ( (const int* ALIGNED( 64 ))a01 )[2];
    d.i[1] = ( (const int* ALIGNED( 64 ))a01 )[3];
    e.i[1] = ( (const int* ALIGNED( 64 ))a01 )[4];
    f.i[1] = ( (const int* ALIGNED( 64 ))a01 )[5];
    g.i[1] = ( (const int* ALIGNED( 64 ))a01 )[6];
    h.i[1] = ( (const int* ALIGNED( 64 ))a01 )[7];

    a.i[2] = ( (const int* ALIGNED( 64 ))a02 )[0];
    b.i[2] = ( (const int* ALIGNED( 64 ))a02 )[1];
    c.i[2] = ( (const int* ALIGNED( 64 ))a02 )[2];
    d.i[2] = ( (const int* ALIGNED( 64 ))a02 )[3];
    e.i[2] = ( (const int* ALIGNED( 64 ))a02 )[4];
    f.i[2] = ( (const int* ALIGNED( 64 ))a02 )[5];
    g.i[2] = ( (const int* ALIGNED( 64 ))a02 )[6];
    h.i[2] = ( (const int* ALIGNED( 64 ))a02 )[7];

    a.i[3] = ( (const int* ALIGNED( 64 ))a03 )[0];
    b.i[3] = ( (const int* ALIGNED( 64 ))a03 )[1];
    c.i[3] = ( (const int* ALIGNED( 64 ))a03 )[2];
    d.i[3] = ( (const int* ALIGNED( 64 ))a03 )[3];
    e.i[3] = ( (const int* ALIGNED( 64 ))a03 )[4];
    f.i[3] = ( (const int* ALIGNED( 64 ))a03 )[5];
    g.i[3] = ( (const int* ALIGNED( 64 ))a03 )[6];
    h.i[3] = ( (const int* ALIGNED( 64 ))a03 )[7];

    a.i[4] = ( (const int* ALIGNED( 64 ))a04 )[0];
    b.i[4] = ( (const int* ALIGNED( 64 ))a04 )[1];
    c.i[4] = ( (const int* ALIGNED( 64 ))a04 )[2];
    d.i[4] = ( (const int* ALIGNED( 64 ))a04 )[3];
    e.i[4] = ( (const int* ALIGNED( 64 ))a04 )[4];
    f.i[4] = ( (const int* ALIGNED( 64 ))a04 )[5];
    g.i[4] = ( (const int* ALIGNED( 64 ))a04 )[6];
    h.i[4] = ( (const int* ALIGNED( 64 ))a04 )[7];

    a.i[5] = ( (const int* ALIGNED( 64 ))a05 )[0];
    b.i[5] = ( (const int* ALIGNED( 64 ))a05 )[1];
    c.i[5] = ( (const int* ALIGNED( 64 ))a05 )[2];
    d.i[5] = ( (const int* ALIGNED( 64 ))a05 )[3];
    e.i[5] = ( (const int* ALIGNED( 64 ))a05 )[4];
    f.i[5] = ( (const int* ALIGNED( 64 ))a05 )[5];
    g.i[5] = ( (const int* ALIGNED( 64 ))a05 )[6];
    h.i[5] = ( (const int* ALIGNED( 64 ))a05 )[7];

    a.i[6] = ( (const int* ALIGNED( 64 ))a06 )[0];
    b.i[6] = ( (const int* ALIGNED( 64 ))a06 )[1];
    c.i[6] = ( (const int* ALIGNED( 64 ))a06 )[2];
    d.i[6] = ( (const int* ALIGNED( 64 ))a06 )[3];
    e.i[6] = ( (const int* ALIGNED( 64 ))a06 )[4];
    f.i[6] = ( (const int* ALIGNED( 64 ))a06 )[5];
    g.i[6] = ( (const int* ALIGNED( 64 ))a06 )[6];
    h.i[6] = ( (const int* ALIGNED( 64 ))a06 )[7];

    a.i[7] = ( (const int* ALIGNED( 64 ))a07 )[0];
    b.i[7] = ( (const int* ALIGNED( 64 ))a07 )[1];
    c.i[7] = ( (const int* ALIGNED( 64 ))a07 )[2];
    d.i[7] = ( (const int* ALIGNED( 64 ))a07 )[3];
    e.i[7] = ( (const int* ALIGNED( 64 ))a07 )[4];
    f.i[7] = ( (const int* ALIGNED( 64 ))a07 )[5];
    g.i[7] = ( (const int* ALIGNED( 64 ))a07 )[6];
    h.i[7] = ( (const int* ALIGNED( 64 ))a07 )[7];

    a.i[8] = ( (const int* ALIGNED( 64 ))a08 )[0];
    b.i[8] = ( (const int* ALIGNED( 64 ))a08 )[1];
    c.i[8] = ( (const int* ALIGNED( 64 ))a08 )[2];
    d.i[8] = ( (const int* ALIGNED( 64 ))a08 )[3];
    e.i[8] = ( (const int* ALIGNED( 64 ))a08 )[4];
    f.i[8] = ( (const int* ALIGNED( 64 ))a08 )[5];
    g.i[8] = ( (const int* ALIGNED( 64 ))a08 )[6];
    h.i[8] = ( (const int* ALIGNED( 64 ))a08 )[7];

    a.i[9] = ( (const int* ALIGNED( 64 ))a09 )[0];
    b.i[9] = ( (const int* ALIGNED( 64 ))a09 )[1];
    c.i[9] = ( (const int* ALIGNED( 64 ))a09 )[2];
    d.i[9] = ( (const int* ALIGNED( 64 ))a09 )[3];
    e.i[9] = ( (const int* ALIGNED( 64 ))a09 )[4];
    f.i[9] = ( (const int* ALIGNED( 64 ))a09 )[5];
    g.i[9] = ( (const int* ALIGNED( 64 ))a09 )[6];
    h.i[9] = ( (const int* ALIGNED( 64 ))a09 )[7];

    a.i[10] = ( (const int* ALIGNED( 64 ))a10 )[0];
    b.i[10] = ( (const int* ALIGNED( 64 ))a10 )[1];
    c.i[10] = ( (const int* ALIGNED( 64 ))a10 )[2];
    d.i[10] = ( (const int* ALIGNED( 64 ))a10 )[3];
    e.i[10] = ( (const int* ALIGNED( 64 ))a10 )[4];
    f.i[10] = ( (const int* ALIGNED( 64 ))a10 )[5];
    g.i[10] = ( (const int* ALIGNED( 64 ))a10 )[6];
    h.i[10] = ( (const int* ALIGNED( 64 ))a10 )[7];

    a.i[11] = ( (const int* ALIGNED( 64 ))a11 )[0];
    b.i[11] = ( (const int* ALIGNED( 64 ))a11 )[1];
    c.i[11] = ( (const int* ALIGNED( 64 ))a11 )[2];
    d.i[11] = ( (const int* ALIGNED( 64 ))a11 )[3];
    e.i[11] = ( (const int* ALIGNED( 64 ))a11 )[4];
    f.i[11] = ( (const int* ALIGNED( 64 ))a11 )[5];
    g.i[11] = ( (const int* ALIGNED( 64 ))a11 )[6];
    h.i[11] = ( (const int* ALIGNED( 64 ))a11 )[7];

    a.i[12] = ( (const int* ALIGNED( 64 ))a12 )[0];
    b.i[12] = ( (const int* ALIGNED( 64 ))a12 )[1];
    c.i[12] = ( (const int* ALIGNED( 64 ))a12 )[2];
    d.i[12] = ( (const int* ALIGNED( 64 ))a12 )[3];
    e.i[12] = ( (const int* ALIGNED( 64 ))a12 )[4];
    f.i[12] = ( (const int* ALIGNED( 64 ))a12 )[5];
    g.i[12] = ( (const int* ALIGNED( 64 ))a12 )[6];
    h.i[12] = ( (const int* ALIGNED( 64 ))a12 )[7];

    a.i[13] = ( (const int* ALIGNED( 64 ))a13 )[0];
    b.i[13] = ( (const int* ALIGNED( 64 ))a13 )[1];
    c.i[13] = ( (const int* ALIGNED( 64 ))a13 )[2];
    d.i[13] = ( (const int* ALIGNED( 64 ))a13 )[3];
    e.i[13] = ( (const int* ALIGNED( 64 ))a13 )[4];
    f.i[13] = ( (const int* ALIGNED( 64 ))a13 )[5];
    g.i[13] = ( (const int* ALIGNED( 64 ))a13 )[6];
    h.i[13] = ( (const int* ALIGNED( 64 ))a13 )[7];

    a.i[14] = ( (const int* ALIGNED( 64 ))a14 )[0];
    b.i[14] = ( (const int* ALIGNED( 64 ))a14 )[1];
    c.i[14] = ( (const int* ALIGNED( 64 ))a14 )[2];
    d.i[14] = ( (const int* ALIGNED( 64 ))a14 )[3];
    e.i[14] = ( (const int* ALIGNED( 64 ))a14 )[4];
    f.i[14] = ( (const int* ALIGNED( 64 ))a14 )[5];
    g.i[14] = ( (const int* ALIGNED( 64 ))a14 )[6];
    h.i[14] = ( (const int* ALIGNED( 64 ))a14 )[7];

    a.i[15] = ( (const int* ALIGNED( 64 ))a15 )[0];
    b.i[15] = ( (const int* ALIGNED( 64 ))a15 )[1];
    c.i[15] = ( (const int* ALIGNED( 64 ))a15 )[2];
    d.i[15] = ( (const int* ALIGNED( 64 ))a15 )[3];
    e.i[15] = ( (const int* ALIGNED( 64 ))a15 )[4];
    f.i[15] = ( (const int* ALIGNED( 64 ))a15 )[5];
    g.i[15] = ( (const int* ALIGNED( 64 ))a15 )[6];
    h.i[15] = ( (const int* ALIGNED( 64 ))a15 )[7];
}

inline void
load_16x16_tr( const void* ALIGNED( 64 ) a00, const void* ALIGNED( 64 ) a01,
               const void* ALIGNED( 64 ) a02, const void* ALIGNED( 64 ) a03,
               const void* ALIGNED( 64 ) a04, const void* ALIGNED( 64 ) a05,
               const void* ALIGNED( 64 ) a06, const void* ALIGNED( 64 ) a07,
               const void* ALIGNED( 64 ) a08, const void* ALIGNED( 64 ) a09,
               const void* ALIGNED( 64 ) a10, const void* ALIGNED( 64 ) a11,
               const void* ALIGNED( 64 ) a12, const void* ALIGNED( 64 ) a13,
               const void* ALIGNED( 64 ) a14, const void* ALIGNED( 64 ) a15,
               v16& b00, v16& b01, v16& b02, v16& b03, v16& b04, v16& b05,
               v16& b06, v16& b07, v16& b08, v16& b09, v16& b10, v16& b11,
               v16& b12, v16& b13, v16& b14, v16& b15 )
{
    b00.i[0] = ( (const int* ALIGNED( 64 ))a00 )[0];
    b01.i[0] = ( (const int* ALIGNED( 64 ))a00 )[1];
    b02.i[0] = ( (const int* ALIGNED( 64 ))a00 )[2];
    b03.i[0] = ( (const int* ALIGNED( 64 ))a00 )[3];
    b04.i[0] = ( (const int* ALIGNED( 64 ))a00 )[4];
    b05.i[0] = ( (const int* ALIGNED( 64 ))a00 )[5];
    b06.i[0] = ( (const int* ALIGNED( 64 ))a00 )[6];
    b07.i[0] = ( (const int* ALIGNED( 64 ))a00 )[7];
    b08.i[0] = ( (const int* ALIGNED( 64 ))a00 )[8];
    b09.i[0] = ( (const int* ALIGNED( 64 ))a00 )[9];
    b10.i[0] = ( (const int* ALIGNED( 64 ))a00 )[10];
    b11.i[0] = ( (const int* ALIGNED( 64 ))a00 )[11];
    b12.i[0] = ( (const int* ALIGNED( 64 ))a00 )[12];
    b13.i[0] = ( (const int* ALIGNED( 64 ))a00 )[13];
    b14.i[0] = ( (const int* ALIGNED( 64 ))a00 )[14];
    b15.i[0] = ( (const int* ALIGNED( 64 ))a00 )[15];

    b00.i[1] = ( (const int* ALIGNED( 64 ))a01 )[0];
    b01.i[1] = ( (const int* ALIGNED( 64 ))a01 )[1];
    b02.i[1] = ( (const int* ALIGNED( 64 ))a01 )[2];
    b03.i[1] = ( (const int* ALIGNED( 64 ))a01 )[3];
    b04.i[1] = ( (const int* ALIGNED( 64 ))a01 )[4];
    b05.i[1] = ( (const int* ALIGNED( 64 ))a01 )[5];
    b06.i[1] = ( (const int* ALIGNED( 64 ))a01 )[6];
    b07.i[1] = ( (const int* ALIGNED( 64 ))a01 )[7];
    b08.i[1] = ( (const int* ALIGNED( 64 ))a01 )[8];
    b09.i[1] = ( (const int* ALIGNED( 64 ))a01 )[9];
    b10.i[1] = ( (const int* ALIGNED( 64 ))a01 )[10];
    b11.i[1] = ( (const int* ALIGNED( 64 ))a01 )[11];
    b12.i[1] = ( (const int* ALIGNED( 64 ))a01 )[12];
    b13.i[1] = ( (const int* ALIGNED( 64 ))a01 )[13];
    b14.i[1] = ( (const int* ALIGNED( 64 ))a01 )[14];
    b15.i[1] = ( (const int* ALIGNED( 64 ))a01 )[15];

    b00.i[2] = ( (const int* ALIGNED( 64 ))a02 )[0];
    b01.i[2] = ( (const int* ALIGNED( 64 ))a02 )[1];
    b02.i[2] = ( (const int* ALIGNED( 64 ))a02 )[2];
    b03.i[2] = ( (const int* ALIGNED( 64 ))a02 )[3];
    b04.i[2] = ( (const int* ALIGNED( 64 ))a02 )[4];
    b05.i[2] = ( (const int* ALIGNED( 64 ))a02 )[5];
    b06.i[2] = ( (const int* ALIGNED( 64 ))a02 )[6];
    b07.i[2] = ( (const int* ALIGNED( 64 ))a02 )[7];
    b08.i[2] = ( (const int* ALIGNED( 64 ))a02 )[8];
    b09.i[2] = ( (const int* ALIGNED( 64 ))a02 )[9];
    b10.i[2] = ( (const int* ALIGNED( 64 ))a02 )[10];
    b11.i[2] = ( (const int* ALIGNED( 64 ))a02 )[11];
    b12.i[2] = ( (const int* ALIGNED( 64 ))a02 )[12];
    b13.i[2] = ( (const int* ALIGNED( 64 ))a02 )[13];
    b14.i[2] = ( (const int* ALIGNED( 64 ))a02 )[14];
    b15.i[2] = ( (const int* ALIGNED( 64 ))a02 )[15];

    b00.i[3] = ( (const int* ALIGNED( 64 ))a03 )[0];
    b01.i[3] = ( (const int* ALIGNED( 64 ))a03 )[1];
    b02.i[3] = ( (const int* ALIGNED( 64 ))a03 )[2];
    b03.i[3] = ( (const int* ALIGNED( 64 ))a03 )[3];
    b04.i[3] = ( (const int* ALIGNED( 64 ))a03 )[4];
    b05.i[3] = ( (const int* ALIGNED( 64 ))a03 )[5];
    b06.i[3] = ( (const int* ALIGNED( 64 ))a03 )[6];
    b07.i[3] = ( (const int* ALIGNED( 64 ))a03 )[7];
    b08.i[3] = ( (const int* ALIGNED( 64 ))a03 )[8];
    b09.i[3] = ( (const int* ALIGNED( 64 ))a03 )[9];
    b10.i[3] = ( (const int* ALIGNED( 64 ))a03 )[10];
    b11.i[3] = ( (const int* ALIGNED( 64 ))a03 )[11];
    b12.i[3] = ( (const int* ALIGNED( 64 ))a03 )[12];
    b13.i[3] = ( (const int* ALIGNED( 64 ))a03 )[13];
    b14.i[3] = ( (const int* ALIGNED( 64 ))a03 )[14];
    b15.i[3] = ( (const int* ALIGNED( 64 ))a03 )[15];

    b00.i[4] = ( (const int* ALIGNED( 64 ))a04 )[0];
    b01.i[4] = ( (const int* ALIGNED( 64 ))a04 )[1];
    b02.i[4] = ( (const int* ALIGNED( 64 ))a04 )[2];
    b03.i[4] = ( (const int* ALIGNED( 64 ))a04 )[3];
    b04.i[4] = ( (const int* ALIGNED( 64 ))a04 )[4];
    b05.i[4] = ( (const int* ALIGNED( 64 ))a04 )[5];
    b06.i[4] = ( (const int* ALIGNED( 64 ))a04 )[6];
    b07.i[4] = ( (const int* ALIGNED( 64 ))a04 )[7];
    b08.i[4] = ( (const int* ALIGNED( 64 ))a04 )[8];
    b09.i[4] = ( (const int* ALIGNED( 64 ))a04 )[9];
    b10.i[4] = ( (const int* ALIGNED( 64 ))a04 )[10];
    b11.i[4] = ( (const int* ALIGNED( 64 ))a04 )[11];
    b12.i[4] = ( (const int* ALIGNED( 64 ))a04 )[12];
    b13.i[4] = ( (const int* ALIGNED( 64 ))a04 )[13];
    b14.i[4] = ( (const int* ALIGNED( 64 ))a04 )[14];
    b15.i[4] = ( (const int* ALIGNED( 64 ))a04 )[15];

    b00.i[5] = ( (const int* ALIGNED( 64 ))a05 )[0];
    b01.i[5] = ( (const int* ALIGNED( 64 ))a05 )[1];
    b02.i[5] = ( (const int* ALIGNED( 64 ))a05 )[2];
    b03.i[5] = ( (const int* ALIGNED( 64 ))a05 )[3];
    b04.i[5] = ( (const int* ALIGNED( 64 ))a05 )[4];
    b05.i[5] = ( (const int* ALIGNED( 64 ))a05 )[5];
    b06.i[5] = ( (const int* ALIGNED( 64 ))a05 )[6];
    b07.i[5] = ( (const int* ALIGNED( 64 ))a05 )[7];
    b08.i[5] = ( (const int* ALIGNED( 64 ))a05 )[8];
    b09.i[5] = ( (const int* ALIGNED( 64 ))a05 )[9];
    b10.i[5] = ( (const int* ALIGNED( 64 ))a05 )[10];
    b11.i[5] = ( (const int* ALIGNED( 64 ))a05 )[11];
    b12.i[5] = ( (const int* ALIGNED( 64 ))a05 )[12];
    b13.i[5] = ( (const int* ALIGNED( 64 ))a05 )[13];
    b14.i[5] = ( (const int* ALIGNED( 64 ))a05 )[14];
    b15.i[5] = ( (const int* ALIGNED( 64 ))a05 )[15];

    b00.i[6] = ( (const int* ALIGNED( 64 ))a06 )[0];
    b01.i[6] = ( (const int* ALIGNED( 64 ))a06 )[1];
    b02.i[6] = ( (const int* ALIGNED( 64 ))a06 )[2];
    b03.i[6] = ( (const int* ALIGNED( 64 ))a06 )[3];
    b04.i[6] = ( (const int* ALIGNED( 64 ))a06 )[4];
    b05.i[6] = ( (const int* ALIGNED( 64 ))a06 )[5];
    b06.i[6] = ( (const int* ALIGNED( 64 ))a06 )[6];
    b07.i[6] = ( (const int* ALIGNED( 64 ))a06 )[7];
    b08.i[6] = ( (const int* ALIGNED( 64 ))a06 )[8];
    b09.i[6] = ( (const int* ALIGNED( 64 ))a06 )[9];
    b10.i[6] = ( (const int* ALIGNED( 64 ))a06 )[10];
    b11.i[6] = ( (const int* ALIGNED( 64 ))a06 )[11];
    b12.i[6] = ( (const int* ALIGNED( 64 ))a06 )[12];
    b13.i[6] = ( (const int* ALIGNED( 64 ))a06 )[13];
    b14.i[6] = ( (const int* ALIGNED( 64 ))a06 )[14];
    b15.i[6] = ( (const int* ALIGNED( 64 ))a06 )[15];

    b00.i[7] = ( (const int* ALIGNED( 64 ))a07 )[0];
    b01.i[7] = ( (const int* ALIGNED( 64 ))a07 )[1];
    b02.i[7] = ( (const int* ALIGNED( 64 ))a07 )[2];
    b03.i[7] = ( (const int* ALIGNED( 64 ))a07 )[3];
    b04.i[7] = ( (const int* ALIGNED( 64 ))a07 )[4];
    b05.i[7] = ( (const int* ALIGNED( 64 ))a07 )[5];
    b06.i[7] = ( (const int* ALIGNED( 64 ))a07 )[6];
    b07.i[7] = ( (const int* ALIGNED( 64 ))a07 )[7];
    b08.i[7] = ( (const int* ALIGNED( 64 ))a07 )[8];
    b09.i[7] = ( (const int* ALIGNED( 64 ))a07 )[9];
    b10.i[7] = ( (const int* ALIGNED( 64 ))a07 )[10];
    b11.i[7] = ( (const int* ALIGNED( 64 ))a07 )[11];
    b12.i[7] = ( (const int* ALIGNED( 64 ))a07 )[12];
    b13.i[7] = ( (const int* ALIGNED( 64 ))a07 )[13];
    b14.i[7] = ( (const int* ALIGNED( 64 ))a07 )[14];
    b15.i[7] = ( (const int* ALIGNED( 64 ))a07 )[15];

    b00.i[8] = ( (const int* ALIGNED( 64 ))a08 )[0];
    b01.i[8] = ( (const int* ALIGNED( 64 ))a08 )[1];
    b02.i[8] = ( (const int* ALIGNED( 64 ))a08 )[2];
    b03.i[8] = ( (const int* ALIGNED( 64 ))a08 )[3];
    b04.i[8] = ( (const int* ALIGNED( 64 ))a08 )[4];
    b05.i[8] = ( (const int* ALIGNED( 64 ))a08 )[5];
    b06.i[8] = ( (const int* ALIGNED( 64 ))a08 )[6];
    b07.i[8] = ( (const int* ALIGNED( 64 ))a08 )[7];
    b08.i[8] = ( (const int* ALIGNED( 64 ))a08 )[8];
    b09.i[8] = ( (const int* ALIGNED( 64 ))a08 )[9];
    b10.i[8] = ( (const int* ALIGNED( 64 ))a08 )[10];
    b11.i[8] = ( (const int* ALIGNED( 64 ))a08 )[11];
    b12.i[8] = ( (const int* ALIGNED( 64 ))a08 )[12];
    b13.i[8] = ( (const int* ALIGNED( 64 ))a08 )[13];
    b14.i[8] = ( (const int* ALIGNED( 64 ))a08 )[14];
    b15.i[8] = ( (const int* ALIGNED( 64 ))a08 )[15];

    b00.i[9] = ( (const int* ALIGNED( 64 ))a09 )[0];
    b01.i[9] = ( (const int* ALIGNED( 64 ))a09 )[1];
    b02.i[9] = ( (const int* ALIGNED( 64 ))a09 )[2];
    b03.i[9] = ( (const int* ALIGNED( 64 ))a09 )[3];
    b04.i[9] = ( (const int* ALIGNED( 64 ))a09 )[4];
    b05.i[9] = ( (const int* ALIGNED( 64 ))a09 )[5];
    b06.i[9] = ( (const int* ALIGNED( 64 ))a09 )[6];
    b07.i[9] = ( (const int* ALIGNED( 64 ))a09 )[7];
    b08.i[9] = ( (const int* ALIGNED( 64 ))a09 )[8];
    b09.i[9] = ( (const int* ALIGNED( 64 ))a09 )[9];
    b10.i[9] = ( (const int* ALIGNED( 64 ))a09 )[10];
    b11.i[9] = ( (const int* ALIGNED( 64 ))a09 )[11];
    b12.i[9] = ( (const int* ALIGNED( 64 ))a09 )[12];
    b13.i[9] = ( (const int* ALIGNED( 64 ))a09 )[13];
    b14.i[9] = ( (const int* ALIGNED( 64 ))a09 )[14];
    b15.i[9] = ( (const int* ALIGNED( 64 ))a09 )[15];

    b00.i[10] = ( (const int* ALIGNED( 64 ))a10 )[0];
    b01.i[10] = ( (const int* ALIGNED( 64 ))a10 )[1];
    b02.i[10] = ( (const int* ALIGNED( 64 ))a10 )[2];
    b03.i[10] = ( (const int* ALIGNED( 64 ))a10 )[3];
    b04.i[10] = ( (const int* ALIGNED( 64 ))a10 )[4];
    b05.i[10] = ( (const int* ALIGNED( 64 ))a10 )[5];
    b06.i[10] = ( (const int* ALIGNED( 64 ))a10 )[6];
    b07.i[10] = ( (const int* ALIGNED( 64 ))a10 )[7];
    b08.i[10] = ( (const int* ALIGNED( 64 ))a10 )[8];
    b09.i[10] = ( (const int* ALIGNED( 64 ))a10 )[9];
    b10.i[10] = ( (const int* ALIGNED( 64 ))a10 )[10];
    b11.i[10] = ( (const int* ALIGNED( 64 ))a10 )[11];
    b12.i[10] = ( (const int* ALIGNED( 64 ))a10 )[12];
    b13.i[10] = ( (const int* ALIGNED( 64 ))a10 )[13];
    b14.i[10] = ( (const int* ALIGNED( 64 ))a10 )[14];
    b15.i[10] = ( (const int* ALIGNED( 64 ))a10 )[15];

    b00.i[11] = ( (const int* ALIGNED( 64 ))a11 )[0];
    b01.i[11] = ( (const int* ALIGNED( 64 ))a11 )[1];
    b02.i[11] = ( (const int* ALIGNED( 64 ))a11 )[2];
    b03.i[11] = ( (const int* ALIGNED( 64 ))a11 )[3];
    b04.i[11] = ( (const int* ALIGNED( 64 ))a11 )[4];
    b05.i[11] = ( (const int* ALIGNED( 64 ))a11 )[5];
    b06.i[11] = ( (const int* ALIGNED( 64 ))a11 )[6];
    b07.i[11] = ( (const int* ALIGNED( 64 ))a11 )[7];
    b08.i[11] = ( (const int* ALIGNED( 64 ))a11 )[8];
    b09.i[11] = ( (const int* ALIGNED( 64 ))a11 )[9];
    b10.i[11] = ( (const int* ALIGNED( 64 ))a11 )[10];
    b11.i[11] = ( (const int* ALIGNED( 64 ))a11 )[11];
    b12.i[11] = ( (const int* ALIGNED( 64 ))a11 )[12];
    b13.i[11] = ( (const int* ALIGNED( 64 ))a11 )[13];
    b14.i[11] = ( (const int* ALIGNED( 64 ))a11 )[14];
    b15.i[11] = ( (const int* ALIGNED( 64 ))a11 )[15];

    b00.i[12] = ( (const int* ALIGNED( 64 ))a12 )[0];
    b01.i[12] = ( (const int* ALIGNED( 64 ))a12 )[1];
    b02.i[12] = ( (const int* ALIGNED( 64 ))a12 )[2];
    b03.i[12] = ( (const int* ALIGNED( 64 ))a12 )[3];
    b04.i[12] = ( (const int* ALIGNED( 64 ))a12 )[4];
    b05.i[12] = ( (const int* ALIGNED( 64 ))a12 )[5];
    b06.i[12] = ( (const int* ALIGNED( 64 ))a12 )[6];
    b07.i[12] = ( (const int* ALIGNED( 64 ))a12 )[7];
    b08.i[12] = ( (const int* ALIGNED( 64 ))a12 )[8];
    b09.i[12] = ( (const int* ALIGNED( 64 ))a12 )[9];
    b10.i[12] = ( (const int* ALIGNED( 64 ))a12 )[10];
    b11.i[12] = ( (const int* ALIGNED( 64 ))a12 )[11];
    b12.i[12] = ( (const int* ALIGNED( 64 ))a12 )[12];
    b13.i[12] = ( (const int* ALIGNED( 64 ))a12 )[13];
    b14.i[12] = ( (const int* ALIGNED( 64 ))a12 )[14];
    b15.i[12] = ( (const int* ALIGNED( 64 ))a12 )[15];

    b00.i[13] = ( (const int* ALIGNED( 64 ))a13 )[0];
    b01.i[13] = ( (const int* ALIGNED( 64 ))a13 )[1];
    b02.i[13] = ( (const int* ALIGNED( 64 ))a13 )[2];
    b03.i[13] = ( (const int* ALIGNED( 64 ))a13 )[3];
    b04.i[13] = ( (const int* ALIGNED( 64 ))a13 )[4];
    b05.i[13] = ( (const int* ALIGNED( 64 ))a13 )[5];
    b06.i[13] = ( (const int* ALIGNED( 64 ))a13 )[6];
    b07.i[13] = ( (const int* ALIGNED( 64 ))a13 )[7];
    b08.i[13] = ( (const int* ALIGNED( 64 ))a13 )[8];
    b09.i[13] = ( (const int* ALIGNED( 64 ))a13 )[9];
    b10.i[13] = ( (const int* ALIGNED( 64 ))a13 )[10];
    b11.i[13] = ( (const int* ALIGNED( 64 ))a13 )[11];
    b12.i[13] = ( (const int* ALIGNED( 64 ))a13 )[12];
    b13.i[13] = ( (const int* ALIGNED( 64 ))a13 )[13];
    b14.i[13] = ( (const int* ALIGNED( 64 ))a13 )[14];
    b15.i[13] = ( (const int* ALIGNED( 64 ))a13 )[15];

    b00.i[14] = ( (const int* ALIGNED( 64 ))a14 )[0];
    b01.i[14] = ( (const int* ALIGNED( 64 ))a14 )[1];
    b02.i[14] = ( (const int* ALIGNED( 64 ))a14 )[2];
    b03.i[14] = ( (const int* ALIGNED( 64 ))a14 )[3];
    b04.i[14] = ( (const int* ALIGNED( 64 ))a14 )[4];
    b05.i[14] = ( (const int* ALIGNED( 64 ))a14 )[5];
    b06.i[14] = ( (const int* ALIGNED( 64 ))a14 )[6];
    b07.i[14] = ( (const int* ALIGNED( 64 ))a14 )[7];
    b08.i[14] = ( (const int* ALIGNED( 64 ))a14 )[8];
    b09.i[14] = ( (const int* ALIGNED( 64 ))a14 )[9];
    b10.i[14] = ( (const int* ALIGNED( 64 ))a14 )[10];
    b11.i[14] = ( (const int* ALIGNED( 64 ))a14 )[11];
    b12.i[14] = ( (const int* ALIGNED( 64 ))a14 )[12];
    b13.i[14] = ( (const int* ALIGNED( 64 ))a14 )[13];
    b14.i[14] = ( (const int* ALIGNED( 64 ))a14 )[14];
    b15.i[14] = ( (const int* ALIGNED( 64 ))a14 )[15];

    b00.i[15] = ( (const int* ALIGNED( 64 ))a15 )[0];
    b01.i[15] = ( (const int* ALIGNED( 64 ))a15 )[1];
    b02.i[15] = ( (const int* ALIGNED( 64 ))a15 )[2];
    b03.i[15] = ( (const int* ALIGNED( 64 ))a15 )[3];
    b04.i[15] = ( (const int* ALIGNED( 64 ))a15 )[4];
    b05.i[15] = ( (const int* ALIGNED( 64 ))a15 )[5];
    b06.i[15] = ( (const int* ALIGNED( 64 ))a15 )[6];
    b07.i[15] = ( (const int* ALIGNED( 64 ))a15 )[7];
    b08.i[15] = ( (const int* ALIGNED( 64 ))a15 )[8];
    b09.i[15] = ( (const int* ALIGNED( 64 ))a15 )[9];
    b10.i[15] = ( (const int* ALIGNED( 64 ))a15 )[10];
    b11.i[15] = ( (const int* ALIGNED( 64 ))a15 )[11];
    b12.i[15] = ( (const int* ALIGNED( 64 ))a15 )[12];
    b13.i[15] = ( (const int* ALIGNED( 64 ))a15 )[13];
    b14.i[15] = ( (const int* ALIGNED( 64 ))a15 )[14];
    b15.i[15] = ( (const int* ALIGNED( 64 ))a15 )[15];
}

inline void
load_16x8_tr_p( const void* ALIGNED( 64 ) a00, const void* ALIGNED( 64 ) a01,
                const void* ALIGNED( 64 ) a02, const void* ALIGNED( 64 ) a03,
                const void* ALIGNED( 64 ) a04, const void* ALIGNED( 64 ) a05,
                const void* ALIGNED( 64 ) a06, const void* ALIGNED( 64 ) a07,
                v16& b00, v16& b01, v16& b02, v16& b03, v16& b04, v16& b05,
                v16& b06, v16& b07 )
{
    b00.i[0] = ( (const int* ALIGNED( 64 ))a00 )[0];
    b01.i[0] = ( (const int* ALIGNED( 64 ))a00 )[1];
    b02.i[0] = ( (const int* ALIGNED( 64 ))a00 )[2];
    b03.i[0] = ( (const int* ALIGNED( 64 ))a00 )[3];
    b04.i[0] = ( (const int* ALIGNED( 64 ))a00 )[4];
    b05.i[0] = ( (const int* ALIGNED( 64 ))a00 )[5];
    b06.i[0] = ( (const int* ALIGNED( 64 ))a00 )[6];
    b07.i[0] = ( (const int* ALIGNED( 64 ))a00 )[7];
    b00.i[1] = ( (const int* ALIGNED( 64 ))a00 )[8];
    b01.i[1] = ( (const int* ALIGNED( 64 ))a00 )[9];
    b02.i[1] = ( (const int* ALIGNED( 64 ))a00 )[10];
    b03.i[1] = ( (const int* ALIGNED( 64 ))a00 )[11];
    b04.i[1] = ( (const int* ALIGNED( 64 ))a00 )[12];
    b05.i[1] = ( (const int* ALIGNED( 64 ))a00 )[13];
    b06.i[1] = ( (const int* ALIGNED( 64 ))a00 )[14];
    b07.i[1] = ( (const int* ALIGNED( 64 ))a00 )[15];

    b00.i[2] = ( (const int* ALIGNED( 64 ))a01 )[0];
    b01.i[2] = ( (const int* ALIGNED( 64 ))a01 )[1];
    b02.i[2] = ( (const int* ALIGNED( 64 ))a01 )[2];
    b03.i[2] = ( (const int* ALIGNED( 64 ))a01 )[3];
    b04.i[2] = ( (const int* ALIGNED( 64 ))a01 )[4];
    b05.i[2] = ( (const int* ALIGNED( 64 ))a01 )[5];
    b06.i[2] = ( (const int* ALIGNED( 64 ))a01 )[6];
    b07.i[2] = ( (const int* ALIGNED( 64 ))a01 )[7];
    b00.i[3] = ( (const int* ALIGNED( 64 ))a01 )[8];
    b01.i[3] = ( (const int* ALIGNED( 64 ))a01 )[9];
    b02.i[3] = ( (const int* ALIGNED( 64 ))a01 )[10];
    b03.i[3] = ( (const int* ALIGNED( 64 ))a01 )[11];
    b04.i[3] = ( (const int* ALIGNED( 64 ))a01 )[12];
    b05.i[3] = ( (const int* ALIGNED( 64 ))a01 )[13];
    b06.i[3] = ( (const int* ALIGNED( 64 ))a01 )[14];
    b07.i[3] = ( (const int* ALIGNED( 64 ))a01 )[15];

    b00.i[4] = ( (const int* ALIGNED( 64 ))a02 )[0];
    b01.i[4] = ( (const int* ALIGNED( 64 ))a02 )[1];
    b02.i[4] = ( (const int* ALIGNED( 64 ))a02 )[2];
    b03.i[4] = ( (const int* ALIGNED( 64 ))a02 )[3];
    b04.i[4] = ( (const int* ALIGNED( 64 ))a02 )[4];
    b05.i[4] = ( (const int* ALIGNED( 64 ))a02 )[5];
    b06.i[4] = ( (const int* ALIGNED( 64 ))a02 )[6];
    b07.i[4] = ( (const int* ALIGNED( 64 ))a02 )[7];
    b00.i[5] = ( (const int* ALIGNED( 64 ))a02 )[8];
    b01.i[5] = ( (const int* ALIGNED( 64 ))a02 )[9];
    b02.i[5] = ( (const int* ALIGNED( 64 ))a02 )[10];
    b03.i[5] = ( (const int* ALIGNED( 64 ))a02 )[11];
    b04.i[5] = ( (const int* ALIGNED( 64 ))a02 )[12];
    b05.i[5] = ( (const int* ALIGNED( 64 ))a02 )[13];
    b06.i[5] = ( (const int* ALIGNED( 64 ))a02 )[14];
    b07.i[5] = ( (const int* ALIGNED( 64 ))a02 )[15];

    b00.i[6] = ( (const int* ALIGNED( 64 ))a03 )[0];
    b01.i[6] = ( (const int* ALIGNED( 64 ))a03 )[1];
    b02.i[6] = ( (const int* ALIGNED( 64 ))a03 )[2];
    b03.i[6] = ( (const int* ALIGNED( 64 ))a03 )[3];
    b04.i[6] = ( (const int* ALIGNED( 64 ))a03 )[4];
    b05.i[6] = ( (const int* ALIGNED( 64 ))a03 )[5];
    b06.i[6] = ( (const int* ALIGNED( 64 ))a03 )[6];
    b07.i[6] = ( (const int* ALIGNED( 64 ))a03 )[7];
    b00.i[7] = ( (const int* ALIGNED( 64 ))a03 )[8];
    b01.i[7] = ( (const int* ALIGNED( 64 ))a03 )[9];
    b02.i[7] = ( (const int* ALIGNED( 64 ))a03 )[10];
    b03.i[7] = ( (const int* ALIGNED( 64 ))a03 )[11];
    b04.i[7] = ( (const int* ALIGNED( 64 ))a03 )[12];
    b05.i[7] = ( (const int* ALIGNED( 64 ))a03 )[13];
    b06.i[7] = ( (const int* ALIGNED( 64 ))a03 )[14];
    b07.i[7] = ( (const int* ALIGNED( 64 ))a03 )[15];

    b00.i[8] = ( (const int* ALIGNED( 64 ))a04 )[0];
    b01.i[8] = ( (const int* ALIGNED( 64 ))a04 )[1];
    b02.i[8] = ( (const int* ALIGNED( 64 ))a04 )[2];
    b03.i[8] = ( (const int* ALIGNED( 64 ))a04 )[3];
    b04.i[8] = ( (const int* ALIGNED( 64 ))a04 )[4];
    b05.i[8] = ( (const int* ALIGNED( 64 ))a04 )[5];
    b06.i[8] = ( (const int* ALIGNED( 64 ))a04 )[6];
    b07.i[8] = ( (const int* ALIGNED( 64 ))a04 )[7];
    b00.i[9] = ( (const int* ALIGNED( 64 ))a04 )[8];
    b01.i[9] = ( (const int* ALIGNED( 64 ))a04 )[9];
    b02.i[9] = ( (const int* ALIGNED( 64 ))a04 )[10];
    b03.i[9] = ( (const int* ALIGNED( 64 ))a04 )[11];
    b04.i[9] = ( (const int* ALIGNED( 64 ))a04 )[12];
    b05.i[9] = ( (const int* ALIGNED( 64 ))a04 )[13];
    b06.i[9] = ( (const int* ALIGNED( 64 ))a04 )[14];
    b07.i[9] = ( (const int* ALIGNED( 64 ))a04 )[15];

    b00.i[10] = ( (const int* ALIGNED( 64 ))a05 )[0];
    b01.i[10] = ( (const int* ALIGNED( 64 ))a05 )[1];
    b02.i[10] = ( (const int* ALIGNED( 64 ))a05 )[2];
    b03.i[10] = ( (const int* ALIGNED( 64 ))a05 )[3];
    b04.i[10] = ( (const int* ALIGNED( 64 ))a05 )[4];
    b05.i[10] = ( (const int* ALIGNED( 64 ))a05 )[5];
    b06.i[10] = ( (const int* ALIGNED( 64 ))a05 )[6];
    b07.i[10] = ( (const int* ALIGNED( 64 ))a05 )[7];
    b00.i[11] = ( (const int* ALIGNED( 64 ))a05 )[8];
    b01.i[11] = ( (const int* ALIGNED( 64 ))a05 )[9];
    b02.i[11] = ( (const int* ALIGNED( 64 ))a05 )[10];
    b03.i[11] = ( (const int* ALIGNED( 64 ))a05 )[11];
    b04.i[11] = ( (const int* ALIGNED( 64 ))a05 )[12];
    b05.i[11] = ( (const int* ALIGNED( 64 ))a05 )[13];
    b06.i[11] = ( (const int* ALIGNED( 64 ))a05 )[14];
    b07.i[11] = ( (const int* ALIGNED( 64 ))a05 )[15];

    b00.i[12] = ( (const int* ALIGNED( 64 ))a06 )[0];
    b01.i[12] = ( (const int* ALIGNED( 64 ))a06 )[1];
    b02.i[12] = ( (const int* ALIGNED( 64 ))a06 )[2];
    b03.i[12] = ( (const int* ALIGNED( 64 ))a06 )[3];
    b04.i[12] = ( (const int* ALIGNED( 64 ))a06 )[4];
    b05.i[12] = ( (const int* ALIGNED( 64 ))a06 )[5];
    b06.i[12] = ( (const int* ALIGNED( 64 ))a06 )[6];
    b07.i[12] = ( (const int* ALIGNED( 64 ))a06 )[7];
    b00.i[13] = ( (const int* ALIGNED( 64 ))a06 )[8];
    b01.i[13] = ( (const int* ALIGNED( 64 ))a06 )[9];
    b02.i[13] = ( (const int* ALIGNED( 64 ))a06 )[10];
    b03.i[13] = ( (const int* ALIGNED( 64 ))a06 )[11];
    b04.i[13] = ( (const int* ALIGNED( 64 ))a06 )[12];
    b05.i[13] = ( (const int* ALIGNED( 64 ))a06 )[13];
    b06.i[13] = ( (const int* ALIGNED( 64 ))a06 )[14];
    b07.i[13] = ( (const int* ALIGNED( 64 ))a06 )[15];

    b00.i[14] = ( (const int* ALIGNED( 64 ))a07 )[0];
    b01.i[14] = ( (const int* ALIGNED( 64 ))a07 )[1];
    b02.i[14] = ( (const int* ALIGNED( 64 ))a07 )[2];
    b03.i[14] = ( (const int* ALIGNED( 64 ))a07 )[3];
    b04.i[14] = ( (const int* ALIGNED( 64 ))a07 )[4];
    b05.i[14] = ( (const int* ALIGNED( 64 ))a07 )[5];
    b06.i[14] = ( (const int* ALIGNED( 64 ))a07 )[6];
    b07.i[14] = ( (const int* ALIGNED( 64 ))a07 )[7];
    b00.i[15] = ( (const int* ALIGNED( 64 ))a07 )[8];
    b01.i[15] = ( (const int* ALIGNED( 64 ))a07 )[9];
    b02.i[15] = ( (const int* ALIGNED( 64 ))a07 )[10];
    b03.i[15] = ( (const int* ALIGNED( 64 ))a07 )[11];
    b04.i[15] = ( (const int* ALIGNED( 64 ))a07 )[12];
    b05.i[15] = ( (const int* ALIGNED( 64 ))a07 )[13];
    b06.i[15] = ( (const int* ALIGNED( 64 ))a07 )[14];
    b07.i[15] = ( (const int* ALIGNED( 64 ))a07 )[15];
}

inline void
load_16x16_tr_p( const void* ALIGNED( 64 ) a00, const void* ALIGNED( 64 ) a01,
                 const void* ALIGNED( 64 ) a02, const void* ALIGNED( 64 ) a03,
                 const void* ALIGNED( 64 ) a04, const void* ALIGNED( 64 ) a05,
                 const void* ALIGNED( 64 ) a06, const void* ALIGNED( 64 ) a07,
                 const void* ALIGNED( 64 ) a08, const void* ALIGNED( 64 ) a09,
                 const void* ALIGNED( 64 ) a10, const void* ALIGNED( 64 ) a11,
                 const void* ALIGNED( 64 ) a12, const void* ALIGNED( 64 ) a13,
                 const void* ALIGNED( 64 ) a14, const void* ALIGNED( 64 ) a15,
                 v16& b00, v16& b01, v16& b02, v16& b03, v16& b04, v16& b05,
                 v16& b06, v16& b07, v16& b08, v16& b09, v16& b10, v16& b11,
                 v16& b12, v16& b13, v16& b14, v16& b15 )
{
    b00.i[0] = ( (const int* ALIGNED( 64 ))a00 )[0];
    b01.i[0] = ( (const int* ALIGNED( 64 ))a00 )[1];
    b02.i[0] = ( (const int* ALIGNED( 64 ))a00 )[2];
    b03.i[0] = ( (const int* ALIGNED( 64 ))a00 )[3];
    b04.i[0] = ( (const int* ALIGNED( 64 ))a00 )[4];
    b05.i[0] = ( (const int* ALIGNED( 64 ))a00 )[5];
    b06.i[0] = ( (const int* ALIGNED( 64 ))a00 )[6];
    b07.i[0] = ( (const int* ALIGNED( 64 ))a00 )[7];
    b00.i[1] = ( (const int* ALIGNED( 64 ))a00 )[8];
    b01.i[1] = ( (const int* ALIGNED( 64 ))a00 )[9];
    b02.i[1] = ( (const int* ALIGNED( 64 ))a00 )[10];
    b03.i[1] = ( (const int* ALIGNED( 64 ))a00 )[11];
    b04.i[1] = ( (const int* ALIGNED( 64 ))a00 )[12];
    b05.i[1] = ( (const int* ALIGNED( 64 ))a00 )[13];
    b06.i[1] = ( (const int* ALIGNED( 64 ))a00 )[14];
    b07.i[1] = ( (const int* ALIGNED( 64 ))a00 )[15];

    b00.i[2] = ( (const int* ALIGNED( 64 ))a01 )[0];
    b01.i[2] = ( (const int* ALIGNED( 64 ))a01 )[1];
    b02.i[2] = ( (const int* ALIGNED( 64 ))a01 )[2];
    b03.i[2] = ( (const int* ALIGNED( 64 ))a01 )[3];
    b04.i[2] = ( (const int* ALIGNED( 64 ))a01 )[4];
    b05.i[2] = ( (const int* ALIGNED( 64 ))a01 )[5];
    b06.i[2] = ( (const int* ALIGNED( 64 ))a01 )[6];
    b07.i[2] = ( (const int* ALIGNED( 64 ))a01 )[7];
    b00.i[3] = ( (const int* ALIGNED( 64 ))a01 )[8];
    b01.i[3] = ( (const int* ALIGNED( 64 ))a01 )[9];
    b02.i[3] = ( (const int* ALIGNED( 64 ))a01 )[10];
    b03.i[3] = ( (const int* ALIGNED( 64 ))a01 )[11];
    b04.i[3] = ( (const int* ALIGNED( 64 ))a01 )[12];
    b05.i[3] = ( (const int* ALIGNED( 64 ))a01 )[13];
    b06.i[3] = ( (const int* ALIGNED( 64 ))a01 )[14];
    b07.i[3] = ( (const int* ALIGNED( 64 ))a01 )[15];

    b00.i[4] = ( (const int* ALIGNED( 64 ))a02 )[0];
    b01.i[4] = ( (const int* ALIGNED( 64 ))a02 )[1];
    b02.i[4] = ( (const int* ALIGNED( 64 ))a02 )[2];
    b03.i[4] = ( (const int* ALIGNED( 64 ))a02 )[3];
    b04.i[4] = ( (const int* ALIGNED( 64 ))a02 )[4];
    b05.i[4] = ( (const int* ALIGNED( 64 ))a02 )[5];
    b06.i[4] = ( (const int* ALIGNED( 64 ))a02 )[6];
    b07.i[4] = ( (const int* ALIGNED( 64 ))a02 )[7];
    b00.i[5] = ( (const int* ALIGNED( 64 ))a02 )[8];
    b01.i[5] = ( (const int* ALIGNED( 64 ))a02 )[9];
    b02.i[5] = ( (const int* ALIGNED( 64 ))a02 )[10];
    b03.i[5] = ( (const int* ALIGNED( 64 ))a02 )[11];
    b04.i[5] = ( (const int* ALIGNED( 64 ))a02 )[12];
    b05.i[5] = ( (const int* ALIGNED( 64 ))a02 )[13];
    b06.i[5] = ( (const int* ALIGNED( 64 ))a02 )[14];
    b07.i[5] = ( (const int* ALIGNED( 64 ))a02 )[15];

    b00.i[6] = ( (const int* ALIGNED( 64 ))a03 )[0];
    b01.i[6] = ( (const int* ALIGNED( 64 ))a03 )[1];
    b02.i[6] = ( (const int* ALIGNED( 64 ))a03 )[2];
    b03.i[6] = ( (const int* ALIGNED( 64 ))a03 )[3];
    b04.i[6] = ( (const int* ALIGNED( 64 ))a03 )[4];
    b05.i[6] = ( (const int* ALIGNED( 64 ))a03 )[5];
    b06.i[6] = ( (const int* ALIGNED( 64 ))a03 )[6];
    b07.i[6] = ( (const int* ALIGNED( 64 ))a03 )[7];
    b00.i[7] = ( (const int* ALIGNED( 64 ))a03 )[8];
    b01.i[7] = ( (const int* ALIGNED( 64 ))a03 )[9];
    b02.i[7] = ( (const int* ALIGNED( 64 ))a03 )[10];
    b03.i[7] = ( (const int* ALIGNED( 64 ))a03 )[11];
    b04.i[7] = ( (const int* ALIGNED( 64 ))a03 )[12];
    b05.i[7] = ( (const int* ALIGNED( 64 ))a03 )[13];
    b06.i[7] = ( (const int* ALIGNED( 64 ))a03 )[14];
    b07.i[7] = ( (const int* ALIGNED( 64 ))a03 )[15];

    b00.i[8] = ( (const int* ALIGNED( 64 ))a04 )[0];
    b01.i[8] = ( (const int* ALIGNED( 64 ))a04 )[1];
    b02.i[8] = ( (const int* ALIGNED( 64 ))a04 )[2];
    b03.i[8] = ( (const int* ALIGNED( 64 ))a04 )[3];
    b04.i[8] = ( (const int* ALIGNED( 64 ))a04 )[4];
    b05.i[8] = ( (const int* ALIGNED( 64 ))a04 )[5];
    b06.i[8] = ( (const int* ALIGNED( 64 ))a04 )[6];
    b07.i[8] = ( (const int* ALIGNED( 64 ))a04 )[7];
    b00.i[9] = ( (const int* ALIGNED( 64 ))a04 )[8];
    b01.i[9] = ( (const int* ALIGNED( 64 ))a04 )[9];
    b02.i[9] = ( (const int* ALIGNED( 64 ))a04 )[10];
    b03.i[9] = ( (const int* ALIGNED( 64 ))a04 )[11];
    b04.i[9] = ( (const int* ALIGNED( 64 ))a04 )[12];
    b05.i[9] = ( (const int* ALIGNED( 64 ))a04 )[13];
    b06.i[9] = ( (const int* ALIGNED( 64 ))a04 )[14];
    b07.i[9] = ( (const int* ALIGNED( 64 ))a04 )[15];

    b00.i[10] = ( (const int* ALIGNED( 64 ))a05 )[0];
    b01.i[10] = ( (const int* ALIGNED( 64 ))a05 )[1];
    b02.i[10] = ( (const int* ALIGNED( 64 ))a05 )[2];
    b03.i[10] = ( (const int* ALIGNED( 64 ))a05 )[3];
    b04.i[10] = ( (const int* ALIGNED( 64 ))a05 )[4];
    b05.i[10] = ( (const int* ALIGNED( 64 ))a05 )[5];
    b06.i[10] = ( (const int* ALIGNED( 64 ))a05 )[6];
    b07.i[10] = ( (const int* ALIGNED( 64 ))a05 )[7];
    b00.i[11] = ( (const int* ALIGNED( 64 ))a05 )[8];
    b01.i[11] = ( (const int* ALIGNED( 64 ))a05 )[9];
    b02.i[11] = ( (const int* ALIGNED( 64 ))a05 )[10];
    b03.i[11] = ( (const int* ALIGNED( 64 ))a05 )[11];
    b04.i[11] = ( (const int* ALIGNED( 64 ))a05 )[12];
    b05.i[11] = ( (const int* ALIGNED( 64 ))a05 )[13];
    b06.i[11] = ( (const int* ALIGNED( 64 ))a05 )[14];
    b07.i[11] = ( (const int* ALIGNED( 64 ))a05 )[15];

    b00.i[12] = ( (const int* ALIGNED( 64 ))a06 )[0];
    b01.i[12] = ( (const int* ALIGNED( 64 ))a06 )[1];
    b02.i[12] = ( (const int* ALIGNED( 64 ))a06 )[2];
    b03.i[12] = ( (const int* ALIGNED( 64 ))a06 )[3];
    b04.i[12] = ( (const int* ALIGNED( 64 ))a06 )[4];
    b05.i[12] = ( (const int* ALIGNED( 64 ))a06 )[5];
    b06.i[12] = ( (const int* ALIGNED( 64 ))a06 )[6];
    b07.i[12] = ( (const int* ALIGNED( 64 ))a06 )[7];
    b00.i[13] = ( (const int* ALIGNED( 64 ))a06 )[8];
    b01.i[13] = ( (const int* ALIGNED( 64 ))a06 )[9];
    b02.i[13] = ( (const int* ALIGNED( 64 ))a06 )[10];
    b03.i[13] = ( (const int* ALIGNED( 64 ))a06 )[11];
    b04.i[13] = ( (const int* ALIGNED( 64 ))a06 )[12];
    b05.i[13] = ( (const int* ALIGNED( 64 ))a06 )[13];
    b06.i[13] = ( (const int* ALIGNED( 64 ))a06 )[14];
    b07.i[13] = ( (const int* ALIGNED( 64 ))a06 )[15];

    b00.i[14] = ( (const int* ALIGNED( 64 ))a07 )[0];
    b01.i[14] = ( (const int* ALIGNED( 64 ))a07 )[1];
    b02.i[14] = ( (const int* ALIGNED( 64 ))a07 )[2];
    b03.i[14] = ( (const int* ALIGNED( 64 ))a07 )[3];
    b04.i[14] = ( (const int* ALIGNED( 64 ))a07 )[4];
    b05.i[14] = ( (const int* ALIGNED( 64 ))a07 )[5];
    b06.i[14] = ( (const int* ALIGNED( 64 ))a07 )[6];
    b07.i[14] = ( (const int* ALIGNED( 64 ))a07 )[7];
    b00.i[15] = ( (const int* ALIGNED( 64 ))a07 )[8];
    b01.i[15] = ( (const int* ALIGNED( 64 ))a07 )[9];
    b02.i[15] = ( (const int* ALIGNED( 64 ))a07 )[10];
    b03.i[15] = ( (const int* ALIGNED( 64 ))a07 )[11];
    b04.i[15] = ( (const int* ALIGNED( 64 ))a07 )[12];
    b05.i[15] = ( (const int* ALIGNED( 64 ))a07 )[13];
    b06.i[15] = ( (const int* ALIGNED( 64 ))a07 )[14];
    b07.i[15] = ( (const int* ALIGNED( 64 ))a07 )[15];

    b08.i[0] = ( (const int* ALIGNED( 64 ))a08 )[0];
    b09.i[0] = ( (const int* ALIGNED( 64 ))a08 )[1];
    b10.i[0] = ( (const int* ALIGNED( 64 ))a08 )[2];
    b11.i[0] = ( (const int* ALIGNED( 64 ))a08 )[3];
    b12.i[0] = ( (const int* ALIGNED( 64 ))a08 )[4];
    b13.i[0] = ( (const int* ALIGNED( 64 ))a08 )[5];
    b14.i[0] = ( (const int* ALIGNED( 64 ))a08 )[6];
    b15.i[0] = ( (const int* ALIGNED( 64 ))a08 )[7];
    b08.i[1] = ( (const int* ALIGNED( 64 ))a08 )[8];
    b09.i[1] = ( (const int* ALIGNED( 64 ))a08 )[9];
    b10.i[1] = ( (const int* ALIGNED( 64 ))a08 )[10];
    b11.i[1] = ( (const int* ALIGNED( 64 ))a08 )[11];
    b12.i[1] = ( (const int* ALIGNED( 64 ))a08 )[12];
    b13.i[1] = ( (const int* ALIGNED( 64 ))a08 )[13];
    b14.i[1] = ( (const int* ALIGNED( 64 ))a08 )[14];
    b15.i[1] = ( (const int* ALIGNED( 64 ))a08 )[15];

    b08.i[2] = ( (const int* ALIGNED( 64 ))a09 )[0];
    b09.i[2] = ( (const int* ALIGNED( 64 ))a09 )[1];
    b10.i[2] = ( (const int* ALIGNED( 64 ))a09 )[2];
    b11.i[2] = ( (const int* ALIGNED( 64 ))a09 )[3];
    b12.i[2] = ( (const int* ALIGNED( 64 ))a09 )[4];
    b13.i[2] = ( (const int* ALIGNED( 64 ))a09 )[5];
    b14.i[2] = ( (const int* ALIGNED( 64 ))a09 )[6];
    b15.i[2] = ( (const int* ALIGNED( 64 ))a09 )[7];
    b08.i[3] = ( (const int* ALIGNED( 64 ))a09 )[8];
    b09.i[3] = ( (const int* ALIGNED( 64 ))a09 )[9];
    b10.i[3] = ( (const int* ALIGNED( 64 ))a09 )[10];
    b11.i[3] = ( (const int* ALIGNED( 64 ))a09 )[11];
    b12.i[3] = ( (const int* ALIGNED( 64 ))a09 )[12];
    b13.i[3] = ( (const int* ALIGNED( 64 ))a09 )[13];
    b14.i[3] = ( (const int* ALIGNED( 64 ))a09 )[14];
    b15.i[3] = ( (const int* ALIGNED( 64 ))a09 )[15];

    b08.i[4] = ( (const int* ALIGNED( 64 ))a10 )[0];
    b09.i[4] = ( (const int* ALIGNED( 64 ))a10 )[1];
    b10.i[4] = ( (const int* ALIGNED( 64 ))a10 )[2];
    b11.i[4] = ( (const int* ALIGNED( 64 ))a10 )[3];
    b12.i[4] = ( (const int* ALIGNED( 64 ))a10 )[4];
    b13.i[4] = ( (const int* ALIGNED( 64 ))a10 )[5];
    b14.i[4] = ( (const int* ALIGNED( 64 ))a10 )[6];
    b15.i[4] = ( (const int* ALIGNED( 64 ))a10 )[7];
    b08.i[5] = ( (const int* ALIGNED( 64 ))a10 )[8];
    b09.i[5] = ( (const int* ALIGNED( 64 ))a10 )[9];
    b10.i[5] = ( (const int* ALIGNED( 64 ))a10 )[10];
    b11.i[5] = ( (const int* ALIGNED( 64 ))a10 )[11];
    b12.i[5] = ( (const int* ALIGNED( 64 ))a10 )[12];
    b13.i[5] = ( (const int* ALIGNED( 64 ))a10 )[13];
    b14.i[5] = ( (const int* ALIGNED( 64 ))a10 )[14];
    b15.i[5] = ( (const int* ALIGNED( 64 ))a10 )[15];

    b08.i[6] = ( (const int* ALIGNED( 64 ))a11 )[0];
    b09.i[6] = ( (const int* ALIGNED( 64 ))a11 )[1];
    b10.i[6] = ( (const int* ALIGNED( 64 ))a11 )[2];
    b11.i[6] = ( (const int* ALIGNED( 64 ))a11 )[3];
    b12.i[6] = ( (const int* ALIGNED( 64 ))a11 )[4];
    b13.i[6] = ( (const int* ALIGNED( 64 ))a11 )[5];
    b14.i[6] = ( (const int* ALIGNED( 64 ))a11 )[6];
    b15.i[6] = ( (const int* ALIGNED( 64 ))a11 )[7];
    b08.i[7] = ( (const int* ALIGNED( 64 ))a11 )[8];
    b09.i[7] = ( (const int* ALIGNED( 64 ))a11 )[9];
    b10.i[7] = ( (const int* ALIGNED( 64 ))a11 )[10];
    b11.i[7] = ( (const int* ALIGNED( 64 ))a11 )[11];
    b12.i[7] = ( (const int* ALIGNED( 64 ))a11 )[12];
    b13.i[7] = ( (const int* ALIGNED( 64 ))a11 )[13];
    b14.i[7] = ( (const int* ALIGNED( 64 ))a11 )[14];
    b15.i[7] = ( (const int* ALIGNED( 64 ))a11 )[15];

    b08.i[8] = ( (const int* ALIGNED( 64 ))a12 )[0];
    b09.i[8] = ( (const int* ALIGNED( 64 ))a12 )[1];
    b10.i[8] = ( (const int* ALIGNED( 64 ))a12 )[2];
    b11.i[8] = ( (const int* ALIGNED( 64 ))a12 )[3];
    b12.i[8] = ( (const int* ALIGNED( 64 ))a12 )[4];
    b13.i[8] = ( (const int* ALIGNED( 64 ))a12 )[5];
    b14.i[8] = ( (const int* ALIGNED( 64 ))a12 )[6];
    b15.i[8] = ( (const int* ALIGNED( 64 ))a12 )[7];
    b08.i[9] = ( (const int* ALIGNED( 64 ))a12 )[8];
    b09.i[9] = ( (const int* ALIGNED( 64 ))a12 )[9];
    b10.i[9] = ( (const int* ALIGNED( 64 ))a12 )[10];
    b11.i[9] = ( (const int* ALIGNED( 64 ))a12 )[11];
    b12.i[9] = ( (const int* ALIGNED( 64 ))a12 )[12];
    b13.i[9] = ( (const int* ALIGNED( 64 ))a12 )[13];
    b14.i[9] = ( (const int* ALIGNED( 64 ))a12 )[14];
    b15.i[9] = ( (const int* ALIGNED( 64 ))a12 )[15];

    b08.i[10] = ( (const int* ALIGNED( 64 ))a13 )[0];
    b09.i[10] = ( (const int* ALIGNED( 64 ))a13 )[1];
    b10.i[10] = ( (const int* ALIGNED( 64 ))a13 )[2];
    b11.i[10] = ( (const int* ALIGNED( 64 ))a13 )[3];
    b12.i[10] = ( (const int* ALIGNED( 64 ))a13 )[4];
    b13.i[10] = ( (const int* ALIGNED( 64 ))a13 )[5];
    b14.i[10] = ( (const int* ALIGNED( 64 ))a13 )[6];
    b15.i[10] = ( (const int* ALIGNED( 64 ))a13 )[7];
    b08.i[11] = ( (const int* ALIGNED( 64 ))a13 )[8];
    b09.i[11] = ( (const int* ALIGNED( 64 ))a13 )[9];
    b10.i[11] = ( (const int* ALIGNED( 64 ))a13 )[10];
    b11.i[11] = ( (const int* ALIGNED( 64 ))a13 )[11];
    b12.i[11] = ( (const int* ALIGNED( 64 ))a13 )[12];
    b13.i[11] = ( (const int* ALIGNED( 64 ))a13 )[13];
    b14.i[11] = ( (const int* ALIGNED( 64 ))a13 )[14];
    b15.i[11] = ( (const int* ALIGNED( 64 ))a13 )[15];

    b08.i[12] = ( (const int* ALIGNED( 64 ))a14 )[0];
    b09.i[12] = ( (const int* ALIGNED( 64 ))a14 )[1];
    b10.i[12] = ( (const int* ALIGNED( 64 ))a14 )[2];
    b11.i[12] = ( (const int* ALIGNED( 64 ))a14 )[3];
    b12.i[12] = ( (const int* ALIGNED( 64 ))a14 )[4];
    b13.i[12] = ( (const int* ALIGNED( 64 ))a14 )[5];
    b14.i[12] = ( (const int* ALIGNED( 64 ))a14 )[6];
    b15.i[12] = ( (const int* ALIGNED( 64 ))a14 )[7];
    b08.i[13] = ( (const int* ALIGNED( 64 ))a14 )[8];
    b09.i[13] = ( (const int* ALIGNED( 64 ))a14 )[9];
    b10.i[13] = ( (const int* ALIGNED( 64 ))a14 )[10];
    b11.i[13] = ( (const int* ALIGNED( 64 ))a14 )[11];
    b12.i[13] = ( (const int* ALIGNED( 64 ))a14 )[12];
    b13.i[13] = ( (const int* ALIGNED( 64 ))a14 )[13];
    b14.i[13] = ( (const int* ALIGNED( 64 ))a14 )[14];
    b15.i[13] = ( (const int* ALIGNED( 64 ))a14 )[15];

    b08.i[14] = ( (const int* ALIGNED( 64 ))a15 )[0];
    b09.i[14] = ( (const int* ALIGNED( 64 ))a15 )[1];
    b10.i[14] = ( (const int* ALIGNED( 64 ))a15 )[2];
    b11.i[14] = ( (const int* ALIGNED( 64 ))a15 )[3];
    b12.i[14] = ( (const int* ALIGNED( 64 ))a15 )[4];
    b13.i[14] = ( (const int* ALIGNED( 64 ))a15 )[5];
    b14.i[14] = ( (const int* ALIGNED( 64 ))a15 )[6];
    b15.i[14] = ( (const int* ALIGNED( 64 ))a15 )[7];
    b08.i[15] = ( (const int* ALIGNED( 64 ))a15 )[8];
    b09.i[15] = ( (const int* ALIGNED( 64 ))a15 )[9];
    b10.i[15] = ( (const int* ALIGNED( 64 ))a15 )[10];
    b11.i[15] = ( (const int* ALIGNED( 64 ))a15 )[11];
    b12.i[15] = ( (const int* ALIGNED( 64 ))a15 )[12];
    b13.i[15] = ( (const int* ALIGNED( 64 ))a15 )[13];
    b14.i[15] = ( (const int* ALIGNED( 64 ))a15 )[14];
    b15.i[15] = ( (const int* ALIGNED( 64 ))a15 )[15];
}

inline void store_16x1_tr( const v16& a, void* a00, void* a01, void* a02,
                           void* a03, void* a04, void* a05, void* a06,
                           void* a07, void* a08, void* a09, void* a10,
                           void* a11, void* a12, void* a13, void* a14,
                           void* a15 )
{
    ( (int*)a00 )[0] = a.i[0];
    ( (int*)a01 )[0] = a.i[1];
    ( (int*)a02 )[0] = a.i[2];
    ( (int*)a03 )[0] = a.i[3];
    ( (int*)a04 )[0] = a.i[4];
    ( (int*)a05 )[0] = a.i[5];
    ( (int*)a06 )[0] = a.i[6];
    ( (int*)a07 )[0] = a.i[7];
    ( (int*)a08 )[0] = a.i[8];
    ( (int*)a09 )[0] = a.i[9];
    ( (int*)a10 )[0] = a.i[10];
    ( (int*)a11 )[0] = a.i[11];
    ( (int*)a12 )[0] = a.i[12];
    ( (int*)a13 )[0] = a.i[13];
    ( (int*)a14 )[0] = a.i[14];
    ( (int*)a15 )[0] = a.i[15];
}

inline void store_16x2_tr( const v16& a, const v16& b, void* ALIGNED( 8 ) a00,
                           void* ALIGNED( 8 ) a01, void* ALIGNED( 8 ) a02,
                           void* ALIGNED( 8 ) a03, void* ALIGNED( 8 ) a04,
                           void* ALIGNED( 8 ) a05, void* ALIGNED( 8 ) a06,
                           void* ALIGNED( 8 ) a07, void* ALIGNED( 8 ) a08,
                           void* ALIGNED( 8 ) a09, void* ALIGNED( 8 ) a10,
                           void* ALIGNED( 8 ) a11, void* ALIGNED( 8 ) a12,
                           void* ALIGNED( 8 ) a13, void* ALIGNED( 8 ) a14,
                           void* ALIGNED( 8 ) a15 )
{
    ( (int* ALIGNED( 8 ))a00 )[0] = a.i[0];
    ( (int* ALIGNED( 8 ))a00 )[1] = b.i[0];

    ( (int* ALIGNED( 8 ))a01 )[0] = a.i[1];
    ( (int* ALIGNED( 8 ))a01 )[1] = b.i[1];

    ( (int* ALIGNED( 8 ))a02 )[0] = a.i[2];
    ( (int* ALIGNED( 8 ))a02 )[1] = b.i[2];

    ( (int* ALIGNED( 8 ))a03 )[0] = a.i[3];
    ( (int* ALIGNED( 8 ))a03 )[1] = b.i[3];

    ( (int* ALIGNED( 8 ))a04 )[0] = a.i[4];
    ( (int* ALIGNED( 8 ))a04 )[1] = b.i[4];

    ( (int* ALIGNED( 8 ))a05 )[0] = a.i[5];
    ( (int* ALIGNED( 8 ))a05 )[1] = b.i[5];

    ( (int* ALIGNED( 8 ))a06 )[0] = a.i[6];
    ( (int* ALIGNED( 8 ))a06 )[1] = b.i[6];

    ( (int* ALIGNED( 8 ))a07 )[0] = a.i[7];
    ( (int* ALIGNED( 8 ))a07 )[1] = b.i[7];

    ( (int* ALIGNED( 8 ))a08 )[0] = a.i[8];
    ( (int* ALIGNED( 8 ))a08 )[1] = b.i[8];

    ( (int* ALIGNED( 8 ))a09 )[0] = a.i[9];
    ( (int* ALIGNED( 8 ))a09 )[1] = b.i[9];

    ( (int* ALIGNED( 8 ))a10 )[0] = a.i[10];
    ( (int* ALIGNED( 8 ))a10 )[1] = b.i[10];

    ( (int* ALIGNED( 8 ))a11 )[0] = a.i[11];
    ( (int* ALIGNED( 8 ))a11 )[1] = b.i[11];

    ( (int* ALIGNED( 8 ))a12 )[0] = a.i[12];
    ( (int* ALIGNED( 8 ))a12 )[1] = b.i[12];

    ( (int* ALIGNED( 8 ))a13 )[0] = a.i[13];
    ( (int* ALIGNED( 8 ))a13 )[1] = b.i[13];

    ( (int* ALIGNED( 8 ))a14 )[0] = a.i[14];
    ( (int* ALIGNED( 8 ))a14 )[1] = b.i[14];

    ( (int* ALIGNED( 8 ))a15 )[0] = a.i[15];
    ( (int* ALIGNED( 8 ))a15 )[1] = b.i[15];
}

inline void store_16x3_tr( const v16& a, const v16& b, const v16& c,
                           void* ALIGNED( 64 ) a00, void* ALIGNED( 64 ) a01,
                           void* ALIGNED( 64 ) a02, void* ALIGNED( 64 ) a03,
                           void* ALIGNED( 64 ) a04, void* ALIGNED( 64 ) a05,
                           void* ALIGNED( 64 ) a06, void* ALIGNED( 64 ) a07,
                           void* ALIGNED( 64 ) a08, void* ALIGNED( 64 ) a09,
                           void* ALIGNED( 64 ) a10, void* ALIGNED( 64 ) a11,
                           void* ALIGNED( 64 ) a12, void* ALIGNED( 64 ) a13,
                           void* ALIGNED( 64 ) a14, void* ALIGNED( 64 ) a15 )
{
    ( (int* ALIGNED( 64 ))a00 )[0] = a.i[0];
    ( (int* ALIGNED( 64 ))a00 )[1] = b.i[0];
    ( (int* ALIGNED( 64 ))a00 )[2] = c.i[0];

    ( (int* ALIGNED( 64 ))a01 )[0] = a.i[1];
    ( (int* ALIGNED( 64 ))a01 )[1] = b.i[1];
    ( (int* ALIGNED( 64 ))a01 )[2] = c.i[1];

    ( (int* ALIGNED( 64 ))a02 )[0] = a.i[2];
    ( (int* ALIGNED( 64 ))a02 )[1] = b.i[2];
    ( (int* ALIGNED( 64 ))a02 )[2] = c.i[2];

    ( (int* ALIGNED( 64 ))a03 )[0] = a.i[3];
    ( (int* ALIGNED( 64 ))a03 )[1] = b.i[3];
    ( (int* ALIGNED( 64 ))a03 )[2] = c.i[3];

    ( (int* ALIGNED( 64 ))a04 )[0] = a.i[4];
    ( (int* ALIGNED( 64 ))a04 )[1] = b.i[4];
    ( (int* ALIGNED( 64 ))a04 )[2] = c.i[4];

    ( (int* ALIGNED( 64 ))a05 )[0] = a.i[5];
    ( (int* ALIGNED( 64 ))a05 )[1] = b.i[5];
    ( (int* ALIGNED( 64 ))a05 )[2] = c.i[5];

    ( (int* ALIGNED( 64 ))a06 )[0] = a.i[6];
    ( (int* ALIGNED( 64 ))a06 )[1] = b.i[6];
    ( (int* ALIGNED( 64 ))a06 )[2] = c.i[6];

    ( (int* ALIGNED( 64 ))a07 )[0] = a.i[7];
    ( (int* ALIGNED( 64 ))a07 )[1] = b.i[7];
    ( (int* ALIGNED( 64 ))a07 )[2] = c.i[7];

    ( (int* ALIGNED( 64 ))a08 )[0] = a.i[8];
    ( (int* ALIGNED( 64 ))a08 )[1] = b.i[8];
    ( (int* ALIGNED( 64 ))a08 )[2] = c.i[8];

    ( (int* ALIGNED( 64 ))a09 )[0] = a.i[9];
    ( (int* ALIGNED( 64 ))a09 )[1] = b.i[9];
    ( (int* ALIGNED( 64 ))a09 )[2] = c.i[9];

    ( (int* ALIGNED( 64 ))a10 )[0] = a.i[10];
    ( (int* ALIGNED( 64 ))a10 )[1] = b.i[10];
    ( (int* ALIGNED( 64 ))a10 )[2] = c.i[10];

    ( (int* ALIGNED( 64 ))a11 )[0] = a.i[11];
    ( (int* ALIGNED( 64 ))a11 )[1] = b.i[11];
    ( (int* ALIGNED( 64 ))a11 )[2] = c.i[11];

    ( (int* ALIGNED( 64 ))a12 )[0] = a.i[12];
    ( (int* ALIGNED( 64 ))a12 )[1] = b.i[12];
    ( (int* ALIGNED( 64 ))a12 )[2] = c.i[12];

    ( (int* ALIGNED( 64 ))a13 )[0] = a.i[13];
    ( (int* ALIGNED( 64 ))a13 )[1] = b.i[13];
    ( (int* ALIGNED( 64 ))a13 )[2] = c.i[13];

    ( (int* ALIGNED( 64 ))a14 )[0] = a.i[14];
    ( (int* ALIGNED( 64 ))a14 )[1] = b.i[14];
    ( (int* ALIGNED( 64 ))a14 )[2] = c.i[14];

    ( (int* ALIGNED( 64 ))a15 )[0] = a.i[15];
    ( (int* ALIGNED( 64 ))a15 )[1] = b.i[15];
    ( (int* ALIGNED( 64 ))a15 )[2] = c.i[15];
}

inline void store_16x4_tr( const v16& a, const v16& b, const v16& c,
                           const v16& d, void* ALIGNED( 64 ) a00,
                           void* ALIGNED( 64 ) a01, void* ALIGNED( 64 ) a02,
                           void* ALIGNED( 64 ) a03, void* ALIGNED( 64 ) a04,
                           void* ALIGNED( 64 ) a05, void* ALIGNED( 64 ) a06,
                           void* ALIGNED( 64 ) a07, void* ALIGNED( 64 ) a08,
                           void* ALIGNED( 64 ) a09, void* ALIGNED( 64 ) a10,
                           void* ALIGNED( 64 ) a11, void* ALIGNED( 64 ) a12,
                           void* ALIGNED( 64 ) a13, void* ALIGNED( 64 ) a14,
                           void* ALIGNED( 64 ) a15 )
{
    ( (int* ALIGNED( 64 ))a00 )[0] = a.i[0];
    ( (int* ALIGNED( 64 ))a00 )[1] = b.i[0];
    ( (int* ALIGNED( 64 ))a00 )[2] = c.i[0];
    ( (int* ALIGNED( 64 ))a00 )[3] = d.i[0];

    ( (int* ALIGNED( 64 ))a01 )[0] = a.i[1];
    ( (int* ALIGNED( 64 ))a01 )[1] = b.i[1];
    ( (int* ALIGNED( 64 ))a01 )[2] = c.i[1];
    ( (int* ALIGNED( 64 ))a01 )[3] = d.i[1];

    ( (int* ALIGNED( 64 ))a02 )[0] = a.i[2];
    ( (int* ALIGNED( 64 ))a02 )[1] = b.i[2];
    ( (int* ALIGNED( 64 ))a02 )[2] = c.i[2];
    ( (int* ALIGNED( 64 ))a02 )[3] = d.i[2];

    ( (int* ALIGNED( 64 ))a03 )[0] = a.i[3];
    ( (int* ALIGNED( 64 ))a03 )[1] = b.i[3];
    ( (int* ALIGNED( 64 ))a03 )[2] = c.i[3];
    ( (int* ALIGNED( 64 ))a03 )[3] = d.i[3];

    ( (int* ALIGNED( 64 ))a04 )[0] = a.i[4];
    ( (int* ALIGNED( 64 ))a04 )[1] = b.i[4];
    ( (int* ALIGNED( 64 ))a04 )[2] = c.i[4];
    ( (int* ALIGNED( 64 ))a04 )[3] = d.i[4];

    ( (int* ALIGNED( 64 ))a05 )[0] = a.i[5];
    ( (int* ALIGNED( 64 ))a05 )[1] = b.i[5];
    ( (int* ALIGNED( 64 ))a05 )[2] = c.i[5];
    ( (int* ALIGNED( 64 ))a05 )[3] = d.i[5];

    ( (int* ALIGNED( 64 ))a06 )[0] = a.i[6];
    ( (int* ALIGNED( 64 ))a06 )[1] = b.i[6];
    ( (int* ALIGNED( 64 ))a06 )[2] = c.i[6];
    ( (int* ALIGNED( 64 ))a06 )[3] = d.i[6];

    ( (int* ALIGNED( 64 ))a07 )[0] = a.i[7];
    ( (int* ALIGNED( 64 ))a07 )[1] = b.i[7];
    ( (int* ALIGNED( 64 ))a07 )[2] = c.i[7];
    ( (int* ALIGNED( 64 ))a07 )[3] = d.i[7];

    ( (int* ALIGNED( 64 ))a08 )[0] = a.i[8];
    ( (int* ALIGNED( 64 ))a08 )[1] = b.i[8];
    ( (int* ALIGNED( 64 ))a08 )[2] = c.i[8];
    ( (int* ALIGNED( 64 ))a08 )[3] = d.i[8];

    ( (int* ALIGNED( 64 ))a09 )[0] = a.i[9];
    ( (int* ALIGNED( 64 ))a09 )[1] = b.i[9];
    ( (int* ALIGNED( 64 ))a09 )[2] = c.i[9];
    ( (int* ALIGNED( 64 ))a09 )[3] = d.i[9];

    ( (int* ALIGNED( 64 ))a10 )[0] = a.i[10];
    ( (int* ALIGNED( 64 ))a10 )[1] = b.i[10];
    ( (int* ALIGNED( 64 ))a10 )[2] = c.i[10];
    ( (int* ALIGNED( 64 ))a10 )[3] = d.i[10];

    ( (int* ALIGNED( 64 ))a11 )[0] = a.i[11];
    ( (int* ALIGNED( 64 ))a11 )[1] = b.i[11];
    ( (int* ALIGNED( 64 ))a11 )[2] = c.i[11];
    ( (int* ALIGNED( 64 ))a11 )[3] = d.i[11];

    ( (int* ALIGNED( 64 ))a12 )[0] = a.i[12];
    ( (int* ALIGNED( 64 ))a12 )[1] = b.i[12];
    ( (int* ALIGNED( 64 ))a12 )[2] = c.i[12];
    ( (int* ALIGNED( 64 ))a12 )[3] = d.i[12];

    ( (int* ALIGNED( 64 ))a13 )[0] = a.i[13];
    ( (int* ALIGNED( 64 ))a13 )[1] = b.i[13];
    ( (int* ALIGNED( 64 ))a13 )[2] = c.i[13];
    ( (int* ALIGNED( 64 ))a13 )[3] = d.i[13];

    ( (int* ALIGNED( 64 ))a14 )[0] = a.i[14];
    ( (int* ALIGNED( 64 ))a14 )[1] = b.i[14];
    ( (int* ALIGNED( 64 ))a14 )[2] = c.i[14];
    ( (int* ALIGNED( 64 ))a14 )[3] = d.i[14];

    ( (int* ALIGNED( 64 ))a15 )[0] = a.i[15];
    ( (int* ALIGNED( 64 ))a15 )[1] = b.i[15];
    ( (int* ALIGNED( 64 ))a15 )[2] = c.i[15];
    ( (int* ALIGNED( 64 ))a15 )[3] = d.i[15];
}

inline void store_16x8_tr(
    const v16& a, const v16& b, const v16& c, const v16& d, const v16& e,
    const v16& f, const v16& g, const v16& h, void* ALIGNED( 64 ) a00,
    void* ALIGNED( 64 ) a01, void* ALIGNED( 64 ) a02, void* ALIGNED( 64 ) a03,
    void* ALIGNED( 64 ) a04, void* ALIGNED( 64 ) a05, void* ALIGNED( 64 ) a06,
    void* ALIGNED( 64 ) a07, void* ALIGNED( 64 ) a08, void* ALIGNED( 64 ) a09,
    void* ALIGNED( 64 ) a10, void* ALIGNED( 64 ) a11, void* ALIGNED( 64 ) a12,
    void* ALIGNED( 64 ) a13, void* ALIGNED( 64 ) a14, void* ALIGNED( 64 ) a15 )
{
    ( (int* ALIGNED( 64 ))a00 )[0] = a.i[0];
    ( (int* ALIGNED( 64 ))a00 )[1] = b.i[0];
    ( (int* ALIGNED( 64 ))a00 )[2] = c.i[0];
    ( (int* ALIGNED( 64 ))a00 )[3] = d.i[0];
    ( (int* ALIGNED( 64 ))a00 )[4] = e.i[0];
    ( (int* ALIGNED( 64 ))a00 )[5] = f.i[0];
    ( (int* ALIGNED( 64 ))a00 )[6] = g.i[0];
    ( (int* ALIGNED( 64 ))a00 )[7] = h.i[0];

    ( (int* ALIGNED( 64 ))a01 )[0] = a.i[1];
    ( (int* ALIGNED( 64 ))a01 )[1] = b.i[1];
    ( (int* ALIGNED( 64 ))a01 )[2] = c.i[1];
    ( (int* ALIGNED( 64 ))a01 )[3] = d.i[1];
    ( (int* ALIGNED( 64 ))a01 )[4] = e.i[1];
    ( (int* ALIGNED( 64 ))a01 )[5] = f.i[1];
    ( (int* ALIGNED( 64 ))a01 )[6] = g.i[1];
    ( (int* ALIGNED( 64 ))a01 )[7] = h.i[1];

    ( (int* ALIGNED( 64 ))a02 )[0] = a.i[2];
    ( (int* ALIGNED( 64 ))a02 )[1] = b.i[2];
    ( (int* ALIGNED( 64 ))a02 )[2] = c.i[2];
    ( (int* ALIGNED( 64 ))a02 )[3] = d.i[2];
    ( (int* ALIGNED( 64 ))a02 )[4] = e.i[2];
    ( (int* ALIGNED( 64 ))a02 )[5] = f.i[2];
    ( (int* ALIGNED( 64 ))a02 )[6] = g.i[2];
    ( (int* ALIGNED( 64 ))a02 )[7] = h.i[2];

    ( (int* ALIGNED( 64 ))a03 )[0] = a.i[3];
    ( (int* ALIGNED( 64 ))a03 )[1] = b.i[3];
    ( (int* ALIGNED( 64 ))a03 )[2] = c.i[3];
    ( (int* ALIGNED( 64 ))a03 )[3] = d.i[3];
    ( (int* ALIGNED( 64 ))a03 )[4] = e.i[3];
    ( (int* ALIGNED( 64 ))a03 )[5] = f.i[3];
    ( (int* ALIGNED( 64 ))a03 )[6] = g.i[3];
    ( (int* ALIGNED( 64 ))a03 )[7] = h.i[3];

    ( (int* ALIGNED( 64 ))a04 )[0] = a.i[4];
    ( (int* ALIGNED( 64 ))a04 )[1] = b.i[4];
    ( (int* ALIGNED( 64 ))a04 )[2] = c.i[4];
    ( (int* ALIGNED( 64 ))a04 )[3] = d.i[4];
    ( (int* ALIGNED( 64 ))a04 )[4] = e.i[4];
    ( (int* ALIGNED( 64 ))a04 )[5] = f.i[4];
    ( (int* ALIGNED( 64 ))a04 )[6] = g.i[4];
    ( (int* ALIGNED( 64 ))a04 )[7] = h.i[4];

    ( (int* ALIGNED( 64 ))a05 )[0] = a.i[5];
    ( (int* ALIGNED( 64 ))a05 )[1] = b.i[5];
    ( (int* ALIGNED( 64 ))a05 )[2] = c.i[5];
    ( (int* ALIGNED( 64 ))a05 )[3] = d.i[5];
    ( (int* ALIGNED( 64 ))a05 )[4] = e.i[5];
    ( (int* ALIGNED( 64 ))a05 )[5] = f.i[5];
    ( (int* ALIGNED( 64 ))a05 )[6] = g.i[5];
    ( (int* ALIGNED( 64 ))a05 )[7] = h.i[5];

    ( (int* ALIGNED( 64 ))a06 )[0] = a.i[6];
    ( (int* ALIGNED( 64 ))a06 )[1] = b.i[6];
    ( (int* ALIGNED( 64 ))a06 )[2] = c.i[6];
    ( (int* ALIGNED( 64 ))a06 )[3] = d.i[6];
    ( (int* ALIGNED( 64 ))a06 )[4] = e.i[6];
    ( (int* ALIGNED( 64 ))a06 )[5] = f.i[6];
    ( (int* ALIGNED( 64 ))a06 )[6] = g.i[6];
    ( (int* ALIGNED( 64 ))a06 )[7] = h.i[6];

    ( (int* ALIGNED( 64 ))a07 )[0] = a.i[7];
    ( (int* ALIGNED( 64 ))a07 )[1] = b.i[7];
    ( (int* ALIGNED( 64 ))a07 )[2] = c.i[7];
    ( (int* ALIGNED( 64 ))a07 )[3] = d.i[7];
    ( (int* ALIGNED( 64 ))a07 )[4] = e.i[7];
    ( (int* ALIGNED( 64 ))a07 )[5] = f.i[7];
    ( (int* ALIGNED( 64 ))a07 )[6] = g.i[7];
    ( (int* ALIGNED( 64 ))a07 )[7] = h.i[7];

    ( (int* ALIGNED( 64 ))a08 )[0] = a.i[8];
    ( (int* ALIGNED( 64 ))a08 )[1] = b.i[8];
    ( (int* ALIGNED( 64 ))a08 )[2] = c.i[8];
    ( (int* ALIGNED( 64 ))a08 )[3] = d.i[8];
    ( (int* ALIGNED( 64 ))a08 )[4] = e.i[8];
    ( (int* ALIGNED( 64 ))a08 )[5] = f.i[8];
    ( (int* ALIGNED( 64 ))a08 )[6] = g.i[8];
    ( (int* ALIGNED( 64 ))a08 )[7] = h.i[8];

    ( (int* ALIGNED( 64 ))a09 )[0] = a.i[9];
    ( (int* ALIGNED( 64 ))a09 )[1] = b.i[9];
    ( (int* ALIGNED( 64 ))a09 )[2] = c.i[9];
    ( (int* ALIGNED( 64 ))a09 )[3] = d.i[9];
    ( (int* ALIGNED( 64 ))a09 )[4] = e.i[9];
    ( (int* ALIGNED( 64 ))a09 )[5] = f.i[9];
    ( (int* ALIGNED( 64 ))a09 )[6] = g.i[9];
    ( (int* ALIGNED( 64 ))a09 )[7] = h.i[9];

    ( (int* ALIGNED( 64 ))a10 )[0] = a.i[10];
    ( (int* ALIGNED( 64 ))a10 )[1] = b.i[10];
    ( (int* ALIGNED( 64 ))a10 )[2] = c.i[10];
    ( (int* ALIGNED( 64 ))a10 )[3] = d.i[10];
    ( (int* ALIGNED( 64 ))a10 )[4] = e.i[10];
    ( (int* ALIGNED( 64 ))a10 )[5] = f.i[10];
    ( (int* ALIGNED( 64 ))a10 )[6] = g.i[10];
    ( (int* ALIGNED( 64 ))a10 )[7] = h.i[10];

    ( (int* ALIGNED( 64 ))a11 )[0] = a.i[11];
    ( (int* ALIGNED( 64 ))a11 )[1] = b.i[11];
    ( (int* ALIGNED( 64 ))a11 )[2] = c.i[11];
    ( (int* ALIGNED( 64 ))a11 )[3] = d.i[11];
    ( (int* ALIGNED( 64 ))a11 )[4] = e.i[11];
    ( (int* ALIGNED( 64 ))a11 )[5] = f.i[11];
    ( (int* ALIGNED( 64 ))a11 )[6] = g.i[11];
    ( (int* ALIGNED( 64 ))a11 )[7] = h.i[11];

    ( (int* ALIGNED( 64 ))a12 )[0] = a.i[12];
    ( (int* ALIGNED( 64 ))a12 )[1] = b.i[12];
    ( (int* ALIGNED( 64 ))a12 )[2] = c.i[12];
    ( (int* ALIGNED( 64 ))a12 )[3] = d.i[12];
    ( (int* ALIGNED( 64 ))a12 )[4] = e.i[12];
    ( (int* ALIGNED( 64 ))a12 )[5] = f.i[12];
    ( (int* ALIGNED( 64 ))a12 )[6] = g.i[12];
    ( (int* ALIGNED( 64 ))a12 )[7] = h.i[12];

    ( (int* ALIGNED( 64 ))a13 )[0] = a.i[13];
    ( (int* ALIGNED( 64 ))a13 )[1] = b.i[13];
    ( (int* ALIGNED( 64 ))a13 )[2] = c.i[13];
    ( (int* ALIGNED( 64 ))a13 )[3] = d.i[13];
    ( (int* ALIGNED( 64 ))a13 )[4] = e.i[13];
    ( (int* ALIGNED( 64 ))a13 )[5] = f.i[13];
    ( (int* ALIGNED( 64 ))a13 )[6] = g.i[13];
    ( (int* ALIGNED( 64 ))a13 )[7] = h.i[13];

    ( (int* ALIGNED( 64 ))a14 )[0] = a.i[14];
    ( (int* ALIGNED( 64 ))a14 )[1] = b.i[14];
    ( (int* ALIGNED( 64 ))a14 )[2] = c.i[14];
    ( (int* ALIGNED( 64 ))a14 )[3] = d.i[14];
    ( (int* ALIGNED( 64 ))a14 )[4] = e.i[14];
    ( (int* ALIGNED( 64 ))a14 )[5] = f.i[14];
    ( (int* ALIGNED( 64 ))a14 )[6] = g.i[14];
    ( (int* ALIGNED( 64 ))a14 )[7] = h.i[14];

    ( (int* ALIGNED( 64 ))a15 )[0] = a.i[15];
    ( (int* ALIGNED( 64 ))a15 )[1] = b.i[15];
    ( (int* ALIGNED( 64 ))a15 )[2] = c.i[15];
    ( (int* ALIGNED( 64 ))a15 )[3] = d.i[15];
    ( (int* ALIGNED( 64 ))a15 )[4] = e.i[15];
    ( (int* ALIGNED( 64 ))a15 )[5] = f.i[15];
    ( (int* ALIGNED( 64 ))a15 )[6] = g.i[15];
    ( (int* ALIGNED( 64 ))a15 )[7] = h.i[15];
}

inline void store_16x16_tr(
    const v16& b00, const v16& b01, const v16& b02, const v16& b03,
    const v16& b04, const v16& b05, const v16& b06, const v16& b07,
    const v16& b08, const v16& b09, const v16& b10, const v16& b11,
    const v16& b12, const v16& b13, const v16& b14, const v16& b15,
    void* ALIGNED( 64 ) a00, void* ALIGNED( 64 ) a01, void* ALIGNED( 64 ) a02,
    void* ALIGNED( 64 ) a03, void* ALIGNED( 64 ) a04, void* ALIGNED( 64 ) a05,
    void* ALIGNED( 64 ) a06, void* ALIGNED( 64 ) a07, void* ALIGNED( 64 ) a08,
    void* ALIGNED( 64 ) a09, void* ALIGNED( 64 ) a10, void* ALIGNED( 64 ) a11,
    void* ALIGNED( 64 ) a12, void* ALIGNED( 64 ) a13, void* ALIGNED( 64 ) a14,
    void* ALIGNED( 64 ) a15 )
{
    ( (int* ALIGNED( 64 ))a00 )[0] = b00.i[0];
    ( (int* ALIGNED( 64 ))a00 )[1] = b01.i[0];
    ( (int* ALIGNED( 64 ))a00 )[2] = b02.i[0];
    ( (int* ALIGNED( 64 ))a00 )[3] = b03.i[0];
    ( (int* ALIGNED( 64 ))a00 )[4] = b04.i[0];
    ( (int* ALIGNED( 64 ))a00 )[5] = b05.i[0];
    ( (int* ALIGNED( 64 ))a00 )[6] = b06.i[0];
    ( (int* ALIGNED( 64 ))a00 )[7] = b07.i[0];
    ( (int* ALIGNED( 64 ))a00 )[8] = b08.i[0];
    ( (int* ALIGNED( 64 ))a00 )[9] = b09.i[0];
    ( (int* ALIGNED( 64 ))a00 )[10] = b10.i[0];
    ( (int* ALIGNED( 64 ))a00 )[11] = b11.i[0];
    ( (int* ALIGNED( 64 ))a00 )[12] = b12.i[0];
    ( (int* ALIGNED( 64 ))a00 )[13] = b13.i[0];
    ( (int* ALIGNED( 64 ))a00 )[14] = b14.i[0];
    ( (int* ALIGNED( 64 ))a00 )[15] = b15.i[0];

    ( (int* ALIGNED( 64 ))a01 )[0] = b00.i[1];
    ( (int* ALIGNED( 64 ))a01 )[1] = b01.i[1];
    ( (int* ALIGNED( 64 ))a01 )[2] = b02.i[1];
    ( (int* ALIGNED( 64 ))a01 )[3] = b03.i[1];
    ( (int* ALIGNED( 64 ))a01 )[4] = b04.i[1];
    ( (int* ALIGNED( 64 ))a01 )[5] = b05.i[1];
    ( (int* ALIGNED( 64 ))a01 )[6] = b06.i[1];
    ( (int* ALIGNED( 64 ))a01 )[7] = b07.i[1];
    ( (int* ALIGNED( 64 ))a01 )[8] = b08.i[1];
    ( (int* ALIGNED( 64 ))a01 )[9] = b09.i[1];
    ( (int* ALIGNED( 64 ))a01 )[10] = b10.i[1];
    ( (int* ALIGNED( 64 ))a01 )[11] = b11.i[1];
    ( (int* ALIGNED( 64 ))a01 )[12] = b12.i[1];
    ( (int* ALIGNED( 64 ))a01 )[13] = b13.i[1];
    ( (int* ALIGNED( 64 ))a01 )[14] = b14.i[1];
    ( (int* ALIGNED( 64 ))a01 )[15] = b15.i[1];

    ( (int* ALIGNED( 64 ))a02 )[0] = b00.i[2];
    ( (int* ALIGNED( 64 ))a02 )[1] = b01.i[2];
    ( (int* ALIGNED( 64 ))a02 )[2] = b02.i[2];
    ( (int* ALIGNED( 64 ))a02 )[3] = b03.i[2];
    ( (int* ALIGNED( 64 ))a02 )[4] = b04.i[2];
    ( (int* ALIGNED( 64 ))a02 )[5] = b05.i[2];
    ( (int* ALIGNED( 64 ))a02 )[6] = b06.i[2];
    ( (int* ALIGNED( 64 ))a02 )[7] = b07.i[2];
    ( (int* ALIGNED( 64 ))a02 )[8] = b08.i[2];
    ( (int* ALIGNED( 64 ))a02 )[9] = b09.i[2];
    ( (int* ALIGNED( 64 ))a02 )[10] = b10.i[2];
    ( (int* ALIGNED( 64 ))a02 )[11] = b11.i[2];
    ( (int* ALIGNED( 64 ))a02 )[12] = b12.i[2];
    ( (int* ALIGNED( 64 ))a02 )[13] = b13.i[2];
    ( (int* ALIGNED( 64 ))a02 )[14] = b14.i[2];
    ( (int* ALIGNED( 64 ))a02 )[15] = b15.i[2];

    ( (int* ALIGNED( 64 ))a03 )[0] = b00.i[3];
    ( (int* ALIGNED( 64 ))a03 )[1] = b01.i[3];
    ( (int* ALIGNED( 64 ))a03 )[2] = b02.i[3];
    ( (int* ALIGNED( 64 ))a03 )[3] = b03.i[3];
    ( (int* ALIGNED( 64 ))a03 )[4] = b04.i[3];
    ( (int* ALIGNED( 64 ))a03 )[5] = b05.i[3];
    ( (int* ALIGNED( 64 ))a03 )[6] = b06.i[3];
    ( (int* ALIGNED( 64 ))a03 )[7] = b07.i[3];
    ( (int* ALIGNED( 64 ))a03 )[8] = b08.i[3];
    ( (int* ALIGNED( 64 ))a03 )[9] = b09.i[3];
    ( (int* ALIGNED( 64 ))a03 )[10] = b10.i[3];
    ( (int* ALIGNED( 64 ))a03 )[11] = b11.i[3];
    ( (int* ALIGNED( 64 ))a03 )[12] = b12.i[3];
    ( (int* ALIGNED( 64 ))a03 )[13] = b13.i[3];
    ( (int* ALIGNED( 64 ))a03 )[14] = b14.i[3];
    ( (int* ALIGNED( 64 ))a03 )[15] = b15.i[3];

    ( (int* ALIGNED( 64 ))a04 )[0] = b00.i[4];
    ( (int* ALIGNED( 64 ))a04 )[1] = b01.i[4];
    ( (int* ALIGNED( 64 ))a04 )[2] = b02.i[4];
    ( (int* ALIGNED( 64 ))a04 )[3] = b03.i[4];
    ( (int* ALIGNED( 64 ))a04 )[4] = b04.i[4];
    ( (int* ALIGNED( 64 ))a04 )[5] = b05.i[4];
    ( (int* ALIGNED( 64 ))a04 )[6] = b06.i[4];
    ( (int* ALIGNED( 64 ))a04 )[7] = b07.i[4];
    ( (int* ALIGNED( 64 ))a04 )[8] = b08.i[4];
    ( (int* ALIGNED( 64 ))a04 )[9] = b09.i[4];
    ( (int* ALIGNED( 64 ))a04 )[10] = b10.i[4];
    ( (int* ALIGNED( 64 ))a04 )[11] = b11.i[4];
    ( (int* ALIGNED( 64 ))a04 )[12] = b12.i[4];
    ( (int* ALIGNED( 64 ))a04 )[13] = b13.i[4];
    ( (int* ALIGNED( 64 ))a04 )[14] = b14.i[4];
    ( (int* ALIGNED( 64 ))a04 )[15] = b15.i[4];

    ( (int* ALIGNED( 64 ))a05 )[0] = b00.i[5];
    ( (int* ALIGNED( 64 ))a05 )[1] = b01.i[5];
    ( (int* ALIGNED( 64 ))a05 )[2] = b02.i[5];
    ( (int* ALIGNED( 64 ))a05 )[3] = b03.i[5];
    ( (int* ALIGNED( 64 ))a05 )[4] = b04.i[5];
    ( (int* ALIGNED( 64 ))a05 )[5] = b05.i[5];
    ( (int* ALIGNED( 64 ))a05 )[6] = b06.i[5];
    ( (int* ALIGNED( 64 ))a05 )[7] = b07.i[5];
    ( (int* ALIGNED( 64 ))a05 )[8] = b08.i[5];
    ( (int* ALIGNED( 64 ))a05 )[9] = b09.i[5];
    ( (int* ALIGNED( 64 ))a05 )[10] = b10.i[5];
    ( (int* ALIGNED( 64 ))a05 )[11] = b11.i[5];
    ( (int* ALIGNED( 64 ))a05 )[12] = b12.i[5];
    ( (int* ALIGNED( 64 ))a05 )[13] = b13.i[5];
    ( (int* ALIGNED( 64 ))a05 )[14] = b14.i[5];
    ( (int* ALIGNED( 64 ))a05 )[15] = b15.i[5];

    ( (int* ALIGNED( 64 ))a06 )[0] = b00.i[6];
    ( (int* ALIGNED( 64 ))a06 )[1] = b01.i[6];
    ( (int* ALIGNED( 64 ))a06 )[2] = b02.i[6];
    ( (int* ALIGNED( 64 ))a06 )[3] = b03.i[6];
    ( (int* ALIGNED( 64 ))a06 )[4] = b04.i[6];
    ( (int* ALIGNED( 64 ))a06 )[5] = b05.i[6];
    ( (int* ALIGNED( 64 ))a06 )[6] = b06.i[6];
    ( (int* ALIGNED( 64 ))a06 )[7] = b07.i[6];
    ( (int* ALIGNED( 64 ))a06 )[8] = b08.i[6];
    ( (int* ALIGNED( 64 ))a06 )[9] = b09.i[6];
    ( (int* ALIGNED( 64 ))a06 )[10] = b10.i[6];
    ( (int* ALIGNED( 64 ))a06 )[11] = b11.i[6];
    ( (int* ALIGNED( 64 ))a06 )[12] = b12.i[6];
    ( (int* ALIGNED( 64 ))a06 )[13] = b13.i[6];
    ( (int* ALIGNED( 64 ))a06 )[14] = b14.i[6];
    ( (int* ALIGNED( 64 ))a06 )[15] = b15.i[6];

    ( (int* ALIGNED( 64 ))a07 )[0] = b00.i[7];
    ( (int* ALIGNED( 64 ))a07 )[1] = b01.i[7];
    ( (int* ALIGNED( 64 ))a07 )[2] = b02.i[7];
    ( (int* ALIGNED( 64 ))a07 )[3] = b03.i[7];
    ( (int* ALIGNED( 64 ))a07 )[4] = b04.i[7];
    ( (int* ALIGNED( 64 ))a07 )[5] = b05.i[7];
    ( (int* ALIGNED( 64 ))a07 )[6] = b06.i[7];
    ( (int* ALIGNED( 64 ))a07 )[7] = b07.i[7];
    ( (int* ALIGNED( 64 ))a07 )[8] = b08.i[7];
    ( (int* ALIGNED( 64 ))a07 )[9] = b09.i[7];
    ( (int* ALIGNED( 64 ))a07 )[10] = b10.i[7];
    ( (int* ALIGNED( 64 ))a07 )[11] = b11.i[7];
    ( (int* ALIGNED( 64 ))a07 )[12] = b12.i[7];
    ( (int* ALIGNED( 64 ))a07 )[13] = b13.i[7];
    ( (int* ALIGNED( 64 ))a07 )[14] = b14.i[7];
    ( (int* ALIGNED( 64 ))a07 )[15] = b15.i[7];

    ( (int* ALIGNED( 64 ))a08 )[0] = b00.i[8];
    ( (int* ALIGNED( 64 ))a08 )[1] = b01.i[8];
    ( (int* ALIGNED( 64 ))a08 )[2] = b02.i[8];
    ( (int* ALIGNED( 64 ))a08 )[3] = b03.i[8];
    ( (int* ALIGNED( 64 ))a08 )[4] = b04.i[8];
    ( (int* ALIGNED( 64 ))a08 )[5] = b05.i[8];
    ( (int* ALIGNED( 64 ))a08 )[6] = b06.i[8];
    ( (int* ALIGNED( 64 ))a08 )[7] = b07.i[8];
    ( (int* ALIGNED( 64 ))a08 )[8] = b08.i[8];
    ( (int* ALIGNED( 64 ))a08 )[9] = b09.i[8];
    ( (int* ALIGNED( 64 ))a08 )[10] = b10.i[8];
    ( (int* ALIGNED( 64 ))a08 )[11] = b11.i[8];
    ( (int* ALIGNED( 64 ))a08 )[12] = b12.i[8];
    ( (int* ALIGNED( 64 ))a08 )[13] = b13.i[8];
    ( (int* ALIGNED( 64 ))a08 )[14] = b14.i[8];
    ( (int* ALIGNED( 64 ))a08 )[15] = b15.i[8];

    ( (int* ALIGNED( 64 ))a09 )[0] = b00.i[9];
    ( (int* ALIGNED( 64 ))a09 )[1] = b01.i[9];
    ( (int* ALIGNED( 64 ))a09 )[2] = b02.i[9];
    ( (int* ALIGNED( 64 ))a09 )[3] = b03.i[9];
    ( (int* ALIGNED( 64 ))a09 )[4] = b04.i[9];
    ( (int* ALIGNED( 64 ))a09 )[5] = b05.i[9];
    ( (int* ALIGNED( 64 ))a09 )[6] = b06.i[9];
    ( (int* ALIGNED( 64 ))a09 )[7] = b07.i[9];
    ( (int* ALIGNED( 64 ))a09 )[8] = b08.i[9];
    ( (int* ALIGNED( 64 ))a09 )[9] = b09.i[9];
    ( (int* ALIGNED( 64 ))a09 )[10] = b10.i[9];
    ( (int* ALIGNED( 64 ))a09 )[11] = b11.i[9];
    ( (int* ALIGNED( 64 ))a09 )[12] = b12.i[9];
    ( (int* ALIGNED( 64 ))a09 )[13] = b13.i[9];
    ( (int* ALIGNED( 64 ))a09 )[14] = b14.i[9];
    ( (int* ALIGNED( 64 ))a09 )[15] = b15.i[9];

    ( (int* ALIGNED( 64 ))a10 )[0] = b00.i[10];
    ( (int* ALIGNED( 64 ))a10 )[1] = b01.i[10];
    ( (int* ALIGNED( 64 ))a10 )[2] = b02.i[10];
    ( (int* ALIGNED( 64 ))a10 )[3] = b03.i[10];
    ( (int* ALIGNED( 64 ))a10 )[4] = b04.i[10];
    ( (int* ALIGNED( 64 ))a10 )[5] = b05.i[10];
    ( (int* ALIGNED( 64 ))a10 )[6] = b06.i[10];
    ( (int* ALIGNED( 64 ))a10 )[7] = b07.i[10];
    ( (int* ALIGNED( 64 ))a10 )[8] = b08.i[10];
    ( (int* ALIGNED( 64 ))a10 )[9] = b09.i[10];
    ( (int* ALIGNED( 64 ))a10 )[10] = b10.i[10];
    ( (int* ALIGNED( 64 ))a10 )[11] = b11.i[10];
    ( (int* ALIGNED( 64 ))a10 )[12] = b12.i[10];
    ( (int* ALIGNED( 64 ))a10 )[13] = b13.i[10];
    ( (int* ALIGNED( 64 ))a10 )[14] = b14.i[10];
    ( (int* ALIGNED( 64 ))a10 )[15] = b15.i[10];

    ( (int* ALIGNED( 64 ))a11 )[0] = b00.i[11];
    ( (int* ALIGNED( 64 ))a11 )[1] = b01.i[11];
    ( (int* ALIGNED( 64 ))a11 )[2] = b02.i[11];
    ( (int* ALIGNED( 64 ))a11 )[3] = b03.i[11];
    ( (int* ALIGNED( 64 ))a11 )[4] = b04.i[11];
    ( (int* ALIGNED( 64 ))a11 )[5] = b05.i[11];
    ( (int* ALIGNED( 64 ))a11 )[6] = b06.i[11];
    ( (int* ALIGNED( 64 ))a11 )[7] = b07.i[11];
    ( (int* ALIGNED( 64 ))a11 )[8] = b08.i[11];
    ( (int* ALIGNED( 64 ))a11 )[9] = b09.i[11];
    ( (int* ALIGNED( 64 ))a11 )[10] = b10.i[11];
    ( (int* ALIGNED( 64 ))a11 )[11] = b11.i[11];
    ( (int* ALIGNED( 64 ))a11 )[12] = b12.i[11];
    ( (int* ALIGNED( 64 ))a11 )[13] = b13.i[11];
    ( (int* ALIGNED( 64 ))a11 )[14] = b14.i[11];
    ( (int* ALIGNED( 64 ))a11 )[15] = b15.i[11];

    ( (int* ALIGNED( 64 ))a12 )[0] = b00.i[12];
    ( (int* ALIGNED( 64 ))a12 )[1] = b01.i[12];
    ( (int* ALIGNED( 64 ))a12 )[2] = b02.i[12];
    ( (int* ALIGNED( 64 ))a12 )[3] = b03.i[12];
    ( (int* ALIGNED( 64 ))a12 )[4] = b04.i[12];
    ( (int* ALIGNED( 64 ))a12 )[5] = b05.i[12];
    ( (int* ALIGNED( 64 ))a12 )[6] = b06.i[12];
    ( (int* ALIGNED( 64 ))a12 )[7] = b07.i[12];
    ( (int* ALIGNED( 64 ))a12 )[8] = b08.i[12];
    ( (int* ALIGNED( 64 ))a12 )[9] = b09.i[12];
    ( (int* ALIGNED( 64 ))a12 )[10] = b10.i[12];
    ( (int* ALIGNED( 64 ))a12 )[11] = b11.i[12];
    ( (int* ALIGNED( 64 ))a12 )[12] = b12.i[12];
    ( (int* ALIGNED( 64 ))a12 )[13] = b13.i[12];
    ( (int* ALIGNED( 64 ))a12 )[14] = b14.i[12];
    ( (int* ALIGNED( 64 ))a12 )[15] = b15.i[12];

    ( (int* ALIGNED( 64 ))a13 )[0] = b00.i[13];
    ( (int* ALIGNED( 64 ))a13 )[1] = b01.i[13];
    ( (int* ALIGNED( 64 ))a13 )[2] = b02.i[13];
    ( (int* ALIGNED( 64 ))a13 )[3] = b03.i[13];
    ( (int* ALIGNED( 64 ))a13 )[4] = b04.i[13];
    ( (int* ALIGNED( 64 ))a13 )[5] = b05.i[13];
    ( (int* ALIGNED( 64 ))a13 )[6] = b06.i[13];
    ( (int* ALIGNED( 64 ))a13 )[7] = b07.i[13];
    ( (int* ALIGNED( 64 ))a13 )[8] = b08.i[13];
    ( (int* ALIGNED( 64 ))a13 )[9] = b09.i[13];
    ( (int* ALIGNED( 64 ))a13 )[10] = b10.i[13];
    ( (int* ALIGNED( 64 ))a13 )[11] = b11.i[13];
    ( (int* ALIGNED( 64 ))a13 )[12] = b12.i[13];
    ( (int* ALIGNED( 64 ))a13 )[13] = b13.i[13];
    ( (int* ALIGNED( 64 ))a13 )[14] = b14.i[13];
    ( (int* ALIGNED( 64 ))a13 )[15] = b15.i[13];

    ( (int* ALIGNED( 64 ))a14 )[0] = b00.i[14];
    ( (int* ALIGNED( 64 ))a14 )[1] = b01.i[14];
    ( (int* ALIGNED( 64 ))a14 )[2] = b02.i[14];
    ( (int* ALIGNED( 64 ))a14 )[3] = b03.i[14];
    ( (int* ALIGNED( 64 ))a14 )[4] = b04.i[14];
    ( (int* ALIGNED( 64 ))a14 )[5] = b05.i[14];
    ( (int* ALIGNED( 64 ))a14 )[6] = b06.i[14];
    ( (int* ALIGNED( 64 ))a14 )[7] = b07.i[14];
    ( (int* ALIGNED( 64 ))a14 )[8] = b08.i[14];
    ( (int* ALIGNED( 64 ))a14 )[9] = b09.i[14];
    ( (int* ALIGNED( 64 ))a14 )[10] = b10.i[14];
    ( (int* ALIGNED( 64 ))a14 )[11] = b11.i[14];
    ( (int* ALIGNED( 64 ))a14 )[12] = b12.i[14];
    ( (int* ALIGNED( 64 ))a14 )[13] = b13.i[14];
    ( (int* ALIGNED( 64 ))a14 )[14] = b14.i[14];
    ( (int* ALIGNED( 64 ))a14 )[15] = b15.i[14];

    ( (int* ALIGNED( 64 ))a15 )[0] = b00.i[15];
    ( (int* ALIGNED( 64 ))a15 )[1] = b01.i[15];
    ( (int* ALIGNED( 64 ))a15 )[2] = b02.i[15];
    ( (int* ALIGNED( 64 ))a15 )[3] = b03.i[15];
    ( (int* ALIGNED( 64 ))a15 )[4] = b04.i[15];
    ( (int* ALIGNED( 64 ))a15 )[5] = b05.i[15];
    ( (int* ALIGNED( 64 ))a15 )[6] = b06.i[15];
    ( (int* ALIGNED( 64 ))a15 )[7] = b07.i[15];
    ( (int* ALIGNED( 64 ))a15 )[8] = b08.i[15];
    ( (int* ALIGNED( 64 ))a15 )[9] = b09.i[15];
    ( (int* ALIGNED( 64 ))a15 )[10] = b10.i[15];
    ( (int* ALIGNED( 64 ))a15 )[11] = b11.i[15];
    ( (int* ALIGNED( 64 ))a15 )[12] = b12.i[15];
    ( (int* ALIGNED( 64 ))a15 )[13] = b13.i[15];
    ( (int* ALIGNED( 64 ))a15 )[14] = b14.i[15];
    ( (int* ALIGNED( 64 ))a15 )[15] = b15.i[15];
}

inline void store_16x8_tr_p( const v16& b00, const v16& b01, const v16& b02,
                             const v16& b03, const v16& b04, const v16& b05,
                             const v16& b06, const v16& b07,
                             void* ALIGNED( 64 ) a00, void* ALIGNED( 64 ) a01,
                             void* ALIGNED( 64 ) a02, void* ALIGNED( 64 ) a03,
                             void* ALIGNED( 64 ) a04, void* ALIGNED( 64 ) a05,
                             void* ALIGNED( 64 ) a06, void* ALIGNED( 64 ) a07 )
{
    ( (int* ALIGNED( 64 ))a00 )[0] = b00.i[0];
    ( (int* ALIGNED( 64 ))a00 )[1] = b01.i[0];
    ( (int* ALIGNED( 64 ))a00 )[2] = b02.i[0];
    ( (int* ALIGNED( 64 ))a00 )[3] = b03.i[0];
    ( (int* ALIGNED( 64 ))a00 )[4] = b04.i[0];
    ( (int* ALIGNED( 64 ))a00 )[5] = b05.i[0];
    ( (int* ALIGNED( 64 ))a00 )[6] = b06.i[0];
    ( (int* ALIGNED( 64 ))a00 )[7] = b07.i[0];
    ( (int* ALIGNED( 64 ))a00 )[8] = b00.i[1];
    ( (int* ALIGNED( 64 ))a00 )[9] = b01.i[1];
    ( (int* ALIGNED( 64 ))a00 )[10] = b02.i[1];
    ( (int* ALIGNED( 64 ))a00 )[11] = b03.i[1];
    ( (int* ALIGNED( 64 ))a00 )[12] = b04.i[1];
    ( (int* ALIGNED( 64 ))a00 )[13] = b05.i[1];
    ( (int* ALIGNED( 64 ))a00 )[14] = b06.i[1];
    ( (int* ALIGNED( 64 ))a00 )[15] = b07.i[1];

    ( (int* ALIGNED( 64 ))a01 )[0] = b00.i[2];
    ( (int* ALIGNED( 64 ))a01 )[1] = b01.i[2];
    ( (int* ALIGNED( 64 ))a01 )[2] = b02.i[2];
    ( (int* ALIGNED( 64 ))a01 )[3] = b03.i[2];
    ( (int* ALIGNED( 64 ))a01 )[4] = b04.i[2];
    ( (int* ALIGNED( 64 ))a01 )[5] = b05.i[2];
    ( (int* ALIGNED( 64 ))a01 )[6] = b06.i[2];
    ( (int* ALIGNED( 64 ))a01 )[7] = b07.i[2];
    ( (int* ALIGNED( 64 ))a01 )[8] = b00.i[3];
    ( (int* ALIGNED( 64 ))a01 )[9] = b01.i[3];
    ( (int* ALIGNED( 64 ))a01 )[10] = b02.i[3];
    ( (int* ALIGNED( 64 ))a01 )[11] = b03.i[3];
    ( (int* ALIGNED( 64 ))a01 )[12] = b04.i[3];
    ( (int* ALIGNED( 64 ))a01 )[13] = b05.i[3];
    ( (int* ALIGNED( 64 ))a01 )[14] = b06.i[3];
    ( (int* ALIGNED( 64 ))a01 )[15] = b07.i[3];

    ( (int* ALIGNED( 64 ))a02 )[0] = b00.i[4];
    ( (int* ALIGNED( 64 ))a02 )[1] = b01.i[4];
    ( (int* ALIGNED( 64 ))a02 )[2] = b02.i[4];
    ( (int* ALIGNED( 64 ))a02 )[3] = b03.i[4];
    ( (int* ALIGNED( 64 ))a02 )[4] = b04.i[4];
    ( (int* ALIGNED( 64 ))a02 )[5] = b05.i[4];
    ( (int* ALIGNED( 64 ))a02 )[6] = b06.i[4];
    ( (int* ALIGNED( 64 ))a02 )[7] = b07.i[4];
    ( (int* ALIGNED( 64 ))a02 )[8] = b00.i[5];
    ( (int* ALIGNED( 64 ))a02 )[9] = b01.i[5];
    ( (int* ALIGNED( 64 ))a02 )[10] = b02.i[5];
    ( (int* ALIGNED( 64 ))a02 )[11] = b03.i[5];
    ( (int* ALIGNED( 64 ))a02 )[12] = b04.i[5];
    ( (int* ALIGNED( 64 ))a02 )[13] = b05.i[5];
    ( (int* ALIGNED( 64 ))a02 )[14] = b06.i[5];
    ( (int* ALIGNED( 64 ))a02 )[15] = b07.i[5];

    ( (int* ALIGNED( 64 ))a03 )[0] = b00.i[6];
    ( (int* ALIGNED( 64 ))a03 )[1] = b01.i[6];
    ( (int* ALIGNED( 64 ))a03 )[2] = b02.i[6];
    ( (int* ALIGNED( 64 ))a03 )[3] = b03.i[6];
    ( (int* ALIGNED( 64 ))a03 )[4] = b04.i[6];
    ( (int* ALIGNED( 64 ))a03 )[5] = b05.i[6];
    ( (int* ALIGNED( 64 ))a03 )[6] = b06.i[6];
    ( (int* ALIGNED( 64 ))a03 )[7] = b07.i[6];
    ( (int* ALIGNED( 64 ))a03 )[8] = b00.i[7];
    ( (int* ALIGNED( 64 ))a03 )[9] = b01.i[7];
    ( (int* ALIGNED( 64 ))a03 )[10] = b02.i[7];
    ( (int* ALIGNED( 64 ))a03 )[11] = b03.i[7];
    ( (int* ALIGNED( 64 ))a03 )[12] = b04.i[7];
    ( (int* ALIGNED( 64 ))a03 )[13] = b05.i[7];
    ( (int* ALIGNED( 64 ))a03 )[14] = b06.i[7];
    ( (int* ALIGNED( 64 ))a03 )[15] = b07.i[7];

    ( (int* ALIGNED( 64 ))a04 )[0] = b00.i[8];
    ( (int* ALIGNED( 64 ))a04 )[1] = b01.i[8];
    ( (int* ALIGNED( 64 ))a04 )[2] = b02.i[8];
    ( (int* ALIGNED( 64 ))a04 )[3] = b03.i[8];
    ( (int* ALIGNED( 64 ))a04 )[4] = b04.i[8];
    ( (int* ALIGNED( 64 ))a04 )[5] = b05.i[8];
    ( (int* ALIGNED( 64 ))a04 )[6] = b06.i[8];
    ( (int* ALIGNED( 64 ))a04 )[7] = b07.i[8];
    ( (int* ALIGNED( 64 ))a04 )[8] = b00.i[9];
    ( (int* ALIGNED( 64 ))a04 )[9] = b01.i[9];
    ( (int* ALIGNED( 64 ))a04 )[10] = b02.i[9];
    ( (int* ALIGNED( 64 ))a04 )[11] = b03.i[9];
    ( (int* ALIGNED( 64 ))a04 )[12] = b04.i[9];
    ( (int* ALIGNED( 64 ))a04 )[13] = b05.i[9];
    ( (int* ALIGNED( 64 ))a04 )[14] = b06.i[9];
    ( (int* ALIGNED( 64 ))a04 )[15] = b07.i[9];

    ( (int* ALIGNED( 64 ))a05 )[0] = b00.i[10];
    ( (int* ALIGNED( 64 ))a05 )[1] = b01.i[10];
    ( (int* ALIGNED( 64 ))a05 )[2] = b02.i[10];
    ( (int* ALIGNED( 64 ))a05 )[3] = b03.i[10];
    ( (int* ALIGNED( 64 ))a05 )[4] = b04.i[10];
    ( (int* ALIGNED( 64 ))a05 )[5] = b05.i[10];
    ( (int* ALIGNED( 64 ))a05 )[6] = b06.i[10];
    ( (int* ALIGNED( 64 ))a05 )[7] = b07.i[10];
    ( (int* ALIGNED( 64 ))a05 )[8] = b00.i[11];
    ( (int* ALIGNED( 64 ))a05 )[9] = b01.i[11];
    ( (int* ALIGNED( 64 ))a05 )[10] = b02.i[11];
    ( (int* ALIGNED( 64 ))a05 )[11] = b03.i[11];
    ( (int* ALIGNED( 64 ))a05 )[12] = b04.i[11];
    ( (int* ALIGNED( 64 ))a05 )[13] = b05.i[11];
    ( (int* ALIGNED( 64 ))a05 )[14] = b06.i[11];
    ( (int* ALIGNED( 64 ))a05 )[15] = b07.i[11];

    ( (int* ALIGNED( 64 ))a06 )[0] = b00.i[12];
    ( (int* ALIGNED( 64 ))a06 )[1] = b01.i[12];
    ( (int* ALIGNED( 64 ))a06 )[2] = b02.i[12];
    ( (int* ALIGNED( 64 ))a06 )[3] = b03.i[12];
    ( (int* ALIGNED( 64 ))a06 )[4] = b04.i[12];
    ( (int* ALIGNED( 64 ))a06 )[5] = b05.i[12];
    ( (int* ALIGNED( 64 ))a06 )[6] = b06.i[12];
    ( (int* ALIGNED( 64 ))a06 )[7] = b07.i[12];
    ( (int* ALIGNED( 64 ))a06 )[8] = b00.i[13];
    ( (int* ALIGNED( 64 ))a06 )[9] = b01.i[13];
    ( (int* ALIGNED( 64 ))a06 )[10] = b02.i[13];
    ( (int* ALIGNED( 64 ))a06 )[11] = b03.i[13];
    ( (int* ALIGNED( 64 ))a06 )[12] = b04.i[13];
    ( (int* ALIGNED( 64 ))a06 )[13] = b05.i[13];
    ( (int* ALIGNED( 64 ))a06 )[14] = b06.i[13];
    ( (int* ALIGNED( 64 ))a06 )[15] = b07.i[13];

    ( (int* ALIGNED( 64 ))a07 )[0] = b00.i[14];
    ( (int* ALIGNED( 64 ))a07 )[1] = b01.i[14];
    ( (int* ALIGNED( 64 ))a07 )[2] = b02.i[14];
    ( (int* ALIGNED( 64 ))a07 )[3] = b03.i[14];
    ( (int* ALIGNED( 64 ))a07 )[4] = b04.i[14];
    ( (int* ALIGNED( 64 ))a07 )[5] = b05.i[14];
    ( (int* ALIGNED( 64 ))a07 )[6] = b06.i[14];
    ( (int* ALIGNED( 64 ))a07 )[7] = b07.i[14];
    ( (int* ALIGNED( 64 ))a07 )[8] = b00.i[15];
    ( (int* ALIGNED( 64 ))a07 )[9] = b01.i[15];
    ( (int* ALIGNED( 64 ))a07 )[10] = b02.i[15];
    ( (int* ALIGNED( 64 ))a07 )[11] = b03.i[15];
    ( (int* ALIGNED( 64 ))a07 )[12] = b04.i[15];
    ( (int* ALIGNED( 64 ))a07 )[13] = b05.i[15];
    ( (int* ALIGNED( 64 ))a07 )[14] = b06.i[15];
    ( (int* ALIGNED( 64 ))a07 )[15] = b07.i[15];
}

inline void store_16x16_tr_p(
    const v16& b00, const v16& b01, const v16& b02, const v16& b03,
    const v16& b04, const v16& b05, const v16& b06, const v16& b07,
    const v16& b08, const v16& b09, const v16& b10, const v16& b11,
    const v16& b12, const v16& b13, const v16& b14, const v16& b15,
    void* ALIGNED( 64 ) a00, void* ALIGNED( 64 ) a01, void* ALIGNED( 64 ) a02,
    void* ALIGNED( 64 ) a03, void* ALIGNED( 64 ) a04, void* ALIGNED( 64 ) a05,
    void* ALIGNED( 64 ) a06, void* ALIGNED( 64 ) a07, void* ALIGNED( 64 ) a08,
    void* ALIGNED( 64 ) a09, void* ALIGNED( 64 ) a10, void* ALIGNED( 64 ) a11,
    void* ALIGNED( 64 ) a12, void* ALIGNED( 64 ) a13, void* ALIGNED( 64 ) a14,
    void* ALIGNED( 64 ) a15 )
{
    ( (int* ALIGNED( 64 ))a00 )[0] = b00.i[0];
    ( (int* ALIGNED( 64 ))a00 )[1] = b01.i[0];
    ( (int* ALIGNED( 64 ))a00 )[2] = b02.i[0];
    ( (int* ALIGNED( 64 ))a00 )[3] = b03.i[0];
    ( (int* ALIGNED( 64 ))a00 )[4] = b04.i[0];
    ( (int* ALIGNED( 64 ))a00 )[5] = b05.i[0];
    ( (int* ALIGNED( 64 ))a00 )[6] = b06.i[0];
    ( (int* ALIGNED( 64 ))a00 )[7] = b07.i[0];
    ( (int* ALIGNED( 64 ))a00 )[8] = b00.i[1];
    ( (int* ALIGNED( 64 ))a00 )[9] = b01.i[1];
    ( (int* ALIGNED( 64 ))a00 )[10] = b02.i[1];
    ( (int* ALIGNED( 64 ))a00 )[11] = b03.i[1];
    ( (int* ALIGNED( 64 ))a00 )[12] = b04.i[1];
    ( (int* ALIGNED( 64 ))a00 )[13] = b05.i[1];
    ( (int* ALIGNED( 64 ))a00 )[14] = b06.i[1];
    ( (int* ALIGNED( 64 ))a00 )[15] = b07.i[1];

    ( (int* ALIGNED( 64 ))a01 )[0] = b00.i[2];
    ( (int* ALIGNED( 64 ))a01 )[1] = b01.i[2];
    ( (int* ALIGNED( 64 ))a01 )[2] = b02.i[2];
    ( (int* ALIGNED( 64 ))a01 )[3] = b03.i[2];
    ( (int* ALIGNED( 64 ))a01 )[4] = b04.i[2];
    ( (int* ALIGNED( 64 ))a01 )[5] = b05.i[2];
    ( (int* ALIGNED( 64 ))a01 )[6] = b06.i[2];
    ( (int* ALIGNED( 64 ))a01 )[7] = b07.i[2];
    ( (int* ALIGNED( 64 ))a01 )[8] = b00.i[3];
    ( (int* ALIGNED( 64 ))a01 )[9] = b01.i[3];
    ( (int* ALIGNED( 64 ))a01 )[10] = b02.i[3];
    ( (int* ALIGNED( 64 ))a01 )[11] = b03.i[3];
    ( (int* ALIGNED( 64 ))a01 )[12] = b04.i[3];
    ( (int* ALIGNED( 64 ))a01 )[13] = b05.i[3];
    ( (int* ALIGNED( 64 ))a01 )[14] = b06.i[3];
    ( (int* ALIGNED( 64 ))a01 )[15] = b07.i[3];

    ( (int* ALIGNED( 64 ))a02 )[0] = b00.i[4];
    ( (int* ALIGNED( 64 ))a02 )[1] = b01.i[4];
    ( (int* ALIGNED( 64 ))a02 )[2] = b02.i[4];
    ( (int* ALIGNED( 64 ))a02 )[3] = b03.i[4];
    ( (int* ALIGNED( 64 ))a02 )[4] = b04.i[4];
    ( (int* ALIGNED( 64 ))a02 )[5] = b05.i[4];
    ( (int* ALIGNED( 64 ))a02 )[6] = b06.i[4];
    ( (int* ALIGNED( 64 ))a02 )[7] = b07.i[4];
    ( (int* ALIGNED( 64 ))a02 )[8] = b00.i[5];
    ( (int* ALIGNED( 64 ))a02 )[9] = b01.i[5];
    ( (int* ALIGNED( 64 ))a02 )[10] = b02.i[5];
    ( (int* ALIGNED( 64 ))a02 )[11] = b03.i[5];
    ( (int* ALIGNED( 64 ))a02 )[12] = b04.i[5];
    ( (int* ALIGNED( 64 ))a02 )[13] = b05.i[5];
    ( (int* ALIGNED( 64 ))a02 )[14] = b06.i[5];
    ( (int* ALIGNED( 64 ))a02 )[15] = b07.i[5];

    ( (int* ALIGNED( 64 ))a03 )[0] = b00.i[6];
    ( (int* ALIGNED( 64 ))a03 )[1] = b01.i[6];
    ( (int* ALIGNED( 64 ))a03 )[2] = b02.i[6];
    ( (int* ALIGNED( 64 ))a03 )[3] = b03.i[6];
    ( (int* ALIGNED( 64 ))a03 )[4] = b04.i[6];
    ( (int* ALIGNED( 64 ))a03 )[5] = b05.i[6];
    ( (int* ALIGNED( 64 ))a03 )[6] = b06.i[6];
    ( (int* ALIGNED( 64 ))a03 )[7] = b07.i[6];
    ( (int* ALIGNED( 64 ))a03 )[8] = b00.i[7];
    ( (int* ALIGNED( 64 ))a03 )[9] = b01.i[7];
    ( (int* ALIGNED( 64 ))a03 )[10] = b02.i[7];
    ( (int* ALIGNED( 64 ))a03 )[11] = b03.i[7];
    ( (int* ALIGNED( 64 ))a03 )[12] = b04.i[7];
    ( (int* ALIGNED( 64 ))a03 )[13] = b05.i[7];
    ( (int* ALIGNED( 64 ))a03 )[14] = b06.i[7];
    ( (int* ALIGNED( 64 ))a03 )[15] = b07.i[7];

    ( (int* ALIGNED( 64 ))a04 )[0] = b00.i[8];
    ( (int* ALIGNED( 64 ))a04 )[1] = b01.i[8];
    ( (int* ALIGNED( 64 ))a04 )[2] = b02.i[8];
    ( (int* ALIGNED( 64 ))a04 )[3] = b03.i[8];
    ( (int* ALIGNED( 64 ))a04 )[4] = b04.i[8];
    ( (int* ALIGNED( 64 ))a04 )[5] = b05.i[8];
    ( (int* ALIGNED( 64 ))a04 )[6] = b06.i[8];
    ( (int* ALIGNED( 64 ))a04 )[7] = b07.i[8];
    ( (int* ALIGNED( 64 ))a04 )[8] = b00.i[9];
    ( (int* ALIGNED( 64 ))a04 )[9] = b01.i[9];
    ( (int* ALIGNED( 64 ))a04 )[10] = b02.i[9];
    ( (int* ALIGNED( 64 ))a04 )[11] = b03.i[9];
    ( (int* ALIGNED( 64 ))a04 )[12] = b04.i[9];
    ( (int* ALIGNED( 64 ))a04 )[13] = b05.i[9];
    ( (int* ALIGNED( 64 ))a04 )[14] = b06.i[9];
    ( (int* ALIGNED( 64 ))a04 )[15] = b07.i[9];

    ( (int* ALIGNED( 64 ))a05 )[0] = b00.i[10];
    ( (int* ALIGNED( 64 ))a05 )[1] = b01.i[10];
    ( (int* ALIGNED( 64 ))a05 )[2] = b02.i[10];
    ( (int* ALIGNED( 64 ))a05 )[3] = b03.i[10];
    ( (int* ALIGNED( 64 ))a05 )[4] = b04.i[10];
    ( (int* ALIGNED( 64 ))a05 )[5] = b05.i[10];
    ( (int* ALIGNED( 64 ))a05 )[6] = b06.i[10];
    ( (int* ALIGNED( 64 ))a05 )[7] = b07.i[10];
    ( (int* ALIGNED( 64 ))a05 )[8] = b00.i[11];
    ( (int* ALIGNED( 64 ))a05 )[9] = b01.i[11];
    ( (int* ALIGNED( 64 ))a05 )[10] = b02.i[11];
    ( (int* ALIGNED( 64 ))a05 )[11] = b03.i[11];
    ( (int* ALIGNED( 64 ))a05 )[12] = b04.i[11];
    ( (int* ALIGNED( 64 ))a05 )[13] = b05.i[11];
    ( (int* ALIGNED( 64 ))a05 )[14] = b06.i[11];
    ( (int* ALIGNED( 64 ))a05 )[15] = b07.i[11];

    ( (int* ALIGNED( 64 ))a06 )[0] = b00.i[12];
    ( (int* ALIGNED( 64 ))a06 )[1] = b01.i[12];
    ( (int* ALIGNED( 64 ))a06 )[2] = b02.i[12];
    ( (int* ALIGNED( 64 ))a06 )[3] = b03.i[12];
    ( (int* ALIGNED( 64 ))a06 )[4] = b04.i[12];
    ( (int* ALIGNED( 64 ))a06 )[5] = b05.i[12];
    ( (int* ALIGNED( 64 ))a06 )[6] = b06.i[12];
    ( (int* ALIGNED( 64 ))a06 )[7] = b07.i[12];
    ( (int* ALIGNED( 64 ))a06 )[8] = b00.i[13];
    ( (int* ALIGNED( 64 ))a06 )[9] = b01.i[13];
    ( (int* ALIGNED( 64 ))a06 )[10] = b02.i[13];
    ( (int* ALIGNED( 64 ))a06 )[11] = b03.i[13];
    ( (int* ALIGNED( 64 ))a06 )[12] = b04.i[13];
    ( (int* ALIGNED( 64 ))a06 )[13] = b05.i[13];
    ( (int* ALIGNED( 64 ))a06 )[14] = b06.i[13];
    ( (int* ALIGNED( 64 ))a06 )[15] = b07.i[13];

    ( (int* ALIGNED( 64 ))a07 )[0] = b00.i[14];
    ( (int* ALIGNED( 64 ))a07 )[1] = b01.i[14];
    ( (int* ALIGNED( 64 ))a07 )[2] = b02.i[14];
    ( (int* ALIGNED( 64 ))a07 )[3] = b03.i[14];
    ( (int* ALIGNED( 64 ))a07 )[4] = b04.i[14];
    ( (int* ALIGNED( 64 ))a07 )[5] = b05.i[14];
    ( (int* ALIGNED( 64 ))a07 )[6] = b06.i[14];
    ( (int* ALIGNED( 64 ))a07 )[7] = b07.i[14];
    ( (int* ALIGNED( 64 ))a07 )[8] = b00.i[15];
    ( (int* ALIGNED( 64 ))a07 )[9] = b01.i[15];
    ( (int* ALIGNED( 64 ))a07 )[10] = b02.i[15];
    ( (int* ALIGNED( 64 ))a07 )[11] = b03.i[15];
    ( (int* ALIGNED( 64 ))a07 )[12] = b04.i[15];
    ( (int* ALIGNED( 64 ))a07 )[13] = b05.i[15];
    ( (int* ALIGNED( 64 ))a07 )[14] = b06.i[15];
    ( (int* ALIGNED( 64 ))a07 )[15] = b07.i[15];

    ( (int* ALIGNED( 64 ))a08 )[0] = b08.i[0];
    ( (int* ALIGNED( 64 ))a08 )[1] = b09.i[0];
    ( (int* ALIGNED( 64 ))a08 )[2] = b10.i[0];
    ( (int* ALIGNED( 64 ))a08 )[3] = b11.i[0];
    ( (int* ALIGNED( 64 ))a08 )[4] = b12.i[0];
    ( (int* ALIGNED( 64 ))a08 )[5] = b13.i[0];
    ( (int* ALIGNED( 64 ))a08 )[6] = b14.i[0];
    ( (int* ALIGNED( 64 ))a08 )[7] = b15.i[0];
    ( (int* ALIGNED( 64 ))a08 )[8] = b08.i[1];
    ( (int* ALIGNED( 64 ))a08 )[9] = b09.i[1];
    ( (int* ALIGNED( 64 ))a08 )[10] = b10.i[1];
    ( (int* ALIGNED( 64 ))a08 )[11] = b11.i[1];
    ( (int* ALIGNED( 64 ))a08 )[12] = b12.i[1];
    ( (int* ALIGNED( 64 ))a08 )[13] = b13.i[1];
    ( (int* ALIGNED( 64 ))a08 )[14] = b14.i[1];
    ( (int* ALIGNED( 64 ))a08 )[15] = b15.i[1];

    ( (int* ALIGNED( 64 ))a09 )[0] = b08.i[2];
    ( (int* ALIGNED( 64 ))a09 )[1] = b09.i[2];
    ( (int* ALIGNED( 64 ))a09 )[2] = b10.i[2];
    ( (int* ALIGNED( 64 ))a09 )[3] = b11.i[2];
    ( (int* ALIGNED( 64 ))a09 )[4] = b12.i[2];
    ( (int* ALIGNED( 64 ))a09 )[5] = b13.i[2];
    ( (int* ALIGNED( 64 ))a09 )[6] = b14.i[2];
    ( (int* ALIGNED( 64 ))a09 )[7] = b15.i[2];
    ( (int* ALIGNED( 64 ))a09 )[8] = b08.i[3];
    ( (int* ALIGNED( 64 ))a09 )[9] = b09.i[3];
    ( (int* ALIGNED( 64 ))a09 )[10] = b10.i[3];
    ( (int* ALIGNED( 64 ))a09 )[11] = b11.i[3];
    ( (int* ALIGNED( 64 ))a09 )[12] = b12.i[3];
    ( (int* ALIGNED( 64 ))a09 )[13] = b13.i[3];
    ( (int* ALIGNED( 64 ))a09 )[14] = b14.i[3];
    ( (int* ALIGNED( 64 ))a09 )[15] = b15.i[3];

    ( (int* ALIGNED( 64 ))a10 )[0] = b08.i[4];
    ( (int* ALIGNED( 64 ))a10 )[1] = b09.i[4];
    ( (int* ALIGNED( 64 ))a10 )[2] = b10.i[4];
    ( (int* ALIGNED( 64 ))a10 )[3] = b11.i[4];
    ( (int* ALIGNED( 64 ))a10 )[4] = b12.i[4];
    ( (int* ALIGNED( 64 ))a10 )[5] = b13.i[4];
    ( (int* ALIGNED( 64 ))a10 )[6] = b14.i[4];
    ( (int* ALIGNED( 64 ))a10 )[7] = b15.i[4];
    ( (int* ALIGNED( 64 ))a10 )[8] = b08.i[5];
    ( (int* ALIGNED( 64 ))a10 )[9] = b09.i[5];
    ( (int* ALIGNED( 64 ))a10 )[10] = b10.i[5];
    ( (int* ALIGNED( 64 ))a10 )[11] = b11.i[5];
    ( (int* ALIGNED( 64 ))a10 )[12] = b12.i[5];
    ( (int* ALIGNED( 64 ))a10 )[13] = b13.i[5];
    ( (int* ALIGNED( 64 ))a10 )[14] = b14.i[5];
    ( (int* ALIGNED( 64 ))a10 )[15] = b15.i[5];

    ( (int* ALIGNED( 64 ))a11 )[0] = b08.i[6];
    ( (int* ALIGNED( 64 ))a11 )[1] = b09.i[6];
    ( (int* ALIGNED( 64 ))a11 )[2] = b10.i[6];
    ( (int* ALIGNED( 64 ))a11 )[3] = b11.i[6];
    ( (int* ALIGNED( 64 ))a11 )[4] = b12.i[6];
    ( (int* ALIGNED( 64 ))a11 )[5] = b13.i[6];
    ( (int* ALIGNED( 64 ))a11 )[6] = b14.i[6];
    ( (int* ALIGNED( 64 ))a11 )[7] = b15.i[6];
    ( (int* ALIGNED( 64 ))a11 )[8] = b08.i[7];
    ( (int* ALIGNED( 64 ))a11 )[9] = b09.i[7];
    ( (int* ALIGNED( 64 ))a11 )[10] = b10.i[7];
    ( (int* ALIGNED( 64 ))a11 )[11] = b11.i[7];
    ( (int* ALIGNED( 64 ))a11 )[12] = b12.i[7];
    ( (int* ALIGNED( 64 ))a11 )[13] = b13.i[7];
    ( (int* ALIGNED( 64 ))a11 )[14] = b14.i[7];
    ( (int* ALIGNED( 64 ))a11 )[15] = b15.i[7];

    ( (int* ALIGNED( 64 ))a12 )[0] = b08.i[8];
    ( (int* ALIGNED( 64 ))a12 )[1] = b09.i[8];
    ( (int* ALIGNED( 64 ))a12 )[2] = b10.i[8];
    ( (int* ALIGNED( 64 ))a12 )[3] = b11.i[8];
    ( (int* ALIGNED( 64 ))a12 )[4] = b12.i[8];
    ( (int* ALIGNED( 64 ))a12 )[5] = b13.i[8];
    ( (int* ALIGNED( 64 ))a12 )[6] = b14.i[8];
    ( (int* ALIGNED( 64 ))a12 )[7] = b15.i[8];
    ( (int* ALIGNED( 64 ))a12 )[8] = b08.i[9];
    ( (int* ALIGNED( 64 ))a12 )[9] = b09.i[9];
    ( (int* ALIGNED( 64 ))a12 )[10] = b10.i[9];
    ( (int* ALIGNED( 64 ))a12 )[11] = b11.i[9];
    ( (int* ALIGNED( 64 ))a12 )[12] = b12.i[9];
    ( (int* ALIGNED( 64 ))a12 )[13] = b13.i[9];
    ( (int* ALIGNED( 64 ))a12 )[14] = b14.i[9];
    ( (int* ALIGNED( 64 ))a12 )[15] = b15.i[9];

    ( (int* ALIGNED( 64 ))a13 )[0] = b08.i[10];
    ( (int* ALIGNED( 64 ))a13 )[1] = b09.i[10];
    ( (int* ALIGNED( 64 ))a13 )[2] = b10.i[10];
    ( (int* ALIGNED( 64 ))a13 )[3] = b11.i[10];
    ( (int* ALIGNED( 64 ))a13 )[4] = b12.i[10];
    ( (int* ALIGNED( 64 ))a13 )[5] = b13.i[10];
    ( (int* ALIGNED( 64 ))a13 )[6] = b14.i[10];
    ( (int* ALIGNED( 64 ))a13 )[7] = b15.i[10];
    ( (int* ALIGNED( 64 ))a13 )[8] = b08.i[11];
    ( (int* ALIGNED( 64 ))a13 )[9] = b09.i[11];
    ( (int* ALIGNED( 64 ))a13 )[10] = b10.i[11];
    ( (int* ALIGNED( 64 ))a13 )[11] = b11.i[11];
    ( (int* ALIGNED( 64 ))a13 )[12] = b12.i[11];
    ( (int* ALIGNED( 64 ))a13 )[13] = b13.i[11];
    ( (int* ALIGNED( 64 ))a13 )[14] = b14.i[11];
    ( (int* ALIGNED( 64 ))a13 )[15] = b15.i[11];

    ( (int* ALIGNED( 64 ))a14 )[0] = b08.i[12];
    ( (int* ALIGNED( 64 ))a14 )[1] = b09.i[12];
    ( (int* ALIGNED( 64 ))a14 )[2] = b10.i[12];
    ( (int* ALIGNED( 64 ))a14 )[3] = b11.i[12];
    ( (int* ALIGNED( 64 ))a14 )[4] = b12.i[12];
    ( (int* ALIGNED( 64 ))a14 )[5] = b13.i[12];
    ( (int* ALIGNED( 64 ))a14 )[6] = b14.i[12];
    ( (int* ALIGNED( 64 ))a14 )[7] = b15.i[12];
    ( (int* ALIGNED( 64 ))a14 )[8] = b08.i[13];
    ( (int* ALIGNED( 64 ))a14 )[9] = b09.i[13];
    ( (int* ALIGNED( 64 ))a14 )[10] = b10.i[13];
    ( (int* ALIGNED( 64 ))a14 )[11] = b11.i[13];
    ( (int* ALIGNED( 64 ))a14 )[12] = b12.i[13];
    ( (int* ALIGNED( 64 ))a14 )[13] = b13.i[13];
    ( (int* ALIGNED( 64 ))a14 )[14] = b14.i[13];
    ( (int* ALIGNED( 64 ))a14 )[15] = b15.i[13];

    ( (int* ALIGNED( 64 ))a15 )[0] = b08.i[14];
    ( (int* ALIGNED( 64 ))a15 )[1] = b09.i[14];
    ( (int* ALIGNED( 64 ))a15 )[2] = b10.i[14];
    ( (int* ALIGNED( 64 ))a15 )[3] = b11.i[14];
    ( (int* ALIGNED( 64 ))a15 )[4] = b12.i[14];
    ( (int* ALIGNED( 64 ))a15 )[5] = b13.i[14];
    ( (int* ALIGNED( 64 ))a15 )[6] = b14.i[14];
    ( (int* ALIGNED( 64 ))a15 )[7] = b15.i[14];
    ( (int* ALIGNED( 64 ))a15 )[8] = b08.i[15];
    ( (int* ALIGNED( 64 ))a15 )[9] = b09.i[15];
    ( (int* ALIGNED( 64 ))a15 )[10] = b10.i[15];
    ( (int* ALIGNED( 64 ))a15 )[11] = b11.i[15];
    ( (int* ALIGNED( 64 ))a15 )[12] = b12.i[15];
    ( (int* ALIGNED( 64 ))a15 )[13] = b13.i[15];
    ( (int* ALIGNED( 64 ))a15 )[14] = b14.i[15];
    ( (int* ALIGNED( 64 ))a15 )[15] = b15.i[15];
}

//////////////
// v16int class

class v16int : public v16
{
    // v16int prefix unary operator friends

    friend inline v16int operator+( const v16int& a ) ALWAYS_INLINE;
    friend inline v16int operator-( const v16int& a ) ALWAYS_INLINE;
    friend inline v16int operator~( const v16int& a ) ALWAYS_INLINE;
    friend inline v16int operator!( const v16int& a ) ALWAYS_INLINE;
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector

    // v16int prefix increment / decrement operator friends

    friend inline v16int operator++( v16int& a ) ALWAYS_INLINE;
    friend inline v16int operator--( v16int& a ) ALWAYS_INLINE;

    // v16int postfix increment / decrement operator friends

    friend inline v16int operator++( v16int& a, int ) ALWAYS_INLINE;
    friend inline v16int operator--( v16int& a, int ) ALWAYS_INLINE;

    // v16int binary operator friends

    friend inline v16int operator+( const v16int& a,
                                    const v16int& b ) ALWAYS_INLINE;
    friend inline v16int operator-( const v16int& a,
                                    const v16int& b ) ALWAYS_INLINE;
    friend inline v16int operator*(const v16int& a,
                                   const v16int& b)ALWAYS_INLINE;
    friend inline v16int operator/( const v16int& a,
                                    const v16int& b ) ALWAYS_INLINE;
    friend inline v16int operator%( const v16int& a,
                                    const v16int& b ) ALWAYS_INLINE;
    friend inline v16int operator^( const v16int& a,
                                    const v16int& b ) ALWAYS_INLINE;
    friend inline v16int operator&(const v16int& a,
                                   const v16int& b)ALWAYS_INLINE;
    friend inline v16int operator|( const v16int& a,
                                    const v16int& b ) ALWAYS_INLINE;
    friend inline v16int operator<<( const v16int& a,
                                     const v16int& b ) ALWAYS_INLINE;
    friend inline v16int operator>>( const v16int& a,
                                     const v16int& b ) ALWAYS_INLINE;

    // v16int logical operator friends

    friend inline v16int operator<( const v16int& a,
                                    const v16int& b ) ALWAYS_INLINE;
    friend inline v16int operator>( const v16int& a,
                                    const v16int& b ) ALWAYS_INLINE;
    friend inline v16int operator==( const v16int& a,
                                     const v16int& b ) ALWAYS_INLINE;
    friend inline v16int operator!=( const v16int& a,
                                     const v16int& b ) ALWAYS_INLINE;
    friend inline v16int operator<=( const v16int& a,
                                     const v16int& b ) ALWAYS_INLINE;
    friend inline v16int operator>=( const v16int& a,
                                     const v16int& b ) ALWAYS_INLINE;
    friend inline v16int operator&&( const v16int& a,
                                     const v16int& b ) ALWAYS_INLINE;
    friend inline v16int operator||( const v16int& a,
                                     const v16int& b ) ALWAYS_INLINE;

    // v16int miscellaneous friends

    friend inline v16int abs( const v16int& a ) ALWAYS_INLINE;
    friend inline v16 czero( const v16int& c, const v16& a ) ALWAYS_INLINE;
    friend inline v16 notczero( const v16int& c, const v16& a ) ALWAYS_INLINE;
    // FIXME: cswap, notcswap!
    friend inline v16 merge( const v16int& c, const v16& t,
                             const v16& f ) ALWAYS_INLINE;

    // v16float unary operator friends

    friend inline v16int operator!( const v16float& a ) ALWAYS_INLINE;

    // v16float logical operator friends

    friend inline v16int operator<( const v16float& a,
                                    const v16float& b ) ALWAYS_INLINE;
    friend inline v16int operator>( const v16float& a,
                                    const v16float& b ) ALWAYS_INLINE;
    friend inline v16int operator==( const v16float& a,
                                     const v16float& b ) ALWAYS_INLINE;
    friend inline v16int operator!=( const v16float& a,
                                     const v16float& b ) ALWAYS_INLINE;
    friend inline v16int operator<=( const v16float& a,
                                     const v16float& b ) ALWAYS_INLINE;
    friend inline v16int operator>=( const v16float& a,
                                     const v16float& b ) ALWAYS_INLINE;
    friend inline v16int operator&&( const v16float& a,
                                     const v16float& b ) ALWAYS_INLINE;
    friend inline v16int operator||( const v16float& a,
                                     const v16float& b ) ALWAYS_INLINE;

    // v16float miscellaneous friends

    friend inline v16float clear_bits( const v16int& m,
                                       const v16float& a ) ALWAYS_INLINE;
    friend inline v16float set_bits( const v16int& m,
                                     const v16float& a ) ALWAYS_INLINE;
    friend inline v16float toggle_bits( const v16int& m,
                                        const v16float& a ) ALWAYS_INLINE;

  public:
    // v16int constructors / destructors

    v16int() {} // Default constructor

    v16int( const v16int& a ) // Copy constructor
    {
        ALWAYS_VECTORIZE
        for ( int j = 0; j < 16; j++ )
            i[j] = a.i[j];
    }

    v16int( const v16& a ) // Init from mixed
    {
        ALWAYS_VECTORIZE
        for ( int j = 0; j < 16; j++ )
            i[j] = a.i[j];
    }

    v16int( int a ) // Init from scalar
    {
        ALWAYS_VECTORIZE
        for ( int j = 0; j < 16; j++ )
            i[j] = a;
    }

    v16int( int i00, int i01, int i02, int i03, int i04, int i05, int i06,
            int i07, int i08, int i09, int i10, int i11, int i12, int i13,
            int i14, int i15 ) // Init from scalars
    {
        i[0] = i00;
        i[1] = i01;
        i[2] = i02;
        i[3] = i03;
        i[4] = i04;
        i[5] = i05;
        i[6] = i06;
        i[7] = i07;
        i[8] = i08;
        i[9] = i09;
        i[10] = i10;
        i[11] = i11;
        i[12] = i12;
        i[13] = i13;
        i[14] = i14;
        i[15] = i15;
    }

    ~v16int() {} // Destructor

    // v16int assignment operators

#define ASSIGN( op )                                                           \
    inline v16int& operator op( const v16int& b )                              \
    {                                                                          \
        ALWAYS_VECTORIZE                                                       \
        for ( int j = 0; j < 16; j++ )                                         \
            i[j] op b.i[j];                                                    \
        return *this;                                                          \
    }

    ASSIGN( = )
    ASSIGN( += )
    ASSIGN( -= )
    ASSIGN( *= )
    ASSIGN( /= )
    ASSIGN( %= )
    ASSIGN( ^= )
    ASSIGN( &= )
    ASSIGN( |= )
    ASSIGN( <<= )
    ASSIGN( >>= )

#undef ASSIGN

    // v16int member access operator

    inline int& operator[]( int n ) { return i[n]; }

    inline int operator()( int n ) { return i[n]; }
};

// v16int prefix unary operators

#define PREFIX_UNARY( op )                                                     \
    inline v16int operator op( const v16int& a )                               \
    {                                                                          \
        v16int b;                                                              \
        ALWAYS_VECTORIZE                                                       \
        for ( int j = 0; j < 16; j++ )                                         \
            b.i[j] = ( op a.i[j] );                                            \
        return b;                                                              \
    }

PREFIX_UNARY( +)
PREFIX_UNARY( -)

inline v16int operator!( const v16int& a )
{
    v16int b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.i[j] = -( !a.i[j] );

    return b;
}

PREFIX_UNARY( ~)

#undef PREFIX_UNARY

// v16int prefix increment / decrement

#define PREFIX_INCDEC( op )                                                    \
    inline v16int operator op( v16int& a )                                     \
    {                                                                          \
        v16int b;                                                              \
        ALWAYS_VECTORIZE                                                       \
        for ( int j = 0; j < 16; j++ )                                         \
            b.i[j] = ( op a.i[j] );                                            \
        return b;                                                              \
    }

PREFIX_INCDEC( ++)
PREFIX_INCDEC( --)

#undef PREFIX_INCDEC

// v16int postfix increment / decrement

#define POSTFIX_INCDEC( op )                                                   \
    inline v16int operator op( v16int& a, int )                                \
    {                                                                          \
        v16int b;                                                              \
        ALWAYS_VECTORIZE                                                       \
        for ( int j = 0; j < 16; j++ )                                         \
            b.i[j] = ( a.i[j] op );                                            \
        return b;                                                              \
    }

POSTFIX_INCDEC( ++)
POSTFIX_INCDEC( --)

#undef POSTFIX_INCDEC

// v16int binary operators

#define BINARY( op )                                                           \
    inline v16int operator op( const v16int& a, const v16int& b )              \
    {                                                                          \
        v16int c;                                                              \
        ALWAYS_VECTORIZE                                                       \
        for ( int j = 0; j < 16; j++ )                                         \
            c.i[j] = a.i[j] op b.i[j];                                         \
        return c;                                                              \
    }

BINARY( +)
BINARY( -)
BINARY( * )
BINARY( / )
BINARY( % )
BINARY( ^)
BINARY( & )
BINARY( | )
BINARY( << )
BINARY( >> )

#undef BINARY

// v16int logical operators

#define LOGICAL( op )                                                          \
    inline v16int operator op( const v16int& a, const v16int& b )              \
    {                                                                          \
        v16int c;                                                              \
        ALWAYS_VECTORIZE                                                       \
        for ( int j = 0; j < 16; j++ )                                         \
            c.i[j] = -( a.i[j] op b.i[j] );                                    \
        return c;                                                              \
    }

LOGICAL( < )
LOGICAL( > )
LOGICAL( == )
LOGICAL( != )
LOGICAL( <= )
LOGICAL( >= )
LOGICAL( &&)
LOGICAL( || )

#undef LOGICAL

// v16int miscellaneous functions

inline v16int abs( const v16int& a )
{
    v16int b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.i[j] = ( a.i[j] >= 0 ) ? a.i[j] : -a.i[j];

    return b;
}

inline v16 czero( const v16int& c, const v16& a )
{
    v16 b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.i[j] = a.i[j] & ~c.i[j];

    return b;
}

inline v16 notczero( const v16int& c, const v16& a )
{
    v16 b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.i[j] = a.i[j] & c.i[j];

    return b;
}

inline v16 merge( const v16int& c, const v16& t, const v16& f )
{
    v16 m;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        m.i[j] = ( f.i[j] & ~c.i[j] ) | ( t.i[j] & c.i[j] );

    return m;
}

////////////////
// v16float class

class v16float : public v16
{
    // v16float prefix unary operator friends

    friend inline v16float operator+( const v16float& a ) ALWAYS_INLINE;
    friend inline v16float operator-( const v16float& a ) ALWAYS_INLINE;
    friend inline v16float operator~( const v16float& a ) ALWAYS_INLINE;
    friend inline v16int operator!( const v16float& a ) ALWAYS_INLINE;
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector

    // v16float prefix increment / decrement operator friends

    friend inline v16float operator++( v16float& a ) ALWAYS_INLINE;
    friend inline v16float operator--( v16float& a ) ALWAYS_INLINE;

    // v16float postfix increment / decrement operator friends

    friend inline v16float operator++( v16float& a, int ) ALWAYS_INLINE;
    friend inline v16float operator--( v16float& a, int ) ALWAYS_INLINE;

    // v16float binary operator friends

    friend inline v16float operator+( const v16float& a,
                                      const v16float& b ) ALWAYS_INLINE;
    friend inline v16float operator-( const v16float& a,
                                      const v16float& b ) ALWAYS_INLINE;
    friend inline v16float operator*(const v16float& a,
                                     const v16float& b)ALWAYS_INLINE;
    friend inline v16float operator/( const v16float& a,
                                      const v16float& b ) ALWAYS_INLINE;

    // v16float logical operator friends

    friend inline v16int operator<( const v16float& a,
                                    const v16float& b ) ALWAYS_INLINE;
    friend inline v16int operator>( const v16float& a,
                                    const v16float& b ) ALWAYS_INLINE;
    friend inline v16int operator==( const v16float& a,
                                     const v16float& b ) ALWAYS_INLINE;
    friend inline v16int operator!=( const v16float& a,
                                     const v16float& b ) ALWAYS_INLINE;
    friend inline v16int operator<=( const v16float& a,
                                     const v16float& b ) ALWAYS_INLINE;
    friend inline v16int operator>=( const v16float& a,
                                     const v16float& b ) ALWAYS_INLINE;
    friend inline v16int operator&&( const v16float& a,
                                     const v16float& b ) ALWAYS_INLINE;
    friend inline v16int operator||( const v16float& a,
                                     const v16float& b ) ALWAYS_INLINE;

    // v16float math library friends

#define CMATH_FR1( fn )                                                        \
    friend inline v16float fn( const v16float& a ) ALWAYS_INLINE
#define CMATH_FR2( fn )                                                        \
    friend inline v16float fn( const v16float& a, const v16float& b )          \
        ALWAYS_INLINE

    CMATH_FR1( acos );
    CMATH_FR1( asin );
    CMATH_FR1( atan );
    CMATH_FR2( atan2 );
    CMATH_FR1( ceil );
    CMATH_FR1( cos );
    CMATH_FR1( cosh );
    CMATH_FR1( exp );
    CMATH_FR1( fabs );
    CMATH_FR1( floor );
    CMATH_FR2( fmod );
    CMATH_FR1( log );
    CMATH_FR1( log10 );
    CMATH_FR2( pow );
    CMATH_FR1( sin );
    CMATH_FR1( sinh );
    CMATH_FR1( sqrt );
    CMATH_FR1( tan );
    CMATH_FR1( tanh );

    CMATH_FR2( copysign );

#undef CMATH_FR1
#undef CMATH_FR2

    // v16float miscellaneous friends

    friend inline v16float rsqrt_approx( const v16float& a ) ALWAYS_INLINE;
    friend inline v16float rsqrt( const v16float& a ) ALWAYS_INLINE;
    friend inline v16float rcp_approx( const v16float& a ) ALWAYS_INLINE;
    friend inline v16float rcp( const v16float& a ) ALWAYS_INLINE;
    friend inline v16float fma( const v16float& a, const v16float& b,
                                const v16float& c ) ALWAYS_INLINE;
    friend inline v16float fms( const v16float& a, const v16float& b,
                                const v16float& c ) ALWAYS_INLINE;
    friend inline v16float fnms( const v16float& a, const v16float& b,
                                 const v16float& c ) ALWAYS_INLINE;
    friend inline v16float clear_bits( const v16int& m,
                                       const v16float& a ) ALWAYS_INLINE;
    friend inline v16float set_bits( const v16int& m,
                                     const v16float& a ) ALWAYS_INLINE;
    friend inline v16float toggle_bits( const v16int& m,
                                        const v16float& a ) ALWAYS_INLINE;
    friend inline void increment_16x1( float* ALIGNED( 64 ) p,
                                       const v16float& a ) ALWAYS_INLINE;
    friend inline void decrement_16x1( float* ALIGNED( 64 ) p,
                                       const v16float& a ) ALWAYS_INLINE;
    friend inline void scale_16x1( float* ALIGNED( 64 ) p,
                                   const v16float& a ) ALWAYS_INLINE;

  public:
    // v16float constructors / destructors

    v16float() {} // Default constructor

    v16float( const v16float& a ) // Copy constructor
    {
        ALWAYS_VECTORIZE
        for ( int j = 0; j < 16; j++ )
            f[j] = a.f[j];
    }

    v16float( const v16& a ) // Init from mixed
    {
        ALWAYS_VECTORIZE
        for ( int j = 0; j < 16; j++ )
            f[j] = a.f[j];
    }

    v16float( float a ) // Init from scalar
    {
        ALWAYS_VECTORIZE
        for ( int j = 0; j < 16; j++ )
            f[j] = a;
    }

    v16float( float f00, float f01, float f02, float f03, float f04, float f05,
              float f06, float f07, float f08, float f09, float f10, float f11,
              float f12, float f13, float f14, float f15 ) // Init from scalars
    {
        f[0] = f00;
        f[1] = f01;
        f[2] = f02;
        f[3] = f03;
        f[4] = f04;
        f[5] = f05;
        f[6] = f06;
        f[7] = f07;
        f[8] = f08;
        f[9] = f09;
        f[10] = f10;
        f[11] = f11;
        f[12] = f12;
        f[13] = f13;
        f[14] = f14;
        f[15] = f15;
    }

    ~v16float() {} // Destructor

    // v16float assignment operators

#define ASSIGN( op )                                                           \
    inline v16float& operator op( const v16float& b )                          \
    {                                                                          \
        ALWAYS_VECTORIZE                                                       \
        for ( int j = 0; j < 16; j++ )                                         \
            f[j] op b.f[j];                                                    \
        return *this;                                                          \
    }

    ASSIGN( = )
    ASSIGN( += )
    ASSIGN( -= )
    ASSIGN( *= )
    ASSIGN( /= )

#undef ASSIGN

    // v16float member access operator

    inline float& operator[]( int n ) { return f[n]; }

    inline float operator()( int n ) { return f[n]; }
};

// v16float prefix unary operators

inline v16float operator+( const v16float& a )
{
    v16float b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.f[j] = +a.f[j];

    return b;
}

inline v16float operator-( const v16float& a )
{
    v16float b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.f[j] = -a.f[j];

    return b;
}

inline v16int operator!( const v16float& a )
{
    v16int b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.i[j] = a.i[j] ? 0 : -1;

    return b;
}

// v16float prefix increment / decrement operators

inline v16float operator++( v16float& a )
{
    v16float b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.f[j] = ++a.f[j];

    return b;
}

inline v16float operator--( v16float& a )
{
    v16float b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.f[j] = --a.f[j];

    return b;
}

// v16float postfix increment / decrement operators

inline v16float operator++( v16float& a, int )
{
    v16float b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.f[j] = a.f[j]++;

    return b;
}

inline v16float operator--( v16float& a, int )
{
    v16float b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.f[j] = a.f[j]--;

    return b;
}

// v16float binary operators

#define BINARY( op )                                                           \
    inline v16float operator op( const v16float& a, const v16float& b )        \
    {                                                                          \
        v16float c;                                                            \
        ALWAYS_VECTORIZE                                                       \
        for ( int j = 0; j < 16; j++ )                                         \
            c.f[j] = a.f[j] op b.f[j];                                         \
        return c;                                                              \
    }

BINARY( +)
BINARY( -)
BINARY( * )
BINARY( / )

#undef BINARY

// v16float logical operators

#define LOGICAL( op )                                                          \
    inline v16int operator op( const v16float& a, const v16float& b )          \
    {                                                                          \
        v16int c;                                                              \
        ALWAYS_VECTORIZE                                                       \
        for ( int j = 0; j < 16; j++ )                                         \
            c.i[j] = -( a.f[j] op b.f[j] );                                    \
        return c;                                                              \
    }

LOGICAL( < )
LOGICAL( > )
LOGICAL( == )
LOGICAL( != )
LOGICAL( <= )
LOGICAL( >= )
LOGICAL( &&)
LOGICAL( || )

#undef LOGICAL

// v16float math library functions

#define CMATH_FR1( fn )                                                        \
    inline v16float fn( const v16float& a )                                    \
    {                                                                          \
        v16float b;                                                            \
        ALWAYS_VECTORIZE                                                       \
        for ( int j = 0; j < 16; j++ )                                         \
            b.f[j] = ::fn( a.f[j] );                                           \
        return b;                                                              \
    }

#define CMATH_FR2( fn )                                                        \
    inline v16float fn( const v16float& a, const v16float& b )                 \
    {                                                                          \
        v16float c;                                                            \
        ALWAYS_VECTORIZE                                                       \
        for ( int j = 0; j < 16; j++ )                                         \
            c.f[j] = ::fn( a.f[j], b.f[j] );                                   \
        return c;                                                              \
    }

CMATH_FR1( acos )
CMATH_FR1( asin ) CMATH_FR1( atan ) CMATH_FR2( atan2 ) CMATH_FR1( ceil )
    CMATH_FR1( cos ) CMATH_FR1( cosh ) CMATH_FR1( exp ) CMATH_FR1( fabs )
        CMATH_FR1( floor ) CMATH_FR2( fmod ) CMATH_FR1( log ) CMATH_FR1( log10 )
            CMATH_FR2( pow ) CMATH_FR1( sin ) CMATH_FR1( sinh )
                CMATH_FR1( sqrt ) CMATH_FR1( tan ) CMATH_FR1( tanh )

                    inline v16float
    copysign( const v16float& a, const v16float& b )
{
    v16float c;
    float t;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
    {
        t = ::fabs( a.f[j] );
        if ( b.f[j] < 0 )
            t = -t;
        c.f[j] = t;
    }

    return c;
}

#undef CMATH_FR1
#undef CMATH_FR2

// v16float miscellaneous functions

inline v16float rsqrt_approx( const v16float& a )
{
    v16float b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.f[j] = ::sqrt( 1.0f / a.f[j] );

    return b;
}

inline v16float rsqrt( const v16float& a )
{
    v16float b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.f[j] = ::sqrt( 1.0f / a.f[j] );

    return b;
}

inline v16float rcp_approx( const v16float& a )
{
    v16float b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.f[j] = 1.0f / a.f[j];

    return b;
}

inline v16float rcp( const v16float& a )
{
    v16float b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.f[j] = 1.0f / a.f[j];

    return b;
}

inline v16float fma( const v16float& a, const v16float& b, const v16float& c )
{
    v16float d;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        d.f[j] = a.f[j] * b.f[j] + c.f[j];

    return d;
}

inline v16float fms( const v16float& a, const v16float& b, const v16float& c )
{
    v16float d;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        d.f[j] = a.f[j] * b.f[j] - c.f[j];

    return d;
}

inline v16float fnms( const v16float& a, const v16float& b, const v16float& c )
{
    v16float d;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        d.f[j] = c.f[j] - a.f[j] * b.f[j];

    return d;
}

inline v16float clear_bits( const v16int& m, const v16float& a )
{
    v16float b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.i[j] = ( ~m.i[j] ) & a.i[j];

    return b;
}

inline v16float set_bits( const v16int& m, const v16float& a )
{
    v16float b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.i[j] = m.i[j] | a.i[j];

    return b;
}

inline v16float toggle_bits( const v16int& m, const v16float& a )
{
    v16float b;

    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        b.i[j] = m.i[j] ^ a.i[j];

    return b;
}

inline void increment_16x1( float* ALIGNED( 64 ) p, const v16float& a )
{
    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        p[j] += a.f[j];
}

inline void decrement_16x1( float* ALIGNED( 64 ) p, const v16float& a )
{
    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        p[j] -= a.f[j];
}

inline void scale_16x1( float* ALIGNED( 64 ) p, const v16float& a )
{
    ALWAYS_VECTORIZE
    for ( int j = 0; j < 16; j++ )
        p[j] *= a.f[j];
}

} // namespace v16

#endif // _v16_portable_h_
