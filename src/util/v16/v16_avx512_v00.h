#ifndef _v16_avx512_h_
#define _v16_avx512_h_

#ifndef IN_v16_h
#error "Do not include v16_avx512.h directly; use v16.h"
#endif

#define V16_ACCELERATION
#define V16_AVX512_ACCELERATION

#include <immintrin.h>
#include <math.h>

#ifndef ALIGNED
#define ALIGNED(n)
#endif

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

    friend inline int any( const v16 &a );
    friend inline int all( const v16 &a );

    template<int n>
    friend inline v16 splat( const v16 &a );

    template<int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, int i10, int i11, int i12, int i13, int i14, int i15>
    friend inline v16 shuffle( const v16 &a );

    friend inline void swap( v16 &a, v16 &b );
    friend inline void transpose( v16 &a00, v16 &a01, v16 &a02, v16 &a03,
				  v16 &a04, v16 &a05, v16 &a06, v16 &a07,
				  v16 &a08, v16 &a09, v16 &a10, v16 &a11,
				  v16 &a12, v16 &a13, v16 &a14, v16 &a15 );

    // v16int miscellaneous friends

    friend inline v16    czero( const v16int &c, const v16 &a );
    friend inline v16 notczero( const v16int &c, const v16 &a );
    friend inline v16    merge( const v16int &c, const v16 &a, const v16 &b );

    // v16 memory manipulation friends

    friend inline void   load_16x1( const void * ALIGNED(64) p, v16 &a );
    friend inline void  store_16x1( const v16 &a, void * ALIGNED(64) p );
    friend inline void stream_16x1( const v16 &a, void * ALIGNED(64) p );
    friend inline void  clear_16x1( void * ALIGNED(64) dst );
    friend inline void   copy_16x1( void * ALIGNED(64) dst,
				    const void * ALIGNED(64) src );
    friend inline void   swap_16x1( void * ALIGNED(64) a, void * ALIGNED(64) b );

    // v16 transposed memory manipulation friends
    // Note: Half aligned values are permissible in the 16x2_tr variants.

    friend inline void load_16x1_tr( const void *a00, const void *a01,
				     const void *a02, const void *a03,
				     const void *a04, const void *a05,
				     const void *a06, const void *a07,
				     const void *a08, const void *a09,
				     const void *a10, const void *a11,
				     const void *a12, const void *a13,
				     const void *a14, const void *a15,
				     v16 &a );
    friend inline void load_16x2_tr( const void * ALIGNED(8) a00,
				     const void * ALIGNED(8) a01,
				     const void * ALIGNED(8) a02,
				     const void * ALIGNED(8) a03,
				     const void * ALIGNED(8) a04,
				     const void * ALIGNED(8) a05,
				     const void * ALIGNED(8) a06,
				     const void * ALIGNED(8) a07,
				     const void * ALIGNED(8) a08,
				     const void * ALIGNED(8) a09,
				     const void * ALIGNED(8) a10,
				     const void * ALIGNED(8) a11,
				     const void * ALIGNED(8) a12,
				     const void * ALIGNED(8) a13,
				     const void * ALIGNED(8) a14,
				     const void * ALIGNED(8) a15,
				     v16 &a, v16 &b );
    friend inline void load_16x3_tr( const void * ALIGNED(64) a00,
				     const void * ALIGNED(64) a01,
				     const void * ALIGNED(64) a02,
				     const void * ALIGNED(64) a03,
				     const void * ALIGNED(64) a04,
				     const void * ALIGNED(64) a05,
				     const void * ALIGNED(64) a06,
				     const void * ALIGNED(64) a07,
				     const void * ALIGNED(64) a08,
				     const void * ALIGNED(64) a09,
				     const void * ALIGNED(64) a10,
				     const void * ALIGNED(64) a11,
				     const void * ALIGNED(64) a12,
				     const void * ALIGNED(64) a13,
				     const void * ALIGNED(64) a14,
				     const void * ALIGNED(64) a15,
				     v16 &a, v16 &b, v16 &c );
    friend inline void load_16x4_tr( const void * ALIGNED(64) a00,
				     const void * ALIGNED(64) a01,
				     const void * ALIGNED(64) a02,
				     const void * ALIGNED(64) a03,
				     const void * ALIGNED(64) a04,
				     const void * ALIGNED(64) a05,
				     const void * ALIGNED(64) a06,
				     const void * ALIGNED(64) a07,
				     const void * ALIGNED(64) a08,
				     const void * ALIGNED(64) a09,
				     const void * ALIGNED(64) a10,
				     const void * ALIGNED(64) a11,
				     const void * ALIGNED(64) a12,
				     const void * ALIGNED(64) a13,
				     const void * ALIGNED(64) a14,
				     const void * ALIGNED(64) a15,
				     v16 &a, v16 &b, v16 &c, v16 &d );
    friend inline void load_16x8_tr( const void * ALIGNED(64) a00,
				     const void * ALIGNED(64) a01,
				     const void * ALIGNED(64) a02,
				     const void * ALIGNED(64) a03,
				     const void * ALIGNED(64) a04,
				     const void * ALIGNED(64) a05,
				     const void * ALIGNED(64) a06,
				     const void * ALIGNED(64) a07,
				     const void * ALIGNED(64) a08,
				     const void * ALIGNED(64) a09,
				     const void * ALIGNED(64) a10,
				     const void * ALIGNED(64) a11,
				     const void * ALIGNED(64) a12,
				     const void * ALIGNED(64) a13,
				     const void * ALIGNED(64) a14,
				     const void * ALIGNED(64) a15,
				     v16 &a, v16 &b, v16 &c, v16 &d,
				     v16 &e, v16 &f, v16 &g, v16 &h );
    friend inline void load_16x16_tr( const void * ALIGNED(64) a00,
				      const void * ALIGNED(64) a01,
				      const void * ALIGNED(64) a02,
				      const void * ALIGNED(64) a03,
				      const void * ALIGNED(64) a04,
				      const void * ALIGNED(64) a05,
				      const void * ALIGNED(64) a06,
				      const void * ALIGNED(64) a07,
				      const void * ALIGNED(64) a08,
				      const void * ALIGNED(64) a09,
				      const void * ALIGNED(64) a10,
				      const void * ALIGNED(64) a11,
				      const void * ALIGNED(64) a12,
				      const void * ALIGNED(64) a13,
				      const void * ALIGNED(64) a14,
				      const void * ALIGNED(64) a15,
				      v16 &b00, v16 &b01, v16 &b02, v16 &b03,
				      v16 &b04, v16 &b05, v16 &b06, v16 &b07,
				      v16 &b08, v16 &b09, v16 &b10, v16 &b11,
				      v16 &b12, v16 &b13, v16 &b14, v16 &b15 );
    friend inline void load_16x16_tr_p( const void * ALIGNED(64) a00,
					const void * ALIGNED(64) a01,
					const void * ALIGNED(64) a02,
					const void * ALIGNED(64) a03,
					const void * ALIGNED(64) a04,
					const void * ALIGNED(64) a05,
					const void * ALIGNED(64) a06,
					const void * ALIGNED(64) a07,
					const void * ALIGNED(64) a08,
					const void * ALIGNED(64) a09,
					const void * ALIGNED(64) a10,
					const void * ALIGNED(64) a11,
					const void * ALIGNED(64) a12,
					const void * ALIGNED(64) a13,
					const void * ALIGNED(64) a14,
					const void * ALIGNED(64) a15,
					v16 &b00, v16 &b01, v16 &b02, v16 &b03,
					v16 &b04, v16 &b05, v16 &b06, v16 &b07,
					v16 &b08, v16 &b09, v16 &b10, v16 &b11,
					v16 &b12, v16 &b13, v16 &b14, v16 &b15 );

    friend inline void store_16x1_tr( const v16 &a,
				      void *a00, void *a01, void *a02, void *a03,
				      void *a04, void *a05, void *a06, void *a07,
				      void *a08, void *a09, void *a10, void *a11,
				      void *a12, void *a13, void *a14, void *a15 );
    friend inline void store_16x2_tr( const v16 &a, const v16 &b,
				      void * ALIGNED(8) a00,
				      void * ALIGNED(8) a01,
				      void * ALIGNED(8) a02,
				      void * ALIGNED(8) a03,
				      void * ALIGNED(8) a04,
				      void * ALIGNED(8) a05,
				      void * ALIGNED(8) a06,
				      void * ALIGNED(8) a07,
				      void * ALIGNED(8) a08,
				      void * ALIGNED(8) a09,
				      void * ALIGNED(8) a10,
				      void * ALIGNED(8) a11,
				      void * ALIGNED(8) a12,
				      void * ALIGNED(8) a13,
				      void * ALIGNED(8) a14,
				      void * ALIGNED(8) a15 );
    friend inline void store_16x3_tr( const v16 &a, const v16 &b, const v16 &c,
				      void * ALIGNED(64) a00,
				      void * ALIGNED(64) a01,
				      void * ALIGNED(64) a02,
				      void * ALIGNED(64) a03,
				      void * ALIGNED(64) a04,
				      void * ALIGNED(64) a05,
				      void * ALIGNED(64) a06,
				      void * ALIGNED(64) a07,
				      void * ALIGNED(64) a08,
				      void * ALIGNED(64) a09,
				      void * ALIGNED(64) a10,
				      void * ALIGNED(64) a11,
				      void * ALIGNED(64) a12,
				      void * ALIGNED(64) a13,
				      void * ALIGNED(64) a14,
				      void * ALIGNED(64) a15 );
    friend inline void store_16x4_tr( const v16 &a, const v16 &b,
				      const v16 &c, const v16 &d,
				      void * ALIGNED(64) a00,
				      void * ALIGNED(64) a01,
				      void * ALIGNED(64) a02,
				      void * ALIGNED(64) a03,
				      void * ALIGNED(64) a04,
				      void * ALIGNED(64) a05,
				      void * ALIGNED(64) a06,
				      void * ALIGNED(64) a07,
				      void * ALIGNED(64) a08,
				      void * ALIGNED(64) a09,
				      void * ALIGNED(64) a10,
				      void * ALIGNED(64) a11,
				      void * ALIGNED(64) a12,
				      void * ALIGNED(64) a13,
				      void * ALIGNED(64) a14,
				      void * ALIGNED(64) a15 );
    friend inline void store_16x8_tr( const v16 &a, const v16 &b,
				      const v16 &c, const v16 &d,
				      const v16 &e, const v16 &f,
				      const v16 &g, const v16 &h,
				      void * ALIGNED(64) a00,
				      void * ALIGNED(64) a01,
				      void * ALIGNED(64) a02,
				      void * ALIGNED(64) a03,
				      void * ALIGNED(64) a04,
				      void * ALIGNED(64) a05,
				      void * ALIGNED(64) a06,
				      void * ALIGNED(64) a07,
				      void * ALIGNED(64) a08,
				      void * ALIGNED(64) a09,
				      void * ALIGNED(64) a10,
				      void * ALIGNED(64) a11,
				      void * ALIGNED(64) a12,
				      void * ALIGNED(64) a13,
				      void * ALIGNED(64) a14,
				      void * ALIGNED(64) a15 );
    friend inline void store_16x16_tr( const v16 &b00, const v16 &b01,
				       const v16 &b02, const v16 &b03,
				       const v16 &b04, const v16 &b05,
				       const v16 &b06, const v16 &b07,
				       const v16 &b08, const v16 &b09,
				       const v16 &b10, const v16 &b11,
				       const v16 &b12, const v16 &b13,
				       const v16 &b14, const v16 &b15,
				       void * ALIGNED(64) a00,
				       void * ALIGNED(64) a01,
				       void * ALIGNED(64) a02,
				       void * ALIGNED(64) a03,
				       void * ALIGNED(64) a04,
				       void * ALIGNED(64) a05,
				       void * ALIGNED(64) a06,
				       void * ALIGNED(64) a07,
				       void * ALIGNED(64) a08,
				       void * ALIGNED(64) a09,
				       void * ALIGNED(64) a10,
				       void * ALIGNED(64) a11,
				       void * ALIGNED(64) a12,
				       void * ALIGNED(64) a13,
				       void * ALIGNED(64) a14,
				       void * ALIGNED(64) a15 );
    friend inline void store_16x16_tr_p( const v16 &b00, const v16 &b01,
					 const v16 &b02, const v16 &b03,
					 const v16 &b04, const v16 &b05,
					 const v16 &b06, const v16 &b07,
					 const v16 &b08, const v16 &b09,
					 const v16 &b10, const v16 &b11,
					 const v16 &b12, const v16 &b13,
					 const v16 &b14, const v16 &b15,
					 void * ALIGNED(64) a00,
					 void * ALIGNED(64) a01,
					 void * ALIGNED(64) a02,
					 void * ALIGNED(64) a03,
					 void * ALIGNED(64) a04,
					 void * ALIGNED(64) a05,
					 void * ALIGNED(64) a06,
					 void * ALIGNED(64) a07,
					 void * ALIGNED(64) a08,
					 void * ALIGNED(64) a09,
					 void * ALIGNED(64) a10,
					 void * ALIGNED(64) a11,
					 void * ALIGNED(64) a12,
					 void * ALIGNED(64) a13,
					 void * ALIGNED(64) a14,
					 void * ALIGNED(64) a15 );

  protected:

    union
    {
      int   i[16];
      float f[16];
      __m512 v;
    };

  public:

    v16() {}                    // Default constructor

    v16( const v16 &a )         // Copy constructor
    {
      i[ 0]=a.i[ 0]; i[ 1]=a.i[ 1]; i[ 2]=a.i[ 2]; i[ 3]=a.i[ 3];
      i[ 4]=a.i[ 4]; i[ 5]=a.i[ 5]; i[ 6]=a.i[ 6]; i[ 7]=a.i[ 7];
      i[ 8]=a.i[ 8]; i[ 9]=a.i[ 9]; i[10]=a.i[10]; i[11]=a.i[11];
      i[12]=a.i[12]; i[13]=a.i[13]; i[14]=a.i[14]; i[15]=a.i[15];
    }

    ~v16() {}                   // Default destructor
  };

  // v16 miscellaneous functions

  inline int any( const v16 &a )
  {
    return a.i[ 0] || a.i[ 1] || a.i[ 2] || a.i[ 3] ||
           a.i[ 4] || a.i[ 5] || a.i[ 6] || a.i[ 7] ||
           a.i[ 8] || a.i[ 9] || a.i[10] || a.i[11] ||
           a.i[12] || a.i[13] || a.i[14] || a.i[15];
  }

  inline int all( const v16 &a )
  {
    return a.i[ 0] && a.i[ 1] && a.i[ 2] && a.i[ 3] &&
           a.i[ 4] && a.i[ 5] && a.i[ 6] && a.i[ 7] &&
           a.i[ 8] && a.i[ 9] && a.i[10] && a.i[11] &&
           a.i[12] && a.i[13] && a.i[14] && a.i[15];
  }

  template<int n>
  inline v16 splat( const v16 & a )
  {
    v16 b;

    b.i[ 0] = a.i[n];
    b.i[ 1] = a.i[n];
    b.i[ 2] = a.i[n];
    b.i[ 3] = a.i[n];
    b.i[ 4] = a.i[n];
    b.i[ 5] = a.i[n];
    b.i[ 6] = a.i[n];
    b.i[ 7] = a.i[n];
    b.i[ 8] = a.i[n];
    b.i[ 9] = a.i[n];
    b.i[10] = a.i[n];
    b.i[11] = a.i[n];
    b.i[12] = a.i[n];
    b.i[13] = a.i[n];
    b.i[14] = a.i[n];
    b.i[15] = a.i[n];

    return b;
  }

  template<int i00, int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11, int i12, int i13, int i14, int i16>
  inline v16 shuffle( const v16 & a )
  {
    v16 b;

    b.i[ 0] = a.i[i00];
    b.i[ 1] = a.i[i01];
    b.i[ 2] = a.i[i02];
    b.i[ 3] = a.i[i03];
    b.i[ 4] = a.i[i04];
    b.i[ 5] = a.i[i05];
    b.i[ 6] = a.i[i06];
    b.i[ 7] = a.i[i07];
    b.i[ 8] = a.i[i08];
    b.i[ 9] = a.i[i09];
    b.i[10] = a.i[i10];
    b.i[11] = a.i[i11];
    b.i[12] = a.i[i12];
    b.i[13] = a.i[i13];
    b.i[14] = a.i[i14];
    b.i[15] = a.i[i15];

    return b;
  }

# define sw(x,y) x^=y, y^=x, x^=y

  inline void swap( v16 &a, v16 &b )
  {
    sw( a.i[ 0], b.i[ 0] );
    sw( a.i[ 1], b.i[ 1] );
    sw( a.i[ 2], b.i[ 2] );
    sw( a.i[ 3], b.i[ 3] );
    sw( a.i[ 4], b.i[ 4] );
    sw( a.i[ 5], b.i[ 5] );
    sw( a.i[ 6], b.i[ 6] );
    sw( a.i[ 7], b.i[ 7] );
    sw( a.i[ 8], b.i[ 8] );
    sw( a.i[ 9], b.i[ 9] );
    sw( a.i[10], b.i[10] );
    sw( a.i[11], b.i[11] );
    sw( a.i[12], b.i[12] );
    sw( a.i[13], b.i[13] );
    sw( a.i[14], b.i[14] );
    sw( a.i[15], b.i[15] );
  }

  inline void transpose( v16 &a00, v16 &a01, v16 &a02, v16 &a03,
			 v16 &a04, v16 &a05, v16 &a06, v16 &a07,
			 v16 &a08, v16 &a09, v16 &a10, v16 &a11,
			 v16 &a12, v16 &a13, v16 &a14, v16 &a15 )
  {
    sw( a00.i[1],a01.i[0] ); sw( a00.i[2],a02.i[0] ); sw( a00.i[3],a03.i[0] ); sw( a00.i[4],a04.i[0] ); sw( a00.i[5],a05.i[0] ); sw( a00.i[6],a06.i[0] ); sw( a00.i[7],a07.i[0] ); sw( a00.i[8],a08.i[0] ); sw( a00.i[9],a09.i[0] ); sw( a00.i[10],a10.i[0] ); sw( a00.i[11],a11.i[ 0] ); sw( a00.i[12],a12.i[ 0] ); sw( a00.i[13],a13.i[ 0] ); sw( a00.i[14],a14.i[ 0] ); sw( a00.i[15],a15.i[ 0] );
                             sw( a01.i[2],a02.i[1] ); sw( a01.i[3],a03.i[1] ); sw( a01.i[4],a04.i[1] ); sw( a01.i[5],a05.i[1] ); sw( a01.i[6],a06.i[1] ); sw( a01.i[7],a07.i[1] ); sw( a01.i[8],a08.i[1] ); sw( a01.i[9],a09.i[1] ); sw( a01.i[10],a10.i[1] ); sw( a01.i[11],a11.i[ 1] ); sw( a01.i[12],a12.i[ 1] ); sw( a01.i[13],a13.i[ 1] ); sw( a01.i[14],a14.i[ 1] ); sw( a01.i[15],a15.i[ 1] );
                                                      sw( a02.i[3],a03.i[2] ); sw( a02.i[4],a04.i[2] ); sw( a02.i[5],a05.i[2] ); sw( a02.i[6],a06.i[2] ); sw( a02.i[7],a07.i[2] ); sw( a02.i[8],a08.i[2] ); sw( a02.i[9],a09.i[2] ); sw( a02.i[10],a10.i[2] ); sw( a02.i[11],a11.i[ 2] ); sw( a02.i[12],a12.i[ 2] ); sw( a02.i[13],a13.i[ 2] ); sw( a02.i[14],a14.i[ 2] ); sw( a02.i[15],a15.i[ 2] );
                                                                               sw( a03.i[4],a04.i[3] ); sw( a03.i[5],a05.i[3] ); sw( a03.i[6],a06.i[3] ); sw( a03.i[7],a07.i[3] ); sw( a03.i[8],a08.i[3] ); sw( a03.i[9],a09.i[3] ); sw( a03.i[10],a10.i[3] ); sw( a03.i[11],a11.i[ 3] ); sw( a03.i[12],a12.i[ 3] ); sw( a03.i[13],a13.i[ 3] ); sw( a03.i[14],a14.i[ 3] ); sw( a03.i[15],a15.i[ 3] );
                                                                                                        sw( a04.i[5],a05.i[4] ); sw( a04.i[6],a06.i[4] ); sw( a04.i[7],a07.i[4] ); sw( a04.i[8],a08.i[4] ); sw( a04.i[9],a09.i[4] ); sw( a04.i[10],a10.i[4] ); sw( a04.i[11],a11.i[ 4] ); sw( a04.i[12],a12.i[ 4] ); sw( a04.i[13],a13.i[ 4] ); sw( a04.i[14],a14.i[ 4] ); sw( a04.i[15],a15.i[ 4] );
                                                                                                                                 sw( a05.i[6],a06.i[5] ); sw( a05.i[7],a07.i[5] ); sw( a05.i[8],a08.i[5] ); sw( a05.i[9],a09.i[5] ); sw( a05.i[10],a10.i[5] ); sw( a05.i[11],a11.i[ 5] ); sw( a05.i[12],a12.i[ 5] ); sw( a05.i[13],a13.i[ 5] ); sw( a05.i[14],a14.i[ 5] ); sw( a05.i[15],a15.i[ 5] );
                                                                                                                                                          sw( a06.i[7],a07.i[6] ); sw( a06.i[8],a08.i[6] ); sw( a06.i[9],a09.i[6] ); sw( a06.i[10],a10.i[6] ); sw( a06.i[11],a11.i[ 6] ); sw( a06.i[12],a12.i[ 6] ); sw( a06.i[13],a13.i[ 6] ); sw( a06.i[14],a14.i[ 6] ); sw( a06.i[15],a15.i[ 6] );
                                                                                                                                                                                   sw( a07.i[8],a08.i[7] ); sw( a07.i[9],a09.i[7] ); sw( a07.i[10],a10.i[7] ); sw( a07.i[11],a11.i[ 7] ); sw( a07.i[12],a12.i[ 7] ); sw( a07.i[13],a13.i[ 7] ); sw( a07.i[14],a14.i[ 7] ); sw( a07.i[15],a15.i[ 7] );
                                                                                                                                                                                                            sw( a08.i[9],a09.i[8] ); sw( a08.i[10],a10.i[8] ); sw( a08.i[11],a11.i[ 8] ); sw( a08.i[12],a12.i[ 8] ); sw( a08.i[13],a13.i[ 8] ); sw( a08.i[14],a14.i[ 8] ); sw( a08.i[15],a15.i[ 8] );
                                                                                                                                                                                                                                     sw( a09.i[10],a10.i[9] ); sw( a09.i[11],a11.i[ 9] ); sw( a09.i[12],a12.i[ 9] ); sw( a09.i[13],a13.i[ 9] ); sw( a09.i[14],a14.i[ 9] ); sw( a09.i[15],a15.i[ 9] );
                                                                                                                                                                                                                                                               sw( a10.i[11],a11.i[10] ); sw( a10.i[12],a12.i[10] ); sw( a10.i[13],a13.i[10] ); sw( a10.i[14],a14.i[10] ); sw( a10.i[15],a15.i[10] );
                                                                                                                                                                                                                                                                                          sw( a11.i[12],a12.i[11] ); sw( a11.i[13],a13.i[11] ); sw( a11.i[14],a14.i[11] ); sw( a11.i[15],a15.i[11] );
                                                                                                                                                                                                                                                                                                                     sw( a12.i[13],a13.i[12] ); sw( a12.i[14],a14.i[12] ); sw( a12.i[15],a15.i[12] );
                                                                                                                                                                                                                                                                                                                                                sw( a13.i[14],a14.i[13] ); sw( a13.i[15],a15.i[13] );
                                                                                                                                                                                                                                                                                                                                                                           sw( a14.i[15],a15.i[14] );
  }

# undef sw

  // v16 memory manipulation functions

  inline void load_16x1( const void * ALIGNED(64) p,
			 v16 &a )
  {
    a.i[ 0] = ((const int * ALIGNED(64))p)[ 0];
    a.i[ 1] = ((const int * ALIGNED(64))p)[ 1];
    a.i[ 2] = ((const int * ALIGNED(64))p)[ 2];
    a.i[ 3] = ((const int * ALIGNED(64))p)[ 3];
    a.i[ 4] = ((const int * ALIGNED(64))p)[ 4];
    a.i[ 5] = ((const int * ALIGNED(64))p)[ 5];
    a.i[ 6] = ((const int * ALIGNED(64))p)[ 6];
    a.i[ 7] = ((const int * ALIGNED(64))p)[ 7];
    a.i[ 8] = ((const int * ALIGNED(64))p)[ 8];
    a.i[ 9] = ((const int * ALIGNED(64))p)[ 9];
    a.i[10] = ((const int * ALIGNED(64))p)[10];
    a.i[11] = ((const int * ALIGNED(64))p)[11];
    a.i[12] = ((const int * ALIGNED(64))p)[12];
    a.i[13] = ((const int * ALIGNED(64))p)[13];
    a.i[14] = ((const int * ALIGNED(64))p)[14];
    a.i[15] = ((const int * ALIGNED(64))p)[15];
  }

  inline void store_16x1( const v16 &a,
			  void * ALIGNED(64) p )
  {
    ((int * ALIGNED(64))p)[ 0] = a.i[ 0];
    ((int * ALIGNED(64))p)[ 1] = a.i[ 1];
    ((int * ALIGNED(64))p)[ 2] = a.i[ 2];
    ((int * ALIGNED(64))p)[ 3] = a.i[ 3];
    ((int * ALIGNED(64))p)[ 4] = a.i[ 4];
    ((int * ALIGNED(64))p)[ 5] = a.i[ 5];
    ((int * ALIGNED(64))p)[ 6] = a.i[ 6];
    ((int * ALIGNED(64))p)[ 7] = a.i[ 7];
    ((int * ALIGNED(64))p)[ 8] = a.i[ 8];
    ((int * ALIGNED(64))p)[ 9] = a.i[ 9];
    ((int * ALIGNED(64))p)[10] = a.i[10];
    ((int * ALIGNED(64))p)[11] = a.i[11];
    ((int * ALIGNED(64))p)[12] = a.i[12];
    ((int * ALIGNED(64))p)[13] = a.i[13];
    ((int * ALIGNED(64))p)[14] = a.i[14];
    ((int * ALIGNED(64))p)[15] = a.i[15];
  }

  inline void stream_16x1( const v16 &a,
			   void * ALIGNED(64) p )
  {
    ((int * ALIGNED(64))p)[ 0] = a.i[ 0];
    ((int * ALIGNED(64))p)[ 1] = a.i[ 1];
    ((int * ALIGNED(64))p)[ 2] = a.i[ 2];
    ((int * ALIGNED(64))p)[ 3] = a.i[ 3];
    ((int * ALIGNED(64))p)[ 4] = a.i[ 4];
    ((int * ALIGNED(64))p)[ 5] = a.i[ 5];
    ((int * ALIGNED(64))p)[ 6] = a.i[ 6];
    ((int * ALIGNED(64))p)[ 7] = a.i[ 7];
    ((int * ALIGNED(64))p)[ 8] = a.i[ 8];
    ((int * ALIGNED(64))p)[ 9] = a.i[ 9];
    ((int * ALIGNED(64))p)[10] = a.i[10];
    ((int * ALIGNED(64))p)[11] = a.i[11];
    ((int * ALIGNED(64))p)[12] = a.i[12];
    ((int * ALIGNED(64))p)[13] = a.i[13];
    ((int * ALIGNED(64))p)[14] = a.i[14];
    ((int * ALIGNED(64))p)[15] = a.i[15];
  }

  inline void clear_16x1( void * ALIGNED(64) p )
  {
    ((int * ALIGNED(64))p)[ 0] = 0;
    ((int * ALIGNED(64))p)[ 1] = 0;
    ((int * ALIGNED(64))p)[ 2] = 0;
    ((int * ALIGNED(64))p)[ 3] = 0;
    ((int * ALIGNED(64))p)[ 4] = 0;
    ((int * ALIGNED(64))p)[ 5] = 0;
    ((int * ALIGNED(64))p)[ 6] = 0;
    ((int * ALIGNED(64))p)[ 7] = 0;
    ((int * ALIGNED(64))p)[ 8] = 0;
    ((int * ALIGNED(64))p)[ 9] = 0;
    ((int * ALIGNED(64))p)[10] = 0;
    ((int * ALIGNED(64))p)[11] = 0;
    ((int * ALIGNED(64))p)[12] = 0;
    ((int * ALIGNED(64))p)[13] = 0;
    ((int * ALIGNED(64))p)[14] = 0;
    ((int * ALIGNED(64))p)[15] = 0;
  }

  // FIXME: Ordering semantics
  inline void copy_16x1( void * ALIGNED(64) dst,
			 const void * ALIGNED(64) src )
  {
    ((int * ALIGNED(64))dst)[ 0] = ((const int * ALIGNED(64))src)[ 0];
    ((int * ALIGNED(64))dst)[ 1] = ((const int * ALIGNED(64))src)[ 1];
    ((int * ALIGNED(64))dst)[ 2] = ((const int * ALIGNED(64))src)[ 2];
    ((int * ALIGNED(64))dst)[ 3] = ((const int * ALIGNED(64))src)[ 3];
    ((int * ALIGNED(64))dst)[ 4] = ((const int * ALIGNED(64))src)[ 4];
    ((int * ALIGNED(64))dst)[ 5] = ((const int * ALIGNED(64))src)[ 5];
    ((int * ALIGNED(64))dst)[ 6] = ((const int * ALIGNED(64))src)[ 6];
    ((int * ALIGNED(64))dst)[ 7] = ((const int * ALIGNED(64))src)[ 7];
    ((int * ALIGNED(64))dst)[ 8] = ((const int * ALIGNED(64))src)[ 8];
    ((int * ALIGNED(64))dst)[ 9] = ((const int * ALIGNED(64))src)[ 9];
    ((int * ALIGNED(64))dst)[10] = ((const int * ALIGNED(64))src)[10];
    ((int * ALIGNED(64))dst)[11] = ((const int * ALIGNED(64))src)[11];
    ((int * ALIGNED(64))dst)[12] = ((const int * ALIGNED(64))src)[12];
    ((int * ALIGNED(64))dst)[13] = ((const int * ALIGNED(64))src)[13];
    ((int * ALIGNED(64))dst)[14] = ((const int * ALIGNED(64))src)[14];
    ((int * ALIGNED(64))dst)[15] = ((const int * ALIGNED(64))src)[15];
  }

  inline void swap_16x1( void * ALIGNED(64) a,
			 void * ALIGNED(64) b )
  {
    int t;

    t = ((int * ALIGNED(64))a)[ 0];
    ((int * ALIGNED(64))a)[ 0] = ((int * ALIGNED(64))b)[ 0];
    ((int * ALIGNED(64))b)[ 0] = t;

    t = ((int * ALIGNED(64))a)[ 1];
    ((int * ALIGNED(64))a)[ 1] = ((int * ALIGNED(64))b)[ 1];
    ((int * ALIGNED(64))b)[ 1] = t;

    t = ((int * ALIGNED(64))a)[ 2];
    ((int * ALIGNED(64))a)[ 2] = ((int * ALIGNED(64))b)[ 2];
    ((int * ALIGNED(64))b)[ 2] = t;

    t = ((int * ALIGNED(64))a)[ 3];
    ((int * ALIGNED(64))a)[ 3] = ((int * ALIGNED(64))b)[ 3];
    ((int * ALIGNED(64))b)[ 3] = t;

    t = ((int * ALIGNED(64))a)[ 4];
    ((int * ALIGNED(64))a)[ 4] = ((int * ALIGNED(64))b)[ 4];
    ((int * ALIGNED(64))b)[ 4] = t;

    t = ((int * ALIGNED(64))a)[ 5];
    ((int * ALIGNED(64))a)[ 5] = ((int * ALIGNED(64))b)[ 5];
    ((int * ALIGNED(64))b)[ 5] = t;

    t = ((int * ALIGNED(64))a)[ 6];
    ((int * ALIGNED(64))a)[ 6] = ((int * ALIGNED(64))b)[ 6];
    ((int * ALIGNED(64))b)[ 6] = t;

    t = ((int * ALIGNED(64))a)[ 7];
    ((int * ALIGNED(64))a)[ 7] = ((int * ALIGNED(64))b)[ 7];
    ((int * ALIGNED(64))b)[ 7] = t;

    t = ((int * ALIGNED(64))a)[ 8];
    ((int * ALIGNED(64))a)[ 8] = ((int * ALIGNED(64))b)[ 8];
    ((int * ALIGNED(64))b)[ 8] = t;

    t = ((int * ALIGNED(64))a)[ 9];
    ((int * ALIGNED(64))a)[ 9] = ((int * ALIGNED(64))b)[ 9];
    ((int * ALIGNED(64))b)[ 9] = t;

    t = ((int * ALIGNED(64))a)[10];
    ((int * ALIGNED(64))a)[10] = ((int * ALIGNED(64))b)[10];
    ((int * ALIGNED(64))b)[10] = t;

    t = ((int * ALIGNED(64))a)[11];
    ((int * ALIGNED(64))a)[11] = ((int * ALIGNED(64))b)[11];
    ((int * ALIGNED(64))b)[11] = t;

    t = ((int * ALIGNED(64))a)[12];
    ((int * ALIGNED(64))a)[12] = ((int * ALIGNED(64))b)[12];
    ((int * ALIGNED(64))b)[12] = t;

    t = ((int * ALIGNED(64))a)[13];
    ((int * ALIGNED(64))a)[13] = ((int * ALIGNED(64))b)[13];
    ((int * ALIGNED(64))b)[13] = t;

    t = ((int * ALIGNED(64))a)[14];
    ((int * ALIGNED(64))a)[14] = ((int * ALIGNED(64))b)[14];
    ((int * ALIGNED(64))b)[14] = t;

    t = ((int * ALIGNED(64))a)[15];
    ((int * ALIGNED(64))a)[15] = ((int * ALIGNED(64))b)[15];
    ((int * ALIGNED(64))b)[15] = t;
  }

  // v16 transposed memory manipulation functions

  inline void load_16x1_tr( const void *a00, const void *a01,
                            const void *a02, const void *a03,
                            const void *a04, const void *a05,
                            const void *a06, const void *a07,
			    const void *a08, const void *a09,
                            const void *a10, const void *a11,
                            const void *a12, const void *a13,
                            const void *a14, const void *a15,
			    v16 &a )
  {
    a.i[ 0] = ((const int *)a00)[0];
    a.i[ 1] = ((const int *)a01)[0];
    a.i[ 2] = ((const int *)a02)[0];
    a.i[ 3] = ((const int *)a03)[0];
    a.i[ 4] = ((const int *)a04)[0];
    a.i[ 5] = ((const int *)a05)[0];
    a.i[ 6] = ((const int *)a06)[0];
    a.i[ 7] = ((const int *)a07)[0];
    a.i[ 8] = ((const int *)a08)[0];
    a.i[ 9] = ((const int *)a09)[0];
    a.i[10] = ((const int *)a10)[0];
    a.i[11] = ((const int *)a11)[0];
    a.i[12] = ((const int *)a12)[0];
    a.i[13] = ((const int *)a13)[0];
    a.i[14] = ((const int *)a14)[0];
    a.i[15] = ((const int *)a15)[0];
  }

  inline void load_16x2_tr( const void * ALIGNED(8) a00,
			    const void * ALIGNED(8) a01,
			    const void * ALIGNED(8) a02,
			    const void * ALIGNED(8) a03,
			    const void * ALIGNED(8) a04,
			    const void * ALIGNED(8) a05,
			    const void * ALIGNED(8) a06,
			    const void * ALIGNED(8) a07,
			    const void * ALIGNED(8) a08,
			    const void * ALIGNED(8) a09,
			    const void * ALIGNED(8) a10,
			    const void * ALIGNED(8) a11,
			    const void * ALIGNED(8) a12,
			    const void * ALIGNED(8) a13,
			    const void * ALIGNED(8) a14,
			    const void * ALIGNED(8) a15,
			    v16 &a, v16 &b )
  {
    a.i[ 0] = ((const int * ALIGNED(8))a00)[0];
    b.i[ 0] = ((const int * ALIGNED(8))a00)[1];

    a.i[ 1] = ((const int * ALIGNED(8))a01)[0];
    b.i[ 1] = ((const int * ALIGNED(8))a01)[1];

    a.i[ 2] = ((const int * ALIGNED(8))a02)[0];
    b.i[ 2] = ((const int * ALIGNED(8))a02)[1];

    a.i[ 3] = ((const int * ALIGNED(8))a03)[0];
    b.i[ 3] = ((const int * ALIGNED(8))a03)[1];

    a.i[ 4] = ((const int * ALIGNED(8))a04)[0];
    b.i[ 4] = ((const int * ALIGNED(8))a04)[1];

    a.i[ 5] = ((const int * ALIGNED(8))a05)[0];
    b.i[ 5] = ((const int * ALIGNED(8))a05)[1];

    a.i[ 6] = ((const int * ALIGNED(8))a06)[0];
    b.i[ 6] = ((const int * ALIGNED(8))a06)[1];

    a.i[ 7] = ((const int * ALIGNED(8))a07)[0];
    b.i[ 7] = ((const int * ALIGNED(8))a07)[1];

    a.i[ 8] = ((const int * ALIGNED(8))a08)[0];
    b.i[ 8] = ((const int * ALIGNED(8))a08)[1];

    a.i[ 9] = ((const int * ALIGNED(8))a09)[0];
    b.i[ 9] = ((const int * ALIGNED(8))a09)[1];

    a.i[10] = ((const int * ALIGNED(8))a10)[0];
    b.i[10] = ((const int * ALIGNED(8))a10)[1];

    a.i[11] = ((const int * ALIGNED(8))a11)[0];
    b.i[11] = ((const int * ALIGNED(8))a11)[1];

    a.i[12] = ((const int * ALIGNED(8))a12)[0];
    b.i[12] = ((const int * ALIGNED(8))a12)[1];

    a.i[13] = ((const int * ALIGNED(8))a13)[0];
    b.i[13] = ((const int * ALIGNED(8))a13)[1];

    a.i[14] = ((const int * ALIGNED(8))a14)[0];
    b.i[14] = ((const int * ALIGNED(8))a14)[1];

    a.i[15] = ((const int * ALIGNED(8))a15)[0];
    b.i[15] = ((const int * ALIGNED(8))a15)[1];
  }

  inline void load_16x3_tr( const void * ALIGNED(64) a00,
                            const void * ALIGNED(64) a01,
                            const void * ALIGNED(64) a02,
                            const void * ALIGNED(64) a03,
			    const void * ALIGNED(64) a04,
			    const void * ALIGNED(64) a05,
			    const void * ALIGNED(64) a06,
			    const void * ALIGNED(64) a07,
			    const void * ALIGNED(64) a08,
                            const void * ALIGNED(64) a09,
                            const void * ALIGNED(64) a10,
                            const void * ALIGNED(64) a11,
			    const void * ALIGNED(64) a12,
			    const void * ALIGNED(64) a13,
			    const void * ALIGNED(64) a14,
			    const void * ALIGNED(64) a15,
			    v16 &a, v16 &b, v16 &c )
  {
    a.i[ 0] = ((const int * ALIGNED(64))a00)[0];
    b.i[ 0] = ((const int * ALIGNED(64))a00)[1];
    c.i[ 0] = ((const int * ALIGNED(64))a00)[2];

    a.i[ 1] = ((const int * ALIGNED(64))a01)[0];
    b.i[ 1] = ((const int * ALIGNED(64))a01)[1];
    c.i[ 1] = ((const int * ALIGNED(64))a01)[2];

    a.i[ 2] = ((const int * ALIGNED(64))a02)[0];
    b.i[ 2] = ((const int * ALIGNED(64))a02)[1];
    c.i[ 2] = ((const int * ALIGNED(64))a02)[2];

    a.i[ 3] = ((const int * ALIGNED(64))a03)[0];
    b.i[ 3] = ((const int * ALIGNED(64))a03)[1];
    c.i[ 3] = ((const int * ALIGNED(64))a03)[2]; 

    a.i[ 4] = ((const int * ALIGNED(64))a04)[0];
    b.i[ 4] = ((const int * ALIGNED(64))a04)[1];
    c.i[ 4] = ((const int * ALIGNED(64))a04)[2];

    a.i[ 5] = ((const int * ALIGNED(64))a05)[0];
    b.i[ 5] = ((const int * ALIGNED(64))a05)[1];
    c.i[ 5] = ((const int * ALIGNED(64))a05)[2];

    a.i[ 6] = ((const int * ALIGNED(64))a06)[0];
    b.i[ 6] = ((const int * ALIGNED(64))a06)[1];
    c.i[ 6] = ((const int * ALIGNED(64))a06)[2];

    a.i[ 7] = ((const int * ALIGNED(64))a07)[0];
    b.i[ 7] = ((const int * ALIGNED(64))a07)[1];
    c.i[ 7] = ((const int * ALIGNED(64))a07)[2]; 

    a.i[ 8] = ((const int * ALIGNED(64))a08)[0];
    b.i[ 8] = ((const int * ALIGNED(64))a08)[1];
    c.i[ 8] = ((const int * ALIGNED(64))a08)[2];

    a.i[ 9] = ((const int * ALIGNED(64))a09)[0];
    b.i[ 9] = ((const int * ALIGNED(64))a09)[1];
    c.i[ 9] = ((const int * ALIGNED(64))a09)[2];

    a.i[10] = ((const int * ALIGNED(64))a10)[0];
    b.i[10] = ((const int * ALIGNED(64))a10)[1];
    c.i[10] = ((const int * ALIGNED(64))a10)[2];

    a.i[11] = ((const int * ALIGNED(64))a11)[0];
    b.i[11] = ((const int * ALIGNED(64))a11)[1];
    c.i[11] = ((const int * ALIGNED(64))a11)[2]; 

    a.i[12] = ((const int * ALIGNED(64))a12)[0];
    b.i[12] = ((const int * ALIGNED(64))a12)[1];
    c.i[12] = ((const int * ALIGNED(64))a12)[2];

    a.i[13] = ((const int * ALIGNED(64))a13)[0];
    b.i[13] = ((const int * ALIGNED(64))a13)[1];
    c.i[13] = ((const int * ALIGNED(64))a13)[2];

    a.i[14] = ((const int * ALIGNED(64))a14)[0];
    b.i[14] = ((const int * ALIGNED(64))a14)[1];
    c.i[14] = ((const int * ALIGNED(64))a14)[2];

    a.i[15] = ((const int * ALIGNED(64))a15)[0];
    b.i[15] = ((const int * ALIGNED(64))a15)[1];
    c.i[15] = ((const int * ALIGNED(64))a15)[2]; 
  }

  inline void load_16x4_tr( const void * ALIGNED(64) a00,
			    const void * ALIGNED(64) a01,
			    const void * ALIGNED(64) a02,
			    const void * ALIGNED(64) a03,
			    const void * ALIGNED(64) a04,
			    const void * ALIGNED(64) a05,
			    const void * ALIGNED(64) a06,
			    const void * ALIGNED(64) a07,
			    const void * ALIGNED(64) a08,
			    const void * ALIGNED(64) a09,
			    const void * ALIGNED(64) a10,
			    const void * ALIGNED(64) a11,
			    const void * ALIGNED(64) a12,
			    const void * ALIGNED(64) a13,
			    const void * ALIGNED(64) a14,
			    const void * ALIGNED(64) a15,
			    v16 &a, v16 &b, v16 &c, v16 &d )
  {
    a.i[ 0] = ((const int * ALIGNED(64))a00)[0];
    b.i[ 0] = ((const int * ALIGNED(64))a00)[1];
    c.i[ 0] = ((const int * ALIGNED(64))a00)[2];
    d.i[ 0] = ((const int * ALIGNED(64))a00)[3];

    a.i[ 1] = ((const int * ALIGNED(64))a01)[0];
    b.i[ 1] = ((const int * ALIGNED(64))a01)[1];
    c.i[ 1] = ((const int * ALIGNED(64))a01)[2];
    d.i[ 1] = ((const int * ALIGNED(64))a01)[3];

    a.i[ 2] = ((const int * ALIGNED(64))a02)[0];
    b.i[ 2] = ((const int * ALIGNED(64))a02)[1];
    c.i[ 2] = ((const int * ALIGNED(64))a02)[2];
    d.i[ 2] = ((const int * ALIGNED(64))a02)[3];

    a.i[ 3] = ((const int * ALIGNED(64))a03)[0];
    b.i[ 3] = ((const int * ALIGNED(64))a03)[1];
    c.i[ 3] = ((const int * ALIGNED(64))a03)[2];
    d.i[ 3] = ((const int * ALIGNED(64))a03)[3];

    a.i[ 4] = ((const int * ALIGNED(64))a04)[0];
    b.i[ 4] = ((const int * ALIGNED(64))a04)[1];
    c.i[ 4] = ((const int * ALIGNED(64))a04)[2];
    d.i[ 4] = ((const int * ALIGNED(64))a04)[3];

    a.i[ 5] = ((const int * ALIGNED(64))a05)[0];
    b.i[ 5] = ((const int * ALIGNED(64))a05)[1];
    c.i[ 5] = ((const int * ALIGNED(64))a05)[2];
    d.i[ 5] = ((const int * ALIGNED(64))a05)[3];

    a.i[ 6] = ((const int * ALIGNED(64))a06)[0];
    b.i[ 6] = ((const int * ALIGNED(64))a06)[1];
    c.i[ 6] = ((const int * ALIGNED(64))a06)[2];
    d.i[ 6] = ((const int * ALIGNED(64))a06)[3];

    a.i[ 7] = ((const int * ALIGNED(64))a07)[0];
    b.i[ 7] = ((const int * ALIGNED(64))a07)[1];
    c.i[ 7] = ((const int * ALIGNED(64))a07)[2];
    d.i[ 7] = ((const int * ALIGNED(64))a07)[3];

    a.i[ 8] = ((const int * ALIGNED(64))a08)[0];
    b.i[ 8] = ((const int * ALIGNED(64))a08)[1];
    c.i[ 8] = ((const int * ALIGNED(64))a08)[2];
    d.i[ 8] = ((const int * ALIGNED(64))a08)[3];

    a.i[ 9] = ((const int * ALIGNED(64))a09)[0];
    b.i[ 9] = ((const int * ALIGNED(64))a09)[1];
    c.i[ 9] = ((const int * ALIGNED(64))a09)[2];
    d.i[ 9] = ((const int * ALIGNED(64))a09)[3];

    a.i[10] = ((const int * ALIGNED(64))a10)[0];
    b.i[10] = ((const int * ALIGNED(64))a10)[1];
    c.i[10] = ((const int * ALIGNED(64))a10)[2];
    d.i[10] = ((const int * ALIGNED(64))a10)[3];

    a.i[11] = ((const int * ALIGNED(64))a11)[0];
    b.i[11] = ((const int * ALIGNED(64))a11)[1];
    c.i[11] = ((const int * ALIGNED(64))a11)[2];
    d.i[11] = ((const int * ALIGNED(64))a11)[3];

    a.i[12] = ((const int * ALIGNED(64))a12)[0];
    b.i[12] = ((const int * ALIGNED(64))a12)[1];
    c.i[12] = ((const int * ALIGNED(64))a12)[2];
    d.i[12] = ((const int * ALIGNED(64))a12)[3];

    a.i[13] = ((const int * ALIGNED(64))a13)[0];
    b.i[13] = ((const int * ALIGNED(64))a13)[1];
    c.i[13] = ((const int * ALIGNED(64))a13)[2];
    d.i[13] = ((const int * ALIGNED(64))a13)[3];

    a.i[14] = ((const int * ALIGNED(64))a14)[0];
    b.i[14] = ((const int * ALIGNED(64))a14)[1];
    c.i[14] = ((const int * ALIGNED(64))a14)[2];
    d.i[14] = ((const int * ALIGNED(64))a14)[3];

    a.i[15] = ((const int * ALIGNED(64))a15)[0];
    b.i[15] = ((const int * ALIGNED(64))a15)[1];
    c.i[15] = ((const int * ALIGNED(64))a15)[2];
    d.i[15] = ((const int * ALIGNED(64))a15)[3];
  }

  inline void load_16x8_tr( const void * ALIGNED(64) a00,
			    const void * ALIGNED(64) a01,
			    const void * ALIGNED(64) a02,
			    const void * ALIGNED(64) a03,
			    const void * ALIGNED(64) a04,
			    const void * ALIGNED(64) a05,
			    const void * ALIGNED(64) a06,
			    const void * ALIGNED(64) a07,
			    const void * ALIGNED(64) a08,
			    const void * ALIGNED(64) a09,
			    const void * ALIGNED(64) a10,
			    const void * ALIGNED(64) a11,
			    const void * ALIGNED(64) a12,
			    const void * ALIGNED(64) a13,
			    const void * ALIGNED(64) a14,
			    const void * ALIGNED(64) a15,
			    v16 &a, v16 &b, v16 &c, v16 &d,
			    v16 &e, v16 &f, v16 &g, v16 &h )
  {
    a.i[ 0] = ((const int * ALIGNED(64))a00)[0];
    b.i[ 0] = ((const int * ALIGNED(64))a00)[1];
    c.i[ 0] = ((const int * ALIGNED(64))a00)[2];
    d.i[ 0] = ((const int * ALIGNED(64))a00)[3];
    e.i[ 0] = ((const int * ALIGNED(64))a00)[4];
    f.i[ 0] = ((const int * ALIGNED(64))a00)[5];
    g.i[ 0] = ((const int * ALIGNED(64))a00)[6];
    h.i[ 0] = ((const int * ALIGNED(64))a00)[7];

    a.i[ 1] = ((const int * ALIGNED(64))a01)[0];
    b.i[ 1] = ((const int * ALIGNED(64))a01)[1];
    c.i[ 1] = ((const int * ALIGNED(64))a01)[2];
    d.i[ 1] = ((const int * ALIGNED(64))a01)[3];
    e.i[ 1] = ((const int * ALIGNED(64))a01)[4];
    f.i[ 1] = ((const int * ALIGNED(64))a01)[5];
    g.i[ 1] = ((const int * ALIGNED(64))a01)[6];
    h.i[ 1] = ((const int * ALIGNED(64))a01)[7];

    a.i[ 2] = ((const int * ALIGNED(64))a02)[0];
    b.i[ 2] = ((const int * ALIGNED(64))a02)[1];
    c.i[ 2] = ((const int * ALIGNED(64))a02)[2];
    d.i[ 2] = ((const int * ALIGNED(64))a02)[3];
    e.i[ 2] = ((const int * ALIGNED(64))a02)[4];
    f.i[ 2] = ((const int * ALIGNED(64))a02)[5];
    g.i[ 2] = ((const int * ALIGNED(64))a02)[6];
    h.i[ 2] = ((const int * ALIGNED(64))a02)[7];

    a.i[ 3] = ((const int * ALIGNED(64))a03)[0];
    b.i[ 3] = ((const int * ALIGNED(64))a03)[1];
    c.i[ 3] = ((const int * ALIGNED(64))a03)[2];
    d.i[ 3] = ((const int * ALIGNED(64))a03)[3];
    e.i[ 3] = ((const int * ALIGNED(64))a03)[4];
    f.i[ 3] = ((const int * ALIGNED(64))a03)[5];
    g.i[ 3] = ((const int * ALIGNED(64))a03)[6];
    h.i[ 3] = ((const int * ALIGNED(64))a03)[7];

    a.i[ 4] = ((const int * ALIGNED(64))a04)[0];
    b.i[ 4] = ((const int * ALIGNED(64))a04)[1];
    c.i[ 4] = ((const int * ALIGNED(64))a04)[2];
    d.i[ 4] = ((const int * ALIGNED(64))a04)[3];
    e.i[ 4] = ((const int * ALIGNED(64))a04)[4];
    f.i[ 4] = ((const int * ALIGNED(64))a04)[5];
    g.i[ 4] = ((const int * ALIGNED(64))a04)[6];
    h.i[ 4] = ((const int * ALIGNED(64))a04)[7];

    a.i[ 5] = ((const int * ALIGNED(64))a05)[0];
    b.i[ 5] = ((const int * ALIGNED(64))a05)[1];
    c.i[ 5] = ((const int * ALIGNED(64))a05)[2];
    d.i[ 5] = ((const int * ALIGNED(64))a05)[3];
    e.i[ 5] = ((const int * ALIGNED(64))a05)[4];
    f.i[ 5] = ((const int * ALIGNED(64))a05)[5];
    g.i[ 5] = ((const int * ALIGNED(64))a05)[6];
    h.i[ 5] = ((const int * ALIGNED(64))a05)[7];

    a.i[ 6] = ((const int * ALIGNED(64))a06)[0];
    b.i[ 6] = ((const int * ALIGNED(64))a06)[1];
    c.i[ 6] = ((const int * ALIGNED(64))a06)[2];
    d.i[ 6] = ((const int * ALIGNED(64))a06)[3];
    e.i[ 6] = ((const int * ALIGNED(64))a06)[4];
    f.i[ 6] = ((const int * ALIGNED(64))a06)[5];
    g.i[ 6] = ((const int * ALIGNED(64))a06)[6];
    h.i[ 6] = ((const int * ALIGNED(64))a06)[7];

    a.i[ 7] = ((const int * ALIGNED(64))a07)[0];
    b.i[ 7] = ((const int * ALIGNED(64))a07)[1];
    c.i[ 7] = ((const int * ALIGNED(64))a07)[2];
    d.i[ 7] = ((const int * ALIGNED(64))a07)[3];
    e.i[ 7] = ((const int * ALIGNED(64))a07)[4];
    f.i[ 7] = ((const int * ALIGNED(64))a07)[5];
    g.i[ 7] = ((const int * ALIGNED(64))a07)[6];
    h.i[ 7] = ((const int * ALIGNED(64))a07)[7];

    a.i[ 8] = ((const int * ALIGNED(64))a08)[0];
    b.i[ 8] = ((const int * ALIGNED(64))a08)[1];
    c.i[ 8] = ((const int * ALIGNED(64))a08)[2];
    d.i[ 8] = ((const int * ALIGNED(64))a08)[3];
    e.i[ 8] = ((const int * ALIGNED(64))a08)[4];
    f.i[ 8] = ((const int * ALIGNED(64))a08)[5];
    g.i[ 8] = ((const int * ALIGNED(64))a08)[6];
    h.i[ 8] = ((const int * ALIGNED(64))a08)[7];

    a.i[ 9] = ((const int * ALIGNED(64))a09)[0];
    b.i[ 9] = ((const int * ALIGNED(64))a09)[1];
    c.i[ 9] = ((const int * ALIGNED(64))a09)[2];
    d.i[ 9] = ((const int * ALIGNED(64))a09)[3];
    e.i[ 9] = ((const int * ALIGNED(64))a09)[4];
    f.i[ 9] = ((const int * ALIGNED(64))a09)[5];
    g.i[ 9] = ((const int * ALIGNED(64))a09)[6];
    h.i[ 9] = ((const int * ALIGNED(64))a09)[7];

    a.i[10] = ((const int * ALIGNED(64))a10)[0];
    b.i[10] = ((const int * ALIGNED(64))a10)[1];
    c.i[10] = ((const int * ALIGNED(64))a10)[2];
    d.i[10] = ((const int * ALIGNED(64))a10)[3];
    e.i[10] = ((const int * ALIGNED(64))a10)[4];
    f.i[10] = ((const int * ALIGNED(64))a10)[5];
    g.i[10] = ((const int * ALIGNED(64))a10)[6];
    h.i[10] = ((const int * ALIGNED(64))a10)[7];

    a.i[11] = ((const int * ALIGNED(64))a11)[0];
    b.i[11] = ((const int * ALIGNED(64))a11)[1];
    c.i[11] = ((const int * ALIGNED(64))a11)[2];
    d.i[11] = ((const int * ALIGNED(64))a11)[3];
    e.i[11] = ((const int * ALIGNED(64))a11)[4];
    f.i[11] = ((const int * ALIGNED(64))a11)[5];
    g.i[11] = ((const int * ALIGNED(64))a11)[6];
    h.i[11] = ((const int * ALIGNED(64))a11)[7];

    a.i[12] = ((const int * ALIGNED(64))a12)[0];
    b.i[12] = ((const int * ALIGNED(64))a12)[1];
    c.i[12] = ((const int * ALIGNED(64))a12)[2];
    d.i[12] = ((const int * ALIGNED(64))a12)[3];
    e.i[12] = ((const int * ALIGNED(64))a12)[4];
    f.i[12] = ((const int * ALIGNED(64))a12)[5];
    g.i[12] = ((const int * ALIGNED(64))a12)[6];
    h.i[12] = ((const int * ALIGNED(64))a12)[7];

    a.i[13] = ((const int * ALIGNED(64))a13)[0];
    b.i[13] = ((const int * ALIGNED(64))a13)[1];
    c.i[13] = ((const int * ALIGNED(64))a13)[2];
    d.i[13] = ((const int * ALIGNED(64))a13)[3];
    e.i[13] = ((const int * ALIGNED(64))a13)[4];
    f.i[13] = ((const int * ALIGNED(64))a13)[5];
    g.i[13] = ((const int * ALIGNED(64))a13)[6];
    h.i[13] = ((const int * ALIGNED(64))a13)[7];

    a.i[14] = ((const int * ALIGNED(64))a14)[0];
    b.i[14] = ((const int * ALIGNED(64))a14)[1];
    c.i[14] = ((const int * ALIGNED(64))a14)[2];
    d.i[14] = ((const int * ALIGNED(64))a14)[3];
    e.i[14] = ((const int * ALIGNED(64))a14)[4];
    f.i[14] = ((const int * ALIGNED(64))a14)[5];
    g.i[14] = ((const int * ALIGNED(64))a14)[6];
    h.i[14] = ((const int * ALIGNED(64))a14)[7];

    a.i[15] = ((const int * ALIGNED(64))a15)[0];
    b.i[15] = ((const int * ALIGNED(64))a15)[1];
    c.i[15] = ((const int * ALIGNED(64))a15)[2];
    d.i[15] = ((const int * ALIGNED(64))a15)[3];
    e.i[15] = ((const int * ALIGNED(64))a15)[4];
    f.i[15] = ((const int * ALIGNED(64))a15)[5];
    g.i[15] = ((const int * ALIGNED(64))a15)[6];
    h.i[15] = ((const int * ALIGNED(64))a15)[7];
  }

  inline void load_16x16_tr( const void * ALIGNED(64) a00,
			     const void * ALIGNED(64) a01,
			     const void * ALIGNED(64) a02,
			     const void * ALIGNED(64) a03,
			     const void * ALIGNED(64) a04,
			     const void * ALIGNED(64) a05,
			     const void * ALIGNED(64) a06,
			     const void * ALIGNED(64) a07,
			     const void * ALIGNED(64) a08,
			     const void * ALIGNED(64) a09,
			     const void * ALIGNED(64) a10,
			     const void * ALIGNED(64) a11,
			     const void * ALIGNED(64) a12,
			     const void * ALIGNED(64) a13,
			     const void * ALIGNED(64) a14,
			     const void * ALIGNED(64) a15,
			     v16 &b00, v16 &b01, v16 &b02, v16 &b03,
			     v16 &b04, v16 &b05, v16 &b06, v16 &b07,
			     v16 &b08, v16 &b09, v16 &b10, v16 &b11,
			     v16 &b12, v16 &b13, v16 &b14, v16 &b15 )
  {
    b00.i[ 0] = ((const int * ALIGNED(64))a00)[ 0];
    b01.i[ 0] = ((const int * ALIGNED(64))a00)[ 1];
    b02.i[ 0] = ((const int * ALIGNED(64))a00)[ 2];
    b03.i[ 0] = ((const int * ALIGNED(64))a00)[ 3];
    b04.i[ 0] = ((const int * ALIGNED(64))a00)[ 4];
    b05.i[ 0] = ((const int * ALIGNED(64))a00)[ 5];
    b06.i[ 0] = ((const int * ALIGNED(64))a00)[ 6];
    b07.i[ 0] = ((const int * ALIGNED(64))a00)[ 7];
    b08.i[ 0] = ((const int * ALIGNED(64))a00)[ 8];
    b09.i[ 0] = ((const int * ALIGNED(64))a00)[ 9];
    b10.i[ 0] = ((const int * ALIGNED(64))a00)[10];
    b11.i[ 0] = ((const int * ALIGNED(64))a00)[11];
    b12.i[ 0] = ((const int * ALIGNED(64))a00)[12];
    b13.i[ 0] = ((const int * ALIGNED(64))a00)[13];
    b14.i[ 0] = ((const int * ALIGNED(64))a00)[14];
    b15.i[ 0] = ((const int * ALIGNED(64))a00)[15];

    b00.i[ 1] = ((const int * ALIGNED(64))a01)[ 0];
    b01.i[ 1] = ((const int * ALIGNED(64))a01)[ 1];
    b02.i[ 1] = ((const int * ALIGNED(64))a01)[ 2];
    b03.i[ 1] = ((const int * ALIGNED(64))a01)[ 3];
    b04.i[ 1] = ((const int * ALIGNED(64))a01)[ 4];
    b05.i[ 1] = ((const int * ALIGNED(64))a01)[ 5];
    b06.i[ 1] = ((const int * ALIGNED(64))a01)[ 6];
    b07.i[ 1] = ((const int * ALIGNED(64))a01)[ 7];
    b08.i[ 1] = ((const int * ALIGNED(64))a01)[ 8];
    b09.i[ 1] = ((const int * ALIGNED(64))a01)[ 9];
    b10.i[ 1] = ((const int * ALIGNED(64))a01)[10];
    b11.i[ 1] = ((const int * ALIGNED(64))a01)[11];
    b12.i[ 1] = ((const int * ALIGNED(64))a01)[12];
    b13.i[ 1] = ((const int * ALIGNED(64))a01)[13];
    b14.i[ 1] = ((const int * ALIGNED(64))a01)[14];
    b15.i[ 1] = ((const int * ALIGNED(64))a01)[15];

    b00.i[ 2] = ((const int * ALIGNED(64))a02)[ 0];
    b01.i[ 2] = ((const int * ALIGNED(64))a02)[ 1];
    b02.i[ 2] = ((const int * ALIGNED(64))a02)[ 2];
    b03.i[ 2] = ((const int * ALIGNED(64))a02)[ 3];
    b04.i[ 2] = ((const int * ALIGNED(64))a02)[ 4];
    b05.i[ 2] = ((const int * ALIGNED(64))a02)[ 5];
    b06.i[ 2] = ((const int * ALIGNED(64))a02)[ 6];
    b07.i[ 2] = ((const int * ALIGNED(64))a02)[ 7];
    b08.i[ 2] = ((const int * ALIGNED(64))a02)[ 8];
    b09.i[ 2] = ((const int * ALIGNED(64))a02)[ 9];
    b10.i[ 2] = ((const int * ALIGNED(64))a02)[10];
    b11.i[ 2] = ((const int * ALIGNED(64))a02)[11];
    b12.i[ 2] = ((const int * ALIGNED(64))a02)[12];
    b13.i[ 2] = ((const int * ALIGNED(64))a02)[13];
    b14.i[ 2] = ((const int * ALIGNED(64))a02)[14];
    b15.i[ 2] = ((const int * ALIGNED(64))a02)[15];

    b00.i[ 3] = ((const int * ALIGNED(64))a03)[ 0];
    b01.i[ 3] = ((const int * ALIGNED(64))a03)[ 1];
    b02.i[ 3] = ((const int * ALIGNED(64))a03)[ 2];
    b03.i[ 3] = ((const int * ALIGNED(64))a03)[ 3];
    b04.i[ 3] = ((const int * ALIGNED(64))a03)[ 4];
    b05.i[ 3] = ((const int * ALIGNED(64))a03)[ 5];
    b06.i[ 3] = ((const int * ALIGNED(64))a03)[ 6];
    b07.i[ 3] = ((const int * ALIGNED(64))a03)[ 7];
    b08.i[ 3] = ((const int * ALIGNED(64))a03)[ 8];
    b09.i[ 3] = ((const int * ALIGNED(64))a03)[ 9];
    b10.i[ 3] = ((const int * ALIGNED(64))a03)[10];
    b11.i[ 3] = ((const int * ALIGNED(64))a03)[11];
    b12.i[ 3] = ((const int * ALIGNED(64))a03)[12];
    b13.i[ 3] = ((const int * ALIGNED(64))a03)[13];
    b14.i[ 3] = ((const int * ALIGNED(64))a03)[14];
    b15.i[ 3] = ((const int * ALIGNED(64))a03)[15];

    b00.i[ 4] = ((const int * ALIGNED(64))a04)[ 0];
    b01.i[ 4] = ((const int * ALIGNED(64))a04)[ 1];
    b02.i[ 4] = ((const int * ALIGNED(64))a04)[ 2];
    b03.i[ 4] = ((const int * ALIGNED(64))a04)[ 3];
    b04.i[ 4] = ((const int * ALIGNED(64))a04)[ 4];
    b05.i[ 4] = ((const int * ALIGNED(64))a04)[ 5];
    b06.i[ 4] = ((const int * ALIGNED(64))a04)[ 6];
    b07.i[ 4] = ((const int * ALIGNED(64))a04)[ 7];
    b08.i[ 4] = ((const int * ALIGNED(64))a04)[ 8];
    b09.i[ 4] = ((const int * ALIGNED(64))a04)[ 9];
    b10.i[ 4] = ((const int * ALIGNED(64))a04)[10];
    b11.i[ 4] = ((const int * ALIGNED(64))a04)[11];
    b12.i[ 4] = ((const int * ALIGNED(64))a04)[12];
    b13.i[ 4] = ((const int * ALIGNED(64))a04)[13];
    b14.i[ 4] = ((const int * ALIGNED(64))a04)[14];
    b15.i[ 4] = ((const int * ALIGNED(64))a04)[15];

    b00.i[ 5] = ((const int * ALIGNED(64))a05)[ 0];
    b01.i[ 5] = ((const int * ALIGNED(64))a05)[ 1];
    b02.i[ 5] = ((const int * ALIGNED(64))a05)[ 2];
    b03.i[ 5] = ((const int * ALIGNED(64))a05)[ 3];
    b04.i[ 5] = ((const int * ALIGNED(64))a05)[ 4];
    b05.i[ 5] = ((const int * ALIGNED(64))a05)[ 5];
    b06.i[ 5] = ((const int * ALIGNED(64))a05)[ 6];
    b07.i[ 5] = ((const int * ALIGNED(64))a05)[ 7];
    b08.i[ 5] = ((const int * ALIGNED(64))a05)[ 8];
    b09.i[ 5] = ((const int * ALIGNED(64))a05)[ 9];
    b10.i[ 5] = ((const int * ALIGNED(64))a05)[10];
    b11.i[ 5] = ((const int * ALIGNED(64))a05)[11];
    b12.i[ 5] = ((const int * ALIGNED(64))a05)[12];
    b13.i[ 5] = ((const int * ALIGNED(64))a05)[13];
    b14.i[ 5] = ((const int * ALIGNED(64))a05)[14];
    b15.i[ 5] = ((const int * ALIGNED(64))a05)[15];

    b00.i[ 6] = ((const int * ALIGNED(64))a06)[ 0];
    b01.i[ 6] = ((const int * ALIGNED(64))a06)[ 1];
    b02.i[ 6] = ((const int * ALIGNED(64))a06)[ 2];
    b03.i[ 6] = ((const int * ALIGNED(64))a06)[ 3];
    b04.i[ 6] = ((const int * ALIGNED(64))a06)[ 4];
    b05.i[ 6] = ((const int * ALIGNED(64))a06)[ 5];
    b06.i[ 6] = ((const int * ALIGNED(64))a06)[ 6];
    b07.i[ 6] = ((const int * ALIGNED(64))a06)[ 7];
    b08.i[ 6] = ((const int * ALIGNED(64))a06)[ 8];
    b09.i[ 6] = ((const int * ALIGNED(64))a06)[ 9];
    b10.i[ 6] = ((const int * ALIGNED(64))a06)[10];
    b11.i[ 6] = ((const int * ALIGNED(64))a06)[11];
    b12.i[ 6] = ((const int * ALIGNED(64))a06)[12];
    b13.i[ 6] = ((const int * ALIGNED(64))a06)[13];
    b14.i[ 6] = ((const int * ALIGNED(64))a06)[14];
    b15.i[ 6] = ((const int * ALIGNED(64))a06)[15];

    b00.i[ 7] = ((const int * ALIGNED(64))a07)[ 0];
    b01.i[ 7] = ((const int * ALIGNED(64))a07)[ 1];
    b02.i[ 7] = ((const int * ALIGNED(64))a07)[ 2];
    b03.i[ 7] = ((const int * ALIGNED(64))a07)[ 3];
    b04.i[ 7] = ((const int * ALIGNED(64))a07)[ 4];
    b05.i[ 7] = ((const int * ALIGNED(64))a07)[ 5];
    b06.i[ 7] = ((const int * ALIGNED(64))a07)[ 6];
    b07.i[ 7] = ((const int * ALIGNED(64))a07)[ 7];
    b08.i[ 7] = ((const int * ALIGNED(64))a07)[ 8];
    b09.i[ 7] = ((const int * ALIGNED(64))a07)[ 9];
    b10.i[ 7] = ((const int * ALIGNED(64))a07)[10];
    b11.i[ 7] = ((const int * ALIGNED(64))a07)[11];
    b12.i[ 7] = ((const int * ALIGNED(64))a07)[12];
    b13.i[ 7] = ((const int * ALIGNED(64))a07)[13];
    b14.i[ 7] = ((const int * ALIGNED(64))a07)[14];
    b15.i[ 7] = ((const int * ALIGNED(64))a07)[15];

    b00.i[ 8] = ((const int * ALIGNED(64))a08)[ 0];
    b01.i[ 8] = ((const int * ALIGNED(64))a08)[ 1];
    b02.i[ 8] = ((const int * ALIGNED(64))a08)[ 2];
    b03.i[ 8] = ((const int * ALIGNED(64))a08)[ 3];
    b04.i[ 8] = ((const int * ALIGNED(64))a08)[ 4];
    b05.i[ 8] = ((const int * ALIGNED(64))a08)[ 5];
    b06.i[ 8] = ((const int * ALIGNED(64))a08)[ 6];
    b07.i[ 8] = ((const int * ALIGNED(64))a08)[ 7];
    b08.i[ 8] = ((const int * ALIGNED(64))a08)[ 8];
    b09.i[ 8] = ((const int * ALIGNED(64))a08)[ 9];
    b10.i[ 8] = ((const int * ALIGNED(64))a08)[10];
    b11.i[ 8] = ((const int * ALIGNED(64))a08)[11];
    b12.i[ 8] = ((const int * ALIGNED(64))a08)[12];
    b13.i[ 8] = ((const int * ALIGNED(64))a08)[13];
    b14.i[ 8] = ((const int * ALIGNED(64))a08)[14];
    b15.i[ 8] = ((const int * ALIGNED(64))a08)[15];

    b00.i[ 9] = ((const int * ALIGNED(64))a09)[ 0];
    b01.i[ 9] = ((const int * ALIGNED(64))a09)[ 1];
    b02.i[ 9] = ((const int * ALIGNED(64))a09)[ 2];
    b03.i[ 9] = ((const int * ALIGNED(64))a09)[ 3];
    b04.i[ 9] = ((const int * ALIGNED(64))a09)[ 4];
    b05.i[ 9] = ((const int * ALIGNED(64))a09)[ 5];
    b06.i[ 9] = ((const int * ALIGNED(64))a09)[ 6];
    b07.i[ 9] = ((const int * ALIGNED(64))a09)[ 7];
    b08.i[ 9] = ((const int * ALIGNED(64))a09)[ 8];
    b09.i[ 9] = ((const int * ALIGNED(64))a09)[ 9];
    b10.i[ 9] = ((const int * ALIGNED(64))a09)[10];
    b11.i[ 9] = ((const int * ALIGNED(64))a09)[11];
    b12.i[ 9] = ((const int * ALIGNED(64))a09)[12];
    b13.i[ 9] = ((const int * ALIGNED(64))a09)[13];
    b14.i[ 9] = ((const int * ALIGNED(64))a09)[14];
    b15.i[ 9] = ((const int * ALIGNED(64))a09)[15];

    b00.i[10] = ((const int * ALIGNED(64))a10)[ 0];
    b01.i[10] = ((const int * ALIGNED(64))a10)[ 1];
    b02.i[10] = ((const int * ALIGNED(64))a10)[ 2];
    b03.i[10] = ((const int * ALIGNED(64))a10)[ 3];
    b04.i[10] = ((const int * ALIGNED(64))a10)[ 4];
    b05.i[10] = ((const int * ALIGNED(64))a10)[ 5];
    b06.i[10] = ((const int * ALIGNED(64))a10)[ 6];
    b07.i[10] = ((const int * ALIGNED(64))a10)[ 7];
    b08.i[10] = ((const int * ALIGNED(64))a10)[ 8];
    b09.i[10] = ((const int * ALIGNED(64))a10)[ 9];
    b10.i[10] = ((const int * ALIGNED(64))a10)[10];
    b11.i[10] = ((const int * ALIGNED(64))a10)[11];
    b12.i[10] = ((const int * ALIGNED(64))a10)[12];
    b13.i[10] = ((const int * ALIGNED(64))a10)[13];
    b14.i[10] = ((const int * ALIGNED(64))a10)[14];
    b15.i[10] = ((const int * ALIGNED(64))a10)[15];

    b00.i[11] = ((const int * ALIGNED(64))a11)[ 0];
    b01.i[11] = ((const int * ALIGNED(64))a11)[ 1];
    b02.i[11] = ((const int * ALIGNED(64))a11)[ 2];
    b03.i[11] = ((const int * ALIGNED(64))a11)[ 3];
    b04.i[11] = ((const int * ALIGNED(64))a11)[ 4];
    b05.i[11] = ((const int * ALIGNED(64))a11)[ 5];
    b06.i[11] = ((const int * ALIGNED(64))a11)[ 6];
    b07.i[11] = ((const int * ALIGNED(64))a11)[ 7];
    b08.i[11] = ((const int * ALIGNED(64))a11)[ 8];
    b09.i[11] = ((const int * ALIGNED(64))a11)[ 9];
    b10.i[11] = ((const int * ALIGNED(64))a11)[10];
    b11.i[11] = ((const int * ALIGNED(64))a11)[11];
    b12.i[11] = ((const int * ALIGNED(64))a11)[12];
    b13.i[11] = ((const int * ALIGNED(64))a11)[13];
    b14.i[11] = ((const int * ALIGNED(64))a11)[14];
    b15.i[11] = ((const int * ALIGNED(64))a11)[15];

    b00.i[12] = ((const int * ALIGNED(64))a12)[ 0];
    b01.i[12] = ((const int * ALIGNED(64))a12)[ 1];
    b02.i[12] = ((const int * ALIGNED(64))a12)[ 2];
    b03.i[12] = ((const int * ALIGNED(64))a12)[ 3];
    b04.i[12] = ((const int * ALIGNED(64))a12)[ 4];
    b05.i[12] = ((const int * ALIGNED(64))a12)[ 5];
    b06.i[12] = ((const int * ALIGNED(64))a12)[ 6];
    b07.i[12] = ((const int * ALIGNED(64))a12)[ 7];
    b08.i[12] = ((const int * ALIGNED(64))a12)[ 8];
    b09.i[12] = ((const int * ALIGNED(64))a12)[ 9];
    b10.i[12] = ((const int * ALIGNED(64))a12)[10];
    b11.i[12] = ((const int * ALIGNED(64))a12)[11];
    b12.i[12] = ((const int * ALIGNED(64))a12)[12];
    b13.i[12] = ((const int * ALIGNED(64))a12)[13];
    b14.i[12] = ((const int * ALIGNED(64))a12)[14];
    b15.i[12] = ((const int * ALIGNED(64))a12)[15];

    b00.i[13] = ((const int * ALIGNED(64))a13)[ 0];
    b01.i[13] = ((const int * ALIGNED(64))a13)[ 1];
    b02.i[13] = ((const int * ALIGNED(64))a13)[ 2];
    b03.i[13] = ((const int * ALIGNED(64))a13)[ 3];
    b04.i[13] = ((const int * ALIGNED(64))a13)[ 4];
    b05.i[13] = ((const int * ALIGNED(64))a13)[ 5];
    b06.i[13] = ((const int * ALIGNED(64))a13)[ 6];
    b07.i[13] = ((const int * ALIGNED(64))a13)[ 7];
    b08.i[13] = ((const int * ALIGNED(64))a13)[ 8];
    b09.i[13] = ((const int * ALIGNED(64))a13)[ 9];
    b10.i[13] = ((const int * ALIGNED(64))a13)[10];
    b11.i[13] = ((const int * ALIGNED(64))a13)[11];
    b12.i[13] = ((const int * ALIGNED(64))a13)[12];
    b13.i[13] = ((const int * ALIGNED(64))a13)[13];
    b14.i[13] = ((const int * ALIGNED(64))a13)[14];
    b15.i[13] = ((const int * ALIGNED(64))a13)[15];

    b00.i[14] = ((const int * ALIGNED(64))a14)[ 0];
    b01.i[14] = ((const int * ALIGNED(64))a14)[ 1];
    b02.i[14] = ((const int * ALIGNED(64))a14)[ 2];
    b03.i[14] = ((const int * ALIGNED(64))a14)[ 3];
    b04.i[14] = ((const int * ALIGNED(64))a14)[ 4];
    b05.i[14] = ((const int * ALIGNED(64))a14)[ 5];
    b06.i[14] = ((const int * ALIGNED(64))a14)[ 6];
    b07.i[14] = ((const int * ALIGNED(64))a14)[ 7];
    b08.i[14] = ((const int * ALIGNED(64))a14)[ 8];
    b09.i[14] = ((const int * ALIGNED(64))a14)[ 9];
    b10.i[14] = ((const int * ALIGNED(64))a14)[10];
    b11.i[14] = ((const int * ALIGNED(64))a14)[11];
    b12.i[14] = ((const int * ALIGNED(64))a14)[12];
    b13.i[14] = ((const int * ALIGNED(64))a14)[13];
    b14.i[14] = ((const int * ALIGNED(64))a14)[14];
    b15.i[14] = ((const int * ALIGNED(64))a14)[15];

    b00.i[15] = ((const int * ALIGNED(64))a15)[ 0];
    b01.i[15] = ((const int * ALIGNED(64))a15)[ 1];
    b02.i[15] = ((const int * ALIGNED(64))a15)[ 2];
    b03.i[15] = ((const int * ALIGNED(64))a15)[ 3];
    b04.i[15] = ((const int * ALIGNED(64))a15)[ 4];
    b05.i[15] = ((const int * ALIGNED(64))a15)[ 5];
    b06.i[15] = ((const int * ALIGNED(64))a15)[ 6];
    b07.i[15] = ((const int * ALIGNED(64))a15)[ 7];
    b08.i[15] = ((const int * ALIGNED(64))a15)[ 8];
    b09.i[15] = ((const int * ALIGNED(64))a15)[ 9];
    b10.i[15] = ((const int * ALIGNED(64))a15)[10];
    b11.i[15] = ((const int * ALIGNED(64))a15)[11];
    b12.i[15] = ((const int * ALIGNED(64))a15)[12];
    b13.i[15] = ((const int * ALIGNED(64))a15)[13];
    b14.i[15] = ((const int * ALIGNED(64))a15)[14];
    b15.i[15] = ((const int * ALIGNED(64))a15)[15];
  }

  inline void load_16x16_tr_p( const void * ALIGNED(64) a00,
			       const void * ALIGNED(64) a01,
			       const void * ALIGNED(64) a02,
			       const void * ALIGNED(64) a03,
			       const void * ALIGNED(64) a04,
			       const void * ALIGNED(64) a05,
			       const void * ALIGNED(64) a06,
			       const void * ALIGNED(64) a07,
			       const void * ALIGNED(64) a08,
			       const void * ALIGNED(64) a09,
			       const void * ALIGNED(64) a10,
			       const void * ALIGNED(64) a11,
			       const void * ALIGNED(64) a12,
			       const void * ALIGNED(64) a13,
			       const void * ALIGNED(64) a14,
			       const void * ALIGNED(64) a15,
			       v16 &b00, v16 &b01, v16 &b02, v16 &b03,
			       v16 &b04, v16 &b05, v16 &b06, v16 &b07,
			       v16 &b08, v16 &b09, v16 &b10, v16 &b11,
			       v16 &b12, v16 &b13, v16 &b14, v16 &b15 )
  {
    b00.i[ 0] = ((const int * ALIGNED(64))a00)[ 0];
    b01.i[ 0] = ((const int * ALIGNED(64))a00)[ 1];
    b02.i[ 0] = ((const int * ALIGNED(64))a00)[ 2];
    b03.i[ 0] = ((const int * ALIGNED(64))a00)[ 3];
    b04.i[ 0] = ((const int * ALIGNED(64))a00)[ 4];
    b05.i[ 0] = ((const int * ALIGNED(64))a00)[ 5];
    b06.i[ 0] = ((const int * ALIGNED(64))a00)[ 6];
    b07.i[ 0] = ((const int * ALIGNED(64))a00)[ 7];
    b00.i[ 1] = ((const int * ALIGNED(64))a00)[ 8];
    b01.i[ 1] = ((const int * ALIGNED(64))a00)[ 9];
    b02.i[ 1] = ((const int * ALIGNED(64))a00)[10];
    b03.i[ 1] = ((const int * ALIGNED(64))a00)[11];
    b04.i[ 1] = ((const int * ALIGNED(64))a00)[12];
    b05.i[ 1] = ((const int * ALIGNED(64))a00)[13];
    b06.i[ 1] = ((const int * ALIGNED(64))a00)[14];
    b07.i[ 1] = ((const int * ALIGNED(64))a00)[15];

    b00.i[ 2] = ((const int * ALIGNED(64))a01)[ 0];
    b01.i[ 2] = ((const int * ALIGNED(64))a01)[ 1];
    b02.i[ 2] = ((const int * ALIGNED(64))a01)[ 2];
    b03.i[ 2] = ((const int * ALIGNED(64))a01)[ 3];
    b04.i[ 2] = ((const int * ALIGNED(64))a01)[ 4];
    b05.i[ 2] = ((const int * ALIGNED(64))a01)[ 5];
    b06.i[ 2] = ((const int * ALIGNED(64))a01)[ 6];
    b07.i[ 2] = ((const int * ALIGNED(64))a01)[ 7];
    b00.i[ 3] = ((const int * ALIGNED(64))a01)[ 8];
    b01.i[ 3] = ((const int * ALIGNED(64))a01)[ 9];
    b02.i[ 3] = ((const int * ALIGNED(64))a01)[10];
    b03.i[ 3] = ((const int * ALIGNED(64))a01)[11];
    b04.i[ 3] = ((const int * ALIGNED(64))a01)[12];
    b05.i[ 3] = ((const int * ALIGNED(64))a01)[13];
    b06.i[ 3] = ((const int * ALIGNED(64))a01)[14];
    b07.i[ 3] = ((const int * ALIGNED(64))a01)[15];

    b00.i[ 4] = ((const int * ALIGNED(64))a02)[ 0];
    b01.i[ 4] = ((const int * ALIGNED(64))a02)[ 1];
    b02.i[ 4] = ((const int * ALIGNED(64))a02)[ 2];
    b03.i[ 4] = ((const int * ALIGNED(64))a02)[ 3];
    b04.i[ 4] = ((const int * ALIGNED(64))a02)[ 4];
    b05.i[ 4] = ((const int * ALIGNED(64))a02)[ 5];
    b06.i[ 4] = ((const int * ALIGNED(64))a02)[ 6];
    b07.i[ 4] = ((const int * ALIGNED(64))a02)[ 7];
    b00.i[ 5] = ((const int * ALIGNED(64))a02)[ 8];
    b01.i[ 5] = ((const int * ALIGNED(64))a02)[ 9];
    b02.i[ 5] = ((const int * ALIGNED(64))a02)[10];
    b03.i[ 5] = ((const int * ALIGNED(64))a02)[11];
    b04.i[ 5] = ((const int * ALIGNED(64))a02)[12];
    b05.i[ 5] = ((const int * ALIGNED(64))a02)[13];
    b06.i[ 5] = ((const int * ALIGNED(64))a02)[14];
    b07.i[ 5] = ((const int * ALIGNED(64))a02)[15];

    b00.i[ 6] = ((const int * ALIGNED(64))a03)[ 0];
    b01.i[ 6] = ((const int * ALIGNED(64))a03)[ 1];
    b02.i[ 6] = ((const int * ALIGNED(64))a03)[ 2];
    b03.i[ 6] = ((const int * ALIGNED(64))a03)[ 3];
    b04.i[ 6] = ((const int * ALIGNED(64))a03)[ 4];
    b05.i[ 6] = ((const int * ALIGNED(64))a03)[ 5];
    b06.i[ 6] = ((const int * ALIGNED(64))a03)[ 6];
    b07.i[ 6] = ((const int * ALIGNED(64))a03)[ 7];
    b00.i[ 7] = ((const int * ALIGNED(64))a03)[ 8];
    b01.i[ 7] = ((const int * ALIGNED(64))a03)[ 9];
    b02.i[ 7] = ((const int * ALIGNED(64))a03)[10];
    b03.i[ 7] = ((const int * ALIGNED(64))a03)[11];
    b04.i[ 7] = ((const int * ALIGNED(64))a03)[12];
    b05.i[ 7] = ((const int * ALIGNED(64))a03)[13];
    b06.i[ 7] = ((const int * ALIGNED(64))a03)[14];
    b07.i[ 7] = ((const int * ALIGNED(64))a03)[15];

    b00.i[ 8] = ((const int * ALIGNED(64))a04)[ 0];
    b01.i[ 8] = ((const int * ALIGNED(64))a04)[ 1];
    b02.i[ 8] = ((const int * ALIGNED(64))a04)[ 2];
    b03.i[ 8] = ((const int * ALIGNED(64))a04)[ 3];
    b04.i[ 8] = ((const int * ALIGNED(64))a04)[ 4];
    b05.i[ 8] = ((const int * ALIGNED(64))a04)[ 5];
    b06.i[ 8] = ((const int * ALIGNED(64))a04)[ 6];
    b07.i[ 8] = ((const int * ALIGNED(64))a04)[ 7];
    b00.i[ 9] = ((const int * ALIGNED(64))a04)[ 8];
    b01.i[ 9] = ((const int * ALIGNED(64))a04)[ 9];
    b02.i[ 9] = ((const int * ALIGNED(64))a04)[10];
    b03.i[ 9] = ((const int * ALIGNED(64))a04)[11];
    b04.i[ 9] = ((const int * ALIGNED(64))a04)[12];
    b05.i[ 9] = ((const int * ALIGNED(64))a04)[13];
    b06.i[ 9] = ((const int * ALIGNED(64))a04)[14];
    b07.i[ 9] = ((const int * ALIGNED(64))a04)[15];

    b00.i[10] = ((const int * ALIGNED(64))a05)[ 0];
    b01.i[10] = ((const int * ALIGNED(64))a05)[ 1];
    b02.i[10] = ((const int * ALIGNED(64))a05)[ 2];
    b03.i[10] = ((const int * ALIGNED(64))a05)[ 3];
    b04.i[10] = ((const int * ALIGNED(64))a05)[ 4];
    b05.i[10] = ((const int * ALIGNED(64))a05)[ 5];
    b06.i[10] = ((const int * ALIGNED(64))a05)[ 6];
    b07.i[10] = ((const int * ALIGNED(64))a05)[ 7];
    b00.i[11] = ((const int * ALIGNED(64))a05)[ 8];
    b01.i[11] = ((const int * ALIGNED(64))a05)[ 9];
    b02.i[11] = ((const int * ALIGNED(64))a05)[10];
    b03.i[11] = ((const int * ALIGNED(64))a05)[11];
    b04.i[11] = ((const int * ALIGNED(64))a05)[12];
    b05.i[11] = ((const int * ALIGNED(64))a05)[13];
    b06.i[11] = ((const int * ALIGNED(64))a05)[14];
    b07.i[11] = ((const int * ALIGNED(64))a05)[15];

    b00.i[12] = ((const int * ALIGNED(64))a06)[ 0];
    b01.i[12] = ((const int * ALIGNED(64))a06)[ 1];
    b02.i[12] = ((const int * ALIGNED(64))a06)[ 2];
    b03.i[12] = ((const int * ALIGNED(64))a06)[ 3];
    b04.i[12] = ((const int * ALIGNED(64))a06)[ 4];
    b05.i[12] = ((const int * ALIGNED(64))a06)[ 5];
    b06.i[12] = ((const int * ALIGNED(64))a06)[ 6];
    b07.i[12] = ((const int * ALIGNED(64))a06)[ 7];
    b00.i[13] = ((const int * ALIGNED(64))a06)[ 8];
    b01.i[13] = ((const int * ALIGNED(64))a06)[ 9];
    b02.i[13] = ((const int * ALIGNED(64))a06)[10];
    b03.i[13] = ((const int * ALIGNED(64))a06)[11];
    b04.i[13] = ((const int * ALIGNED(64))a06)[12];
    b05.i[13] = ((const int * ALIGNED(64))a06)[13];
    b06.i[13] = ((const int * ALIGNED(64))a06)[14];
    b07.i[13] = ((const int * ALIGNED(64))a06)[15];

    b00.i[14] = ((const int * ALIGNED(64))a07)[ 0];
    b01.i[14] = ((const int * ALIGNED(64))a07)[ 1];
    b02.i[14] = ((const int * ALIGNED(64))a07)[ 2];
    b03.i[14] = ((const int * ALIGNED(64))a07)[ 3];
    b04.i[14] = ((const int * ALIGNED(64))a07)[ 4];
    b05.i[14] = ((const int * ALIGNED(64))a07)[ 5];
    b06.i[14] = ((const int * ALIGNED(64))a07)[ 6];
    b07.i[14] = ((const int * ALIGNED(64))a07)[ 7];
    b00.i[15] = ((const int * ALIGNED(64))a07)[ 8];
    b01.i[15] = ((const int * ALIGNED(64))a07)[ 9];
    b02.i[15] = ((const int * ALIGNED(64))a07)[10];
    b03.i[15] = ((const int * ALIGNED(64))a07)[11];
    b04.i[15] = ((const int * ALIGNED(64))a07)[12];
    b05.i[15] = ((const int * ALIGNED(64))a07)[13];
    b06.i[15] = ((const int * ALIGNED(64))a07)[14];
    b07.i[15] = ((const int * ALIGNED(64))a07)[15];

    b08.i[ 0] = ((const int * ALIGNED(64))a08)[ 0];
    b09.i[ 0] = ((const int * ALIGNED(64))a08)[ 1];
    b10.i[ 0] = ((const int * ALIGNED(64))a08)[ 2];
    b11.i[ 0] = ((const int * ALIGNED(64))a08)[ 3];
    b12.i[ 0] = ((const int * ALIGNED(64))a08)[ 4];
    b13.i[ 0] = ((const int * ALIGNED(64))a08)[ 5];
    b14.i[ 0] = ((const int * ALIGNED(64))a08)[ 6];
    b15.i[ 0] = ((const int * ALIGNED(64))a08)[ 7];
    b08.i[ 1] = ((const int * ALIGNED(64))a08)[ 8];
    b09.i[ 1] = ((const int * ALIGNED(64))a08)[ 9];
    b10.i[ 1] = ((const int * ALIGNED(64))a08)[10];
    b11.i[ 1] = ((const int * ALIGNED(64))a08)[11];
    b12.i[ 1] = ((const int * ALIGNED(64))a08)[12];
    b13.i[ 1] = ((const int * ALIGNED(64))a08)[13];
    b14.i[ 1] = ((const int * ALIGNED(64))a08)[14];
    b15.i[ 1] = ((const int * ALIGNED(64))a08)[15];

    b08.i[ 2] = ((const int * ALIGNED(64))a09)[ 0];
    b09.i[ 2] = ((const int * ALIGNED(64))a09)[ 1];
    b10.i[ 2] = ((const int * ALIGNED(64))a09)[ 2];
    b11.i[ 2] = ((const int * ALIGNED(64))a09)[ 3];
    b12.i[ 2] = ((const int * ALIGNED(64))a09)[ 4];
    b13.i[ 2] = ((const int * ALIGNED(64))a09)[ 5];
    b14.i[ 2] = ((const int * ALIGNED(64))a09)[ 6];
    b15.i[ 2] = ((const int * ALIGNED(64))a09)[ 7];
    b08.i[ 3] = ((const int * ALIGNED(64))a09)[ 8];
    b09.i[ 3] = ((const int * ALIGNED(64))a09)[ 9];
    b10.i[ 3] = ((const int * ALIGNED(64))a09)[10];
    b11.i[ 3] = ((const int * ALIGNED(64))a09)[11];
    b12.i[ 3] = ((const int * ALIGNED(64))a09)[12];
    b13.i[ 3] = ((const int * ALIGNED(64))a09)[13];
    b14.i[ 3] = ((const int * ALIGNED(64))a09)[14];
    b15.i[ 3] = ((const int * ALIGNED(64))a09)[15];

    b08.i[ 4] = ((const int * ALIGNED(64))a10)[ 0];
    b09.i[ 4] = ((const int * ALIGNED(64))a10)[ 1];
    b10.i[ 4] = ((const int * ALIGNED(64))a10)[ 2];
    b11.i[ 4] = ((const int * ALIGNED(64))a10)[ 3];
    b12.i[ 4] = ((const int * ALIGNED(64))a10)[ 4];
    b13.i[ 4] = ((const int * ALIGNED(64))a10)[ 5];
    b14.i[ 4] = ((const int * ALIGNED(64))a10)[ 6];
    b15.i[ 4] = ((const int * ALIGNED(64))a10)[ 7];
    b08.i[ 5] = ((const int * ALIGNED(64))a10)[ 8];
    b09.i[ 5] = ((const int * ALIGNED(64))a10)[ 9];
    b10.i[ 5] = ((const int * ALIGNED(64))a10)[10];
    b11.i[ 5] = ((const int * ALIGNED(64))a10)[11];
    b12.i[ 5] = ((const int * ALIGNED(64))a10)[12];
    b13.i[ 5] = ((const int * ALIGNED(64))a10)[13];
    b14.i[ 5] = ((const int * ALIGNED(64))a10)[14];
    b15.i[ 5] = ((const int * ALIGNED(64))a10)[15];

    b08.i[ 6] = ((const int * ALIGNED(64))a11)[ 0];
    b09.i[ 6] = ((const int * ALIGNED(64))a11)[ 1];
    b10.i[ 6] = ((const int * ALIGNED(64))a11)[ 2];
    b11.i[ 6] = ((const int * ALIGNED(64))a11)[ 3];
    b12.i[ 6] = ((const int * ALIGNED(64))a11)[ 4];
    b13.i[ 6] = ((const int * ALIGNED(64))a11)[ 5];
    b14.i[ 6] = ((const int * ALIGNED(64))a11)[ 6];
    b15.i[ 6] = ((const int * ALIGNED(64))a11)[ 7];
    b08.i[ 7] = ((const int * ALIGNED(64))a11)[ 8];
    b09.i[ 7] = ((const int * ALIGNED(64))a11)[ 9];
    b10.i[ 7] = ((const int * ALIGNED(64))a11)[10];
    b11.i[ 7] = ((const int * ALIGNED(64))a11)[11];
    b12.i[ 7] = ((const int * ALIGNED(64))a11)[12];
    b13.i[ 7] = ((const int * ALIGNED(64))a11)[13];
    b14.i[ 7] = ((const int * ALIGNED(64))a11)[14];
    b15.i[ 7] = ((const int * ALIGNED(64))a11)[15];

    b08.i[ 8] = ((const int * ALIGNED(64))a12)[ 0];
    b09.i[ 8] = ((const int * ALIGNED(64))a12)[ 1];
    b10.i[ 8] = ((const int * ALIGNED(64))a12)[ 2];
    b11.i[ 8] = ((const int * ALIGNED(64))a12)[ 3];
    b12.i[ 8] = ((const int * ALIGNED(64))a12)[ 4];
    b13.i[ 8] = ((const int * ALIGNED(64))a12)[ 5];
    b14.i[ 8] = ((const int * ALIGNED(64))a12)[ 6];
    b15.i[ 8] = ((const int * ALIGNED(64))a12)[ 7];
    b08.i[ 9] = ((const int * ALIGNED(64))a12)[ 8];
    b09.i[ 9] = ((const int * ALIGNED(64))a12)[ 9];
    b10.i[ 9] = ((const int * ALIGNED(64))a12)[10];
    b11.i[ 9] = ((const int * ALIGNED(64))a12)[11];
    b12.i[ 9] = ((const int * ALIGNED(64))a12)[12];
    b13.i[ 9] = ((const int * ALIGNED(64))a12)[13];
    b14.i[ 9] = ((const int * ALIGNED(64))a12)[14];
    b15.i[ 9] = ((const int * ALIGNED(64))a12)[15];

    b08.i[10] = ((const int * ALIGNED(64))a13)[ 0];
    b09.i[10] = ((const int * ALIGNED(64))a13)[ 1];
    b10.i[10] = ((const int * ALIGNED(64))a13)[ 2];
    b11.i[10] = ((const int * ALIGNED(64))a13)[ 3];
    b12.i[10] = ((const int * ALIGNED(64))a13)[ 4];
    b13.i[10] = ((const int * ALIGNED(64))a13)[ 5];
    b14.i[10] = ((const int * ALIGNED(64))a13)[ 6];
    b15.i[10] = ((const int * ALIGNED(64))a13)[ 7];
    b08.i[11] = ((const int * ALIGNED(64))a13)[ 8];
    b09.i[11] = ((const int * ALIGNED(64))a13)[ 9];
    b10.i[11] = ((const int * ALIGNED(64))a13)[10];
    b11.i[11] = ((const int * ALIGNED(64))a13)[11];
    b12.i[11] = ((const int * ALIGNED(64))a13)[12];
    b13.i[11] = ((const int * ALIGNED(64))a13)[13];
    b14.i[11] = ((const int * ALIGNED(64))a13)[14];
    b15.i[11] = ((const int * ALIGNED(64))a13)[15];

    b08.i[12] = ((const int * ALIGNED(64))a14)[ 0];
    b09.i[12] = ((const int * ALIGNED(64))a14)[ 1];
    b10.i[12] = ((const int * ALIGNED(64))a14)[ 2];
    b11.i[12] = ((const int * ALIGNED(64))a14)[ 3];
    b12.i[12] = ((const int * ALIGNED(64))a14)[ 4];
    b13.i[12] = ((const int * ALIGNED(64))a14)[ 5];
    b14.i[12] = ((const int * ALIGNED(64))a14)[ 6];
    b15.i[12] = ((const int * ALIGNED(64))a14)[ 7];
    b08.i[13] = ((const int * ALIGNED(64))a14)[ 8];
    b09.i[13] = ((const int * ALIGNED(64))a14)[ 9];
    b10.i[13] = ((const int * ALIGNED(64))a14)[10];
    b11.i[13] = ((const int * ALIGNED(64))a14)[11];
    b12.i[13] = ((const int * ALIGNED(64))a14)[12];
    b13.i[13] = ((const int * ALIGNED(64))a14)[13];
    b14.i[13] = ((const int * ALIGNED(64))a14)[14];
    b15.i[13] = ((const int * ALIGNED(64))a14)[15];

    b08.i[14] = ((const int * ALIGNED(64))a15)[ 0];
    b09.i[14] = ((const int * ALIGNED(64))a15)[ 1];
    b10.i[14] = ((const int * ALIGNED(64))a15)[ 2];
    b11.i[14] = ((const int * ALIGNED(64))a15)[ 3];
    b12.i[14] = ((const int * ALIGNED(64))a15)[ 4];
    b13.i[14] = ((const int * ALIGNED(64))a15)[ 5];
    b14.i[14] = ((const int * ALIGNED(64))a15)[ 6];
    b15.i[14] = ((const int * ALIGNED(64))a15)[ 7];
    b08.i[15] = ((const int * ALIGNED(64))a15)[ 8];
    b09.i[15] = ((const int * ALIGNED(64))a15)[ 9];
    b10.i[15] = ((const int * ALIGNED(64))a15)[10];
    b11.i[15] = ((const int * ALIGNED(64))a15)[11];
    b12.i[15] = ((const int * ALIGNED(64))a15)[12];
    b13.i[15] = ((const int * ALIGNED(64))a15)[13];
    b14.i[15] = ((const int * ALIGNED(64))a15)[14];
    b15.i[15] = ((const int * ALIGNED(64))a15)[15];
  }

  inline void store_16x1_tr( const v16 &a,
			     void *a00, void *a01, void *a02, void *a03,
			     void *a04, void *a05, void *a06, void *a07,
			     void *a08, void *a09, void *a10, void *a11,
			     void *a12, void *a13, void *a14, void *a15 )
  {
    ((int *)a00)[0] = a.i[ 0];
    ((int *)a01)[0] = a.i[ 1];
    ((int *)a02)[0] = a.i[ 2];
    ((int *)a03)[0] = a.i[ 3];
    ((int *)a04)[0] = a.i[ 4];
    ((int *)a05)[0] = a.i[ 5];
    ((int *)a06)[0] = a.i[ 6];
    ((int *)a07)[0] = a.i[ 7];
    ((int *)a08)[0] = a.i[ 8];
    ((int *)a09)[0] = a.i[ 9];
    ((int *)a10)[0] = a.i[10];
    ((int *)a11)[0] = a.i[11];
    ((int *)a12)[0] = a.i[12];
    ((int *)a13)[0] = a.i[13];
    ((int *)a14)[0] = a.i[14];
    ((int *)a15)[0] = a.i[15];
  }

  inline void store_16x2_tr( const v16 &a, const v16 &b,
			     void * ALIGNED(8) a00, void * ALIGNED(8) a01,
			     void * ALIGNED(8) a02, void * ALIGNED(8) a03,
			     void * ALIGNED(8) a04, void * ALIGNED(8) a05,
			     void * ALIGNED(8) a06, void * ALIGNED(8) a07,
			     void * ALIGNED(8) a08, void * ALIGNED(8) a09,
			     void * ALIGNED(8) a10, void * ALIGNED(8) a11,
			     void * ALIGNED(8) a12, void * ALIGNED(8) a13,
			     void * ALIGNED(8) a14, void * ALIGNED(8) a15 )
  {
    ((int * ALIGNED(8))a00)[0] = a.i[ 0];
    ((int * ALIGNED(8))a00)[1] = b.i[ 0];

    ((int * ALIGNED(8))a01)[0] = a.i[ 1];
    ((int * ALIGNED(8))a01)[1] = b.i[ 1];

    ((int * ALIGNED(8))a02)[0] = a.i[ 2];
    ((int * ALIGNED(8))a02)[1] = b.i[ 2];

    ((int * ALIGNED(8))a03)[0] = a.i[ 3];
    ((int * ALIGNED(8))a03)[1] = b.i[ 3];

    ((int * ALIGNED(8))a04)[0] = a.i[ 4];
    ((int * ALIGNED(8))a04)[1] = b.i[ 4];

    ((int * ALIGNED(8))a05)[0] = a.i[ 5];
    ((int * ALIGNED(8))a05)[1] = b.i[ 5];

    ((int * ALIGNED(8))a06)[0] = a.i[ 6];
    ((int * ALIGNED(8))a06)[1] = b.i[ 6];

    ((int * ALIGNED(8))a07)[0] = a.i[ 7];
    ((int * ALIGNED(8))a07)[1] = b.i[ 7];

    ((int * ALIGNED(8))a08)[0] = a.i[ 8];
    ((int * ALIGNED(8))a08)[1] = b.i[ 8];

    ((int * ALIGNED(8))a09)[0] = a.i[ 9];
    ((int * ALIGNED(8))a09)[1] = b.i[ 9];

    ((int * ALIGNED(8))a10)[0] = a.i[10];
    ((int * ALIGNED(8))a10)[1] = b.i[10];

    ((int * ALIGNED(8))a11)[0] = a.i[11];
    ((int * ALIGNED(8))a11)[1] = b.i[11];

    ((int * ALIGNED(8))a12)[0] = a.i[12];
    ((int * ALIGNED(8))a12)[1] = b.i[12];

    ((int * ALIGNED(8))a13)[0] = a.i[13];
    ((int * ALIGNED(8))a13)[1] = b.i[13];

    ((int * ALIGNED(8))a14)[0] = a.i[14];
    ((int * ALIGNED(8))a14)[1] = b.i[14];

    ((int * ALIGNED(8))a15)[0] = a.i[15];
    ((int * ALIGNED(8))a15)[1] = b.i[15];
  }

  inline void store_16x3_tr( const v16 &a, const v16 &b, const v16 &c,
			     void * ALIGNED(64) a00, void * ALIGNED(64) a01,
			     void * ALIGNED(64) a02, void * ALIGNED(64) a03,
			     void * ALIGNED(64) a04, void * ALIGNED(64) a05,
			     void * ALIGNED(64) a06, void * ALIGNED(64) a07,
			     void * ALIGNED(64) a08, void * ALIGNED(64) a09,
			     void * ALIGNED(64) a10, void * ALIGNED(64) a11,
			     void * ALIGNED(64) a12, void * ALIGNED(64) a13,
			     void * ALIGNED(64) a14, void * ALIGNED(64) a15 )
  {
    ((int * ALIGNED(64))a00)[0] = a.i[ 0];
    ((int * ALIGNED(64))a00)[1] = b.i[ 0];
    ((int * ALIGNED(64))a00)[2] = c.i[ 0];

    ((int * ALIGNED(64))a01)[0] = a.i[ 1];
    ((int * ALIGNED(64))a01)[1] = b.i[ 1];
    ((int * ALIGNED(64))a01)[2] = c.i[ 1];

    ((int * ALIGNED(64))a02)[0] = a.i[ 2];
    ((int * ALIGNED(64))a02)[1] = b.i[ 2];
    ((int * ALIGNED(64))a02)[2] = c.i[ 2];

    ((int * ALIGNED(64))a03)[0] = a.i[ 3];
    ((int * ALIGNED(64))a03)[1] = b.i[ 3];
    ((int * ALIGNED(64))a03)[2] = c.i[ 3];

    ((int * ALIGNED(64))a04)[0] = a.i[ 4];
    ((int * ALIGNED(64))a04)[1] = b.i[ 4];
    ((int * ALIGNED(64))a04)[2] = c.i[ 4];

    ((int * ALIGNED(64))a05)[0] = a.i[ 5];
    ((int * ALIGNED(64))a05)[1] = b.i[ 5];
    ((int * ALIGNED(64))a05)[2] = c.i[ 5];

    ((int * ALIGNED(64))a06)[0] = a.i[ 6];
    ((int * ALIGNED(64))a06)[1] = b.i[ 6];
    ((int * ALIGNED(64))a06)[2] = c.i[ 6];

    ((int * ALIGNED(64))a07)[0] = a.i[ 7];
    ((int * ALIGNED(64))a07)[1] = b.i[ 7];
    ((int * ALIGNED(64))a07)[2] = c.i[ 7];

    ((int * ALIGNED(64))a08)[0] = a.i[ 8];
    ((int * ALIGNED(64))a08)[1] = b.i[ 8];
    ((int * ALIGNED(64))a08)[2] = c.i[ 8];

    ((int * ALIGNED(64))a09)[0] = a.i[ 9];
    ((int * ALIGNED(64))a09)[1] = b.i[ 9];
    ((int * ALIGNED(64))a09)[2] = c.i[ 9];

    ((int * ALIGNED(64))a10)[0] = a.i[10];
    ((int * ALIGNED(64))a10)[1] = b.i[10];
    ((int * ALIGNED(64))a10)[2] = c.i[10];

    ((int * ALIGNED(64))a11)[0] = a.i[11];
    ((int * ALIGNED(64))a11)[1] = b.i[11];
    ((int * ALIGNED(64))a11)[2] = c.i[11];

    ((int * ALIGNED(64))a12)[0] = a.i[12];
    ((int * ALIGNED(64))a12)[1] = b.i[12];
    ((int * ALIGNED(64))a12)[2] = c.i[12];

    ((int * ALIGNED(64))a13)[0] = a.i[13];
    ((int * ALIGNED(64))a13)[1] = b.i[13];
    ((int * ALIGNED(64))a13)[2] = c.i[13];

    ((int * ALIGNED(64))a14)[0] = a.i[14];
    ((int * ALIGNED(64))a14)[1] = b.i[14];
    ((int * ALIGNED(64))a14)[2] = c.i[14];

    ((int * ALIGNED(64))a15)[0] = a.i[15];
    ((int * ALIGNED(64))a15)[1] = b.i[15];
    ((int * ALIGNED(64))a15)[2] = c.i[15];
  }

  inline void store_16x4_tr( const v16 &a, const v16 &b, const v16 &c, const v16 &d,
			     void * ALIGNED(64) a00, void * ALIGNED(64) a01,
			     void * ALIGNED(64) a02, void * ALIGNED(64) a03,
			     void * ALIGNED(64) a04, void * ALIGNED(64) a05,
			     void * ALIGNED(64) a06, void * ALIGNED(64) a07,
			     void * ALIGNED(64) a08, void * ALIGNED(64) a09,
			     void * ALIGNED(64) a10, void * ALIGNED(64) a11,
			     void * ALIGNED(64) a12, void * ALIGNED(64) a13,
			     void * ALIGNED(64) a14, void * ALIGNED(64) a15 )
  {
    ((int * ALIGNED(64))a00)[0] = a.i[ 0];
    ((int * ALIGNED(64))a00)[1] = b.i[ 0];
    ((int * ALIGNED(64))a00)[2] = c.i[ 0];
    ((int * ALIGNED(64))a00)[3] = d.i[ 0];

    ((int * ALIGNED(64))a01)[0] = a.i[ 1];
    ((int * ALIGNED(64))a01)[1] = b.i[ 1];
    ((int * ALIGNED(64))a01)[2] = c.i[ 1];
    ((int * ALIGNED(64))a01)[3] = d.i[ 1];

    ((int * ALIGNED(64))a02)[0] = a.i[ 2];
    ((int * ALIGNED(64))a02)[1] = b.i[ 2];
    ((int * ALIGNED(64))a02)[2] = c.i[ 2];
    ((int * ALIGNED(64))a02)[3] = d.i[ 2];

    ((int * ALIGNED(64))a03)[0] = a.i[ 3];
    ((int * ALIGNED(64))a03)[1] = b.i[ 3];
    ((int * ALIGNED(64))a03)[2] = c.i[ 3];
    ((int * ALIGNED(64))a03)[3] = d.i[ 3];

    ((int * ALIGNED(64))a04)[0] = a.i[ 4];
    ((int * ALIGNED(64))a04)[1] = b.i[ 4];
    ((int * ALIGNED(64))a04)[2] = c.i[ 4];
    ((int * ALIGNED(64))a04)[3] = d.i[ 4];

    ((int * ALIGNED(64))a05)[0] = a.i[ 5];
    ((int * ALIGNED(64))a05)[1] = b.i[ 5];
    ((int * ALIGNED(64))a05)[2] = c.i[ 5];
    ((int * ALIGNED(64))a05)[3] = d.i[ 5];

    ((int * ALIGNED(64))a06)[0] = a.i[ 6];
    ((int * ALIGNED(64))a06)[1] = b.i[ 6];
    ((int * ALIGNED(64))a06)[2] = c.i[ 6];
    ((int * ALIGNED(64))a06)[3] = d.i[ 6];

    ((int * ALIGNED(64))a07)[0] = a.i[ 7];
    ((int * ALIGNED(64))a07)[1] = b.i[ 7];
    ((int * ALIGNED(64))a07)[2] = c.i[ 7];
    ((int * ALIGNED(64))a07)[3] = d.i[ 7];

    ((int * ALIGNED(64))a08)[0] = a.i[ 8];
    ((int * ALIGNED(64))a08)[1] = b.i[ 8];
    ((int * ALIGNED(64))a08)[2] = c.i[ 8];
    ((int * ALIGNED(64))a08)[3] = d.i[ 8];

    ((int * ALIGNED(64))a09)[0] = a.i[ 9];
    ((int * ALIGNED(64))a09)[1] = b.i[ 9];
    ((int * ALIGNED(64))a09)[2] = c.i[ 9];
    ((int * ALIGNED(64))a09)[3] = d.i[ 9];

    ((int * ALIGNED(64))a10)[0] = a.i[10];
    ((int * ALIGNED(64))a10)[1] = b.i[10];
    ((int * ALIGNED(64))a10)[2] = c.i[10];
    ((int * ALIGNED(64))a10)[3] = d.i[10];

    ((int * ALIGNED(64))a11)[0] = a.i[11];
    ((int * ALIGNED(64))a11)[1] = b.i[11];
    ((int * ALIGNED(64))a11)[2] = c.i[11];
    ((int * ALIGNED(64))a11)[3] = d.i[11];

    ((int * ALIGNED(64))a12)[0] = a.i[12];
    ((int * ALIGNED(64))a12)[1] = b.i[12];
    ((int * ALIGNED(64))a12)[2] = c.i[12];
    ((int * ALIGNED(64))a12)[3] = d.i[12];

    ((int * ALIGNED(64))a13)[0] = a.i[13];
    ((int * ALIGNED(64))a13)[1] = b.i[13];
    ((int * ALIGNED(64))a13)[2] = c.i[13];
    ((int * ALIGNED(64))a13)[3] = d.i[13];

    ((int * ALIGNED(64))a14)[0] = a.i[14];
    ((int * ALIGNED(64))a14)[1] = b.i[14];
    ((int * ALIGNED(64))a14)[2] = c.i[14];
    ((int * ALIGNED(64))a14)[3] = d.i[14];

    ((int * ALIGNED(64))a15)[0] = a.i[15];
    ((int * ALIGNED(64))a15)[1] = b.i[15];
    ((int * ALIGNED(64))a15)[2] = c.i[15];
    ((int * ALIGNED(64))a15)[3] = d.i[15];
  }

  inline void store_16x8_tr( const v16 &a, const v16 &b, const v16 &c, const v16 &d,
			     const v16 &e, const v16 &f, const v16 &g, const v16 &h,
			     void * ALIGNED(64) a00, void * ALIGNED(64) a01,
			     void * ALIGNED(64) a02, void * ALIGNED(64) a03,
			     void * ALIGNED(64) a04, void * ALIGNED(64) a05,
			     void * ALIGNED(64) a06, void * ALIGNED(64) a07,
			     void * ALIGNED(64) a08, void * ALIGNED(64) a09,
			     void * ALIGNED(64) a10, void * ALIGNED(64) a11,
			     void * ALIGNED(64) a12, void * ALIGNED(64) a13,
			     void * ALIGNED(64) a14, void * ALIGNED(64) a15 )
  {
    ((int * ALIGNED(64))a00)[0] = a.i[ 0];
    ((int * ALIGNED(64))a00)[1] = b.i[ 0];
    ((int * ALIGNED(64))a00)[2] = c.i[ 0];
    ((int * ALIGNED(64))a00)[3] = d.i[ 0];
    ((int * ALIGNED(64))a00)[4] = e.i[ 0];
    ((int * ALIGNED(64))a00)[5] = f.i[ 0];
    ((int * ALIGNED(64))a00)[6] = g.i[ 0];
    ((int * ALIGNED(64))a00)[7] = h.i[ 0];

    ((int * ALIGNED(64))a01)[0] = a.i[ 1];
    ((int * ALIGNED(64))a01)[1] = b.i[ 1];
    ((int * ALIGNED(64))a01)[2] = c.i[ 1];
    ((int * ALIGNED(64))a01)[3] = d.i[ 1];
    ((int * ALIGNED(64))a01)[4] = e.i[ 1];
    ((int * ALIGNED(64))a01)[5] = f.i[ 1];
    ((int * ALIGNED(64))a01)[6] = g.i[ 1];
    ((int * ALIGNED(64))a01)[7] = h.i[ 1];

    ((int * ALIGNED(64))a02)[0] = a.i[ 2];
    ((int * ALIGNED(64))a02)[1] = b.i[ 2];
    ((int * ALIGNED(64))a02)[2] = c.i[ 2];
    ((int * ALIGNED(64))a02)[3] = d.i[ 2];
    ((int * ALIGNED(64))a02)[4] = e.i[ 2];
    ((int * ALIGNED(64))a02)[5] = f.i[ 2];
    ((int * ALIGNED(64))a02)[6] = g.i[ 2];
    ((int * ALIGNED(64))a02)[7] = h.i[ 2];

    ((int * ALIGNED(64))a03)[0] = a.i[ 3];
    ((int * ALIGNED(64))a03)[1] = b.i[ 3];
    ((int * ALIGNED(64))a03)[2] = c.i[ 3];
    ((int * ALIGNED(64))a03)[3] = d.i[ 3];
    ((int * ALIGNED(64))a03)[4] = e.i[ 3];
    ((int * ALIGNED(64))a03)[5] = f.i[ 3];
    ((int * ALIGNED(64))a03)[6] = g.i[ 3];
    ((int * ALIGNED(64))a03)[7] = h.i[ 3];

    ((int * ALIGNED(64))a04)[0] = a.i[ 4];
    ((int * ALIGNED(64))a04)[1] = b.i[ 4];
    ((int * ALIGNED(64))a04)[2] = c.i[ 4];
    ((int * ALIGNED(64))a04)[3] = d.i[ 4];
    ((int * ALIGNED(64))a04)[4] = e.i[ 4];
    ((int * ALIGNED(64))a04)[5] = f.i[ 4];
    ((int * ALIGNED(64))a04)[6] = g.i[ 4];
    ((int * ALIGNED(64))a04)[7] = h.i[ 4];

    ((int * ALIGNED(64))a05)[0] = a.i[ 5];
    ((int * ALIGNED(64))a05)[1] = b.i[ 5];
    ((int * ALIGNED(64))a05)[2] = c.i[ 5];
    ((int * ALIGNED(64))a05)[3] = d.i[ 5];
    ((int * ALIGNED(64))a05)[4] = e.i[ 5];
    ((int * ALIGNED(64))a05)[5] = f.i[ 5];
    ((int * ALIGNED(64))a05)[6] = g.i[ 5];
    ((int * ALIGNED(64))a05)[7] = h.i[ 5];

    ((int * ALIGNED(64))a06)[0] = a.i[ 6];
    ((int * ALIGNED(64))a06)[1] = b.i[ 6];
    ((int * ALIGNED(64))a06)[2] = c.i[ 6];
    ((int * ALIGNED(64))a06)[3] = d.i[ 6];
    ((int * ALIGNED(64))a06)[4] = e.i[ 6];
    ((int * ALIGNED(64))a06)[5] = f.i[ 6];
    ((int * ALIGNED(64))a06)[6] = g.i[ 6];
    ((int * ALIGNED(64))a06)[7] = h.i[ 6];

    ((int * ALIGNED(64))a07)[0] = a.i[ 7];
    ((int * ALIGNED(64))a07)[1] = b.i[ 7];
    ((int * ALIGNED(64))a07)[2] = c.i[ 7];
    ((int * ALIGNED(64))a07)[3] = d.i[ 7];
    ((int * ALIGNED(64))a07)[4] = e.i[ 7];
    ((int * ALIGNED(64))a07)[5] = f.i[ 7];
    ((int * ALIGNED(64))a07)[6] = g.i[ 7];
    ((int * ALIGNED(64))a07)[7] = h.i[ 7];

    ((int * ALIGNED(64))a08)[0] = a.i[ 8];
    ((int * ALIGNED(64))a08)[1] = b.i[ 8];
    ((int * ALIGNED(64))a08)[2] = c.i[ 8];
    ((int * ALIGNED(64))a08)[3] = d.i[ 8];
    ((int * ALIGNED(64))a08)[4] = e.i[ 8];
    ((int * ALIGNED(64))a08)[5] = f.i[ 8];
    ((int * ALIGNED(64))a08)[6] = g.i[ 8];
    ((int * ALIGNED(64))a08)[7] = h.i[ 8];

    ((int * ALIGNED(64))a09)[0] = a.i[ 9];
    ((int * ALIGNED(64))a09)[1] = b.i[ 9];
    ((int * ALIGNED(64))a09)[2] = c.i[ 9];
    ((int * ALIGNED(64))a09)[3] = d.i[ 9];
    ((int * ALIGNED(64))a09)[4] = e.i[ 9];
    ((int * ALIGNED(64))a09)[5] = f.i[ 9];
    ((int * ALIGNED(64))a09)[6] = g.i[ 9];
    ((int * ALIGNED(64))a09)[7] = h.i[ 9];

    ((int * ALIGNED(64))a10)[0] = a.i[10];
    ((int * ALIGNED(64))a10)[1] = b.i[10];
    ((int * ALIGNED(64))a10)[2] = c.i[10];
    ((int * ALIGNED(64))a10)[3] = d.i[10];
    ((int * ALIGNED(64))a10)[4] = e.i[10];
    ((int * ALIGNED(64))a10)[5] = f.i[10];
    ((int * ALIGNED(64))a10)[6] = g.i[10];
    ((int * ALIGNED(64))a10)[7] = h.i[10];

    ((int * ALIGNED(64))a11)[0] = a.i[11];
    ((int * ALIGNED(64))a11)[1] = b.i[11];
    ((int * ALIGNED(64))a11)[2] = c.i[11];
    ((int * ALIGNED(64))a11)[3] = d.i[11];
    ((int * ALIGNED(64))a11)[4] = e.i[11];
    ((int * ALIGNED(64))a11)[5] = f.i[11];
    ((int * ALIGNED(64))a11)[6] = g.i[11];
    ((int * ALIGNED(64))a11)[7] = h.i[11];

    ((int * ALIGNED(64))a12)[0] = a.i[12];
    ((int * ALIGNED(64))a12)[1] = b.i[12];
    ((int * ALIGNED(64))a12)[2] = c.i[12];
    ((int * ALIGNED(64))a12)[3] = d.i[12];
    ((int * ALIGNED(64))a12)[4] = e.i[12];
    ((int * ALIGNED(64))a12)[5] = f.i[12];
    ((int * ALIGNED(64))a12)[6] = g.i[12];
    ((int * ALIGNED(64))a12)[7] = h.i[12];

    ((int * ALIGNED(64))a13)[0] = a.i[13];
    ((int * ALIGNED(64))a13)[1] = b.i[13];
    ((int * ALIGNED(64))a13)[2] = c.i[13];
    ((int * ALIGNED(64))a13)[3] = d.i[13];
    ((int * ALIGNED(64))a13)[4] = e.i[13];
    ((int * ALIGNED(64))a13)[5] = f.i[13];
    ((int * ALIGNED(64))a13)[6] = g.i[13];
    ((int * ALIGNED(64))a13)[7] = h.i[13];

    ((int * ALIGNED(64))a14)[0] = a.i[14];
    ((int * ALIGNED(64))a14)[1] = b.i[14];
    ((int * ALIGNED(64))a14)[2] = c.i[14];
    ((int * ALIGNED(64))a14)[3] = d.i[14];
    ((int * ALIGNED(64))a14)[4] = e.i[14];
    ((int * ALIGNED(64))a14)[5] = f.i[14];
    ((int * ALIGNED(64))a14)[6] = g.i[14];
    ((int * ALIGNED(64))a14)[7] = h.i[14];

    ((int * ALIGNED(64))a15)[0] = a.i[15];
    ((int * ALIGNED(64))a15)[1] = b.i[15];
    ((int * ALIGNED(64))a15)[2] = c.i[15];
    ((int * ALIGNED(64))a15)[3] = d.i[15];
    ((int * ALIGNED(64))a15)[4] = e.i[15];
    ((int * ALIGNED(64))a15)[5] = f.i[15];
    ((int * ALIGNED(64))a15)[6] = g.i[15];
    ((int * ALIGNED(64))a15)[7] = h.i[15];
  }

  inline void store_16x16_tr( const v16 &b00, const v16 &b01, const v16 &b02, const v16 &b03,
			      const v16 &b04, const v16 &b05, const v16 &b06, const v16 &b07,
			      const v16 &b08, const v16 &b09, const v16 &b10, const v16 &b11,
			      const v16 &b12, const v16 &b13, const v16 &b14, const v16 &b15,
			      void * ALIGNED(64) a00, void * ALIGNED(64) a01,
			      void * ALIGNED(64) a02, void * ALIGNED(64) a03,
			      void * ALIGNED(64) a04, void * ALIGNED(64) a05,
			      void * ALIGNED(64) a06, void * ALIGNED(64) a07,
			      void * ALIGNED(64) a08, void * ALIGNED(64) a09,
			      void * ALIGNED(64) a10, void * ALIGNED(64) a11,
			      void * ALIGNED(64) a12, void * ALIGNED(64) a13,
			      void * ALIGNED(64) a14, void * ALIGNED(64) a15 )
  {
    ((int * ALIGNED(64))a00)[ 0] = b00.i[ 0];
    ((int * ALIGNED(64))a00)[ 1] = b01.i[ 0];
    ((int * ALIGNED(64))a00)[ 2] = b02.i[ 0];
    ((int * ALIGNED(64))a00)[ 3] = b03.i[ 0];
    ((int * ALIGNED(64))a00)[ 4] = b04.i[ 0];
    ((int * ALIGNED(64))a00)[ 5] = b05.i[ 0];
    ((int * ALIGNED(64))a00)[ 6] = b06.i[ 0];
    ((int * ALIGNED(64))a00)[ 7] = b07.i[ 0];
    ((int * ALIGNED(64))a00)[ 8] = b08.i[ 0];
    ((int * ALIGNED(64))a00)[ 9] = b09.i[ 0];
    ((int * ALIGNED(64))a00)[10] = b10.i[ 0];
    ((int * ALIGNED(64))a00)[11] = b11.i[ 0];
    ((int * ALIGNED(64))a00)[12] = b12.i[ 0];
    ((int * ALIGNED(64))a00)[13] = b13.i[ 0];
    ((int * ALIGNED(64))a00)[14] = b14.i[ 0];
    ((int * ALIGNED(64))a00)[15] = b15.i[ 0];

    ((int * ALIGNED(64))a01)[ 0] = b00.i[ 1];
    ((int * ALIGNED(64))a01)[ 1] = b01.i[ 1];
    ((int * ALIGNED(64))a01)[ 2] = b02.i[ 1];
    ((int * ALIGNED(64))a01)[ 3] = b03.i[ 1];
    ((int * ALIGNED(64))a01)[ 4] = b04.i[ 1];
    ((int * ALIGNED(64))a01)[ 5] = b05.i[ 1];
    ((int * ALIGNED(64))a01)[ 6] = b06.i[ 1];
    ((int * ALIGNED(64))a01)[ 7] = b07.i[ 1];
    ((int * ALIGNED(64))a01)[ 8] = b08.i[ 1];
    ((int * ALIGNED(64))a01)[ 9] = b09.i[ 1];
    ((int * ALIGNED(64))a01)[10] = b10.i[ 1];
    ((int * ALIGNED(64))a01)[11] = b11.i[ 1];
    ((int * ALIGNED(64))a01)[12] = b12.i[ 1];
    ((int * ALIGNED(64))a01)[13] = b13.i[ 1];
    ((int * ALIGNED(64))a01)[14] = b14.i[ 1];
    ((int * ALIGNED(64))a01)[15] = b15.i[ 1];

    ((int * ALIGNED(64))a02)[ 0] = b00.i[ 2];
    ((int * ALIGNED(64))a02)[ 1] = b01.i[ 2];
    ((int * ALIGNED(64))a02)[ 2] = b02.i[ 2];
    ((int * ALIGNED(64))a02)[ 3] = b03.i[ 2];
    ((int * ALIGNED(64))a02)[ 4] = b04.i[ 2];
    ((int * ALIGNED(64))a02)[ 5] = b05.i[ 2];
    ((int * ALIGNED(64))a02)[ 6] = b06.i[ 2];
    ((int * ALIGNED(64))a02)[ 7] = b07.i[ 2];
    ((int * ALIGNED(64))a02)[ 8] = b08.i[ 2];
    ((int * ALIGNED(64))a02)[ 9] = b09.i[ 2];
    ((int * ALIGNED(64))a02)[10] = b10.i[ 2];
    ((int * ALIGNED(64))a02)[11] = b11.i[ 2];
    ((int * ALIGNED(64))a02)[12] = b12.i[ 2];
    ((int * ALIGNED(64))a02)[13] = b13.i[ 2];
    ((int * ALIGNED(64))a02)[14] = b14.i[ 2];
    ((int * ALIGNED(64))a02)[15] = b15.i[ 2];

    ((int * ALIGNED(64))a03)[ 0] = b00.i[ 3];
    ((int * ALIGNED(64))a03)[ 1] = b01.i[ 3];
    ((int * ALIGNED(64))a03)[ 2] = b02.i[ 3];
    ((int * ALIGNED(64))a03)[ 3] = b03.i[ 3];
    ((int * ALIGNED(64))a03)[ 4] = b04.i[ 3];
    ((int * ALIGNED(64))a03)[ 5] = b05.i[ 3];
    ((int * ALIGNED(64))a03)[ 6] = b06.i[ 3];
    ((int * ALIGNED(64))a03)[ 7] = b07.i[ 3];
    ((int * ALIGNED(64))a03)[ 8] = b08.i[ 3];
    ((int * ALIGNED(64))a03)[ 9] = b09.i[ 3];
    ((int * ALIGNED(64))a03)[10] = b10.i[ 3];
    ((int * ALIGNED(64))a03)[11] = b11.i[ 3];
    ((int * ALIGNED(64))a03)[12] = b12.i[ 3];
    ((int * ALIGNED(64))a03)[13] = b13.i[ 3];
    ((int * ALIGNED(64))a03)[14] = b14.i[ 3];
    ((int * ALIGNED(64))a03)[15] = b15.i[ 3];

    ((int * ALIGNED(64))a04)[ 0] = b00.i[ 4];
    ((int * ALIGNED(64))a04)[ 1] = b01.i[ 4];
    ((int * ALIGNED(64))a04)[ 2] = b02.i[ 4];
    ((int * ALIGNED(64))a04)[ 3] = b03.i[ 4];
    ((int * ALIGNED(64))a04)[ 4] = b04.i[ 4];
    ((int * ALIGNED(64))a04)[ 5] = b05.i[ 4];
    ((int * ALIGNED(64))a04)[ 6] = b06.i[ 4];
    ((int * ALIGNED(64))a04)[ 7] = b07.i[ 4];
    ((int * ALIGNED(64))a04)[ 8] = b08.i[ 4];
    ((int * ALIGNED(64))a04)[ 9] = b09.i[ 4];
    ((int * ALIGNED(64))a04)[10] = b10.i[ 4];
    ((int * ALIGNED(64))a04)[11] = b11.i[ 4];
    ((int * ALIGNED(64))a04)[12] = b12.i[ 4];
    ((int * ALIGNED(64))a04)[13] = b13.i[ 4];
    ((int * ALIGNED(64))a04)[14] = b14.i[ 4];
    ((int * ALIGNED(64))a04)[15] = b15.i[ 4];

    ((int * ALIGNED(64))a05)[ 0] = b00.i[ 5];
    ((int * ALIGNED(64))a05)[ 1] = b01.i[ 5];
    ((int * ALIGNED(64))a05)[ 2] = b02.i[ 5];
    ((int * ALIGNED(64))a05)[ 3] = b03.i[ 5];
    ((int * ALIGNED(64))a05)[ 4] = b04.i[ 5];
    ((int * ALIGNED(64))a05)[ 5] = b05.i[ 5];
    ((int * ALIGNED(64))a05)[ 6] = b06.i[ 5];
    ((int * ALIGNED(64))a05)[ 7] = b07.i[ 5];
    ((int * ALIGNED(64))a05)[ 8] = b08.i[ 5];
    ((int * ALIGNED(64))a05)[ 9] = b09.i[ 5];
    ((int * ALIGNED(64))a05)[10] = b10.i[ 5];
    ((int * ALIGNED(64))a05)[11] = b11.i[ 5];
    ((int * ALIGNED(64))a05)[12] = b12.i[ 5];
    ((int * ALIGNED(64))a05)[13] = b13.i[ 5];
    ((int * ALIGNED(64))a05)[14] = b14.i[ 5];
    ((int * ALIGNED(64))a05)[15] = b15.i[ 5];

    ((int * ALIGNED(64))a06)[ 0] = b00.i[ 6];
    ((int * ALIGNED(64))a06)[ 1] = b01.i[ 6];
    ((int * ALIGNED(64))a06)[ 2] = b02.i[ 6];
    ((int * ALIGNED(64))a06)[ 3] = b03.i[ 6];
    ((int * ALIGNED(64))a06)[ 4] = b04.i[ 6];
    ((int * ALIGNED(64))a06)[ 5] = b05.i[ 6];
    ((int * ALIGNED(64))a06)[ 6] = b06.i[ 6];
    ((int * ALIGNED(64))a06)[ 7] = b07.i[ 6];
    ((int * ALIGNED(64))a06)[ 8] = b08.i[ 6];
    ((int * ALIGNED(64))a06)[ 9] = b09.i[ 6];
    ((int * ALIGNED(64))a06)[10] = b10.i[ 6];
    ((int * ALIGNED(64))a06)[11] = b11.i[ 6];
    ((int * ALIGNED(64))a06)[12] = b12.i[ 6];
    ((int * ALIGNED(64))a06)[13] = b13.i[ 6];
    ((int * ALIGNED(64))a06)[14] = b14.i[ 6];
    ((int * ALIGNED(64))a06)[15] = b15.i[ 6];

    ((int * ALIGNED(64))a07)[ 0] = b00.i[ 7];
    ((int * ALIGNED(64))a07)[ 1] = b01.i[ 7];
    ((int * ALIGNED(64))a07)[ 2] = b02.i[ 7];
    ((int * ALIGNED(64))a07)[ 3] = b03.i[ 7];
    ((int * ALIGNED(64))a07)[ 4] = b04.i[ 7];
    ((int * ALIGNED(64))a07)[ 5] = b05.i[ 7];
    ((int * ALIGNED(64))a07)[ 6] = b06.i[ 7];
    ((int * ALIGNED(64))a07)[ 7] = b07.i[ 7];
    ((int * ALIGNED(64))a07)[ 8] = b08.i[ 7];
    ((int * ALIGNED(64))a07)[ 9] = b09.i[ 7];
    ((int * ALIGNED(64))a07)[10] = b10.i[ 7];
    ((int * ALIGNED(64))a07)[11] = b11.i[ 7];
    ((int * ALIGNED(64))a07)[12] = b12.i[ 7];
    ((int * ALIGNED(64))a07)[13] = b13.i[ 7];
    ((int * ALIGNED(64))a07)[14] = b14.i[ 7];
    ((int * ALIGNED(64))a07)[15] = b15.i[ 7];

    ((int * ALIGNED(64))a08)[ 0] = b00.i[ 8];
    ((int * ALIGNED(64))a08)[ 1] = b01.i[ 8];
    ((int * ALIGNED(64))a08)[ 2] = b02.i[ 8];
    ((int * ALIGNED(64))a08)[ 3] = b03.i[ 8];
    ((int * ALIGNED(64))a08)[ 4] = b04.i[ 8];
    ((int * ALIGNED(64))a08)[ 5] = b05.i[ 8];
    ((int * ALIGNED(64))a08)[ 6] = b06.i[ 8];
    ((int * ALIGNED(64))a08)[ 7] = b07.i[ 8];
    ((int * ALIGNED(64))a08)[ 8] = b08.i[ 8];
    ((int * ALIGNED(64))a08)[ 9] = b09.i[ 8];
    ((int * ALIGNED(64))a08)[10] = b10.i[ 8];
    ((int * ALIGNED(64))a08)[11] = b11.i[ 8];
    ((int * ALIGNED(64))a08)[12] = b12.i[ 8];
    ((int * ALIGNED(64))a08)[13] = b13.i[ 8];
    ((int * ALIGNED(64))a08)[14] = b14.i[ 8];
    ((int * ALIGNED(64))a08)[15] = b15.i[ 8];

    ((int * ALIGNED(64))a09)[ 0] = b00.i[ 9];
    ((int * ALIGNED(64))a09)[ 1] = b01.i[ 9];
    ((int * ALIGNED(64))a09)[ 2] = b02.i[ 9];
    ((int * ALIGNED(64))a09)[ 3] = b03.i[ 9];
    ((int * ALIGNED(64))a09)[ 4] = b04.i[ 9];
    ((int * ALIGNED(64))a09)[ 5] = b05.i[ 9];
    ((int * ALIGNED(64))a09)[ 6] = b06.i[ 9];
    ((int * ALIGNED(64))a09)[ 7] = b07.i[ 9];
    ((int * ALIGNED(64))a09)[ 8] = b08.i[ 9];
    ((int * ALIGNED(64))a09)[ 9] = b09.i[ 9];
    ((int * ALIGNED(64))a09)[10] = b10.i[ 9];
    ((int * ALIGNED(64))a09)[11] = b11.i[ 9];
    ((int * ALIGNED(64))a09)[12] = b12.i[ 9];
    ((int * ALIGNED(64))a09)[13] = b13.i[ 9];
    ((int * ALIGNED(64))a09)[14] = b14.i[ 9];
    ((int * ALIGNED(64))a09)[15] = b15.i[ 9];

    ((int * ALIGNED(64))a10)[ 0] = b00.i[10];
    ((int * ALIGNED(64))a10)[ 1] = b01.i[10];
    ((int * ALIGNED(64))a10)[ 2] = b02.i[10];
    ((int * ALIGNED(64))a10)[ 3] = b03.i[10];
    ((int * ALIGNED(64))a10)[ 4] = b04.i[10];
    ((int * ALIGNED(64))a10)[ 5] = b05.i[10];
    ((int * ALIGNED(64))a10)[ 6] = b06.i[10];
    ((int * ALIGNED(64))a10)[ 7] = b07.i[10];
    ((int * ALIGNED(64))a10)[ 8] = b08.i[10];
    ((int * ALIGNED(64))a10)[ 9] = b09.i[10];
    ((int * ALIGNED(64))a10)[10] = b10.i[10];
    ((int * ALIGNED(64))a10)[11] = b11.i[10];
    ((int * ALIGNED(64))a10)[12] = b12.i[10];
    ((int * ALIGNED(64))a10)[13] = b13.i[10];
    ((int * ALIGNED(64))a10)[14] = b14.i[10];
    ((int * ALIGNED(64))a10)[15] = b15.i[10];

    ((int * ALIGNED(64))a11)[ 0] = b00.i[11];
    ((int * ALIGNED(64))a11)[ 1] = b01.i[11];
    ((int * ALIGNED(64))a11)[ 2] = b02.i[11];
    ((int * ALIGNED(64))a11)[ 3] = b03.i[11];
    ((int * ALIGNED(64))a11)[ 4] = b04.i[11];
    ((int * ALIGNED(64))a11)[ 5] = b05.i[11];
    ((int * ALIGNED(64))a11)[ 6] = b06.i[11];
    ((int * ALIGNED(64))a11)[ 7] = b07.i[11];
    ((int * ALIGNED(64))a11)[ 8] = b08.i[11];
    ((int * ALIGNED(64))a11)[ 9] = b09.i[11];
    ((int * ALIGNED(64))a11)[10] = b10.i[11];
    ((int * ALIGNED(64))a11)[11] = b11.i[11];
    ((int * ALIGNED(64))a11)[12] = b12.i[11];
    ((int * ALIGNED(64))a11)[13] = b13.i[11];
    ((int * ALIGNED(64))a11)[14] = b14.i[11];
    ((int * ALIGNED(64))a11)[15] = b15.i[11];

    ((int * ALIGNED(64))a12)[ 0] = b00.i[12];
    ((int * ALIGNED(64))a12)[ 1] = b01.i[12];
    ((int * ALIGNED(64))a12)[ 2] = b02.i[12];
    ((int * ALIGNED(64))a12)[ 3] = b03.i[12];
    ((int * ALIGNED(64))a12)[ 4] = b04.i[12];
    ((int * ALIGNED(64))a12)[ 5] = b05.i[12];
    ((int * ALIGNED(64))a12)[ 6] = b06.i[12];
    ((int * ALIGNED(64))a12)[ 7] = b07.i[12];
    ((int * ALIGNED(64))a12)[ 8] = b08.i[12];
    ((int * ALIGNED(64))a12)[ 9] = b09.i[12];
    ((int * ALIGNED(64))a12)[10] = b10.i[12];
    ((int * ALIGNED(64))a12)[11] = b11.i[12];
    ((int * ALIGNED(64))a12)[12] = b12.i[12];
    ((int * ALIGNED(64))a12)[13] = b13.i[12];
    ((int * ALIGNED(64))a12)[14] = b14.i[12];
    ((int * ALIGNED(64))a12)[15] = b15.i[12];

    ((int * ALIGNED(64))a13)[ 0] = b00.i[13];
    ((int * ALIGNED(64))a13)[ 1] = b01.i[13];
    ((int * ALIGNED(64))a13)[ 2] = b02.i[13];
    ((int * ALIGNED(64))a13)[ 3] = b03.i[13];
    ((int * ALIGNED(64))a13)[ 4] = b04.i[13];
    ((int * ALIGNED(64))a13)[ 5] = b05.i[13];
    ((int * ALIGNED(64))a13)[ 6] = b06.i[13];
    ((int * ALIGNED(64))a13)[ 7] = b07.i[13];
    ((int * ALIGNED(64))a13)[ 8] = b08.i[13];
    ((int * ALIGNED(64))a13)[ 9] = b09.i[13];
    ((int * ALIGNED(64))a13)[10] = b10.i[13];
    ((int * ALIGNED(64))a13)[11] = b11.i[13];
    ((int * ALIGNED(64))a13)[12] = b12.i[13];
    ((int * ALIGNED(64))a13)[13] = b13.i[13];
    ((int * ALIGNED(64))a13)[14] = b14.i[13];
    ((int * ALIGNED(64))a13)[15] = b15.i[13];

    ((int * ALIGNED(64))a14)[ 0] = b00.i[14];
    ((int * ALIGNED(64))a14)[ 1] = b01.i[14];
    ((int * ALIGNED(64))a14)[ 2] = b02.i[14];
    ((int * ALIGNED(64))a14)[ 3] = b03.i[14];
    ((int * ALIGNED(64))a14)[ 4] = b04.i[14];
    ((int * ALIGNED(64))a14)[ 5] = b05.i[14];
    ((int * ALIGNED(64))a14)[ 6] = b06.i[14];
    ((int * ALIGNED(64))a14)[ 7] = b07.i[14];
    ((int * ALIGNED(64))a14)[ 8] = b08.i[14];
    ((int * ALIGNED(64))a14)[ 9] = b09.i[14];
    ((int * ALIGNED(64))a14)[10] = b10.i[14];
    ((int * ALIGNED(64))a14)[11] = b11.i[14];
    ((int * ALIGNED(64))a14)[12] = b12.i[14];
    ((int * ALIGNED(64))a14)[13] = b13.i[14];
    ((int * ALIGNED(64))a14)[14] = b14.i[14];
    ((int * ALIGNED(64))a14)[15] = b15.i[14];

    ((int * ALIGNED(64))a15)[ 0] = b00.i[15];
    ((int * ALIGNED(64))a15)[ 1] = b01.i[15];
    ((int * ALIGNED(64))a15)[ 2] = b02.i[15];
    ((int * ALIGNED(64))a15)[ 3] = b03.i[15];
    ((int * ALIGNED(64))a15)[ 4] = b04.i[15];
    ((int * ALIGNED(64))a15)[ 5] = b05.i[15];
    ((int * ALIGNED(64))a15)[ 6] = b06.i[15];
    ((int * ALIGNED(64))a15)[ 7] = b07.i[15];
    ((int * ALIGNED(64))a15)[ 8] = b08.i[15];
    ((int * ALIGNED(64))a15)[ 9] = b09.i[15];
    ((int * ALIGNED(64))a15)[10] = b10.i[15];
    ((int * ALIGNED(64))a15)[11] = b11.i[15];
    ((int * ALIGNED(64))a15)[12] = b12.i[15];
    ((int * ALIGNED(64))a15)[13] = b13.i[15];
    ((int * ALIGNED(64))a15)[14] = b14.i[15];
    ((int * ALIGNED(64))a15)[15] = b15.i[15];
  }

  inline void store_16x16_tr_p( const v16 &b00, const v16 &b01, const v16 &b02, const v16 &b03,
				const v16 &b04, const v16 &b05, const v16 &b06, const v16 &b07,
				const v16 &b08, const v16 &b09, const v16 &b10, const v16 &b11,
				const v16 &b12, const v16 &b13, const v16 &b14, const v16 &b15,
				void * ALIGNED(64) a00, void * ALIGNED(64) a01,
				void * ALIGNED(64) a02, void * ALIGNED(64) a03,
				void * ALIGNED(64) a04, void * ALIGNED(64) a05,
				void * ALIGNED(64) a06, void * ALIGNED(64) a07,
				void * ALIGNED(64) a08, void * ALIGNED(64) a09,
				void * ALIGNED(64) a10, void * ALIGNED(64) a11,
				void * ALIGNED(64) a12, void * ALIGNED(64) a13,
				void * ALIGNED(64) a14, void * ALIGNED(64) a15 )
  {
    ((int * ALIGNED(64))a00)[ 0] = b00.i[ 0];
    ((int * ALIGNED(64))a00)[ 1] = b01.i[ 0];
    ((int * ALIGNED(64))a00)[ 2] = b02.i[ 0];
    ((int * ALIGNED(64))a00)[ 3] = b03.i[ 0];
    ((int * ALIGNED(64))a00)[ 4] = b04.i[ 0];
    ((int * ALIGNED(64))a00)[ 5] = b05.i[ 0];
    ((int * ALIGNED(64))a00)[ 6] = b06.i[ 0];
    ((int * ALIGNED(64))a00)[ 7] = b07.i[ 0];
    ((int * ALIGNED(64))a00)[ 8] = b00.i[ 1];
    ((int * ALIGNED(64))a00)[ 9] = b01.i[ 1];
    ((int * ALIGNED(64))a00)[10] = b02.i[ 1];
    ((int * ALIGNED(64))a00)[11] = b03.i[ 1];
    ((int * ALIGNED(64))a00)[12] = b04.i[ 1];
    ((int * ALIGNED(64))a00)[13] = b05.i[ 1];
    ((int * ALIGNED(64))a00)[14] = b06.i[ 1];
    ((int * ALIGNED(64))a00)[15] = b07.i[ 1];

    ((int * ALIGNED(64))a01)[ 0] = b00.i[ 2];
    ((int * ALIGNED(64))a01)[ 1] = b01.i[ 2];
    ((int * ALIGNED(64))a01)[ 2] = b02.i[ 2];
    ((int * ALIGNED(64))a01)[ 3] = b03.i[ 2];
    ((int * ALIGNED(64))a01)[ 4] = b04.i[ 2];
    ((int * ALIGNED(64))a01)[ 5] = b05.i[ 2];
    ((int * ALIGNED(64))a01)[ 6] = b06.i[ 2];
    ((int * ALIGNED(64))a01)[ 7] = b07.i[ 2];
    ((int * ALIGNED(64))a01)[ 8] = b00.i[ 3];
    ((int * ALIGNED(64))a01)[ 9] = b01.i[ 3];
    ((int * ALIGNED(64))a01)[10] = b02.i[ 3];
    ((int * ALIGNED(64))a01)[11] = b03.i[ 3];
    ((int * ALIGNED(64))a01)[12] = b04.i[ 3];
    ((int * ALIGNED(64))a01)[13] = b05.i[ 3];
    ((int * ALIGNED(64))a01)[14] = b06.i[ 3];
    ((int * ALIGNED(64))a01)[15] = b07.i[ 3];

    ((int * ALIGNED(64))a02)[ 0] = b00.i[ 4];
    ((int * ALIGNED(64))a02)[ 1] = b01.i[ 4];
    ((int * ALIGNED(64))a02)[ 2] = b02.i[ 4];
    ((int * ALIGNED(64))a02)[ 3] = b03.i[ 4];
    ((int * ALIGNED(64))a02)[ 4] = b04.i[ 4];
    ((int * ALIGNED(64))a02)[ 5] = b05.i[ 4];
    ((int * ALIGNED(64))a02)[ 6] = b06.i[ 4];
    ((int * ALIGNED(64))a02)[ 7] = b07.i[ 4];
    ((int * ALIGNED(64))a02)[ 8] = b00.i[ 5];
    ((int * ALIGNED(64))a02)[ 9] = b01.i[ 5];
    ((int * ALIGNED(64))a02)[10] = b02.i[ 5];
    ((int * ALIGNED(64))a02)[11] = b03.i[ 5];
    ((int * ALIGNED(64))a02)[12] = b04.i[ 5];
    ((int * ALIGNED(64))a02)[13] = b05.i[ 5];
    ((int * ALIGNED(64))a02)[14] = b06.i[ 5];
    ((int * ALIGNED(64))a02)[15] = b07.i[ 5];

    ((int * ALIGNED(64))a03)[ 0] = b00.i[ 6];
    ((int * ALIGNED(64))a03)[ 1] = b01.i[ 6];
    ((int * ALIGNED(64))a03)[ 2] = b02.i[ 6];
    ((int * ALIGNED(64))a03)[ 3] = b03.i[ 6];
    ((int * ALIGNED(64))a03)[ 4] = b04.i[ 6];
    ((int * ALIGNED(64))a03)[ 5] = b05.i[ 6];
    ((int * ALIGNED(64))a03)[ 6] = b06.i[ 6];
    ((int * ALIGNED(64))a03)[ 7] = b07.i[ 6];
    ((int * ALIGNED(64))a03)[ 8] = b00.i[ 7];
    ((int * ALIGNED(64))a03)[ 9] = b01.i[ 7];
    ((int * ALIGNED(64))a03)[10] = b02.i[ 7];
    ((int * ALIGNED(64))a03)[11] = b03.i[ 7];
    ((int * ALIGNED(64))a03)[12] = b04.i[ 7];
    ((int * ALIGNED(64))a03)[13] = b05.i[ 7];
    ((int * ALIGNED(64))a03)[14] = b06.i[ 7];
    ((int * ALIGNED(64))a03)[15] = b07.i[ 7];

    ((int * ALIGNED(64))a04)[ 0] = b00.i[ 8];
    ((int * ALIGNED(64))a04)[ 1] = b01.i[ 8];
    ((int * ALIGNED(64))a04)[ 2] = b02.i[ 8];
    ((int * ALIGNED(64))a04)[ 3] = b03.i[ 8];
    ((int * ALIGNED(64))a04)[ 4] = b04.i[ 8];
    ((int * ALIGNED(64))a04)[ 5] = b05.i[ 8];
    ((int * ALIGNED(64))a04)[ 6] = b06.i[ 8];
    ((int * ALIGNED(64))a04)[ 7] = b07.i[ 8];
    ((int * ALIGNED(64))a04)[ 8] = b00.i[ 9];
    ((int * ALIGNED(64))a04)[ 9] = b01.i[ 9];
    ((int * ALIGNED(64))a04)[10] = b02.i[ 9];
    ((int * ALIGNED(64))a04)[11] = b03.i[ 9];
    ((int * ALIGNED(64))a04)[12] = b04.i[ 9];
    ((int * ALIGNED(64))a04)[13] = b05.i[ 9];
    ((int * ALIGNED(64))a04)[14] = b06.i[ 9];
    ((int * ALIGNED(64))a04)[15] = b07.i[ 9];

    ((int * ALIGNED(64))a05)[ 0] = b00.i[10];
    ((int * ALIGNED(64))a05)[ 1] = b01.i[10];
    ((int * ALIGNED(64))a05)[ 2] = b02.i[10];
    ((int * ALIGNED(64))a05)[ 3] = b03.i[10];
    ((int * ALIGNED(64))a05)[ 4] = b04.i[10];
    ((int * ALIGNED(64))a05)[ 5] = b05.i[10];
    ((int * ALIGNED(64))a05)[ 6] = b06.i[10];
    ((int * ALIGNED(64))a05)[ 7] = b07.i[10];
    ((int * ALIGNED(64))a05)[ 8] = b00.i[11];
    ((int * ALIGNED(64))a05)[ 9] = b01.i[11];
    ((int * ALIGNED(64))a05)[10] = b02.i[11];
    ((int * ALIGNED(64))a05)[11] = b03.i[11];
    ((int * ALIGNED(64))a05)[12] = b04.i[11];
    ((int * ALIGNED(64))a05)[13] = b05.i[11];
    ((int * ALIGNED(64))a05)[14] = b06.i[11];
    ((int * ALIGNED(64))a05)[15] = b07.i[11];

    ((int * ALIGNED(64))a06)[ 0] = b00.i[12];
    ((int * ALIGNED(64))a06)[ 1] = b01.i[12];
    ((int * ALIGNED(64))a06)[ 2] = b02.i[12];
    ((int * ALIGNED(64))a06)[ 3] = b03.i[12];
    ((int * ALIGNED(64))a06)[ 4] = b04.i[12];
    ((int * ALIGNED(64))a06)[ 5] = b05.i[12];
    ((int * ALIGNED(64))a06)[ 6] = b06.i[12];
    ((int * ALIGNED(64))a06)[ 7] = b07.i[12];
    ((int * ALIGNED(64))a06)[ 8] = b00.i[13];
    ((int * ALIGNED(64))a06)[ 9] = b01.i[13];
    ((int * ALIGNED(64))a06)[10] = b02.i[13];
    ((int * ALIGNED(64))a06)[11] = b03.i[13];
    ((int * ALIGNED(64))a06)[12] = b04.i[13];
    ((int * ALIGNED(64))a06)[13] = b05.i[13];
    ((int * ALIGNED(64))a06)[14] = b06.i[13];
    ((int * ALIGNED(64))a06)[15] = b07.i[13];

    ((int * ALIGNED(64))a07)[ 0] = b00.i[14];
    ((int * ALIGNED(64))a07)[ 1] = b01.i[14];
    ((int * ALIGNED(64))a07)[ 2] = b02.i[14];
    ((int * ALIGNED(64))a07)[ 3] = b03.i[14];
    ((int * ALIGNED(64))a07)[ 4] = b04.i[14];
    ((int * ALIGNED(64))a07)[ 5] = b05.i[14];
    ((int * ALIGNED(64))a07)[ 6] = b06.i[14];
    ((int * ALIGNED(64))a07)[ 7] = b07.i[14];
    ((int * ALIGNED(64))a07)[ 8] = b00.i[15];
    ((int * ALIGNED(64))a07)[ 9] = b01.i[15];
    ((int * ALIGNED(64))a07)[10] = b02.i[15];
    ((int * ALIGNED(64))a07)[11] = b03.i[15];
    ((int * ALIGNED(64))a07)[12] = b04.i[15];
    ((int * ALIGNED(64))a07)[13] = b05.i[15];
    ((int * ALIGNED(64))a07)[14] = b06.i[15];
    ((int * ALIGNED(64))a07)[15] = b07.i[15];

    ((int * ALIGNED(64))a08)[ 0] = b08.i[ 0];
    ((int * ALIGNED(64))a08)[ 1] = b09.i[ 0];
    ((int * ALIGNED(64))a08)[ 2] = b10.i[ 0];
    ((int * ALIGNED(64))a08)[ 3] = b11.i[ 0];
    ((int * ALIGNED(64))a08)[ 4] = b12.i[ 0];
    ((int * ALIGNED(64))a08)[ 5] = b13.i[ 0];
    ((int * ALIGNED(64))a08)[ 6] = b14.i[ 0];
    ((int * ALIGNED(64))a08)[ 7] = b15.i[ 0];
    ((int * ALIGNED(64))a08)[ 8] = b08.i[ 1];
    ((int * ALIGNED(64))a08)[ 9] = b09.i[ 1];
    ((int * ALIGNED(64))a08)[10] = b10.i[ 1];
    ((int * ALIGNED(64))a08)[11] = b11.i[ 1];
    ((int * ALIGNED(64))a08)[12] = b12.i[ 1];
    ((int * ALIGNED(64))a08)[13] = b13.i[ 1];
    ((int * ALIGNED(64))a08)[14] = b14.i[ 1];
    ((int * ALIGNED(64))a08)[15] = b15.i[ 1];

    ((int * ALIGNED(64))a09)[ 0] = b08.i[ 2];
    ((int * ALIGNED(64))a09)[ 1] = b09.i[ 2];
    ((int * ALIGNED(64))a09)[ 2] = b10.i[ 2];
    ((int * ALIGNED(64))a09)[ 3] = b11.i[ 2];
    ((int * ALIGNED(64))a09)[ 4] = b12.i[ 2];
    ((int * ALIGNED(64))a09)[ 5] = b13.i[ 2];
    ((int * ALIGNED(64))a09)[ 6] = b14.i[ 2];
    ((int * ALIGNED(64))a09)[ 7] = b15.i[ 2];
    ((int * ALIGNED(64))a09)[ 8] = b08.i[ 3];
    ((int * ALIGNED(64))a09)[ 9] = b09.i[ 3];
    ((int * ALIGNED(64))a09)[10] = b10.i[ 3];
    ((int * ALIGNED(64))a09)[11] = b11.i[ 3];
    ((int * ALIGNED(64))a09)[12] = b12.i[ 3];
    ((int * ALIGNED(64))a09)[13] = b13.i[ 3];
    ((int * ALIGNED(64))a09)[14] = b14.i[ 3];
    ((int * ALIGNED(64))a09)[15] = b15.i[ 3];

    ((int * ALIGNED(64))a10)[ 0] = b08.i[ 4];
    ((int * ALIGNED(64))a10)[ 1] = b09.i[ 4];
    ((int * ALIGNED(64))a10)[ 2] = b10.i[ 4];
    ((int * ALIGNED(64))a10)[ 3] = b11.i[ 4];
    ((int * ALIGNED(64))a10)[ 4] = b12.i[ 4];
    ((int * ALIGNED(64))a10)[ 5] = b13.i[ 4];
    ((int * ALIGNED(64))a10)[ 6] = b14.i[ 4];
    ((int * ALIGNED(64))a10)[ 7] = b15.i[ 4];
    ((int * ALIGNED(64))a10)[ 8] = b08.i[ 5];
    ((int * ALIGNED(64))a10)[ 9] = b09.i[ 5];
    ((int * ALIGNED(64))a10)[10] = b10.i[ 5];
    ((int * ALIGNED(64))a10)[11] = b11.i[ 5];
    ((int * ALIGNED(64))a10)[12] = b12.i[ 5];
    ((int * ALIGNED(64))a10)[13] = b13.i[ 5];
    ((int * ALIGNED(64))a10)[14] = b14.i[ 5];
    ((int * ALIGNED(64))a10)[15] = b15.i[ 5];

    ((int * ALIGNED(64))a11)[ 0] = b08.i[ 6];
    ((int * ALIGNED(64))a11)[ 1] = b09.i[ 6];
    ((int * ALIGNED(64))a11)[ 2] = b10.i[ 6];
    ((int * ALIGNED(64))a11)[ 3] = b11.i[ 6];
    ((int * ALIGNED(64))a11)[ 4] = b12.i[ 6];
    ((int * ALIGNED(64))a11)[ 5] = b13.i[ 6];
    ((int * ALIGNED(64))a11)[ 6] = b14.i[ 6];
    ((int * ALIGNED(64))a11)[ 7] = b15.i[ 6];
    ((int * ALIGNED(64))a11)[ 8] = b08.i[ 7];
    ((int * ALIGNED(64))a11)[ 9] = b09.i[ 7];
    ((int * ALIGNED(64))a11)[10] = b10.i[ 7];
    ((int * ALIGNED(64))a11)[11] = b11.i[ 7];
    ((int * ALIGNED(64))a11)[12] = b12.i[ 7];
    ((int * ALIGNED(64))a11)[13] = b13.i[ 7];
    ((int * ALIGNED(64))a11)[14] = b14.i[ 7];
    ((int * ALIGNED(64))a11)[15] = b15.i[ 7];

    ((int * ALIGNED(64))a12)[ 0] = b08.i[ 8];
    ((int * ALIGNED(64))a12)[ 1] = b09.i[ 8];
    ((int * ALIGNED(64))a12)[ 2] = b10.i[ 8];
    ((int * ALIGNED(64))a12)[ 3] = b11.i[ 8];
    ((int * ALIGNED(64))a12)[ 4] = b12.i[ 8];
    ((int * ALIGNED(64))a12)[ 5] = b13.i[ 8];
    ((int * ALIGNED(64))a12)[ 6] = b14.i[ 8];
    ((int * ALIGNED(64))a12)[ 7] = b15.i[ 8];
    ((int * ALIGNED(64))a12)[ 8] = b08.i[ 9];
    ((int * ALIGNED(64))a12)[ 9] = b09.i[ 9];
    ((int * ALIGNED(64))a12)[10] = b10.i[ 9];
    ((int * ALIGNED(64))a12)[11] = b11.i[ 9];
    ((int * ALIGNED(64))a12)[12] = b12.i[ 9];
    ((int * ALIGNED(64))a12)[13] = b13.i[ 9];
    ((int * ALIGNED(64))a12)[14] = b14.i[ 9];
    ((int * ALIGNED(64))a12)[15] = b15.i[ 9];

    ((int * ALIGNED(64))a13)[ 0] = b08.i[10];
    ((int * ALIGNED(64))a13)[ 1] = b09.i[10];
    ((int * ALIGNED(64))a13)[ 2] = b10.i[10];
    ((int * ALIGNED(64))a13)[ 3] = b11.i[10];
    ((int * ALIGNED(64))a13)[ 4] = b12.i[10];
    ((int * ALIGNED(64))a13)[ 5] = b13.i[10];
    ((int * ALIGNED(64))a13)[ 6] = b14.i[10];
    ((int * ALIGNED(64))a13)[ 7] = b15.i[10];
    ((int * ALIGNED(64))a13)[ 8] = b08.i[11];
    ((int * ALIGNED(64))a13)[ 9] = b09.i[11];
    ((int * ALIGNED(64))a13)[10] = b10.i[11];
    ((int * ALIGNED(64))a13)[11] = b11.i[11];
    ((int * ALIGNED(64))a13)[12] = b12.i[11];
    ((int * ALIGNED(64))a13)[13] = b13.i[11];
    ((int * ALIGNED(64))a13)[14] = b14.i[11];
    ((int * ALIGNED(64))a13)[15] = b15.i[11];

    ((int * ALIGNED(64))a14)[ 0] = b08.i[12];
    ((int * ALIGNED(64))a14)[ 1] = b09.i[12];
    ((int * ALIGNED(64))a14)[ 2] = b10.i[12];
    ((int * ALIGNED(64))a14)[ 3] = b11.i[12];
    ((int * ALIGNED(64))a14)[ 4] = b12.i[12];
    ((int * ALIGNED(64))a14)[ 5] = b13.i[12];
    ((int * ALIGNED(64))a14)[ 6] = b14.i[12];
    ((int * ALIGNED(64))a14)[ 7] = b15.i[12];
    ((int * ALIGNED(64))a14)[ 8] = b08.i[13];
    ((int * ALIGNED(64))a14)[ 9] = b09.i[13];
    ((int * ALIGNED(64))a14)[10] = b10.i[13];
    ((int * ALIGNED(64))a14)[11] = b11.i[13];
    ((int * ALIGNED(64))a14)[12] = b12.i[13];
    ((int * ALIGNED(64))a14)[13] = b13.i[13];
    ((int * ALIGNED(64))a14)[14] = b14.i[13];
    ((int * ALIGNED(64))a14)[15] = b15.i[13];

    ((int * ALIGNED(64))a15)[ 0] = b08.i[14];
    ((int * ALIGNED(64))a15)[ 1] = b09.i[14];
    ((int * ALIGNED(64))a15)[ 2] = b10.i[14];
    ((int * ALIGNED(64))a15)[ 3] = b11.i[14];
    ((int * ALIGNED(64))a15)[ 4] = b12.i[14];
    ((int * ALIGNED(64))a15)[ 5] = b13.i[14];
    ((int * ALIGNED(64))a15)[ 6] = b14.i[14];
    ((int * ALIGNED(64))a15)[ 7] = b15.i[14];
    ((int * ALIGNED(64))a15)[ 8] = b08.i[15];
    ((int * ALIGNED(64))a15)[ 9] = b09.i[15];
    ((int * ALIGNED(64))a15)[10] = b10.i[15];
    ((int * ALIGNED(64))a15)[11] = b11.i[15];
    ((int * ALIGNED(64))a15)[12] = b12.i[15];
    ((int * ALIGNED(64))a15)[13] = b13.i[15];
    ((int * ALIGNED(64))a15)[14] = b14.i[15];
    ((int * ALIGNED(64))a15)[15] = b15.i[15];
  }

  //////////////
  // v16int class

  class v16int : public v16
  {
    // v16int prefix unary operator friends

    friend inline v16int operator  +( const v16int & a );
    friend inline v16int operator  -( const v16int & a );
    friend inline v16int operator  ~( const v16int & a );
    friend inline v16int operator  !( const v16int & a );
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector

    // v16int prefix increment / decrement operator friends

    friend inline v16int operator ++( v16int & a );
    friend inline v16int operator --( v16int & a );

    // v16int postfix increment / decrement operator friends

    friend inline v16int operator ++( v16int & a, int );
    friend inline v16int operator --( v16int & a, int );

    // v16int binary operator friends

    friend inline v16int operator  +( const v16int &a, const v16int &b );
    friend inline v16int operator  -( const v16int &a, const v16int &b );
    friend inline v16int operator  *( const v16int &a, const v16int &b );
    friend inline v16int operator  /( const v16int &a, const v16int &b );
    friend inline v16int operator  %( const v16int &a, const v16int &b );
    friend inline v16int operator  ^( const v16int &a, const v16int &b );
    friend inline v16int operator  &( const v16int &a, const v16int &b );
    friend inline v16int operator  |( const v16int &a, const v16int &b );
    friend inline v16int operator <<( const v16int &a, const v16int &b );
    friend inline v16int operator >>( const v16int &a, const v16int &b );

    // v16int logical operator friends

    friend inline v16int operator  <( const v16int &a, const v16int &b );
    friend inline v16int operator  >( const v16int &a, const v16int &b );
    friend inline v16int operator ==( const v16int &a, const v16int &b );
    friend inline v16int operator !=( const v16int &a, const v16int &b );
    friend inline v16int operator <=( const v16int &a, const v16int &b );
    friend inline v16int operator >=( const v16int &a, const v16int &b );
    friend inline v16int operator &&( const v16int &a, const v16int &b );
    friend inline v16int operator ||( const v16int &a, const v16int &b );

    // v16int miscellaneous friends

    friend inline v16int abs( const v16int &a );
    friend inline v16    czero( const v16int &c, const v16 &a );
    friend inline v16 notczero( const v16int &c, const v16 &a );
    // FIXME: cswap, notcswap!
    friend inline v16 merge( const v16int &c, const v16 &t, const v16 &f );

    // v16float unary operator friends

    friend inline v16int operator  !( const v16float & a );

    // v16float logical operator friends

    friend inline v16int operator  <( const v16float &a, const v16float &b );
    friend inline v16int operator  >( const v16float &a, const v16float &b );
    friend inline v16int operator ==( const v16float &a, const v16float &b );
    friend inline v16int operator !=( const v16float &a, const v16float &b );
    friend inline v16int operator <=( const v16float &a, const v16float &b );
    friend inline v16int operator >=( const v16float &a, const v16float &b );
    friend inline v16int operator &&( const v16float &a, const v16float &b );
    friend inline v16int operator ||( const v16float &a, const v16float &b );

    // v16float miscellaneous friends

    friend inline v16float clear_bits(  const v16int &m, const v16float &a );
    friend inline v16float set_bits(    const v16int &m, const v16float &a );
    friend inline v16float toggle_bits( const v16int &m, const v16float &a );

  public:

    // v16int constructors / destructors

    v16int() {}                                  // Default constructor

    v16int( const v16int &a )                    // Copy constructor
    {
      i[ 0] = a.i[ 0]; i[ 1] = a.i[ 1]; i[ 2] = a.i[ 2]; i[ 3] = a.i[ 3];
      i[ 4] = a.i[ 4]; i[ 5] = a.i[ 5]; i[ 6] = a.i[ 6]; i[ 7] = a.i[ 7];
      i[ 8] = a.i[ 8]; i[ 9] = a.i[ 9]; i[10] = a.i[10]; i[11] = a.i[11];
      i[12] = a.i[12]; i[13] = a.i[13]; i[14] = a.i[14]; i[15] = a.i[15];
    }

    v16int( const v16 &a )                       // Init from mixed
    {
      i[ 0] = a.i[ 0]; i[ 1] = a.i[ 1]; i[ 2] = a.i[ 2]; i[ 3] = a.i[ 3];
      i[ 4] = a.i[ 4]; i[ 5] = a.i[ 5]; i[ 6] = a.i[ 6]; i[ 7] = a.i[ 7];
      i[ 8] = a.i[ 8]; i[ 9] = a.i[ 9]; i[10] = a.i[10]; i[11] = a.i[11];
      i[12] = a.i[12]; i[13] = a.i[13]; i[14] = a.i[14]; i[15] = a.i[15];
    }

    v16int( int a )                              // Init from scalar
    {
      i[ 0] = a; i[ 1] = a; i[ 2] = a; i[ 3] = a;
      i[ 4] = a; i[ 5] = a; i[ 6] = a; i[ 7] = a;
      i[ 8] = a; i[ 9] = a; i[10] = a; i[11] = a;
      i[12] = a; i[13] = a; i[14] = a; i[15] = a;
    }

    v16int( int i00, int i01, int i02, int i03,
	    int i04, int i05, int i06, int i07,
	    int i08, int i09, int i10, int i11,
	    int i12, int i13, int i14, int i15 ) // Init from scalars
    {
      i[ 0] = i00; i[ 1] = i01; i[ 2] = i02; i[ 3] = i03;
      i[ 4] = i04; i[ 5] = i05; i[ 6] = i06; i[ 7] = i07;
      i[ 8] = i08; i[ 9] = i09; i[10] = i10; i[11] = i11;
      i[12] = i12; i[13] = i13; i[14] = i14; i[15] = i15;
    }

    ~v16int() {}                                 // Destructor

    // v16int assignment operators

#   define ASSIGN(op)			          \
    inline v16int &operator op( const v16int &b ) \
    {						  \
      i[ 0] op b.i[ 0];                           \
      i[ 1] op b.i[ 1];                           \
      i[ 2] op b.i[ 2];                           \
      i[ 3] op b.i[ 3];                           \
      i[ 4] op b.i[ 4];                           \
      i[ 5] op b.i[ 5];                           \
      i[ 6] op b.i[ 6];                           \
      i[ 7] op b.i[ 7];                           \
      i[ 8] op b.i[ 8];                           \
      i[ 9] op b.i[ 9];                           \
      i[10] op b.i[10];                           \
      i[11] op b.i[11];                           \
      i[12] op b.i[12];                           \
      i[13] op b.i[13];                           \
      i[14] op b.i[14];                           \
      i[15] op b.i[15];                           \
      return *this;                               \
    }

    ASSIGN( =)
    ASSIGN(+=)
    ASSIGN(-=)
    ASSIGN(*=)
    ASSIGN(/=)
    ASSIGN(%=)
    ASSIGN(^=)
    ASSIGN(&=)
    ASSIGN(|=)
    ASSIGN(<<=)
    ASSIGN(>>=)

#   undef ASSIGN

    // v16int member access operator

    inline int &operator []( int n )
    {
      return i[n];
    }

    inline int  operator ()( int n )
    {
      return i[n];
    }
  };

  // v16int prefix unary operators

# define PREFIX_UNARY(op)                       \
  inline v16int operator op( const v16int & a ) \
  {						\
    v16int b;                                   \
    b.i[ 0] = (op a.i[ 0]);                     \
    b.i[ 1] = (op a.i[ 1]);                     \
    b.i[ 2] = (op a.i[ 2]);                     \
    b.i[ 3] = (op a.i[ 3]);                     \
    b.i[ 4] = (op a.i[ 4]);                     \
    b.i[ 5] = (op a.i[ 5]);                     \
    b.i[ 6] = (op a.i[ 6]);                     \
    b.i[ 7] = (op a.i[ 7]);                     \
    b.i[ 8] = (op a.i[ 8]);                     \
    b.i[ 9] = (op a.i[ 9]);                     \
    b.i[10] = (op a.i[10]);                     \
    b.i[11] = (op a.i[11]);                     \
    b.i[12] = (op a.i[12]);                     \
    b.i[13] = (op a.i[13]);                     \
    b.i[14] = (op a.i[14]);                     \
    b.i[15] = (op a.i[15]);                     \
    return b;                                   \
  }

  PREFIX_UNARY(+)
  PREFIX_UNARY(-)

  inline v16int operator !( const v16int & a )
  {
    v16int b;
    b.i[ 0] = - ( !a.i[ 0] );
    b.i[ 1] = - ( !a.i[ 1] );
    b.i[ 2] = - ( !a.i[ 2] );
    b.i[ 3] = - ( !a.i[ 3] );
    b.i[ 4] = - ( !a.i[ 4] );
    b.i[ 5] = - ( !a.i[ 5] );
    b.i[ 6] = - ( !a.i[ 6] );
    b.i[ 7] = - ( !a.i[ 7] );
    b.i[ 8] = - ( !a.i[ 8] );
    b.i[ 9] = - ( !a.i[ 9] );
    b.i[10] = - ( !a.i[10] );
    b.i[11] = - ( !a.i[11] );
    b.i[12] = - ( !a.i[12] );
    b.i[13] = - ( !a.i[13] );
    b.i[14] = - ( !a.i[14] );
    b.i[15] = - ( !a.i[15] );
    return b;
  }

  PREFIX_UNARY(~)

# undef PREFIX_UNARY

  // v16int prefix increment / decrement

# define PREFIX_INCDEC(op)                      \
  inline v16int operator op( v16int & a )       \
  {						\
    v16int b;                                   \
    b.i[ 0] = ( op a.i[ 0] );                   \
    b.i[ 1] = ( op a.i[ 1] );                   \
    b.i[ 2] = ( op a.i[ 2] );                   \
    b.i[ 3] = ( op a.i[ 3] );                   \
    b.i[ 4] = ( op a.i[ 4] );                   \
    b.i[ 5] = ( op a.i[ 5] );                   \
    b.i[ 6] = ( op a.i[ 6] );                   \
    b.i[ 7] = ( op a.i[ 7] );                   \
    b.i[ 8] = ( op a.i[ 8] );                   \
    b.i[ 9] = ( op a.i[ 9] );                   \
    b.i[10] = ( op a.i[10] );                   \
    b.i[11] = ( op a.i[11] );                   \
    b.i[12] = ( op a.i[12] );                   \
    b.i[13] = ( op a.i[13] );                   \
    b.i[14] = ( op a.i[14] );                   \
    b.i[15] = ( op a.i[15] );                   \
    return b;                                   \
  }

  PREFIX_INCDEC(++)
  PREFIX_INCDEC(--)

# undef PREFIX_INCDEC

  // v16int postfix increment / decrement

# define POSTFIX_INCDEC(op)                    \
  inline v16int operator op( v16int & a, int ) \
  {					       \
    v16int b;                                  \
    b.i[ 0] = ( a.i[ 0] op );                  \
    b.i[ 1] = ( a.i[ 1] op );                  \
    b.i[ 2] = ( a.i[ 2] op );                  \
    b.i[ 3] = ( a.i[ 3] op );                  \
    b.i[ 4] = ( a.i[ 4] op );                  \
    b.i[ 5] = ( a.i[ 5] op );                  \
    b.i[ 6] = ( a.i[ 6] op );                  \
    b.i[ 7] = ( a.i[ 7] op );                  \
    b.i[ 8] = ( a.i[ 8] op );                  \
    b.i[ 9] = ( a.i[ 9] op );                  \
    b.i[10] = ( a.i[10] op );                  \
    b.i[11] = ( a.i[11] op );                  \
    b.i[12] = ( a.i[12] op );                  \
    b.i[13] = ( a.i[13] op );                  \
    b.i[14] = ( a.i[14] op );                  \
    b.i[15] = ( a.i[15] op );                  \
    return b;                                  \
  }

  POSTFIX_INCDEC(++)
  POSTFIX_INCDEC(--)

# undef POSTFIX_INCDEC

  // v16int binary operators

# define BINARY(op)                                             \
  inline v16int operator op( const v16int &a, const v16int &b ) \
  {								\
    v16int c;                                                   \
    c.i[ 0] = a.i[ 0] op b.i[ 0];                               \
    c.i[ 1] = a.i[ 1] op b.i[ 1];                               \
    c.i[ 2] = a.i[ 2] op b.i[ 2];                               \
    c.i[ 3] = a.i[ 3] op b.i[ 3];                               \
    c.i[ 4] = a.i[ 4] op b.i[ 4];                               \
    c.i[ 5] = a.i[ 5] op b.i[ 5];                               \
    c.i[ 6] = a.i[ 6] op b.i[ 6];                               \
    c.i[ 7] = a.i[ 7] op b.i[ 7];                               \
    c.i[ 8] = a.i[ 8] op b.i[ 8];                               \
    c.i[ 9] = a.i[ 9] op b.i[ 9];                               \
    c.i[10] = a.i[10] op b.i[10];                               \
    c.i[11] = a.i[11] op b.i[11];                               \
    c.i[12] = a.i[12] op b.i[12];                               \
    c.i[13] = a.i[13] op b.i[13];                               \
    c.i[14] = a.i[14] op b.i[14];                               \
    c.i[15] = a.i[15] op b.i[15];                               \
    return c;                                                   \
  }

  BINARY(+)
  BINARY(-)
  BINARY(*)
  BINARY(/)
  BINARY(%)
  BINARY(^)
  BINARY(&)
  BINARY(|)
  BINARY(<<)
  BINARY(>>)

# undef BINARY

  // v16int logical operators

# define LOGICAL(op)                                            \
  inline v16int operator op( const v16int &a, const v16int &b ) \
  {                                                             \
    v16int c;                                                   \
    c.i[ 0] = - ( a.i[ 0] op b.i[ 0] );                         \
    c.i[ 1] = - ( a.i[ 1] op b.i[ 1] );                         \
    c.i[ 2] = - ( a.i[ 2] op b.i[ 2] );                         \
    c.i[ 3] = - ( a.i[ 3] op b.i[ 3] );                         \
    c.i[ 4] = - ( a.i[ 4] op b.i[ 4] );                         \
    c.i[ 5] = - ( a.i[ 5] op b.i[ 5] );                         \
    c.i[ 6] = - ( a.i[ 6] op b.i[ 6] );                         \
    c.i[ 7] = - ( a.i[ 7] op b.i[ 7] );                         \
    c.i[ 8] = - ( a.i[ 8] op b.i[ 8] );                         \
    c.i[ 9] = - ( a.i[ 9] op b.i[ 9] );                         \
    c.i[10] = - ( a.i[10] op b.i[10] );                         \
    c.i[11] = - ( a.i[11] op b.i[11] );                         \
    c.i[12] = - ( a.i[12] op b.i[12] );                         \
    c.i[13] = - ( a.i[13] op b.i[13] );                         \
    c.i[14] = - ( a.i[14] op b.i[14] );                         \
    c.i[15] = - ( a.i[15] op b.i[15] );                         \
    return c;                                                   \
  }

  LOGICAL(<)
  LOGICAL(>)
  LOGICAL(==)
  LOGICAL(!=)
  LOGICAL(<=)
  LOGICAL(>=)
  LOGICAL(&&)
  LOGICAL(||)

# undef LOGICAL

  // v16int miscellaneous functions

  inline v16int abs( const v16int &a )
  {
    v16int b;

    b.i[ 0] = ( a.i[ 0] >= 0 ) ? a.i[ 0] : - a.i[ 0];
    b.i[ 1] = ( a.i[ 1] >= 0 ) ? a.i[ 1] : - a.i[ 1];
    b.i[ 2] = ( a.i[ 2] >= 0 ) ? a.i[ 2] : - a.i[ 2];
    b.i[ 3] = ( a.i[ 3] >= 0 ) ? a.i[ 3] : - a.i[ 3];
    b.i[ 4] = ( a.i[ 4] >= 0 ) ? a.i[ 4] : - a.i[ 4];
    b.i[ 5] = ( a.i[ 5] >= 0 ) ? a.i[ 5] : - a.i[ 5];
    b.i[ 6] = ( a.i[ 6] >= 0 ) ? a.i[ 6] : - a.i[ 6];
    b.i[ 7] = ( a.i[ 7] >= 0 ) ? a.i[ 7] : - a.i[ 7];
    b.i[ 8] = ( a.i[ 8] >= 0 ) ? a.i[ 8] : - a.i[ 8];
    b.i[ 9] = ( a.i[ 9] >= 0 ) ? a.i[ 9] : - a.i[ 9];
    b.i[10] = ( a.i[10] >= 0 ) ? a.i[10] : - a.i[10];
    b.i[11] = ( a.i[11] >= 0 ) ? a.i[11] : - a.i[11];
    b.i[12] = ( a.i[12] >= 0 ) ? a.i[12] : - a.i[12];
    b.i[13] = ( a.i[13] >= 0 ) ? a.i[13] : - a.i[13];
    b.i[14] = ( a.i[14] >= 0 ) ? a.i[14] : - a.i[14];
    b.i[15] = ( a.i[15] >= 0 ) ? a.i[15] : - a.i[15];

    return b;
  }

  inline v16 czero( const v16int &c, const v16 &a )
  {
    v16 b;

    b.i[ 0] = a.i[ 0] & ~c.i[ 0];
    b.i[ 1] = a.i[ 1] & ~c.i[ 1];
    b.i[ 2] = a.i[ 2] & ~c.i[ 2];
    b.i[ 3] = a.i[ 3] & ~c.i[ 3];
    b.i[ 4] = a.i[ 4] & ~c.i[ 4];
    b.i[ 5] = a.i[ 5] & ~c.i[ 5];
    b.i[ 6] = a.i[ 6] & ~c.i[ 6];
    b.i[ 7] = a.i[ 7] & ~c.i[ 7];
    b.i[ 8] = a.i[ 8] & ~c.i[ 8];
    b.i[ 9] = a.i[ 9] & ~c.i[ 9];
    b.i[10] = a.i[10] & ~c.i[10];
    b.i[11] = a.i[11] & ~c.i[11];
    b.i[12] = a.i[12] & ~c.i[12];
    b.i[13] = a.i[13] & ~c.i[13];
    b.i[14] = a.i[14] & ~c.i[14];
    b.i[15] = a.i[15] & ~c.i[15];

    return b;
  }

  inline v16 notczero( const v16int &c, const v16 &a )
  {
    v16 b;

    b.i[ 0] = a.i[ 0] & c.i[ 0];
    b.i[ 1] = a.i[ 1] & c.i[ 1];
    b.i[ 2] = a.i[ 2] & c.i[ 2];
    b.i[ 3] = a.i[ 3] & c.i[ 3];
    b.i[ 4] = a.i[ 4] & c.i[ 4];
    b.i[ 5] = a.i[ 5] & c.i[ 5];
    b.i[ 6] = a.i[ 6] & c.i[ 6];
    b.i[ 7] = a.i[ 7] & c.i[ 7];
    b.i[ 8] = a.i[ 8] & c.i[ 8];
    b.i[ 9] = a.i[ 9] & c.i[ 9];
    b.i[10] = a.i[10] & c.i[10];
    b.i[11] = a.i[11] & c.i[11];
    b.i[12] = a.i[12] & c.i[12];
    b.i[13] = a.i[13] & c.i[13];
    b.i[14] = a.i[14] & c.i[14];
    b.i[15] = a.i[15] & c.i[15];

    return b;
  }

  inline v16 merge( const v16int &c, const v16 &t, const v16 &f )
  {
    v16 m;

    m.i[ 0] = ( f.i[ 0] & ~c.i[ 0] ) | ( t.i[ 0] & c.i[ 0] );
    m.i[ 1] = ( f.i[ 1] & ~c.i[ 1] ) | ( t.i[ 1] & c.i[ 1] );
    m.i[ 2] = ( f.i[ 2] & ~c.i[ 2] ) | ( t.i[ 2] & c.i[ 2] );
    m.i[ 3] = ( f.i[ 3] & ~c.i[ 3] ) | ( t.i[ 3] & c.i[ 3] );
    m.i[ 4] = ( f.i[ 4] & ~c.i[ 4] ) | ( t.i[ 4] & c.i[ 4] );
    m.i[ 5] = ( f.i[ 5] & ~c.i[ 5] ) | ( t.i[ 5] & c.i[ 5] );
    m.i[ 6] = ( f.i[ 6] & ~c.i[ 6] ) | ( t.i[ 6] & c.i[ 6] );
    m.i[ 7] = ( f.i[ 7] & ~c.i[ 7] ) | ( t.i[ 7] & c.i[ 7] );
    m.i[ 8] = ( f.i[ 8] & ~c.i[ 8] ) | ( t.i[ 8] & c.i[ 8] );
    m.i[ 9] = ( f.i[ 9] & ~c.i[ 9] ) | ( t.i[ 9] & c.i[ 9] );
    m.i[10] = ( f.i[10] & ~c.i[10] ) | ( t.i[10] & c.i[10] );
    m.i[11] = ( f.i[11] & ~c.i[11] ) | ( t.i[11] & c.i[11] );
    m.i[12] = ( f.i[12] & ~c.i[12] ) | ( t.i[12] & c.i[12] );
    m.i[13] = ( f.i[13] & ~c.i[13] ) | ( t.i[13] & c.i[13] );
    m.i[14] = ( f.i[14] & ~c.i[14] ) | ( t.i[14] & c.i[14] );
    m.i[15] = ( f.i[15] & ~c.i[15] ) | ( t.i[15] & c.i[15] );

    return m;
  }

  ////////////////
  // v16float class

  class v16float : public v16
  {
    // v16float prefix unary operator friends

    friend inline v16float operator  +( const v16float &a );
    friend inline v16float operator  -( const v16float &a );
    friend inline v16float operator  ~( const v16float &a );
    friend inline v16int   operator  !( const v16float &a );
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector

    // v16float prefix increment / decrement operator friends

    friend inline v16float operator ++( v16float &a );
    friend inline v16float operator --( v16float &a );

    // v16float postfix increment / decrement operator friends

    friend inline v16float operator ++( v16float &a, int );
    friend inline v16float operator --( v16float &a, int );

    // v16float binary operator friends

    friend inline v16float operator  +( const v16float &a, const v16float &b );
    friend inline v16float operator  -( const v16float &a, const v16float &b );
    friend inline v16float operator  *( const v16float &a, const v16float &b );
    friend inline v16float operator  /( const v16float &a, const v16float &b );

    // v16float logical operator friends

    friend inline v16int operator  <( const v16float &a, const v16float &b );
    friend inline v16int operator  >( const v16float &a, const v16float &b );
    friend inline v16int operator ==( const v16float &a, const v16float &b );
    friend inline v16int operator !=( const v16float &a, const v16float &b );
    friend inline v16int operator <=( const v16float &a, const v16float &b );
    friend inline v16int operator >=( const v16float &a, const v16float &b );
    friend inline v16int operator &&( const v16float &a, const v16float &b );
    friend inline v16int operator ||( const v16float &a, const v16float &b );

    // v16float math library friends

#   define CMATH_FR1(fn) friend inline v16float fn( const v16float &a )
#   define CMATH_FR2(fn) friend inline v16float fn( const v16float &a,  \
                                                    const v16float &b )

    CMATH_FR1(acos);  CMATH_FR1(asin);  CMATH_FR1(atan); CMATH_FR2(atan2);
    CMATH_FR1(ceil);  CMATH_FR1(cos);   CMATH_FR1(cosh); CMATH_FR1(exp);
    CMATH_FR1(fabs);  CMATH_FR1(floor); CMATH_FR2(fmod); CMATH_FR1(log);
    CMATH_FR1(log10); CMATH_FR2(pow);   CMATH_FR1(sin);  CMATH_FR1(sinh);
    CMATH_FR1(sqrt);  CMATH_FR1(tan);   CMATH_FR1(tanh);

    CMATH_FR2(copysign);

#   undef CMATH_FR1
#   undef CMATH_FR2

    // v16float miscellaneous friends

    friend inline v16float rsqrt_approx( const v16float &a );
    friend inline v16float rsqrt       ( const v16float &a );
    friend inline v16float rcp_approx( const v16float &a );
    friend inline v16float rcp       ( const v16float &a );
    friend inline v16float fma ( const v16float &a, const v16float &b, const v16float &c );
    friend inline v16float fms ( const v16float &a, const v16float &b, const v16float &c );
    friend inline v16float fnms( const v16float &a, const v16float &b, const v16float &c );
    friend inline v16float  clear_bits( const v16int &m, const v16float &a );
    friend inline v16float    set_bits( const v16int &m, const v16float &a );
    friend inline v16float toggle_bits( const v16int &m, const v16float &a );
    friend inline void increment_16x1( float * ALIGNED(64) p, const v16float &a );
    friend inline void decrement_16x1( float * ALIGNED(64) p, const v16float &a );
    friend inline void     scale_16x1( float * ALIGNED(64) p, const v16float &a );

  public:

    // v16float constructors / destructors

    v16float() {}                                          // Default constructor

    v16float( const v16float &a )                          // Copy constructor
    {
      f[ 0] = a.f[ 0]; f[ 1] = a.f[ 1]; f[ 2] = a.f[ 2]; f[ 3] = a.f[ 3];
      f[ 4] = a.f[ 4]; f[ 5] = a.f[ 5]; f[ 6] = a.f[ 6]; f[ 7] = a.f[ 7];
      f[ 8] = a.f[ 8]; f[ 9] = a.f[ 9]; f[10] = a.f[10]; f[11] = a.f[11];
      f[12] = a.f[12]; f[13] = a.f[13]; f[14] = a.f[14]; f[15] = a.f[15];
    }

    v16float( const v16 &a )                               // Init from mixed
    {
      f[ 0] = a.f[ 0]; f[ 1] = a.f[ 1]; f[ 2] = a.f[ 2]; f[ 3] = a.f[ 3];
      f[ 4] = a.f[ 4]; f[ 5] = a.f[ 5]; f[ 6] = a.f[ 6]; f[ 7] = a.f[ 7];
      f[ 8] = a.f[ 8]; f[ 9] = a.f[ 9]; f[10] = a.f[10]; f[11] = a.f[11];
      f[12] = a.f[12]; f[13] = a.f[13]; f[14] = a.f[14]; f[15] = a.f[15];
    }

    v16float( float a )                                    // Init from scalar
    {
      f[ 0] = a; f[ 1] = a; f[ 2] = a; f[ 3] = a;
      f[ 4] = a; f[ 5] = a; f[ 6] = a; f[ 7] = a;
      f[ 8] = a; f[ 9] = a; f[10] = a; f[11] = a;
      f[12] = a; f[13] = a; f[14] = a; f[15] = a;
    }

    v16float( float f00, float f01, float f02, float f03,
	      float f04, float f05, float f06, float f07,
	      float f08, float f09, float f10, float f11,
	      float f12, float f13, float f14, float f15 ) // Init from scalars
    {
      f[ 0] = f00; f[ 1] = f01; f[ 2] = f02; f[ 3] = f03;
      f[ 4] = f04; f[ 5] = f05; f[ 6] = f06; f[ 7] = f07;
      f[ 8] = f08; f[ 9] = f09; f[10] = f10; f[11] = f11;
      f[12] = f12; f[13] = f13; f[14] = f14; f[15] = f15;
    }

    ~v16float() {}                                         // Destructor

    // v16float assignment operators

#   define ASSIGN(op)                                   \
    inline v16float &operator op( const v16float &b )   \
    {							\
      f[ 0] op b.f[ 0];		             		\
      f[ 1] op b.f[ 1];                                 \
      f[ 2] op b.f[ 2];                                 \
      f[ 3] op b.f[ 3];                                 \
      f[ 4] op b.f[ 4];                                 \
      f[ 5] op b.f[ 5];                                 \
      f[ 6] op b.f[ 6];                                 \
      f[ 7] op b.f[ 7];                                 \
      f[ 8] op b.f[ 8];                                 \
      f[ 9] op b.f[ 9];                                 \
      f[10] op b.f[10];                                 \
      f[11] op b.f[11];                                 \
      f[12] op b.f[12];                                 \
      f[13] op b.f[13];                                 \
      f[14] op b.f[14];                                 \
      f[15] op b.f[15];                                 \
      return *this;                                     \
    }

    ASSIGN(=)
    ASSIGN(+=)
    ASSIGN(-=)
    ASSIGN(*=)
    ASSIGN(/=)

#   undef ASSIGN

    // v16float member access operator

    inline float &operator []( int n )
    {
      return f[n];
    }

    inline float  operator ()( int n )
    {
      return f[n];
    }
  };

  // v16float prefix unary operators

  inline v16float operator +( const v16float &a )
  {
    v16float b;

    b.f[ 0] = +a.f[ 0];
    b.f[ 1] = +a.f[ 1];
    b.f[ 2] = +a.f[ 2];
    b.f[ 3] = +a.f[ 3];
    b.f[ 4] = +a.f[ 4];
    b.f[ 5] = +a.f[ 5];
    b.f[ 6] = +a.f[ 6];
    b.f[ 7] = +a.f[ 7];
    b.f[ 8] = +a.f[ 8];
    b.f[ 9] = +a.f[ 9];
    b.f[10] = +a.f[10];
    b.f[11] = +a.f[11];
    b.f[12] = +a.f[12];
    b.f[13] = +a.f[13];
    b.f[14] = +a.f[14];
    b.f[15] = +a.f[15];

    return b;
  }

  inline v16float operator -( const v16float &a )
  {
    v16float b;

    b.f[ 0] = -a.f[ 0];
    b.f[ 1] = -a.f[ 1];
    b.f[ 2] = -a.f[ 2];
    b.f[ 3] = -a.f[ 3];
    b.f[ 4] = -a.f[ 4];
    b.f[ 5] = -a.f[ 5];
    b.f[ 6] = -a.f[ 6];
    b.f[ 7] = -a.f[ 7];
    b.f[ 8] = -a.f[ 8];
    b.f[ 9] = -a.f[ 9];
    b.f[10] = -a.f[10];
    b.f[11] = -a.f[11];
    b.f[12] = -a.f[12];
    b.f[13] = -a.f[13];
    b.f[14] = -a.f[14];
    b.f[15] = -a.f[15];

    return b;
  }

  inline v16int operator !( const v16float &a )
  {
    v16int b;

    b.i[ 0] = a.i[ 0] ? 0 : -1;
    b.i[ 1] = a.i[ 1] ? 0 : -1;
    b.i[ 2] = a.i[ 2] ? 0 : -1;
    b.i[ 3] = a.i[ 3] ? 0 : -1;
    b.i[ 4] = a.i[ 4] ? 0 : -1;
    b.i[ 5] = a.i[ 5] ? 0 : -1;
    b.i[ 6] = a.i[ 6] ? 0 : -1;
    b.i[ 7] = a.i[ 7] ? 0 : -1;
    b.i[ 8] = a.i[ 8] ? 0 : -1;
    b.i[ 9] = a.i[ 9] ? 0 : -1;
    b.i[10] = a.i[10] ? 0 : -1;
    b.i[11] = a.i[11] ? 0 : -1;
    b.i[12] = a.i[12] ? 0 : -1;
    b.i[13] = a.i[13] ? 0 : -1;
    b.i[14] = a.i[14] ? 0 : -1;
    b.i[15] = a.i[15] ? 0 : -1;

    return b;
  }

  // v16float prefix increment / decrement operators

  inline v16float operator ++( v16float &a )
  {
    v16float b;

    b.f[ 0] = ++a.f[ 0];
    b.f[ 1] = ++a.f[ 1];
    b.f[ 2] = ++a.f[ 2];
    b.f[ 3] = ++a.f[ 3];
    b.f[ 4] = ++a.f[ 4];
    b.f[ 5] = ++a.f[ 5];
    b.f[ 6] = ++a.f[ 6];
    b.f[ 7] = ++a.f[ 7];
    b.f[ 8] = ++a.f[ 8];
    b.f[ 9] = ++a.f[ 9];
    b.f[10] = ++a.f[10];
    b.f[11] = ++a.f[11];
    b.f[12] = ++a.f[12];
    b.f[13] = ++a.f[13];
    b.f[14] = ++a.f[14];
    b.f[15] = ++a.f[15];

    return b;
  }

  inline v16float operator --( v16float &a )
  {
    v16float b;

    b.f[ 0] = --a.f[ 0];
    b.f[ 1] = --a.f[ 1];
    b.f[ 2] = --a.f[ 2];
    b.f[ 3] = --a.f[ 3];
    b.f[ 4] = --a.f[ 4];
    b.f[ 5] = --a.f[ 5];
    b.f[ 6] = --a.f[ 6];
    b.f[ 7] = --a.f[ 7];
    b.f[ 8] = --a.f[ 8];
    b.f[ 9] = --a.f[ 9];
    b.f[10] = --a.f[10];
    b.f[11] = --a.f[11];
    b.f[12] = --a.f[12];
    b.f[13] = --a.f[13];
    b.f[14] = --a.f[14];
    b.f[15] = --a.f[15];

    return b;
  }

  // v16float postfix increment / decrement operators

  inline v16float operator ++( v16float &a, int )
  {
    v16float b;

    b.f[ 0] = a.f[ 0]++;
    b.f[ 1] = a.f[ 1]++;
    b.f[ 2] = a.f[ 2]++;
    b.f[ 3] = a.f[ 3]++;
    b.f[ 4] = a.f[ 4]++;
    b.f[ 5] = a.f[ 5]++;
    b.f[ 6] = a.f[ 6]++;
    b.f[ 7] = a.f[ 7]++;
    b.f[ 8] = a.f[ 8]++;
    b.f[ 9] = a.f[ 9]++;
    b.f[10] = a.f[10]++;
    b.f[11] = a.f[11]++;
    b.f[12] = a.f[12]++;
    b.f[13] = a.f[13]++;
    b.f[14] = a.f[14]++;
    b.f[15] = a.f[15]++;

    return b;
  }

  inline v16float operator --( v16float &a, int )
  {
    v16float b;

    b.f[ 0] = a.f[ 0]--;
    b.f[ 1] = a.f[ 1]--;
    b.f[ 2] = a.f[ 2]--;
    b.f[ 3] = a.f[ 3]--;
    b.f[ 4] = a.f[ 4]--;
    b.f[ 5] = a.f[ 5]--;
    b.f[ 6] = a.f[ 6]--;
    b.f[ 7] = a.f[ 7]--;
    b.f[ 8] = a.f[ 8]--;
    b.f[ 9] = a.f[ 9]--;
    b.f[10] = a.f[10]--;
    b.f[11] = a.f[11]--;
    b.f[12] = a.f[12]--;
    b.f[13] = a.f[13]--;
    b.f[14] = a.f[14]--;
    b.f[15] = a.f[15]--;

    return b;
  }

  // v16float binary operators

# define BINARY(op)                                                   \
  inline v16float operator op( const v16float &a, const v16float &b ) \
  {								      \
    v16float c;                                                       \
    c.f[ 0] = a.f[ 0] op b.f[ 0];                                     \
    c.f[ 1] = a.f[ 1] op b.f[ 1];                                     \
    c.f[ 2] = a.f[ 2] op b.f[ 2];                                     \
    c.f[ 3] = a.f[ 3] op b.f[ 3];                                     \
    c.f[ 4] = a.f[ 4] op b.f[ 4];                                     \
    c.f[ 5] = a.f[ 5] op b.f[ 5];                                     \
    c.f[ 6] = a.f[ 6] op b.f[ 6];                                     \
    c.f[ 7] = a.f[ 7] op b.f[ 7];                                     \
    c.f[ 8] = a.f[ 8] op b.f[ 8];                                     \
    c.f[ 9] = a.f[ 9] op b.f[ 9];                                     \
    c.f[10] = a.f[10] op b.f[10];                                     \
    c.f[11] = a.f[11] op b.f[11];                                     \
    c.f[12] = a.f[12] op b.f[12];                                     \
    c.f[13] = a.f[13] op b.f[13];                                     \
    c.f[14] = a.f[14] op b.f[14];                                     \
    c.f[15] = a.f[15] op b.f[15];                                     \
    return c;                                                         \
  }

  BINARY(+)
  BINARY(-)
  BINARY(*)
  BINARY(/)

# undef BINARY

  // v16float logical operators

# define LOGICAL(op)                                                \
  inline v16int operator op( const v16float &a, const v16float &b ) \
  {								    \
    v16int c;                                                       \
    c.i[ 0] = -( a.f[ 0] op b.f[ 0] );                              \
    c.i[ 1] = -( a.f[ 1] op b.f[ 1] );                              \
    c.i[ 2] = -( a.f[ 2] op b.f[ 2] );                              \
    c.i[ 3] = -( a.f[ 3] op b.f[ 3] );                              \
    c.i[ 4] = -( a.f[ 4] op b.f[ 4] );                              \
    c.i[ 5] = -( a.f[ 5] op b.f[ 5] );                              \
    c.i[ 6] = -( a.f[ 6] op b.f[ 6] );                              \
    c.i[ 7] = -( a.f[ 7] op b.f[ 7] );                              \
    c.i[ 8] = -( a.f[ 8] op b.f[ 8] );                              \
    c.i[ 9] = -( a.f[ 9] op b.f[ 9] );                              \
    c.i[10] = -( a.f[10] op b.f[10] );                              \
    c.i[11] = -( a.f[11] op b.f[11] );                              \
    c.i[12] = -( a.f[12] op b.f[12] );                              \
    c.i[13] = -( a.f[13] op b.f[13] );                              \
    c.i[14] = -( a.f[14] op b.f[14] );                              \
    c.i[15] = -( a.f[15] op b.f[15] );                              \
    return c;                                                       \
  }

  LOGICAL(< )
  LOGICAL(> )
  LOGICAL(==)
  LOGICAL(!=)
  LOGICAL(<=)
  LOGICAL(>=)
  LOGICAL(&&)
  LOGICAL(||)

# undef LOGICAL

  // v16float math library functions

# define CMATH_FR1(fn)                          \
  inline v16float fn( const v16float &a )       \
  {						\
    v16float b;                                 \
    b.f[ 0] = ::fn( a.f[ 0] );                  \
    b.f[ 1] = ::fn( a.f[ 1] );                  \
    b.f[ 2] = ::fn( a.f[ 2] );                  \
    b.f[ 3] = ::fn( a.f[ 3] );                  \
    b.f[ 4] = ::fn( a.f[ 4] );                  \
    b.f[ 5] = ::fn( a.f[ 5] );                  \
    b.f[ 6] = ::fn( a.f[ 6] );                  \
    b.f[ 7] = ::fn( a.f[ 7] );                  \
    b.f[ 8] = ::fn( a.f[ 8] );                  \
    b.f[ 9] = ::fn( a.f[ 9] );                  \
    b.f[10] = ::fn( a.f[10] );                  \
    b.f[11] = ::fn( a.f[11] );                  \
    b.f[12] = ::fn( a.f[12] );                  \
    b.f[13] = ::fn( a.f[13] );                  \
    b.f[14] = ::fn( a.f[14] );                  \
    b.f[15] = ::fn( a.f[15] );                  \
    return b;                                   \
  }

# define CMATH_FR2(fn)                                          \
  inline v16float fn( const v16float &a, const v16float &b )    \
  {								\
    v16float c;                                                 \
    c.f[ 0] = ::fn( a.f[ 0], b.f[ 0] );                         \
    c.f[ 1] = ::fn( a.f[ 1], b.f[ 1] );                         \
    c.f[ 2] = ::fn( a.f[ 2], b.f[ 2] );                         \
    c.f[ 3] = ::fn( a.f[ 3], b.f[ 3] );                         \
    c.f[ 4] = ::fn( a.f[ 4], b.f[ 4] );                         \
    c.f[ 5] = ::fn( a.f[ 5], b.f[ 5] );                         \
    c.f[ 6] = ::fn( a.f[ 6], b.f[ 6] );                         \
    c.f[ 7] = ::fn( a.f[ 7], b.f[ 7] );                         \
    c.f[ 8] = ::fn( a.f[ 8], b.f[ 8] );                         \
    c.f[ 9] = ::fn( a.f[ 9], b.f[ 9] );                         \
    c.f[10] = ::fn( a.f[10], b.f[10] );                         \
    c.f[11] = ::fn( a.f[11], b.f[11] );                         \
    c.f[12] = ::fn( a.f[12], b.f[12] );                         \
    c.f[13] = ::fn( a.f[13], b.f[13] );                         \
    c.f[14] = ::fn( a.f[14], b.f[14] );                         \
    c.f[15] = ::fn( a.f[15], b.f[15] );                         \
    return c;                                                   \
  }

  CMATH_FR1(acos)     CMATH_FR1(asin)  CMATH_FR1(atan) CMATH_FR2(atan2)
  CMATH_FR1(ceil)     CMATH_FR1(cos)   CMATH_FR1(cosh) CMATH_FR1(exp)
  CMATH_FR1(fabs)     CMATH_FR1(floor) CMATH_FR2(fmod) CMATH_FR1(log)
  CMATH_FR1(log10)    CMATH_FR2(pow)   CMATH_FR1(sin)  CMATH_FR1(sinh)
  CMATH_FR1(sqrt)     CMATH_FR1(tan)   CMATH_FR1(tanh)

  inline v16float copysign( const v16float &a, const v16float &b )
  {
    v16float c;
    float t;

    t = ::fabs( a.f[ 0] );
    if( b.f[ 0] < 0 ) t = -t;
    c.f[ 0] = t;

    t = ::fabs( a.f[ 1] );
    if( b.f[ 1] < 0 ) t = -t;
    c.f[ 1] = t;

    t = ::fabs( a.f[ 2] );
    if( b.f[ 2] < 0 ) t = -t;
    c.f[ 2] = t;

    t = ::fabs( a.f[ 3] );
    if( b.f[ 3] < 0 ) t = -t;
    c.f[ 3] = t;

    t = ::fabs( a.f[ 4] );
    if( b.f[ 4] < 0 ) t = -t;
    c.f[ 4] = t;

    t = ::fabs( a.f[ 5] );
    if( b.f[ 5] < 0 ) t = -t;
    c.f[ 5] = t;

    t = ::fabs( a.f[ 6] );
    if( b.f[ 6] < 0 ) t = -t;
    c.f[ 6] = t;

    t = ::fabs( a.f[ 7] );
    if( b.f[ 7] < 0 ) t = -t;
    c.f[ 7] = t;

    t = ::fabs( a.f[ 8] );
    if( b.f[ 8] < 0 ) t = -t;
    c.f[ 8] = t;

    t = ::fabs( a.f[ 9] );
    if( b.f[ 9] < 0 ) t = -t;
    c.f[ 9] = t;

    t = ::fabs( a.f[10] );
    if( b.f[10] < 0 ) t = -t;
    c.f[10] = t;

    t = ::fabs( a.f[11] );
    if( b.f[11] < 0 ) t = -t;
    c.f[11] = t;

    t = ::fabs( a.f[12] );
    if( b.f[12] < 0 ) t = -t;
    c.f[12] = t;

    t = ::fabs( a.f[13] );
    if( b.f[13] < 0 ) t = -t;
    c.f[13] = t;

    t = ::fabs( a.f[14] );
    if( b.f[14] < 0 ) t = -t;
    c.f[14] = t;

    t = ::fabs( a.f[15] );
    if( b.f[15] < 0 ) t = -t;
    c.f[15] = t;

    return c;
  }

# undef CMATH_FR1
# undef CMATH_FR2

  // v16float miscellaneous functions

  inline v16float rsqrt_approx( const v16float &a )
  {
    v16float b;

    b.f[ 0] = ::sqrt( 1.0f/a.f[ 0] );
    b.f[ 1] = ::sqrt( 1.0f/a.f[ 1] );
    b.f[ 2] = ::sqrt( 1.0f/a.f[ 2] );
    b.f[ 3] = ::sqrt( 1.0f/a.f[ 3] );
    b.f[ 4] = ::sqrt( 1.0f/a.f[ 4] );
    b.f[ 5] = ::sqrt( 1.0f/a.f[ 5] );
    b.f[ 6] = ::sqrt( 1.0f/a.f[ 6] );
    b.f[ 7] = ::sqrt( 1.0f/a.f[ 7] );
    b.f[ 8] = ::sqrt( 1.0f/a.f[ 8] );
    b.f[ 9] = ::sqrt( 1.0f/a.f[ 9] );
    b.f[10] = ::sqrt( 1.0f/a.f[10] );
    b.f[11] = ::sqrt( 1.0f/a.f[11] );
    b.f[12] = ::sqrt( 1.0f/a.f[12] );
    b.f[13] = ::sqrt( 1.0f/a.f[13] );
    b.f[14] = ::sqrt( 1.0f/a.f[14] );
    b.f[15] = ::sqrt( 1.0f/a.f[15] );

    return b;
  }

  inline v16float rsqrt( const v16float &a )
  {
    v16float b;

    b.f[ 0] = ::sqrt( 1.0f/a.f[ 0] );
    b.f[ 1] = ::sqrt( 1.0f/a.f[ 1] );
    b.f[ 2] = ::sqrt( 1.0f/a.f[ 2] );
    b.f[ 3] = ::sqrt( 1.0f/a.f[ 3] );
    b.f[ 4] = ::sqrt( 1.0f/a.f[ 4] );
    b.f[ 5] = ::sqrt( 1.0f/a.f[ 5] );
    b.f[ 6] = ::sqrt( 1.0f/a.f[ 6] );
    b.f[ 7] = ::sqrt( 1.0f/a.f[ 7] );
    b.f[ 8] = ::sqrt( 1.0f/a.f[ 8] );
    b.f[ 9] = ::sqrt( 1.0f/a.f[ 9] );
    b.f[10] = ::sqrt( 1.0f/a.f[10] );
    b.f[11] = ::sqrt( 1.0f/a.f[11] );
    b.f[12] = ::sqrt( 1.0f/a.f[12] );
    b.f[13] = ::sqrt( 1.0f/a.f[13] );
    b.f[14] = ::sqrt( 1.0f/a.f[14] );
    b.f[15] = ::sqrt( 1.0f/a.f[15] );

    return b;
  }

  inline v16float rcp_approx( const v16float &a )
  {
    v16float b;

    b.f[ 0] = 1.0f/a.f[ 0];
    b.f[ 1] = 1.0f/a.f[ 1];
    b.f[ 2] = 1.0f/a.f[ 2];
    b.f[ 3] = 1.0f/a.f[ 3];
    b.f[ 4] = 1.0f/a.f[ 4];
    b.f[ 5] = 1.0f/a.f[ 5];
    b.f[ 6] = 1.0f/a.f[ 6];
    b.f[ 7] = 1.0f/a.f[ 7];
    b.f[ 8] = 1.0f/a.f[ 8];
    b.f[ 9] = 1.0f/a.f[ 9];
    b.f[10] = 1.0f/a.f[10];
    b.f[11] = 1.0f/a.f[11];
    b.f[12] = 1.0f/a.f[12];
    b.f[13] = 1.0f/a.f[13];
    b.f[14] = 1.0f/a.f[14];
    b.f[15] = 1.0f/a.f[15];

    return b;
  }

  inline v16float rcp( const v16float &a )
  {
    v16float b;

    b.f[ 0] = 1.0f/a.f[ 0];
    b.f[ 1] = 1.0f/a.f[ 1];
    b.f[ 2] = 1.0f/a.f[ 2];
    b.f[ 3] = 1.0f/a.f[ 3];
    b.f[ 4] = 1.0f/a.f[ 4];
    b.f[ 5] = 1.0f/a.f[ 5];
    b.f[ 6] = 1.0f/a.f[ 6];
    b.f[ 7] = 1.0f/a.f[ 7];
    b.f[ 8] = 1.0f/a.f[ 8];
    b.f[ 9] = 1.0f/a.f[ 9];
    b.f[10] = 1.0f/a.f[10];
    b.f[11] = 1.0f/a.f[11];
    b.f[12] = 1.0f/a.f[12];
    b.f[13] = 1.0f/a.f[13];
    b.f[14] = 1.0f/a.f[14];
    b.f[15] = 1.0f/a.f[15];

    return b;
  }

  inline v16float fma( const v16float &a, const v16float &b, const v16float &c )
  {
    v16float d;

    d.f[ 0] = a.f[ 0] * b.f[ 0] + c.f[ 0];
    d.f[ 1] = a.f[ 1] * b.f[ 1] + c.f[ 1];
    d.f[ 2] = a.f[ 2] * b.f[ 2] + c.f[ 2];
    d.f[ 3] = a.f[ 3] * b.f[ 3] + c.f[ 3];
    d.f[ 4] = a.f[ 4] * b.f[ 4] + c.f[ 4];
    d.f[ 5] = a.f[ 5] * b.f[ 5] + c.f[ 5];
    d.f[ 6] = a.f[ 6] * b.f[ 6] + c.f[ 6];
    d.f[ 7] = a.f[ 7] * b.f[ 7] + c.f[ 7];
    d.f[ 8] = a.f[ 8] * b.f[ 8] + c.f[ 8];
    d.f[ 9] = a.f[ 9] * b.f[ 9] + c.f[ 9];
    d.f[10] = a.f[10] * b.f[10] + c.f[10];
    d.f[11] = a.f[11] * b.f[11] + c.f[11];
    d.f[12] = a.f[12] * b.f[12] + c.f[12];
    d.f[13] = a.f[13] * b.f[13] + c.f[13];
    d.f[14] = a.f[14] * b.f[14] + c.f[14];
    d.f[15] = a.f[15] * b.f[15] + c.f[15];

    return d;
  }

  inline v16float fms( const v16float &a, const v16float &b, const v16float &c )
  {
    v16float d;

    d.f[0] = a.f[0] * b.f[0] - c.f[0];
    d.f[1] = a.f[1] * b.f[1] - c.f[1];
    d.f[2] = a.f[2] * b.f[2] - c.f[2];
    d.f[3] = a.f[3] * b.f[3] - c.f[3];
    d.f[4] = a.f[4] * b.f[4] - c.f[4];
    d.f[5] = a.f[5] * b.f[5] - c.f[5];
    d.f[6] = a.f[6] * b.f[6] - c.f[6];
    d.f[7] = a.f[7] * b.f[7] - c.f[7];
    d.f[0] = a.f[0] * b.f[0] - c.f[0];
    d.f[1] = a.f[1] * b.f[1] - c.f[1];
    d.f[2] = a.f[2] * b.f[2] - c.f[2];
    d.f[3] = a.f[3] * b.f[3] - c.f[3];
    d.f[4] = a.f[4] * b.f[4] - c.f[4];
    d.f[5] = a.f[5] * b.f[5] - c.f[5];
    d.f[6] = a.f[6] * b.f[6] - c.f[6];
    d.f[7] = a.f[7] * b.f[7] - c.f[7];

    return d;
  }

  inline v16float fnms( const v16float &a, const v16float &b, const v16float &c )
  {
    v16float d;

    d.f[ 0] = c.f[ 0] - a.f[ 0] * b.f[ 0];
    d.f[ 1] = c.f[ 1] - a.f[ 1] * b.f[ 1];
    d.f[ 2] = c.f[ 2] - a.f[ 2] * b.f[ 2];
    d.f[ 3] = c.f[ 3] - a.f[ 3] * b.f[ 3];
    d.f[ 4] = c.f[ 4] - a.f[ 4] * b.f[ 4];
    d.f[ 5] = c.f[ 5] - a.f[ 5] * b.f[ 5];
    d.f[ 6] = c.f[ 6] - a.f[ 6] * b.f[ 6];
    d.f[ 7] = c.f[ 7] - a.f[ 7] * b.f[ 7];
    d.f[ 8] = c.f[ 8] - a.f[ 8] * b.f[ 8];
    d.f[ 9] = c.f[ 9] - a.f[ 9] * b.f[ 9];
    d.f[10] = c.f[10] - a.f[10] * b.f[10];
    d.f[11] = c.f[11] - a.f[11] * b.f[11];
    d.f[12] = c.f[12] - a.f[12] * b.f[12];
    d.f[13] = c.f[13] - a.f[13] * b.f[13];
    d.f[14] = c.f[14] - a.f[14] * b.f[14];
    d.f[15] = c.f[15] - a.f[15] * b.f[15];

    return d;
  }

  inline v16float clear_bits( const v16int &m, const v16float &a )
  {
    v16float b;

    b.i[ 0] = ( ~m.i[ 0] ) & a.i[ 0];
    b.i[ 1] = ( ~m.i[ 1] ) & a.i[ 1];
    b.i[ 2] = ( ~m.i[ 2] ) & a.i[ 2];
    b.i[ 3] = ( ~m.i[ 3] ) & a.i[ 3];
    b.i[ 4] = ( ~m.i[ 4] ) & a.i[ 4];
    b.i[ 5] = ( ~m.i[ 5] ) & a.i[ 5];
    b.i[ 6] = ( ~m.i[ 6] ) & a.i[ 6];
    b.i[ 7] = ( ~m.i[ 7] ) & a.i[ 7];
    b.i[ 8] = ( ~m.i[ 8] ) & a.i[ 8];
    b.i[ 9] = ( ~m.i[ 9] ) & a.i[ 9];
    b.i[10] = ( ~m.i[10] ) & a.i[10];
    b.i[11] = ( ~m.i[11] ) & a.i[11];
    b.i[12] = ( ~m.i[12] ) & a.i[12];
    b.i[13] = ( ~m.i[13] ) & a.i[13];
    b.i[14] = ( ~m.i[14] ) & a.i[14];
    b.i[15] = ( ~m.i[15] ) & a.i[15];

    return b;
  }

  inline v16float set_bits( const v16int &m, const v16float &a )
  {
    v16float b;

    b.i[ 0] = m.i[ 0] | a.i[ 0];
    b.i[ 1] = m.i[ 1] | a.i[ 1];
    b.i[ 2] = m.i[ 2] | a.i[ 2];
    b.i[ 3] = m.i[ 3] | a.i[ 3];
    b.i[ 4] = m.i[ 4] | a.i[ 4];
    b.i[ 5] = m.i[ 5] | a.i[ 5];
    b.i[ 6] = m.i[ 6] | a.i[ 6];
    b.i[ 7] = m.i[ 7] | a.i[ 7];
    b.i[ 8] = m.i[ 8] | a.i[ 8];
    b.i[ 9] = m.i[ 9] | a.i[ 9];
    b.i[10] = m.i[10] | a.i[10];
    b.i[11] = m.i[11] | a.i[11];
    b.i[12] = m.i[12] | a.i[12];
    b.i[13] = m.i[13] | a.i[13];
    b.i[14] = m.i[14] | a.i[14];
    b.i[15] = m.i[15] | a.i[15];

    return b;
  }

  inline v16float toggle_bits( const v16int &m, const v16float &a )
  {
    v16float b;

    b.i[ 0] = m.i[ 0] ^ a.i[ 0];
    b.i[ 1] = m.i[ 1] ^ a.i[ 1];
    b.i[ 2] = m.i[ 2] ^ a.i[ 2];
    b.i[ 3] = m.i[ 3] ^ a.i[ 3];
    b.i[ 4] = m.i[ 4] ^ a.i[ 4];
    b.i[ 5] = m.i[ 5] ^ a.i[ 5];
    b.i[ 6] = m.i[ 6] ^ a.i[ 6];
    b.i[ 7] = m.i[ 7] ^ a.i[ 7];
    b.i[ 8] = m.i[ 8] ^ a.i[ 8];
    b.i[ 9] = m.i[ 9] ^ a.i[ 9];
    b.i[10] = m.i[10] ^ a.i[10];
    b.i[11] = m.i[11] ^ a.i[11];
    b.i[12] = m.i[12] ^ a.i[12];
    b.i[13] = m.i[13] ^ a.i[13];
    b.i[14] = m.i[14] ^ a.i[14];
    b.i[15] = m.i[15] ^ a.i[15];

    return b;
  }

  inline void increment_16x1( float * ALIGNED(64) p, const v16float &a )
  {
    p[ 0] += a.f[ 0];
    p[ 1] += a.f[ 1];
    p[ 2] += a.f[ 2];
    p[ 3] += a.f[ 3];
    p[ 4] += a.f[ 4];
    p[ 5] += a.f[ 5];
    p[ 6] += a.f[ 6];
    p[ 7] += a.f[ 7];
    p[ 8] += a.f[ 8];
    p[ 9] += a.f[ 9];
    p[10] += a.f[10];
    p[11] += a.f[11];
    p[12] += a.f[12];
    p[13] += a.f[13];
    p[14] += a.f[14];
    p[15] += a.f[15];
  }

  inline void decrement_16x1( float * ALIGNED(64) p, const v16float &a )
  {
    p[ 0] -= a.f[ 0];
    p[ 1] -= a.f[ 1];
    p[ 2] -= a.f[ 2];
    p[ 3] -= a.f[ 3];
    p[ 4] -= a.f[ 4];
    p[ 5] -= a.f[ 5];
    p[ 6] -= a.f[ 6];
    p[ 7] -= a.f[ 7];
    p[ 8] -= a.f[ 8];
    p[ 9] -= a.f[ 9];
    p[10] -= a.f[10];
    p[11] -= a.f[11];
    p[12] -= a.f[12];
    p[13] -= a.f[13];
    p[14] -= a.f[14];
    p[15] -= a.f[15];
  }

  inline void scale_16x1( float * ALIGNED(64) p, const v16float &a )
  {
    p[ 0] *= a.f[ 0];
    p[ 1] *= a.f[ 1];
    p[ 2] *= a.f[ 2];
    p[ 3] *= a.f[ 3];
    p[ 4] *= a.f[ 4];
    p[ 5] *= a.f[ 5];
    p[ 6] *= a.f[ 6];
    p[ 7] *= a.f[ 7];
    p[ 8] *= a.f[ 8];
    p[ 9] *= a.f[ 9];
    p[10] *= a.f[10];
    p[11] *= a.f[11];
    p[12] *= a.f[12];
    p[13] *= a.f[13];
    p[14] *= a.f[14];
    p[15] *= a.f[15];
  }

} // namespace v16

#endif // _v16_avx512_h_
