#ifndef _v8_portable_h_
#define _v8_portable_h_

#ifndef IN_v8_h
#error "Do not include v8_portable.h directly; use v8.h"
#endif

#include <math.h>

#define V8_ACCELERATION
#define V8_PORTABLE_ACCELERATION

#ifndef ALIGNED
#define ALIGNED(n)
#endif

#define ALWAYS_INLINE __attribute__((always_inline))

namespace v8
{
  class v8;
  class v8int;
  class v8float;

  ////////////////
  // v8 base class

  class v8
  {
    friend class v8int;
    friend class v8float;

    // v8 miscellaneous friends

    friend inline int any( const v8 &a ) ALWAYS_INLINE;
    friend inline int all( const v8 &a ) ALWAYS_INLINE;

    template<int n>
    friend inline v8 splat( const v8 &a ) ALWAYS_INLINE;

    template<int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7>
    friend inline v8 shuffle( const v8 &a ) ALWAYS_INLINE;

    friend inline void swap( v8 &a, v8 &b ) ALWAYS_INLINE;
    friend inline void transpose( v8 &a0, v8 &a1, v8 &a2, v8 &a3,
				  v8 &a4, v8 &a5, v8 &a6, v8 &a7 ) ALWAYS_INLINE;

    // v8int miscellaneous friends

    friend inline v8    czero( const v8int &c, const v8 &a ) ALWAYS_INLINE;
    friend inline v8 notczero( const v8int &c, const v8 &a ) ALWAYS_INLINE;
    friend inline v8 merge( const v8int &c, const v8 &a, const v8 &b ) ALWAYS_INLINE;

    // v8 memory manipulation friends

    friend inline void   load_8x1( const void * ALIGNED(16) p, v8 &a ) ALWAYS_INLINE;
    friend inline void  store_8x1( const v8 &a, void * ALIGNED(16) p ) ALWAYS_INLINE;
    friend inline void stream_8x1( const v8 &a, void * ALIGNED(16) p ) ALWAYS_INLINE;
    friend inline void  clear_8x1( void * ALIGNED(16) dst ) ALWAYS_INLINE;
    friend inline void   copy_8x1( void * ALIGNED(16) dst,
                                   const void * ALIGNED(16) src ) ALWAYS_INLINE;
    friend inline void   swap_8x1( void * ALIGNED(16) a, void * ALIGNED(16) b ) ALWAYS_INLINE;

    // v8 transposed memory manipulation friends
    // Note: Half aligned values are permissible in the 8x2_tr variants.

    friend inline void load_8x1_tr( const void *a0, const void *a1,
                                    const void *a2, const void *a3,
				    const void *a4, const void *a5,
                                    const void *a6, const void *a7,
                                    v8 &a ) ALWAYS_INLINE;

    friend inline void load_8x2_tr( const void * ALIGNED(8) a0,
                                    const void * ALIGNED(8) a1,
                                    const void * ALIGNED(8) a2,
                                    const void * ALIGNED(8) a3,
				    const void * ALIGNED(8) a4,
                                    const void * ALIGNED(8) a5,
                                    const void * ALIGNED(8) a6,
                                    const void * ALIGNED(8) a7,
                                    v8 &a, v8 &b ) ALWAYS_INLINE;

    friend inline void load_8x3_tr( const void * ALIGNED(16) a0,
                                    const void * ALIGNED(16) a1,
                                    const void * ALIGNED(16) a2,
                                    const void * ALIGNED(16) a3,
				    const void * ALIGNED(16) a4,
                                    const void * ALIGNED(16) a5,
                                    const void * ALIGNED(16) a6,
                                    const void * ALIGNED(16) a7,
                                    v8 &a, v8 &b, v8 &c ) ALWAYS_INLINE;

    friend inline void load_8x4_tr( const void * ALIGNED(16) a0,
                                    const void * ALIGNED(16) a1,
                                    const void * ALIGNED(16) a2,
                                    const void * ALIGNED(16) a3,
				    const void * ALIGNED(16) a4,
                                    const void * ALIGNED(16) a5,
                                    const void * ALIGNED(16) a6,
                                    const void * ALIGNED(16) a7,
                                    v8 &a, v8 &b, v8 &c, v8 &d ) ALWAYS_INLINE;

    friend inline void load_8x8_tr( const void * ALIGNED(16) a0,
                                    const void * ALIGNED(16) a1,
                                    const void * ALIGNED(16) a2,
                                    const void * ALIGNED(16) a3,
				    const void * ALIGNED(16) a4,
                                    const void * ALIGNED(16) a5,
                                    const void * ALIGNED(16) a6,
                                    const void * ALIGNED(16) a7,
                                    v8 &a, v8 &b, v8 &c, v8 &d,
                                    v8 &e, v8 &f, v8 &g, v8 &h ) ALWAYS_INLINE;

    friend inline void store_8x1_tr( const v8 &a,
                                     void *a0, void *a1, void *a2, void *a3,
                                     void *a4, void *a5, void *a6, void *a7 ) ALWAYS_INLINE;

    friend inline void store_8x2_tr( const v8 &a, const v8 &b,
                                     void * ALIGNED(8) a0,
                                     void * ALIGNED(8) a1,
                                     void * ALIGNED(8) a2,
                                     void * ALIGNED(8) a3,
                                     void * ALIGNED(8) a4,
                                     void * ALIGNED(8) a5,
                                     void * ALIGNED(8) a6,
                                     void * ALIGNED(8) a7 ) ALWAYS_INLINE;

    friend inline void store_8x3_tr( const v8 &a, const v8 &b, const v8 &c,
                                     void * ALIGNED(16) a0,
                                     void * ALIGNED(16) a1,
                                     void * ALIGNED(16) a2,
                                     void * ALIGNED(16) a3,
                                     void * ALIGNED(16) a4,
                                     void * ALIGNED(16) a5,
                                     void * ALIGNED(16) a6,
                                     void * ALIGNED(16) a7 ) ALWAYS_INLINE;

    friend inline void store_8x4_tr( const v8 &a, const v8 &b,
                                     const v8 &c, const v8 &d,
                                     void * ALIGNED(16) a0,
                                     void * ALIGNED(16) a1,
                                     void * ALIGNED(16) a2,
                                     void * ALIGNED(16) a3,
                                     void * ALIGNED(16) a4,
                                     void * ALIGNED(16) a5,
                                     void * ALIGNED(16) a6,
                                     void * ALIGNED(16) a7 ) ALWAYS_INLINE;

    friend inline void store_8x8_tr( const v8 &a, const v8 &b,
                                     const v8 &c, const v8 &d,
                                     const v8 &e, const v8 &f,
                                     const v8 &g, const v8 &h,
                                     void * ALIGNED(16) a0,
                                     void * ALIGNED(16) a1,
                                     void * ALIGNED(16) a2,
                                     void * ALIGNED(16) a3,
                                     void * ALIGNED(16) a4,
                                     void * ALIGNED(16) a5,
                                     void * ALIGNED(16) a6,
                                     void * ALIGNED(16) a7 ) ALWAYS_INLINE;

  protected:

    union
    {
      int i[8];
      float f[8];
    };

  public:

    v8() {}                    // Default constructor

    v8( const v8 &a )          // Copy constructor
    {
      i[0]=a.i[0]; i[1]=a.i[1]; i[2]=a.i[2]; i[3]=a.i[3];
      i[4]=a.i[4]; i[5]=a.i[5]; i[6]=a.i[6]; i[7]=a.i[7];
    }

    ~v8() {}                   // Default destructor
  };

  // v8 miscellaneous functions

  inline int any( const v8 &a )
  {
    return a.i[0] || a.i[1] || a.i[2] || a.i[3] ||
           a.i[4] || a.i[5] || a.i[6] || a.i[7];
  }

  inline int all( const v8 &a )
  {
    return a.i[0] && a.i[1] && a.i[2] && a.i[3] &&
           a.i[4] && a.i[5] && a.i[6] && a.i[7];
  }

  template<int n>
  inline v8 splat( const v8 & a )
  {
    v8 b;

    b.i[0] = a.i[n];
    b.i[1] = a.i[n];
    b.i[2] = a.i[n];
    b.i[3] = a.i[n];
    b.i[4] = a.i[n];
    b.i[5] = a.i[n];
    b.i[6] = a.i[n];
    b.i[7] = a.i[n];

    return b;
  }

  template<int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7>
  inline v8 shuffle( const v8 & a )
  {
    v8 b;

    b.i[0] = a.i[i0];
    b.i[1] = a.i[i1];
    b.i[2] = a.i[i2];
    b.i[3] = a.i[i3];
    b.i[4] = a.i[i4];
    b.i[5] = a.i[i5];
    b.i[6] = a.i[i6];
    b.i[7] = a.i[i7];

    return b;
  }

# define sw(x,y) x^=y, y^=x, x^=y

  inline void swap( v8 &a, v8 &b )
  {
    sw( a.i[0], b.i[0] );
    sw( a.i[1], b.i[1] );
    sw( a.i[2], b.i[2] );
    sw( a.i[3], b.i[3] );
    sw( a.i[4], b.i[4] );
    sw( a.i[5], b.i[5] );
    sw( a.i[6], b.i[6] );
    sw( a.i[7], b.i[7] );
  }

  inline void transpose( v8 &a0, v8 &a1, v8 &a2, v8 &a3,
			 v8 &a4, v8 &a5, v8 &a6, v8 &a7 )
  {
    sw( a0.i[1],a1.i[0] ); sw( a0.i[2],a2.i[0] ); sw( a0.i[3],a3.i[0] ); sw( a0.i[4],a4.i[0] ); sw( a0.i[5],a5.i[0] ); sw( a0.i[6],a6.i[0] ); sw( a0.i[7],a7.i[0] );
                           sw( a1.i[2],a2.i[1] ); sw( a1.i[3],a3.i[1] ); sw( a1.i[4],a4.i[1] ); sw( a1.i[5],a5.i[1] ); sw( a1.i[6],a6.i[1] ); sw( a1.i[7],a7.i[1] );
                                                  sw( a2.i[3],a3.i[2] ); sw( a2.i[4],a4.i[2] ); sw( a2.i[5],a5.i[2] ); sw( a2.i[6],a6.i[2] ); sw( a2.i[7],a7.i[2] );
                                                                         sw( a3.i[4],a4.i[3] ); sw( a3.i[5],a5.i[3] ); sw( a3.i[6],a6.i[3] ); sw( a3.i[7],a7.i[3] );
                                                                                                sw( a4.i[5],a5.i[4] ); sw( a4.i[6],a6.i[4] ); sw( a4.i[7],a7.i[4] );
                                                                                                                       sw( a5.i[6],a6.i[5] ); sw( a5.i[7],a7.i[5] );
                                                                                                                                              sw( a6.i[7],a7.i[6] );
  }

# undef sw

  // v8 memory manipulation functions

  inline void load_8x1( const void * ALIGNED(16) p,
			v8 &a )
  {
    a.i[0] = ((const int * ALIGNED(16))p)[0];
    a.i[1] = ((const int * ALIGNED(16))p)[1];
    a.i[2] = ((const int * ALIGNED(16))p)[2];
    a.i[3] = ((const int * ALIGNED(16))p)[3];
    a.i[4] = ((const int * ALIGNED(16))p)[4];
    a.i[5] = ((const int * ALIGNED(16))p)[5];
    a.i[6] = ((const int * ALIGNED(16))p)[6];
    a.i[7] = ((const int * ALIGNED(16))p)[7];
  }

  inline void store_8x1( const v8 &a,
			 void * ALIGNED(16) p )
  {
    ((int * ALIGNED(16))p)[0] = a.i[0];
    ((int * ALIGNED(16))p)[1] = a.i[1];
    ((int * ALIGNED(16))p)[2] = a.i[2];
    ((int * ALIGNED(16))p)[3] = a.i[3];
    ((int * ALIGNED(16))p)[4] = a.i[4];
    ((int * ALIGNED(16))p)[5] = a.i[5];
    ((int * ALIGNED(16))p)[6] = a.i[6];
    ((int * ALIGNED(16))p)[7] = a.i[7];
  }

  inline void stream_8x1( const v8 &a,
			  void * ALIGNED(16) p )
  {
    ((int * ALIGNED(16))p)[0] = a.i[0];
    ((int * ALIGNED(16))p)[1] = a.i[1];
    ((int * ALIGNED(16))p)[2] = a.i[2];
    ((int * ALIGNED(16))p)[3] = a.i[3];
    ((int * ALIGNED(16))p)[4] = a.i[4];
    ((int * ALIGNED(16))p)[5] = a.i[5];
    ((int * ALIGNED(16))p)[6] = a.i[6];
    ((int * ALIGNED(16))p)[7] = a.i[7];
  }

  inline void clear_8x1( void * ALIGNED(16) p )
  {
    ((int * ALIGNED(16))p)[0] = 0;
    ((int * ALIGNED(16))p)[1] = 0;
    ((int * ALIGNED(16))p)[2] = 0;
    ((int * ALIGNED(16))p)[3] = 0;
    ((int * ALIGNED(16))p)[4] = 0;
    ((int * ALIGNED(16))p)[5] = 0;
    ((int * ALIGNED(16))p)[6] = 0;
    ((int * ALIGNED(16))p)[7] = 0;
  }

  // FIXME: Ordering semantics
  inline void copy_8x1( void * ALIGNED(16) dst,
                        const void * ALIGNED(16) src )
  {
    ((int * ALIGNED(16))dst)[0] = ((const int * ALIGNED(16))src)[0];
    ((int * ALIGNED(16))dst)[1] = ((const int * ALIGNED(16))src)[1];
    ((int * ALIGNED(16))dst)[2] = ((const int * ALIGNED(16))src)[2];
    ((int * ALIGNED(16))dst)[3] = ((const int * ALIGNED(16))src)[3];
    ((int * ALIGNED(16))dst)[4] = ((const int * ALIGNED(16))src)[4];
    ((int * ALIGNED(16))dst)[5] = ((const int * ALIGNED(16))src)[5];
    ((int * ALIGNED(16))dst)[6] = ((const int * ALIGNED(16))src)[6];
    ((int * ALIGNED(16))dst)[7] = ((const int * ALIGNED(16))src)[7];
  }

  inline void swap_8x1( void * ALIGNED(16) a,
			void * ALIGNED(16) b )
  {
    int t;

    t = ((int * ALIGNED(16))a)[0];
    ((int * ALIGNED(16))a)[0] = ((int * ALIGNED(16))b)[0];
    ((int * ALIGNED(16))b)[0] = t;

    t = ((int * ALIGNED(16))a)[1];
    ((int * ALIGNED(16))a)[1] = ((int * ALIGNED(16))b)[1];
    ((int * ALIGNED(16))b)[1] = t;

    t = ((int * ALIGNED(16))a)[2];
    ((int * ALIGNED(16))a)[2] = ((int * ALIGNED(16))b)[2];
    ((int * ALIGNED(16))b)[2] = t;

    t = ((int * ALIGNED(16))a)[3];
    ((int * ALIGNED(16))a)[3] = ((int * ALIGNED(16))b)[3];
    ((int * ALIGNED(16))b)[3] = t;

    t = ((int * ALIGNED(16))a)[4];
    ((int * ALIGNED(16))a)[4] = ((int * ALIGNED(16))b)[4];
    ((int * ALIGNED(16))b)[4] = t;

    t = ((int * ALIGNED(16))a)[5];
    ((int * ALIGNED(16))a)[5] = ((int * ALIGNED(16))b)[5];
    ((int * ALIGNED(16))b)[5] = t;

    t = ((int * ALIGNED(16))a)[6];
    ((int * ALIGNED(16))a)[6] = ((int * ALIGNED(16))b)[6];
    ((int * ALIGNED(16))b)[6] = t;

    t = ((int * ALIGNED(16))a)[7];
    ((int * ALIGNED(16))a)[7] = ((int * ALIGNED(16))b)[7];
    ((int * ALIGNED(16))b)[7] = t;
  }

  // v8 transposed memory manipulation functions

  inline void load_8x1_tr( const void *a0, const void *a1,
                           const void *a2, const void *a3,
                           const void *a4, const void *a5,
                           const void *a6, const void *a7,
			   v8 &a )
  {
    a.i[0] = ((const int *)a0)[0];
    a.i[1] = ((const int *)a1)[0];
    a.i[2] = ((const int *)a2)[0];
    a.i[3] = ((const int *)a3)[0];
    a.i[4] = ((const int *)a4)[0];
    a.i[5] = ((const int *)a5)[0];
    a.i[6] = ((const int *)a6)[0];
    a.i[7] = ((const int *)a7)[0];
  }

  inline void load_8x2_tr( const void * ALIGNED(8) a0,
                           const void * ALIGNED(8) a1,
                           const void * ALIGNED(8) a2,
                           const void * ALIGNED(8) a3,
			   const void * ALIGNED(8) a4,
                           const void * ALIGNED(8) a5,
                           const void * ALIGNED(8) a6,
                           const void * ALIGNED(8) a7,
                           v8 &a, v8 &b )
  {
    a.i[0] = ((const int * ALIGNED(8))a0)[0];
    b.i[0] = ((const int * ALIGNED(8))a0)[1];

    a.i[1] = ((const int * ALIGNED(8))a1)[0];
    b.i[1] = ((const int * ALIGNED(8))a1)[1];

    a.i[2] = ((const int * ALIGNED(8))a2)[0];
    b.i[2] = ((const int * ALIGNED(8))a2)[1];

    a.i[3] = ((const int * ALIGNED(8))a3)[0];
    b.i[3] = ((const int * ALIGNED(8))a3)[1];

    a.i[4] = ((const int * ALIGNED(8))a4)[0];
    b.i[4] = ((const int * ALIGNED(8))a4)[1];

    a.i[5] = ((const int * ALIGNED(8))a5)[0];
    b.i[5] = ((const int * ALIGNED(8))a5)[1];

    a.i[6] = ((const int * ALIGNED(8))a6)[0];
    b.i[6] = ((const int * ALIGNED(8))a6)[1];

    a.i[7] = ((const int * ALIGNED(8))a7)[0];
    b.i[7] = ((const int * ALIGNED(8))a7)[1];
  }

  inline void load_8x3_tr( const void * ALIGNED(16) a0,
                           const void * ALIGNED(16) a1,
                           const void * ALIGNED(16) a2,
                           const void * ALIGNED(16) a3,
 			   const void * ALIGNED(16) a4,
                           const void * ALIGNED(16) a5,
                           const void * ALIGNED(16) a6,
                           const void * ALIGNED(16) a7,
                           v8 &a, v8 &b, v8 &c )
  {
    a.i[0] = ((const int * ALIGNED(16))a0)[0];
    b.i[0] = ((const int * ALIGNED(16))a0)[1];
    c.i[0] = ((const int * ALIGNED(16))a0)[2];

    a.i[1] = ((const int * ALIGNED(16))a1)[0];
    b.i[1] = ((const int * ALIGNED(16))a1)[1];
    c.i[1] = ((const int * ALIGNED(16))a1)[2];

    a.i[2] = ((const int * ALIGNED(16))a2)[0];
    b.i[2] = ((const int * ALIGNED(16))a2)[1];
    c.i[2] = ((const int * ALIGNED(16))a2)[2];

    a.i[3] = ((const int * ALIGNED(16))a3)[0];
    b.i[3] = ((const int * ALIGNED(16))a3)[1];
    c.i[3] = ((const int * ALIGNED(16))a3)[2]; 

    a.i[4] = ((const int * ALIGNED(16))a4)[0];
    b.i[4] = ((const int * ALIGNED(16))a4)[1];
    c.i[4] = ((const int * ALIGNED(16))a4)[2];

    a.i[5] = ((const int * ALIGNED(16))a5)[0];
    b.i[5] = ((const int * ALIGNED(16))a5)[1];
    c.i[5] = ((const int * ALIGNED(16))a5)[2];

    a.i[6] = ((const int * ALIGNED(16))a6)[0];
    b.i[6] = ((const int * ALIGNED(16))a6)[1];
    c.i[6] = ((const int * ALIGNED(16))a6)[2];

    a.i[7] = ((const int * ALIGNED(16))a7)[0];
    b.i[7] = ((const int * ALIGNED(16))a7)[1];
    c.i[7] = ((const int * ALIGNED(16))a7)[2]; 
   }

  inline void load_8x4_tr( const void * ALIGNED(16) a0,
                           const void * ALIGNED(16) a1,
                           const void * ALIGNED(16) a2,
                           const void * ALIGNED(16) a3,
			   const void * ALIGNED(16) a4,
                           const void * ALIGNED(16) a5,
                           const void * ALIGNED(16) a6,
                           const void * ALIGNED(16) a7,
                           v8 &a, v8 &b, v8 &c, v8 &d )
  {
    a.i[0] = ((const int * ALIGNED(16))a0)[0];
    b.i[0] = ((const int * ALIGNED(16))a0)[1];
    c.i[0] = ((const int * ALIGNED(16))a0)[2];
    d.i[0] = ((const int * ALIGNED(16))a0)[3];

    a.i[1] = ((const int * ALIGNED(16))a1)[0];
    b.i[1] = ((const int * ALIGNED(16))a1)[1];
    c.i[1] = ((const int * ALIGNED(16))a1)[2];
    d.i[1] = ((const int * ALIGNED(16))a1)[3];

    a.i[2] = ((const int * ALIGNED(16))a2)[0];
    b.i[2] = ((const int * ALIGNED(16))a2)[1];
    c.i[2] = ((const int * ALIGNED(16))a2)[2];
    d.i[2] = ((const int * ALIGNED(16))a2)[3];

    a.i[3] = ((const int * ALIGNED(16))a3)[0];
    b.i[3] = ((const int * ALIGNED(16))a3)[1];
    c.i[3] = ((const int * ALIGNED(16))a3)[2];
    d.i[3] = ((const int * ALIGNED(16))a3)[3];

    a.i[4] = ((const int * ALIGNED(16))a4)[0];
    b.i[4] = ((const int * ALIGNED(16))a4)[1];
    c.i[4] = ((const int * ALIGNED(16))a4)[2];
    d.i[4] = ((const int * ALIGNED(16))a4)[3];

    a.i[5] = ((const int * ALIGNED(16))a5)[0];
    b.i[5] = ((const int * ALIGNED(16))a5)[1];
    c.i[5] = ((const int * ALIGNED(16))a5)[2];
    d.i[5] = ((const int * ALIGNED(16))a5)[3];

    a.i[6] = ((const int * ALIGNED(16))a6)[0];
    b.i[6] = ((const int * ALIGNED(16))a6)[1];
    c.i[6] = ((const int * ALIGNED(16))a6)[2];
    d.i[6] = ((const int * ALIGNED(16))a6)[3];

    a.i[7] = ((const int * ALIGNED(16))a7)[0];
    b.i[7] = ((const int * ALIGNED(16))a7)[1];
    c.i[7] = ((const int * ALIGNED(16))a7)[2];
    d.i[7] = ((const int * ALIGNED(16))a7)[3];
  }

  inline void load_8x8_tr( const void * ALIGNED(16) a0,
                           const void * ALIGNED(16) a1,
                           const void * ALIGNED(16) a2,
                           const void * ALIGNED(16) a3,
			   const void * ALIGNED(16) a4,
                           const void * ALIGNED(16) a5,
                           const void * ALIGNED(16) a6,
                           const void * ALIGNED(16) a7,
                           v8 &a, v8 &b, v8 &c, v8 &d,
                           v8 &e, v8 &f, v8 &g, v8 &h )
  {
    a.i[0] = ((const int * ALIGNED(16))a0)[0];
    b.i[0] = ((const int * ALIGNED(16))a0)[1];
    c.i[0] = ((const int * ALIGNED(16))a0)[2];
    d.i[0] = ((const int * ALIGNED(16))a0)[3];
    e.i[0] = ((const int * ALIGNED(16))a0)[4];
    f.i[0] = ((const int * ALIGNED(16))a0)[5];
    g.i[0] = ((const int * ALIGNED(16))a0)[6];
    h.i[0] = ((const int * ALIGNED(16))a0)[7];

    a.i[1] = ((const int * ALIGNED(16))a1)[0];
    b.i[1] = ((const int * ALIGNED(16))a1)[1];
    c.i[1] = ((const int * ALIGNED(16))a1)[2];
    d.i[1] = ((const int * ALIGNED(16))a1)[3];
    e.i[1] = ((const int * ALIGNED(16))a1)[4];
    f.i[1] = ((const int * ALIGNED(16))a1)[5];
    g.i[1] = ((const int * ALIGNED(16))a1)[6];
    h.i[1] = ((const int * ALIGNED(16))a1)[7];

    a.i[2] = ((const int * ALIGNED(16))a2)[0];
    b.i[2] = ((const int * ALIGNED(16))a2)[1];
    c.i[2] = ((const int * ALIGNED(16))a2)[2];
    d.i[2] = ((const int * ALIGNED(16))a2)[3];
    e.i[2] = ((const int * ALIGNED(16))a2)[4];
    f.i[2] = ((const int * ALIGNED(16))a2)[5];
    g.i[2] = ((const int * ALIGNED(16))a2)[6];
    h.i[2] = ((const int * ALIGNED(16))a2)[7];

    a.i[3] = ((const int * ALIGNED(16))a3)[0];
    b.i[3] = ((const int * ALIGNED(16))a3)[1];
    c.i[3] = ((const int * ALIGNED(16))a3)[2];
    d.i[3] = ((const int * ALIGNED(16))a3)[3];
    e.i[3] = ((const int * ALIGNED(16))a3)[4];
    f.i[3] = ((const int * ALIGNED(16))a3)[5];
    g.i[3] = ((const int * ALIGNED(16))a3)[6];
    h.i[3] = ((const int * ALIGNED(16))a3)[7];

    a.i[4] = ((const int * ALIGNED(16))a4)[0];
    b.i[4] = ((const int * ALIGNED(16))a4)[1];
    c.i[4] = ((const int * ALIGNED(16))a4)[2];
    d.i[4] = ((const int * ALIGNED(16))a4)[3];
    e.i[4] = ((const int * ALIGNED(16))a4)[4];
    f.i[4] = ((const int * ALIGNED(16))a4)[5];
    g.i[4] = ((const int * ALIGNED(16))a4)[6];
    h.i[4] = ((const int * ALIGNED(16))a4)[7];

    a.i[5] = ((const int * ALIGNED(16))a5)[0];
    b.i[5] = ((const int * ALIGNED(16))a5)[1];
    c.i[5] = ((const int * ALIGNED(16))a5)[2];
    d.i[5] = ((const int * ALIGNED(16))a5)[3];
    e.i[5] = ((const int * ALIGNED(16))a5)[4];
    f.i[5] = ((const int * ALIGNED(16))a5)[5];
    g.i[5] = ((const int * ALIGNED(16))a5)[6];
    h.i[5] = ((const int * ALIGNED(16))a5)[7];

    a.i[6] = ((const int * ALIGNED(16))a6)[0];
    b.i[6] = ((const int * ALIGNED(16))a6)[1];
    c.i[6] = ((const int * ALIGNED(16))a6)[2];
    d.i[6] = ((const int * ALIGNED(16))a6)[3];
    e.i[6] = ((const int * ALIGNED(16))a6)[4];
    f.i[6] = ((const int * ALIGNED(16))a6)[5];
    g.i[6] = ((const int * ALIGNED(16))a6)[6];
    h.i[6] = ((const int * ALIGNED(16))a6)[7];

    a.i[7] = ((const int * ALIGNED(16))a7)[0];
    b.i[7] = ((const int * ALIGNED(16))a7)[1];
    c.i[7] = ((const int * ALIGNED(16))a7)[2];
    d.i[7] = ((const int * ALIGNED(16))a7)[3];
    e.i[7] = ((const int * ALIGNED(16))a7)[4];
    f.i[7] = ((const int * ALIGNED(16))a7)[5];
    g.i[7] = ((const int * ALIGNED(16))a7)[6];
    h.i[7] = ((const int * ALIGNED(16))a7)[7];
  }

  inline void store_8x1_tr( const v8 &a,
                            void *a0, void *a1, void *a2, void *a3,
                            void *a4, void *a5, void *a6, void *a7 )
  {
    ((int *)a0)[0] = a.i[0];
    ((int *)a1)[0] = a.i[1];
    ((int *)a2)[0] = a.i[2];
    ((int *)a3)[0] = a.i[3];
    ((int *)a4)[0] = a.i[4];
    ((int *)a5)[0] = a.i[5];
    ((int *)a6)[0] = a.i[6];
    ((int *)a7)[0] = a.i[7];
  }

  inline void store_8x2_tr( const v8 &a, const v8 &b,
                            void * ALIGNED(8) a0, void * ALIGNED(8) a1,
                            void * ALIGNED(8) a2, void * ALIGNED(8) a3,
                            void * ALIGNED(8) a4, void * ALIGNED(8) a5,
                            void * ALIGNED(8) a6, void * ALIGNED(8) a7 )
  {
    ((int * ALIGNED(8))a0)[0] = a.i[0];
    ((int * ALIGNED(8))a0)[1] = b.i[0];

    ((int * ALIGNED(8))a1)[0] = a.i[1];
    ((int * ALIGNED(8))a1)[1] = b.i[1];

    ((int * ALIGNED(8))a2)[0] = a.i[2];
    ((int * ALIGNED(8))a2)[1] = b.i[2];

    ((int * ALIGNED(8))a3)[0] = a.i[3];
    ((int * ALIGNED(8))a3)[1] = b.i[3];

    ((int * ALIGNED(8))a4)[0] = a.i[4];
    ((int * ALIGNED(8))a4)[1] = b.i[4];

    ((int * ALIGNED(8))a5)[0] = a.i[5];
    ((int * ALIGNED(8))a5)[1] = b.i[5];

    ((int * ALIGNED(8))a6)[0] = a.i[6];
    ((int * ALIGNED(8))a6)[1] = b.i[6];

    ((int * ALIGNED(8))a7)[0] = a.i[7];
    ((int * ALIGNED(8))a7)[1] = b.i[7];
  }

  inline void store_8x3_tr( const v8 &a, const v8 &b, const v8 &c,
                            void * ALIGNED(16) a0, void * ALIGNED(16) a1,
                            void * ALIGNED(16) a2, void * ALIGNED(16) a3,
                            void * ALIGNED(16) a4, void * ALIGNED(16) a5,
                            void * ALIGNED(16) a6, void * ALIGNED(16) a7 )
  {
    ((int * ALIGNED(16))a0)[0] = a.i[0];
    ((int * ALIGNED(16))a0)[1] = b.i[0];
    ((int * ALIGNED(16))a0)[2] = c.i[0];

    ((int * ALIGNED(16))a1)[0] = a.i[1];
    ((int * ALIGNED(16))a1)[1] = b.i[1];
    ((int * ALIGNED(16))a1)[2] = c.i[1];

    ((int * ALIGNED(16))a2)[0] = a.i[2];
    ((int * ALIGNED(16))a2)[1] = b.i[2];
    ((int * ALIGNED(16))a2)[2] = c.i[2];

    ((int * ALIGNED(16))a3)[0] = a.i[3];
    ((int * ALIGNED(16))a3)[1] = b.i[3];
    ((int * ALIGNED(16))a3)[2] = c.i[3];

    ((int * ALIGNED(16))a4)[0] = a.i[4];
    ((int * ALIGNED(16))a4)[1] = b.i[4];
    ((int * ALIGNED(16))a4)[2] = c.i[4];

    ((int * ALIGNED(16))a5)[0] = a.i[5];
    ((int * ALIGNED(16))a5)[1] = b.i[5];
    ((int * ALIGNED(16))a5)[2] = c.i[5];

    ((int * ALIGNED(16))a6)[0] = a.i[6];
    ((int * ALIGNED(16))a6)[1] = b.i[6];
    ((int * ALIGNED(16))a6)[2] = c.i[6];

    ((int * ALIGNED(16))a7)[0] = a.i[7];
    ((int * ALIGNED(16))a7)[1] = b.i[7];
    ((int * ALIGNED(16))a7)[2] = c.i[7];
  }

  inline void store_8x4_tr( const v8 &a, const v8 &b, const v8 &c, const v8 &d,
                            void * ALIGNED(16) a0, void * ALIGNED(16) a1,
                            void * ALIGNED(16) a2, void * ALIGNED(16) a3,
                            void * ALIGNED(16) a4, void * ALIGNED(16) a5,
                            void * ALIGNED(16) a6, void * ALIGNED(16) a7 )
  {
    ((int * ALIGNED(16))a0)[0] = a.i[0];
    ((int * ALIGNED(16))a0)[1] = b.i[0];
    ((int * ALIGNED(16))a0)[2] = c.i[0];
    ((int * ALIGNED(16))a0)[3] = d.i[0];

    ((int * ALIGNED(16))a1)[0] = a.i[1];
    ((int * ALIGNED(16))a1)[1] = b.i[1];
    ((int * ALIGNED(16))a1)[2] = c.i[1];
    ((int * ALIGNED(16))a1)[3] = d.i[1];

    ((int * ALIGNED(16))a2)[0] = a.i[2];
    ((int * ALIGNED(16))a2)[1] = b.i[2];
    ((int * ALIGNED(16))a2)[2] = c.i[2];
    ((int * ALIGNED(16))a2)[3] = d.i[2];

    ((int * ALIGNED(16))a3)[0] = a.i[3];
    ((int * ALIGNED(16))a3)[1] = b.i[3];
    ((int * ALIGNED(16))a3)[2] = c.i[3];
    ((int * ALIGNED(16))a3)[3] = d.i[3];

    ((int * ALIGNED(16))a4)[0] = a.i[4];
    ((int * ALIGNED(16))a4)[1] = b.i[4];
    ((int * ALIGNED(16))a4)[2] = c.i[4];
    ((int * ALIGNED(16))a4)[3] = d.i[4];

    ((int * ALIGNED(16))a5)[0] = a.i[5];
    ((int * ALIGNED(16))a5)[1] = b.i[5];
    ((int * ALIGNED(16))a5)[2] = c.i[5];
    ((int * ALIGNED(16))a5)[3] = d.i[5];

    ((int * ALIGNED(16))a6)[0] = a.i[6];
    ((int * ALIGNED(16))a6)[1] = b.i[6];
    ((int * ALIGNED(16))a6)[2] = c.i[6];
    ((int * ALIGNED(16))a6)[3] = d.i[6];

    ((int * ALIGNED(16))a7)[0] = a.i[7];
    ((int * ALIGNED(16))a7)[1] = b.i[7];
    ((int * ALIGNED(16))a7)[2] = c.i[7];
    ((int * ALIGNED(16))a7)[3] = d.i[7];
  }

  inline void store_8x8_tr( const v8 &a, const v8 &b, const v8 &c, const v8 &d,
			    const v8 &e, const v8 &f, const v8 &g, const v8 &h,
                            void * ALIGNED(16) a0, void * ALIGNED(16) a1,
                            void * ALIGNED(16) a2, void * ALIGNED(16) a3,
                            void * ALIGNED(16) a4, void * ALIGNED(16) a5,
                            void * ALIGNED(16) a6, void * ALIGNED(16) a7 )
  {
    ((int * ALIGNED(16))a0)[0] = a.i[0];
    ((int * ALIGNED(16))a0)[1] = b.i[0];
    ((int * ALIGNED(16))a0)[2] = c.i[0];
    ((int * ALIGNED(16))a0)[3] = d.i[0];
    ((int * ALIGNED(16))a0)[4] = e.i[0];
    ((int * ALIGNED(16))a0)[5] = f.i[0];
    ((int * ALIGNED(16))a0)[6] = g.i[0];
    ((int * ALIGNED(16))a0)[7] = h.i[0];

    ((int * ALIGNED(16))a1)[0] = a.i[1];
    ((int * ALIGNED(16))a1)[1] = b.i[1];
    ((int * ALIGNED(16))a1)[2] = c.i[1];
    ((int * ALIGNED(16))a1)[3] = d.i[1];
    ((int * ALIGNED(16))a1)[4] = e.i[1];
    ((int * ALIGNED(16))a1)[5] = f.i[1];
    ((int * ALIGNED(16))a1)[6] = g.i[1];
    ((int * ALIGNED(16))a1)[7] = h.i[1];

    ((int * ALIGNED(16))a2)[0] = a.i[2];
    ((int * ALIGNED(16))a2)[1] = b.i[2];
    ((int * ALIGNED(16))a2)[2] = c.i[2];
    ((int * ALIGNED(16))a2)[3] = d.i[2];
    ((int * ALIGNED(16))a2)[4] = e.i[2];
    ((int * ALIGNED(16))a2)[5] = f.i[2];
    ((int * ALIGNED(16))a2)[6] = g.i[2];
    ((int * ALIGNED(16))a2)[7] = h.i[2];

    ((int * ALIGNED(16))a3)[0] = a.i[3];
    ((int * ALIGNED(16))a3)[1] = b.i[3];
    ((int * ALIGNED(16))a3)[2] = c.i[3];
    ((int * ALIGNED(16))a3)[3] = d.i[3];
    ((int * ALIGNED(16))a3)[4] = e.i[3];
    ((int * ALIGNED(16))a3)[5] = f.i[3];
    ((int * ALIGNED(16))a3)[6] = g.i[3];
    ((int * ALIGNED(16))a3)[7] = h.i[3];

    ((int * ALIGNED(16))a4)[0] = a.i[4];
    ((int * ALIGNED(16))a4)[1] = b.i[4];
    ((int * ALIGNED(16))a4)[2] = c.i[4];
    ((int * ALIGNED(16))a4)[3] = d.i[4];
    ((int * ALIGNED(16))a4)[4] = e.i[4];
    ((int * ALIGNED(16))a4)[5] = f.i[4];
    ((int * ALIGNED(16))a4)[6] = g.i[4];
    ((int * ALIGNED(16))a4)[7] = h.i[4];

    ((int * ALIGNED(16))a5)[0] = a.i[5];
    ((int * ALIGNED(16))a5)[1] = b.i[5];
    ((int * ALIGNED(16))a5)[2] = c.i[5];
    ((int * ALIGNED(16))a5)[3] = d.i[5];
    ((int * ALIGNED(16))a5)[4] = e.i[5];
    ((int * ALIGNED(16))a5)[5] = f.i[5];
    ((int * ALIGNED(16))a5)[6] = g.i[5];
    ((int * ALIGNED(16))a5)[7] = h.i[5];

    ((int * ALIGNED(16))a6)[0] = a.i[6];
    ((int * ALIGNED(16))a6)[1] = b.i[6];
    ((int * ALIGNED(16))a6)[2] = c.i[6];
    ((int * ALIGNED(16))a6)[3] = d.i[6];
    ((int * ALIGNED(16))a6)[4] = e.i[6];
    ((int * ALIGNED(16))a6)[5] = f.i[6];
    ((int * ALIGNED(16))a6)[6] = g.i[6];
    ((int * ALIGNED(16))a6)[7] = h.i[6];

    ((int * ALIGNED(16))a7)[0] = a.i[7];
    ((int * ALIGNED(16))a7)[1] = b.i[7];
    ((int * ALIGNED(16))a7)[2] = c.i[7];
    ((int * ALIGNED(16))a7)[3] = d.i[7];
    ((int * ALIGNED(16))a7)[4] = e.i[7];
    ((int * ALIGNED(16))a7)[5] = f.i[7];
    ((int * ALIGNED(16))a7)[6] = g.i[7];
    ((int * ALIGNED(16))a7)[7] = h.i[7];
  }

  //////////////
  // v8int class

  class v8int : public v8
  {
    // v8int prefix unary operator friends

    friend inline v8int operator  +( const v8int & a ) ALWAYS_INLINE;
    friend inline v8int operator  -( const v8int & a ) ALWAYS_INLINE;
    friend inline v8int operator  ~( const v8int & a ) ALWAYS_INLINE;
    friend inline v8int operator  !( const v8int & a ) ALWAYS_INLINE;
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector

    // v8int prefix increment / decrement operator friends

    friend inline v8int operator ++( v8int & a ) ALWAYS_INLINE;
    friend inline v8int operator --( v8int & a ) ALWAYS_INLINE;

    // v8int postfix increment / decrement operator friends

    friend inline v8int operator ++( v8int & a, int ) ALWAYS_INLINE;
    friend inline v8int operator --( v8int & a, int ) ALWAYS_INLINE;

    // v8int binary operator friends

    friend inline v8int operator  +( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator  -( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator  *( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator  /( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator  %( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator  ^( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator  &( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator  |( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator <<( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator >>( const v8int &a, const v8int &b ) ALWAYS_INLINE;

    // v8int logical operator friends

    friend inline v8int operator  <( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator  >( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator ==( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator !=( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator <=( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator >=( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator &&( const v8int &a, const v8int &b ) ALWAYS_INLINE;
    friend inline v8int operator ||( const v8int &a, const v8int &b ) ALWAYS_INLINE;

    // v8int miscellaneous friends

    friend inline v8int abs( const v8int &a ) ALWAYS_INLINE;
    friend inline v8    czero( const v8int &c, const v8 &a ) ALWAYS_INLINE;
    friend inline v8 notczero( const v8int &c, const v8 &a ) ALWAYS_INLINE;
    // FIXME: cswap, notcswap!
    friend inline v8 merge( const v8int &c, const v8 &t, const v8 &f ) ALWAYS_INLINE;

    // v8float unary operator friends

    friend inline v8int operator  !( const v8float & a ) ALWAYS_INLINE;

    // v8float logical operator friends

    friend inline v8int operator  <( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8int operator  >( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8int operator ==( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8int operator !=( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8int operator <=( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8int operator >=( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8int operator &&( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8int operator ||( const v8float &a, const v8float &b ) ALWAYS_INLINE;

    // v8float miscellaneous friends

    friend inline v8float clear_bits(  const v8int &m, const v8float &a ) ALWAYS_INLINE;
    friend inline v8float set_bits(    const v8int &m, const v8float &a ) ALWAYS_INLINE;
    friend inline v8float toggle_bits( const v8int &m, const v8float &a ) ALWAYS_INLINE;

  public:

    // v8int constructors / destructors

    v8int() {}                                // Default constructor

    v8int( const v8int &a ) : v8(a) {}        // Copy constructor from v8

    v8int( const v8 &a )                      // Init from mixed
    {
      i[0] = a.i[0]; i[1] = a.i[1]; i[2] = a.i[2]; i[3] = a.i[3];
      i[4] = a.i[4]; i[5] = a.i[5]; i[6] = a.i[6]; i[7] = a.i[7];
    }

    v8int( int a )                            // Init from scalar
    {
      i[0] = a; i[1] = a; i[2] = a; i[3] = a;
      i[4] = a; i[5] = a; i[6] = a; i[7] = a;
    }

    v8int( int i0, int i1, int i2, int i3,
	   int i4, int i5, int i6, int i7 )   // Init from scalars
    {
      i[0] = i0; i[1] = i1; i[2] = i2; i[3] = i3;
      i[4] = i4; i[5] = i5; i[6] = i6; i[7] = i7;
    }

    ~v8int() {}                               // Destructor

    // v8int assignment operators

#   define ASSIGN(op)			          \
    inline v8int &operator op( const v8int &b )   \
    {						  \
      i[0] op b.i[0];                             \
      i[1] op b.i[1];                             \
      i[2] op b.i[2];                             \
      i[3] op b.i[3];                             \
      i[4] op b.i[4];                             \
      i[5] op b.i[5];                             \
      i[6] op b.i[6];                             \
      i[7] op b.i[7];                             \
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

    // v8int member access operator

    inline int &operator []( int n )
    {
      return i[n];
    }

    inline int  operator ()( int n )
    {
      return i[n];
    }
  };

  // v8int prefix unary operators

# define PREFIX_UNARY(op)                       \
  inline v8int operator op( const v8int & a )   \
  {						\
    v8int b;                                    \
    b.i[0] = ( op a.i[0] );                     \
    b.i[1] = ( op a.i[1] );                     \
    b.i[2] = ( op a.i[2] );                     \
    b.i[3] = ( op a.i[3] );                     \
    b.i[4] = ( op a.i[4] );                     \
    b.i[5] = ( op a.i[5] );                     \
    b.i[6] = ( op a.i[6] );                     \
    b.i[7] = ( op a.i[7] );                     \
    return b;                                   \
  }

  PREFIX_UNARY(+)
  PREFIX_UNARY(-)

  inline v8int operator !( const v8int & a )
  {
    v8int b;

    b.i[0] = - ( !a.i[0] );
    b.i[1] = - ( !a.i[1] );
    b.i[2] = - ( !a.i[2] );
    b.i[3] = - ( !a.i[3] );
    b.i[4] = - ( !a.i[4] );
    b.i[5] = - ( !a.i[5] );
    b.i[6] = - ( !a.i[6] );
    b.i[7] = - ( !a.i[7] );

    return b;
  }

  PREFIX_UNARY(~)

# undef PREFIX_UNARY

  // v8int prefix increment / decrement

# define PREFIX_INCDEC(op)                      \
  inline v8int operator op( v8int & a )         \
  {						\
    v8int b;                                    \
    b.i[0] = ( op a.i[0] );                     \
    b.i[1] = ( op a.i[1] );                     \
    b.i[2] = ( op a.i[2] );                     \
    b.i[3] = ( op a.i[3] );                     \
    b.i[4] = ( op a.i[4] );                     \
    b.i[5] = ( op a.i[5] );                     \
    b.i[6] = ( op a.i[6] );                     \
    b.i[7] = ( op a.i[7] );                     \
    return b;                                   \
  }

  PREFIX_INCDEC(++)
  PREFIX_INCDEC(--)

# undef PREFIX_INCDEC

  // v8int postfix increment / decrement

# define POSTFIX_INCDEC(op)                    \
  inline v8int operator op( v8int & a, int )   \
  {					       \
    v8int b;                                   \
    b.i[0] = ( a.i[0] op );                    \
    b.i[1] = ( a.i[1] op );                    \
    b.i[2] = ( a.i[2] op );                    \
    b.i[3] = ( a.i[3] op );                    \
    b.i[4] = ( a.i[4] op );                    \
    b.i[5] = ( a.i[5] op );                    \
    b.i[6] = ( a.i[6] op );                    \
    b.i[7] = ( a.i[7] op );                    \
    return b;                                  \
  }

  POSTFIX_INCDEC(++)
  POSTFIX_INCDEC(--)

# undef POSTFIX_INCDEC

  // v8int binary operators

# define BINARY(op)                                             \
  inline v8int operator op( const v8int &a, const v8int &b )    \
  {								\
    v8int c;                                                    \
    c.i[0] = a.i[0] op b.i[0];                                  \
    c.i[1] = a.i[1] op b.i[1];                                  \
    c.i[2] = a.i[2] op b.i[2];                                  \
    c.i[3] = a.i[3] op b.i[3];                                  \
    c.i[4] = a.i[4] op b.i[4];                                  \
    c.i[5] = a.i[5] op b.i[5];                                  \
    c.i[6] = a.i[6] op b.i[6];                                  \
    c.i[7] = a.i[7] op b.i[7];                                  \
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

  // v8int logical operators

# define LOGICAL(op)                                           \
  inline v8int operator op( const v8int &a, const v8int &b )   \
  {							       \
    v8int c;                                                   \
    c.i[0] = - ( a.i[0] op b.i[0] );                           \
    c.i[1] = - ( a.i[1] op b.i[1] );                           \
    c.i[2] = - ( a.i[2] op b.i[2] );                           \
    c.i[3] = - ( a.i[3] op b.i[3] );                           \
    c.i[4] = - ( a.i[4] op b.i[4] );                           \
    c.i[5] = - ( a.i[5] op b.i[5] );                           \
    c.i[6] = - ( a.i[6] op b.i[6] );                           \
    c.i[7] = - ( a.i[7] op b.i[7] );                           \
    return c;                                                  \
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

  // v8int miscellaneous functions

  inline v8int abs( const v8int &a )
  {
    v8int b;

    b.i[0] = ( a.i[0] >= 0 ) ? a.i[0] : -a.i[0];
    b.i[1] = ( a.i[1] >= 0 ) ? a.i[1] : -a.i[1];
    b.i[2] = ( a.i[2] >= 0 ) ? a.i[2] : -a.i[2];
    b.i[3] = ( a.i[3] >= 0 ) ? a.i[3] : -a.i[3];
    b.i[4] = ( a.i[4] >= 0 ) ? a.i[4] : -a.i[4];
    b.i[5] = ( a.i[5] >= 0 ) ? a.i[5] : -a.i[5];
    b.i[6] = ( a.i[6] >= 0 ) ? a.i[6] : -a.i[6];
    b.i[7] = ( a.i[7] >= 0 ) ? a.i[7] : -a.i[7];

    return b;
  }

  inline v8 czero( const v8int &c, const v8 &a )
  {
    v8 b;

    b.i[0] = a.i[0] & ~c.i[0];
    b.i[1] = a.i[1] & ~c.i[1];
    b.i[2] = a.i[2] & ~c.i[2];
    b.i[3] = a.i[3] & ~c.i[3];
    b.i[4] = a.i[4] & ~c.i[4];
    b.i[5] = a.i[5] & ~c.i[5];
    b.i[6] = a.i[6] & ~c.i[6];
    b.i[7] = a.i[7] & ~c.i[7];

    return b;
  }

  inline v8 notczero( const v8int &c, const v8 &a )
  {
    v8 b;

    b.i[0] = a.i[0] & c.i[0];
    b.i[1] = a.i[1] & c.i[1];
    b.i[2] = a.i[2] & c.i[2];
    b.i[3] = a.i[3] & c.i[3];
    b.i[4] = a.i[4] & c.i[4];
    b.i[5] = a.i[5] & c.i[5];
    b.i[6] = a.i[6] & c.i[6];
    b.i[7] = a.i[7] & c.i[7];

    return b;
  }

  inline v8 merge( const v8int &c, const v8 &t, const v8 &f )
  {
    v8 m;

    m.i[0] = ( f.i[0] & ~c.i[0] ) | ( t.i[0] & c.i[0] );
    m.i[1] = ( f.i[1] & ~c.i[1] ) | ( t.i[1] & c.i[1] );
    m.i[2] = ( f.i[2] & ~c.i[2] ) | ( t.i[2] & c.i[2] );
    m.i[3] = ( f.i[3] & ~c.i[3] ) | ( t.i[3] & c.i[3] );
    m.i[4] = ( f.i[4] & ~c.i[4] ) | ( t.i[4] & c.i[4] );
    m.i[5] = ( f.i[5] & ~c.i[5] ) | ( t.i[5] & c.i[5] );
    m.i[6] = ( f.i[6] & ~c.i[6] ) | ( t.i[6] & c.i[6] );
    m.i[7] = ( f.i[7] & ~c.i[7] ) | ( t.i[7] & c.i[7] );

    return m;
  }

  ////////////////
  // v8float class

  class v8float : public v8
  {
    // v8float prefix unary operator friends

    friend inline v8float operator  +( const v8float &a ) ALWAYS_INLINE;
    friend inline v8float operator  -( const v8float &a ) ALWAYS_INLINE;
    friend inline v8float operator  ~( const v8float &a ) ALWAYS_INLINE;
    friend inline v8int   operator  !( const v8float &a ) ALWAYS_INLINE;
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector

    // v8float prefix increment / decrement operator friends

    friend inline v8float operator ++( v8float &a ) ALWAYS_INLINE;
    friend inline v8float operator --( v8float &a ) ALWAYS_INLINE;

    // v8float postfix increment / decrement operator friends

    friend inline v8float operator ++( v8float &a, int ) ALWAYS_INLINE;
    friend inline v8float operator --( v8float &a, int ) ALWAYS_INLINE;

    // v8float binary operator friends

    friend inline v8float operator  +( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8float operator  -( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8float operator  *( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8float operator  /( const v8float &a, const v8float &b ) ALWAYS_INLINE;

    // v8float logical operator friends

    friend inline v8int operator  <( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8int operator  >( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8int operator ==( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8int operator !=( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8int operator <=( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8int operator >=( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8int operator &&( const v8float &a, const v8float &b ) ALWAYS_INLINE;
    friend inline v8int operator ||( const v8float &a, const v8float &b ) ALWAYS_INLINE;

    // v8float math library friends

#   define CMATH_FR1(fn) friend inline v8float fn( const v8float &a ) ALWAYS_INLINE
#   define CMATH_FR2(fn) friend inline v8float fn( const v8float &a,  \
                                                   const v8float &b ) ALWAYS_INLINE

    CMATH_FR1(acos);  CMATH_FR1(asin);  CMATH_FR1(atan); CMATH_FR2(atan2);
    CMATH_FR1(ceil);  CMATH_FR1(cos);   CMATH_FR1(cosh); CMATH_FR1(exp);
    CMATH_FR1(fabs);  CMATH_FR1(floor); CMATH_FR2(fmod); CMATH_FR1(log);
    CMATH_FR1(log10); CMATH_FR2(pow);   CMATH_FR1(sin);  CMATH_FR1(sinh);
    CMATH_FR1(sqrt);  CMATH_FR1(tan);   CMATH_FR1(tanh);

    CMATH_FR2(copysign);

#   undef CMATH_FR1
#   undef CMATH_FR2

    // v8float miscellaneous friends

    friend inline v8float rsqrt_approx( const v8float &a ) ALWAYS_INLINE;
    friend inline v8float rsqrt       ( const v8float &a ) ALWAYS_INLINE;
    friend inline v8float rcp_approx( const v8float &a ) ALWAYS_INLINE;
    friend inline v8float rcp       ( const v8float &a ) ALWAYS_INLINE;
    friend inline v8float fma ( const v8float &a, const v8float &b, const v8float &c ) ALWAYS_INLINE;
    friend inline v8float fms ( const v8float &a, const v8float &b, const v8float &c ) ALWAYS_INLINE;
    friend inline v8float fnms( const v8float &a, const v8float &b, const v8float &c ) ALWAYS_INLINE;
    friend inline v8float  clear_bits( const v8int &m, const v8float &a ) ALWAYS_INLINE;
    friend inline v8float    set_bits( const v8int &m, const v8float &a ) ALWAYS_INLINE;
    friend inline v8float toggle_bits( const v8int &m, const v8float &a ) ALWAYS_INLINE;
    friend inline void increment_8x1( float * ALIGNED(16) p, const v8float &a ) ALWAYS_INLINE;
    friend inline void decrement_8x1( float * ALIGNED(16) p, const v8float &a ) ALWAYS_INLINE;
    friend inline void     scale_8x1( float * ALIGNED(16) p, const v8float &a ) ALWAYS_INLINE;

  public:

    // v8float constructors / destructors

    v8float() {}                                        // Default constructor

    v8float( const v8float &a ) : v8()                  // Copy constructor
    {
      f[0] = a.f[0];
      f[1] = a.f[1];
      f[2] = a.f[2];
      f[3] = a.f[3];
      f[4] = a.f[4];
      f[5] = a.f[5];
      f[6] = a.f[6];
      f[7] = a.f[7];
    }

    v8float( const v8 &a )                              // Init from mixed
    {
      f[0] = a.f[0];
      f[1] = a.f[1];
      f[2] = a.f[2];
      f[3] = a.f[3];
      f[4] = a.f[4];
      f[5] = a.f[5];
      f[6] = a.f[6];
      f[7] = a.f[7];
    }

    v8float( float a )                                  // Init from scalar
    {
      f[0] = a;
      f[1] = a;
      f[2] = a;
      f[3] = a;
      f[4] = a;
      f[5] = a;
      f[6] = a;
      f[7] = a;
    }

    v8float( float f0, float f1, float f2, float f3,
	     float f4, float f5, float f6, float f7 )   // Init from scalars
    {
      f[0] = f0;
      f[1] = f1;
      f[2] = f2;
      f[3] = f3;
      f[4] = f4;
      f[5] = f5;
      f[6] = f6;
      f[7] = f7;
    }

    ~v8float() {}                                       // Destructor

    // v8float assignment operators

#   define ASSIGN(op)                                   \
    inline v8float &operator op( const v8float &b )     \
    {							\
      f[0] op b.f[0];		             		\
      f[1] op b.f[1];                                   \
      f[2] op b.f[2];                                   \
      f[3] op b.f[3];                                   \
      f[4] op b.f[4];		             		\
      f[5] op b.f[5];                                   \
      f[6] op b.f[6];                                   \
      f[7] op b.f[7];                                   \
      return *this;                                     \
    }

    ASSIGN(=)
    ASSIGN(+=)
    ASSIGN(-=)
    ASSIGN(*=)
    ASSIGN(/=)

#   undef ASSIGN

    // v8float member access operator

    inline float &operator []( int n )
    {
      return f[n];
    }

    inline float  operator ()( int n )
    {
      return f[n];
    }
  };

  // v8float prefix unary operators

  inline v8float operator +( const v8float &a )
  {
    v8float b;

    b.f[0] = +a.f[0];
    b.f[1] = +a.f[1];
    b.f[2] = +a.f[2];
    b.f[3] = +a.f[3];
    b.f[4] = +a.f[4];
    b.f[5] = +a.f[5];
    b.f[6] = +a.f[6];
    b.f[7] = +a.f[7];

    return b;
  }

  inline v8float operator -( const v8float &a )
  {
    v8float b;

    b.f[0] = -a.f[0];
    b.f[1] = -a.f[1];
    b.f[2] = -a.f[2];
    b.f[3] = -a.f[3];
    b.f[4] = -a.f[4];
    b.f[5] = -a.f[5];
    b.f[6] = -a.f[6];
    b.f[7] = -a.f[7];

    return b;
  }

  inline v8int operator !( const v8float &a )
  {
    v8int b;

    b.i[0] = a.i[0] ? 0 : -1;
    b.i[1] = a.i[1] ? 0 : -1;
    b.i[2] = a.i[2] ? 0 : -1;
    b.i[3] = a.i[3] ? 0 : -1;
    b.i[4] = a.i[4] ? 0 : -1;
    b.i[5] = a.i[5] ? 0 : -1;
    b.i[6] = a.i[6] ? 0 : -1;
    b.i[7] = a.i[7] ? 0 : -1;

    return b;
  }

  // v8float prefix increment / decrement operators

  inline v8float operator ++( v8float &a )
  {
    v8float b;

    b.f[0] = ++a.f[0];
    b.f[1] = ++a.f[1];
    b.f[2] = ++a.f[2];
    b.f[3] = ++a.f[3];
    b.f[4] = ++a.f[4];
    b.f[5] = ++a.f[5];
    b.f[6] = ++a.f[6];
    b.f[7] = ++a.f[7];

    return b;
  }

  inline v8float operator --( v8float &a )
  {
    v8float b;

    b.f[0] = --a.f[0];
    b.f[1] = --a.f[1];
    b.f[2] = --a.f[2];
    b.f[3] = --a.f[3];
    b.f[4] = --a.f[4];
    b.f[5] = --a.f[5];
    b.f[6] = --a.f[6];
    b.f[7] = --a.f[7];

    return b;
  }

  // v8float postfix increment / decrement operators

  inline v8float operator ++( v8float &a, int )
  {
    v8float b;

    b.f[0] = a.f[0]++;
    b.f[1] = a.f[1]++;
    b.f[2] = a.f[2]++;
    b.f[3] = a.f[3]++;
    b.f[4] = a.f[4]++;
    b.f[5] = a.f[5]++;
    b.f[6] = a.f[6]++;
    b.f[7] = a.f[7]++;

    return b;
  }

  inline v8float operator --( v8float &a, int )
  {
    v8float b;

    b.f[0] = a.f[0]--;
    b.f[1] = a.f[1]--;
    b.f[2] = a.f[2]--;
    b.f[3] = a.f[3]--;
    b.f[4] = a.f[4]--;
    b.f[5] = a.f[5]--;
    b.f[6] = a.f[6]--;
    b.f[7] = a.f[7]--;

    return b;
  }

  // v8float binary operators

# define BINARY(op)                                                  \
  inline v8float operator op( const v8float &a, const v8float &b )   \
  {								     \
    v8float c;                                                       \
    c.f[0] = a.f[0] op b.f[0];                                       \
    c.f[1] = a.f[1] op b.f[1];                                       \
    c.f[2] = a.f[2] op b.f[2];                                       \
    c.f[3] = a.f[3] op b.f[3];                                       \
    c.f[4] = a.f[4] op b.f[4];                                       \
    c.f[5] = a.f[5] op b.f[5];                                       \
    c.f[6] = a.f[6] op b.f[6];                                       \
    c.f[7] = a.f[7] op b.f[7];                                       \
    return c;                                                        \
  }

  BINARY(+)
  BINARY(-)
  BINARY(*)
  BINARY(/)

# undef BINARY

  // v8float logical operators

# define LOGICAL(op)                                               \
  inline v8int operator op( const v8float &a, const v8float &b )   \
  {								   \
    v8int c;                                                       \
    c.i[0] = -( a.f[0] op b.f[0] );                                \
    c.i[1] = -( a.f[1] op b.f[1] );                                \
    c.i[2] = -( a.f[2] op b.f[2] );                                \
    c.i[3] = -( a.f[3] op b.f[3] );                                \
    c.i[4] = -( a.f[4] op b.f[4] );                                \
    c.i[5] = -( a.f[5] op b.f[5] );                                \
    c.i[6] = -( a.f[6] op b.f[6] );                                \
    c.i[7] = -( a.f[7] op b.f[7] );                                \
    return c;                                                      \
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

  // v8float math library functions

# define CMATH_FR1(fn)                          \
  inline v8float fn( const v8float &a )         \
  {						\
    v8float b;                                  \
    b.f[0] = ::fn( a.f[0] );                    \
    b.f[1] = ::fn( a.f[1] );                    \
    b.f[2] = ::fn( a.f[2] );                    \
    b.f[3] = ::fn( a.f[3] );                    \
    b.f[4] = ::fn( a.f[4] );                    \
    b.f[5] = ::fn( a.f[5] );                    \
    b.f[6] = ::fn( a.f[6] );                    \
    b.f[7] = ::fn( a.f[7] );                    \
    return b;                                   \
  }

# define CMATH_FR2(fn)                                          \
  inline v8float fn( const v8float &a, const v8float &b )       \
  {								\
    v8float c;                                                  \
    c.f[0] = ::fn( a.f[0], b.f[0] );                            \
    c.f[1] = ::fn( a.f[1], b.f[1] );                            \
    c.f[2] = ::fn( a.f[2], b.f[2] );                            \
    c.f[3] = ::fn( a.f[3], b.f[3] );                            \
    c.f[4] = ::fn( a.f[4], b.f[4] );                            \
    c.f[5] = ::fn( a.f[5], b.f[5] );                            \
    c.f[6] = ::fn( a.f[6], b.f[6] );                            \
    c.f[7] = ::fn( a.f[7], b.f[7] );                            \
    return c;                                                   \
  }

  CMATH_FR1(acos)     CMATH_FR1(asin)  CMATH_FR1(atan) CMATH_FR2(atan2)
  CMATH_FR1(ceil)     CMATH_FR1(cos)   CMATH_FR1(cosh) CMATH_FR1(exp)
  CMATH_FR1(fabs)     CMATH_FR1(floor) CMATH_FR2(fmod) CMATH_FR1(log)
  CMATH_FR1(log10)    CMATH_FR2(pow)   CMATH_FR1(sin)  CMATH_FR1(sinh)
  CMATH_FR1(sqrt)     CMATH_FR1(tan)   CMATH_FR1(tanh)

  inline v8float copysign( const v8float &a, const v8float &b )
  {
    v8float c;
    float t;

    t = ::fabs( a.f[0] );
    if( b.f[0] < 0 ) t = -t;
    c.f[0] = t;

    t = ::fabs( a.f[1] );
    if( b.f[1] < 0 ) t = -t;
    c.f[1] = t;

    t = ::fabs( a.f[2] );
    if( b.f[2] < 0 ) t = -t;
    c.f[2] = t;

    t = ::fabs( a.f[3] );
    if( b.f[3] < 0 ) t = -t;
    c.f[3] = t;

    t = ::fabs( a.f[4] );
    if( b.f[4] < 0 ) t = -t;
    c.f[4] = t;

    t = ::fabs( a.f[5] );
    if( b.f[5] < 0 ) t = -t;
    c.f[5] = t;

    t = ::fabs( a.f[6] );
    if( b.f[6] < 0 ) t = -t;
    c.f[6] = t;

    t = ::fabs( a.f[7] );
    if( b.f[7] < 0 ) t = -t;
    c.f[7] = t;

    return c;
  }

# undef CMATH_FR1
# undef CMATH_FR2

  // v8float miscellaneous functions

  inline v8float rsqrt_approx( const v8float &a )
  {
    v8float b;

    b.f[0] = ::sqrt( 1.0f / a.f[0] );
    b.f[1] = ::sqrt( 1.0f / a.f[1] );
    b.f[2] = ::sqrt( 1.0f / a.f[2] );
    b.f[3] = ::sqrt( 1.0f / a.f[3] );
    b.f[4] = ::sqrt( 1.0f / a.f[4] );
    b.f[5] = ::sqrt( 1.0f / a.f[5] );
    b.f[6] = ::sqrt( 1.0f / a.f[6] );
    b.f[7] = ::sqrt( 1.0f / a.f[7] );

    return b;
  }

  inline v8float rsqrt( const v8float &a )
  {
    v8float b;

    b.f[0] = ::sqrt( 1.0f / a.f[0] );
    b.f[1] = ::sqrt( 1.0f / a.f[1] );
    b.f[2] = ::sqrt( 1.0f / a.f[2] );
    b.f[3] = ::sqrt( 1.0f / a.f[3] );
    b.f[4] = ::sqrt( 1.0f / a.f[4] );
    b.f[5] = ::sqrt( 1.0f / a.f[5] );
    b.f[6] = ::sqrt( 1.0f / a.f[6] );
    b.f[7] = ::sqrt( 1.0f / a.f[7] );

    return b;
  }

  inline v8float rcp_approx( const v8float &a )
  {
    v8float b;

    b.f[0] = 1.0f / a.f[0];
    b.f[1] = 1.0f / a.f[1];
    b.f[2] = 1.0f / a.f[2];
    b.f[3] = 1.0f / a.f[3];
    b.f[4] = 1.0f / a.f[4];
    b.f[5] = 1.0f / a.f[5];
    b.f[6] = 1.0f / a.f[6];
    b.f[7] = 1.0f / a.f[7];

    return b;
  }

  inline v8float rcp( const v8float &a )
  {
    v8float b;

    b.f[0] = 1.0f / a.f[0];
    b.f[1] = 1.0f / a.f[1];
    b.f[2] = 1.0f / a.f[2];
    b.f[3] = 1.0f / a.f[3];
    b.f[4] = 1.0f / a.f[4];
    b.f[5] = 1.0f / a.f[5];
    b.f[6] = 1.0f / a.f[6];
    b.f[7] = 1.0f / a.f[7];

    return b;
  }

  inline v8float fma( const v8float &a, const v8float &b, const v8float &c )
  {
    v8float d;

    d.f[0] = a.f[0] * b.f[0] + c.f[0];
    d.f[1] = a.f[1] * b.f[1] + c.f[1];
    d.f[2] = a.f[2] * b.f[2] + c.f[2];
    d.f[3] = a.f[3] * b.f[3] + c.f[3];
    d.f[4] = a.f[4] * b.f[4] + c.f[4];
    d.f[5] = a.f[5] * b.f[5] + c.f[5];
    d.f[6] = a.f[6] * b.f[6] + c.f[6];
    d.f[7] = a.f[7] * b.f[7] + c.f[7];

    return d;
  }

  inline v8float fms( const v8float &a, const v8float &b, const v8float &c )
  {
    v8float d;

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

  inline v8float fnms( const v8float &a, const v8float &b, const v8float &c )
  {
    v8float d;

    d.f[0] = c.f[0] - a.f[0] * b.f[0];
    d.f[1] = c.f[1] - a.f[1] * b.f[1];
    d.f[2] = c.f[2] - a.f[2] * b.f[2];
    d.f[3] = c.f[3] - a.f[3] * b.f[3];
    d.f[4] = c.f[4] - a.f[4] * b.f[4];
    d.f[5] = c.f[5] - a.f[5] * b.f[5];
    d.f[6] = c.f[6] - a.f[6] * b.f[6];
    d.f[7] = c.f[7] - a.f[7] * b.f[7];

    return d;
  }

  inline v8float clear_bits( const v8int &m, const v8float &a )
  {
    v8float b;

    b.i[0] = ( ~m.i[0] ) & a.i[0];
    b.i[1] = ( ~m.i[1] ) & a.i[1];
    b.i[2] = ( ~m.i[2] ) & a.i[2];
    b.i[3] = ( ~m.i[3] ) & a.i[3];
    b.i[4] = ( ~m.i[4] ) & a.i[4];
    b.i[5] = ( ~m.i[5] ) & a.i[5];
    b.i[6] = ( ~m.i[6] ) & a.i[6];
    b.i[7] = ( ~m.i[7] ) & a.i[7];

    return b;
  }

  inline v8float set_bits( const v8int &m, const v8float &a )
  {
    v8float b;

    b.i[0] = m.i[0] | a.i[0];
    b.i[1] = m.i[1] | a.i[1];
    b.i[2] = m.i[2] | a.i[2];
    b.i[3] = m.i[3] | a.i[3];
    b.i[4] = m.i[4] | a.i[4];
    b.i[5] = m.i[5] | a.i[5];
    b.i[6] = m.i[6] | a.i[6];
    b.i[7] = m.i[7] | a.i[7];

    return b;
  }

  inline v8float toggle_bits( const v8int &m, const v8float &a )
  {
    v8float b;

    b.i[0] = m.i[0] ^ a.i[0];
    b.i[1] = m.i[1] ^ a.i[1];
    b.i[2] = m.i[2] ^ a.i[2];
    b.i[3] = m.i[3] ^ a.i[3];
    b.i[4] = m.i[4] ^ a.i[4];
    b.i[5] = m.i[5] ^ a.i[5];
    b.i[6] = m.i[6] ^ a.i[6];
    b.i[7] = m.i[7] ^ a.i[7];

    return b;
  }

  inline void increment_8x1( float * ALIGNED(16) p, const v8float &a )
  {
    p[0] += a.f[0];
    p[1] += a.f[1];
    p[2] += a.f[2];
    p[3] += a.f[3];
    p[4] += a.f[4];
    p[5] += a.f[5];
    p[6] += a.f[6];
    p[7] += a.f[7];
  }

  inline void decrement_8x1( float * ALIGNED(16) p, const v8float &a )
  {
    p[0] -= a.f[0];
    p[1] -= a.f[1];
    p[2] -= a.f[2];
    p[3] -= a.f[3];
    p[4] -= a.f[4];
    p[5] -= a.f[5];
    p[6] -= a.f[6];
    p[7] -= a.f[7];
  }

  inline void scale_8x1( float * ALIGNED(16) p, const v8float &a )
  {
    p[0] *= a.f[0];
    p[1] *= a.f[1];
    p[2] *= a.f[2];
    p[3] *= a.f[3];
    p[4] *= a.f[4];
    p[5] *= a.f[5];
    p[6] *= a.f[6];
    p[7] *= a.f[7];
  }

} // namespace v8

#endif // _v8_portable_h_
