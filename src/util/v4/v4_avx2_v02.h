#ifndef _v4_avx2_h_
#define _v4_avx2_h_

#ifndef IN_v4_h
#error "Do not include v4_avx2.h directly; use v4.h"
#endif

#define V4_ACCELERATION
#define V4_AVX2_ACCELERATION

#include <immintrin.h>
#include <math.h>

#ifndef ALIGNED
#define ALIGNED(n)
#endif

namespace v4
{
  class v4;
  class v4int;
  class v4float;

  ////////////////
  // v4 base class

  class v4
  {
    friend class v4int;
    friend class v4float;

    // v4 miscellenous friends

    friend inline int any( const v4 &a );
    friend inline int all( const v4 &a );

    template<int n>
    friend inline v4 splat( const v4 &a );

    template<int i0, int i1, int i2, int i3>
    friend inline v4 shuffle( const v4 &a );

    friend inline void swap( v4 &a, v4 &b );
    friend inline void transpose( v4 &a0, v4 &a1, v4 &a2, v4 &a3 );

    // v4int miscellaneous friends

    friend inline v4    czero( const v4int &c, const v4 &a );
    friend inline v4 notczero( const v4int &c, const v4 &a );
    friend inline v4 merge( const v4int &c, const v4 &a, const v4 &b );

    // v4 memory manipulation friends

    friend inline void   load_4x1( const void * ALIGNED(16) p, v4 &a );
    friend inline void  store_4x1( const v4 &a, void * ALIGNED(16) p );
    friend inline void stream_4x1( const v4 &a, void * ALIGNED(16) p );
    friend inline void   copy_4x1( void * ALIGNED(16) dst,
                                   const void * ALIGNED(16) src );
    friend inline void   swap_4x1( void * ALIGNED(16) a, void * ALIGNED(16) b );

    // v4 transposed memory manipulation friends
    // Note: Half aligned values are permissible in the 4x2_tr variants!

    friend inline void load_4x1_tr( const void *a0, const void *a1,
                                    const void *a2, const void *a3,
                                    v4 &a );
    friend inline void load_4x2_tr( const void * ALIGNED(8) a0,
                                    const void * ALIGNED(8) a1,
                                    const void * ALIGNED(8) a2,
                                    const void * ALIGNED(8) a3,
                                    v4 &a, v4 &b );
    friend inline void load_4x3_tr( const void * ALIGNED(16) a0,
                                    const void * ALIGNED(16) a1,
                                    const void * ALIGNED(16) a2,
                                    const void * ALIGNED(16) a3,
                                    v4 &a, v4 &b, v4 &c );
    friend inline void load_4x4_tr( const void * ALIGNED(16) a0,
                                    const void * ALIGNED(16) a1,
                                    const void * ALIGNED(16) a2,
                                    const void * ALIGNED(16) a3,
                                    v4 &a, v4 &b, v4 &c, v4 &d );

    friend inline void store_4x1_tr( const v4 &a,
                                     void *a0, void *a1, void *a2, void *a3 );
    friend inline void store_4x2_tr( const v4 &a, const v4 &b,
                                     void * ALIGNED(8) a0,
                                     void * ALIGNED(8) a1,
                                     void * ALIGNED(8) a2,
                                     void * ALIGNED(8) a3 );
    friend inline void store_4x3_tr( const v4 &a, const v4 &b, const v4 &c,
                                     void * ALIGNED(16) a0,
                                     void * ALIGNED(16) a1,
                                     void * ALIGNED(16) a2,
                                     void * ALIGNED(16) a3 );
    friend inline void store_4x4_tr( const v4 &a, const v4 &b,
                                     const v4 &c, const v4 &d,
                                     void * ALIGNED(16) a0,
                                     void * ALIGNED(16) a1,
                                     void * ALIGNED(16) a2,
                                     void * ALIGNED(16) a3 );

  protected:

    union
    {
      int i[4];
      float f[4];
      __m128 v;
    };

  public:

    v4() {}                    // Default constructor

    v4( const v4 &a )          // Copy constructor
    {
      v=a.v;
    }

    ~v4() {}                   // Default destructor
  };

  // v4 miscellaneous functions

  inline int any( const v4 &a )
  {
    return a.i[0] || a.i[1] || a.i[2] || a.i[3];
  }

  inline int all( const v4 &a )
  {
    return a.i[0] && a.i[1] && a.i[2] && a.i[3];
  }

  template<int n>
  inline v4 splat( const v4 & a )
  {
    v4 b;
    __m128 a_v = a.v;

    b.v = _mm_shuffle_ps( a_v, a_v, ( n*permute<1,1,1,1>::value ) );

    return b;
  }

  template<int i0, int i1, int i2, int i3>
  inline v4 shuffle( const v4 & a )
  {
    v4 b;
    __m128 a_v = a.v;

    b.v = _mm_shuffle_ps( a_v, a_v, ( permute<i0,i1,i2,i3>::value ) );

    return b;
  }

# define sw(x,y) x^=y, y^=x, x^=y

  inline void swap( v4 &a, v4 &b )
  { 
    __m128 a_v = a.v;

    a.v = b.v;

    b.v = a_v;
  }

  inline void transpose( v4 &a0, v4 &a1, v4 &a2, v4 &a3 )
  {
    sw( a0.i[1],a1.i[0] ); sw( a0.i[2],a2.i[0] ); sw( a0.i[3],a3.i[0] );
                           sw( a1.i[2],a2.i[1] ); sw( a1.i[3],a3.i[1] );
                                                  sw( a2.i[3],a3.i[2] );
  }

# undef sw

  // v4 memory manipulation functions
  
  inline void load_4x1( const void * ALIGNED(16) p, v4 &a )
  {
    a.i[0] = ((const int * ALIGNED(16))p)[0];
    a.i[1] = ((const int * ALIGNED(16))p)[1];
    a.i[2] = ((const int * ALIGNED(16))p)[2];
    a.i[3] = ((const int * ALIGNED(16))p)[3];
  }

  inline void store_4x1( const v4 &a, void * ALIGNED(16) p )
  {
    ((int * ALIGNED(16))p)[0] = a.i[0];
    ((int * ALIGNED(16))p)[1] = a.i[1];
    ((int * ALIGNED(16))p)[2] = a.i[2];
    ((int * ALIGNED(16))p)[3] = a.i[3];
  }

  inline void stream_4x1( const v4 &a, void * ALIGNED(16) p )
  {
    ((int * ALIGNED(16))p)[0] = a.i[0];
    ((int * ALIGNED(16))p)[1] = a.i[1];
    ((int * ALIGNED(16))p)[2] = a.i[2];
    ((int * ALIGNED(16))p)[3] = a.i[3];
  }

  inline void clear_4x1( void * ALIGNED(16) p )
  {
    ((int * ALIGNED(16))p)[0] = 0;
    ((int * ALIGNED(16))p)[1] = 0;
    ((int * ALIGNED(16))p)[2] = 0;
    ((int * ALIGNED(16))p)[3] = 0;
  }

  // FIXME: Ordering semantics
  inline void copy_4x1( void * ALIGNED(16) dst,
                        const void * ALIGNED(16) src )
  {
    ((int * ALIGNED(16))dst)[0] = ((const int * ALIGNED(16))src)[0];
    ((int * ALIGNED(16))dst)[1] = ((const int * ALIGNED(16))src)[1];
    ((int * ALIGNED(16))dst)[2] = ((const int * ALIGNED(16))src)[2];
    ((int * ALIGNED(16))dst)[3] = ((const int * ALIGNED(16))src)[3];
  }

  inline void swap_4x1( void * ALIGNED(16) a, void * ALIGNED(16) b )
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
  }

  // v4 transposed memory manipulation functions

  inline void load_4x1_tr( const void *a0, const void *a1,
                           const void *a2, const void *a3,
			   v4 &a )
  {
    a.i[0] = ((const int *)a0)[0];
    a.i[1] = ((const int *)a1)[0];
    a.i[2] = ((const int *)a2)[0];
    a.i[3] = ((const int *)a3)[0];
  }

  inline void load_4x2_tr( const void * ALIGNED(8) a0,
                           const void * ALIGNED(8) a1,
                           const void * ALIGNED(8) a2,
                           const void * ALIGNED(8) a3,
                           v4 &a, v4 &b )
  {
    a.i[0] = ((const int * ALIGNED(8))a0)[0];
    b.i[0] = ((const int * ALIGNED(8))a0)[1];

    a.i[1] = ((const int * ALIGNED(8))a1)[0];
    b.i[1] = ((const int * ALIGNED(8))a1)[1];

    a.i[2] = ((const int * ALIGNED(8))a2)[0];
    b.i[2] = ((const int * ALIGNED(8))a2)[1];

    a.i[3] = ((const int * ALIGNED(8))a3)[0];
    b.i[3] = ((const int * ALIGNED(8))a3)[1];
  }

  inline void load_4x3_tr( const void * ALIGNED(16) a0,
                           const void * ALIGNED(16) a1,
                           const void * ALIGNED(16) a2,
                           const void * ALIGNED(16) a3,
                           v4 &a, v4 &b, v4 &c )
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
  }

  inline void load_4x4_tr( const void * ALIGNED(16) a0,
                           const void * ALIGNED(16) a1,
                           const void * ALIGNED(16) a2,
                           const void * ALIGNED(16) a3,
                           v4 &a, v4 &b, v4 &c, v4 &d )
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
  }

  inline void store_4x1_tr( const v4 &a,
                            void *a0, void *a1, void *a2, void *a3 )
  {
    ((int *)a0)[0] = a.i[0];
    ((int *)a1)[0] = a.i[1];
    ((int *)a2)[0] = a.i[2];
    ((int *)a3)[0] = a.i[3];
  }

  inline void store_4x2_tr( const v4 &a, const v4 &b,
                            void * ALIGNED(8) a0, void * ALIGNED(8) a1,
                            void * ALIGNED(8) a2, void * ALIGNED(8) a3 )
  {
    ((int * ALIGNED(8))a0)[0] = a.i[0];
    ((int * ALIGNED(8))a0)[1] = b.i[0];

    ((int * ALIGNED(8))a1)[0] = a.i[1];
    ((int * ALIGNED(8))a1)[1] = b.i[1];

    ((int * ALIGNED(8))a2)[0] = a.i[2];
    ((int * ALIGNED(8))a2)[1] = b.i[2];

    ((int * ALIGNED(8))a3)[0] = a.i[3];
    ((int * ALIGNED(8))a3)[1] = b.i[3];
  }

  inline void store_4x3_tr( const v4 &a, const v4 &b, const v4 &c,
                            void * ALIGNED(16) a0, void * ALIGNED(16) a1,
                            void * ALIGNED(16) a2, void * ALIGNED(16) a3 )
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
  }

  inline void store_4x4_tr( const v4 &a, const v4 &b, const v4 &c, const v4 &d,
                            void * ALIGNED(16) a0, void * ALIGNED(16) a1,
                            void * ALIGNED(16) a2, void * ALIGNED(16) a3 )
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
  }

  //////////////
  // v4int class

  class v4int : public v4
  {
    // v4int prefix unary operator friends

    friend inline v4int operator  +( const v4int & a );
    friend inline v4int operator  -( const v4int & a );
    friend inline v4int operator  ~( const v4int & a );
    friend inline v4int operator  !( const v4int & a );
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector

    // v4int prefix increment / decrement operator friends

    friend inline v4int operator ++( v4int & a );
    friend inline v4int operator --( v4int & a );

    // v4int postfix increment / decrement operator friends

    friend inline v4int operator ++( v4int & a, int );
    friend inline v4int operator --( v4int & a, int );

    // v4int binary operator friends

    friend inline v4int operator  +( const v4int &a, const v4int &b );
    friend inline v4int operator  -( const v4int &a, const v4int &b );
    friend inline v4int operator  *( const v4int &a, const v4int &b );
    friend inline v4int operator  /( const v4int &a, const v4int &b );
    friend inline v4int operator  %( const v4int &a, const v4int &b );
    friend inline v4int operator  ^( const v4int &a, const v4int &b );
    friend inline v4int operator  &( const v4int &a, const v4int &b );
    friend inline v4int operator  |( const v4int &a, const v4int &b );
    friend inline v4int operator <<( const v4int &a, const v4int &b );
    friend inline v4int operator >>( const v4int &a, const v4int &b );

    // v4int logical operator friends

    friend inline v4int operator  <( const v4int &a, const v4int &b );
    friend inline v4int operator  >( const v4int &a, const v4int &b );
    friend inline v4int operator ==( const v4int &a, const v4int &b );
    friend inline v4int operator !=( const v4int &a, const v4int &b );
    friend inline v4int operator <=( const v4int &a, const v4int &b );
    friend inline v4int operator >=( const v4int &a, const v4int &b );
    friend inline v4int operator &&( const v4int &a, const v4int &b );
    friend inline v4int operator ||( const v4int &a, const v4int &b );

    // v4int miscellaneous friends

    friend inline v4int abs( const v4int &a );
    friend inline v4    czero( const v4int &c, const v4 &a );
    friend inline v4 notczero( const v4int &c, const v4 &a );
    // FIXME: cswap, notcswap!
    friend inline v4 merge( const v4int &c, const v4 &t, const v4 &f );

    // v4float unary operator friends

    friend inline v4int operator  !( const v4float & a ); 

    // v4float logical operator friends

    friend inline v4int operator  <( const v4float &a, const v4float &b );
    friend inline v4int operator  >( const v4float &a, const v4float &b );
    friend inline v4int operator ==( const v4float &a, const v4float &b );
    friend inline v4int operator !=( const v4float &a, const v4float &b );
    friend inline v4int operator <=( const v4float &a, const v4float &b );
    friend inline v4int operator >=( const v4float &a, const v4float &b );
    friend inline v4int operator &&( const v4float &a, const v4float &b );
    friend inline v4int operator ||( const v4float &a, const v4float &b );

    // v4float miscellaneous friends

    friend inline v4float clear_bits(  const v4int &m, const v4float &a );
    friend inline v4float set_bits(    const v4int &m, const v4float &a );
    friend inline v4float toggle_bits( const v4int &m, const v4float &a );

  public:

    // v4int constructors / destructors

    v4int() {}                                // Default constructor

    v4int( const v4int &a )                   // Copy constructor
    {
      v = a.v;
    }

    v4int( const v4 &a )                      // Init from mixed
    {
      v = a.v;
    }

    v4int( int a )                            // Init from scalar
    {
      union
      {
	int i;
	float f;
      } u;

      u.i = a;
      v   = _mm_set1_ps( u.f );
    }

    v4int( int i0, int i1, int i2, int i3 )   // Init from scalars
    {
      union
      {
	int i;
	float f;
      } u0, u1, u2, u3;

      u0.i = i0;
      u1.i = i1;
      u2.i = i2;
      u3.i = i3;

      v = _mm_setr_ps( u0.f, u1.f, u2.f, u3.f );
    }

    ~v4int() {}                               // Destructor

    // v4int assignment operators
  
#   define ASSIGN(op)			          \
    inline v4int &operator op( const v4int &b )   \
    {						  \
      i[0] op b.i[0];                             \
      i[1] op b.i[1];                             \
      i[2] op b.i[2];                             \
      i[3] op b.i[3];                             \
      return *this;                               \
    }

    inline v4int &operator =( const v4int &b )
    {
      v = b.v;

      return *this;
    }

    ASSIGN(+=)
    ASSIGN(-=)
    ASSIGN(*=)
    ASSIGN(/=)
    ASSIGN(%=)

    inline v4int &operator ^=( const v4int &b )
    {
      v = _mm_xor_ps( v, b.v );

      return *this;
    }

    inline v4int &operator &=( const v4int &b )
    {
      v = _mm_and_ps( v, b.v );

      return *this;
    }

    inline v4int &operator |=( const v4int &b )
    {
      v = _mm_or_ps( v, b.v );

      return *this;
    }

    ASSIGN(<<=)
    ASSIGN(>>=)

#   undef ASSIGN

    // v4int member access operator

    inline int &operator []( int n )
    {
      return i[n];
    }

    inline int  operator ()( int n )
    {
      return i[n];
    }
  };

  // v4int prefix unary operators

# define PREFIX_UNARY(op)                       \
  inline v4int operator op( const v4int & a )   \
  {						\
    v4int b;                                    \
    b.i[0] = ( op a.i[0] );                     \
    b.i[1] = ( op a.i[1] );                     \
    b.i[2] = ( op a.i[2] );                     \
    b.i[3] = ( op a.i[3] );                     \
    return b;                                   \
  }

  inline v4int operator +( const v4int & a )
  {
    v4int b;

    b.v = a.v;

    return b;
  }

  PREFIX_UNARY(-)

  inline v4int operator !( const v4int & a )
  {
    v4int b;

    b.i[0] = - ( !a.i[0] );
    b.i[1] = - ( !a.i[1] );
    b.i[2] = - ( !a.i[2] );
    b.i[3] = - ( !a.i[3] );

    return b;
  }

  inline v4int operator ~( const v4int & a )
  {
    v4int b;

    union
    {
      int i;
      float f;
    } u;

    u.i = -1;
    b.v = _mm_xor_ps( a.v, _mm_set1_ps( u.f ) );

    return b;
  }

# undef PREFIX_UNARY

  // v4int prefix increment / decrement

# define PREFIX_INCDEC(op)                      \
  inline v4int operator op( v4int & a )         \
  {						\
    v4int b;                                    \
    b.i[0] = ( op a.i[0] );                     \
    b.i[1] = ( op a.i[1] );                     \
    b.i[2] = ( op a.i[2] );                     \
    b.i[3] = ( op a.i[3] );                     \
    return b;                                   \
  }

  PREFIX_INCDEC(++)
  PREFIX_INCDEC(--)

# undef PREFIX_INCDEC

  // v4int postfix increment / decrement

# define POSTFIX_INCDEC(op)                    \
  inline v4int operator op( v4int & a, int )   \
  {					       \
    v4int b;                                   \
    b.i[0] = ( a.i[0] op );                    \
    b.i[1] = ( a.i[1] op );                    \
    b.i[2] = ( a.i[2] op );                    \
    b.i[3] = ( a.i[3] op );                    \
    return b;                                  \
  }

  POSTFIX_INCDEC(++)
  POSTFIX_INCDEC(--)

# undef POSTFIX_INCDEC

  // v4int binary operators
  
# define BINARY(op)                                             \
  inline v4int operator op( const v4int &a, const v4int &b )    \
  {								\
    v4int c;                                                    \
    c.i[0] = a.i[0] op b.i[0];                                  \
    c.i[1] = a.i[1] op b.i[1];                                  \
    c.i[2] = a.i[2] op b.i[2];                                  \
    c.i[3] = a.i[3] op b.i[3];                                  \
    return c;                                                   \
  }

  BINARY(+)
  BINARY(-)
  BINARY(*)
  BINARY(/)
  BINARY(%)

  inline v4int operator ^( const v4int &a, const v4int &b )
  {
    v4int c;

    c.v = _mm_xor_ps( a.v, b.v );

    return c;
  }

  inline v4int operator &( const v4int &a, const v4int &b )
  {
    v4int c;

    c.v = _mm_and_ps( a.v, b.v );

    return c;
  }

  inline v4int operator |( const v4int &a, const v4int &b )
  {
    v4int c;

    c.v = _mm_or_ps( a.v, b.v );

    return c;
  }

  BINARY(<<)
  BINARY(>>)

# undef BINARY

  // v4int logical operators

# define LOGICAL(op)                                           \
  inline v4int operator op( const v4int &a, const v4int &b )   \
  {							       \
    v4int c;                                                   \
    c.i[0] = -(a.i[0] op b.i[0]);                              \
    c.i[1] = -(a.i[1] op b.i[1]);                              \
    c.i[2] = -(a.i[2] op b.i[2]);                              \
    c.i[3] = -(a.i[3] op b.i[3]);                              \
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

  // v4int miscellaneous functions

  inline v4int abs( const v4int &a )
  {
    v4int b;

    b.i[0] = ( a.i[0] >= 0 ) ? a.i[0] : -a.i[0];
    b.i[1] = ( a.i[1] >= 0 ) ? a.i[1] : -a.i[1];
    b.i[2] = ( a.i[2] >= 0 ) ? a.i[2] : -a.i[2];
    b.i[3] = ( a.i[3] >= 0 ) ? a.i[3] : -a.i[3];

    return b;
  }

  inline v4 czero( const v4int &c, const v4 &a )
  {
    v4 b;

    b.v = _mm_andnot_ps( c.v, a.v );

    return b;
  }

  inline v4 notczero( const v4int &c, const v4 &a )
  {
    v4 b;

    b.v = _mm_and_ps( c.v, a.v );

    return b;
  }

  inline v4 merge( const v4int &c, const v4 &t, const v4 &f )
  {
    __m128 c_v = c.v;
    v4 tf;

    tf.v = _mm_or_ps( _mm_andnot_ps( c_v, f.v ),
		      _mm_and_ps( c_v, t.v ) );

    return tf;
  }

  ////////////////
  // v4float class

  class v4float : public v4
  {
    // v4float prefix unary operator friends

    friend inline v4float operator  +( const v4float &a );
    friend inline v4float operator  -( const v4float &a );
    friend inline v4float operator  ~( const v4float &a );
    friend inline v4int   operator  !( const v4float &a );
    // Note: Referencing (*) and dereferencing (&) apply to the whole vector

    // v4float prefix increment / decrement operator friends

    friend inline v4float operator ++( v4float &a );
    friend inline v4float operator --( v4float &a );

    // v4float postfix increment / decrement operator friends

    friend inline v4float operator ++( v4float &a, int );
    friend inline v4float operator --( v4float &a, int );

    // v4float binary operator friends

    friend inline v4float operator  +( const v4float &a, const v4float &b );
    friend inline v4float operator  -( const v4float &a, const v4float &b );
    friend inline v4float operator  *( const v4float &a, const v4float &b );
    friend inline v4float operator  /( const v4float &a, const v4float &b );

    // v4float logical operator friends

    friend inline v4int operator  <( const v4float &a, const v4float &b );
    friend inline v4int operator  >( const v4float &a, const v4float &b );
    friend inline v4int operator ==( const v4float &a, const v4float &b );
    friend inline v4int operator !=( const v4float &a, const v4float &b );
    friend inline v4int operator <=( const v4float &a, const v4float &b );
    friend inline v4int operator >=( const v4float &a, const v4float &b );
    friend inline v4int operator &&( const v4float &a, const v4float &b );
    friend inline v4int operator ||( const v4float &a, const v4float &b );

    // v4float math library friends

#   define CMATH_FR1(fn) friend inline v4float fn( const v4float &a )
#   define CMATH_FR2(fn) friend inline v4float fn( const v4float &a,  \
                                                   const v4float &b )

    CMATH_FR1(acos);  CMATH_FR1(asin);  CMATH_FR1(atan); CMATH_FR2(atan2);
    CMATH_FR1(ceil);  CMATH_FR1(cos);   CMATH_FR1(cosh); CMATH_FR1(exp);
    CMATH_FR1(fabs);  CMATH_FR1(floor); CMATH_FR2(fmod); CMATH_FR1(log);
    CMATH_FR1(log10); CMATH_FR2(pow);   CMATH_FR1(sin);  CMATH_FR1(sinh);
    CMATH_FR1(sqrt);  CMATH_FR1(tan);   CMATH_FR1(tanh);

    CMATH_FR2(copysign);

#   undef CMATH_FR1
#   undef CMATH_FR2

    // v4float miscellaneous friends

    friend inline v4float rsqrt_approx( const v4float &a );
    friend inline v4float rsqrt       ( const v4float &a );
    friend inline v4float rcp_approx( const v4float &a );
    friend inline v4float rcp       ( const v4float &a );
    friend inline v4float fma ( const v4float &a, const v4float &b, const v4float &c );
    friend inline v4float fms ( const v4float &a, const v4float &b, const v4float &c );
    friend inline v4float fnms( const v4float &a, const v4float &b, const v4float &c );
    friend inline v4float  clear_bits( const v4int &m, const v4float &a );
    friend inline v4float    set_bits( const v4int &m, const v4float &a );
    friend inline v4float toggle_bits( const v4int &m, const v4float &a );
    friend inline void increment_4x1( float * ALIGNED(16) p, const v4float &a );
    friend inline void decrement_4x1( float * ALIGNED(16) p, const v4float &a );
    friend inline void     scale_4x1( float * ALIGNED(16) p, const v4float &a );
    // FIXME: crack
    friend inline void trilinear( v4float & wl, v4float & wh );

  public:

    // v4float constructors / destructors

    v4float() {}                                        // Default constructor

    v4float( const v4float &a )                         // Copy constructor
    {
      v = a.v;
    }

    v4float( const v4 &a )                              // Init from mixed
    {
      v = a.v;
    }

    v4float( float a )                                  // Init from scalar
    {
      v = _mm_set1_ps( a );
    }

    v4float( float f0, float f1, float f2, float f3 )   // Init from scalars
    {
      v = _mm_setr_ps( f0, f1, f2, f3 );
    }

    ~v4float() {}                                       // Destructor

    // v4float assignment operators

#   define ASSIGN(op,intrin)				\
    inline v4float &operator op( const v4float &b )     \
    {							\
      v = intrin( v, b.v );                             \
      return *this;					\
    }

    inline v4float &operator =( const v4float &b )
    {
      v = b.v;

      return *this;
    }

    ASSIGN( +=, _mm_add_ps )
    ASSIGN( -=, _mm_sub_ps )
    ASSIGN( *=, _mm_mul_ps )
    ASSIGN( /=, _mm_div_ps )

#   undef ASSIGN

    // v4float member access operator

    inline float &operator []( int n )
    {
      return f[n];
    }

    inline float  operator ()( int n )
    {
      return f[n];
    }
  };

  // v4float prefix unary operators

  inline v4float operator +( const v4float &a )
  {
    v4float b;

    b.v = a.v;

    return b;
  }

  inline v4float operator -( const v4float &a )
  {
    v4float b;

    b.v = _mm_sub_ps( _mm_setzero_ps(), a.v );

    return b;
  }

  inline v4int operator !( const v4float &a )
  {
    v4int b;

    b.v = _mm_cmpeq_ps( _mm_setzero_ps(), a.v );

    return b;
  }

  // v4float prefix increment / decrement operators

  inline v4float operator ++( v4float &a )
  {
    v4float b;

    __m128 t = _mm_add_ps( a.v, _mm_set1_ps( 1 ) );

    a.v = t;
    b.v = t;

    return b;
  }

  inline v4float operator --( v4float &a )
  {
    v4float b;

    __m128 t = _mm_sub_ps( a.v, _mm_set1_ps( 1 ) );

    a.v = t;
    b.v = t;

    return b;
  }

  // v4float postfix increment / decrement operators

  inline v4float operator ++( v4float &a, int )
  {
    v4float b;

    __m128 a_v = a.v;

    a.v = _mm_add_ps( a_v, _mm_set1_ps( 1 ) );
    b.v = a_v;

    return b;
  }

  inline v4float operator --( v4float &a, int )
  {
    v4float b;

    __m128 a_v = a.v;

    a.v = _mm_sub_ps( a_v, _mm_set1_ps( 1 ) );
    b.v = a_v;

    return b;
  }

  // v4float binary operators

# define BINARY(op,intrin)                                           \
  inline v4float operator op( const v4float &a, const v4float &b )   \
  {								     \
    v4float c;                                                       \
    c.v = intrin( a.v, b.v );                                        \
    return c;                                                        \
  }

  BINARY( +, _mm_add_ps )
  BINARY( -, _mm_sub_ps )
  BINARY( *, _mm_mul_ps )
  BINARY( /, _mm_div_ps )

# undef BINARY

  // v4float logical operators

# define LOGICAL(op,intrin)                                        \
  inline v4int operator op( const v4float &a, const v4float &b )   \
  {								   \
    v4int c;                                                       \
    c.v = intrin( a.v, b.v );                                      \
    return c;                                                      \
  }

  LOGICAL(  <, _mm_cmplt_ps )
  LOGICAL(  >, _mm_cmpgt_ps )
  LOGICAL( ==, _mm_cmpeq_ps )
  LOGICAL( !=, _mm_cmpneq_ps )
  LOGICAL( <=, _mm_cmple_ps )
  LOGICAL( >=, _mm_cmpge_ps )

  inline v4int operator &&( const v4float &a, const v4float &b )
  {
    v4int c;

    __m128 vzero = _mm_setzero_ps();

    c.v = _mm_and_ps( _mm_cmpneq_ps( a.v, vzero ),
		      _mm_cmpneq_ps( b.v, vzero ) );

    return c;
  }

  inline v4int operator ||( const v4float &a, const v4float &b )
  {
    v4int c;

    __m128 vzero = _mm_setzero_ps();

    c.v = _mm_or_ps( _mm_cmpneq_ps( a.v, vzero ),
		     _mm_cmpneq_ps( b.v, vzero ) );

    return c;
  }

# undef LOGICAL

  // v4float math library functions

# define CMATH_FR1(fn)                          \
  inline v4float fn( const v4float &a )         \
  {						\
    v4float b;                                  \
    b.f[0] = ::fn( a.f[0] );                    \
    b.f[1] = ::fn( a.f[1] );                    \
    b.f[2] = ::fn( a.f[2] );                    \
    b.f[3] = ::fn( a.f[3] );                    \
    return b;                                   \
  }

# define CMATH_FR2(fn)                                          \
  inline v4float fn( const v4float &a, const v4float &b )       \
  {								\
    v4float c;                                                  \
    c.f[0] = ::fn( a.f[0], b.f[0] );                            \
    c.f[1] = ::fn( a.f[1], b.f[1] );                            \
    c.f[2] = ::fn( a.f[2], b.f[2] );                            \
    c.f[3] = ::fn( a.f[3], b.f[3] );                            \
    return c;                                                   \
  }

  CMATH_FR1(acos)     CMATH_FR1(asin)  CMATH_FR1(atan) CMATH_FR2(atan2)
  CMATH_FR1(ceil)     CMATH_FR1(cos)   CMATH_FR1(cosh) CMATH_FR1(exp)
  /*CMATH_FR1(fabs)*/ CMATH_FR1(floor) CMATH_FR2(fmod) CMATH_FR1(log)
  CMATH_FR1(log10)    CMATH_FR2(pow)   CMATH_FR1(sin)  CMATH_FR1(sinh)
  /*CMATH_FR1(sqrt)*/ CMATH_FR1(tan)   CMATH_FR1(tanh)

  inline v4float fabs( const v4float &a )
  {
    v4float b;

    b.v = _mm_andnot_ps( _mm_set1_ps( -0.0f ), a.v );

    return b;
  }

  inline v4float sqrt( const v4float &a )
  {
    v4float b;

    b.v = _mm_sqrt_ps( a.v );

    return b;
  }

  inline v4float copysign( const v4float &a, const v4float &b )
  {
    v4float c;

    __m128 t = _mm_set1_ps( -0.0f );

    c.v = _mm_or_ps( _mm_and_ps( t, b.v ),
		     _mm_andnot_ps( t, a.v ) );

    return c;
  }

# undef CMATH_FR1
# undef CMATH_FR2

  // v4float miscelleanous functions

  inline v4float rsqrt_approx( const v4float &a )
  {
    v4float b;

    b.v = _mm_rsqrt_ps( a.v );

    return b;
  }

  inline v4float rsqrt( const v4float &a )
  {
    v4float b;

    __m128 b_v;

    b_v = _mm_rsqrt_ps( a.v );

    // Note: It is quicker to just call div_ps and sqrt_ps if more
    // refinement desired.

    b.v = _mm_fmadd_ps( _mm_set1_ps( 0.5f ),
			_mm_fnmadd_ps( a.v,
				       _mm_mul_ps( b_v,
						   _mm_mul_ps( b_v, b_v ) ),
				       b_v ),
			b_v );

    return b;
  }

  inline v4float rcp_approx( const v4float &a )
  {
    v4float b;

    b.v = _mm_rcp_ps( a.v );

    return b;
  }
  
  inline v4float rcp( const v4float &a )
  {
    v4float b;

    __m128 b_v;

    b_v = _mm_rcp_ps( a.v );

    b.v = _mm_fnmadd_ps( a.v,
			 _mm_mul_ps( b_v, b_v ),
			 _mm_add_ps( b_v, b_v ) );

    return b;
  }

  inline v4float fma( const v4float &a, const v4float &b, const v4float &c )
  {
    v4float d;

    d.v = _mm_fmadd_ps( a.v, b.v, c.v );

    return d;
  }

  inline v4float fms( const v4float &a, const v4float &b, const v4float &c )
  {
    v4float d;

    d.v = _mm_fmsub_ps( a.v, b.v, c.v );

    return d;
  }

  inline v4float fnms( const v4float &a, const v4float &b, const v4float &c )
  {
    v4float d;

    d.v = _mm_fnmadd_ps( a.v, b.v, c.v );

    return d;
  }

  inline v4float clear_bits( const v4int &m, const v4float &a )
  {
    v4float b;

    b.v = _mm_andnot_ps( m.v, a.v );

    return b;
  }

  inline v4float set_bits( const v4int &m, const v4float &a )
  {
    v4float b;

    b.v = _mm_or_ps( m.v, a.v );

    return b;
  }

  inline v4float toggle_bits( const v4int &m, const v4float &a )
  {
    v4float b;

    b.v = _mm_xor_ps( m.v, a.v );

    return b;
  }

  inline void increment_4x1( float * ALIGNED(16) p, const v4float &a )
  {
    _mm_store_ps( p, _mm_add_ps( _mm_load_ps( p ), a.v ) );
  }

  inline void decrement_4x1( float * ALIGNED(16) p, const v4float &a )
  {
    _mm_store_ps( p, _mm_sub_ps( _mm_load_ps( p ), a.v ) );
  }

  inline void scale_4x1( float * ALIGNED(16) p, const v4float &a )
  {
    _mm_store_ps( p, _mm_mul_ps( _mm_load_ps( p ), a.v ) );
  }

  // Given wl = x y z w, compute:
  // wl = (1-x)(1-y)(1-z) (1+x)(1-y)(1-z) (1-x)(1+y)(1-z) (1+x)(1+y)(1-z)
  // wh = (1-x)(1-y)(1+z) (1+x)(1-y)(1+z) (1-x)(1+y)(1+z) (1+x)(1+y)(1+z)
  inline void trilinear( v4float &wl, v4float &wh )
  {
    __m128 l = _mm_set1_ps( 1.0f ), s = _mm_setr_ps( -0.0f, +0.0f, -0.0f, +0.0f );
    __m128 z = wl.v, xy;

    xy = _mm_add_ps( l, _mm_xor_ps( s, _mm_shuffle_ps( z,z, PERM(0,0,1,1) ) ) );

    z  = _mm_add_ps( l, _mm_xor_ps( s, _mm_shuffle_ps( z,z, PERM(2,2,2,2) ) ) );

    xy = _mm_mul_ps( _mm_shuffle_ps( xy,xy, PERM(0,1,0,1) ),
                     _mm_shuffle_ps( xy,xy, PERM(2,2,3,3) ) );

    wl.v = _mm_mul_ps( xy, _mm_shuffle_ps( z,z, PERM(0,0,0,0) ) );

    wh.v = _mm_mul_ps( xy, _mm_shuffle_ps( z,z, PERM(1,1,1,1) ) );
  }

} // namespace v4

#endif // _v4_avx2_h_
