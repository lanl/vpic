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
    friend inline void load_16x8_tr_p( const void * ALIGNED(64) a00,
				       const void * ALIGNED(64) a01,
				       const void * ALIGNED(64) a02,
				       const void * ALIGNED(64) a03,
				       const void * ALIGNED(64) a04,
				       const void * ALIGNED(64) a05,
				       const void * ALIGNED(64) a06,
				       const void * ALIGNED(64) a07,
				       v16 &a, v16 &b, v16 &c, v16 &d,
				       v16 &e, v16 &f, v16 &g, v16 &h );
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
    friend inline void store_16x8_tr_p( const v16 &a, const v16 &b,
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
					void * ALIGNED(64) a07 );
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
      v = a.v;
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

    b.v = _mm512_set1_ps( a.v[n] );

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
    __m512 a_v = a.v;

    a.v = b.v;

    b.v = a_v;
  }

  inline void transpose( v16 &a00, v16 &a01, v16 &a02, v16 &a03,
			 v16 &a04, v16 &a05, v16 &a06, v16 &a07,
			 v16 &a08, v16 &a09, v16 &a10, v16 &a11,
			 v16 &a12, v16 &a13, v16 &a14, v16 &a15 )
  {
    __m512 t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14, t15;

    // Start                                 a00 =   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    //                                       a01 =  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31
    //                                       a02 =  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47
    //                                       a03 =  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63
    //                                       a04 =  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79
    //                                       a05 =  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95
    //                                       a06 =  96  97  98  99 100 101 102 103 104 105 106 107 108 109 110 111
    //                                       a07 = 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127
    //                                       a08 = 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143
    //                                       a09 = 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159
    //                                       a10 = 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175
    //                                       a11 = 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191
    //                                       a12 = 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207
    //                                       a13 = 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223
    //                                       a14 = 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239
    //                                       a15 = 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255

    t00   = _mm512_unpacklo_ps( a00.v, a01.v ); //   0  16   1  17   4  20   5  21   8  24   9  25  12  28  13  29 
    t01   = _mm512_unpackhi_ps( a00.v, a01.v ); //   2  18   3  19   6  22   7  23  10  26  11  27  14  30  15  31
    t02   = _mm512_unpacklo_ps( a02.v, a03.v ); //  32  48  33  49  36  52  37  53  40  56  41  57  44  60  45  61
    t03   = _mm512_unpackhi_ps( a02.v, a03.v ); //  34  50  35  51  38  54  39  55  42  58  43  59  46  62  47  63
    t04   = _mm512_unpacklo_ps( a04.v, a05.v ); //  64  80  65  81  68  84  69  85  72  88  73  89  76  92  77  93
    t05   = _mm512_unpackhi_ps( a04.v, a05.v ); //  66  82  67  83  70  86  71  87  74  90  75  91  78  94  79  95
    t06   = _mm512_unpacklo_ps( a06.v, a07.v ); //  96 112  97 113 100 116 101 117 104 120 105 121 108 124 109 125
    t07   = _mm512_unpackhi_ps( a06.v, a07.v ); //  98 114  99 115 102 118 103 119 106 122 107 123 110 126 111 127
    t08   = _mm512_unpacklo_ps( a08.v, a09.v ); // 128 144 129 145 132 148 133 149 136 152 137 153 140 156 141 157
    t09   = _mm512_unpackhi_ps( a08.v, a09.v ); // 130 146 131 147 134 150 135 151 138 154 139 155 142 158 143 159
    t10   = _mm512_unpacklo_ps( a10.v, a11.v ); // 160 176 161 177 164 180 165 181 168 184 169 185 172 188 173 189
    t11   = _mm512_unpackhi_ps( a10.v, a11.v ); // 162 178 163 179 166 182 167 183 170 186 171 187 174 190 175 191
    t12   = _mm512_unpacklo_ps( a12.v, a13.v ); // 192 208 193 209 196 212 197 213 200 216 201 217 204 220 205 221
    t13   = _mm512_unpackhi_ps( a12.v, a13.v ); // 194 210 195 211 198 214 199 215 202 218 203 219 206 222 207 223
    t14   = _mm512_unpacklo_ps( a14.v, a15.v ); // 224 240 225 241 228 244 229 245 232 248 233 249 236 252 237 253
    t15   = _mm512_unpackhi_ps( a14.v, a15.v ); // 226 242 227 243 230 246 231 247 234 250 235 251 238 254 239 255

    a00.v = _mm512_shuffle_ps( t00, t02, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //   0  16  32  48 
    a01.v = _mm512_shuffle_ps( t00, t02, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //   1  17  33  49 
    a02.v = _mm512_shuffle_ps( t01, t03, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //   2  18  34  50 
    a03.v = _mm512_shuffle_ps( t01, t03, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //   3  19  35  51 
    a04.v = _mm512_shuffle_ps( t04, t06, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //  64  80  96 112 
    a05.v = _mm512_shuffle_ps( t04, t06, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //  65  81  97 113 
    a06.v = _mm512_shuffle_ps( t05, t07, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //  66  82  98 114 
    a07.v = _mm512_shuffle_ps( t05, t07, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //  67  83  99 115 
    a08.v = _mm512_shuffle_ps( t08, t10, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 128 144 160 176 
    a09.v = _mm512_shuffle_ps( t08, t10, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 129 145 161 177 
    a10.v = _mm512_shuffle_ps( t09, t11, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 130 146 162 178 
    a11.v = _mm512_shuffle_ps( t09, t11, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 131 147 163 179 
    a12.v = _mm512_shuffle_ps( t12, t14, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 192 208 228 240 
    a13.v = _mm512_shuffle_ps( t12, t14, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 193 209 229 241 
    a14.v = _mm512_shuffle_ps( t13, t15, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 194 210 230 242 
    a15.v = _mm512_shuffle_ps( t13, t15, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 195 211 231 243 

    t00   = _mm512_shuffle_f32x4( a00.v, a04.v, 0x88 ); //   0  16  32  48   8  24  40  56  64  80  96  112 ...
    t01   = _mm512_shuffle_f32x4( a01.v, a05.v, 0x88 ); //   1  17  33  49 ...
    t02   = _mm512_shuffle_f32x4( a02.v, a06.v, 0x88 ); //   2  18  34  50 ...
    t03   = _mm512_shuffle_f32x4( a03.v, a07.v, 0x88 ); //   3  19  35  51 ...
    t04   = _mm512_shuffle_f32x4( a00.v, a04.v, 0xdd ); //   4  20  36  52 ...
    t05   = _mm512_shuffle_f32x4( a01.v, a05.v, 0xdd ); //   5  21  37  53 ...
    t06   = _mm512_shuffle_f32x4( a02.v, a06.v, 0xdd ); //   6  22  38  54 ...
    t07   = _mm512_shuffle_f32x4( a03.v, a07.v, 0xdd ); //   7  23  39  55 ...
    t08   = _mm512_shuffle_f32x4( a08.v, a12.v, 0x88 ); // 128 144 160 176 ...
    t09   = _mm512_shuffle_f32x4( a09.v, a13.v, 0x88 ); // 129 145 161 177 ...
    t10   = _mm512_shuffle_f32x4( a10.v, a14.v, 0x88 ); // 130 146 162 178 ...
    t11   = _mm512_shuffle_f32x4( a11.v, a15.v, 0x88 ); // 131 147 163 179 ...
    t12   = _mm512_shuffle_f32x4( a08.v, a12.v, 0xdd ); // 132 148 164 180 ...
    t13   = _mm512_shuffle_f32x4( a09.v, a13.v, 0xdd ); // 133 149 165 181 ...
    t14   = _mm512_shuffle_f32x4( a10.v, a14.v, 0xdd ); // 134 150 166 182 ...
    t15   = _mm512_shuffle_f32x4( a11.v, a15.v, 0xdd ); // 135 151 167 183 ...

    a00.v = _mm512_shuffle_f32x4( t00, t08, 0x88 ); //   0  16  32  48  64  80  96 112 ... 240
    a01.v = _mm512_shuffle_f32x4( t01, t09, 0x88 ); //   1  17  33  49  66  81  97 113 ... 241
    a02.v = _mm512_shuffle_f32x4( t02, t10, 0x88 ); //   2  18  34  50  67  82  98 114 ... 242
    a03.v = _mm512_shuffle_f32x4( t03, t11, 0x88 ); //   3  19  35  51  68  83  99 115 ... 243
    a04.v = _mm512_shuffle_f32x4( t04, t12, 0x88 ); //   4 ...
    a05.v = _mm512_shuffle_f32x4( t05, t13, 0x88 ); //   5 ...
    a06.v = _mm512_shuffle_f32x4( t06, t14, 0x88 ); //   6 ...
    a07.v = _mm512_shuffle_f32x4( t07, t15, 0x88 ); //   7 ...
    a08.v = _mm512_shuffle_f32x4( t00, t08, 0xdd ); //   8 ...
    a09.v = _mm512_shuffle_f32x4( t01, t09, 0xdd ); //   9 ...
    a10.v = _mm512_shuffle_f32x4( t02, t10, 0xdd ); //  10 ...
    a11.v = _mm512_shuffle_f32x4( t03, t11, 0xdd ); //  11 ...
    a12.v = _mm512_shuffle_f32x4( t04, t12, 0xdd ); //  12 ...
    a13.v = _mm512_shuffle_f32x4( t05, t13, 0xdd ); //  13 ...
    a14.v = _mm512_shuffle_f32x4( t06, t14, 0xdd ); //  14 ...
    a15.v = _mm512_shuffle_f32x4( t07, t15, 0xdd ); //  15  31  47  63  79  96 111 127 ... 255
  }

#if 0
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
#endif

# undef sw

  // v16 memory manipulation functions

  inline void load_16x1( const void * ALIGNED(64) p,
			 v16 &a )
  {
    for( int j = 0; j < 16; j++ )
      a.i[j] = ((const int * ALIGNED(64))p)[j];
  }

  inline void store_16x1( const v16 &a,
			  void * ALIGNED(64) p )
  {
    for( int j = 0; j < 16; j++ )
      ((int * ALIGNED(64))p)[j] = a.i[j];
  }

  inline void stream_16x1( const v16 &a,
			   void * ALIGNED(64) p )
  {
    for( int j = 0; j < 16; j++ )
      ((int * ALIGNED(64))p)[j] = a.i[j];
  }

  inline void clear_16x1( void * ALIGNED(64) p )
  {
    for( int j = 0; j < 16; j++ )
      ((int * ALIGNED(64))p)[j] = 0;
  }

  // FIXME: Ordering semantics
  inline void copy_16x1( void * ALIGNED(64) dst,
			 const void * ALIGNED(64) src )
  {
    for( int j = 0; j < 16; j++ )
      ((int * ALIGNED(64))dst)[j] = ((const int * ALIGNED(64))src)[j];
  }

  inline void swap_16x1( void * ALIGNED(64) a,
			 void * ALIGNED(64) b )
  {
    int t;

    for( int j = 0; j < 16; j++ )
    {
      t = ((int * ALIGNED(64))a)[j];
      ((int * ALIGNED(64))a)[j] = ((int * ALIGNED(64))b)[j];
      ((int * ALIGNED(64))b)[j] = t;
    }
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
    __m512 t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14, t15;

    b00.v = _mm512_load_ps( (const float *)a00 );                      //   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    b01.v = _mm512_load_ps( (const float *)a01 );                      //  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31
    b02.v = _mm512_load_ps( (const float *)a02 );                      //  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47
    b03.v = _mm512_load_ps( (const float *)a03 );                      //  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63
    b04.v = _mm512_load_ps( (const float *)a04 );                      //  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79
    b05.v = _mm512_load_ps( (const float *)a05 );                      //  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95
    b06.v = _mm512_load_ps( (const float *)a06 );                      //  96  97  98  99 100 101 102 103 104 105 106 107 108 109 110 111
    b07.v = _mm512_load_ps( (const float *)a07 );                      // 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127
    b08.v = _mm512_load_ps( (const float *)a08 );                      // 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143
    b09.v = _mm512_load_ps( (const float *)a09 );                      // 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159
    b10.v = _mm512_load_ps( (const float *)a10 );                      // 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175
    b11.v = _mm512_load_ps( (const float *)a11 );                      // 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191
    b12.v = _mm512_load_ps( (const float *)a12 );                      // 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207
    b13.v = _mm512_load_ps( (const float *)a13 );                      // 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223
    b14.v = _mm512_load_ps( (const float *)a14 );                      // 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239
    b15.v = _mm512_load_ps( (const float *)a15 );                      // 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255

    t00   = _mm512_unpacklo_ps( b00.v, b01.v );                        //   0  16   1  17   4  20   5  21   8  24   9  25  12  28  13  29 
    t01   = _mm512_unpackhi_ps( b00.v, b01.v );                        //   2  18   3  19   6  22   7  23  10  26  11  27  14  30  15  31
    t02   = _mm512_unpacklo_ps( b02.v, b03.v );                        //  32  48  33  49  36  52  37  53  40  56  41  57  44  60  45  61
    t03   = _mm512_unpackhi_ps( b02.v, b03.v );                        //  34  50  35  51  38  54  39  55  42  58  43  59  46  62  47  63
    t04   = _mm512_unpacklo_ps( b04.v, b05.v );                        //  64  80  65  81  68  84  69  85  72  88  73  89  76  92  77  93
    t05   = _mm512_unpackhi_ps( b04.v, b05.v );                        //  66  82  67  83  70  86  71  87  74  90  75  91  78  94  79  95
    t06   = _mm512_unpacklo_ps( b06.v, b07.v );                        //  96 112  97 113 100 116 101 117 104 120 105 121 108 124 109 125
    t07   = _mm512_unpackhi_ps( b06.v, b07.v );                        //  98 114  99 115 102 118 103 119 106 122 107 123 110 126 111 127
    t08   = _mm512_unpacklo_ps( b08.v, b09.v );                        // 128 144 129 145 132 148 133 149 136 152 137 153 140 156 141 157
    t09   = _mm512_unpackhi_ps( b08.v, b09.v );                        // 130 146 131 147 134 150 135 151 138 154 139 155 142 158 143 159
    t10   = _mm512_unpacklo_ps( b10.v, b11.v );                        // 160 176 161 177 164 180 165 181 168 184 169 185 172 188 173 189
    t11   = _mm512_unpackhi_ps( b10.v, b11.v );                        // 162 178 163 179 166 182 167 183 170 186 171 187 174 190 175 191
    t12   = _mm512_unpacklo_ps( b12.v, b13.v );                        // 192 208 193 209 196 212 197 213 200 216 201 217 204 220 205 221
    t13   = _mm512_unpackhi_ps( b12.v, b13.v );                        // 194 210 195 211 198 214 199 215 202 218 203 219 206 222 207 223
    t14   = _mm512_unpacklo_ps( b14.v, b15.v );                        // 224 240 225 241 228 244 229 245 232 248 233 249 236 252 237 253
    t15   = _mm512_unpackhi_ps( b14.v, b15.v );                        // 226 242 227 243 230 246 231 247 234 250 235 251 238 254 239 255

    b00.v = _mm512_shuffle_ps( t00, t02, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //   0  16  32  48   4  20  36  52   8  24  40  56  12  28  44  60
    b01.v = _mm512_shuffle_ps( t00, t02, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //   1  17  33  49   5  21  37  53   9  25  41  57  13  29  45  61
    b02.v = _mm512_shuffle_ps( t01, t03, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //   2  18  34  50   6  22  38  54  10  26  42  58  14  30  46  62
    b03.v = _mm512_shuffle_ps( t01, t03, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //   3  19  35  51   7  23  39  55  11  27  43  59  15  31  47  63
    b04.v = _mm512_shuffle_ps( t04, t06, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //  64  80  96 112  68  84 100 116  72  88 104 120  76  92 108 124
    b05.v = _mm512_shuffle_ps( t04, t06, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //  65  81  97 113  69  85 101 117  73  89 105 121  77  93 109 125
    b06.v = _mm512_shuffle_ps( t05, t07, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //  66  82  98 114  70  86 102 118  74  90 106 122  78  94 110 126
    b07.v = _mm512_shuffle_ps( t05, t07, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //  67  83  99 115  71  87 103 119  75  91 107 123  79  95 111 127
    b08.v = _mm512_shuffle_ps( t08, t10, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 128 144 160 176 132 148 164 180 136 152 168 184 140 156 172 188
    b09.v = _mm512_shuffle_ps( t08, t10, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 129 145 161 177 133 149 165 181 137 153 169 185 141 157 173 189
    b10.v = _mm512_shuffle_ps( t09, t11, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 130 146 162 178 134 150 166 182 138 154 170 186 142 158 174 190
    b11.v = _mm512_shuffle_ps( t09, t11, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 131 147 163 179 135 151 167 183 139 155 171 187 143 159 175 191
    b12.v = _mm512_shuffle_ps( t12, t14, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 192 208 224 240 196 212 228 244 200 216 232 248 204 220 236 252
    b13.v = _mm512_shuffle_ps( t12, t14, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 193 209 225 241 197 213 229 245 201 217 233 249 205 221 237 253
    b14.v = _mm512_shuffle_ps( t13, t15, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 194 210 226 242 198 214 230 246 202 218 234 250 206 222 238 254
    b15.v = _mm512_shuffle_ps( t13, t15, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 195 211 227 243 199 215 231 247 203 219 235 251 207 223 239 255

    t00   = _mm512_shuffle_f32x4( b00.v, b04.v, 0x88 );                //   0  16  32  48   8  24  40  56  64  80  96 112  72  88 104 120
    t01   = _mm512_shuffle_f32x4( b01.v, b05.v, 0x88 );                //   1  17  33  49   9  25  41  57  65  81  97 113  73  89 105 121
    t02   = _mm512_shuffle_f32x4( b02.v, b06.v, 0x88 );                //   2  18  34  50  10  26  42  58  66  82  98 114  74  90 106 122
    t03   = _mm512_shuffle_f32x4( b03.v, b07.v, 0x88 );                //   3  19  35  51  11  27  43  59  67  83  99 115  75  91 107 123
    t04   = _mm512_shuffle_f32x4( b00.v, b04.v, 0xdd );                //   4  20  36  52  12  28  44  60  68  84 100 116  76  92 108 124
    t05   = _mm512_shuffle_f32x4( b01.v, b05.v, 0xdd );                //   5  21  37  53  13  29  45  61  69  85 101 117  77  93 109 125
    t06   = _mm512_shuffle_f32x4( b02.v, b06.v, 0xdd );                //   6  22  38  54  14  30  46  62  70  86 102 118  78  94 110 126
    t07   = _mm512_shuffle_f32x4( b03.v, b07.v, 0xdd );                //   7  23  39  55  15  31  47  63  71  87 103 119  79  95 111 127
    t08   = _mm512_shuffle_f32x4( b08.v, b12.v, 0x88 );                // 128 144 160 176 136 152 168 184 192 208 224 240 200 216 232 248
    t09   = _mm512_shuffle_f32x4( b09.v, b13.v, 0x88 );                // 129 145 161 177 137 153 169 185 193 209 225 241 201 217 233 249
    t10   = _mm512_shuffle_f32x4( b10.v, b14.v, 0x88 );                // 130 146 162 178 138 154 170 186 194 210 226 242 202 218 234 250
    t11   = _mm512_shuffle_f32x4( b11.v, b15.v, 0x88 );                // 131 147 163 179 139 155 171 187 195 211 227 243 203 219 235 251
    t12   = _mm512_shuffle_f32x4( b08.v, b12.v, 0xdd );                // 132 148 164 180 140 156 172 188 196 212 228 244 204 220 236 252
    t13   = _mm512_shuffle_f32x4( b09.v, b13.v, 0xdd );                // 133 149 165 181 141 157 173 189 197 213 229 245 205 221 237 253
    t14   = _mm512_shuffle_f32x4( b10.v, b14.v, 0xdd );                // 134 150 166 182 142 158 174 190 198 214 230 246 206 222 238 254
    t15   = _mm512_shuffle_f32x4( b11.v, b15.v, 0xdd );                // 135 151 167 183 143 159 175 191 199 215 231 247 207 223 239 255

    b00.v = _mm512_shuffle_f32x4( t00, t08, 0x88 );                    //   0  16  32  48  64  80  96 112 128 144 160 176 192 208 224 240
    b01.v = _mm512_shuffle_f32x4( t01, t09, 0x88 );                    //   1  17  33  49  66  81  97 113 129 145 161 177 193 209 225 241
    b02.v = _mm512_shuffle_f32x4( t02, t10, 0x88 );                    //   2  18  34  50  67  82  98 114 130 146 162 178 194 210 226 242
    b03.v = _mm512_shuffle_f32x4( t03, t11, 0x88 );                    //   3  19  35  51  68  83  99 115 131 147 163 179 195 211 227 243
    b04.v = _mm512_shuffle_f32x4( t04, t12, 0x88 );                    //   4  20  36  52  69  84 100 116 132 148 164 180 196 212 228 244
    b05.v = _mm512_shuffle_f32x4( t05, t13, 0x88 );                    //   5  21  37  53  70  85 101 117 133 149 165 181 197 213 229 245
    b06.v = _mm512_shuffle_f32x4( t06, t14, 0x88 );                    //   6  22  38  54  71  86 102 118 134 150 166 182 198 214 230 246
    b07.v = _mm512_shuffle_f32x4( t07, t15, 0x88 );                    //   7  23  39  55  72  87 103 119 135 151 167 183 199 215 231 247
    b08.v = _mm512_shuffle_f32x4( t00, t08, 0xdd );                    //   8  24  40  56  73  88 104 120 136 152 168 184 200 216 232 248
    b09.v = _mm512_shuffle_f32x4( t01, t09, 0xdd );                    //   9  25  41  57  74  89 105 121 137 153 169 185 201 217 233 249
    b10.v = _mm512_shuffle_f32x4( t02, t10, 0xdd );                    //  10  26  42  58  75  90 106 122 138 154 170 186 202 218 234 250
    b11.v = _mm512_shuffle_f32x4( t03, t11, 0xdd );                    //  11  27  43  59  76  91 107 123 139 155 171 187 203 219 235 251
    b12.v = _mm512_shuffle_f32x4( t04, t12, 0xdd );                    //  12  28  44  60  77  92 108 124 140 156 172 188 204 220 236 252
    b13.v = _mm512_shuffle_f32x4( t05, t13, 0xdd );                    //  13  29  45  61  78  93 109 125 141 157 173 189 205 221 237 253
    b14.v = _mm512_shuffle_f32x4( t06, t14, 0xdd );                    //  14  30  46  62  79  94 110 126 142 158 174 190 206 222 238 254
    b15.v = _mm512_shuffle_f32x4( t07, t15, 0xdd );                    //  15  31  47  63  79  95 111 127 143 159 175 191 207 223 239 255
  }

#if 0
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
#endif

  inline void load_16x8_tr_p( const void * ALIGNED(64) a00,
			      const void * ALIGNED(64) a01,
			      const void * ALIGNED(64) a02,
			      const void * ALIGNED(64) a03,
			      const void * ALIGNED(64) a04,
			      const void * ALIGNED(64) a05,
			      const void * ALIGNED(64) a06,
			      const void * ALIGNED(64) a07,
			      v16 &b00, v16 &b01, v16 &b02, v16 &b03,
			      v16 &b04, v16 &b05, v16 &b06, v16 &b07 )
  {
    __m512 t00, t01, t02, t03, t04, t05, t06, t07;

    __m512i idx = _mm512_set_epi32( 15, 11, 14, 10, 13, 9, 12, 8, 7, 3, 6, 2, 5, 1, 4, 0 );

    b00.v = _mm512_load_ps( (const float *)a00 );                      //   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    b01.v = _mm512_load_ps( (const float *)a01 );                      //  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31
    b02.v = _mm512_load_ps( (const float *)a02 );                      //  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47
    b03.v = _mm512_load_ps( (const float *)a03 );                      //  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63
    b04.v = _mm512_load_ps( (const float *)a04 );                      //  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79
    b05.v = _mm512_load_ps( (const float *)a05 );                      //  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95
    b06.v = _mm512_load_ps( (const float *)a06 );                      //  96  97  98  99 100 101 102 103 104 105 106 107 108 109 110 111
    b07.v = _mm512_load_ps( (const float *)a07 );                      // 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127
    // b08.v = _mm512_load_ps( (const float *)a08 );                      // 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143
    // b09.v = _mm512_load_ps( (const float *)a09 );                      // 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159
    // b10.v = _mm512_load_ps( (const float *)a10 );                      // 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175
    // b11.v = _mm512_load_ps( (const float *)a11 );                      // 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191
    // b12.v = _mm512_load_ps( (const float *)a12 );                      // 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207
    // b13.v = _mm512_load_ps( (const float *)a13 );                      // 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223
    // b14.v = _mm512_load_ps( (const float *)a14 );                      // 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239
    // b15.v = _mm512_load_ps( (const float *)a15 );                      // 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255

    t00   = _mm512_unpacklo_ps( b00.v, b01.v );                        //   0  16   1  17   4  20   5  21   8  24   9  25  12  28  13  29 
    t01   = _mm512_unpackhi_ps( b00.v, b01.v );                        //   2  18   3  19   6  22   7  23  10  26  11  27  14  30  15  31
    t02   = _mm512_unpacklo_ps( b02.v, b03.v );                        //  32  48  33  49  36  52  37  53  40  56  41  57  44  60  45  61
    t03   = _mm512_unpackhi_ps( b02.v, b03.v );                        //  34  50  35  51  38  54  39  55  42  58  43  59  46  62  47  63
    t04   = _mm512_unpacklo_ps( b04.v, b05.v );                        //  64  80  65  81  68  84  69  85  72  88  73  89  76  92  77  93
    t05   = _mm512_unpackhi_ps( b04.v, b05.v );                        //  66  82  67  83  70  86  71  87  74  90  75  91  78  94  79  95
    t06   = _mm512_unpacklo_ps( b06.v, b07.v );                        //  96 112  97 113 100 116 101 117 104 120 105 121 108 124 109 125
    t07   = _mm512_unpackhi_ps( b06.v, b07.v );                        //  98 114  99 115 102 118 103 119 106 122 107 123 110 126 111 127
    // t08   = _mm512_unpacklo_ps( b08.v, b09.v );                        // 128 144 129 145 132 148 133 149 136 152 137 153 140 156 141 157
    // t09   = _mm512_unpackhi_ps( b08.v, b09.v );                        // 130 146 131 147 134 150 135 151 138 154 139 155 142 158 143 159
    // t10   = _mm512_unpacklo_ps( b10.v, b11.v );                        // 160 176 161 177 164 180 165 181 168 184 169 185 172 188 173 189
    // t11   = _mm512_unpackhi_ps( b10.v, b11.v );                        // 162 178 163 179 166 182 167 183 170 186 171 187 174 190 175 191
    // t12   = _mm512_unpacklo_ps( b12.v, b13.v );                        // 192 208 193 209 196 212 197 213 200 216 201 217 204 220 205 221
    // t13   = _mm512_unpackhi_ps( b12.v, b13.v );                        // 194 210 195 211 198 214 199 215 202 218 203 219 206 222 207 223
    // t14   = _mm512_unpacklo_ps( b14.v, b15.v );                        // 224 240 225 241 228 244 229 245 232 248 233 249 236 252 237 253
    // t15   = _mm512_unpackhi_ps( b14.v, b15.v );                        // 226 242 227 243 230 246 231 247 234 250 235 251 238 254 239 255

    b00.v = _mm512_shuffle_ps( t00, t02, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //   0  16  32  48   4  20  36  52   8  24  40  56  12  28  44  60
    b01.v = _mm512_shuffle_ps( t00, t02, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //   1  17  33  49   5  21  37  53   9  25  41  57  13  29  45  61
    b02.v = _mm512_shuffle_ps( t01, t03, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //   2  18  34  50   6  22  38  54  10  26  42  58  14  30  46  62
    b03.v = _mm512_shuffle_ps( t01, t03, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //   3  19  35  51   7  23  39  55  11  27  43  59  15  31  47  63
    b04.v = _mm512_shuffle_ps( t04, t06, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //  64  80  96 112  68  84 100 116  72  88 104 120  76  92 108 124
    b05.v = _mm512_shuffle_ps( t04, t06, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //  65  81  97 113  69  85 101 117  73  89 105 121  77  93 109 125
    b06.v = _mm512_shuffle_ps( t05, t07, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //  66  82  98 114  70  86 102 118  74  90 106 122  78  94 110 126
    b07.v = _mm512_shuffle_ps( t05, t07, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //  67  83  99 115  71  87 103 119  75  91 107 123  79  95 111 127
    // b08.v = _mm512_shuffle_ps( t08, t10, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 128 144 160 176 132 148 164 180 136 152 168 184 140 156 172 188
    // b09.v = _mm512_shuffle_ps( t08, t10, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 129 145 161 177 133 149 165 181 137 153 169 185 141 157 173 189
    // b10.v = _mm512_shuffle_ps( t09, t11, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 130 146 162 178 134 150 166 182 138 154 170 186 142 158 174 190
    // b11.v = _mm512_shuffle_ps( t09, t11, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 131 147 163 179 135 151 167 183 139 155 171 187 143 159 175 191
    // b12.v = _mm512_shuffle_ps( t12, t14, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 192 208 224 240 196 212 228 244 200 216 232 248 204 220 236 252
    // b13.v = _mm512_shuffle_ps( t12, t14, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 193 209 225 241 197 213 229 245 201 217 233 249 205 221 237 253
    // b14.v = _mm512_shuffle_ps( t13, t15, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 194 210 226 242 198 214 230 246 202 218 234 250 206 222 238 254
    // b15.v = _mm512_shuffle_ps( t13, t15, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 195 211 227 243 199 215 231 247 203 219 235 251 207 223 239 255

    t00   = _mm512_shuffle_f32x4( b00.v, b04.v, 0x88 );                //   0  16  32  48   8  24  40  56  64  80  96 112  72  88 104 120
    t01   = _mm512_shuffle_f32x4( b01.v, b05.v, 0x88 );                //   1  17  33  49   9  25  41  57  65  81  97 113  73  89 105 121
    t02   = _mm512_shuffle_f32x4( b02.v, b06.v, 0x88 );                //   2  18  34  50  10  26  42  58  66  82  98 114  74  90 106 122
    t03   = _mm512_shuffle_f32x4( b03.v, b07.v, 0x88 );                //   3  19  35  51  11  27  43  59  67  83  99 115  75  91 107 123
    t04   = _mm512_shuffle_f32x4( b00.v, b04.v, 0xdd );                //   4  20  36  52  12  28  44  60  68  84 100 116  76  92 108 124
    t05   = _mm512_shuffle_f32x4( b01.v, b05.v, 0xdd );                //   5  21  37  53  13  29  45  61  69  85 101 117  77  93 109 125
    t06   = _mm512_shuffle_f32x4( b02.v, b06.v, 0xdd );                //   6  22  38  54  14  30  46  62  70  86 102 118  78  94 110 126
    t07   = _mm512_shuffle_f32x4( b03.v, b07.v, 0xdd );                //   7  23  39  55  15  31  47  63  71  87 103 119  79  95 111 127
    // t08   = _mm512_shuffle_f32x4( b08.v, b12.v, 0x88 );                // 128 144 160 176 136 152 168 184 192 208 224 240 200 216 232 248
    // t09   = _mm512_shuffle_f32x4( b09.v, b13.v, 0x88 );                // 129 145 161 177 137 153 169 185 193 209 225 241 201 217 233 249
    // t10   = _mm512_shuffle_f32x4( b10.v, b14.v, 0x88 );                // 130 146 162 178 138 154 170 186 194 210 226 242 202 218 234 250
    // t11   = _mm512_shuffle_f32x4( b11.v, b15.v, 0x88 );                // 131 147 163 179 139 155 171 187 195 211 227 243 203 219 235 251
    // t12   = _mm512_shuffle_f32x4( b08.v, b12.v, 0xdd );                // 132 148 164 180 140 156 172 188 196 212 228 244 204 220 236 252
    // t13   = _mm512_shuffle_f32x4( b09.v, b13.v, 0xdd );                // 133 149 165 181 141 157 173 189 197 213 229 245 205 221 237 253
    // t14   = _mm512_shuffle_f32x4( b10.v, b14.v, 0xdd );                // 134 150 166 182 142 158 174 190 198 214 230 246 206 222 238 254
    // t15   = _mm512_shuffle_f32x4( b11.v, b15.v, 0xdd );                // 135 151 167 183 143 159 175 191 199 215 231 247 207 223 239 255

    b00.v = _mm512_permutexvar_ps( idx, t00 );                         //   0   8  16  24  32  40  48  56  64  72  80  88  96 104 112 120
    b01.v = _mm512_permutexvar_ps( idx, t01 );                         //   1   9  17  25  33  41  49  57  65  73  81  89  97 105 113 121
    b02.v = _mm512_permutexvar_ps( idx, t02 );                         //   2  10  18  26  34  42  50  58  66  74  82  90  98 106 114 122
    b03.v = _mm512_permutexvar_ps( idx, t03 );                         //   3  11  19  27  35  43  51  59  67  75  83  91  99 107 115 123
    b04.v = _mm512_permutexvar_ps( idx, t04 );                         //   4  12  20  28  36  44  52  60  68  76  84  92 100 108 116 124
    b05.v = _mm512_permutexvar_ps( idx, t05 );                         //   5  13  21  29  37  45  53  61  69  77  85  93 101 109 117 125
    b06.v = _mm512_permutexvar_ps( idx, t06 );                         //   6  14  22  30  38  46  54  62  70  78  86  94 102 110 118 126
    b07.v = _mm512_permutexvar_ps( idx, t07 );                         //   7  15  23  31  39  47  55  63  71  79  87  95 103 111 119 127
    // b08.v = _mm512_permutexvar_ps( idx, t08 );                         // 128 136 144 152 160 168 176 184 192 200 208 216 224 232 240 248
    // b09.v = _mm512_permutexvar_ps( idx, t09 );                         // 129 137 145 153 161 169 177 185 193 201 209 217 225 233 241 249
    // b10.v = _mm512_permutexvar_ps( idx, t10 );                         // 130 138 146 154 162 170 178 186 194 202 210 218 226 234 242 250
    // b11.v = _mm512_permutexvar_ps( idx, t11 );                         // 131 139 147 155 163 171 179 187 195 203 211 219 227 235 243 251
    // b12.v = _mm512_permutexvar_ps( idx, t12 );                         // 132 140 148 156 164 172 180 188 196 204 212 220 228 236 244 252
    // b13.v = _mm512_permutexvar_ps( idx, t13 );                         // 133 141 149 157 165 173 181 189 197 205 213 221 229 237 245 253
    // b14.v = _mm512_permutexvar_ps( idx, t14 );                         // 134 142 150 158 166 174 182 190 198 206 214 222 230 238 246 254
    // b15.v = _mm512_permutexvar_ps( idx, t15 );                         // 135 143 151 159 167 175 183 191 199 207 215 223 231 239 247 255
  }

#if 0
  inline void load_16x8_tr_p( const void * ALIGNED(64) a00,
			      const void * ALIGNED(64) a01,
			      const void * ALIGNED(64) a02,
			      const void * ALIGNED(64) a03,
			      const void * ALIGNED(64) a04,
			      const void * ALIGNED(64) a05,
			      const void * ALIGNED(64) a06,
			      const void * ALIGNED(64) a07,
			      v16 &b00, v16 &b01, v16 &b02, v16 &b03,
			      v16 &b04, v16 &b05, v16 &b06, v16 &b07 )
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
  }
#endif

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
    __m512 t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14, t15;

    __m512i idx = _mm512_set_epi32( 15, 11, 14, 10, 13, 9, 12, 8, 7, 3, 6, 2, 5, 1, 4, 0 );

    b00.v = _mm512_load_ps( (const float *)a00 );                      //   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    b01.v = _mm512_load_ps( (const float *)a01 );                      //  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31
    b02.v = _mm512_load_ps( (const float *)a02 );                      //  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47
    b03.v = _mm512_load_ps( (const float *)a03 );                      //  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63
    b04.v = _mm512_load_ps( (const float *)a04 );                      //  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79
    b05.v = _mm512_load_ps( (const float *)a05 );                      //  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95
    b06.v = _mm512_load_ps( (const float *)a06 );                      //  96  97  98  99 100 101 102 103 104 105 106 107 108 109 110 111
    b07.v = _mm512_load_ps( (const float *)a07 );                      // 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127
    b08.v = _mm512_load_ps( (const float *)a08 );                      // 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143
    b09.v = _mm512_load_ps( (const float *)a09 );                      // 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159
    b10.v = _mm512_load_ps( (const float *)a10 );                      // 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175
    b11.v = _mm512_load_ps( (const float *)a11 );                      // 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191
    b12.v = _mm512_load_ps( (const float *)a12 );                      // 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207
    b13.v = _mm512_load_ps( (const float *)a13 );                      // 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223
    b14.v = _mm512_load_ps( (const float *)a14 );                      // 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239
    b15.v = _mm512_load_ps( (const float *)a15 );                      // 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255

    t00   = _mm512_unpacklo_ps( b00.v, b01.v );                        //   0  16   1  17   4  20   5  21   8  24   9  25  12  28  13  29 
    t01   = _mm512_unpackhi_ps( b00.v, b01.v );                        //   2  18   3  19   6  22   7  23  10  26  11  27  14  30  15  31
    t02   = _mm512_unpacklo_ps( b02.v, b03.v );                        //  32  48  33  49  36  52  37  53  40  56  41  57  44  60  45  61
    t03   = _mm512_unpackhi_ps( b02.v, b03.v );                        //  34  50  35  51  38  54  39  55  42  58  43  59  46  62  47  63
    t04   = _mm512_unpacklo_ps( b04.v, b05.v );                        //  64  80  65  81  68  84  69  85  72  88  73  89  76  92  77  93
    t05   = _mm512_unpackhi_ps( b04.v, b05.v );                        //  66  82  67  83  70  86  71  87  74  90  75  91  78  94  79  95
    t06   = _mm512_unpacklo_ps( b06.v, b07.v );                        //  96 112  97 113 100 116 101 117 104 120 105 121 108 124 109 125
    t07   = _mm512_unpackhi_ps( b06.v, b07.v );                        //  98 114  99 115 102 118 103 119 106 122 107 123 110 126 111 127
    t08   = _mm512_unpacklo_ps( b08.v, b09.v );                        // 128 144 129 145 132 148 133 149 136 152 137 153 140 156 141 157
    t09   = _mm512_unpackhi_ps( b08.v, b09.v );                        // 130 146 131 147 134 150 135 151 138 154 139 155 142 158 143 159
    t10   = _mm512_unpacklo_ps( b10.v, b11.v );                        // 160 176 161 177 164 180 165 181 168 184 169 185 172 188 173 189
    t11   = _mm512_unpackhi_ps( b10.v, b11.v );                        // 162 178 163 179 166 182 167 183 170 186 171 187 174 190 175 191
    t12   = _mm512_unpacklo_ps( b12.v, b13.v );                        // 192 208 193 209 196 212 197 213 200 216 201 217 204 220 205 221
    t13   = _mm512_unpackhi_ps( b12.v, b13.v );                        // 194 210 195 211 198 214 199 215 202 218 203 219 206 222 207 223
    t14   = _mm512_unpacklo_ps( b14.v, b15.v );                        // 224 240 225 241 228 244 229 245 232 248 233 249 236 252 237 253
    t15   = _mm512_unpackhi_ps( b14.v, b15.v );                        // 226 242 227 243 230 246 231 247 234 250 235 251 238 254 239 255

    b00.v = _mm512_shuffle_ps( t00, t02, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //   0  16  32  48   4  20  36  52   8  24  40  56  12  28  44  60
    b01.v = _mm512_shuffle_ps( t00, t02, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //   1  17  33  49   5  21  37  53   9  25  41  57  13  29  45  61
    b02.v = _mm512_shuffle_ps( t01, t03, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //   2  18  34  50   6  22  38  54  10  26  42  58  14  30  46  62
    b03.v = _mm512_shuffle_ps( t01, t03, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //   3  19  35  51   7  23  39  55  11  27  43  59  15  31  47  63
    b04.v = _mm512_shuffle_ps( t04, t06, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //  64  80  96 112  68  84 100 116  72  88 104 120  76  92 108 124
    b05.v = _mm512_shuffle_ps( t04, t06, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //  65  81  97 113  69  85 101 117  73  89 105 121  77  93 109 125
    b06.v = _mm512_shuffle_ps( t05, t07, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //  66  82  98 114  70  86 102 118  74  90 106 122  78  94 110 126
    b07.v = _mm512_shuffle_ps( t05, t07, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //  67  83  99 115  71  87 103 119  75  91 107 123  79  95 111 127
    b08.v = _mm512_shuffle_ps( t08, t10, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 128 144 160 176 132 148 164 180 136 152 168 184 140 156 172 188
    b09.v = _mm512_shuffle_ps( t08, t10, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 129 145 161 177 133 149 165 181 137 153 169 185 141 157 173 189
    b10.v = _mm512_shuffle_ps( t09, t11, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 130 146 162 178 134 150 166 182 138 154 170 186 142 158 174 190
    b11.v = _mm512_shuffle_ps( t09, t11, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 131 147 163 179 135 151 167 183 139 155 171 187 143 159 175 191
    b12.v = _mm512_shuffle_ps( t12, t14, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 192 208 224 240 196 212 228 244 200 216 232 248 204 220 236 252
    b13.v = _mm512_shuffle_ps( t12, t14, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 193 209 225 241 197 213 229 245 201 217 233 249 205 221 237 253
    b14.v = _mm512_shuffle_ps( t13, t15, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 194 210 226 242 198 214 230 246 202 218 234 250 206 222 238 254
    b15.v = _mm512_shuffle_ps( t13, t15, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 195 211 227 243 199 215 231 247 203 219 235 251 207 223 239 255

    t00   = _mm512_shuffle_f32x4( b00.v, b04.v, 0x88 );                //   0  16  32  48   8  24  40  56  64  80  96 112  72  88 104 120
    t01   = _mm512_shuffle_f32x4( b01.v, b05.v, 0x88 );                //   1  17  33  49   9  25  41  57  65  81  97 113  73  89 105 121
    t02   = _mm512_shuffle_f32x4( b02.v, b06.v, 0x88 );                //   2  18  34  50  10  26  42  58  66  82  98 114  74  90 106 122
    t03   = _mm512_shuffle_f32x4( b03.v, b07.v, 0x88 );                //   3  19  35  51  11  27  43  59  67  83  99 115  75  91 107 123
    t04   = _mm512_shuffle_f32x4( b00.v, b04.v, 0xdd );                //   4  20  36  52  12  28  44  60  68  84 100 116  76  92 108 124
    t05   = _mm512_shuffle_f32x4( b01.v, b05.v, 0xdd );                //   5  21  37  53  13  29  45  61  69  85 101 117  77  93 109 125
    t06   = _mm512_shuffle_f32x4( b02.v, b06.v, 0xdd );                //   6  22  38  54  14  30  46  62  70  86 102 118  78  94 110 126
    t07   = _mm512_shuffle_f32x4( b03.v, b07.v, 0xdd );                //   7  23  39  55  15  31  47  63  71  87 103 119  79  95 111 127
    t08   = _mm512_shuffle_f32x4( b08.v, b12.v, 0x88 );                // 128 144 160 176 136 152 168 184 192 208 224 240 200 216 232 248
    t09   = _mm512_shuffle_f32x4( b09.v, b13.v, 0x88 );                // 129 145 161 177 137 153 169 185 193 209 225 241 201 217 233 249
    t10   = _mm512_shuffle_f32x4( b10.v, b14.v, 0x88 );                // 130 146 162 178 138 154 170 186 194 210 226 242 202 218 234 250
    t11   = _mm512_shuffle_f32x4( b11.v, b15.v, 0x88 );                // 131 147 163 179 139 155 171 187 195 211 227 243 203 219 235 251
    t12   = _mm512_shuffle_f32x4( b08.v, b12.v, 0xdd );                // 132 148 164 180 140 156 172 188 196 212 228 244 204 220 236 252
    t13   = _mm512_shuffle_f32x4( b09.v, b13.v, 0xdd );                // 133 149 165 181 141 157 173 189 197 213 229 245 205 221 237 253
    t14   = _mm512_shuffle_f32x4( b10.v, b14.v, 0xdd );                // 134 150 166 182 142 158 174 190 198 214 230 246 206 222 238 254
    t15   = _mm512_shuffle_f32x4( b11.v, b15.v, 0xdd );                // 135 151 167 183 143 159 175 191 199 215 231 247 207 223 239 255

    b00.v = _mm512_permutexvar_ps( idx, t00 );                         //   0   8  16  24  32  40  48  56  64  72  80  88  96 104 112 120
    b01.v = _mm512_permutexvar_ps( idx, t01 );                         //   1   9  17  25  33  41  49  57  65  73  81  89  97 105 113 121
    b02.v = _mm512_permutexvar_ps( idx, t02 );                         //   2  10  18  26  34  42  50  58  66  74  82  90  98 106 114 122
    b03.v = _mm512_permutexvar_ps( idx, t03 );                         //   3  11  19  27  35  43  51  59  67  75  83  91  99 107 115 123
    b04.v = _mm512_permutexvar_ps( idx, t04 );                         //   4  12  20  28  36  44  52  60  68  76  84  92 100 108 116 124
    b05.v = _mm512_permutexvar_ps( idx, t05 );                         //   5  13  21  29  37  45  53  61  69  77  85  93 101 109 117 125
    b06.v = _mm512_permutexvar_ps( idx, t06 );                         //   6  14  22  30  38  46  54  62  70  78  86  94 102 110 118 126
    b07.v = _mm512_permutexvar_ps( idx, t07 );                         //   7  15  23  31  39  47  55  63  71  79  87  95 103 111 119 127
    b08.v = _mm512_permutexvar_ps( idx, t08 );                         // 128 136 144 152 160 168 176 184 192 200 208 216 224 232 240 248
    b09.v = _mm512_permutexvar_ps( idx, t09 );                         // 129 137 145 153 161 169 177 185 193 201 209 217 225 233 241 249
    b10.v = _mm512_permutexvar_ps( idx, t10 );                         // 130 138 146 154 162 170 178 186 194 202 210 218 226 234 242 250
    b11.v = _mm512_permutexvar_ps( idx, t11 );                         // 131 139 147 155 163 171 179 187 195 203 211 219 227 235 243 251
    b12.v = _mm512_permutexvar_ps( idx, t12 );                         // 132 140 148 156 164 172 180 188 196 204 212 220 228 236 244 252
    b13.v = _mm512_permutexvar_ps( idx, t13 );                         // 133 141 149 157 165 173 181 189 197 205 213 221 229 237 245 253
    b14.v = _mm512_permutexvar_ps( idx, t14 );                         // 134 142 150 158 166 174 182 190 198 206 214 222 230 238 246 254
    b15.v = _mm512_permutexvar_ps( idx, t15 );                         // 135 143 151 159 167 175 183 191 199 207 215 223 231 239 247 255
  }

#if 0
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
#endif

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
    __m512 t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14, t15;
    __m512 u00, u01, u02, u03, u04, u05, u06, u07, u08, u09, u10, u11, u12, u13, u14, u15;

    // Start                                 a00 =   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    //                                       a01 =  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31
    //                                       a02 =  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47
    //                                       a03 =  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63
    //                                       a04 =  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79
    //                                       a05 =  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95
    //                                       a06 =  96  97  98  99 100 101 102 103 104 105 106 107 108 109 110 111
    //                                       a07 = 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127
    //                                       a08 = 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143
    //                                       a09 = 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159
    //                                       a10 = 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175
    //                                       a11 = 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191
    //                                       a12 = 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207
    //                                       a13 = 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223
    //                                       a14 = 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239
    //                                       a15 = 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255

    t00 = _mm512_unpacklo_ps( b00.v, b01.v ); //   0  16   1  17   4  20   5  21   8  24   9  25  12  28  13  29 
    t01 = _mm512_unpackhi_ps( b00.v, b01.v ); //   2  18   3  19   6  22   7  23  10  26  11  27  14  30  15  31
    t02 = _mm512_unpacklo_ps( b02.v, b03.v ); //  32  48  33  49  36  52  37  53  40  56  41  57  44  60  45  61
    t03 = _mm512_unpackhi_ps( b02.v, b03.v ); //  34  50  35  51  38  54  39  55  42  58  43  59  46  62  47  63
    t04 = _mm512_unpacklo_ps( b04.v, b05.v ); //  64  80  65  81  68  84  69  85  72  88  73  89  76  92  77  93
    t05 = _mm512_unpackhi_ps( b04.v, b05.v ); //  66  82  67  83  70  86  71  87  74  90  75  91  78  94  79  95
    t06 = _mm512_unpacklo_ps( b06.v, b07.v ); //  96 112  97 113 100 116 101 117 104 120 105 121 108 124 109 125
    t07 = _mm512_unpackhi_ps( b06.v, b07.v ); //  98 114  99 115 102 118 103 119 106 122 107 123 110 126 111 127
    t08 = _mm512_unpacklo_ps( b08.v, b09.v ); // 128 144 129 145 132 148 133 149 136 152 137 153 140 156 141 157
    t09 = _mm512_unpackhi_ps( b08.v, b09.v ); // 130 146 131 147 134 150 135 151 138 154 139 155 142 158 143 159
    t10 = _mm512_unpacklo_ps( b10.v, b11.v ); // 160 176 161 177 164 180 165 181 168 184 169 185 172 188 173 189
    t11 = _mm512_unpackhi_ps( b10.v, b11.v ); // 162 178 163 179 166 182 167 183 170 186 171 187 174 190 175 191
    t12 = _mm512_unpacklo_ps( b12.v, b13.v ); // 192 208 193 209 196 212 197 213 200 216 201 217 204 220 205 221
    t13 = _mm512_unpackhi_ps( b12.v, b13.v ); // 194 210 195 211 198 214 199 215 202 218 203 219 206 222 207 223
    t14 = _mm512_unpacklo_ps( b14.v, b15.v ); // 224 240 225 241 228 244 229 245 232 248 233 249 236 252 237 253
    t15 = _mm512_unpackhi_ps( b14.v, b15.v ); // 226 242 227 243 230 246 231 247 234 250 235 251 238 254 239 255

    u00 = _mm512_shuffle_ps( t00, t02, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //   0  16  32  48 
    u01 = _mm512_shuffle_ps( t00, t02, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //   1  17  33  49 
    u02 = _mm512_shuffle_ps( t01, t03, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //   2  18  34  50 
    u03 = _mm512_shuffle_ps( t01, t03, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //   3  19  35  51 
    u04 = _mm512_shuffle_ps( t04, t06, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //  64  80  96 112 
    u05 = _mm512_shuffle_ps( t04, t06, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //  65  81  97 113 
    u06 = _mm512_shuffle_ps( t05, t07, _MM_SHUFFLE( 1, 0, 1, 0 ) );  //  66  82  98 114 
    u07 = _mm512_shuffle_ps( t05, t07, _MM_SHUFFLE( 3, 2, 3, 2 ) );  //  67  83  99 115 
    u08 = _mm512_shuffle_ps( t08, t10, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 128 144 160 176 
    u09 = _mm512_shuffle_ps( t08, t10, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 129 145 161 177 
    u10 = _mm512_shuffle_ps( t09, t11, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 130 146 162 178 
    u11 = _mm512_shuffle_ps( t09, t11, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 131 147 163 179 
    u12 = _mm512_shuffle_ps( t12, t14, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 192 208 228 240 
    u13 = _mm512_shuffle_ps( t12, t14, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 193 209 229 241 
    u14 = _mm512_shuffle_ps( t13, t15, _MM_SHUFFLE( 1, 0, 1, 0 ) );  // 194 210 230 242 
    u15 = _mm512_shuffle_ps( t13, t15, _MM_SHUFFLE( 3, 2, 3, 2 ) );  // 195 211 231 243 

    t00 = _mm512_shuffle_f32x4( u00, u04, 0x88 ); //   0  16  32  48   8  24  40  56  64  80  96  112 ...
    t01 = _mm512_shuffle_f32x4( u01, u05, 0x88 ); //   1  17  33  49 ...
    t02 = _mm512_shuffle_f32x4( u02, u06, 0x88 ); //   2  18  34  50 ...
    t03 = _mm512_shuffle_f32x4( u03, u07, 0x88 ); //   3  19  35  51 ...
    t04 = _mm512_shuffle_f32x4( u00, u04, 0xdd ); //   4  20  36  52 ...
    t05 = _mm512_shuffle_f32x4( u01, u05, 0xdd ); //   5  21  37  53 ...
    t06 = _mm512_shuffle_f32x4( u02, u06, 0xdd ); //   6  22  38  54 ...
    t07 = _mm512_shuffle_f32x4( u03, u07, 0xdd ); //   7  23  39  55 ...
    t08 = _mm512_shuffle_f32x4( u08, u12, 0x88 ); // 128 144 160 176 ...
    t09 = _mm512_shuffle_f32x4( u09, u13, 0x88 ); // 129 145 161 177 ...
    t10 = _mm512_shuffle_f32x4( u10, u14, 0x88 ); // 130 146 162 178 ...
    t11 = _mm512_shuffle_f32x4( u11, u15, 0x88 ); // 131 147 163 179 ...
    t12 = _mm512_shuffle_f32x4( u08, u12, 0xdd ); // 132 148 164 180 ...
    t13 = _mm512_shuffle_f32x4( u09, u13, 0xdd ); // 133 149 165 181 ...
    t14 = _mm512_shuffle_f32x4( u10, u14, 0xdd ); // 134 150 166 182 ...
    t15 = _mm512_shuffle_f32x4( u11, u15, 0xdd ); // 135 151 167 183 ...

    u00 = _mm512_shuffle_f32x4( t00, t08, 0x88 ); //   0  16  32  48  64  80  96 112 ... 240
    u01 = _mm512_shuffle_f32x4( t01, t09, 0x88 ); //   1  17  33  49  66  81  97 113 ... 241
    u02 = _mm512_shuffle_f32x4( t02, t10, 0x88 ); //   2  18  34  50  67  82  98 114 ... 242
    u03 = _mm512_shuffle_f32x4( t03, t11, 0x88 ); //   3  19  35  51  68  83  99 115 ... 243
    u04 = _mm512_shuffle_f32x4( t04, t12, 0x88 ); //   4 ...
    u05 = _mm512_shuffle_f32x4( t05, t13, 0x88 ); //   5 ...
    u06 = _mm512_shuffle_f32x4( t06, t14, 0x88 ); //   6 ...
    u07 = _mm512_shuffle_f32x4( t07, t15, 0x88 ); //   7 ...
    u08 = _mm512_shuffle_f32x4( t00, t08, 0xdd ); //   8 ...
    u09 = _mm512_shuffle_f32x4( t01, t09, 0xdd ); //   9 ...
    u10 = _mm512_shuffle_f32x4( t02, t10, 0xdd ); //  10 ...
    u11 = _mm512_shuffle_f32x4( t03, t11, 0xdd ); //  11 ...
    u12 = _mm512_shuffle_f32x4( t04, t12, 0xdd ); //  12 ...
    u13 = _mm512_shuffle_f32x4( t05, t13, 0xdd ); //  13 ...
    u14 = _mm512_shuffle_f32x4( t06, t14, 0xdd ); //  14 ...
    u15 = _mm512_shuffle_f32x4( t07, t15, 0xdd ); //  15  31  47  63  79  96 111 127 ... 255

    _mm512_store_ps( (float *)a00, u00 );
    _mm512_store_ps( (float *)a01, u01 );
    _mm512_store_ps( (float *)a02, u02 );
    _mm512_store_ps( (float *)a03, u03 );
    _mm512_store_ps( (float *)a04, u04 );
    _mm512_store_ps( (float *)a05, u05 );
    _mm512_store_ps( (float *)a06, u06 );
    _mm512_store_ps( (float *)a07, u07 );
    _mm512_store_ps( (float *)a08, u08 );
    _mm512_store_ps( (float *)a09, u09 );
    _mm512_store_ps( (float *)a10, u10 );
    _mm512_store_ps( (float *)a11, u11 );
    _mm512_store_ps( (float *)a12, u12 );
    _mm512_store_ps( (float *)a13, u13 );
    _mm512_store_ps( (float *)a14, u14 );
    _mm512_store_ps( (float *)a15, u15 );
  }

#if 0
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
#endif

  inline void store_16x8_tr_p( const v16 &b00,
			       const v16 &b01,
			       const v16 &b02,
			       const v16 &b03,
			       const v16 &b04,
			       const v16 &b05,
			       const v16 &b06,
			       const v16 &b07,
			       void * ALIGNED(64) a00,
			       void * ALIGNED(64) a01,
			       void * ALIGNED(64) a02,
			       void * ALIGNED(64) a03,
			       void * ALIGNED(64) a04,
			       void * ALIGNED(64) a05,
			       void * ALIGNED(64) a06,
			       void * ALIGNED(64) a07 )
  {
    __m512 t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14, t15;
    __m512 u00, u01, u02, u03, u04, u05, u06, u07, u08, u09, u10, u11, u12, u13, u14, u15;

    __m512i idx = _mm512_set_epi32( 15, 13, 11, 9, 14, 12, 10, 8, 7, 5, 3, 1, 6, 4, 2, 0 );

    __m512i idx1, idx2;

    // Start                                                     b00 =   0   8  16  24  32  40  48  56  64  72  80  88  96 104 112 120
    //                                                           b01 =   1   9  17  25  33  41  49  57  65  73  81  89  97 105 113 121
    //                                                           b02 =   2  10  18  26  34  42  50  58  66  74  82  90  98 106 114 122
    //                                                           b03 =   3  11  19  27  35  43  51  59  67  75  83  91  99 107 115 123
    //                                                           b04 =   4  12  20  28  36  44  52  60  68  76  84  92 100 108 116 124
    //                                                           b05 =   5  13  21  29  37  45  53  61  69  77  85  93 101 109 117 125
    //                                                           b06 =   6  14  22  30  38  46  54  62  70  78  86  94 102 110 118 126
    //                                                           b07 =   7  15  23  31  39  47  55  63  71  79  87  95 103 111 119 127
       //                                                           b08 = 128 136 144 152 160 168 176 184 192 200 208 216 224 232 240 248
       //                                                           b09 = 129 137 145 153 161 169 177 185 193 201 209 217 225 233 241 249
       //                                                           b10 = 130 138 146 154 162 170 178 186 194 202 210 218 226 234 242 250
       //                                                           b11 = 131 139 147 155 163 171 179 187 195 203 211 219 227 235 243 251
       //                                                           b12 = 132 140 148 156 164 172 180 188 196 204 212 220 228 236 244 252
       //                                                           b13 = 133 141 149 157 165 173 181 189 197 205 213 221 229 237 245 253
       //                                                           b14 = 134 142 150 158 166 174 182 190 198 206 214 222 230 238 246 254
       //                                                           b15 = 135 143 151 159 167 175 183 191 199 207 215 223 231 239 247 255

    t00 = _mm512_permutexvar_ps( idx, b00.v );                      //   0  16  32  48   8  24  40  56  64  80  96 112  72  88 104 120
    t01 = _mm512_permutexvar_ps( idx, b01.v );                      //   1  17  33  49   9  25  41  57  65  81  97 113  73  89 105 121
    t02 = _mm512_permutexvar_ps( idx, b02.v );                      //   2  18  34  50  10  26  42  58  66  82  98 114  74  90 106 122
    t03 = _mm512_permutexvar_ps( idx, b03.v );                      //   3  19  35  51  11  27  43  59  67  83  99 115  75  91 107 123
    t04 = _mm512_permutexvar_ps( idx, b04.v );                      //   4  20  36  52  12  28  44  60  68  84 100 116  76  92 108 124
    t05 = _mm512_permutexvar_ps( idx, b05.v );                      //   5  21  37  53  13  29  45  61  69  85 101 117  77  93 109 125
    t06 = _mm512_permutexvar_ps( idx, b06.v );                      //   6  22  38  54  14  30  46  62  70  86 102 118  78  94 110 126
    t07 = _mm512_permutexvar_ps( idx, b07.v );                      //   7  23  39  55  15  31  47  63  71  87 103 119  79  95 111 127
    // t08 = _mm512_permutexvar_ps( idx, b08.v );                      // 128 144 160 176 136 152 168 184 192 208 228 240 200 216 232 248
    // t09 = _mm512_permutexvar_ps( idx, b09.v );                      // 129 145 161 177 137 153 169 185 193 209 229 241 201 217 233 249
    // t10 = _mm512_permutexvar_ps( idx, b10.v );                      // 130 146 162 178 138 154 170 186 194 210 230 242 202 218 234 250
    // t11 = _mm512_permutexvar_ps( idx, b11.v );                      // 131 147 163 179 139 155 171 187 195 211 231 243 203 219 235 251
    // t12 = _mm512_permutexvar_ps( idx, b12.v );                      // 132 148 164 180 140 156 172 188 196 212 228 244 204 220 236 252
    // t13 = _mm512_permutexvar_ps( idx, b13.v );                      // 133 149 165 181 141 157 173 189 197 213 229 245 205 221 237 253
    // t14 = _mm512_permutexvar_ps( idx, b14.v );                      // 134 150 166 182 142 158 174 190 198 214 230 246 206 222 238 254
    // t15 = _mm512_permutexvar_ps( idx, b15.v );                      // 135 151 167 183 143 159 175 191 199 215 231 247 207 223 239 255

    idx1 = _mm512_set_epi32(  7+16,  6+16,  5+16,  4+16,  7,  6,  5,  4,  3+16,  2+16, 1+16, 0+16,  3,  2, 1, 0 );
    idx2 = _mm512_set_epi32( 15+16, 14+16, 13+16, 12+16, 15, 14, 13, 12, 11+16, 10+16, 9+16, 8+16, 11, 10, 9, 8 );

    u00 = _mm512_permutex2var_ps( t00, idx1, t04 );                 //   0  16  32  48   4  20  36  52   8  24  40  56  12  28  44  60
    u01 = _mm512_permutex2var_ps( t01, idx1, t05 );                 //   1  17  33  49   5  21  37  53   9  25  41  57  13  29  45  61
    u02 = _mm512_permutex2var_ps( t02, idx1, t06 );                 //   2  18  34  50   6  22  38  54  10  26  42  58  14  30  46  62
    u03 = _mm512_permutex2var_ps( t03, idx1, t07 );                 //   3  19  35  51   7  23  39  55  11  27  43  59  15  31  47  63
    u04 = _mm512_permutex2var_ps( t00, idx2, t04 );                 //  64  80  96 112  68  84 100 116  72  88 104 120  76  92 108 124
    u05 = _mm512_permutex2var_ps( t01, idx2, t05 );                 //  65  81  97 113  69  85 101 117  73  89 105 121  77  93 109 125
    u06 = _mm512_permutex2var_ps( t02, idx2, t06 );                 //  66  82  98 114  70  86 102 118  74  90 106 122  78  94 110 126
    u07 = _mm512_permutex2var_ps( t03, idx2, t07 );                 //  67  83  99 115  71  87 103 119  75  91 107 123  79  95 111 127
    // u08 = _mm512_permutex2var_ps( t08, idx1, t12 );                 // 128 144 160 176 132 148 164 180 136 152 168 184 140 156 172 188
    // u09 = _mm512_permutex2var_ps( t09, idx1, t13 );                 // 129 145 161 177 133 149 165 181 137 153 169 185 141 157 173 189
    // u10 = _mm512_permutex2var_ps( t10, idx1, t14 );                 // 130 146 162 178 134 150 166 182 138 154 170 186 142 158 174 190
    // u11 = _mm512_permutex2var_ps( t11, idx1, t15 );                 // 131 147 163 179 135 151 167 183 139 155 171 187 143 159 175 191
    // u12 = _mm512_permutex2var_ps( t08, idx2, t12 );                 // 192 208 224 240 196 212 228 244 200 216 232 248 204 220 236 252
    // u13 = _mm512_permutex2var_ps( t09, idx2, t13 );                 // 193 209 225 241 197 213 229 245 201 217 233 249 205 221 237 253
    // u14 = _mm512_permutex2var_ps( t10, idx2, t14 );                 // 194 210 226 242 198 214 230 246 202 218 234 250 206 222 238 254
    // u15 = _mm512_permutex2var_ps( t11, idx2, t15 );                 // 195 211 227 243 199 215 231 247 203 219 235 251 207 223 239 255

    t00 = _mm512_shuffle_ps( u00, u01, _MM_SHUFFLE( 1, 0, 1, 0 ) ); //   0  16   1  17   4  20   5  21   8  24   9  25  12  28  13  29 
    t01 = _mm512_shuffle_ps( u02, u03, _MM_SHUFFLE( 1, 0, 1, 0 ) ); //   2  18   3  19   6  22   7  23  10  26  11  27  14  30  15  31
    t02 = _mm512_shuffle_ps( u00, u01, _MM_SHUFFLE( 3, 2, 3, 2 ) ); //  32  48  33  49  36  52  37  53  40  56  41  57  44  60  45  61
    t03 = _mm512_shuffle_ps( u02, u03, _MM_SHUFFLE( 3, 2, 3, 2 ) ); //  34  50  35  51  38  54  39  55  42  58  43  59  46  62  47  63
    t04 = _mm512_shuffle_ps( u04, u05, _MM_SHUFFLE( 1, 0, 1, 0 ) ); //  64  80  65  81  68  84  69  85  72  88  73  89  76  92  77  93
    t05 = _mm512_shuffle_ps( u06, u07, _MM_SHUFFLE( 1, 0, 1, 0 ) ); //  66  82  67  83  70  86  71  87  74  90  75  91  78  94  79  95
    t06 = _mm512_shuffle_ps( u04, u05, _MM_SHUFFLE( 3, 2, 3, 2 ) ); //  96 112  97 113 100 116 101 117 104 120 105 121 108 124 109 125
    t07 = _mm512_shuffle_ps( u06, u07, _MM_SHUFFLE( 3, 2, 3, 2 ) ); //  98 114  99 115 102 118 103 119 106 122 107 123 110 126 111 127
    // t08 = _mm512_shuffle_ps( u08, u09, _MM_SHUFFLE( 1, 0, 1, 0 ) ); // 128 144 129 145 132 148 133 149 136 152 137 153 140 156 141 157
    // t09 = _mm512_shuffle_ps( u10, u11, _MM_SHUFFLE( 1, 0, 1, 0 ) ); // 130 146 131 147 134 150 135 151 138 154 139 155 142 158 143 159
    // t10 = _mm512_shuffle_ps( u08, u09, _MM_SHUFFLE( 3, 2, 3, 2 ) ); // 160 176 161 177 164 180 165 181 168 184 169 185 172 188 173 189
    // t11 = _mm512_shuffle_ps( u10, u11, _MM_SHUFFLE( 3, 2, 3, 2 ) ); // 162 178 163 179 166 182 167 183 170 186 171 187 174 190 175 191
    // t12 = _mm512_shuffle_ps( u12, u13, _MM_SHUFFLE( 1, 0, 1, 0 ) ); // 192 208 193 209 196 212 197 213 200 216 201 217 204 220 205 221
    // t13 = _mm512_shuffle_ps( u14, u15, _MM_SHUFFLE( 1, 0, 1, 0 ) ); // 194 210 195 211 198 214 199 215 202 218 203 219 206 222 207 223
    // t14 = _mm512_shuffle_ps( u12, u13, _MM_SHUFFLE( 3, 2, 3, 2 ) ); // 224 240 225 241 228 244 229 245 232 248 233 249 236 252 237 253
    // t15 = _mm512_shuffle_ps( u14, u15, _MM_SHUFFLE( 3, 2, 3, 2 ) ); // 226 242 227 243 230 246 231 247 234 250 235 251 238 254 239 255

    u00 = _mm512_shuffle_ps( t00, t01, _MM_SHUFFLE( 2, 0, 2, 0 ) ); //   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    u01 = _mm512_shuffle_ps( t00, t01, _MM_SHUFFLE( 3, 1, 3, 1 ) ); //  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31
    u02 = _mm512_shuffle_ps( t02, t03, _MM_SHUFFLE( 2, 0, 2, 0 ) ); //  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47
    u03 = _mm512_shuffle_ps( t02, t03, _MM_SHUFFLE( 3, 1, 3, 1 ) ); //  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63
    u04 = _mm512_shuffle_ps( t04, t05, _MM_SHUFFLE( 2, 0, 2, 0 ) ); //  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79
    u05 = _mm512_shuffle_ps( t04, t05, _MM_SHUFFLE( 3, 1, 3, 1 ) ); //  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95
    u06 = _mm512_shuffle_ps( t06, t07, _MM_SHUFFLE( 2, 0, 2, 0 ) ); //  96  97  98  99 100 101 102 103 104 105 106 107 108 109 110 111
    u07 = _mm512_shuffle_ps( t06, t07, _MM_SHUFFLE( 3, 1, 3, 1 ) ); // 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127
    // u08 = _mm512_shuffle_ps( t08, t09, _MM_SHUFFLE( 2, 0, 2, 0 ) ); // 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143
    // u09 = _mm512_shuffle_ps( t08, t09, _MM_SHUFFLE( 3, 1, 3, 1 ) ); // 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159
    // u10 = _mm512_shuffle_ps( t10, t11, _MM_SHUFFLE( 2, 0, 2, 0 ) ); // 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175
    // u11 = _mm512_shuffle_ps( t10, t11, _MM_SHUFFLE( 3, 1, 3, 1 ) ); // 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191
    // u12 = _mm512_shuffle_ps( t12, t13, _MM_SHUFFLE( 2, 0, 2, 0 ) ); // 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207
    // u13 = _mm512_shuffle_ps( t12, t13, _MM_SHUFFLE( 3, 1, 3, 1 ) ); // 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223
    // u14 = _mm512_shuffle_ps( t14, t15, _MM_SHUFFLE( 2, 0, 2, 0 ) ); // 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239
    // u15 = _mm512_shuffle_ps( t14, t15, _MM_SHUFFLE( 3, 1, 3, 1 ) ); // 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255

    _mm512_store_ps( (float *)a00, u00 );
    _mm512_store_ps( (float *)a01, u01 );
    _mm512_store_ps( (float *)a02, u02 );
    _mm512_store_ps( (float *)a03, u03 );
    _mm512_store_ps( (float *)a04, u04 );
    _mm512_store_ps( (float *)a05, u05 );
    _mm512_store_ps( (float *)a06, u06 );
    _mm512_store_ps( (float *)a07, u07 );
    // _mm512_store_ps( (float *)a08, u08 );
    // _mm512_store_ps( (float *)a09, u09 );
    // _mm512_store_ps( (float *)a10, u10 );
    // _mm512_store_ps( (float *)a11, u11 );
    // _mm512_store_ps( (float *)a12, u12 );
    // _mm512_store_ps( (float *)a13, u13 );
    // _mm512_store_ps( (float *)a14, u14 );
    // _mm512_store_ps( (float *)a15, u15 );
  }

#if 0
  inline void store_16x8_tr_p( const v16 &b00,
			       const v16 &b01,
			       const v16 &b02,
			       const v16 &b03,
			       const v16 &b04,
			       const v16 &b05,
			       const v16 &b06,
			       const v16 &b07,
			       void * ALIGNED(64) a00,
			       void * ALIGNED(64) a01,
			       void * ALIGNED(64) a02,
			       void * ALIGNED(64) a03,
			       void * ALIGNED(64) a04,
			       void * ALIGNED(64) a05,
			       void * ALIGNED(64) a06,
			       void * ALIGNED(64) a07 )
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
  }
#endif

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
    __m512 t00, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, t11, t12, t13, t14, t15;
    __m512 u00, u01, u02, u03, u04, u05, u06, u07, u08, u09, u10, u11, u12, u13, u14, u15;

    __m512i idx = _mm512_set_epi32( 15, 13, 11, 9, 14, 12, 10, 8, 7, 5, 3, 1, 6, 4, 2, 0 );

    __m512i idx1, idx2;

    // Start                                                     b00 =   0   8  16  24  32  40  48  56  64  72  80  88  96 104 112 120
    //                                                           b01 =   1   9  17  25  33  41  49  57  65  73  81  89  97 105 113 121
    //                                                           b02 =   2  10  18  26  34  42  50  58  66  74  82  90  98 106 114 122
    //                                                           b03 =   3  11  19  27  35  43  51  59  67  75  83  91  99 107 115 123
    //                                                           b04 =   4  12  20  28  36  44  52  60  68  76  84  92 100 108 116 124
    //                                                           b05 =   5  13  21  29  37  45  53  61  69  77  85  93 101 109 117 125
    //                                                           b06 =   6  14  22  30  38  46  54  62  70  78  86  94 102 110 118 126
    //                                                           b07 =   7  15  23  31  39  47  55  63  71  79  87  95 103 111 119 127
    //                                                           b08 = 128 136 144 152 160 168 176 184 192 200 208 216 224 232 240 248
    //                                                           b09 = 129 137 145 153 161 169 177 185 193 201 209 217 225 233 241 249
    //                                                           b10 = 130 138 146 154 162 170 178 186 194 202 210 218 226 234 242 250
    //                                                           b11 = 131 139 147 155 163 171 179 187 195 203 211 219 227 235 243 251
    //                                                           b12 = 132 140 148 156 164 172 180 188 196 204 212 220 228 236 244 252
    //                                                           b13 = 133 141 149 157 165 173 181 189 197 205 213 221 229 237 245 253
    //                                                           b14 = 134 142 150 158 166 174 182 190 198 206 214 222 230 238 246 254
    //                                                           b15 = 135 143 151 159 167 175 183 191 199 207 215 223 231 239 247 255

    t00 = _mm512_permutexvar_ps( idx, b00.v );                      //   0  16  32  48   8  24  40  56  64  80  96 112  72  88 104 120
    t01 = _mm512_permutexvar_ps( idx, b01.v );                      //   1  17  33  49   9  25  41  57  65  81  97 113  73  89 105 121
    t02 = _mm512_permutexvar_ps( idx, b02.v );                      //   2  18  34  50  10  26  42  58  66  82  98 114  74  90 106 122
    t03 = _mm512_permutexvar_ps( idx, b03.v );                      //   3  19  35  51  11  27  43  59  67  83  99 115  75  91 107 123
    t04 = _mm512_permutexvar_ps( idx, b04.v );                      //   4  20  36  52  12  28  44  60  68  84 100 116  76  92 108 124
    t05 = _mm512_permutexvar_ps( idx, b05.v );                      //   5  21  37  53  13  29  45  61  69  85 101 117  77  93 109 125
    t06 = _mm512_permutexvar_ps( idx, b06.v );                      //   6  22  38  54  14  30  46  62  70  86 102 118  78  94 110 126
    t07 = _mm512_permutexvar_ps( idx, b07.v );                      //   7  23  39  55  15  31  47  63  71  87 103 119  79  95 111 127
    t08 = _mm512_permutexvar_ps( idx, b08.v );                      // 128 144 160 176 136 152 168 184 192 208 228 240 200 216 232 248
    t09 = _mm512_permutexvar_ps( idx, b09.v );                      // 129 145 161 177 137 153 169 185 193 209 229 241 201 217 233 249
    t10 = _mm512_permutexvar_ps( idx, b10.v );                      // 130 146 162 178 138 154 170 186 194 210 230 242 202 218 234 250
    t11 = _mm512_permutexvar_ps( idx, b11.v );                      // 131 147 163 179 139 155 171 187 195 211 231 243 203 219 235 251
    t12 = _mm512_permutexvar_ps( idx, b12.v );                      // 132 148 164 180 140 156 172 188 196 212 228 244 204 220 236 252
    t13 = _mm512_permutexvar_ps( idx, b13.v );                      // 133 149 165 181 141 157 173 189 197 213 229 245 205 221 237 253
    t14 = _mm512_permutexvar_ps( idx, b14.v );                      // 134 150 166 182 142 158 174 190 198 214 230 246 206 222 238 254
    t15 = _mm512_permutexvar_ps( idx, b15.v );                      // 135 151 167 183 143 159 175 191 199 215 231 247 207 223 239 255

    idx1 = _mm512_set_epi32(  7+16,  6+16,  5+16,  4+16,  7,  6,  5,  4,  3+16,  2+16, 1+16, 0+16,  3,  2, 1, 0 );
    idx2 = _mm512_set_epi32( 15+16, 14+16, 13+16, 12+16, 15, 14, 13, 12, 11+16, 10+16, 9+16, 8+16, 11, 10, 9, 8 );

    u00 = _mm512_permutex2var_ps( t00, idx1, t04 );                 //   0  16  32  48   4  20  36  52   8  24  40  56  12  28  44  60
    u01 = _mm512_permutex2var_ps( t01, idx1, t05 );                 //   1  17  33  49   5  21  37  53   9  25  41  57  13  29  45  61
    u02 = _mm512_permutex2var_ps( t02, idx1, t06 );                 //   2  18  34  50   6  22  38  54  10  26  42  58  14  30  46  62
    u03 = _mm512_permutex2var_ps( t03, idx1, t07 );                 //   3  19  35  51   7  23  39  55  11  27  43  59  15  31  47  63
    u04 = _mm512_permutex2var_ps( t00, idx2, t04 );                 //  64  80  96 112  68  84 100 116  72  88 104 120  76  92 108 124
    u05 = _mm512_permutex2var_ps( t01, idx2, t05 );                 //  65  81  97 113  69  85 101 117  73  89 105 121  77  93 109 125
    u06 = _mm512_permutex2var_ps( t02, idx2, t06 );                 //  66  82  98 114  70  86 102 118  74  90 106 122  78  94 110 126
    u07 = _mm512_permutex2var_ps( t03, idx2, t07 );                 //  67  83  99 115  71  87 103 119  75  91 107 123  79  95 111 127
    u08 = _mm512_permutex2var_ps( t08, idx1, t12 );                 // 128 144 160 176 132 148 164 180 136 152 168 184 140 156 172 188
    u09 = _mm512_permutex2var_ps( t09, idx1, t13 );                 // 129 145 161 177 133 149 165 181 137 153 169 185 141 157 173 189
    u10 = _mm512_permutex2var_ps( t10, idx1, t14 );                 // 130 146 162 178 134 150 166 182 138 154 170 186 142 158 174 190
    u11 = _mm512_permutex2var_ps( t11, idx1, t15 );                 // 131 147 163 179 135 151 167 183 139 155 171 187 143 159 175 191
    u12 = _mm512_permutex2var_ps( t08, idx2, t12 );                 // 192 208 224 240 196 212 228 244 200 216 232 248 204 220 236 252
    u13 = _mm512_permutex2var_ps( t09, idx2, t13 );                 // 193 209 225 241 197 213 229 245 201 217 233 249 205 221 237 253
    u14 = _mm512_permutex2var_ps( t10, idx2, t14 );                 // 194 210 226 242 198 214 230 246 202 218 234 250 206 222 238 254
    u15 = _mm512_permutex2var_ps( t11, idx2, t15 );                 // 195 211 227 243 199 215 231 247 203 219 235 251 207 223 239 255

    t00 = _mm512_shuffle_ps( u00, u01, _MM_SHUFFLE( 1, 0, 1, 0 ) ); //   0  16   1  17   4  20   5  21   8  24   9  25  12  28  13  29 
    t01 = _mm512_shuffle_ps( u02, u03, _MM_SHUFFLE( 1, 0, 1, 0 ) ); //   2  18   3  19   6  22   7  23  10  26  11  27  14  30  15  31
    t02 = _mm512_shuffle_ps( u00, u01, _MM_SHUFFLE( 3, 2, 3, 2 ) ); //  32  48  33  49  36  52  37  53  40  56  41  57  44  60  45  61
    t03 = _mm512_shuffle_ps( u02, u03, _MM_SHUFFLE( 3, 2, 3, 2 ) ); //  34  50  35  51  38  54  39  55  42  58  43  59  46  62  47  63
    t04 = _mm512_shuffle_ps( u04, u05, _MM_SHUFFLE( 1, 0, 1, 0 ) ); //  64  80  65  81  68  84  69  85  72  88  73  89  76  92  77  93
    t05 = _mm512_shuffle_ps( u06, u07, _MM_SHUFFLE( 1, 0, 1, 0 ) ); //  66  82  67  83  70  86  71  87  74  90  75  91  78  94  79  95
    t06 = _mm512_shuffle_ps( u04, u05, _MM_SHUFFLE( 3, 2, 3, 2 ) ); //  96 112  97 113 100 116 101 117 104 120 105 121 108 124 109 125
    t07 = _mm512_shuffle_ps( u06, u07, _MM_SHUFFLE( 3, 2, 3, 2 ) ); //  98 114  99 115 102 118 103 119 106 122 107 123 110 126 111 127
    t08 = _mm512_shuffle_ps( u08, u09, _MM_SHUFFLE( 1, 0, 1, 0 ) ); // 128 144 129 145 132 148 133 149 136 152 137 153 140 156 141 157
    t09 = _mm512_shuffle_ps( u10, u11, _MM_SHUFFLE( 1, 0, 1, 0 ) ); // 130 146 131 147 134 150 135 151 138 154 139 155 142 158 143 159
    t10 = _mm512_shuffle_ps( u08, u09, _MM_SHUFFLE( 3, 2, 3, 2 ) ); // 160 176 161 177 164 180 165 181 168 184 169 185 172 188 173 189
    t11 = _mm512_shuffle_ps( u10, u11, _MM_SHUFFLE( 3, 2, 3, 2 ) ); // 162 178 163 179 166 182 167 183 170 186 171 187 174 190 175 191
    t12 = _mm512_shuffle_ps( u12, u13, _MM_SHUFFLE( 1, 0, 1, 0 ) ); // 192 208 193 209 196 212 197 213 200 216 201 217 204 220 205 221
    t13 = _mm512_shuffle_ps( u14, u15, _MM_SHUFFLE( 1, 0, 1, 0 ) ); // 194 210 195 211 198 214 199 215 202 218 203 219 206 222 207 223
    t14 = _mm512_shuffle_ps( u12, u13, _MM_SHUFFLE( 3, 2, 3, 2 ) ); // 224 240 225 241 228 244 229 245 232 248 233 249 236 252 237 253
    t15 = _mm512_shuffle_ps( u14, u15, _MM_SHUFFLE( 3, 2, 3, 2 ) ); // 226 242 227 243 230 246 231 247 234 250 235 251 238 254 239 255

    u00 = _mm512_shuffle_ps( t00, t01, _MM_SHUFFLE( 2, 0, 2, 0 ) ); //   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    u01 = _mm512_shuffle_ps( t00, t01, _MM_SHUFFLE( 3, 1, 3, 1 ) ); //  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31
    u02 = _mm512_shuffle_ps( t02, t03, _MM_SHUFFLE( 2, 0, 2, 0 ) ); //  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47
    u03 = _mm512_shuffle_ps( t02, t03, _MM_SHUFFLE( 3, 1, 3, 1 ) ); //  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63
    u04 = _mm512_shuffle_ps( t04, t05, _MM_SHUFFLE( 2, 0, 2, 0 ) ); //  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79
    u05 = _mm512_shuffle_ps( t04, t05, _MM_SHUFFLE( 3, 1, 3, 1 ) ); //  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95
    u06 = _mm512_shuffle_ps( t06, t07, _MM_SHUFFLE( 2, 0, 2, 0 ) ); //  96  97  98  99 100 101 102 103 104 105 106 107 108 109 110 111
    u07 = _mm512_shuffle_ps( t06, t07, _MM_SHUFFLE( 3, 1, 3, 1 ) ); // 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127
    u08 = _mm512_shuffle_ps( t08, t09, _MM_SHUFFLE( 2, 0, 2, 0 ) ); // 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143
    u09 = _mm512_shuffle_ps( t08, t09, _MM_SHUFFLE( 3, 1, 3, 1 ) ); // 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159
    u10 = _mm512_shuffle_ps( t10, t11, _MM_SHUFFLE( 2, 0, 2, 0 ) ); // 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175
    u11 = _mm512_shuffle_ps( t10, t11, _MM_SHUFFLE( 3, 1, 3, 1 ) ); // 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191
    u12 = _mm512_shuffle_ps( t12, t13, _MM_SHUFFLE( 2, 0, 2, 0 ) ); // 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207
    u13 = _mm512_shuffle_ps( t12, t13, _MM_SHUFFLE( 3, 1, 3, 1 ) ); // 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223
    u14 = _mm512_shuffle_ps( t14, t15, _MM_SHUFFLE( 2, 0, 2, 0 ) ); // 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239
    u15 = _mm512_shuffle_ps( t14, t15, _MM_SHUFFLE( 3, 1, 3, 1 ) ); // 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255

    _mm512_store_ps( (float *)a00, u00 );
    _mm512_store_ps( (float *)a01, u01 );
    _mm512_store_ps( (float *)a02, u02 );
    _mm512_store_ps( (float *)a03, u03 );
    _mm512_store_ps( (float *)a04, u04 );
    _mm512_store_ps( (float *)a05, u05 );
    _mm512_store_ps( (float *)a06, u06 );
    _mm512_store_ps( (float *)a07, u07 );
    _mm512_store_ps( (float *)a08, u08 );
    _mm512_store_ps( (float *)a09, u09 );
    _mm512_store_ps( (float *)a10, u10 );
    _mm512_store_ps( (float *)a11, u11 );
    _mm512_store_ps( (float *)a12, u12 );
    _mm512_store_ps( (float *)a13, u13 );
    _mm512_store_ps( (float *)a14, u14 );
    _mm512_store_ps( (float *)a15, u15 );
  }

#if 0
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
#endif

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
      v = a.v;
    }

    v16int( const v16 &a )                       // Init from mixed
    {
      v = a.v;
    }

    v16int( int a )                              // Init from scalar
    {
      union
      {
	int i;
	float f;
      } u;
      u.i = a;
      v = _mm512_set1_ps( u.f );
    }

    v16int( int i00, int i01, int i02, int i03,
	    int i04, int i05, int i06, int i07,
	    int i08, int i09, int i10, int i11,
	    int i12, int i13, int i14, int i15 ) // Init from scalars
    {
      union
      {
	int i;
	float f;
      } u00, u01, u02, u03, u04, u05, u06, u07,
	u08, u09, u10, u11, u12, u13, u14, u15;

      u00.i = i00; u01.i = i01; u02.i = i02; u03.i = i03;
      u04.i = i04; u05.i = i05; u06.i = i06; u07.i = i07;
      u08.i = i08; u09.i = i09; u10.i = i10; u11.i = i11;
      u12.i = i12; u13.i = i13; u14.i = i14; u15.i = i15;

      v = _mm512_setr_ps( u00.f, u01.f, u02.f, u03.f,
			  u04.f, u05.f, u06.f, u07.f,
			  u08.f, u09.f, u10.f, u11.f,
			  u12.f, u13.f, u14.f, u15.f );
    }

    ~v16int() {}                                 // Destructor

    // v16int assignment operators

#   define ASSIGN(op)			          \
    inline v16int &operator op( const v16int &b ) \
    {						  \
      for( int j = 0; j < 16; j++ )               \
        i[j] op b.i[j];                           \
      return *this;                               \
    }

    inline v16int &operator =( const v16int &b )
    {
      v = b.v;
      return *this;
    }

    ASSIGN(+=)
    ASSIGN(-=)
    ASSIGN(*=)
    ASSIGN(/=)
    ASSIGN(%=)

    inline v16int &operator ^=( const v16int &b )
    {
      v = _mm512_xor_ps( v, b.v );
      return *this;
    }

    inline v16int &operator &=( const v16int &b )
    {
      v = _mm512_and_ps( v, b.v );
      return *this;
    }

    inline v16int &operator |=( const v16int &b )
    {
      v = _mm512_or_ps( v, b.v );
      return *this;
    }

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
    for( int j = 0; j < 16; j++ )               \
      b.i[j] = ( op a.i[j] );                   \
    return b;                                   \
  }

  inline v16int operator +( const v16int & a )
  {
    v16int b;

    b.v = a.v;

    return b;
  }

  PREFIX_UNARY(-)

  inline v16int operator !( const v16int & a )
  {
    v16int b;

    for( int j = 0; j < 16; j++ )
      b.i[j] = - ( !a.i[j] );

    return b;
  }

  inline v16int operator ~( const v16int & a )
  {
    v16int b;
    union
    {
      int i;
      float f;
    } u;
    u.i = -1;
    b.v = _mm512_xor_ps( a.v, _mm512_set1_ps( u.f ) );
    return b;
  }

# undef PREFIX_UNARY

  // v16int prefix increment / decrement

# define PREFIX_INCDEC(op)                      \
  inline v16int operator op( v16int & a )       \
  {						\
    v16int b;                                   \
    for( int j = 0; j < 16; j++ )               \
      b.i[j] = ( op a.i[j] );                   \
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
    for( int j = 0; j < 16; j++ )              \
      b.i[j] = ( a.i[j] op );                  \
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
    for( int j = 0; j < 16; j++ )                               \
      c.i[j] = a.i[j] op b.i[j];                                \
    return c;                                                   \
  }

  BINARY(+)
  BINARY(-)
  BINARY(*)
  BINARY(/)
  BINARY(%)

  inline v16int operator ^( const v16int &a, const v16int &b )
  {
    v16int c;

    c.v = _mm512_xor_ps( a.v, b.v );

    return c;
  }

  inline v16int operator &( const v16int &a, const v16int &b )
  {
    v16int c;

    c.v = _mm512_and_ps( a.v, b.v );

    return c;
  }

#if 0
  inline v16int operator |( const v16int &a, const v16int &b )
  {
    v16int c;

    c.v = _mm512_or_ps( a.v, b.v );

    return c;
  }
#endif

  BINARY(|)
  BINARY(<<)
  BINARY(>>)

# undef BINARY

  // v16int logical operators

# define LOGICAL(op)                                            \
  inline v16int operator op( const v16int &a, const v16int &b ) \
  {                                                             \
    v16int c;                                                   \
    for( int j = 0; j < 16; j++ )                               \
      c.i[j] = - ( a.i[j] op b.i[j] );                          \
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

    for( int j = 0; j < 16; j++ )
      b.i[j] = ( a.i[j] >= 0 ) ? a.i[j] : -a.i[j];

    return b;
  }

  inline v16 czero( const v16int &c, const v16 &a )
  {
    v16 b;

    for( int j = 0; j < 16; j++ )
      b.i[j] = a.i[j] & ~c.i[j];

    return b;
  }

#if 0
  inline v16 czero( const v16int &c, const v16 &a )
  {
    v16 b;

    b.v = _mm512_andnot_ps( c.v, a.v );

    return b;
  }
#endif

  inline v16 notczero( const v16int &c, const v16 &a )
  {
    v16 b;

    b.v = _mm512_and_ps( c.v, a.v );

    return b;
  }

  inline v16 merge( const v16int &c, const v16 &t, const v16 &f )
  {
    v16 m;

    for( int j = 0; j < 16; j++ )
      m.i[j] = ( f.i[j] & ~c.i[j] ) | ( t.i[j] & c.i[j] );

    return m;
  }

#if 0
  inline v16 merge( const v16int &c, const v16 &t, const v16 &f )
  {
    __m512 c_v = c.v;

    v16 tf;

    tf.v = _mm512_or_ps( _mm512_andnot_ps( c_v, f.v ),
			 _mm512_and_ps( c_v, t.v ) );

    return tf;
  }
#endif

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
      v = a.v;
    }

    v16float( const v16 &a )                               // Init from mixed
    {
      v = a.v;
    }

    v16float( float a )                                    // Init from scalar
    {
      v = _mm512_set1_ps( a );
    }

    v16float( float f00, float f01, float f02, float f03,
	      float f04, float f05, float f06, float f07,
	      float f08, float f09, float f10, float f11,
	      float f12, float f13, float f14, float f15 ) // Init from scalars
    {
      v = _mm512_setr_ps( f00, f01, f02, f03, f04, f05, f06, f07,
			  f08, f09, f10, f11, f12, f13, f14, f15 );
    }

    ~v16float() {}                                         // Destructor

    // v16float assignment operators

#   define ASSIGN(op,intrin)				\
    inline v16float &operator op( const v16float &b )   \
    {							\
      v = intrin( v, b.v );                             \
      return *this;                                     \
    }

    inline v16float &operator =( const v16float &b )
    {
      v = b.v;
      return *this;
    }

    ASSIGN( +=, _mm512_add_ps )
    ASSIGN( -=, _mm512_sub_ps )
    ASSIGN( *=, _mm512_mul_ps )
    ASSIGN( /=, _mm512_div_ps )

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

    b.v = a.v;

    return b;
  }

  inline v16float operator -( const v16float &a )
  {
    v16float b;

    b.v = _mm512_sub_ps( _mm512_setzero_ps(), a.v );

    return b;
  }

  inline v16int operator !( const v16float &a )
  {
    v16int b;

    for( int j = 0; j < 16; j++ )
      b.i[j] = a.i[j] ? 0 : -1;

    return b;
  }

#if 0
  inline v16int operator !( const v16float &a )
  {
    v16int b;

    b.v = _mm512_cmp_ps( _mm512_setzero_ps(), a.v, _CMP_EQ_OS );

    return b;
  }
#endif

  // v16float prefix increment / decrement operators

  inline v16float operator ++( v16float &a )
  {
    v16float b;
    __m512 t = _mm512_add_ps( a.v, _mm512_set1_ps( 1.0f ) );

    a.v = t;
    b.v = t;

    return b;
  }

  inline v16float operator --( v16float &a )
  {
    v16float b;
    __m512 t = _mm512_sub_ps( a.v, _mm512_set1_ps( 1.0f ) );

    a.v = t;
    b.v = t;

    return b;
  }

  // v16float postfix increment / decrement operators

  inline v16float operator ++( v16float &a, int )
  {
    v16float b;
    __m512 a_v = a.v;

    a.v = _mm512_add_ps( a_v, _mm512_set1_ps( 1.0f ) );
    b.v = a_v;

    return b;
  }

  inline v16float operator --( v16float &a, int )
  {
    v16float b;
    __m512 a_v = a.v;

    a.v = _mm512_sub_ps( a_v, _mm512_set1_ps( 1.0f ) );
    b.v = a_v;

    return b;
  }

  // v16float binary operators

# define BINARY(op,intrin)                                            \
  inline v16float operator op( const v16float &a, const v16float &b ) \
  {								      \
    v16float c;                                                       \
    c.v = intrin( a.v, b.v );                                         \
    return c;                                                         \
  }

  BINARY( +, _mm512_add_ps )
  BINARY( -, _mm512_sub_ps )
  BINARY( *, _mm512_mul_ps )
  BINARY( /, _mm512_div_ps )

# undef BINARY

  // v16float logical operators

# define LOGICAL(op)                                                \
  inline v16int operator op( const v16float &a, const v16float &b ) \
  {								    \
    v16int c;                                                       \
    for( int j = 0; j < 16; j++ )                                   \
      c.i[j] = -( a.f[j] op b.f[j] );                               \
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

#if 0
# define LOGICAL(op,intrin,flag)                                    \
  inline v16int operator op( const v16float &a, const v16float &b ) \
  {							            \
    v16int c;                                                       \
    c.v = intrin( a.v, b.v, flag );				    \
    return c;                                                       \
  }

  LOGICAL( <,  _mm512_cmp_ps, _CMP_LT_OS  )
  LOGICAL( >,  _mm512_cmp_ps, _CMP_GT_OS  )
  LOGICAL( ==, _mm512_cmp_ps, _CMP_EQ_OS  )
  LOGICAL( !=, _mm512_cmp_ps, _CMP_NEQ_OS )
  LOGICAL( <=, _mm512_cmp_ps, _CMP_LE_OS  )
  LOGICAL( >=, _mm512_cmp_ps, _CMP_GE_OS  )

  inline v16int operator &&( const v16float &a, const v16float &b )
  {
    v16int c;
    __m512 vzero = _mm512_setzero_ps();
    c.v = _mm512_and_ps( _mm512_cmp_ps( a.v, vzero, _CMP_NEQ_OS ),
			 _mm512_cmp_ps( b.v, vzero, _CMP_NEQ_OS ) );
    return c;
  }

  inline v16int operator ||( const v16float &a, const v16float &b )
  {
    v16int c;
    __m512 vzero = _mm512_setzero_ps();
    c.v = _mm512_or_ps( _mm512_cmp_ps( a.v, vzero, _CMP_NEQ_OS ),
			_mm512_cmp_ps( b.v, vzero, _CMP_NEQ_OS ) );
    return c;
  }

# undef LOGICAL
#endif

  // v16float math library functions

# define CMATH_FR1(fn)                          \
  inline v16float fn( const v16float &a )       \
  {						\
    v16float b;                                 \
    for( int j = 0; j < 16; j++ )               \
      b.f[j] = ::fn( a.f[j] );                  \
    return b;                                   \
  }

# define CMATH_FR2(fn)                                          \
  inline v16float fn( const v16float &a, const v16float &b )    \
  {								\
    v16float c;                                                 \
    for( int j = 0; j < 16; j++ )                               \
      c.f[j] = ::fn( a.f[j], b.f[j] );                          \
    return c;                                                   \
  }

  CMATH_FR1(acos)     CMATH_FR1(asin)  CMATH_FR1(atan) CMATH_FR2(atan2)
  CMATH_FR1(ceil)     CMATH_FR1(cos)   CMATH_FR1(cosh) CMATH_FR1(exp)
  /*CMATH_FR1(fabs)*/ CMATH_FR1(floor) CMATH_FR2(fmod) CMATH_FR1(log)
  CMATH_FR1(log10)    CMATH_FR2(pow)   CMATH_FR1(sin)  CMATH_FR1(sinh)
  /*CMATH_FR1(sqrt)*/ CMATH_FR1(tan)   CMATH_FR1(tanh)

  inline v16float fabs( const v16float &a )
  {
    v16float b;

    b.v = _mm512_andnot_ps( _mm512_set1_ps( -0.0f ), a.v );

    return b;
  }

  inline v16float sqrt( const v16float &a )
  {
    v16float b;

    b.v = _mm512_sqrt_ps( a.v );

    return b;
  }

  inline v16float copysign( const v16float &a, const v16float &b )
  {
    v16float c;
    __m512 t = _mm512_set1_ps( -0.0f );

    c.v = _mm512_or_ps( _mm512_and_ps( t, b.v ), _mm512_andnot_ps( t, a.v ) );

    return c;
  }

# undef CMATH_FR1
# undef CMATH_FR2

  // v16float miscellaneous functions

  inline v16float rsqrt_approx( const v16float &a )
  {
    v16float b;

    b.v = _mm512_rsqrt14_ps(a.v);

    // b.v = _mm512_rsqrt28_ps(a.v);

    return b;
  }

  inline v16float rsqrt( const v16float &a )
  {
    v16float b;
    __m512 a_v = a.v, b_v;

    // b_v = _mm512_rsqrt28_ps(a_v);

    b_v = _mm512_rsqrt14_ps(a_v);

    b.v = _mm512_add_ps( b_v, _mm512_mul_ps( _mm512_set1_ps( 0.5f ),
					     _mm512_sub_ps( b_v,
							    _mm512_mul_ps( a_v,
									   _mm512_mul_ps( b_v,
											  _mm512_mul_ps( b_v, b_v ) ) ) ) ) );

    // Note: It is quicker to just call div_ps and sqrt_ps if more refinement
    // is desired.
    // b.v = _mm512_div_ps( _mm512_set1_ps( 1.0f ), _mm512_sqrt_ps( a.v ) );

    return b;
  }

  inline v16float rcp_approx( const v16float &a )
  {
    v16float b;

    // b.v = _mm512_rcp28_ps( a.v );

    b.v = _mm512_rcp14_ps( a.v );

    return b;
  }

  inline v16float rcp( const v16float &a )
  {
    v16float b;
    __m512 a_v = a.v, b_v;

    // b_v = _mm512_rcp28_ps( a_v );

    b_v = _mm512_rcp14_ps( a_v );

    b.v = _mm512_sub_ps( _mm512_add_ps( b_v, b_v ),
			 _mm512_mul_ps( a_v, _mm512_mul_ps( b_v, b_v ) ) );

    // b.v = _mm512_div_ps( _mm512_set1_ps( 1.0f ), a.v );

    return b;
  }

  inline v16float fma( const v16float &a, const v16float &b, const v16float &c )
  {
    v16float d;

    d.v = _mm512_fmadd_ps( a.v, b.v, c.v );

    return d;
  }

  inline v16float fms( const v16float &a, const v16float &b, const v16float &c )
  {
    v16float d;

    d.v = _mm512_fmsub_ps( a.v, b.v, c.v );

    return d;
  }

  inline v16float fnms( const v16float &a, const v16float &b, const v16float &c )
  {
    v16float d;

    d.v = _mm512_fnmadd_ps( a.v, b.v, c.v );

    return d;
  }

  inline v16float clear_bits( const v16int &m, const v16float &a )
  {
    v16float b;

    b.v = _mm512_andnot_ps( m.v, a.v );

    return b;
  }

  inline v16float set_bits( const v16int &m, const v16float &a )
  {
    v16float b;

    b.v = _mm512_or_ps( m.v, a.v );

    return b;
  }

  inline v16float toggle_bits( const v16int &m, const v16float &a )
  {
    v16float b;

    b.v = _mm512_xor_ps( m.v, a.v );

    return b;
  }

  inline void increment_16x1( float * ALIGNED(64) p, const v16float &a )
  {
    _mm512_store_ps( p, _mm512_add_ps( _mm512_load_ps( p ), a.v ) );
  }

  inline void decrement_16x1( float * ALIGNED(64) p, const v16float &a )
  {
    _mm512_store_ps( p, _mm512_sub_ps( _mm512_load_ps( p ), a.v ) );
  }

  inline void scale_16x1( float * ALIGNED(64) p, const v16float &a )
  {
    _mm512_store_ps( p, _mm512_mul_ps( _mm512_load_ps( p ), a.v ) );
  }

} // namespace v16

#endif // _v16_avx512_h_
