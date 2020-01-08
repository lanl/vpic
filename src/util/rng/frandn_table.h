#ifndef _frandn_table_h_
#define _frandn_table_h_

#ifndef IN_rng
#error "Do not include frand_table.h; use rng.h"
#endif

#define FRANDN_N 64
#define FRANDN_R ( 3.2159292455085228233657712593185351579450e+00f )

extern const float frandn_zig_x[];
extern const float frandn_zig_y[];

#endif /* _frandn_table_h_ */
