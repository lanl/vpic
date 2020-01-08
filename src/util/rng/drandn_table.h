#ifndef _drand_table_h_
#define _drand_table_h_

#ifndef IN_rng
#error "Do not include drand_table.h; use rng.h"
#endif

#define DRANDN_N 256
#define DRANDN_R ( 3.6554204190269413594915892673498092335649e+00 )

extern const double drandn_zig_x[];
extern const double drandn_zig_y[];

#endif /* _drandn_table_h_ */
