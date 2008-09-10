#ifndef _dumpvars_h_
#define _dumpvars_h_
 
/* number of variables in header; total is 48, 192 bytes total */
#define NWHEAD_INT 24
#define NWHEAD_REAL 24
#define NWHEAD_TOT NWHEAD_INT+NWHEAD_REAL
 
#define NVARFLDMX 16
#define NVARHYDMX 14
#define NVARPARMX 8
#define NVARHISMX 250
#define NVARGRDMX 4
 
#define NFILMX 1
#define NFILHISMX 10
#define LENNAM 8
#define LENLAB 16
 
/* (code and version idents, not currently used) */
#define CODE 118112105099
#define VERS 20080404
 
/* integer flag to dump specific variables */
int iffldvar[NVARFLDMX], ifhydvar[NVARHYDMX], ifparvar[NVARPARMX];
 
int ndmphis;
int ifhisvar[NVARHISMX];   /* on-off switch for his not used */
int nfilhis=0;
 
int ifgrdvar[NVARGRDMX];
 
/* total number of variables of each type dumped */
int nvarfld, nvarhyd, nvarpar, nvarhis, nvargrd;
 
/* names of variables (for diagnostics) */
/* lengths are + 1 for null termination */
char namfld[NVARFLDMX][LENNAM+1], namhyd[NVARHYDMX][LENNAM+1],
     nampar[NVARPARMX][LENNAM+1], namhis[NVARHISMX][LENNAM+1];
char namgrd[NVARGRDMX+3][LENNAM+1];
char cenfld[NVARFLDMX][LENNAM+1];
char cenhyd[NVARHYDMX][LENNAM+1];
 
char names_field[NVARFLDMX][LENNAM+1] = {
 "ex      ","ey      ","ez      ","diveerr ",
 "cbx     ","cby     ","cbz     ","divberr ",
 "tcax    ","tcay    ","tcaz    ","rhob    ",
 "jfx     ","jfy     ","jfz     ","rhof    " };
char center_field[NVARFLDMX][LENNAM+1] = {
 "edgex   ","edgey   ","edgez   ","node    ",
 "facex   ","facey   ","facez   ","cell    ",
 "edgex   ","edgey   ","edgez   ","node    ",
 "edgex   ","edgey   ","edgez   ","node    " };
char names_hydro[NVARHYDMX][LENNAM+1] = {
 "jx      ","jy      ","jz      ","rho     ",
 "px      ","py      ","pz      ","ke      ",
 "txx     ","tyy     ","tzz     ","tyz     ",
 "tzx     ","txy     " };
char center_hydro[NVARHYDMX][LENNAM+1] = {
 "node    ","node    ","node    ","node    ",
 "node    ","node    ","node    ","node    ",
 "node    ","node    ","node    ","node    ",
 "node    ","node    " };
char names_particle[NVARPARMX][LENNAM+1] = {
 "pardx   ","pardy   ","pardz   ","parcell ",
 "parux   ","paruy   ","paruz   ","parq    " };
char names_grid[NVARGRDMX][LENNAM+1] = {
 "mata    ","matb    ","matc    ","matd    "};
 
char names_history[NVARHISMX][LENNAM+1];
 
/* labels of variables (for diagnostics) */
// char labfld[NVARFLDMX][LENLAB+1], labhyd[NVARHYDMX][LENLAB+1],
//      labpar[NVARPARMX][LENLAB+1];
// char labhis[NVARHISMX][LENLAB+1];   /* his not yet implemented */
 
/* units conversion factors to be applied for diagnostic variables */
float unitfld[NVARFLDMX], unithyd[NVARHYDMX], unitpar[NVARPARMX];
// float unithis[NVARHISMX];   /* his not yet implemented */
 
/* integer dump types for standard point-structured dump files: */
/* grd=10, fld=11, hyd=12, par=13, his=15                       */
/* for block-structured:                                        */
/* grd=20, fld=21, hyd=22 (not used for par or his dumps)       */
int itypgrd=10, itypfld=11, ityphyd=12, ityppar=13, ityphis=15;
 
char fdmpnam[256], fheadnam[256], fcatnam[256];
 
/* option to be used for strided output of fields and hydro dumps */
/* stropt=1 first index(default), 2 center index, 3 average value */
int stropt=1;
 
// long headers for the strided dumps. at present, there is a two-word
// initial header, followed by 24 integers and 24 reals. this can be
// modified when needed.
//
// two initial header words:
//   ityphead = integer, dump type (fld, hyd, par, his, grd, etc.)
//
// 24 integer variables in each header:
//   nrnk = integer, rank number of this processor.
//   nprc = integer, total number of processors.
//   itopo = integer, number of processors in i (1) dimension.
//   jtopo = integer, number of processors in j (2) dimension.
//   ktopo = integer, number of processors in k (3) dimension.
//   nvar = integer, number of variables in dump.
//   imx = integer, i (1) dimension of processor arrays.
//   jmx = integer, j (2) dimension of processor arrays.
//   kmx = integer, k (3) dimension of processor arrays.
//   istr = integer, stride in i direction.
//   jstr = integer, stride in j direction.
//   kstr = integer, stride in k direction.
//   imxstr = integer, i (1) dimension of strided processor arrays.
//   jmxstr = integer, j (2) dimension of strided processor arrays.
//   kmxstr = integer, k (3) dimension of strided processor arrays.
//   icelglb = integer, i direction number of global cells.
//   jcelblg = integer, j direction number of global cells.
//   kcelglb = integer, k direction number of global cells.
//   ncyc = integer, cycle (timestep) number.
//   ndmp = integer, number of this dump.
//   npar = integer, number of particles this processor.
//   pstr = integer, stride for particles.
//   nparstr = integer, number of strided processor particles.
//   ipspec = integer, particle species number (0 otherwise).
//
// 24 real variables in each header:
//   tim = real, physical time.
//   dt = real, physical timestep.
//   cvac = real, speed of light in vacuum.
//   eps0 = real, epsilon0 dielectric constant.
//   damp = real, damping factor.
//   qom  = real, q/m ratio for particle dump (otherwise 0).
//   xmn = real, x minimum coordinate.
//   xmx = real, x maximum coordinate.
//   ymn = real, y minimum coordinate.
//   ymx = real, y maximum coordinate.
//   zmn = real, z minimum coordinate.
//   zmx = real, z maximum coordinate.
//   xmnglb = real, x global minimum coordinate.
//   xmxglb = real, x global maximum coordinate.
//   ymnglb = real, y global minimum coordinate.
//   ymxglb = real, y global maximum coordinate.
//   zmnglb = real, z global minimum coordinate.
//   zmxglb = real, z global maximum coordinate.
//   dx = real, delta-x cell size, without stride.
//   dy = real, delta-y cell size, without stride.
//   dz = real, delta-z cell size, without stride.
//   dxstr = real, delta-x including stride (dx*istr).
//   dystr = real, delta-y including stride (dy*jstr).
//   dzstr = real, delta-z including stride (dz*kstr).
 
int   ityp;
 
int   ityphead;
 
int   nrnk,    nprc,    itopo,   jtopo,   ktopo,   nvar,
      imx,     jmx,     kmx,     istr,    jstr,    kstr,
      imxstr,  jmxstr,  kmxstr,  icelglb, jcelglb, kcelglb,
      ncyc,    ndmp,    npar,    pstr,    nparstr, ipspec;
 
float tim,     dt,      cvac,    eps0,    damp,    qom,
      xmn,     xmx,     ymn,     ymx,     zmn,     zmx,
      xmnglb,  xmxglb,  ymnglb,  ymxglb,  zmnglb,  zmxglb,
      dx,      dy,      dz,      dxstr,   dystr,   dzstr;
 
/* data types idattyp currently used to specify header blocks are:  */
/* 101=integer, 103=real, 105=names, 106=labels, 107=centerings     */
/*                        (names, labels, and center are strings)   */
/* 199=user defined (arbitrary length, must appear last in header)  */
int itypdat;
 
 
#endif   /* endif _dumpvars_h_ */
