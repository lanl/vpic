/* --------------------------------------------------------------  */
/* (C)Copyright 2006,2008,                                         */
/* International Business Machines Corporation                     */
/* All Rights Reserved.                                            */
/*                                                                 */
/* Redistribution and use in source and binary forms, with or      */
/* without modification, are permitted provided that the           */
/* following conditions are met:                                   */
/*                                                                 */
/* - Redistributions of source code must retain the above copyright*/
/*   notice, this list of conditions and the following disclaimer. */
/*                                                                 */
/* - Redistributions in binary form must reproduce the above       */
/*   copyright notice, this list of conditions and the following   */
/*   disclaimer in the documentation and/or other materials        */
/*   provided with the distribution.                               */
/*                                                                 */
/* - Neither the name of IBM Corporation nor the names of its      */
/*   contributors may be used to endorse or promote products       */
/*   derived from this software without specific prior written     */
/*   permission.                                                   */
/*                                                                 */
/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND          */
/* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,     */
/* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF        */
/* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE        */
/* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR            */
/* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    */
/* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT    */
/* NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)        */
/* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN       */
/* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR    */
/* OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,  */
/* EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              */
/* --------------------------------------------------------------  */
/* PROLOG END TAG zYx                                              */
/*
 * cache-4way.h
 *
 * service definitions for direct mapped and 4-way
 * set associative cache
 *
 * this file gets included by the top-level cache
 * header file if CACHE_LOG2NWAY is 0 or 2
 *
 */
#ifndef _CACHE_4WAY_H_
#define _CACHE_4WAY_H_
#include <assert.h>

/* callable interfaces */
#define __cache_rd			CACHE_FUNC (cache_rd)
#define __cache_wr			CACHE_FUNC (cache_wr)
#define __cache_rw			CACHE_FUNC (cache_rw)
#define __cache_flush			CACHE_FUNC (cache_flush)
#define __cache_wait			CACHE_FUNC (cache_wait)

/* internal interfaces */
#define __cache_ea_base                 CACHE_VAR (cache_ea_base)
#define __cache_replace			CACHE_FUNC (cache_replace)
#define __cache_replace_idx		CACHE_FUNC (cache_replace_idx)
#define __cache_rd_miss			CACHE_FUNC (cache_rd_miss)
#define __cache_set_lookup		CACHE_FUNC (cache_set_lookup)
#define __cache_line_lookup		CACHE_FUNC (cache_line_lookup)
#define __cache_writeback		CACHE_FUNC (cache_writeback)

#ifdef CACHE_READ_X4
#define __cache_rd_x4			CACHE_FUNC (cache_rd_x4)
#define __cache_addr_x4			CACHE_FUNC (cache_addr_x4)
#define __cacheline_num_x4		CACHE_FUNC (cacheline_num_x4)
#define __cacheline_byte_offset_x4 	CACHE_FUNC (cacheline_byte_offset_x4)
#define __cache_set_num_x4		CACHE_FUNC (cache_set_num_x4)
#endif

/* the cache ea base */
static unsigned long long __cache_ea_base;

/*
 * internal utility macros
 */

/* first line number of a set */
#define __cache_dir_idx(set) ((set) << CACHE_NWAY_SHIFT)

/* Get tagid for a given set. */
#define __cache_tagid(set) CACHE_SET_TAGID ((set))

/* Get tagmask for a given set. */
#define __cache_tagmask(set) (1 << CACHE_SET_TAGID ((set)))

/* Get cached_mem line number. */
#define __cacheline_num(set, idx) \
    (((set) << CACHE_NWAY_SHIFT) + (idx))

/* Hash EA and compute its set number. */
#define __cache_set_num(ea)		\
    ((((ea) >> CACHELINE_SHIFT) ^	\
      ((ea) >> (CACHELINE_SHIFT + 1)))	\
     & CACHE_NSETS_MASK)

/* byte offset within cacheline */
#define __cacheline_byte_offset(ea) ((ea) & ~CACHE_ALIGN_MASK)

/* byte offset from start of cache */
#define __cache_offset(lnum, byte) (((lnum) << CACHELINE_SHIFT) + (byte))

/* write a value to the cache and update the directory */
#define __cache_value(ea, lnum, byte, val)	\
{						\
    *CACHE_ADDR((lnum), (byte)) = (val);	\
    CACHELINE_SETTAGMOD(lnum, ea);              \
}

#define CACHE_ADDR(lnum, byte) \
    ((CACHED_TYPE *)&__cache_mem[__cache_offset((lnum), (byte))]) 

/* utility functions */

/* Pack four scalars into a quadword. */
static inline vec_uint4
__load_vec_uint4 (unsigned int ui1,
	unsigned int ui2, unsigned int ui3, unsigned int ui4)
{
    vec_uchar16 shuffle = (vec_uchar16) {
	0x00, 0x01, 0x02, 0x03, 0x10, 0x11, 0x12, 0x13,
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80
    };
    vec_uint4 iv1 = spu_promote ((unsigned int) ui1, 0);
    vec_uint4 iv2 = spu_promote ((unsigned int) ui2, 0);
    vec_uint4 iv3 = spu_promote ((unsigned int) ui3, 0);
    vec_uint4 iv4 = spu_promote ((unsigned int) ui4, 0);
    return spu_or (spu_shuffle (iv1, iv2, shuffle),
	    spu_shuffle (iv3, iv4, spu_rlqwbyte (shuffle, 8)));
}

/* Pack four quadwords into one. */
static inline vec_uint4
__pack_vec_uint4 (vec_uint4 ui1, vec_uint4 ui2, vec_uint4 ui3, vec_uint4 ui4)
{
    vec_uchar16 shuffle = (vec_uchar16) {
	0x00, 0x01, 0x02, 0x03, 0x10, 0x11, 0x12, 0x13,
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80
    };
    return spu_or (spu_shuffle (ui1, ui2, shuffle),
	    spu_shuffle (ui3, ui4, spu_rlqwbyte (shuffle, 8)));
}
#endif /* end of defined once stuff */

#if (CACHE_NWAY == 1)
static inline int
__cache_replace_idx (int set __attribute__ ((unused)))
{
    return 0;
}

/* Lookup EA in set directory. */
static inline int
__cache_line_lookup (unsigned int set, unsigned int ea)
{
    return (CACHELINE_ISTAG(set, ea) ? (int)set : -1);
}

#define __cache_4way_init(name,base) /* Safe outside these headers */   \
  CACHE_SYM(name,cache_ea_base) = (base)

#else

#define __cache_replace_cntr CACHE_VAR (cache_replace_cntr)
static unsigned * RESTRICT __attribute__ ((aligned (16))) __cache_replace_cntr;

/* Instead of clearing this during init, we could use CACHE_NWAY_MASK
   in code which accesses this. */
/* FIXME: DUE TO EVIL HACKERY INVOLVED HERE, ALL CACHES IN A FILE MUST
   BE EITHER 4-WAY OR 1-WAY. */
#define __cache_4way_init(name,base) /* Safe outside these headers */ \
  SPU_MALLOC( CACHE_SYM(name,cache_replace_cntr),                     \
              CACHE_SYM(name,cache_nsets), 16 );                      \
  CLEAR( CACHE_SYM(name,cache_replace_cntr),                          \
         CACHE_SYM(name,cache_nsets) );                               \
  CACHE_SYM(name,cache_ea_base) = (base)

/* Lookup EA in set directory. */
static inline vec_uint4
__cache_set_lookup (int set, vec_uint4 ea4)
{
    vec_uint4 dir = spu_and(*(vec_uint4 *)CACHEDIR_ADDR(__cache_dir_idx(set)),
	    CACHELINE_TAG_MASK|CACHELINE_VALID);
    vec_uint4 ea4_valid = spu_or (ea4, CACHELINE_VALID);

    return spu_gather (spu_cmpeq (dir, ea4_valid));
}

/* Lookup a line given a set */
static inline int
__cache_line_lookup (int set, vec_uint4 ea4)
{
    vec_uint4 exists = __cache_set_lookup (set, ea4);
    if (likely (spu_extract (exists, 0) != 0))
    {
	int idx = CACHE_NWAY_MASK - spu_extract (spu_sub ((unsigned int) 31,
		    spu_cntlz (exists)), 0);
	return __cacheline_num (set, idx);
    }
    return -1;
}

/* Select set entry to be replaced. */
static inline int
__cache_replace_idx (int set)
{
    unsigned int curr = __cache_replace_cntr[set];

    /*
     * if the line we selected happens to be locked,
     * find an unlocked line.  this is not quite RR,
     * but is close enough considering the potential
     * expense of maintaing an RR cursor of unlocked lines.
     */
    if (unlikely (CACHELINE_ISLOCKED(__cacheline_num(set, curr))))
    {
	vec_uint4 unlocked = spu_gather (spu_cmpeq (spu_and (
			*(vec_uint4 *)CACHEDIR_ADDR(__cache_dir_idx(set)),
			CACHELINE_LOCKED), spu_splats ((unsigned)0)));
	curr = CACHE_NWAY_MASK - spu_extract (spu_sub ((unsigned int) 31,
		    spu_cntlz (unlocked)), 0);
	/* make sure we found an unlocked line */
	assert (curr <= CACHE_NWAY_MASK);
    }
    __cache_replace_cntr[set]  = ((curr + 1) & CACHE_NWAY_MASK);
    return curr;
}
#endif

#if (CACHE_TYPE == CACHE_TYPE_RW)
/* write back a modified line and clear mod bit */
static inline void
__cache_writeback(unsigned int ea_aligned, unsigned int lnum,
	unsigned int tagid)
{
    mfc_put (CACHE_ADDR(lnum, 0), __cache_ea_base + ea_aligned,
             CACHELINE_SIZE, tagid, 0, 0);
    CACHELINE_CLEARMOD(lnum);
}
#endif


/* chooses a line to replace and does write-back as needed */
static inline int
__cache_replace (int set)
{
    int idx = __cache_replace_idx (set);
    int lnum = __cacheline_num (set, idx);
#if (CACHE_TYPE == CACHE_TYPE_RW)
    if (CACHELINE_ISMOD(lnum))
    {
	__cache_writeback (CACHELINE_GETTAG(lnum), lnum, __cache_tagid (set));
	CACHE_STATS_WRITEBACK(set);
    }
#endif
    CACHE_STATS_REPLACE(set);
    return lnum;
}

/*
 * Fetch a missing line and update directory.
 * Caller must ensure that the line is not already resident
 */
static inline int
__cache_rd_miss (unsigned int ea_aligned, int set)
{
    int lnum;
#ifdef DEBUG
#if (CACHE_NWAY == 1)
    lnum = __cache_line_lookup (set, ea_aligned);
#else
    lnum = __cache_line_lookup (set, spu_splats (ea_aligned));
#endif
    assert (lnum < 0);
#endif
    lnum = __cache_replace (set);

    mfc_getf (CACHE_ADDR (lnum, 0), __cache_ea_base + ea_aligned,
              CACHELINE_SIZE, __cache_tagid (set), 0, 0);

    CACHELINE_SETTAG (lnum, ea_aligned);

    return lnum;
}

/* Look up EA and return cached value. */
static inline CACHED_TYPE
__cache_rd (unsigned int ea) {
    unsigned int ea_aligned = ea & CACHE_ALIGN_MASK;
    int set  = __cache_set_num (ea);
#if (CACHE_NWAY == 1)
    int lnum = __cache_line_lookup (set, ea_aligned);
#else
    int lnum = __cache_line_lookup (set, spu_splats (ea_aligned));
#endif
    int byte = __cacheline_byte_offset (ea);
    if (unlikely (lnum < 0)) {
	CACHE_STATS_RD_MISS (set);
	lnum = __cache_rd_miss (ea_aligned, set);
	spu_writech(MFC_WrTagMask, __cache_tagmask (set));
	spu_mfcstat(MFC_TAG_UPDATE_ALL);
    }
#ifdef CACHE_STATS
    else {
	CACHE_STATS_RD_HIT (set);
    }
#endif
    return *CACHE_ADDR (lnum, byte);
}

static inline CACHED_TYPE *
__cache_rw (unsigned int ea, unsigned hit_wait, unsigned miss_wait)
{
    unsigned int ea_aligned = ea & CACHE_ALIGN_MASK;
    int set  = __cache_set_num (ea);
#if (CACHE_NWAY == 1)
    int lnum = __cache_line_lookup (set, ea_aligned);
#else
    int lnum = __cache_line_lookup (set, spu_splats (ea_aligned));
#endif
    int byte = __cacheline_byte_offset (ea);

    CACHED_TYPE *ret = CACHE_ADDR (lnum, byte);

    if (unlikely (lnum < 0))
    {
	CACHE_STATS_RD_MISS(set);
	lnum = __cache_rd_miss (ea_aligned, set);
	if( miss_wait ) {
	    spu_writech (MFC_WrTagMask, __cache_tagmask (set));
	    spu_mfcstat (MFC_TAG_UPDATE_ALL);
	}
#if (CACHE_TYPE == CACHE_TYPE_RW)
	CACHELINE_SETMOD (lnum);
#endif
	return CACHE_ADDR (lnum, byte);
    }
    else {
        if( hit_wait ) {
	    spu_writech (MFC_WrTagMask, __cache_tagmask (set));
	    spu_mfcstat (MFC_TAG_UPDATE_ALL);
        }
#ifdef CACHE_STATS
	CACHE_STATS_RD_HIT (set);
#endif
    }
#if (CACHE_TYPE == CACHE_TYPE_RW)
    CACHELINE_SETMOD (lnum);
#endif
    return ret;
}

/* wait for a previous touch to finish */
static inline void
__cache_wait (unsigned lsa __attribute__ ((unused))) {
    spu_writech (MFC_WrTagMask,
            __cache_tagmask (CACHE_ADDR2LNUM(lsa) >> CACHE_NWAY_SHIFT));
    spu_mfcstat (MFC_TAG_UPDATE_ALL);
}

#if (CACHE_TYPE == CACHE_TYPE_RW)
/* write value to the cache */
static inline void
__cache_wr(unsigned int ea, CACHED_TYPE val)
{
    unsigned int ea_aligned = ea & CACHE_ALIGN_MASK;
    int set  = __cache_set_num (ea);
#if (CACHE_NWAY == 1)
    int lnum = __cache_line_lookup (set, ea_aligned);
#else
    int lnum = __cache_line_lookup (set, spu_splats (ea_aligned));
#endif
    int byte = __cacheline_byte_offset (ea);

    if (unlikely (lnum < 0))
    {
	CACHE_STATS_WR_MISS (set);
	lnum = __cache_rd_miss (ea_aligned, set);
	spu_writech (MFC_WrTagMask, __cache_tagmask (set));
	spu_mfcstat (MFC_TAG_UPDATE_ALL);
	__cache_value (ea_aligned, lnum, byte, val);
	return;
    }
#ifdef CACHE_STATS
    else {
	CACHE_STATS_WR_HIT (set);
    }
#endif
    __cache_value (ea_aligned, lnum, byte, val);
    return;
}

/* write back all dirty lines */
static inline void
__cache_flush()
{
    unsigned lnum = 0;
    while (likely (lnum < CACHE_NWAY*CACHE_NSETS)) {
	if ((CACHELINE_ISMOD(lnum)))
	{
	    __cache_writeback (CACHELINE_GETTAG(lnum),
		    lnum, __cache_tagid(lnum >> CACHE_NWAY_SHIFT));
	    CACHE_STATS_FLUSH(lnum >> CACHE_NWAY_SHIFT);
	    spu_writech (MFC_WrTagMask, __cache_tagmask (lnum >> CACHE_NWAY_SHIFT));
	    spu_mfcstat (MFC_TAG_UPDATE_ALL);
	}
	lnum++;
    }
}
#endif /* r/w */

/*
 * If the cache is 4-way or more, and the cached type
 * is an integral type 32-bits or smaller, CACHE_READ_X4
 * can be defined to enable the rd_x4 interface.
 */
#if defined(CACHE_READ_X4) && (CACHE_NWAY >= 4)
static inline vec_uint4
__cache_addr_x4(vec_uint4 lnum, vec_uint4 byte)
{
    return spu_add (spu_add (spu_sl (lnum, CACHELINE_SHIFT), byte),
	    (unsigned)&__cache_mem[0]);
}

/* Get cached_mem line number, x4. */
static inline vec_uint4
__cacheline_num_x4 (vec_uint4 set, vec_int4 idx)
{
    return spu_add (spu_sl (set, CACHE_NWAY_SHIFT), (vec_uint4) idx);
}

/* Get byte offset within cacheline, x4. */
static inline vec_uint4
__cacheline_byte_offset_x4 (vec_uint4 ea)
{
    return spu_and (ea, (vec_uint4) spu_splats (~CACHE_ALIGN_MASK));
}

/* Hash EA and compute its set number, x4. */
static inline vec_uint4
__cache_set_num_x4 (vec_uint4 ea_x4)
{
    vec_uint4 tmp0 = spu_rlmask (ea_x4, -CACHELINE_SHIFT);
    vec_uint4 tmp1 = spu_rlmask (ea_x4, -(CACHELINE_SHIFT + 1));
    return spu_and (spu_xor (tmp0, tmp1),
	    (vec_uint4) spu_splats (CACHE_NSETS_MASK));
}

/* Look up EA and cache/return data, x4. */
static inline vec_uint4
__cache_rd_x4 (vec_uint4 ea_x4)
{
    vec_uint4 ea_aligned_x4 =
	spu_and ((ea_x4), (vec_uint4) spu_splats (CACHE_ALIGN_MASK));

    vec_uint4 s_x4 = __cache_set_num_x4 (ea_x4);

    int s0 = spu_extract (s_x4, 0);
    int s1 = spu_extract (s_x4, 1);
    int s2 = spu_extract (s_x4, 2);
    int s3 = spu_extract (s_x4, 3);

    /*
     * construct a vector for each ea to lookup.
     * this allows the set lookup service to check
     * 4 lines at a time with one spu op.
     */
    vec_uint4 ea_aligned0 = spu_splats (spu_extract (ea_aligned_x4, 0));
    vec_uint4 ea_aligned1 = spu_splats (spu_extract (ea_aligned_x4, 1));
    vec_uint4 ea_aligned2 = spu_splats (spu_extract (ea_aligned_x4, 2));
    vec_uint4 ea_aligned3 = spu_splats (spu_extract (ea_aligned_x4, 3));

    vec_uint4 found0 = __cache_set_lookup (s0, ea_aligned0);
    vec_uint4 found1 = __cache_set_lookup (s1, ea_aligned1);
    vec_uint4 found2 = __cache_set_lookup (s2, ea_aligned2);
    vec_uint4 found3 = __cache_set_lookup (s3, ea_aligned3);

    vec_uint4 found_x4 = __pack_vec_uint4 (found0, found1, found2, found3);

    /*
     * create vector of indices from the results of the set lookup.
     * -lines that were found will have an index of [0-(NWAY-1)]
     * -lines that were not found will have an index of NWAY.
     */
    vec_int4 i_x4 = (vec_int4)spu_sub (spu_add (spu_cntlz (found_x4),
		(unsigned)CACHE_NWAY_MASK), spu_splats((unsigned)31));

    vec_uint4 missing = spu_rlmask ((vec_uint4) i_x4, -CACHE_NWAY_SHIFT);
    unsigned int ms = spu_extract (spu_gather (missing), 0);

    vec_uint4 ibyte = __cacheline_byte_offset_x4 (ea_x4);
    vec_uint4 iline = __cacheline_num_x4 (s_x4, i_x4);
    vec_uint4 iaddr = __cache_addr_x4 (iline, ibyte);

    unsigned int d0 = *(CACHED_TYPE *)spu_extract(iaddr, 0);
    unsigned int d1 = *(CACHED_TYPE *)spu_extract(iaddr, 1);
    unsigned int d2 = *(CACHED_TYPE *)spu_extract(iaddr, 2);
    unsigned int d3 = *(CACHED_TYPE *)spu_extract(iaddr, 3);

    vec_uint4 ret = __load_vec_uint4 (d0, d1, d2, d3);

    /*
     * if we missed on any line, read them in one at a time,
     * locking them as we go so we don't lose any along the
     * way
     */
    if (unlikely (ms))
    {
	unsigned int ea0 = spu_extract (ea_aligned_x4, 0);
	unsigned int ea1 = spu_extract (ea_aligned_x4, 1);
	unsigned int ea2 = spu_extract (ea_aligned_x4, 2);
	unsigned int ea3 = spu_extract (ea_aligned_x4, 3);

	int lnum0 = spu_extract(iline, 0);
	int lnum1 = spu_extract(iline, 1);
	int lnum2 = spu_extract(iline, 2);
	int lnum3 = spu_extract(iline, 3);

	unsigned int tagmask = 0;

	 /*
	 * recheck each line that was missing to
	 * see if it's still missing so we don't
	 * do a "read miss" unless we have to
	 */
	if (unlikely (ms & 0x8)) {
	    CACHE_STATS_RD_MISS(s0);
	    tagmask |= __cache_tagmask (s0);
	    lnum0 = __cache_rd_miss (ea0, s0);
	}
#ifdef CACHE_STATS
	else {
	    CACHE_STATS_RD_HIT(s0);
	}
#endif
	CACHELINE_LOCK (lnum0);
	if (unlikely (ms & 0x4) &&
	    unlikely ((lnum1 = __cache_line_lookup (s1, ea_aligned1)) < 0))
	{
	    CACHE_STATS_RD_MISS (s1);
	    tagmask |= __cache_tagmask (s1);
	    lnum1 = __cache_rd_miss (ea1, s1);
	}
#ifdef CACHE_STATS
	else {
	    CACHE_STATS_RD_HIT (s1);
	}
#endif
	CACHELINE_LOCK (lnum1);
	if (unlikely (ms & 0x2) &&
	    unlikely ((lnum2 = __cache_line_lookup (s2, ea_aligned2)) < 0))
	{
	    CACHE_STATS_RD_MISS(s2);
	    tagmask |= __cache_tagmask (s2);
	    lnum2 = __cache_rd_miss (ea2, s2);
	}
#ifdef CACHE_STATS
	else {
	    CACHE_STATS_RD_HIT (s2);
	}
#endif
	CACHELINE_LOCK (lnum2);
	if (unlikely (ms & 0x1) &&
	    unlikely ((lnum3 = __cache_line_lookup (s3, ea_aligned3)) < 0))
	{
	    CACHE_STATS_RD_MISS (s3);
	    tagmask |= __cache_tagmask (s3);
	    lnum3 = __cache_rd_miss (ea3, s3);
	}
#ifdef CACHE_STATS
	else {
	    CACHE_STATS_RD_HIT (s3);
	}
#endif
	CACHELINE_UNLOCK (lnum0);
	CACHELINE_UNLOCK (lnum1);
	CACHELINE_UNLOCK (lnum2);

	spu_writech (MFC_WrTagMask, tagmask);
	spu_mfcstat (MFC_TAG_UPDATE_ALL);

	/* all 4 lines now in cache */

	d0 = *CACHE_ADDR (lnum0, spu_extract(ibyte, 0));
	d1 = *CACHE_ADDR (lnum1, spu_extract(ibyte, 1));
	d2 = *CACHE_ADDR (lnum2, spu_extract(ibyte, 2));
	d3 = *CACHE_ADDR (lnum3, spu_extract(ibyte, 3));

	ret = __load_vec_uint4 (d0, d1, d2, d3);
	return ret;
    }
#ifdef CACHE_STATS
    else {
	CACHE_STATS_RD_HIT (s0);
	CACHE_STATS_RD_HIT (s1);
	CACHE_STATS_RD_HIT (s2);
	CACHE_STATS_RD_HIT (s3);
    }
#endif
    return ret;
}
#endif /* CACHE_READ_X4 && nway >= 4 */

