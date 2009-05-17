/* --------------------------------------------------------------  */
/* (C)Copyright 2006,2007,                                         */
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
/* cache-stats.h
 *
 * Internal statistics for
 * software managed cache.
 */

#ifndef _CACHE_STATS_H_
#define _CACHE_STATS_H_

#define __cache_stats	 CACHE_VAR(cache_stats)
#define __cache_pr_stats CACHE_FUNC(cache_pr_stats)

typedef unsigned int cache_cntr_t;

struct cache_metrics {
    cache_cntr_t ld_hits;
    cache_cntr_t st_hits;
    cache_cntr_t ld_misses;
    cache_cntr_t st_misses;
    cache_cntr_t flushes;
    cache_cntr_t replace;
    cache_cntr_t writebacks;
};

#define CACHE_STATS_BANNER0 \
		"\nname=%s type=%s nway=%d nsets=%d linesz=%d totsz=%d\n"
#define CACHE_STATS_BANNER1 \
		"\n                 READ                      WRITE" \
		"\n         --------------------      --------------------\n"
#define CACHE_STATS_BANNER2 \
		"SET      HIT     MISS  HITPCT      HIT     MISS  HITPCT" \
		"     REPL     WRBK FLUSH\n"
#define CACHE_STATS_FORMATTER \
		"%3d%9d%9d%7.1f%%%9d%9d%7.1f%%%9d%9d%6d\n"
#define CACHE_STATS_TOT_FORMATTER \
		"TOT%9d%9d%7.1f%%%9d%9d%7.1f%%%9d%9d%6d\n"

#endif /* end defined once stuff */

/*
 * these are undefed and redefined so they can
 * be called without ifdefing the calling code
 */
#undef CACHE_STATS_RD_HIT
#undef CACHE_STATS_WR_HIT
#undef CACHE_STATS_RD_MISS
#undef CACHE_STATS_WR_MISS
#undef CACHE_STATS_REPLACE
#undef CACHE_STATS_WRITEBACK
#undef CACHE_STATS_FLUSH
#undef cache_pr_stats

#ifdef CACHE_STATS

#define CACHE_STATS_RD_HIT(set)     (__cache_stats[(set)].ld_hits++)
#define CACHE_STATS_WR_HIT(set)     (__cache_stats[(set)].st_hits++)
#define CACHE_STATS_RD_MISS(set)    (__cache_stats[(set)].ld_misses++)
#define CACHE_STATS_WR_MISS(set)    (__cache_stats[(set)].st_misses++)
#define CACHE_STATS_REPLACE(set)    (__cache_stats[(set)].replace++)
#define CACHE_STATS_WRITEBACK(set)  (__cache_stats[(set)].writebacks++)
#define CACHE_STATS_FLUSH(set)	    (__cache_stats[(set)].flushes++)

#define cache_pr_stats(name)	    (name ## _cache_pr_stats())

static struct cache_metrics * __cache_stats;

/* FIXME: ALL INDIVIDUAL CACHES MUST BE USE STATS OR NOT USE STATS IN
   A FILE DUE TO HACKERY SURROUNDING DECLARIONS AND SO FORTH. */
#define __cache_stats_init(name) /* Safe outside these headers */       \
  struct cache_metrics CACHE_SYM(name,stack_stats)                      \
    [ CACHE_SYM(name,cache_nsets) ];                                    \
  CACHE_SYM(name,cache_dir) = CACHE_SYM(name,stack_stats)

static inline void
__cache_pr_stats(void)
{
    cache_cntr_t ld_total, ld_hits, ld_misses;
    cache_cntr_t ld_total_tot = 0, ld_hits_tot = 0,
		 ld_misses_tot = 0;
    float pld_hits, pld_hits_tot;
    cache_cntr_t st_total = 0, st_hits = 0, st_misses = 0;
    cache_cntr_t st_total_tot = 0, st_hits_tot = 0,
		 st_misses_tot = 0;
    float pst_hits = 0.0, pst_hits_tot = 0.0;
    cache_cntr_t writebacks, replace, flushes;
    cache_cntr_t writebacks_tot = 0, replace_tot = 0, flushes_tot = 0;
    int i;

    printf(CACHE_STATS_BANNER0, CACHE_NAME_STRING(),
#if (CACHE_TYPE == CACHE_TYPE_RW)
	    "READ_WRITE",
#else
	    "READ_ONLY",
#endif
	    CACHE_NWAY, CACHE_NSETS, CACHELINE_SIZE, CACHE_MEM_SIZE
	);

    printf(CACHE_STATS_BANNER1);
    printf(CACHE_STATS_BANNER2);

    for (i = 0; i < CACHE_NSETS; i++) {
	ld_hits = __cache_stats[i].ld_hits;
	ld_hits_tot += ld_hits;
	ld_misses = __cache_stats[i].ld_misses;
	ld_misses_tot += ld_misses;
	ld_total = ld_hits + ld_misses;
	ld_total_tot += ld_total;
	pld_hits = (ld_total) ? (double) ld_hits / (double) ld_total : 0.0f;
	pld_hits *= 100.0f;
        st_total = 0;
#if (CACHE_TYPE == CACHE_TYPE_RW)
	st_hits = __cache_stats[i].st_hits;
	st_hits_tot += st_hits;
	st_misses = __cache_stats[i].st_misses;
	st_misses_tot += st_misses;
	st_total = st_hits + st_misses;
	st_total_tot += st_total;
	pst_hits = (st_total) ? (double) st_hits / (double) st_total : 0.0f;
	pst_hits *= 100.0f;
#endif
	writebacks = __cache_stats[i].writebacks;
	writebacks_tot += writebacks;
	replace = __cache_stats[i].replace;
	replace_tot += replace;
	flushes = __cache_stats[i].flushes;
	flushes_tot += flushes;

	printf(CACHE_STATS_FORMATTER, i,
		ld_hits, ld_misses, pld_hits,
		st_hits, st_misses, pst_hits,
		writebacks, replace, flushes);
    }

    /* totals */
    pld_hits_tot = (ld_total_tot) ? (double) ld_hits_tot /
		(double) ld_total_tot : 0.0f;
    pld_hits_tot *= 100.0f;
    pst_hits_tot = (st_total_tot) ? (double) st_hits_tot /
		(double) st_total_tot : 0.0f;
    pst_hits_tot *= 100.0f;

    printf(CACHE_STATS_TOT_FORMATTER,
	    ld_hits_tot, ld_misses_tot, pld_hits_tot,
	    st_hits_tot, st_misses_tot, pst_hits_tot,
	    writebacks_tot, replace_tot, flushes_tot);

}
#else
#define __cache_stats_init(name)
#define CACHE_STATS_RD_HIT(set)
#define CACHE_STATS_WR_HIT(set)
#define CACHE_STATS_RD_MISS(set)
#define CACHE_STATS_WR_MISS(set)
#define CACHE_STATS_REPLACE(set)
#define CACHE_STATS_WRITEBACK(set)
#define CACHE_STATS_FLUSH(set)
#define cache_pr_stats(name)
#endif
