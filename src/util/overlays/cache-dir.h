/* --------------------------------------------------------------  */
/* (C)Copyright 2006                                               */
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
 * cache-dir.h
 *
 * Implementation of cache directory
 *
 * Included by top-level cache header file
 */

#ifndef _CACHE_DIR_H_
#define _CACHE_DIR_H_

#define __cache_dir             CACHE_VAR (cache_dir)
#define __cache_lock            CACHE_FUNC (cache_lock)
#define __cache_unlock          CACHE_FUNC (cache_unlock)

/*
 * The smallest supported line size is 2^4
 * so the 4 lowest order bits of the tag
 * are used as state bits for the line.
 */
#define CACHELINE_VALID		0x1
#define CACHELINE_MODIFIED	0x2
#define CACHELINE_LOCKED	0x4

#define CACHELINE_TAG_MASK CACHE_ALIGN_MASK

/* lock bit */
#define CACHELINE_LOCK(lnum)      (__cache_dir[(lnum)] |= CACHELINE_LOCKED)
#define CACHELINE_UNLOCK(lnum)    (__cache_dir[(lnum)] &= ~CACHELINE_LOCKED)
#define CACHELINE_ISLOCKED(lnum)   (__cache_dir[(lnum)] & CACHELINE_LOCKED)

/* mod bit */
#define CACHELINE_SETMOD(lnum)    (__cache_dir[(lnum)] |= CACHELINE_MODIFIED)
#define CACHELINE_CLEARMOD(lnum)  (__cache_dir[(lnum)] &= ~CACHELINE_MODIFIED)
#define CACHELINE_ISMOD(lnum)      \
    ((__cache_dir[(lnum)] & (CACHELINE_MODIFIED | CACHELINE_VALID)) == \
        (CACHELINE_MODIFIED | CACHELINE_VALID))

/* valid bit */
#define CACHELINE_SETVALID(lnum)  (__cache_dir[(lnum)] |= CACHELINE_VALID)
#define CACHELINE_CLEARVALID(lnum) (__cache_dir[(lnum)] &= ~CACHELINE_VALID)

/* tag */
#define CACHELINE_GETTAG(lnum)   (__cache_dir[(lnum)] & CACHELINE_TAG_MASK)
#define CACHELINE_SETTAG(lnum, ea) \
    (__cache_dir[(lnum)] = (((ea) & CACHELINE_TAG_MASK) | CACHELINE_VALID))
#define CACHELINE_SETTAGMOD(lnum, ea) \
    (__cache_dir[(lnum)] = \
     (((ea) & CACHELINE_TAG_MASK) | CACHELINE_VALID | CACHELINE_MODIFIED))
#define CACHELINE_ISTAG(lnum, ea) \
    ((__cache_dir[(lnum)] & (CACHELINE_TAG_MASK|CACHELINE_VALID)) == \
        (((ea) & CACHELINE_TAG_MASK) | CACHELINE_VALID))

#define CACHEDIR_ADDR(lnum) (&__cache_dir[(lnum)])

#define CACHE_ADDR2LNUM(a) \
    (((unsigned)(a) - (unsigned)&__cache_mem[0]) >> CACHELINE_LOG2SIZE)

#endif

/*
 * The cache directory
 */
static unsigned * RESTRICT __attribute__ ((aligned (16))) __cache_dir;

#define __cache_dir_init(name) /* Safe outside these headers */             \
  SPU_MALLOC( CACHE_SYM(name,cache_dir),                                    \
              CACHE_SYM(name,cache_nsets)*CACHE_SYM(name,cache_nway), 16 ); \
  CLEAR( CACHE_SYM(name,cache_dir), /* Mark cache empty */                  \
         CACHE_SYM(name,cache_nsets)*CACHE_SYM(name,cache_nway) )           \

static inline void
 __cache_lock(unsigned lsa)
{
    CACHELINE_LOCK(CACHE_ADDR2LNUM((lsa)));
    return;
}
static inline void
 __cache_unlock(unsigned lsa)
{
    CACHELINE_UNLOCK(CACHE_ADDR2LNUM((lsa)));
    return;
}
