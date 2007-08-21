/* --------------------------------------------------------------- */
/* (C)Copyright 2006                                               */
/* International Business Machines Corporation,                    */
/* All Rights Reserved.                                            */
/*                                                                 */
/* This program is made available under the terms of the           */
/* Common Public License v1.0 which accompanies this distribution. */
/* --------------------------------------------------------------- */
/* PROLOG END TAG zYx                                              */
/*
 * cache-api.h
 *
 * Copyright (C) 2005, 2006 IBM Corp.
 *
 * Common header file for defining a cache.
 * May be included multiple times to declare
 * multiple caches as needed. 
 *
 * The callable interface from this file is:
 *
 *	cache_rd (name, ... 
 *	cache_wr (name, ... 
 *	cache_rw (name, ... 
 *	cache_touch (name, ... 
 *	cache_wait (name, ... 
 *	cache_flush (name, ... 
 *	cache_lock (name, ... 
 *	cache_unlock (name, ... 
 * #ifdef CACHE_READ_X4
 *	cache_rd_x4 (name, ...
 * #endif
 *
 * The following defines control the naming 
 * and size of various cache attributes and
 * may be redefined before each inclusion of
 * this file to define multiple caches with
 * different attributes:
 *
 *	Symbol			Default Value
 *	--------------------------------------
 *	CACHED_TYPE		none (required)
 * 	CACHE_NAME		none (required)
 *	CACHELINE_LOG2SIZE	7 (128 bytes)
 *	CACHE_LOG2NSETS		6 (64 sets)
 *	CACHE_LOG2NWAY		2 (4-way)
 *	CACHE_TYPE		1 (read/write)
 *	CACHE_SET_TAGID(set) 	(set & 0x1f)
 */

#ifndef _CACHE_API_H_
#define _CACHE_API_H_

#include <spu_mfcio.h>

#ifndef likely
#define likely(_c)      __builtin_expect((_c), 1)
#define unlikely(_c)    __builtin_expect((_c), 0)
#endif

/*
 * Per-cache name constructors.
 */
#define CONCAT(a, b)                    a ## b
#define STRINGIFY(name)                 CONCAT(name, _)
#define CACHE_PREFIX                    STRINGIFY(CACHE_NAME)
#define CACHE_APPEND(name, var)         CONCAT(name, var)
#define CACHE_STR(var)                  CACHE_APPEND(CACHE_PREFIX, var)
#define CACHE_FUNC(func)                CACHE_STR(func)
#define CACHE_VAR(var)                  CACHE_STR(var)
#define CACHE_NAME_STR(s)		#s
#define CACHE_GET_NAME(s)		CACHE_NAME_STR(s)
#define CACHE_NAME_STRING()		CACHE_GET_NAME(CACHE_NAME)

/* external interfaces */
#define cache_rd(name, ea)	(name ## _cache_rd((unsigned)(ea)))
#define cache_rw(name, ea)      (name ## _cache_rw((unsigned)(ea), 1))
#define cache_touch(name, ea)   (name ## _cache_rw((unsigned)(ea), 0))
#define cache_wr(name, ea, val) (name ## _cache_wr((unsigned)(ea), (val)))
#define cache_wait(name, ea)	(name ## _cache_wait((unsigned)(ea)))
#define cache_flush(name)	(name ## _cache_flush())
#define cache_lock(name, lsa)	(name ## _cache_lock((unsigned)(lsa)))
#define cache_unlock(name, lsa)	(name ## _cache_unlock((unsigned)(lsa)))
#ifdef CACHE_READ_X4
#define cache_rd_x4(name, ea)	(name ## _cache_rd_x4((ea)))
#endif

#define CACHE_TYPE_RO 0
#define CACHE_TYPE_RW 1

#endif /* end defined once stuff */

#ifndef CACHED_TYPE
#error CACHED_TYPE not defined!!
#endif

#ifndef CACHE_NAME
#error CACHED_NAME not defined!!
#endif

#ifndef CACHE_TYPE
#define CACHE_TYPE CACHE_TYPE_RW
#endif 


#if (CACHE_TYPE != CACHE_TYPE_RO && CACHE_TYPE != CACHE_TYPE_RW)
#error Invalid CACHE_TYPE!
#endif

/*
 * Log2 of cache line size.
 * Default is 7 (128 bytes).
 */
#ifndef CACHELINE_LOG2SIZE
#define CACHELINE_LOG2SIZE 7
#else
#if (CACHELINE_LOG2SIZE < 4 || CACHELINE_LOG2SIZE > 12)
#error Invalid CACHELINE_LOG2SIZE!
#endif
#endif

#define CACHELINE_SIZE	(1 << CACHELINE_LOG2SIZE)
#define CACHELINE_SHIFT	CACHELINE_LOG2SIZE
#define CACHE_ALIGN_MASK  ~(CACHELINE_SIZE - 1)

/*
 * Log2 of number of cached sets.
 * Default is 6 (64 sets).
 */
#ifndef CACHE_LOG2NSETS
#define CACHE_LOG2NSETS 6
#else
#if (CACHE_LOG2NSETS < 0 || CACHE_LOG2NSETS > 12)
#error Invalid CACHE_LOG2NSETS!
#endif
#endif

#define CACHE_NSETS		(1 << CACHE_LOG2NSETS)
#define CACHE_NSETS_MASK	(CACHE_NSETS - 1)
#define CACHE_NSETS_SHIFT	CACHE_LOG2NSETS

/*
 * Log2 of associativity.
 * Default is 2 (4-way).
 */
#ifndef CACHE_LOG2NWAY
#define CACHE_LOG2NWAY 2
#else
#if (CACHE_LOG2NWAY != 0 && CACHE_LOG2NWAY != 2)
#error Invalid CACHE_LOG2NWAY!
#endif
#endif

#define CACHE_NWAY		(1 << CACHE_LOG2NWAY)
#define CACHE_NWAY_MASK		(CACHE_NWAY - 1)
#define CACHE_NWAY_SHIFT	CACHE_LOG2NWAY

/*
 * TAGID definitions 
 */
#ifndef CACHE_SET_TAGID
#define CACHE_SET_TAGID(set)        ((set) & 0x1f)
#endif

/* the size of the cache */
#define CACHE_MEM_SIZE      (CACHE_NSETS*CACHE_NWAY*CACHELINE_SIZE)

/* the cache */
#define __cache_mem	CACHE_VAR (cache_mem)
char __cache_mem[CACHE_MEM_SIZE] __attribute__ ((aligned (128)));

/* the cache directory */
#include "cache-dir.h"

/* internal stats */
#include "cache-stats.h"

/*
 * include the appropriate cache implementation
 * (only 4-way and direct mapped for now)
 */

#if (CACHE_NWAY == 1 || CACHE_NWAY == 4)
#include "cache-4way.h"
#else
#error Invalid CACHE_NWAY!
#endif

#undef CACHED_TYPE
#undef CACHE_TYPE
#undef CACHE_NAME
#undef CACHE_LOG2NWAY
#undef CACHE_LOG2NSETS
#undef CACHELINE_LOG2SIZE

