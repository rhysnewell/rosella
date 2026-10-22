/* easel/esl_config.h.  Generated from esl_config.h.in by configure.  */
/* esl_config.h.in  [input to configure]
 * 
 * System-dependent configuration of Easel, by autoconf.
 * 
 * This file should be included in all Easel .c files before
 * anything else, because it may set #define's that control
 * behaviour of system includes and system libraries. An example
 * is large file support.
 * 
 */
#ifndef eslCONFIG_INCLUDED
#define eslCONFIG_INCLUDED

#if !defined(__x86_64__) && !defined(__aarch64__)
#error "rosella builds HMMER for x86_64 and aarch64 only"
#endif

/* Version info.
 */
#define EASEL_VERSION "0.49"
#define EASEL_DATE "Aug 2023"
#define EASEL_COPYRIGHT "Copyright (C) 2023 Howard Hughes Medical Institute."
#define EASEL_LICENSE "Freely distributed under the BSD open source license."
#define EASEL_URL "http://bioeasel.org/"

/* Control of debugging instrumentation */
#define eslDEBUGLEVEL 0
/* #undef eslENABLE_ASAN */
/* #undef eslENABLE_TSAN */

/* Optional parallel implementation support */
#if defined(__x86_64__)
#define eslENABLE_SSE 1
#endif
/* #undef eslENABLE_SSE4 */
/* #undef eslENABLE_AVX */
/* #undef eslENABLE_AVX512 */
#if defined(__aarch64__)
#define eslENABLE_NEON 1
#endif
/* #undef eslENABLE_VMX */

#if defined(__aarch64__)
#define eslHAVE_NEON_AARCH64 1
#endif

/* #undef eslENABLE_CUDA */

#if defined(__x86_64__)
#define HAVE_FLUSH_ZERO_MODE 1
#endif
/* #undef HAVE_DENORMALS_ZERO_MODE */

/* #undef HAVE_MPI */
#define HAVE_PTHREAD 1

/* Programs */
#define HAVE_GZIP 1

/* Libraries */
/* #undef HAVE_LIBGSL */

/* Headers */
#if defined(__linux__)
#define HAVE_ENDIAN_H 1
#endif
#define HAVE_INTTYPES_H 1
#define HAVE_STDINT_H 1
#define HAVE_UNISTD_H 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_STRINGS_H 1
#define HAVE_NETINET_IN_H 1

#define HAVE_SYS_PARAM_H 1
#if defined(__APPLE__)
#define HAVE_SYS_SYSCTL_H 1
#endif

/* Types */
/* #undef WORDS_BIGENDIAN */
/* #undef int8_t */
/* #undef int16_t */
/* #undef int32_t */
/* #undef int64_t */
/* #undef uint8_t */
/* #undef uint16_t */
/* #undef uint32_t */
/* #undef uint64_t */
/* #undef off_t */

/* Compiler characteristics */
#define HAVE_FUNC_ATTRIBUTE_NORETURN 1
#define HAVE_FUNC_ATTRIBUTE_FORMAT 1

/* Functions */
/* #undef HAVE_ALIGNED_ALLOC */
#if defined(__APPLE__)
#define HAVE_ERFC 1
#endif
#define HAVE_GETCWD 1
#define HAVE_GETPID 1
/* #undef HAVE__MM_MALLOC */
#define HAVE_POPEN 1
/* #undef HAVE_POSIX_MEMALIGN */
#define HAVE_STRCASECMP 1
#define HAVE_STRSEP 1
#define HAVE_SYSCONF 1
#if defined(__APPLE__)
#define HAVE_SYSCTL 1
#endif
#define HAVE_TIMES 1

#define HAVE_FSEEKO 1

/* System services */
/* #undef _FILE_OFFSET_BITS */
/* #undef _LARGE_FILES */
/* #undef _LARGEFILE_SOURCE */

 
/* Function behavior */
#define eslSTOPWATCH_HIGHRES

#endif /*eslCONFIG_INCLUDED*/

