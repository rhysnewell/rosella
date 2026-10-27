#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <immintrin.h>

#include "easel.h"
#include "hmmer.h"
#include "impl_sse.h"

#include "hmmer_avx2.h"

// Knudsen's band layout from HMMER's ssvfilter.c at 32 lanes. Fourteen bands leave two of the
// sixteen ymm registers for the running maximum and the begin vector.
#define LANES 32
#define BANDS 14

struct rosella_ssv32 {
  int M;
  int Q;
  __m256i *mem;
  __m256i **sbv;
};

static int segments(int M) { return ESL_MAX(2, (M - 1) / LANES + 1); }

rosella_ssv32 *rosella_ssv32_new(const P7_OPROFILE *om) {
  int M = om->M, Q = segments(M), Q16 = p7O_NQB(M), Kp = om->abc->Kp, width = Q + BANDS;
  rosella_ssv32 *held = calloc(1, sizeof(*held));
  if (held == NULL) return NULL;
  held->M = M;
  held->Q = Q;
  held->sbv = malloc(sizeof(__m256i *) * Kp);
  held->mem = aligned_alloc(32, sizeof(__m256i) * width * Kp);
  if (held->sbv == NULL || held->mem == NULL) {
    rosella_ssv32_free(held);
    return NULL;
  }
  for (int x = 0; x < Kp; x++) {
    uint8_t striped[16];
    uint8_t lanes[LANES];
    held->sbv[x] = held->mem + x * width;
    for (int q = 0; q < Q; q++) {
      for (int z = 0; z < LANES; z++) {
        int k = 1 + q + z * Q;
        lanes[z] = 127;
        if (k <= M) {
          memcpy(striped, &om->sbv[x][(k - 1) % Q16], 16);
          lanes[z] = striped[(k - 1) / Q16];
        }
      }
      memcpy(&held->sbv[x][q], lanes, LANES);
    }
    for (int q = Q; q < width; q++) held->sbv[x][q] = held->sbv[x][q % Q];
  }
  return held;
}

void rosella_ssv32_free(rosella_ssv32 *held) {
  if (held == NULL) return;
  free(held->mem);
  free(held->sbv);
  free(held);
}

static inline __m256i wrapped(__m256i sv, __m256i beginv) {
  __m256i low = _mm256_permute2x128_si256(sv, sv, 0x08);
  return _mm256_or_si256(_mm256_alignr_epi8(sv, low, 15), beginv);
}

#define STEP(W)                                     \
  _Pragma("GCC unroll 16") for (int b = 0; b < W; b++) { \
    sv[b] = _mm256_subs_epi8(sv[b], rsc[b]);        \
    xEv = _mm256_max_epu8(xEv, sv[b]);              \
  }

#define TURN(W, at, check, label)                                 \
  _Pragma("GCC unroll 16") for (int j = W - 1; j >= 0; j--) {     \
    if (check && i >= L) goto label;                              \
    rsc = held->sbv[dsq[at]] + Q - (j + 1);                       \
    STEP(W)                                                       \
    sv[j] = wrapped(sv[j], beginv);                               \
    at++;                                                         \
  }

#define BAND(W)                                                                        \
  static __m256i band_##W(const ESL_DSQ *dsq, int L, const rosella_ssv32 *held, int q, \
                          __m256i beginv, __m256i xEv) {                               \
    __m256i sv[W];                                                                     \
    const __m256i *rsc;                                                                \
    int i, i2, Q = held->Q;                                                            \
    dsq++;                                                                             \
    _Pragma("GCC unroll 16") for (int b = 0; b < W; b++) sv[b] = beginv;               \
    for (i = 0; i < L && i < Q - q - W; i++) {                                         \
      rsc = held->sbv[dsq[i]] + i + q;                                                 \
      STEP(W)                                                                          \
    }                                                                                  \
    i = Q - q - W;                                                                     \
    TURN(W, i, 1, first)                                                               \
  first:                                                                               \
    for (i2 = Q - q; i2 < L - Q; i2 += Q) {                                            \
      for (i = 0; i < Q - W; i++) {                                                    \
        rsc = held->sbv[dsq[i2 + i]] + i;                                              \
        STEP(W)                                                                        \
      }                                                                                \
      i += i2;                                                                         \
      TURN(W, i, 0, first)                                                             \
    }                                                                                  \
    for (i = 0; i2 + i < L && i < Q - W; i++) {                                        \
      rsc = held->sbv[dsq[i2 + i]] + i;                                                \
      STEP(W)                                                                          \
    }                                                                                  \
    i += i2;                                                                           \
    TURN(W, i, 1, last)                                                                \
  last:                                                                                \
    return xEv;                                                                        \
  }

BAND(1)
BAND(2)
BAND(3)
BAND(4)
BAND(5)
BAND(6)
BAND(7)
BAND(8)
BAND(9)
BAND(10)
BAND(11)
BAND(12)
BAND(13)
BAND(14)

typedef __m256i (*band_fn)(const ESL_DSQ *, int, const rosella_ssv32 *, int, __m256i, __m256i);

static const band_fn bands[BANDS + 1] = {NULL,     band_1,  band_2,  band_3,  band_4,
                                         band_5,   band_6,  band_7,  band_8,  band_9,
                                         band_10,  band_11, band_12, band_13, band_14};

uint8_t rosella_ssv32_xe(const ESL_DSQ *dsq, int L, const rosella_ssv32 *held) {
  __m256i beginv = _mm256_set1_epi8(-128);
  __m256i xEv = beginv;
  int Q = held->Q, rounds = (Q + BANDS - 1) / BANDS, last = 0;
  for (int round = 0; round < rounds; round++) {
    int q = (Q * (round + 1)) / rounds;
    xEv = bands[q - last](dsq, L, held, last, beginv, xEv);
    last = q;
  }
  __m128i half = _mm_max_epu8(_mm256_castsi256_si128(xEv), _mm256_extracti128_si256(xEv, 1));
  half = _mm_max_epu8(half, _mm_srli_si128(half, 8));
  half = _mm_max_epu8(half, _mm_srli_si128(half, 4));
  half = _mm_max_epu8(half, _mm_srli_si128(half, 2));
  half = _mm_max_epu8(half, _mm_srli_si128(half, 1));
  return (uint8_t) _mm_extract_epi8(half, 0);
}
