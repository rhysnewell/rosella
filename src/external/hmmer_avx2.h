#ifndef ROSELLA_HMMER_AVX2_H
#define ROSELLA_HMMER_AVX2_H

#include <stdint.h>

#include "easel.h"
#include "hmmer.h"

typedef struct rosella_ssv32 rosella_ssv32;

rosella_ssv32 *rosella_ssv32_new(const P7_OPROFILE *om);
void rosella_ssv32_free(rosella_ssv32 *held);
uint8_t rosella_ssv32_xe(const ESL_DSQ *dsq, int L, const rosella_ssv32 *held);

#endif
