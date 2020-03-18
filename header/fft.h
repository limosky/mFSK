/*
 * Copyright (C) 2003-2004, Mark Borgerding
 *
 * All rights reserved.
 *
 * Modified by Steve Sampson for C99 complex math
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI    3.14159265358979f
#endif

#ifndef TAU
#define TAU     (2.0f * M_PI)
#endif
    
/* Complex FFT */

struct fft_state {
    int nfft;
    int inverse;
    int factors[64];
    complex float twiddles[1];
};

typedef struct fft_state *fft_cfg;

/* Real FFT */

struct fftr_state {
    fft_cfg substate;
    complex float *tmpbuf;
    complex float *super_twiddles;
};

typedef struct fftr_state *fftr_cfg;

/* Complex Function Calls */

fft_cfg fft_alloc(int, int, void *, size_t *);
void fft(fft_cfg, const complex float *, complex float *);

/* Real Function Calls */

fftr_cfg fftr_alloc(int, int, void *, size_t *);
void fftr(fftr_cfg, const float *, complex float *);
void fftri(fftr_cfg, const complex float *, float *);

#ifdef __cplusplus
}
#endif
