/*
 * Copyright (C) January 2016 David Rowe
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

#include <stdint.h>
#include <stdbool.h>
#include <complex.h>

#include "fft.h"
#include "statistics.h"

#define MODE_2FSK 2
#define MODE_4FSK 4
#define MAX_TONES 4
#define FFT_SIZE  1024
#define FSK_SCALE 16383
    
#define OVERSAMPLE_RATE    8
#define MIN_FREQUENCY    300
#define MAX_FREQUENCY   2200
#define MIN_SPACING      100

#define cmplx(value) (cosf(value) + sinf(value) * I)
#define cmplxconj(value) (cosf(value) + sinf(value) * -I)

struct FSK {
    struct STATS *stats;
    complex float phase[MAX_TONES];
    float f_est[MAX_TONES];
    float norm_rx_timing;
    float EbNodB;
    float ppm;
    int Ndft;
    int Fs;
    int N;
    int Rs;
    int Ts;
    int Nmem;
    int P;
    int Nsym;
    int Nbits;
    int f1;
    int shift;
    int mode;
    int est_min;
    int est_max;
    int est_space;
    int nstash;
    int nin;
    bool burst_mode;
    fft_cfg fftcfg;
    float *hann_table;
    float *fft_est;
    complex float *samp_old;
};

struct FSK *fsk_create(int, int, int, int);
struct FSK *fsk_create_hbr(int, int, int, int, int, int);
void fsk_destroy(struct FSK *);

/* Getters/Setters */

int fsk_get_nin(struct FSK *);
int fsk_get_N(struct FSK *);
int fsk_get_Nmem(struct FSK *);
int fsk_get_Nbits(struct FSK *);
int fsk_get_Ts(struct FSK *);
float fsk_get_f_est(struct FSK *, int);
void fsk_set_nsym(struct FSK *, int);
void fsk_set_est_limits(struct FSK *, int, int);
void fsk_set_estimators(struct FSK *);

#ifdef __cplusplus
}
#endif

