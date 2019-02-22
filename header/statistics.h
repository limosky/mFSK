/*
 * Copyright (C) June 2015 David Rowe
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

#include "fsk.h"
#include "fft.h"

#define STATS_FFTSIZE   1024
#define ASCALE  ((STATS_FFTSIZE / 2) * 1000)

struct STATS {
    fft_cfg fftcfg;
    float snr_est; /* estimated SNR of rx signal in dB (3 kHz noise BW) */
    float foff; /* estimated freq offset in Hz */
    float rx_timing; /* estimated optimum timing offset in samples */
    float clock_offset; /* Estimated tx/rx sample clock offset in ppm */
    float scale_dB;
    float fftbuffer[STATS_FFTSIZE];
};

struct STATS *stats_open(void);
void stats_close(struct STATS *);
void stats_spectrum(struct STATS *, float [], complex float [], int);

float stats_get_snr_est(struct STATS *);
float stats_get_foff(struct STATS *);
float stats_get_rx_timing(struct STATS *);
float stats_get_clock_offset(struct STATS *);

float *stats_get_fft_buf_ptr(struct STATS *);

void stats_set_snr_est(struct STATS *, float);
void stats_set_foff(struct STATS *, float);
void stats_set_rx_timing(struct STATS *, float);
void stats_set_clock_offset(struct STATS *, float);

#ifdef __cplusplus
}
#endif
