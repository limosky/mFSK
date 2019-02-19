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
extern "C"
{
#endif

#include <complex.h>

#define STATS_FFTSIZE   1024
#define ASCALE  ((STATS_FFTSIZE / 2) * 1000)

int stats_open(void);
void stats_close(void);
void stats_spectrum(float [], complex float [], int);

float stats_get_snr_est(void);
float stats_get_foff(void);
float stats_get_rx_timing(void);
float stats_get_clock_offset(void);
float stats_get_f_est(int);
float stats_get_fft_buf(int);

float *stats_get_fft_buf_ptr(void);

void stats_set_snr_est(float);
void stats_set_foff(float);
void stats_set_rx_timing(float);
void stats_set_clock_offset(float);
void stats_set_f_est(int, float);
void stats_set_fft_buf(int, float);

#ifdef __cplusplus
}
#endif
