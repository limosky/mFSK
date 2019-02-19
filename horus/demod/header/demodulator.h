/*
 * Copyright (C) January 2016 David Rowe, Brady O'Brien
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

#include <stdint.h>
#include <complex.h>

#include "fft.h"
#include "statistics.h"

#ifndef M_PI
#define M_PI    3.14159265358979f
#endif

#define TAU     (2.0f * M_PI)

#define MODE_2FSK 2
#define MODE_4FSK 4
#define MAX_TONES 4
#define FFT_SIZE  1024
#define FSK_SCALE 16383

int fsk_create(int Fs, int Rs, int M, int tx_f1);
int fsk_create_hbr(int Fs, int Rs, int P, int M, int tx_f1, int tx_fs);
void fsk_destroy(void);
void fsk_demod(uint8_t rx_bits[], complex float fsk_in[]);
void fsk_demod_sd(float rx_bits[], complex float fsk_in[]);

/* Getters/Setters */

int fsk_get_nin(void);
int fsk_get_N(void);
int fsk_get_Nmem(void);
int fsk_get_Nbits(void);
int fsk_get_Ts(void);
float fsk_get_f_est(int);
void fsk_set_nsym(int nsym);
void fsk_set_est_limits(int fmin, int fmax);
void fsk_set_estimators(void);

#ifdef __cplusplus
}
#endif
