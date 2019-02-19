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

#include <math.h>
#include <complex.h>

#include "statistics.h"
#include "demodulator.h"
#include "fft.h"

/* BSS memory */

static fft_cfg stats_fftcfg;
static float stats_snr_est; /* estimated SNR of rx signal in dB (3 kHz noise BW) */
static float stats_foff; /* estimated freq offset in Hz */
static float stats_rx_timing; /* estimated optimum timing offset in samples */
static float stats_clock_offset; /* Estimated tx/rx sample clock offset in ppm */
static float stats_f_est[MAX_TONES];
static float stats_scale_dB;
static float stats_fftbuffer[STATS_FFTSIZE];

/* Public Functions */

int stats_open() {
    if ((stats_fftcfg = fft_alloc(STATS_FFTSIZE, 0, NULL, NULL)) == NULL) {
        return 0;
    } else {
        int i;

        stats_snr_est = 0.0f;
        stats_foff = 0.0f;
        stats_rx_timing = 0.0f;
        stats_clock_offset = 0.0f;

        stats_scale_dB = 20.0f * log10f(ASCALE);

        for (i = 0; i < MAX_TONES; i++) {
            stats_f_est[i] = 0.0f;
        }

        for (i = 0; i < STATS_FFTSIZE; i++) {
            stats_fftbuffer[i] = 0.0f;
        }
    }

    return 1;
}

void stats_close() {
}

static float normf(complex float val) {
    float realf = crealf(val);
    float imagf = cimagf(val);
    return realf * realf + imagf * imagf;
}

/* Getters/Setters */

float stats_get_snr_est() {
    return stats_snr_est;
}

float stats_get_foff() {
    return stats_foff;
}

float stats_get_rx_timing() {
    return stats_rx_timing;
}

float stats_get_clock_offset() {
    return stats_clock_offset;
}

float stats_get_f_est(int index) {
    return stats_f_est[index];
}

float *stats_get_fft_buf_ptr() {
    return stats_fftbuffer;
}

void stats_set_snr_est(float val) {
    stats_snr_est = val;
}

void stats_set_foff(float val) {
    stats_foff = val;
}

void stats_set_rx_timing(float val) {
    stats_rx_timing = val;
}

void stats_set_clock_offset(float val) {
    stats_clock_offset = val;
}

void stats_set_f_est(int index, float val) {
    stats_f_est[index] = val;
}

void stats_spectrum(float mag_spec_dB[], complex float rx_cplx[], int len) {
    complex float fftin[STATS_FFTSIZE];
    complex float fftout[STATS_FFTSIZE];
    int i, j;

    /* shift buffer samples left by new length */

    for (i = 0; i < (STATS_FFTSIZE - len); i++) {
        stats_fftbuffer[i] = stats_fftbuffer[i + len];
    }

    /* now add in the new length of samples */

    for (j = 0; j < len; j++, i++) {
        stats_fftbuffer[i] = crealf(rx_cplx[j]);
    }

    /* Window and FFT */

    for (i = 0; i < STATS_FFTSIZE; i++) {
        fftin[i] = (stats_fftbuffer[i] * (0.5f - 0.5f * cosf(TAU * (float) i / STATS_FFTSIZE))) + 0.0f * I;
    }

    fft(stats_fftcfg, fftin, fftout);

    /* Scale and convert to dB */

    for (i = 0; i < (STATS_FFTSIZE / 2); i++) {
        mag_spec_dB[i] = 10.0f * log10f(normf(fftout[i]) + 1E-12f) - stats_scale_dB;
    }
}
