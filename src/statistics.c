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
#include "fft.h"
#include "fsk.h"

/* Private Functions */

static float normf(complex float val) {
    float realf = crealf(val);
    float imagf = cimagf(val);
    return realf * realf + imagf * imagf;
}

/* Public Functions */

struct STATS *stats_open() {
    struct STATS *stats;

    if ((stats = (struct STATS *) malloc(sizeof (struct STATS))) == NULL) {
        return NULL;
    }

    if ((stats->fftcfg = fft_alloc(STATS_FFTSIZE, 0, NULL, NULL)) == NULL) {
        return NULL;
    } else {
        int i;

        stats->snr_est = 0.0f;
        stats->foff = 0.0f;
        stats->rx_timing = 0.0f;
        stats->clock_offset = 0.0f;

        stats->scale_dB = 20.0f * log10f(ASCALE);

        for (i = 0; i < STATS_FFTSIZE; i++) {
            stats->fftbuffer[i] = 0.0f;
        }
    }

    return stats;
}

void stats_close(struct STATS *stats) {
    free(stats->fftcfg);
    free(stats);
}

/* Getters/Setters */

float stats_get_snr_est(struct STATS *stats) {
    return stats->snr_est;
}

float stats_get_foff(struct STATS *stats) {
    return stats->foff;
}

float stats_get_rx_timing(struct STATS *stats) {
    return stats->rx_timing;
}

float stats_get_clock_offset(struct STATS *stats) {
    return stats->clock_offset;
}

float *stats_get_fft_buf_ptr(struct STATS *stats) {
    return stats->fftbuffer;
}

void stats_set_snr_est(struct STATS *stats, float val) {
    stats->snr_est = val;
}

void stats_set_foff(struct STATS *stats, float val) {
    stats->foff = val;
}

void stats_set_rx_timing(struct STATS *stats, float val) {
    stats->rx_timing = val;
}

void stats_set_clock_offset(struct STATS *stats, float val) {
    stats->clock_offset = val;
}

void stats_spectrum(struct STATS *stats, float mag_spec_dB[], complex float rx_cplx[], int len) {
    complex float fftin[STATS_FFTSIZE];
    complex float fftout[STATS_FFTSIZE];
    int i, j;

    /* shift buffer samples left by new length */

    for (i = 0; i < (STATS_FFTSIZE - len); i++) {
        stats->fftbuffer[i] = stats->fftbuffer[i + len];
    }

    /* now add in the new length of samples */

    for (j = 0; j < len; j++, i++) {
        stats->fftbuffer[i] = crealf(rx_cplx[j]);
    }

    /* Window and FFT */

    for (i = 0; i < STATS_FFTSIZE; i++) {
        fftin[i] = (stats->fftbuffer[i] * (0.5f - 0.5f * cosf(TAU * (float) i / STATS_FFTSIZE))) + 0.0f * I;
    }

    fft(stats->fftcfg, fftin, fftout);

    /* Scale and convert to dB */

    for (i = 0; i < (STATS_FFTSIZE / 2); i++) {
        mag_spec_dB[i] = 10.0f * log10f(normf(fftout[i]) + 1E-12f) - stats->scale_dB;
    }
}
