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
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <complex.h>
#include <math.h>

#include "fsk.h"
#include "demodulator.h"
#include "statistics.h"
#include "fft.h"

static int generate_hann_table(struct FSK *);

struct FSK *fsk_create(int sample_rate, int symbol_rate, int mfsk, int tone_f1) {
    struct FSK *fsk;
    struct STATS *stats;

    if ((fsk = (struct FSK *) malloc(sizeof (struct FSK))) == NULL) {
        return NULL;
    }
    
    if ((stats = stats_open()) == NULL) {
        return NULL;
    }

    fsk->stats = stats;
    fsk->Fs = sample_rate;
    fsk->Rs = symbol_rate;
    fsk->Ts = sample_rate / symbol_rate;
    fsk->N = sample_rate;
    fsk->burst_mode = false;
    fsk->P = OVERSAMPLE_RATE;
    fsk->Nsym = fsk->N / fsk->Ts;
    fsk->Ndft = FFT_SIZE;
    fsk->Nmem = fsk->N + (fsk->Ts * 2);
    fsk->f1 = tone_f1;
    fsk->shift = (symbol_rate * 2);
    fsk->nin = fsk->N;
    fsk->mode = mfsk;
    fsk->Nbits = (mfsk == MODE_2FSK) ? fsk->Nsym : fsk->Nsym * 2;

    fsk->est_min = MIN_FREQUENCY;
    fsk->est_max = MAX_FREQUENCY;
    fsk->est_space = MIN_SPACING;

    /* Set up rx state */
    for (int i = 0; i < mfsk; i++) {
        fsk->phase[i] = cmplx(0.0f);
    }

    fsk->nstash = (fsk->Ts * 4);

    if (generate_hann_table(fsk) == 0) {
        return NULL;
    }

    if ((fsk->samp_old = (complex float *) malloc(sizeof (complex float) * fsk->nstash)) == NULL) {
        return NULL;
    }

    for (int i = 0; i < fsk->nstash; i++) {
        fsk->samp_old[i] = 0.0f;
    }

    if ((fsk->fftcfg = fft_alloc(fsk->Ndft, 0, NULL, NULL)) == NULL) {
        free(fsk->samp_old);
        
        return NULL;
    }

    if ((fsk->fft_est = (float *) malloc(sizeof (float) * (fsk->Ndft / 2))) == NULL) {
        free(fsk->samp_old);
        free(fsk->fftcfg);
        
        return NULL;
    }

    for (int i = 0; i < (fsk->Ndft / 2); i++) {
        fsk->fft_est[i] = 0.0f;
    }

    for (int i = 0; i < mfsk; i++) {
        fsk->f_est[i] = 0.0f;
    }

    fsk->norm_rx_timing = 0.0f;
    fsk->EbNodB = 0.0f;
    fsk->ppm = 0.0f;

    if (stats_open() == 0) {
        free(fsk->fft_est);
        free(fsk->samp_old);
        free(fsk->fftcfg);
        
        return NULL;
    }

    return fsk;
}

struct FSK *fsk_create_hbr(int sample_rate, int symbol_rate, int oversample_rate, int mfsk, int f1, int shift) {
    struct FSK *fsk;
    struct STATS *stats;

    if ((fsk = (struct FSK *) malloc(sizeof (struct FSK))) == NULL) {
        return NULL;
    }

    if ((stats = stats_open()) == NULL) {
        return NULL;
    }

    fsk->stats = stats;
    fsk->Fs = sample_rate;
    fsk->Rs = symbol_rate;
    fsk->Ts = sample_rate / symbol_rate; /* cycles per symbol */
    fsk->burst_mode = false;
    fsk->Nsym = 48; /* Number of symbols in a processing frame */
    fsk->N = fsk->Ts * fsk->Nsym;
    fsk->P = oversample_rate;
    fsk->Nmem = fsk->N + (fsk->Ts * 2);
    fsk->f1 = f1;
    fsk->shift = shift;
    fsk->nin = fsk->N;
    fsk->mode = mfsk;
    fsk->Nbits = (mfsk == MODE_2FSK) ? fsk->Nsym : fsk->Nsym * 2;

    int dft = 0;

    for (int i = 1; i; i <<= 1)
        if (fsk->N & i)
            dft = i;

    fsk->Ndft = dft;

    fsk->est_min = (fsk->Rs / 2);
    fsk->est_max = (fsk->Fs / 2) - fsk->Rs;
    fsk->est_space = fsk->Rs - (fsk->Rs / 5);

    /* Set up rx state */

    for (int i = 0; i < mfsk; i++) {
        fsk->phase[i] = cmplx(0.0f);
    }

    fsk->nstash = (fsk->Ts * 4);

    if (generate_hann_table(fsk) == 0) {
        return NULL;
    }

    if ((fsk->fftcfg = fft_alloc(fsk->Ndft, 0, NULL, NULL)) == NULL) {
        return NULL;
    }

    if ((fsk->samp_old = (complex float *) malloc(sizeof (complex float) * fsk->nstash)) == NULL) {
        return NULL;
    }

    for (int i = 0; i < fsk->nstash; i++) {
        fsk->samp_old[i] = 0.0f;
    }

    if ((fsk->fft_est = (float *) malloc(sizeof (float) * (fsk->Ndft / 2))) == NULL) {
        free(fsk->samp_old);
        free(fsk->fftcfg);
        
        return NULL;
    }

    for (int i = 0; i < (fsk->Ndft / 2); i++) {
        fsk->fft_est[i] = 0.0f;
    }

    for (int i = 0; i < mfsk; i++) {
        fsk->f_est[i] = 0.0f;
    }

    fsk->norm_rx_timing = 0.0f;
    fsk->EbNodB = 0.0f;
    fsk->ppm = 0.0f;

    if (stats_open() == 0) {
        free(fsk->fft_est);
        free(fsk->samp_old);
        free(fsk->fftcfg);
        
        return NULL;
    }

    return fsk;
}

int fsk_get_N(struct FSK *fsk) {
    return fsk->N;
}

int fsk_get_Nmem(struct FSK *fsk) {
    return fsk->Nmem;
}

int fsk_get_Nbits(struct FSK *fsk) {
    return fsk->Nbits;
}

int fsk_get_Ts(struct FSK *fsk) {
    return fsk->Ts;
}

float fsk_get_f_est(struct FSK *fsk, int index) {
    return fsk->f_est[index];
}

/*
 * Set the minimum and maximum frequencies
 * which the frequency estimator can find tones
 */
void fsk_set_est_limits(struct FSK *fsk, int min, int max) {
    if (min < MIN_FREQUENCY)
        min = MIN_FREQUENCY;

    if (max > MAX_FREQUENCY)
        max = MAX_FREQUENCY;

    fsk->est_min = min;
    fsk->est_max = max;
}

void fsk_set_nsym(struct FSK *fsk, int nsyms) {
    fsk->N = fsk->Ts * nsyms;
    fsk->Nsym = nsyms;
    fsk->Nmem = fsk->N + (2 * fsk->Ts);
    fsk->nin = fsk->N;
    fsk->Nbits = (fsk->mode == MODE_2FSK) ? fsk->Nsym : fsk->Nsym * 2;

    int dft = 0;

    for (int i = 1; i; i <<= 1) {
        if ((fsk->N) & i) {
            dft = i;
        }
    }

    fsk->Ndft = dft;

    free(fsk->fftcfg);
    free(fsk->fft_est);

    fsk->fftcfg = fft_alloc(fsk->Ndft, 0, NULL, NULL);
    fsk->fft_est = (float *) malloc(sizeof (float) * (fsk->Ndft / 2));

    for (int i = 0; i < (fsk->Ndft / 2); i++) {
        fsk->fft_est[i] = 0.0f;
    }
}

/* Set the FSK demodulator into burst mode */

void fsk_enable_burst_mode(struct FSK *fsk, int nsyms) {
    fsk_set_nsym(fsk, nsyms);
    fsk->nin = fsk->N;
    fsk->burst_mode = true;
}

void fsk_set_estimators(struct FSK *fsk) {
    /* Clear freq estimator state */
    for (int i = 0; i < (fsk->Ndft / 2); i++) {
        fsk->fft_est[i] = 0.0f;
    }

    /* Reset timing diff correction */
    fsk->nin = fsk->N;
}

int fsk_get_nin(struct FSK *fsk) {
    return fsk->nin;
}

void fsk_destroy(struct FSK *fsk) {
    stats_close(fsk->stats);
    free(fsk->fftcfg);
    free(fsk->samp_old);
    free(fsk->fft_est);
}

static int generate_hann_table(struct FSK *fsk) {
    if ((fsk->hann_table = (float *) malloc(sizeof (float) * fsk->Ndft)) == NULL) {
        return 0;
    } else {
        complex float dphi = cmplx(TAU / (float) (fsk->Ndft - 1));
        complex float rphi = conjf(dphi) * .5f;

        for (int i = 0; i < fsk->Ndft; i++) {
            rphi = rphi * dphi;
            fsk->hann_table[i] = .5f - crealf(rphi);
        }
    }
    
    return 1;
}
