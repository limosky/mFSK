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

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <complex.h>
#include <math.h>

#include "demodulator.h"
#include "statistics.h"
#include "fft.h"

#define cmplx(value) (cosf(value) + sinf(value) * I)
#define cmplxconj(value) (cosf(value) + sinf(value) * -I)

#define OVERSAMPLE_RATE    8
#define MIN_FREQUENCY    300
#define MAX_FREQUENCY   2200
#define MIN_SPACING      100

/* BSS Memory */

static complex float fsk_phase[MAX_TONES];
static float fsk_f_est[MAX_TONES];
static float fsk_norm_rx_timing;
static float fsk_EbNodB;
static float fsk_ppm;
static int fsk_Ndft;
static int fsk_Fs;
static int fsk_N;
static int fsk_Rs;
static int fsk_Ts;
static int fsk_Nmem;
static int fsk_P;
static int fsk_Nsym;
static int fsk_Nbits;
static int fsk_f1;
static int fsk_shift;
static int fsk_mode;
static int fsk_est_min;
static int fsk_est_max;
static int fsk_est_space;
static int fsk_nstash;
static int fsk_nin;
static bool fsk_burst_mode;
static fft_cfg fsk_fftcfg;
static float *fsk_hann_table;
static float *fsk_fft_est;
static complex float *fsk_samp_old;

/* Static prototypes/Forward declarations */

static void demodulate(uint8_t [], float [], complex float []);
static void frequency_estimate(complex float [], float []);
static float *fsk_generate_hann_table(int);

/* Public Functions */

static float normf(complex float val) {
    float realf = crealf(val);
    float imagf = cimagf(val);
    return realf * realf + imagf * imagf;
}

void fsk_demod(uint8_t rx_bits[], complex float fsk_in[]) {
    demodulate(rx_bits, (float *) NULL, fsk_in);
}

void fsk_demod_sd(float rx_sd[], complex float fsk_in[]) {
    demodulate(NULL, rx_sd, fsk_in);
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: fsk_create
  AUTHOR......: Brady O'Brien
  DATE CREATED: 7 January 2016

  Create and initialize an instance of the FSK modem.

\*---------------------------------------------------------------------------*/

int fsk_create(int sample_rate, int symbol_rate, int mfsk, int tone_f1) {
    int i;

    /* Set constant config parameters */
    fsk_Fs = sample_rate;
    fsk_Rs = symbol_rate;
    fsk_Ts = sample_rate / symbol_rate;
    fsk_N = sample_rate;
    fsk_burst_mode = false;
    fsk_P = OVERSAMPLE_RATE;
    fsk_Nsym = fsk_N / fsk_Ts;
    fsk_Ndft = FFT_SIZE;
    fsk_Nmem = fsk_N + (fsk_Ts * 2);
    fsk_f1 = tone_f1;
    fsk_shift = (symbol_rate * 2);
    fsk_nin = fsk_N;
    fsk_mode = mfsk;
    fsk_Nbits = (mfsk == MODE_2FSK) ? fsk_Nsym : fsk_Nsym * 2;

    fsk_est_min = MIN_FREQUENCY;
    fsk_est_max = MAX_FREQUENCY;
    fsk_est_space = MIN_SPACING;

    /* Set up rx state */
    for (i = 0; i < mfsk; i++) {
        fsk_phase[i] = cmplx(0.0f);
    }

    fsk_nstash = (fsk_Ts * 4);

    if (fsk_generate_hann_table(fsk_Ndft) == NULL) {
        return 0;
    }

    if ((fsk_samp_old = (complex float *) malloc(sizeof (complex float) * fsk_nstash)) == NULL) {
        return 0;
    }

    for (i = 0; i < fsk_nstash; i++) {
        fsk_samp_old[i] = 0.0f;
    }

    if ((fsk_fftcfg = fft_alloc(fsk_Ndft, 0, NULL, NULL)) == NULL) {
        free(fsk_samp_old);
        
        return 0;
    }

    if ((fsk_fft_est = (float *) malloc(sizeof (float) * (fsk_Ndft / 2))) == NULL) {
        free(fsk_samp_old);
        free(fsk_fftcfg);
        
        return 0;
    }

    for (i = 0; i < (fsk_Ndft / 2); i++) {
        fsk_fft_est[i] = 0.0f;
    }

    for (i = 0; i < mfsk; i++) {
        fsk_f_est[i] = 0.0f;
    }

    fsk_norm_rx_timing = 0.0f;
    fsk_EbNodB = 0.0f;
    fsk_ppm = 0.0f;

    if (stats_open() == 0) {
        free(fsk_fft_est);
        free(fsk_samp_old);
        free(fsk_fftcfg);
        
        return 0;
    }

    return 1;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: fsk_create_hbr
  AUTHOR......: Brady O'Brien
  DATE CREATED: 11 February 2016

  Create and initialize a high bit-rate (HBR) instance of the FSK demod.
  Returns a 1 on success, or 0 on failure

\*---------------------------------------------------------------------------*/

int fsk_create_hbr(int sample_rate, int symbol_rate, int oversample_rate, int mfsk, int f1, int shift) {
    int i;

    /* Set constant config parameters */
    fsk_Fs = sample_rate;
    fsk_Rs = symbol_rate;
    fsk_Ts = sample_rate / symbol_rate; /* cycles per symbol */
    fsk_burst_mode = false;
    fsk_Nsym = 48; /* Number of symbols in a processing frame */
    fsk_N = fsk_Ts * fsk_Nsym;
    fsk_P = oversample_rate;
    fsk_Nmem = fsk_N + (fsk_Ts * 2);
    fsk_f1 = f1;
    fsk_shift = shift;
    fsk_nin = fsk_N;
    fsk_mode = mfsk;
    fsk_Nbits = (mfsk == MODE_2FSK) ? fsk_Nsym : fsk_Nsym * 2;

    int dft = 0;

    for (i = 1; i; i <<= 1)
        if (fsk_N & i)
            dft = i;

    fsk_Ndft = dft;

    fsk_est_min = (fsk_Rs / 2);
    fsk_est_max = (fsk_Fs / 2) - fsk_Rs;
    fsk_est_space = fsk_Rs - (fsk_Rs / 5);

    /* Set up rx state */

    for (i = 0; i < mfsk; i++) {
        fsk_phase[i] = cmplx(0.0f);
    }

    fsk_nstash = (fsk_Ts * 4);

    if (fsk_generate_hann_table(fsk_Ndft) == NULL) {
        return 0;
    }

    if ((fsk_fftcfg = fft_alloc(fsk_Ndft, 0, NULL, NULL)) == NULL) {
        return 0;
    }

    if ((fsk_samp_old = (complex float *) malloc(sizeof (complex float) * fsk_nstash)) == NULL) {
        return 0;
    }

    for (i = 0; i < fsk_nstash; i++) {
        fsk_samp_old[i] = 0.0f;
    }

    if ((fsk_fft_est = (float *) malloc(sizeof (float) * (fsk_Ndft / 2))) == NULL) {
        free(fsk_samp_old);
        free(fsk_fftcfg);
        
        return 0;
    }

    for (i = 0; i < (fsk_Ndft / 2); i++) {
        fsk_fft_est[i] = 0.0f;
    }

    for (i = 0; i < mfsk; i++) {
        fsk_f_est[i] = 0.0f;
    }

    fsk_norm_rx_timing = 0.0f;
    fsk_EbNodB = 0.0f;
    fsk_ppm = 0.0f;

    if (stats_open() == 0) {
        free(fsk_fft_est);
        free(fsk_samp_old);
        free(fsk_fftcfg);
        
        return 0;
    }

    return 1;
}

int fsk_get_N() {
    return fsk_N;
}

int fsk_get_Nmem() {
    return fsk_Nmem;
}

int fsk_get_Nbits() {
    return fsk_Nbits;
}

int fsk_get_Ts() {
    return fsk_Ts;
}

float fsk_get_f_est(int index) {
    return fsk_f_est[index];
}

/*
 * Set the minimum and maximum frequencies
 * which the frequency estimator can find tones
 */
void fsk_set_est_limits(int min, int max) {
    if (min < MIN_FREQUENCY)
        min = MIN_FREQUENCY;

    if (max > MAX_FREQUENCY)
        max = MAX_FREQUENCY;

    fsk_est_min = min;
    fsk_est_max = max;
}

void fsk_set_nsym(int nsyms) {
    int i;

    fsk_N = fsk_Ts * nsyms;
    fsk_Nsym = nsyms;
    fsk_Nmem = fsk_N + (2 * fsk_Ts);
    fsk_nin = fsk_N;
    fsk_Nbits = (fsk_mode == MODE_2FSK) ? fsk_Nsym : fsk_Nsym * 2;

    int dft = 0;

    for (i = 1; i; i <<= 1) {
        if ((fsk_N) & i) {
            dft = i;
        }
    }

    fsk_Ndft = dft;

    free(fsk_fftcfg);
    free(fsk_fft_est);

    fsk_fftcfg = fft_alloc(fsk_Ndft, 0, NULL, NULL);
    fsk_fft_est = (float *) malloc(sizeof (float) * (fsk_Ndft / 2));

    for (i = 0; i < (fsk_Ndft / 2); i++) {
        fsk_fft_est[i] = 0.0f;
    }
}

/* Set the FSK demodulator into burst mode */

void fsk_enable_burst_mode(int nsyms) {
    fsk_set_nsym(nsyms);
    fsk_nin = fsk_N;
    fsk_burst_mode = true;
}

void fsk_set_estimators() {
    int i;

    /* Clear freq estimator state */
    for (i = 0; i < (fsk_Ndft / 2); i++) {
        fsk_fft_est[i] = 0.0f;
    }

    /* Reset timing diff correction */
    fsk_nin = fsk_N;
}

int fsk_get_nin() {
    return fsk_nin;
}

void fsk_destroy() {
    stats_close();
    free(fsk_fftcfg);
    free(fsk_samp_old);
    free(fsk_fft_est);
}

/* Local Functions */

static float *fsk_generate_hann_table(int ndft) {
    if ((fsk_hann_table = (float *) malloc(sizeof (float) * ndft)) == NULL) {
        return NULL;
    } else {
        complex float dphi = cmplx(TAU / (float) (ndft - 1));
        complex float rphi = conjf(dphi) * .5f;
        int i;

        for (i = 0; i < ndft; i++) {
            rphi = rphi * dphi;
            fsk_hann_table[i] = .5f - crealf(rphi);
        }
    }

    return fsk_hann_table;
}

/*
 * Estimate frequencies of the tones within a block of samples.
 *
 * Parameters:
 * fsk_in - block of samples in this demod cycles, must be nin long
 * freqs - Array for the estimated frequencies
 */
static void frequency_estimate(complex float fsk_in[], float freqs[]) {
    complex float fftin[fsk_Ndft];
    complex float fftout[fsk_Ndft];
    int i, j;

    int f_min = (fsk_est_min * fsk_Ndft) / fsk_Fs;
    int f_max = (fsk_est_max * fsk_Ndft) / fsk_Fs;
    int f_zero = (fsk_est_space * fsk_Ndft) / fsk_Fs;

    /* scale averaging time constant based on number of samples */
    float tc = 0.95f * fsk_Ndft / fsk_Fs;

    for (j = 0; j < (fsk_nin / fsk_Ndft); j++) {
        /* 48000 sample rate (for example) will have a spare */
        /* 896 samples besides the 46 "Ndft" samples, so adjust */

        int samps = (fsk_nin - ((j + 1) * fsk_Ndft));
        int fft_samps = (samps >= fsk_Ndft) ? fsk_Ndft : samps;

        /* Copy FSK buffer into reals of FFT buffer and apply window */
        for (i = 0; i < fft_samps; i++) {
            fftin[i] = fsk_in[i + fsk_Ndft * j] * fsk_hann_table[i];
        }

        /* Zero out the remaining slots */
        for (; i < fsk_Ndft; i++) {
            fftin[i] = 0.0f;
        }

        /* Do the FFT */
        fft(fsk_fftcfg, fftin, fftout);

        /* Find the magnitude^2 of each freq slot and stash away in the real
         * value, so this only has to be done once. Since we're only comparing
         * these values and only need the mag of 2 points, we don't need to do
         * a sqrt to each value */
        for (i = 0; i < (fsk_Ndft / 2); i++) {
            fftout[i] = normf(fftout[i]) + cimagf(fftout[i]) * I;
        }

        /* Zero out the minimum and maximum ends */
        for (i = 0; i < f_min; i++) {
            fftout[i] = 0.0f + cimagf(fftout[i]) * I;
        }

        for (i = f_max - 1; i < (fsk_Ndft / 2); i++) {
            fftout[i] = 0.0f + cimagf(fftout[i]) * I;
        }

        /* Mix back in with the previous fft block */

        for (i = 0; i < (fsk_Ndft / 2); i++) {
            fsk_fft_est[i] = (fsk_fft_est[i] * (1.0f - tc)) + (sqrtf(crealf(fftout[i])) * tc);

            /* Copy new fft est into imag of fftout for frequency divination below */
            fftout[i] = crealf(fftout[i]) + fsk_fft_est[i] * I;
        }
    }

    int freqi[fsk_mode];
    float max;
    int imax;

    /* Find the M frequency peaks here */
    for (i = 0; i < fsk_mode; i++) {
        imax = 0;
        max = 0.0f;

        for (j = 0; j < (fsk_Ndft / 2); j++) {
            if (cimagf(fftout[j]) > max) {
                max = cimagf(fftout[j]);
                imax = j;
            }
        }

        /* Blank out FMax +/- Fspace/2 */
        f_min = imax - f_zero;
        f_min = (f_min < 0) ? 0 : f_min;
        f_max = imax + f_zero;
        f_max = (f_max > fsk_Ndft) ? fsk_Ndft : f_max;

        for (j = f_min; j < f_max; j++) {
            fftout[j] = crealf(fftout[j]);
        }

        /* Stick the freq index on the list */
        freqi[i] = imax;
    }

    /* Gnome sort the freq list */
    /* My favorite sort of sort*/
    i = 1;
    while (i < fsk_mode) {
        if (freqi[i] >= freqi[i - 1]) {
            i++;
        } else {
            j = freqi[i];
            freqi[i] = freqi[i - 1];
            freqi[i - 1] = j;

            if (i > 1) {
                i--;
            }
        }
    }

    /* Convert freqs from indices to frequencies */

    for (i = 0; i < fsk_mode; i++) {
        freqs[i] = (float) freqi[i] * ((float) fsk_Fs / (float) fsk_Ndft);
    }
}

static void demodulate(uint8_t rx_bits[], float rx_sd[], complex float fsk_in[]) {
    float freqs[fsk_mode];
    int nold = (fsk_Nmem - fsk_nin);
    int using_old_samps;
    size_t i, j, m, dc_i, cbuf_i;

    complex float *f_int[fsk_mode]; /* Filtered and downsampled symbol tones */
    complex float t[fsk_mode]; /* complex number temps */
    complex float dphi[fsk_mode];
    complex float *sample_src;
    complex float f_intbuf[fsk_Ts];

    /* Estimate tone frequencies */
    frequency_estimate(fsk_in, freqs);

    /* allocate memory for the integrated samples */
    for (m = 0; m < fsk_mode; m++) {
        f_int[m] = (complex float *) malloc(sizeof (complex float) * (fsk_Nsym + 1) * fsk_P);
    }

    if (fsk_f_est[0] < 1.0f) {
        for (m = 0; m < fsk_mode; m++) {
            fsk_f_est[m] = freqs[m];
        }
    }

    /* Initialize downmixers for each symbol tone */
    for (m = 0; m < fsk_mode; m++) {
        /* Back the stored phase off to account for re-integraton of old samples */
        dphi[m] = cexpf(I * (-TAU * (fsk_Nmem - fsk_nin - (fsk_Ts / fsk_P)) * (freqs[m] / (float) fsk_Fs)));
        fsk_phase[m] = fsk_phase[m] * dphi[m];

        /* Figure out how much to nudge each sample downmixer for every sample */
        dphi[m] = cexpf(I * (TAU * (freqs[m] / (float) fsk_Fs)));
    }

    /* Integrate and downsample for symbol tones */
    for (m = 0; m < fsk_mode; m++) {
        /* Copy buffer pointers in to avoid second buffer indirection */
        float f_est_m = freqs[m];
        complex float *f_int_m = &(f_int[m][0]);
        complex float dphi_m = dphi[m];

        dc_i = 0;
        cbuf_i = 0;
        sample_src = &(fsk_samp_old[fsk_nstash - nold]);
        using_old_samps = 1;

        /* Pre-fill integration buffer */
        for (dc_i = 0; dc_i < fsk_Ts - (fsk_Ts / fsk_P); dc_i++) {
            /* Switch sample source to new samples when we run out of old ones */
            if (dc_i >= nold && using_old_samps) {
                sample_src = &fsk_in[0];
                dc_i = 0;
                using_old_samps = 0;

                /* Recalculate delta-phi after switching to new sample source */
                fsk_phase[m] = fsk_phase[m] / cabsf(fsk_phase[m]);
                dphi_m = cmplx(TAU * (f_est_m / (float) fsk_Fs));
            }

            /* Downconvert and place into integration buffer */
            f_intbuf[dc_i] = sample_src[dc_i] * conjf(fsk_phase[m]);

            /* Spin downconversion phases */
            fsk_phase[m] = fsk_phase[m] * dphi_m;
        }

        cbuf_i = dc_i;

        /* Integrate over Ts at offsets of Ts/P */
        for (i = 0; i < ((fsk_Nsym + 1) * fsk_P); i++) {
            /* Downconvert and Place Ts/P samples in the integration buffers */
            for (j = 0; j < (fsk_Ts / fsk_P); j++, dc_i++) {
                /* Switch sample source to new samples when we run out of old ones */
                if (dc_i >= nold && using_old_samps) {
                    sample_src = &fsk_in[0];
                    dc_i = 0;
                    using_old_samps = 0;

                    /* Recalculate delta-phi after switching to new sample source */
                    fsk_phase[m] = fsk_phase[m] / cabsf(fsk_phase[m]);
                    dphi_m = cmplx(TAU * (f_est_m / (float) fsk_Fs));
                }

                /* Downconvert and place into integration buffer */
                f_intbuf[cbuf_i + j] = sample_src[dc_i] * conjf(fsk_phase[m]);

                /* Spin downconversion phases */
                fsk_phase[m] = fsk_phase[m] * dphi_m;
            }

            /* Dump internal samples */
            cbuf_i = cbuf_i + (fsk_Ts / fsk_P);

            if (cbuf_i >= fsk_Ts)
                cbuf_i = 0;

            /* Integrate over the integration buffers, save samples */
            complex float it = 0.0f;

            for (j = 0; j < fsk_Ts; j++) {
                it = it + f_intbuf[j];
            }

            f_int_m[i] = it;
        }
    }

    /* Stash samples away in the old sample buffer for the next round of bit getting */
    memcpy(fsk_samp_old, &fsk_in[fsk_nin - fsk_nstash], sizeof (complex float) * fsk_nstash);

    /* Fine Timing Estimation */
    /* Apply magic nonlinearity to f1_int and f2_int, shift down to 0,
     * extract angle */

    /* Figure out how much to spin the oscillator to extract magic spectral line */
    complex float dphift = cmplx(TAU * ((float) fsk_Rs / (float) (fsk_P * fsk_Rs)));
    complex float phi_ft = cmplx(0.0f);
    complex float t_c = cmplx(0.0f);

    for (i = 0; i < ((fsk_Nsym + 1) * fsk_P); i++) {
        /* Get abs^2 of fx_int[i], and add 'em */
        float ft1 = 0.0f;

        for (m = 0; m < fsk_mode; m++) {
            ft1 += normf(f_int[m][i]);
        }

        /* Down shift and accumulate magic line */
        t_c = t_c + (ft1 * phi_ft);

        /* Spin the oscillator for the magic line shift */
        phi_ft = phi_ft * dphift;
    }

    /* Check for NaNs in the fine timing estimate, return if found */
    /* otherwise segfaults happen */
    if (isnan(crealf(t_c)) || isnan(cimagf(t_c))) {
        return;
    }

    /* Get the magic angle */
    float norm_rx_timing = cargf(t_c) / TAU;
    float rx_timing = norm_rx_timing * (float) fsk_P;

    float old_norm_rx_timing = fsk_norm_rx_timing;

    /* update the global */
    fsk_norm_rx_timing = norm_rx_timing;

    /* Estimate sample clock offset */
    float d_norm_rx_timing = norm_rx_timing - old_norm_rx_timing;

    /* Filter out big jumps in due to nin change */
    if (fabsf(d_norm_rx_timing) < .2f) {
        float appm = 1e6f * d_norm_rx_timing / (float) fsk_Nsym;

        /* Filter and update the global */
        fsk_ppm = .9f * fsk_ppm + .1f * appm;
    }

    /* Figure out how many samples are needed the next modem cycle */
    /* Unless we're in burst mode */
    if (fsk_burst_mode == false) {
        if (norm_rx_timing > 0.25f) {
            fsk_nin = fsk_N + (fsk_Ts / 2);
        } else if (norm_rx_timing < -0.25f) {
            fsk_nin = fsk_N - (fsk_Ts / 2);
        } else {
            fsk_nin = fsk_N;
        }
    }

    /* Re-sample the integrators with linear interpolation magic */
    int low_sample = (int) floorf(rx_timing);
    float fract = rx_timing - (float) low_sample;
    int high_sample = (int) ceilf(rx_timing);

    /* Vars for finding the max-of-4 for each bit */
    float tmax[fsk_mode];
    float mean_ebno = 0.0f;
    float std_ebno = 0.0f;

    /* FINALLY, THE BITS */
    /* also, resample fx_int */
    for (i = 0; i < fsk_Nsym; i++) {
        int st = (i + 1) * fsk_P;

        for (m = 0; m < fsk_mode; m++) {
            t[m] = (1.0f - fract) * f_int[m][st + low_sample];
            t[m] = t[m] + (fract * f_int[m][st + high_sample]);

            /* Figure mag^2 of each resampled fx_int */
            tmax[m] = normf(t[m]);
        }

        float max = tmax[0]; /* Maximum for figuring correct symbol */
        int sym = 0; /* Index of maximum */

        for (m = 0; m < fsk_mode; m++) {
            if (tmax[m] > max) {
                max = tmax[m];
                sym = m;
            }
        }

        /* Get the actual bit */
        if (rx_bits != NULL) {
            /* Get bits for 2FSK and 4FSK */
            /* TODO: Replace this with something more generic maybe */
            if (fsk_mode == MODE_2FSK) {
                rx_bits[i] = sym == 1; /* 2FSK. 1 bit per symbol */
            } else if (fsk_mode == MODE_4FSK) {
                rx_bits[(i * 2) + 1] = (sym & 0x1); /* 4FSK. 2 bits per symbol */
                rx_bits[(i * 2) ] = (sym & 0x2) >> 1;
            }
        }

        /* Produce soft decision symbols */
        if (rx_sd != NULL) {
            /* Convert symbols from max^2 into max */
            for (m = 0; m < fsk_mode; m++) {
                tmax[m] = sqrtf(tmax[m]);
            }

            if (fsk_mode == MODE_2FSK) {
                rx_sd[i] = tmax[0] - tmax[1];
            } else if (fsk_mode == MODE_4FSK) {
                rx_sd[(i * 2) + 1] = -tmax[0]; /* Bits=00 */
                rx_sd[(i * 2) ] = -tmax[0];

                rx_sd[(i * 2) + 1] += tmax[1]; /* Bits=01 */
                rx_sd[(i * 2) ] += -tmax[1];
                rx_sd[(i * 2) + 1] += -tmax[2]; /* Bits=10 */
                rx_sd[(i * 2) ] += tmax[2];
                rx_sd[(i * 2) + 1] += tmax[3]; /* Bits=11 */
                rx_sd[(i * 2) ] += tmax[3];
            }
        }

        /* Accumulate resampled int magnitude for EbNodB estimation */
        /* Standard deviation is calculated by algorithm devised by crafty soviets */
        /* Accumulate the square of the sampled value */
        std_ebno += max;

        /* Figure the abs value of the max tone */
        mean_ebno += sqrtf(max);

        /* Soft output goes here */
    }

    /* Calculate mean for EbNo dB estimation */
    mean_ebno = mean_ebno / (float) fsk_Nsym;

    /* Calculate the std. dev for EbNo dB estimate */
    std_ebno = (std_ebno / (float) fsk_Nsym) - (mean_ebno * mean_ebno);
    std_ebno = sqrt(std_ebno);

    fsk_EbNodB = -6.0f + (20.0f * log10f((1e-6f + mean_ebno) / (1e-6f + std_ebno)));

    /* Save clock offset in ppm */
    stats_set_clock_offset(fsk_ppm);

    /* Calculate and save SNR from EbNo dB estimate */
    stats_set_snr_est((stats_get_snr_est() * .5f) + (fsk_EbNodB * .5f));

    stats_set_rx_timing((float) rx_timing);

    /* Update frequency offset statistics */
    float rx_est_avg = (fsk_f_est[0] + fsk_f_est[1]) / 2.0f;
    float rx_center = ((fsk_f1 * 2) + fsk_shift) / 2.0f;

    stats_set_foff(rx_center - rx_est_avg);

    for (m = 0; m < fsk_mode; m++) {
        free(f_int[m]);
    }
}
