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

#include "fft.h"
#include "fsk.h"
#include "statistics.h"
#include "demodulator.h"

/* Static prototypes/Forward declarations */

static void demodulate(struct FSK *, uint8_t [], float [], complex float []);
static void frequency_estimate(struct FSK *, complex float [], float []);

static float normf(complex float val) {
    float realf = crealf(val);
    float imagf = cimagf(val);
    return realf * realf + imagf * imagf;
}

/* Public Functions */

void fsk_demod(struct FSK *fsk, uint8_t rx_bits[], complex float fsk_in[]) {
    demodulate(fsk, rx_bits, (float *) NULL, fsk_in);
}

void fsk_demod_sd(struct FSK *fsk, float rx_sd[], complex float fsk_in[]) {
    demodulate(fsk, NULL, rx_sd, fsk_in);
}

/* Private Functions */

static void demodulate(struct FSK *fsk, uint8_t rx_bits[], float rx_sd[], complex float fsk_in[]) {
    float freqs[fsk->mode];
    int nold = (fsk->Nmem - fsk->nin);
    int using_old_samps;
    size_t i, j, m, dc_i, cbuf_i;

    complex float *f_int[fsk->mode]; /* Filtered and downsampled symbol tones */
    complex float t[fsk->mode]; /* complex number temps */
    complex float dphi[fsk->mode];
    complex float *sample_src;
    complex float f_intbuf[fsk->Ts];

    /* Estimate tone frequencies */
    frequency_estimate(fsk, fsk_in, freqs);

    /* allocate memory for the integrated samples */
    for (m = 0; m < fsk->mode; m++) {
        f_int[m] = (complex float *) malloc(sizeof (complex float) * (fsk->Nsym + 1) * fsk->P);
    }

    if (fsk->f_est[0] < 1.0f) {
        for (m = 0; m < fsk->mode; m++) {
            fsk->f_est[m] = freqs[m];
        }
    }

    /* Initialize downmixers for each symbol tone */
    for (m = 0; m < fsk->mode; m++) {
        /* Back the stored phase off to account for re-integraton of old samples */
        dphi[m] = cmplx(-TAU * (fsk->Nmem - fsk->nin - (fsk->Ts / fsk->P)) * (freqs[m] / (float) fsk->Fs));
        fsk->phase[m] = fsk->phase[m] * dphi[m];

        /* Figure out how much to nudge each sample downmixer for every sample */
        dphi[m] = cmplx(TAU * (freqs[m] / (float) fsk->Fs));
    }

    /* Integrate and downsample for symbol tones */
    for (m = 0; m < fsk->mode; m++) {
        /* Copy buffer pointers in to avoid second buffer indirection */
        float f_est_m = freqs[m];
        complex float *f_int_m = &(f_int[m][0]);
        complex float dphi_m = dphi[m];

        dc_i = 0;
        cbuf_i = 0;
        sample_src = &(fsk->samp_old[fsk->nstash - nold]);
        using_old_samps = 1;

        /* Pre-fill integration buffer */
        for (dc_i = 0; dc_i < fsk->Ts - (fsk->Ts / fsk->P); dc_i++) {
            /* Switch sample source to new samples when we run out of old ones */
            if (dc_i >= nold && using_old_samps) {
                sample_src = &fsk_in[0];
                dc_i = 0;
                using_old_samps = 0;

                /* Recalculate delta-phi after switching to new sample source */
                fsk->phase[m] = fsk->phase[m] / cabsf(fsk->phase[m]);
                dphi_m = cmplx(TAU * (f_est_m / (float) fsk->Fs));
            }

            /* Downconvert and place into integration buffer */
            f_intbuf[dc_i] = sample_src[dc_i] * conjf(fsk->phase[m]);

            /* Spin downconversion phases */
            fsk->phase[m] = fsk->phase[m] * dphi_m;
        }

        cbuf_i = dc_i;

        /* Integrate over Ts at offsets of Ts/P */
        for (i = 0; i < ((fsk->Nsym + 1) * fsk->P); i++) {
            /* Downconvert and Place Ts/P samples in the integration buffers */
            for (j = 0; j < (fsk->Ts / fsk->P); j++, dc_i++) {
                /* Switch sample source to new samples when we run out of old ones */
                if (dc_i >= nold && using_old_samps) {
                    sample_src = &fsk_in[0];
                    dc_i = 0;
                    using_old_samps = 0;

                    /* Recalculate delta-phi after switching to new sample source */
                    fsk->phase[m] = fsk->phase[m] / cabsf(fsk->phase[m]);
                    dphi_m = cmplx(TAU * (f_est_m / (float) fsk->Fs));
                }

                /* Downconvert and place into integration buffer */
                f_intbuf[cbuf_i + j] = sample_src[dc_i] * conjf(fsk->phase[m]);

                /* Spin downconversion phases */
                fsk->phase[m] = fsk->phase[m] * dphi_m;
            }

            /* Dump internal samples */
            cbuf_i = cbuf_i + (fsk->Ts / fsk->P);

            if (cbuf_i >= fsk->Ts)
                cbuf_i = 0;

            /* Integrate over the integration buffers, save samples */
            complex float it = 0.0f;

            for (j = 0; j < fsk->Ts; j++) {
                it = it + f_intbuf[j];
            }

            f_int_m[i] = it;
        }
    }

    /* Stash samples away in the old sample buffer for the next round of bit getting */
    memcpy(fsk->samp_old, &fsk_in[fsk->nin - fsk->nstash], sizeof (complex float) * fsk->nstash);

    /* Fine Timing Estimation */
    /* Apply magic nonlinearity to f1_int and f2_int, shift down to 0,
     * extract angle */

    /* Figure out how much to spin the oscillator to extract magic spectral line */
    complex float dphift = cmplx(TAU * ((float) fsk->Rs / (float) (fsk->P * fsk->Rs)));
    complex float phi_ft = cmplx(0.0f);
    complex float t_c = 0.0f;

    for (i = 0; i < ((fsk->Nsym + 1) * fsk->P); i++) {
        /* Get abs^2 of fx_int[i], and add 'em */
        float ft1 = 0.0f;

        for (m = 0; m < fsk->mode; m++) {
            ft1 += normf(f_int[m][i]);
        }

        /* Down shift and accumulate magic line */
        t_c += (ft1 * phi_ft);

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
    float rx_timing = norm_rx_timing * (float) fsk->P;

    float old_norm_rx_timing = fsk->norm_rx_timing;

    /* update the global */
    fsk->norm_rx_timing = norm_rx_timing;

    /* Estimate sample clock offset */
    float d_norm_rx_timing = norm_rx_timing - old_norm_rx_timing;

    /* Filter out big jumps in due to nin change */
    if (fabsf(d_norm_rx_timing) < .2f) {
        float appm = 1e6f * d_norm_rx_timing / (float) fsk->Nsym;

        /* Filter and update the global */
        fsk->ppm = .9f * fsk->ppm + .1f * appm;
    }

    /* Figure out how many samples are needed the next modem cycle */
    /* Unless we're in burst mode */
    if (fsk->burst_mode == false) {
        if (norm_rx_timing > 0.25f) {
            fsk->nin = fsk->N + (fsk->Ts / 2);
        } else if (norm_rx_timing < -0.25f) {
            fsk->nin = fsk->N - (fsk->Ts / 2);
        } else {
            fsk->nin = fsk->N;
        }
    }

    /* Re-sample the integrators with linear interpolation magic */
    int low_sample = (int) floorf(rx_timing);
    float fract = rx_timing - (float) low_sample;
    int high_sample = (int) ceilf(rx_timing);

    /* Vars for finding the max-of-4 for each bit */
    float tmax[fsk->mode];
    float mean_ebno = 0.0f;
    float std_ebno = 0.0f;

    /* FINALLY, THE BITS */
    /* also, resample fx_int */
    for (i = 0; i < fsk->Nsym; i++) {
        int st = (i + 1) * fsk->P;

        for (m = 0; m < fsk->mode; m++) {
            t[m] = (1.0f - fract) * f_int[m][st + low_sample];
            t[m] = t[m] + (fract * f_int[m][st + high_sample]);

            /* Figure mag^2 of each resampled fx_int */
            tmax[m] = normf(t[m]);
        }

        float max = tmax[0]; /* Maximum for figuring correct symbol */
        int sym = 0; /* Index of maximum */

        for (m = 0; m < fsk->mode; m++) {
            if (tmax[m] > max) {
                max = tmax[m];
                sym = m;
            }
        }

        /* Get the actual bit */
        if (rx_bits != NULL) {
            /* Get bits for 2FSK and 4FSK */
            /* TODO: Replace this with something more generic maybe */
            if (fsk->mode == MODE_2FSK) {
                rx_bits[i] = sym == 1; /* 2FSK. 1 bit per symbol */
            } else if (fsk->mode == MODE_4FSK) {
                rx_bits[(i * 2) + 1] = (sym & 0x1); /* 4FSK. 2 bits per symbol */
                rx_bits[(i * 2) ] = (sym & 0x2) >> 1;
            }
        }

        /* Produce soft decision symbols */
        if (rx_sd != NULL) {
            /* Convert symbols from max^2 into max */
            for (m = 0; m < fsk->mode; m++) {
                tmax[m] = sqrtf(tmax[m]);
            }

            if (fsk->mode == MODE_2FSK) {
                rx_sd[i] = tmax[0] - tmax[1];
            } else if (fsk->mode == MODE_4FSK) {
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
    mean_ebno = mean_ebno / (float) fsk->Nsym;

    /* Calculate the std. dev for EbNo dB estimate */
    std_ebno = (std_ebno / (float) fsk->Nsym) - (mean_ebno * mean_ebno);
    std_ebno = sqrt(std_ebno);

    fsk->EbNodB = -6.0f + (20.0f * log10f((1e-6f + mean_ebno) / (1e-6f + std_ebno)));

    /* Save clock offset in ppm */
    stats_set_clock_offset(fsk->stats, fsk->ppm);

    /* Calculate and save SNR from EbNo dB estimate */
    stats_set_snr_est(fsk->stats, (stats_get_snr_est(fsk->stats) * .5f) + (fsk->EbNodB * .5f));

    stats_set_rx_timing(fsk->stats, (float) rx_timing);

    /* Update frequency offset statistics */
    float rx_est_avg = (fsk->f_est[0] + fsk->f_est[1]) / 2.0f;
    float rx_center = ((fsk->f1 * 2) + fsk->shift) / 2.0f;

    stats_set_foff(fsk->stats, rx_center - rx_est_avg);

    for (m = 0; m < fsk->mode; m++) {
        free(f_int[m]);
    }
}

/*
 * Estimate frequencies of the tones within a block of samples.
 *
 * Parameters:
 * fsk_in - block of samples in this demod cycles, must be nin long
 * freqs - Array for the estimated frequencies
 */
static void frequency_estimate(struct FSK *fsk, complex float fsk_in[], float freqs[]) {
    complex float fftin[fsk->Ndft];
    complex float fftout[fsk->Ndft];
    int i, j;

    int f_min = (fsk->est_min * fsk->Ndft) / fsk->Fs;
    int f_max = (fsk->est_max * fsk->Ndft) / fsk->Fs;
    int f_zero = (fsk->est_space * fsk->Ndft) / fsk->Fs;

    /* scale averaging time constant based on number of samples */
    float tc = 0.95f * fsk->Ndft / fsk->Fs;

    for (j = 0; j < (fsk->nin / fsk->Ndft); j++) {
        /* 48000 sample rate (for example) will have a spare */
        /* 896 samples besides the 46 "Ndft" samples, so adjust */

        int samps = (fsk->nin - ((j + 1) * fsk->Ndft));
        int fft_samps = (samps >= fsk->Ndft) ? fsk->Ndft : samps;

        /* Copy FSK buffer into reals of FFT buffer and apply window */
        for (i = 0; i < fft_samps; i++) {
            fftin[i] = fsk_in[i + fsk->Ndft * j] * fsk->hann_table[i];
        }

        /* Zero out the remaining slots */
        for (; i < fsk->Ndft; i++) {
            fftin[i] = 0.0f;
        }

        /* Do the FFT */
        fft(fsk->fftcfg, fftin, fftout);

        /* Find the magnitude^2 of each freq slot and stash away in the real
         * value, so this only has to be done once. Since we're only comparing
         * these values and only need the mag of 2 points, we don't need to do
         * a sqrt to each value */
        for (i = 0; i < (fsk->Ndft / 2); i++) {
            fftout[i] = normf(fftout[i]) + cimagf(fftout[i]) * I;
        }

        /* Zero out the minimum and maximum ends */
        for (i = 0; i < f_min; i++) {
            fftout[i] = 0.0f + cimagf(fftout[i]) * I;
        }

        for (i = f_max - 1; i < (fsk->Ndft / 2); i++) {
            fftout[i] = 0.0f + cimagf(fftout[i]) * I;
        }

        /* Mix back in with the previous fft block */

        for (i = 0; i < (fsk->Ndft / 2); i++) {
            fsk->fft_est[i] = (fsk->fft_est[i] * (1.0f - tc)) + (sqrtf(crealf(fftout[i])) * tc);

            /* Copy new fft est into imag of fftout for frequency divination below */
            fftout[i] = crealf(fftout[i]) + fsk->fft_est[i] * I;
        }
    }

    int freqi[fsk->mode];
    float max;
    int imax;

    /* Find the M frequency peaks here */
    for (i = 0; i < fsk->mode; i++) {
        imax = 0;
        max = 0.0f;

        for (j = 0; j < (fsk->Ndft / 2); j++) {
            if (cimagf(fftout[j]) > max) {
                max = cimagf(fftout[j]);
                imax = j;
            }
        }

        /* Blank out FMax +/- Fspace/2 */
        f_min = imax - f_zero;
        f_min = (f_min < 0) ? 0 : f_min;
        f_max = imax + f_zero;
        f_max = (f_max > fsk->Ndft) ? fsk->Ndft : f_max;

        for (j = f_min; j < f_max; j++) {
            fftout[j] = crealf(fftout[j]);
        }

        /* Stick the freq index on the list */
        freqi[i] = imax;
    }

    /* Gnome sort the freq list */
    /* My favorite sort of sort*/
    i = 1;
    while (i < fsk->mode) {
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

    for (i = 0; i < fsk->mode; i++) {
        freqs[i] = (float) freqi[i] * ((float) fsk->Fs / (float) fsk->Ndft);
    }
}
