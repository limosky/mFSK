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
#include <complex.h>
#include <math.h>

#include "demodulator.h"
#include "modulator.h"

#define cmplx(value) (cosf(value) + sinf(value) * I)
#define cmplxconj(value) (cosf(value) + sinf(value) * -I)

/* BSS memory */

static complex float mod_oscillator[MAX_TONES];
static complex float mod_phase;
static float mod_fs;
static float mod_rs;
static float mod_f1;
static float mod_shift;
static int mod_mode;
static int mod_cycles;

int mod_create(int mode, int sample_rate, int symbol_rate, int first_freq, int shift) {
    float tone;
    int i;

    mod_mode = mode;
    mod_fs = (float) sample_rate; /* the audio sample rate Fs */
    mod_rs = (float) symbol_rate; /* the audio symbol rate Rs */
    mod_f1 = (float) first_freq; /* frequency of first tone */
    mod_shift = (float) shift; /* tone separation */
    mod_cycles = (sample_rate / symbol_rate); /* resulting number of sample cycles */

    /* Check for calling parameter errors */

    if ((sample_rate % symbol_rate) != 0) {
        return 0;
    }

    /* Initialize the FSK oscillators */

    for (i = 0, tone = mod_f1; i < mod_mode; i++, tone += mod_shift) {
        mod_oscillator[i] = cmplx(TAU * (tone / mod_fs));
    }

    /* Initialize modulator phase */

    mod_phase = cmplx(0.0f);

    return 1;
}

/* Nothing yet */
void mod_destroy() {
}

void modulate(complex float baseband[], int bits) {
    int i, tone;

    /* limit the bit parameter for safety */

    if (mod_mode == MODE_2FSK) {
        tone = bits & 0x1;
    } else if (mod_mode == MODE_4FSK) {
        tone = bits & 0x3;
    }

    /*
     * Grab one cycle of continuous phase
     * from the selected oscillator
     */
    for (i = 0; i < mod_cycles; i++) {
        mod_phase = mod_phase * mod_oscillator[tone];
        baseband[i] = mod_phase;
    }

    /* Normalize phase */
    mod_phase = mod_phase / cabsf(mod_phase);
}

/*
 * Manchester encoding is twice the bit frequency.
 * So we cut the Ts in half and alternate the signal.
 * 
 * There is no demodulator for it yet.
 */
void manchester_modulate(complex float baseband[], int bit) {
    int i;

    /* limit the bit parameter for safety */

    if ((bit & 0x1) == 0) {
        for (i = 0; i < (mod_cycles / 2); i++) {
            baseband[i] = -1.0f;
        }

        for (i = (mod_cycles / 2); i < mod_cycles; i++) {
            baseband[i] = 1.0f;
        }
    } else {
        for (i = 0; i < (mod_cycles / 2); i++) {
            baseband[i] = 1.0f;
        }

        for (i = (mod_cycles / 2); i < mod_cycles; i++) {
            baseband[i] = -1.0f;
        }
    }
}
