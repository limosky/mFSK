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

#include "modulator.h"

struct MODULATE *mod_create(int mode, int sample_rate, int symbol_rate, int first_freq, int shift) {
    struct MODULATE *mod;
    float tone;
    int i;

    if ((mod = (struct MODULATE *) malloc(sizeof (struct MODULATE))) == NULL) {
        return NULL;
    }

    mod->mode = mode;
    mod->fs = (float) sample_rate; /* the audio sample rate Fs */
    mod->rs = (float) symbol_rate; /* the audio symbol rate Rs */
    mod->f1 = (float) first_freq; /* frequency of first tone */
    mod->shift = (float) shift; /* tone separation */
    mod->cycles = (sample_rate / symbol_rate); /* resulting number of sample cycles */

    /* Check for calling parameter errors */

    if ((sample_rate % symbol_rate) != 0) {
        return 0;
    }

    /* Initialize the FSK oscillators */

    for (i = 0, tone = mod->f1; i < mod->mode; i++, tone += mod->shift) {
        mod->oscillator[i] = cmplx(TAU * (tone / mod->fs));
    }

    /* Initialize modulator phase */

    mod->phase = cmplx(0.0f);

    return mod;
}

/* Nothing yet */
void mod_destroy(struct MODULATE *mod) {
    free(mod);
}

void modulate(struct MODULATE *mod, complex float baseband[], int bits) {
    int i, tone;

    /* limit the bit parameter for safety */

    if (mod->mode == MODE_2FSK) {
        tone = bits & 0x1;
    } else if (mod->mode == MODE_4FSK) {
        tone = bits & 0x3;
    }

    /*
     * Grab one cycle of continuous phase
     * from the selected oscillator
     */
    for (i = 0; i < mod->cycles; i++) {
        mod->phase = mod->phase * mod->oscillator[tone];
        baseband[i] = mod->phase;
    }

    /* Normalize phase */
    mod->phase = mod->phase / cabsf(mod->phase);
}

/*
 * Manchester encoding is twice the bit frequency.
 * So we cut the Ts in half and alternate the signal.
 * 
 * There is no demodulator for it yet.
 */
void manchester_modulate(struct MODULATE *mod, complex float baseband[], int bit) {
    int i;

    /* limit the bit parameter for safety */

    if ((bit & 0x1) == 0) {
        for (i = 0; i < (mod->cycles / 2); i++) {
            baseband[i] = -1.0f;
        }

        for (i = (mod->cycles / 2); i < mod->cycles; i++) {
            baseband[i] = 1.0f;
        }
    } else {
        for (i = 0; i < (mod->cycles / 2); i++) {
            baseband[i] = 1.0f;
        }

        for (i = (mod->cycles / 2); i < mod->cycles; i++) {
            baseband[i] = -1.0f;
        }
    }
}
