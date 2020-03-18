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

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>

#include "fsk.h"

struct MODULATE {
    complex float oscillator[MAX_TONES];
    complex float phase;
    float fs;
    float rs;
    float f1;
    float shift;
    int mode;
    int cycles;
};

struct MODULATE *mod_create(int, int, int, int, int);
void mod_destroy(struct MODULATE *);
void modulate(struct MODULATE *, complex float [], int);
void manchester_modulate(struct MODULATE *, complex float [], int);

#ifdef __cplusplus
}
#endif
