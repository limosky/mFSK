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

#include <stdint.h>
#include <complex.h>
    
void fsk_demod(struct FSK *, uint8_t [], complex float []);

#ifdef __cplusplus
}
#endif
