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
extern "C" {
#endif

#include <complex.h>
    
int mod_create(int, int, int, int, int);
void mod_destroy(void);
void modulate(complex float [], int);
void manchester_modulate(complex float [], int);

#ifdef __cplusplus
}
#endif
