/*
 * Copyright (C) March 2013 David Rowe, Tomas Härdin
 *
 * All rights reserved.
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void golay23_init(void);
int golay23_encode(int);
int golay23_decode(int);
int golay23_count_errors(int, int);

#ifdef __cplusplus
}
#endif
