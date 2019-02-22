/*
 * Copyright (C) March 2013 David Rowe, Tomas Härdin
 *
 * All rights reserved.
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */

#include <stdint.h>

#include "golayenctable.h"
#include "golaydectable.h"

#define POLY    0xC75   //AE3 reversed
#define popcount __builtin_popcount

static int syndrome(int c) {
    int x;

    for (x = 11; x >= 0; x--) {
        if (c & ((1 << 11) << x)) {
            c ^= (POLY << x);
        }
    }

    return c;
}

static int golay23_encode_no_tables(int c) {
    c <<= 11;
    
    return syndrome(c) | c;
}

static int unrotate(unsigned int c, int x) {
    return ((c << x) & 0x7FFFFF) | (c >> (23 - x));
}

static int golay23_decode_no_tables(int c) {
    int c2, x, t, s;
    
    c = unrotate(c, 12);

    for (x = 0; x < 23; x++) {
        s = syndrome(c);

        if (popcount(s) <= 3) {
            return unrotate(c ^ s, x) & 0xFFF;
        }

        for (t = 0; t < 23; t++) {
            c2 = c ^ (1 << t);
            s = syndrome(c2);

            if (popcount(s) <= 2) {
                return unrotate(c2 ^ s, x) & 0xFFF;
            }
        }

        //rotate
        c = (c >> 1) | ((c & 1) << 22);
    }

    return c & 0xFFF;
}

void golay23_init(void) {
}

int  golay23_encode(int c) {
    return golay23_encode_no_tables(c);
}

int  golay23_decode(int c) {
    return unrotate(golay23_decode_no_tables(c), 11);
}

int  golay23_count_errors(int recd_codeword, int corrected_codeword) {
    return popcount(recd_codeword ^ corrected_codeword);
}
