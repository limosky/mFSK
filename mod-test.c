#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "modulator.h"

#define MODE_2FSK 2
#define MODE_4FSK 4

#define MODE MODE_4FSK
#define FIRST_TONE_FREQ 1100
#define SHIFT_FREQ 200
#define FS 8000
#define RS 100
#define CYCLES (FS / RS)

int main() {
    FILE *fout;
    int i, j;

    if ((fout = fopen("/tmp/fsksignal.raw", "wb")) == NULL) {
        fprintf(stderr, "Error opening output waveform file\n");
        return 1;
    }

    if (mod_create(MODE, FS, RS, FIRST_TONE_FREQ, SHIFT_FREQ) == 0) {
        fprintf(stderr, "Unable to create modulator\n");
        return 2;
    }

    complex float baseband[CYCLES];
    short word[CYCLES];

    for (i = 0; i < 512; i++) {
        int bitpair = (rand() % MODE);

        modulate(baseband, bitpair);

        for (j = 0; j < CYCLES; j++) {
            // 50% modulation

            word[j] = (short) (crealf(baseband[j]) * 16384.0) & 0xFFFF;
        }

        fwrite(word, sizeof (short), CYCLES, fout);
    }

    mod_destroy();
    fclose(fout);

    return 0;
}

