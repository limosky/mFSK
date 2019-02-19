/*
 * Program to test the CRC, Scrambler, and Interleaver
 *
 * The number of bits in the frame should be less
 * than the maximum prime value. For example 22 bytes
 * is 176 bits, so you should have a prime number
 * close to that value. Larger frames means larger
 * prime table.
 * 
 * If your frames are always the same size, then you
 * only need one prime number.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define INTERLEAVE                       0
#define DEINTERLEAVE                     1

static const uint16_t primes[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
    233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347
};

static uint16_t crc16(const uint8_t *data, int length) {
    uint8_t x;
    uint16_t crc = 0xFFFF;

    while (length--) {
        x = (crc >> 8) ^ *data++;
        x ^= (x >> 4);
        crc = (crc << 8) ^ ((uint16_t) (x << 12)) ^
                ((uint16_t) (x << 5)) ^ ((uint16_t) x);
    }

    return crc;
}

static void interleave(uint8_t *inout, int nbytes, int dir) {
    uint32_t i, j, n, ibit, ibyte, ishift, jbyte, jshift;
    uint32_t b, tmp;
    uint8_t out[nbytes];

    memset(out, 0, nbytes);

    uint16_t imax = sizeof (primes) / sizeof (uint16_t);
    uint16_t nbits = (uint16_t) (nbytes * 8);

    i = 1;
    while ((primes[i] < nbits) && (i < imax))
        i++;

    b = primes[i - 1]; /* b = nearest prime to length of nbits */
    
    for (n = 0; n < nbits; n++) {
        i = n;
        j = (b * i) % nbits;

        if (dir == DEINTERLEAVE) {
            tmp = j;
            j = i;
            i = tmp;
        }

        ibyte = (i / 8);
        ishift = (i % 8);
        ibit = (inout[ibyte] >> ishift) & 0x1;

        jbyte = (j / 8);
        jshift = (j % 8);

        out[jbyte] |= (ibit << jshift);
    }

    memcpy(inout, out, nbytes);
}

static void scramble(unsigned char *inout, int nbytes) {
    int i, ibit, ibits, ibyte, ishift;

    int nbits = (nbytes * 8);

    uint16_t scrambler = 0x4a80; /* initialize scrambler seed every frame */
    uint16_t scrambler_out;

    /* in place modification of each bit */

    for (i = 0; i < nbits; i++) {
        scrambler_out = ((scrambler & 0x2) >> 1) ^ (scrambler & 0x1);

        ibyte = (i / 8);
        ishift = (i % 8);

        ibit = (inout[ibyte] >> ishift) & 0x1;
        ibits = ibit ^ scrambler_out;

        inout[ibyte] &= ~(1 << ishift);
        inout[ibyte] |= (ibits << ishift);

        /* update scrambler */

        scrambler >>= 1;
        scrambler |= (scrambler_out << 14);
    }
}

int main() {
    int nbytes = 22;
    uint8_t data[nbytes];
    uint8_t crc_data[] = { 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39 };
    int i;

    fprintf(stdout, "Testing CRC, Scrambler, and Interleaver\n\n");

    /* Test CRC */

    uint16_t val = crc16(crc_data, sizeof (crc_data));
    
    if (val == 0x29B1) {
        fprintf(stdout, "CRC Passes !\n\n");
    } else {
        fprintf(stdout, "CRC Failed = 0x%04X should be 0x29B1\n\n", val);
    }
    
    /* Test Scrambler */
    
    for (i = 0; i < nbytes; i++) {
        data[i] = rand() % 256;
    }

    for (i = 0; i < nbytes; i++) {
        fprintf(stdout, "%02X ", data[i]);
    }

    fprintf(stdout, "\n");

    scramble(data, nbytes);

    for (i = 0; i < nbytes; i++) {
        fprintf(stdout, "%02X ", data[i]);
    }

    fprintf(stdout, "\n");

    scramble(data, nbytes);

    for (i = 0; i < nbytes; i++) {
        fprintf(stdout, "%02X ", data[i]);
    }

    fprintf(stdout, "\n\nShould print the same numbers in second and third rows above\n\n");

    /* Test Interleaver */

    unsigned char inter[nbytes];
    unsigned char incopy[nbytes];

    /* copy of input for later compare   */

    for (i = 0; i < nbytes; i++)
        incopy[i] = data[i];

    interleave(data, nbytes, INTERLEAVE); /* interleave */

    memcpy(inter, data, nbytes); /* snap shot of interleaved bytes */
    scramble(data, nbytes);      /* Scramble just for giggles */
    scramble(data, nbytes);
    
    interleave(data, nbytes, DEINTERLEAVE); /* de-interleave */

    /* all ones in last column means it worked! */

    for (i = 0; i < nbytes; i++) {
        printf("%02d %02X %02X %02X %d\n",
                i+1, incopy[i], inter[i], data[i], (incopy[i] == data[i]));
    }

    printf("\nInterleaver tested OK !\n");

    return (EXIT_SUCCESS);
}

