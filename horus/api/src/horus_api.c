/*
 * Copyright (C) March 2018 David Rowe
 *
 * All rights reserved.
 * 
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 * 
 * Modified API to combine l2 functions
 */

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "fsk.h"
#include "modulator.h"
#include "demodulator.h"
#include "golay23.h"
#include "statistics.h"
#include "horus_api.h"

/* Unique word for Horus RTTY 7 bit '$' character, 3 sync bits,
   repeated 5 times */

static const int8_t uw_horus_rtty[] = {
    0, 0, 1, 0, 0, 1, 0, 1, 1, 0,
    0, 0, 1, 0, 0, 1, 0, 1, 1, 0,
    0, 0, 1, 0, 0, 1, 0, 1, 1, 0,
    0, 0, 1, 0, 0, 1, 0, 1, 1, 0,
    0, 0, 1, 0, 0, 1, 0, 1, 1, 0
};

/* Unique word for Horus Binary 0x24 = $ */

static const int8_t uw_horus_binary[] = {
    0, 0, 1, 0, 0, 1, 0, 0,
    0, 0, 1, 0, 0, 1, 0, 0
};

static const char l2_uw[] = {
    '$', '$'
};

static const uint16_t primes[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
    233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347
};

/* BSS Memory */

struct FSK *api_fsk;
struct MODULATE *api_mod;
struct STATS *api_stats;

static int api_verbose;
static int api_mFSK;                    /* number of FSK tones */
static int api_mode;
static int api_Fs;                      /* sample rate in Hz */
static int api_Rs;                      /* symbol rate in Hz */
static int api_uw_thresh;               /* threshold for UW detection */
static int api_uw_len;                  /* length of unique word */
static int api_max_packet_len;          /* max length of a telemetry packet */
static int api_rx_bits_len;             /* length of rx_bits buffer */
static int api_crc_ok;                  /* most recent packet checksum results */
static int8_t api_uw[MAX_UW_LENGTH];    /* unique word bits mapped to +/-1 */
static uint8_t *api_rx_bits;            /* buffer of received bits */

/* Static prototypes/Forward declarations */

static uint16_t crc16(const uint8_t *, int);
static int hex2int(char);
static int horus_find_uw(int);
static int horus_get_num_tx_data_bytes(int);
static int horus_encode_packet(uint8_t *, uint8_t *, int);
static void horus_decode_packet(uint8_t *, uint8_t *, int);
static int horus_extract_rtty(char [], int);
static int horus_extract_binary(char [], int);

/* Public Functions */

/**
 * Initialize API Variables and open and FSK session
 * 
 * hmode - An integer to specify either RTTY or Binary mode.
 * 
 * Returns 1 on success and 0 on error.
 */
int horus_open(int hmode) {
    int i;

    if ((hmode != HORUS_MODE_RTTY) && (hmode != HORUS_MODE_BINARY)) {
        return 0;
    }

    api_Fs = 48000;
    api_Rs = 100; /* 100 baud */
    api_verbose = 0;
    api_mode = hmode;

    if (hmode == HORUS_MODE_RTTY) {
        api_mFSK = MODE_2FSK;
        api_max_packet_len = 1000;  /* bits */

        /* map UW to make it easier to search for */

        for (i = 0; i < sizeof (uw_horus_rtty); i++) {
            api_uw[i] = (uw_horus_rtty[i] * 2) - 1;
        }

        api_uw_len = sizeof (uw_horus_rtty);
        api_uw_thresh = (sizeof (uw_horus_rtty) - 2); /* allow a few bit errors in UW detection */
        api_rx_bits_len = api_max_packet_len;
    } else if (hmode == HORUS_MODE_BINARY) {
        api_mFSK = MODE_4FSK;
        api_max_packet_len = HORUS_BINARY_NUM_BITS;

        for (i = 0; i < sizeof (uw_horus_binary); i++) {
            api_uw[i] = (uw_horus_binary[i] * 2) - 1;
        }

        api_uw_len = sizeof (uw_horus_binary);
        api_uw_thresh = (sizeof (uw_horus_binary) - 2); /* allow a few bit errors in UW detection */

        golay23_init();

        api_rx_bits_len = api_max_packet_len;
    }

    /* 1100 Hz FSK F1 frequency */

    if ((api_fsk = fsk_create(api_Fs, api_Rs, api_mFSK, 1100)) == NULL) {
        return 0;
    }

    /* shift frequency 2 * Rs = 200 Hz */

    if ((api_mod = mod_create(api_mFSK, api_Fs, api_Rs, 1100, api_Rs * 2)) == NULL) {
        return 0;
    }
    
    if ((api_stats = stats_open()) == NULL) {
        return 0;
    }
    
    /* allocate enough room for two packets so we know there will be
       one complete packet if we find a UW at start */

    api_rx_bits_len += fsk_get_Nbits(api_fsk);

    if ((api_rx_bits = (uint8_t *) malloc(api_rx_bits_len)) == NULL) {
        return 0;
    }

    for (i = 0; i < api_rx_bits_len; i++) {
        api_rx_bits[i] = 0;
    }

    api_crc_ok = 0;

    return 1;
}

/**
 * Close the FSK session and free memory.
 */
void horus_close() {
    free(api_rx_bits);
    stats_close(api_stats);
    mod_destroy(api_mod);
    fsk_destroy(api_fsk);

}

/**
 * data_out  - ASCII or HEX data output signal
 * demod_in  - 16-bit short for PCM RAW input signal
 *
 * Returns 1 for RTTY checksum good, or 1 for binary data detected,
 * or else returns 0.
 */

int horus_rx(char data_out[], short demod_in[]) {
    int i, j, uw_loc;
    int nin = fsk_get_nin(api_fsk);
    int packet_detected = 0;
    int Nbits = fsk_get_Nbits(api_fsk);

    if (api_verbose) {
        fprintf(stderr, "max_packet_len: %d rx_bits_len: %d Nbits: %d\n",
                api_max_packet_len, api_rx_bits_len, Nbits);
    }

    /* shift buffer of bits to make room for new bits */

    for (i = 0, j = Nbits; j < api_rx_bits_len; i++, j++) {
        api_rx_bits[i] = api_rx_bits[j];
    }

    /* demodulate latest bits */

    complex float demod_in_comp[nin];

    for (i = 0; i < nin; i++) {
        demod_in_comp[i] = (float) demod_in[i] + 0.0f * I;
    }

    fsk_demod(api_fsk, &api_rx_bits[api_rx_bits_len - Nbits], demod_in_comp);

    /* UW search to see if we can find the start of a packet in the buffer */

    if ((uw_loc = horus_find_uw(Nbits)) != -1) {
        if (api_verbose) {
            fprintf(stderr, "uw_loc: %d mode: %d\n", uw_loc, api_mode);
        }

        /* OK we have found a unique word, and therefore the start of
           a packet, so lets try to extract valid packets */

        if (api_mode == HORUS_MODE_RTTY) {
            packet_detected = horus_extract_rtty(data_out, uw_loc);
        } else if (api_mode == HORUS_MODE_BINARY) {
            packet_detected = horus_extract_binary(data_out, uw_loc);
        }
    }

    return packet_detected;
}

/**
 * This function is for Binary packets
 * 
 * demod_out - 16 bit short for PCM RAW output signal
 * percent   - integer between 10% and 88% for signal amplitude
 * payload   - a pointer to unsigned character data to be transmitted
 *
 * Returns an integer of PCM demodulator words in demod_out
 */

int horus_binary_tx(short demod_out[], int percent, uint8_t *payload) {
    float percent_modulation;

    if ((percent >= 10) && (percent < 100)) {
        percent_modulation = (32768.0f * (((float) percent) / 100.0));
    } else {
        percent_modulation = 16384.0f; /* 50% default */
    }

    int nbytes = horus_get_max_packet_len();            /* payload bytes (22) */

    if (nbytes == 0) {
        return 0;
    }

    int txbytes = horus_get_num_tx_data_bytes(nbytes);
    uint8_t tx[txbytes];    /* payload + overhead (45) */

    ((struct BinaryPacket *) payload)->Checksum = crc16(payload, nbytes - 2); /* no uw */
    horus_encode_packet(tx, payload, nbytes); /* scrambled/interleaved */

    int cycles = api_Fs / api_Rs;
    complex float segment[cycles];
    int i, j, k, m, dibit;

    /*
     * mode is 4fsk in binary
     * This algorithm is brain damaged [SRS]
     */
    for (i = 0; i < txbytes; i++) {
        uint8_t word = tx[i];
        
        for (j = (i * (cycles * 4)), k = 6; k >= 0; j += cycles, k -= 2) {
            dibit = (word >> k) & 0x3; /* MSB pair first */
            
            modulate(api_mod, segment, dibit);

            for (m = 0; m < cycles; m++) {
                demod_out[j + m] = (short) (crealf(segment[m]) * percent_modulation) & 0xFFFF;
            }
        }
    }

    return txbytes * (4 * cycles); /* 4 dibit segments per TX byte */
}

int horus_get_version() {
    return HORUS_API_VERSION;
}

int horus_get_mode() {
    return api_mode;
}

int horus_get_Fs() {
    return api_Fs;
}

int horus_get_Rs() {
    return api_Rs;
}

int horus_get_mFSK() {
    return api_mFSK;
}

int horus_get_nin() {
    return fsk_get_nin(api_fsk);
}

int horus_get_max_demod_in() {
    return fsk_get_Nmem(api_fsk);
}

int horus_get_max_demod_out() {
    int nbytes = horus_get_max_packet_len();            /* payload bytes (22) */
    int txbytes = horus_get_num_tx_data_bytes(nbytes);  /* payload + overhead (45) */
    int cycles = api_Fs / api_Rs;                       /* 480 cycles */

    return txbytes * (4 * cycles);  /* 4 dibit segments per txbyte */
}

int horus_get_max_packet_len() {
    if (api_mode == HORUS_MODE_RTTY) {
        return api_max_packet_len / 10; /* 7 bit ASCII, plus 3 sync bits */
    } else if (api_mode == HORUS_MODE_BINARY) {
        return HORUS_BINARY_NUM_PAYLOAD_BYTES;
    } else {
        return 0; /* bad news */
    }
}

void horus_get_stats(float *snr_est) {
    /* SNR scaled from Eb/No est returned by FSK to SNR in 3000 Hz */

    *snr_est = stats_get_snr_est(api_stats) + 10.0f * log10f(api_Rs / 3000.0f);
}

void horus_get_extended_stats() {
    int i;

    stats_set_snr_est(api_stats, stats_get_snr_est(api_stats) + 10.0f * log10f(api_Rs / 3000.0f));
}

int horus_get_crc_ok() {
    return api_crc_ok;
}

void horus_set_verbose(int verbose) {
    api_verbose = verbose;
}

/**
 * In-Place Algebraic Golden Prime Interleaver
 * 
 * inout  - The byte data to be operated on
 * nbytes - The number of data bytes in the input/output
 * dir    - Direction is 0 to Interleave and 1 to De-Interleave
 */
void horus_interleave(uint8_t *inout, int nbytes, int dir) {
    uint32_t i, j, n, ibit, ibyte, ishift, jbyte, jshift;
    uint32_t b, tmp;
    uint8_t out[nbytes];

    memset(out, 0, nbytes);

    uint16_t imax = sizeof (primes) / sizeof (uint16_t);
    uint16_t nbits = (uint16_t) (nbytes * 8);

    i = 1;
    while ((primes[i] < nbits) && (i < imax))
        i++;

    b = primes[i - 1];  /* b = nearest prime to length of nbits */

    for (n = 0; n < nbits; n++) {
        i = n;
        j = (b * i) % nbits;

        if (dir == DEINTERLEAVE) {
            tmp = j;
            j = i;
            i = tmp;
        }

        ibyte  = (i / 8);
        ishift = (i % 8);
        ibit = (inout[ibyte] >> ishift) & 0x1;

        jbyte  = (j / 8);
        jshift = (j % 8);

        out[jbyte] |= (ibit << jshift);
    }

    memcpy(inout, out, nbytes);
}

/**
 *          In-Place 15 bit additive scrambler with 0x4a80 Frame Sync
 * 
 *  Sync    1   0   0   1   0   1   0   1   0   0   0   0   0   0   0
 *        +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 *        | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11| 12| 13| 14| 15|
 *        +-+-+---+---+---+---+---+---+---+---+---+---+---+---+-+-+-+-+
 *          ^                                                   |   |
 *          |                                                   v   v
 *          |                                                 +-------+
 *          +<------------------------------------------------|   +   |
 *          |                                                 +-------+ 
 *          v
 *        +---+
 * in --->| + |---> out
 *        +---+
 * 
 * inout  - 8-bit unsigned byte input and output data
 * nbytes - integer number of bytes to process
 * 
 * Returns void.
 */

void horus_scramble(uint8_t *inout, int nbytes) {
    int i, ibit, ibits, ibyte, ishift;

    int nbits = (nbytes * 8);

    uint16_t scrambler = HORUS_SCRAMBLE_SEED; /* initialize scrambler every frame */
    uint16_t scrambler_out;

    /* in place modification of each bit */

    for (i = 0; i < nbits; i++) {
        scrambler_out = ((scrambler & 0x2) >> 1) ^ (scrambler & 0x1);

        ibyte  = (i / 8);
        ishift = (i % 8);

        ibit = (inout[ibyte] >> ishift) & 0x1;
        ibits = ibit ^ scrambler_out;

        inout[ibyte] &= ~(1     << ishift);
        inout[ibyte] |=  (ibits << ishift);

        /* update scrambler */

        scrambler >>= 1;
        scrambler  |= (scrambler_out << 14);
    }
}

/* Local Functions */

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

static int horus_find_uw(int n) {
    int i, j, corr, mx, mx_ind;
    int rx_bits_mapped[n + api_uw_len];

    /* map rx_bits to +/-1 for UW search */

    for (i = 0; i < (n + api_uw_len); i++) {
        rx_bits_mapped[i] = 2 * api_rx_bits[i] - 1;
    }

    /* look for UW  */

    mx = 0;
    mx_ind = 0;

    for (i = 0; i < n; i++) {

        /* calculate correlation between bit stream and UW */

        corr = 0;

        for (j = 0; j < api_uw_len; j++) {
            corr += rx_bits_mapped[i + j] * api_uw[j];
        }

        /* peak pick maximum */

        if (corr > mx) {
            mx = corr;
            mx_ind = i;
        }
    }

    if (mx >= api_uw_thresh) {
        return mx_ind;
    } else {
        return -1;
    }
}

static int hex2int(char ch) {
    if (ch >= '0' && ch <= '9')
        return ch - '0';

    if (ch >= 'A' && ch <= 'F')
        return ch - 'A' + 10;

    if (ch >= 'a' && ch <= 'f')
        return ch - 'a' + 10;

    return -1;
}

/*
 * The binary packet format is:
 *
 *    | Unique Word | payload data bits | parity bits |
 *
 * This function works out how much storage the caller of
 * horus_encode_tx_packet() will need to store the tx packet
 */

static int horus_get_num_tx_data_bytes(int num_payload_data_bytes) {
    int num_payload_data_bits = num_payload_data_bytes * 8;
    int num_golay_codewords = num_payload_data_bits / 12;

    if (num_payload_data_bits % 12) /* round up to 12 bits, may mean some unused bits */
        num_golay_codewords++;

    int num_tx_data_bits = (sizeof (l2_uw) * 8) + num_payload_data_bits +
                                   (num_golay_codewords * 11);
    
    int num_tx_data_bytes = (num_tx_data_bits / 8);

    if (num_tx_data_bits % 8) /* round up to nearest byte, may mean some unused bits */
        num_tx_data_bytes++;

    return num_tx_data_bytes;
}

/*
 * Takes an array of payload data bytes, prepends a unique word and appends
 * Golay FEC bits.
 */

static int horus_encode_packet(uint8_t *tx, uint8_t *payload, int nbytes) {
    int inbit, golayparity;
    int ninbyte, shift, golayparitybit;
    int i;

    int txbytes = horus_get_num_tx_data_bytes(nbytes);
    uint8_t *pout = tx;

    memcpy(pout, l2_uw, sizeof (l2_uw));    /* 2 bytes */
    pout += sizeof (l2_uw);

    memcpy(pout, payload, nbytes);
    pout += nbytes;

    /*
     * Bits are written MSB first.
     */

    int num_payload_data_bits = (nbytes * 8);
    int paritybyte = 0;
    int ingolay = 0;
    int ningolay = 0;
    int ninbit = 0;
    int nparitybits = 0;

    while (ninbit < num_payload_data_bits) {

        /* extract input data bit */

        ninbyte = (ninbit / 8);
        shift = 7 - (ninbit % 8);
        inbit = (payload[ninbyte] >> shift) & 0x1;
        ninbit++;

        /* build up input golay codeword */

        ingolay |= inbit;
        ningolay++;

        /* when we get 12 bits do a Golay encode */

        if (ningolay % 12) {
            ingolay <<= 1;
        } else {
            golayparity = golay23_encode(ingolay);
            ingolay = 0;

            /* write parity bits to output data */

            for (i = 0; i < 11; i++) {
                golayparitybit = (golayparity >> (10 - i)) & 0x1;
                paritybyte |= golayparitybit;
                nparitybits++;

                if (nparitybits % 8) {
                    paritybyte <<= 1;
                } else {
                    /* OK we have a full byte ready */
                    *pout = paritybyte;
                    pout++;
                    paritybyte = 0;
                }
            }
        }
    }

    /* Complete final Golay encode, we may have partially finished ingolay, paritybyte */

    if (ningolay % 12) {
        ingolay >>= 1;
        golayparity = golay23_encode(ingolay << 1);
        
        /* write parity bits to output data */

        for (i = 0; i < 11; i++) {
            golayparitybit = (golayparity >> (10 - i)) & 0x1;
            paritybyte = paritybyte | golayparitybit;
            nparitybits++;

            if (nparitybits % 8) {
                paritybyte <<= 1;
            } else {
                /* OK we have a full byte ready */
                *pout++ = (uint8_t) paritybyte;
                paritybyte = 0;
            }
        }
    }

    /* and final, partially complete, parity byte */

    if (nparitybits % 8) {
        paritybyte <<= 7 - (nparitybits % 8); // use MS bits first
        *pout++ = (uint8_t) paritybyte;
    }

    /* we don't interleave or scramble UW */

    horus_interleave(&tx[sizeof (l2_uw)], txbytes - 2, INTERLEAVE);

    /* scrambler to prevent long strings of the same symbol
       which upsets the modem */

    horus_scramble(&tx[sizeof (l2_uw)], txbytes - 2);

    return txbytes;
}

static void horus_decode_packet(uint8_t *output_payload_data,
        uint8_t *input_rx_data,
        int num_payload_data_bytes) {
    int inbit, outbit, outbyte, outdata, ingolay, golayparitybit;
    int ninbyte, shift, i;

    /* we don't interleave UW */

    int num_tx_data_bytes = horus_get_num_tx_data_bytes(num_payload_data_bytes);

    horus_scramble(&input_rx_data[sizeof (l2_uw)], num_tx_data_bytes - 2);
    horus_interleave(&input_rx_data[sizeof (l2_uw)], num_tx_data_bytes - 2, DEINTERLEAVE);

    uint8_t *pin = input_rx_data + sizeof (l2_uw) + num_payload_data_bytes;
    uint8_t *pout = output_payload_data;

    /* Read input data bits one at a time.  When we have 12 read 11 parity bits. Golay decode.
       Write decoded (output data) bits every time we have 8 of them. */

    int num_payload_data_bits = (num_payload_data_bytes * 8);

    int ninbit = 0;
    int ningolay = 0;
    int nparitybits = 0;
    int paritybyte = *pin++;
    int noutbits = 0;
    
    ingolay = 0;
    outbyte = 0;

    while (ninbit < num_payload_data_bits) {

        /* extract input data bit */

        ninbyte = (ninbit / 8) + sizeof (l2_uw);
        shift = 7 - (ninbit % 8);
        inbit = (input_rx_data[ninbyte] >> shift) & 0x1;
        ninbit++;

        /* build up golay codeword */

        ingolay |= inbit;
        ningolay++;
        ingolay <<= 1;

        /* when we get 12 data bits start reading parity bits */

        if ((ningolay % 12) == 0) {
            for (i = 0; i < 11; i++) {
                shift = 7 - (nparitybits % 8);
                golayparitybit = (paritybyte >> shift) & 0x1;
                ingolay |= golayparitybit;

                if (i != 10)
                    ingolay <<= 1;

                nparitybits++;

                if ((nparitybits % 8) == 0) {
                    /* OK grab a new byte */
                    paritybyte = *pin++;
                }
            }

            /* write decoded/error corrected bits to output payload data */

            outdata = golay23_decode(ingolay) >> 11;

            for (i = 0; i < 12; i++) {
                shift = 11 - i;
                outbit = (outdata >> shift) & 0x1;
                outbyte |= outbit;
                noutbits++;

                if (noutbits % 8) {
                    outbyte <<= 1;
                } else {
                    *pout++ = outbyte;
                    outbyte = 0;
                }
            }

            ingolay = 0;
        }
    }

    /* Complete final Golay decode  */

    int golayparity = 0;

    if (ningolay % 12) {
        for (i = 0; i < 11; i++) {
            shift = 7 - (nparitybits % 8);
            golayparitybit = (paritybyte >> shift) & 0x1;
            golayparity |= golayparitybit;

            if (i != 10)
                golayparity <<= 1;

            nparitybits++;

            if ((nparitybits % 8) == 0) {
                /* OK grab a new byte */
                paritybyte = *pin++;
            }
        }

        ingolay >>= 1;

        int codeword = (ingolay << 12) + golayparity;

        outdata = golay23_decode(codeword) >> 11;

        /* write final byte */

        int ntogo = num_payload_data_bits - noutbits;

        for (i = 0; i < ntogo; i++) {
            shift = ntogo - i;
            outbit = (outdata >> shift) & 0x1;
            outbyte |= outbit;
            noutbits++;

            if (noutbits % 8) {
                outbyte <<= 1;
            } else {
                *pout++ = outbyte;
                outbyte = 0;
            }
        }
    }
}

static int horus_extract_rtty(char data_out[], int uw_loc) {
    const int nfield = 7; /* 7 bit ASCII                    */
    const int npad = 3; /* 3 sync bits between characters */
    int st = uw_loc; /* first bit of first char        */
    int en = api_max_packet_len - nfield; /* last bit of max length packet  */
    int i, j;
    char *ptx_crc;

    char *pout = data_out;
    int nout = 0;
    int crc_ok = 0;
    int endpacket = 0;
    uint16_t rx_crc = 0;
    uint16_t tx_crc = 0;

    for (i = st; i < en; i += nfield + npad) {

        /* assemble char LSB to MSB */

        uint8_t char_dec = 0;

        for (j = 0; j < nfield; j++) {
            assert(api_rx_bits[i + j] <= 1);
            char_dec |= api_rx_bits[i + j] * (1 << j);
        }

        if (api_verbose) {
            fprintf(stderr, "i: %4d 0x%02x %c ", i, char_dec, char_dec);

            if ((nout % 6) == 0) {
                fprintf(stderr, "\n");
            }
        }

        /*  if we find a '*' that's the end of the packet for RX CRC calculations */

        if (!endpacket && (char_dec == '*')) {
            endpacket = 1;
            rx_crc = crc16((uint8_t *) & data_out[5], nout - 5);
            ptx_crc = pout + 1; /* start of tx CRC */
        }

        /* build up output array, really only need up to tx crc but
           may end up going further */

        *pout++ = (char) char_dec;
        nout++;
    }

    /* if we found the end of packet flag and have enough chars to compute checksum ... */

    if (endpacket && (pout > (ptx_crc + 3))) {
        tx_crc = 0;

        for (i = 0; i < 4; i++) {
            tx_crc <<= 4;
            tx_crc |= hex2int(ptx_crc[i]);
        }

        crc_ok = (tx_crc == rx_crc);
        *(ptx_crc + 4) = 0; /* terminate ASCII string */
    } else {
        *data_out = '\0';
    }

    if (api_verbose) {
        fprintf(stderr, "\n endpacket: %d nout: %d tx_crc: 0x%04x rx_crc: 0x%04x\n",
                endpacket, nout, tx_crc, rx_crc);
    }

    /* make sure we don't overrun storage */

    assert(nout <= horus_get_max_packet_len());

    api_crc_ok = crc_ok;

    return crc_ok;
}

static int horus_extract_binary(char data_out[], int uw_loc) {
    int b, j;

    int st = uw_loc; /* first bit of first char */
    int en = uw_loc + api_max_packet_len; /* last bit of max length packet  */
    uint8_t rxpacket[api_max_packet_len];

    uint8_t *pout = rxpacket;
    int nout = 0;

    /* convert bits to a packet of bytes */

    for (b = st; b < en; b += 8) { /* 8 bit binary */

        /* assemble bytes MSB to LSB */

        uint8_t rxbyte = 0;

        for (j = 0; j < 8; j++) {
            assert(api_rx_bits[b + j] <= 1);
            rxbyte <<= 1;
            rxbyte |= api_rx_bits[b + j];
        }

        /* build up output array */

        *pout++ = rxbyte;
        nout++;
    }

    if (api_verbose) {
        fprintf(stderr, "nout: %d\nReceived Packet before decoding:\n", nout);

        for (j = 0; j < nout; j++) {
            fprintf(stderr, "%02X", rxpacket[j]);
        }

        fprintf(stderr, "\n");
    }

    uint8_t payload_bytes[HORUS_BINARY_NUM_PAYLOAD_BYTES];
    horus_decode_packet(payload_bytes, rxpacket, HORUS_BINARY_NUM_PAYLOAD_BYTES);

    uint16_t crc_rx = crc16(payload_bytes, HORUS_BINARY_NUM_PAYLOAD_BYTES - 2);
    uint16_t crc_tx = *(uint16_t *) &payload_bytes[HORUS_BINARY_NUM_PAYLOAD_BYTES - 2];

    if (api_verbose) {
        fprintf(stderr, "crc_tx: %04X crc_rx: %04X\n", crc_tx, crc_rx);
    }

    /* convert to ASCII string of hex characters */

    data_out[0] = 0;
    char hex[3];

    for (j = 0; j < HORUS_BINARY_NUM_PAYLOAD_BYTES; j++) {
        sprintf(hex, "%02X", payload_bytes[j]);
        strcat(data_out, hex);
    }

    if (api_verbose) {
        fprintf(stderr, "nout: %d\nDecoded Payload bytes:\n%s", nout, data_out);
    }

    api_crc_ok = (crc_tx == crc_rx);

    /* binary packets always marked as OK, as next layer determines validity */

    return 1;
}
