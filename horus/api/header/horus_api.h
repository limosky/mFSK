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

#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdint.h>

#define HORUS_API_VERSION                1    /* unique number that is bumped if API changes */

#define HORUS_MODE_BINARY                0
#define HORUS_MODE_RTTY                  1
    
#define MAX_UW_LENGTH                  100

#define HORUS_BINARY_NUM_BITS          360    /* fixed number of bits in binary frame    */
#define HORUS_BINARY_NUM_PAYLOAD_BYTES  22    /* fixed number of bytes in binary payload */

#define HORUS_SCRAMBLE_SEED         0x4a80    /* LFSR initial value */

#define INTERLEAVE                       0
#define DEINTERLEAVE                     1

/* Horus binary packet (22 bytes) */

struct BinaryPacket
{
    uint8_t     PayloadID;   /* User Defined */
    uint16_t	Counter;     /* User Defined */
    uint8_t	Hours;       /* GPS Hours */
    uint8_t	Minutes;     /* GPS Minutes */
    uint8_t	Seconds;     /* GPS Seconds */
    float	Latitude;    /* GPS Latitude (4 byte float) */
    float	Longitude;   /* GPS Longitude (4 byte float) */
    uint16_t  	Altitude;    /* GPS Altitude in Feet */
    uint8_t     Speed;       /* Ground Speed in Knots (1-255 knots) */
    uint8_t     Sats;        /* Number of GPS satellites in view */
    int8_t      Temp;        /* Twos Complement Temperature value. */
    uint8_t     BattVoltage; /* 0 = 0.5v, 255 = 2.0V, linear steps in-between */
    uint16_t    Checksum;    /* CRC16-CCITT Checksum */
}  __attribute__ ((packed));

/**
 * Initialize API Variables and open and FSK session
 * 
 * hmode - An integer to specify either RTTY or Binary mode.
 * 
 * Returns 1 on success and 0 on error.
 */
int horus_open(int);

/**
 * Close the FSK session and free memory.
 */
void horus_close(void);

/**
 * data_out  - ASCII or HEX data output signal
 * demod_in  - 16-bit short for PCM RAW input signal
 * 
 * Returns 1 for RTTY checksum good, or 1 for binary data detected,
 * or else returns 0.
 */
int horus_rx(char [], short []);

/**
 * This function is for Binary packets
 * 
 * demod_out - 16 bit short for PCM RAW output signal
 * percent   - integer between 10 and 100 for signal amplitude
 * payload   - a pointer to unsigned character data to be transmitted
 * 
 * Returns an integer of PCM demodulator words in demod_out
 */
int horus_binary_tx(short [], int, uint8_t *);

void horus_set_verbose(int);

int horus_get_version(void);
int horus_get_mode(void);
int horus_get_Fs(void);
int horus_get_Rs(void);
int horus_get_mFSK(void);
int horus_get_nin(void);

void horus_get_stats(float *);
void horus_get_extended_stats(void);

int horus_get_crc_ok(void);
int horus_get_max_demod_in(void);
int horus_get_max_demod_out(void);
int horus_get_max_packet_len(void);

/**
 * In-Place Algebraic Golden Prime Interleaver
 * 
 * inout  - The byte data to be operated on
 * nbytes - The number of data bytes in the input/output
 * dir    - Direction is 0 to Interleave and 1 to De-Interleave
 */
void horus_interleave(uint8_t *, int, int);

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
void horus_scramble(uint8_t *, int);

#ifdef __cplusplus
}
#endif
