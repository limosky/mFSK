/*
 * Copyright (C) April 2018 David Rowe
 *
 * All rights reserved.
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include "fsk.h"
#include "demodulator.h"
#include "horus_api.h"

static struct option long_opts[] = {
    {"help", no_argument, 0, 'h'},
    {"mode", required_argument, 0, 'm'},
    {"stats", optional_argument, 0, 't'},
    {0, 0, 0, 0}
};

int main(int argc, char *argv[]) {
    struct STATS *stats;
    FILE *fin = NULL, *fout = NULL;
    int i, j, Ndft, mode;
    int stats_ctr, stats_loop, stats_rate, verbose, crc_results;
    float loop_time;
    int enable_stats = 0;

    if ((stats = stats_open()) == NULL) {
        return 0;
    }

    stats_loop = 0;
    stats_rate = 8;
    mode = -1;
    verbose = crc_results = 0;

    int opt = 0;
    int opt_idx = 0;

    while (opt != -1) {
        opt = getopt_long(argc, argv, "hvcm:t::", long_opts, &opt_idx);

        switch (opt) {
            case 'm':
                if ((strcmp(optarg, "RTTY") == 0) || (strcmp(optarg, "rtty") == 0)) {
                    mode = HORUS_MODE_RTTY;
                }

                if ((strcmp(optarg, "BINARY") == 0) || (strcmp(optarg, "binary") == 0)) {
                    mode = HORUS_MODE_BINARY;
                }

                if (mode == -1) {
                    fprintf(stderr, "use --mode RTTY or --mode binary\n");
                    exit(1);
                }
                break;
            case 't':
                enable_stats = 1;

                if (optarg != NULL) {
                    stats_rate = atoi(optarg);

                    if (stats_rate == 0) {
                        stats_rate = 8;
                    }
                }
                break;
            case 'v':
                verbose = 1;
                break;
            case 'c':
                crc_results = 1;
                break;
            case 'h':
            case '?':
                goto helpmsg;
        }
    }

    int dx = optind;

    if ((argc - dx) < 1) {
        fprintf(stderr, "Too few arguments\n");
        goto helpmsg;
    }

    if ((argc - dx) > 5) {
        fprintf(stderr, "Too many arguments\n");
helpmsg:
        fprintf(stderr, "usage: %s -m RTTY|binary [-v] [-c] [-t [r]] InputModemRawFile OutputAsciiFile\n", argv[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "InputModemRawFile      48 kHz 16 bit shorts real modem signal from radio\n");
        fprintf(stderr, " -m RTTY|binary\n");
        fprintf(stderr, "--mode=RTTY|binary[r]  RTTY or binary Horus protcols\n");
        fprintf(stderr, " -t[r] --stats=[r]     Print out modem statistics to stderr in JSON.\n");
        fprintf(stderr, "                       r, if provided, sets the number of modem frames\n"
                "                       between statistic printouts\n");
        fprintf(stderr, " -v                    verbose debug info\n");
        fprintf(stderr, " -c                    display CRC results for each packet\n");
        exit(1);
    }

    /* Open files */

    if (verbose) {
        fprintf(stderr, "mode: %d verbose: %d stats_loop: %d stats_rate: %d\n",
                mode, verbose, stats_loop, stats_rate);
    }

    if (strcmp(argv[dx], "-") == 0) {
        fin = stdin;
    } else {
        fin = fopen(argv[dx], "rb");
    }

    if (strcmp(argv[dx + 1], "-") == 0) {
        fout = stdout;
    } else {
        fout = fopen(argv[dx + 1], "w");
    }

    if ((fin == NULL) || (fout == NULL)) {
        fprintf(stderr, "Couldn't open test vector files\n");
        exit(1);
    }

    /* end command line processing */

    if (horus_open(mode) == 0) {
        fprintf(stderr, "Couldn't open Horus API\n");
        exit(1);
    }

    horus_set_verbose(verbose);

    if (enable_stats) {
        loop_time = (float) horus_get_nin() / horus_get_Fs();
        stats_loop = (int) (1.0f / (stats_rate * loop_time));
        stats_ctr = 0;
    }

    int max_demod_in = horus_get_max_demod_in();
    short demod_in[max_demod_in];
    int max_ascii_out = horus_get_max_packet_len();
    char ascii_out[max_ascii_out];

    /* Main loop ----------------------------------------------------------------------- */

    int nin = horus_get_nin();

    while (fread(demod_in, sizeof (short), nin, fin) == nin) {
        if (verbose) {
            fprintf(stderr, "read nin %d\n", nin);
        }

        if (horus_rx(ascii_out, demod_in)) {
            fprintf(fout, "%s", ascii_out);

            if (crc_results) {
                if (horus_get_crc_ok()) {
                    fprintf(fout, "  CRC OK");
                } else {
                    fprintf(fout, "  CRC BAD");
                }
            }

            fprintf(fout, "\n");
        }

        if (enable_stats && stats_ctr <= 0) {
            horus_get_extended_stats();

            /* Print standard 2FSK stats */

            fprintf(stderr, "{\"EbNodB\": %2.2f,\t\"ppm\": %d,", stats_get_snr_est(stats), (int) stats_get_clock_offset(stats));

            stats_ctr = stats_loop;
        }

        stats_ctr--;

        if (fin == stdin || fout == stdin) {
            fflush(fin);
            fflush(fout);
        }

        nin = horus_get_nin();
    }

    horus_close();

    return 0;
}
