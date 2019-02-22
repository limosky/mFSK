#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <complex.h>

#include "horus_api.h"

static struct BinaryPacket packet;
static uint8_t *payload;

int main() {
    FILE *fout;
    int i;

    packet.PayloadID = 25;
    packet.Counter = 0;
    packet.Hours = 0;
    packet.Minutes = 0;
    packet.Seconds = 0;
    packet.Latitude = 35.0;
    packet.Longitude = -97.0;
    packet.Altitude = 400;
    packet.Speed = 0;
    packet.Sats = 3;
    packet.Temp = 128;
    packet.BattVoltage = 128;
    packet.Checksum = 0;

    payload = (uint8_t *) &packet;

    if (horus_open(HORUS_MODE_BINARY) == 0) {
        fprintf(stderr, "Unable to create Horus instance\n");
        return 1;
    }

    if ((fout = fopen("/tmp/horus-test.raw", "wb")) == NULL) {
        fprintf(stderr, "Error opening output waveform file\n");
        return 2;
    }

    int length = horus_get_max_demod_out(); /* 4 * 45 * 480 = 86400 */
    short word[length];

    for (i = 0; i < 10; i++) {
        horus_binary_tx(word, 50, payload);

        fwrite(word, sizeof (short), length, fout);
        packet.Counter = packet.Counter + 1;
    }

    fflush(fout);
    fclose(fout);

    horus_close();

    return 0;
}
