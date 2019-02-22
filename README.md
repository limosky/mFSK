### mFSK Modem
A 2FSK and 4FSK Modem based on the Horus FSK Modem designed by David Rowe and Brady O'Brien. I have changed things around a bit, by using a Dynamic Library for the fsk and horus API's. I'm using the Netbeans IDE and have included the export so you can import and compile any changes. You can also just burst the ZIP's (for example on a Raspberry Pi) and type ```make``` for each of the packages.

The KISS FFT has been highly modified to use the C99 complex syntax instead of all the macros and artificial complex structure.

I also added in the "real" forward, and inverse KISS FFT functions that are also converted to the complex syntax. 

#### Build Instructions
First copy the ```libfsk.so``` and ```libhorus_api.so``` files to ```/usr/local/lib``` Make sure the permissions are set correctly (644).
```
sudo cp libfsk.so /usr/local/lib
sudo cp libhorus_api.so /usr/local/lib
sudo ldconfig
```
Then copy all the header files to /usr/local/include
```
sudo cp header/* /usr/local/include
```
Then compile the test file
```
c99 mod-test.c -o mod-test -lfsk -lm
c99 horus_demod.c -o horus_demod -lhorus_api -lfsk -lm
c99 bin-test.c -o bin-test -lhorus_api -lfsk -lm
```
Executing the ```mod-test``` program will create a RAW audio file ```/tmp/fsksignal.raw``` and you can use ```audacity``` to view the spectrum, or play with the audio.

Executing the ```horus_demod``` test program reads in a RAW audio file and attempts to decode it. You should get something like this:
```
$ horus_demod -m binary -c binary.raw -
0112000000230000000000000000000000001C9A9545  CRC OK
266548A2EB3C960E62156D114BD97BD285C7041D2C19  CRC BAD
0112000000230000000000000000000000001C9A9545  CRC OK
```

Also included is a ```bin-test``` program to exercise the modulator code. It sends the packet data structure, and increments its counter:

```
PayloadID = 25;
Counter = 0;
Hours = 0;
Minutes = 0;
Seconds = 0;
Latitude = 35.0;
Longitude = -97.0;
Altitude = 400;
Speed = 0;
Sats = 3;
Temp = 128;
BattVoltage = 128;
Checksum = 0;
```

Executing the ```bin-test``` program will write 10 packets to ```/tmp/horus-test.raw``` and you can use ```horus_demod``` to see if you can decode it. This is what I get, and it is weird because I increment the Counter word, and you can see the Counter word incrementing (second and third bytes). The line with the CRC Error seems to be buffer garbage, or I have a bug.

```
$ ./horus_demod -m binary -c /tmp/horus-test.raw -
19000000000000000C420000C2C29001000380804FC2  CRC OK
74436DB4DF771CAEC8370107FD6319535659A1D5A408  CRC BAD
19010000000000000C420000C2C2900100038080501C  CRC OK
19020000000000000C420000C2C2900100038080506E  CRC OK
19030000000000000C420000C2C29001000380804FB0  CRC OK
74436DB4DF771E27D0370405FD6319535655CBD5A4CA  CRC BAD
19040000000000000C420000C2C2900100038080508A  CRC OK
19050000000000000C420000C2C29001000380804F54  CRC OK
74436DB4DF3D8EAD4076C107FD631942F65D09D563CA  CRC BAD
19060000000000000C420000C2C29001000380804F26  CRC OK
74436DB4DF771E20D876C107FD6319535659F9D5D648  CRC BAD
19070000000000000C420000C2C290010003808050F8  CRC OK
19080000000000000C420000C2C29001000380807152  CRC OK
```

#### Examples using Audacity
Included is a couple of PNG graphics from the modulator tests, and the resulting RAW audio samples. The audio samples are 16-bit signed PCM, little-endian, at an 8000 sample rate.

The test signal has an F1 frequency of 1100 Hz, a shift frequency of 200 Hz, and a symbol rate of 100 Baud.

In 2FSK this results in a 0-bit producing 1100 Hz, and a 1-bit producing 1300 Hz.

<img src="https://raw.githubusercontent.com/srsampson/mFSK/master/2fsk.png" width="400">

In 4FSK this results in a 00-bits producing 1100 Hz, 01-bits producing 1300 Hz, 10-bits producing 1500 Hz, and 11-bits producing 1700 Hz

<img src="https://raw.githubusercontent.com/srsampson/mFSK/master/4fsk.png" width="400">

Just for kicks I added a Manchester Modulator option. I don't have a demodulator, but I was interested in seeing what the spectrum looked like.

#### Project Development
Currently the modem compiles without error, and the modulator/demodulator seems to work. I am converting it to use internal data structures rather than static BSS storage. So maybe the docs have strange things that don't add up.

There is a little testing program I used to make sure the CRC, Scrambler, and Interleaver actually work. It showed that they produced the right data, and they were bidirectional.

#### Library Calls
In order to use the FSK library, the following functions are available, and the internal data structure is shown in structures:
```
struct FSK {
    complex float phase[MAX_TONES];
    float f_est[MAX_TONES];
    float norm_rx_timing;
    float EbNodB;
    float ppm;
    int Ndft;
    int Fs;
    int N;
    int Rs;
    int Ts;
    int Nmem;
    int P;
    int Nsym;
    int Nbits;
    int f1;
    int shift;
    int mode;
    int est_min;
    int est_max;
    int est_space;
    int nstash;
    int nin;
    bool burst_mode;
    fft_cfg fftcfg;
    float *hann_table;
    float *fft_est;
    complex float *samp_old;
};

struct MODULATE {
    complex float oscillator[MAX_TONES];
    complex float phase;
    float fs;
    float rs;
    float f1;
    float shift;
    int mode;
    int cycles;
};

struct STATS {
    fft_cfg fftcfg;
    float snr_est; /* estimated SNR of rx signal in dB (3 kHz noise BW) */
    float foff; /* estimated freq offset in Hz */
    float rx_timing; /* estimated optimum timing offset in samples */
    float clock_offset; /* Estimated tx/rx sample clock offset in ppm */
    float scale_dB;
    float fftbuffer[STATS_FFTSIZE];
};
```
#### Stats
```
struct STATS *stats_open(void);
void stats_close(struct STATS *);
void stats_spectrum(struct STATS *, float [], complex float [], int);

float stats_get_snr_est(struct STATS *);
float stats_get_foff(struct STATS *);
float stats_get_rx_timing(struct STATS *);
float stats_get_clock_offset(struct STATS *);

float *stats_get_fft_buf_ptr(struct STATS *);

void stats_set_snr_est(struct STATS *, float);
void stats_set_foff(struct STATS *, float);
void stats_set_rx_timing(struct STATS *, float);
void stats_set_clock_offset(struct STATS *, float);
```
#### Fsk
```
struct FSK *fsk_create(int, int, int, int);
struct FSK *fsk_create_hbr(int, int, int, int, int, int);
void fsk_destroy(struct FSK *);

int fsk_get_nin(struct FSK *);
int fsk_get_N(struct FSK *);
int fsk_get_Nmem(struct FSK *);
int fsk_get_Nbits(struct FSK *);
int fsk_get_Ts(struct FSK *);
float fsk_get_f_est(struct FSK *, int);
void fsk_set_nsym(struct FSK *, int);
void fsk_set_est_limits(struct FSK *, int, int);
void fsk_set_estimators(struct FSK *);
```
#### Modulate
```
struct MODULATE *mod_create(int, int, int, int, int);
void mod_destroy(struct MODULATE *);
void modulate(struct MODULATE *, complex float [], int);
void manchester_modulate(struct MODULATE *, complex float [], int);
```
#### Demodulate
```
void fsk_demod(struct FSK *, uint8_t [], complex float []);
void fsk_demod_sd(struct FSK *, float [], complex float []);
```
