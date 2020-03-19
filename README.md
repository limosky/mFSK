### mFSK Modem
A 2FSK and 4FSK Modem based on the FSK Modem designed by David Rowe. I have changed it to a Dynamic Library. I'm using the Netbeans IDE and have included the export so you can import and compile any changes. You can also just burst the ZIP's (for example on a Raspberry Pi) and type ```make```.

#### Build Instructions
First copy the ```libfsk.so``` file to ```/usr/local/lib``` Make sure the permissions are set correctly (644).
```
sudo cp libfsk.so /usr/local/lib
sudo ldconfig
```
Then copy all the header files to /usr/local/include
```
sudo cp header/* /usr/local/include
```
To compile the test file
```
c99 mod-test.c -o mod-test -lfsk -lm
```
Executing the ```mod-test``` program will create a RAW audio file ```/tmp/fsksignal.raw``` and you can use ```audacity``` to view the spectrum, or play with the audio.


#### Examples using Audacity
Included is a couple of PNG graphics from the modulator tests, and the resulting RAW audio samples. The audio samples are 16-bit signed PCM, little-endian, at an 8000 sample rate.

The test signal has an F1 frequency of 1100 Hz, a shift frequency of 200 Hz, and a symbol rate of 100 Baud.

In 2FSK this results in a 0-bit producing 1100 Hz, and a 1-bit producing 1300 Hz.

<img src="https://raw.githubusercontent.com/srsampson/mFSK/master/2fsk.png" width="400">

In 4FSK this results in a 00-bits producing 1100 Hz, 01-bits producing 1300 Hz, 10-bits producing 1500 Hz, and 11-bits producing 1700 Hz

<img src="https://raw.githubusercontent.com/srsampson/mFSK/master/4fsk.png" width="400">

Just for kicks I added a Manchester Modulator option. I don't have a demodulator, but I was interested in seeing what the spectrum looked like.

#### Project Development
Currently the modem compiles without error, and the modulator/demodulator seems to work.

There is a little testing program I used for the CRC, Scrambler, and Interleaver actually work. It shows they produce the right data, and is bidirectional.

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
