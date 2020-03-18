/*
 * Copyright (c) 2003-2004, Mark Borgerding
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or other
 * materials provided with the distribution.
 *
 * Neither the author nor the names of any contributors may be used to endorse or
 * promote products derived from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Modified March 2020 by Steve Sampson
 */

#include "fft.h"
#include "fsk.h"

/* Static prototypes/Forward declarations */

static void kf_bfly2(complex float *, const size_t, const fft_cfg, int);
static void kf_bfly4(complex float *, const size_t, const fft_cfg, const size_t);
static void kf_bfly3(complex float *, const size_t, const fft_cfg, size_t);
static void kf_bfly5(complex float *, const size_t, const fft_cfg, int);
static void kf_bfly_generic(complex float *, const size_t, const fft_cfg, int, int);
static void kf_work(complex float *, const complex float *, const size_t, int *, const fft_cfg);
static void kf_factor(int, int *);
static void fft_stride(fft_cfg, const complex float *, complex float *);

/* Public Functions */

fft_cfg fft_alloc(int nfft, int inverse_fft, void *mem, size_t *lenmem) {
    fft_cfg st = NULL;
    size_t memneeded = sizeof (struct fft_state)
            + sizeof (complex float) * (nfft - 1); /* twiddle factors*/

    if (lenmem == NULL) {
        st = (fft_cfg) malloc(memneeded);
    } else {
        if (mem != NULL && *lenmem >= memneeded)
            st = (fft_cfg) mem;

        *lenmem = memneeded;
    }

    if (st) {
        st->nfft = nfft;
        st->inverse = inverse_fft;

        for (int i = 0; i < nfft; i++) {
            float phase = -TAU * (float) i / (float) nfft;

            if (inverse_fft)
                phase *= -1.0f;

            *(st->twiddles + i) = cmplx(phase);
        }

        kf_factor(nfft, st->factors);
    }

    return st;
}

fftr_cfg fftr_alloc(int nfft, int inverse_fft, void *mem, size_t *lenmem) {
    fftr_cfg st = NULL;
    size_t subsize;

    if (nfft & 1) {
        return NULL;
    }

    nfft >>= 1;

    fft_alloc(nfft, inverse_fft, NULL, &subsize);
    size_t memneeded = sizeof (struct fftr_state) + subsize + sizeof (complex float) * (nfft * 3 / 2);

    if (lenmem == NULL) {
        st = (fftr_cfg) malloc(memneeded);
    } else {
        if (*lenmem >= memneeded) {
            st = (fftr_cfg) mem;
        }

        *lenmem = memneeded;
    }

    if (!st) {
        return NULL;
    }

    st->substate = (fft_cfg) (st + 1); /*just beyond fftr_state struct */
    st->tmpbuf = (complex float *) (((char *) st->substate) + subsize);
    st->super_twiddles = st->tmpbuf + nfft;

    fft_alloc(nfft, inverse_fft, st->substate, &subsize);

    for (int i = 0; i < nfft / 2; ++i) {
        float phase = -M_PI * ((float) (i + 1) / (float) nfft + .5f);

        if (inverse_fft) {
            phase *= -1.0f;
        }

        *(st->super_twiddles + i) = cmplx(phase);
    }

    return st;
}

/* Complex FFT */

void fft(fft_cfg cfg, const complex float *fin, complex float *fout) {
    fft_stride(cfg, fin, fout);
}

/* Real FFT Forward */

void encode_fftr(fftr_cfg st, const float *timedata, complex float *freqdata) {
    complex float fpnk, fpk, f1k, f2k, tw, tdc;

    fft(st->substate, (const complex float *) timedata, st->tmpbuf);

    int ncfft = st->substate->nfft;

    tdc = st->tmpbuf[0];

    freqdata[0] = (crealf(tdc) + cimagf(tdc)) + 0.0f * I;
    freqdata[ncfft] = (crealf(tdc) - cimagf(tdc)) + 0.0f * I;

    for (int k = 1; k <= (ncfft / 2); k++) {
        fpk = st->tmpbuf[k];
        fpnk = conjf(st->tmpbuf[ncfft - k]);

        f1k = fpk + fpnk;
        f2k = fpk - fpnk;
        tw = f2k * st->super_twiddles[k - 1];

        freqdata[k] = ((crealf(f1k) + crealf(tw)) * .5f) + ((cimagf(f1k) + cimagf(tw)) * .5f) * I;
        freqdata[ncfft - k] = ((crealf(f1k) - crealf(tw)) * .5f) + ((cimagf(tw) - cimagf(f1k)) * .5f) * I;
    }
}

/* Real FFT Inverse */

void encode_fftri(fftr_cfg st, const complex float *freqdata, float *timedata) {
    complex float fk, fnkc, fek, fok;

    int ncfft = st->substate->nfft;

    st->tmpbuf[0] = (crealf(freqdata[0]) + crealf(freqdata[ncfft])) +
            (crealf(freqdata[0]) - crealf(freqdata[ncfft])) * I;

    for (int k = 1; k <= (ncfft / 2); k++) {
        fk = freqdata[k];
        fnkc = conjf(freqdata[ncfft - k]);

        fek = fk + fnkc;
        fok = (fk - fnkc) * st->super_twiddles[k - 1];

        st->tmpbuf[k] = fek + fok;
        st->tmpbuf[ncfft - k] = conjf(fek - fok);
    }

    fft_stride(st->substate, st->tmpbuf, (complex float *) timedata);
}

/* Local Functions */

static void kf_bfly2(
        complex float *Fout,
        const size_t fstride,
        const fft_cfg st,
        int m) {
    complex float *Fout2;
    complex float *tw1 = st->twiddles;
    complex float t;

    Fout2 = Fout + m;

    do {
        t = *Fout2 * *tw1;
        tw1 += fstride;

        *Fout2 = *Fout - t;
        *Fout = *Fout + t;

        Fout2++;
        Fout++;
    } while (--m);
}

static void kf_bfly4(
        complex float *Fout,
        const size_t fstride,
        const fft_cfg st,
        const size_t m) {
    complex float *tw1, *tw2, *tw3;
    complex float scratch[6];
    size_t k = m;
    const size_t m2 = 2 * m;
    const size_t m3 = 3 * m;

    tw3 = tw2 = tw1 = st->twiddles;

    do {
        scratch[0] = Fout[m] * *tw1;
        scratch[1] = Fout[m2] * *tw2;
        scratch[2] = Fout[m3] * *tw3;

        scratch[5] = *Fout - scratch[1];
        *Fout = *Fout + scratch[1];

        scratch[3] = scratch[0] + scratch[2];
        scratch[4] = scratch[0] - scratch[2];

        Fout[m2] = *Fout - scratch[3];

        tw1 += fstride;
        tw2 += fstride * 2;
        tw3 += fstride * 3;

        *Fout = *Fout + scratch[3];

        if (st->inverse) {
            Fout[m] = (crealf(scratch[5]) - cimagf(scratch[4])) + (cimagf(scratch[5]) + crealf(scratch[4])) * I;
            Fout[m3] = (crealf(scratch[5]) + cimagf(scratch[4])) + (cimagf(scratch[5]) - crealf(scratch[4])) * I;
        } else {
            Fout[m] = (crealf(scratch[5]) + cimagf(scratch[4])) + (cimagf(scratch[5]) - crealf(scratch[4])) * I;
            Fout[m3] = (crealf(scratch[5]) - cimagf(scratch[4])) + (cimagf(scratch[5]) + crealf(scratch[4])) * I;
        }

        Fout++;
    } while (--k);
}

static void kf_bfly3(
        complex float *Fout,
        const size_t fstride,
        const fft_cfg st,
        size_t m) {
    size_t k = m;
    const size_t m2 = 2 * m;
    complex float scratch[5];
    complex float epi3 = st->twiddles[fstride * m];

    complex float *tw1 = st->twiddles;
    complex float *tw2 = st->twiddles;

    do {
        scratch[1] = Fout[m] * *tw1;
        scratch[2] = Fout[m2] * *tw2;

        scratch[3] = scratch[1] + scratch[2];
        scratch[0] = scratch[1] - scratch[2];

        tw1 += fstride;
        tw2 += fstride * 2;

        Fout[m] = (crealf(*Fout) - (crealf(scratch[3]) * .5f)) + (cimagf(*Fout) - (cimagf(scratch[3]) * .5f)) * I;

        scratch[0] = scratch[0] * cimagf(epi3);
        *Fout = *Fout + scratch[3];

        Fout[m2] = (crealf(Fout[m]) + cimagf(scratch[0])) + (cimagf(Fout[m]) - crealf(scratch[0])) * I;
        Fout[m] = (crealf(Fout[m]) - cimagf(scratch[0])) + (cimagf(Fout[m]) + crealf(scratch[0])) * I;

        Fout++;
    } while (--k);
}

static void kf_bfly5(
        complex float *Fout,
        const size_t fstride,
        const fft_cfg st,
        int m) {
    complex float *Fout0, *Fout1, *Fout2, *Fout3, *Fout4;
    complex float scratch[13];
    complex float *tw = st->twiddles;
    complex float ya = tw[fstride * m];
    complex float yb = tw[fstride * m * 2];

    Fout0 = Fout;
    Fout1 = Fout0 + m;
    Fout2 = Fout0 + 2 * m;
    Fout3 = Fout0 + 3 * m;
    Fout4 = Fout0 + 4 * m;

    for (int u = 0; u < m; u++) {
        scratch[0] = *Fout0;

        scratch[1] = *Fout1 * tw[fstride * u];
        scratch[2] = *Fout2 * tw[fstride * 2 * u];
        scratch[3] = *Fout3 * tw[fstride * 3 * u];
        scratch[4] = *Fout4 * tw[fstride * 4 * u];

        scratch[7] = scratch[1] + scratch[4];
        scratch[10] = scratch[1] - scratch[4];
        scratch[8] = scratch[2] + scratch[3];
        scratch[9] = scratch[2] - scratch[3];

        *Fout0 += (crealf(scratch[7]) + crealf(scratch[8])) + (cimagf(scratch[7]) + cimagf(scratch[8])) * I;

        scratch[5] = (crealf(scratch[0]) + (crealf(scratch[7]) * crealf(ya)) + (crealf(scratch[8]) * crealf(yb))) + (cimagf(scratch[0]) + (cimagf(scratch[7]) * crealf(ya)) + (cimagf(scratch[8]) * crealf(yb))) * I;
        scratch[6] = -((cimagf(scratch[10]) * cimagf(ya)) + (cimagf(scratch[9]) * cimagf(yb))) - ((crealf(scratch[10]) * cimagf(ya)) - (crealf(scratch[9]) * cimagf(yb))) * I;

        *Fout1 = scratch[5] - scratch[6];
        *Fout4 = scratch[5] + scratch[6];

        scratch[11] = (crealf(scratch[0]) + (crealf(scratch[7]) * crealf(yb)) + (crealf(scratch[8]) * crealf(ya))) + (cimagf(scratch[0]) + (cimagf(scratch[7]) * crealf(yb)) + (cimagf(scratch[8]) * crealf(ya))) * I;
        scratch[12] = -((cimagf(scratch[10]) * cimagf(yb)) + (cimagf(scratch[9]) * cimagf(ya))) + ((crealf(scratch[10]) * cimagf(yb)) - (crealf(scratch[9]) * cimagf(ya))) * I;

        *Fout2 = scratch[11] + scratch[12];
        *Fout3 = scratch[11] - scratch[12];

        Fout0++;
        Fout1++;
        Fout2++;
        Fout3++;
        Fout4++;
    }
}

/* perform the butterfly for one stage of a mixed radix FFT */

static void kf_bfly_generic(
        complex float *Fout,
        const size_t fstride,
        const fft_cfg st,
        int m,
        int p) {
    complex float *tw = st->twiddles;
    complex float t;
    int Norig = st->nfft;

    complex float *scratch = (complex float *) malloc(sizeof (complex float) * p);

    for (int u = 0; u < m; u++) {
        int k = u;

        for (int q1 = 0; q1 < p; q1++) {
            scratch[q1] = Fout[ k ];
            k += m;
        }

        k = u;
        for (int q1 = 0; q1 < p; q1++) {
            int twidx = 0;
            Fout[ k ] = scratch[0];

            for (int q = 1; q < p; q++) {
                twidx += fstride * k;

                if (twidx >= Norig)
                    twidx -= Norig;

                t = scratch[q] * tw[twidx];
                Fout[ k ] = Fout[ k ] + t;
            }

            k += m;
        }
    }

    free(scratch);
}

static void kf_work(
        complex float *Fout,
        const complex float *f,
        const size_t fstride,
        int *factors,
        const fft_cfg st) {
    complex float *Fout_beg = Fout;
    const int p = *factors++; /* the radix  */
    const int m = *factors++; /* stage's fft length/p */
    const complex float *Fout_end = Fout + p * m;

    if (m == 1) {
        do {
            *Fout = *f;
            f += fstride;
        } while (++Fout != Fout_end);
    } else {
        do {
            kf_work(Fout, f, fstride*p, factors, st);

            f += fstride;
        } while ((Fout += m) != Fout_end);
    }

    Fout = Fout_beg;

    // recombine the p smaller DFTs
    switch (p) {
        case 2:
            kf_bfly2(Fout, fstride, st, m);
            break;
        case 3:
            kf_bfly3(Fout, fstride, st, m);
            break;
        case 4:
            kf_bfly4(Fout, fstride, st, m);
            break;
        case 5:
            kf_bfly5(Fout, fstride, st, m);
            break;
        default:
            kf_bfly_generic(Fout, fstride, st, m, p);
    }
}

static void kf_factor(int n, int *facbuf) {
    int p = 4;
    float floor_sqrt = floorf(sqrtf((float) n));

    /*factor out powers of 4, powers of 2, then any remaining primes */
    do {
        while (n % p) {
            switch (p) {
                case 4:
                    p = 2;
                    break;
                case 2:
                    p = 3;
                    break;
                default:
                    p += 2;
            }

            if (p > floor_sqrt)
                p = n; /* no more factors, skip to end */
        }

        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    } while (n > 1);
}

static void fft_stride(fft_cfg st, const complex float *fin, complex float *fout) {
    if (fin == fout) {
        complex float *tmpbuf = (complex float *) malloc(sizeof (complex float) * st->nfft);

        kf_work(tmpbuf, fin, 1, st->factors, st);
        memcpy(fout, tmpbuf, sizeof (complex float) * st->nfft);

        free(tmpbuf);
    } else {
        kf_work(fout, fin, 1, st->factors, st);
    }
}
