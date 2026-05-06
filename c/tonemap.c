#include "tonemap.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* BT.601 luma weights for RGB */
#define LUM_R 0.299f
#define LUM_G 0.587f
#define LUM_B 0.114f

static int cmp_float(const void *a, const void *b)
{
    float fa = *(const float *)a;
    float fb = *(const float *)b;
    return (fa > fb) - (fa < fb);
}

void reinhard_tonemap(const float *irr, int W, int H,
                      float alpha, float gamma, float wb_strength,
                      uint8_t *out)
{
    int N = W * H;

    float *E     = (float *)malloc((size_t)N * 3 * sizeof(float));
    float *L     = (float *)malloc((size_t)N *     sizeof(float));
    float *T     = (float *)malloc((size_t)N *     sizeof(float));
    float *Ltone = (float *)malloc((size_t)N *     sizeof(float));

    /* Optional gray-world WB.  Keep wb_strength=0 for exact Python behavior. */
    float gain_r = 1.0f, gain_g = 1.0f, gain_b = 1.0f;
    if (wb_strength > 0.0f) {
        double mean_r = 0.0, mean_g = 0.0, mean_b = 0.0;
        for (int p = 0; p < N; p++) {
            mean_r += irr[p * 3 + 0];
            mean_g += irr[p * 3 + 1];
            mean_b += irr[p * 3 + 2];
        }
        mean_r /= N;
        mean_g /= N;
        mean_b /= N;
        double mean_all = (mean_r + mean_g + mean_b) / 3.0;
        float full_r = (mean_r > 1e-12) ? (float)(mean_all / mean_r) : 1.0f;
        float full_g = (mean_g > 1e-12) ? (float)(mean_all / mean_g) : 1.0f;
        float full_b = (mean_b > 1e-12) ? (float)(mean_all / mean_b) : 1.0f;

        if (wb_strength < 0.0f) wb_strength = 0.0f;
        if (wb_strength > 1.0f) wb_strength = 1.0f;
        gain_r = 1.0f + wb_strength * (full_r - 1.0f);
        gain_g = 1.0f + wb_strength * (full_g - 1.0f);
        gain_b = 1.0f + wb_strength * (full_b - 1.0f);
    }

    /* --- Step 1: per-channel linear normalisation to [0, 1] --- */
    for (int c = 0; c < 3; c++) {
        float emin =  FLT_MAX;
        float emax = -FLT_MAX;
        for (int p = 0; p < N; p++) {
            float v = irr[p * 3 + c];
            if (c == 0) v *= gain_r;
            if (c == 1) v *= gain_g;
            if (c == 2) v *= gain_b;
            if (v < emin) emin = v;
            if (v > emax) emax = v;
        }
        float range = (emax > emin) ? (emax - emin) : 1.0f;
        for (int p = 0; p < N; p++) {
            float v = irr[p * 3 + c];
            if (c == 0) v *= gain_r;
            if (c == 1) v *= gain_g;
            if (c == 2) v *= gain_b;
            E[p * 3 + c] = (v - emin) / range;
        }
    }

    /* --- Step 2: gamma correction --- */
    for (int i = 0; i < N * 3; i++)
        E[i] = powf(E[i], gamma);

    /* --- Step 3: luminance --- */
    for (int p = 0; p < N; p++)
        L[p] = LUM_R * E[p*3+0] + LUM_G * E[p*3+1] + LUM_B * E[p*3+2];

    /* --- Step 4: geometric mean of luminance --- */
    double log_sum = 0.0;
    for (int p = 0; p < N; p++)
        log_sum += log((double)L[p] + 1e-6);
    double L_avg = exp(log_sum / N);

    /* --- Step 5: key-scaled luminance T --- */
    float T_max = 0.0f;
    float scale = (float)(alpha / L_avg);
    for (int p = 0; p < N; p++) {
        T[p] = scale * L[p];
        if (T[p] > T_max) T_max = T[p];
    }

    /* --- Step 6: Reinhard tone curve --- */
    float Tmax2 = T_max * T_max;
    if (Tmax2 < 1e-10f) Tmax2 = 1e-10f;
    for (int p = 0; p < N; p++)
        Ltone[p] = T[p] * (1.0f + T[p] / Tmax2) / (1.0f + T[p]);

    /* --- Step 7: scale each channel and clip to [0,1] --- */
    for (int p = 0; p < N; p++) {
        float M = (L[p] > 1e-7f) ? (Ltone[p] / L[p]) : 0.0f;
        for (int c = 0; c < 3; c++) {
            float v = E[p * 3 + c] * M;
            if (v < 0.0f) v = 0.0f;
            if (v > 1.0f) v = 1.0f;
            E[p * 3 + c] = v;
        }
    }

    /* Optional shadow anchoring when WB is active.
     * WB can lift the shadow floor slightly; anchor black to the 1st-percentile
     * luminance and remap back to [0,1]. */
    if (wb_strength > 0.0f) {
        float *Lsorted = (float *)malloc((size_t)N * sizeof(float));
        if (Lsorted) {
            for (int p = 0; p < N; p++) {
                float r = E[p * 3 + 0];
                float g = E[p * 3 + 1];
                float b = E[p * 3 + 2];
                Lsorted[p] = LUM_R * r + LUM_G * g + LUM_B * b;
            }
            qsort(Lsorted, (size_t)N, sizeof(float), cmp_float);

            int idx = (int)(0.01f * (N - 1));
            float black = Lsorted[idx];
            if (black < 0.0f) black = 0.0f;
            if (black > 0.08f) black = 0.08f;

            if (black > 1e-6f) {
                float inv = 1.0f / (1.0f - black);
                for (int i = 0; i < N * 3; i++) {
                    float v = (E[i] - black) * inv;
                    if (v < 0.0f) v = 0.0f;
                    if (v > 1.0f) v = 1.0f;
                    E[i] = v;
                }
            }
            free(Lsorted);
        }
    }

    /* --- Step 8: quantise to uint8 --- */
    for (int i = 0; i < N * 3; i++)
        out[i] = (uint8_t)(E[i] * 255.0f + 0.5f);

    free(E);
    free(L);
    free(T);
    free(Ltone);
}
