#include "hdr_debevec.h"
#include "gsolve.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/*
 * Sample num_px unique indices from [0, total-1] without replacement.
 * Uses a simple LCG with a fixed seed for reproducible results across runs
 * (matches using a fixed numpy random seed on the Python side).
 */
static void sample_indices(int *out, int num_px, int total)
{
    int *pool = (int *)malloc((size_t)total * sizeof(int));
    for (int i = 0; i < total; i++)
        pool[i] = i;

    /* LCG parameters from Numerical Recipes (32-bit) */
    unsigned int rng = 42u;
    for (int i = 0; i < num_px; i++)
    {
        rng = rng * 1664525u + 1013904223u;
        int j = i + (int)((rng >> 1) % (unsigned int)(total - i));
        int tmp = pool[i];
        pool[i] = pool[j];
        pool[j] = tmp;
    }
    memcpy(out, pool, (size_t)num_px * sizeof(int));
    free(pool);
}

void hdr_debevec(uint8_t **images, int num_images, int W, int H,
                 const float *B, float lambda_, int num_px,
                 float *crf_out, float *w_out)
{
    /* w[z] = min(z, 255-z)  — triangular weighting, mirrors Python:
         np.concatenate((np.arange(128)-0, 255-np.arange(128,256)))
       which gives [0,1,...,127, 127,126,...,0]                       */
    for (int z = 0; z < 256; z++)
        w_out[z] = (float)((z < 128) ? z : (255 - z));

    /* Sample pixel positions */
    int *px_idx = (int *)malloc((size_t)num_px * sizeof(int));
    sample_indices(px_idx, num_px, W * H);

    /* Z[i*num_images + j] = pixel i, image j (per-channel) */
    int *Z = (int *)malloc((size_t)num_px * num_images * sizeof(int));
    float *lE = (float *)malloc((size_t)num_px * sizeof(float));

    printf("  Building CRF for 3 channels (%d pixels x %d images)...\n",
           num_px, num_images);

    for (int ch = 0; ch < 3; ch++)
    {
        /* Fill Z for this channel */
        for (int i = 0; i < num_px; i++)
        {
            int px = px_idx[i];
            for (int j = 0; j < num_images; j++)
                Z[i * num_images + j] = images[j][px * 3 + ch];
        }

        float *g = crf_out + ch * 256;
        gsolve(Z, num_px, num_images, B, lambda_, w_out, g, lE);
        printf("    Channel %d: g[0]=%.4f  g[128]=%.6f  g[255]=%.4f\n",
               ch, g[0], g[128], g[255]);
    }

    free(px_idx);
    free(Z);
    free(lE);
}
