#include "load_images.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

/* Parse exposure time (seconds) from filename like "32_1.png" or "1_4.jpg".
   Looks at the basename, strips extension, splits on the LAST underscore. */
static float parse_exposure(const char *path)
{
    /* Find basename (last / or \ separator) */
    const char *base = path;
    for (const char *p = path; *p; p++)
        if (*p == '/' || *p == '\\') base = p + 1;

    char buf[256];
    strncpy(buf, base, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';

    /* Strip extension */
    char *dot = strrchr(buf, '.');
    if (dot) *dot = '\0';

    /* Split numerator/denominator on last '_' */
    char *under = strrchr(buf, '_');
    if (!under) return 1.0f;
    *under = '\0';

    int num = atoi(buf);
    int den = atoi(under + 1);
    return (den != 0) ? (float)num / (float)den : 1.0f;
}

typedef struct { uint8_t *pixels; float exposure; } ImgEntry;

static int cmp_desc_exposure(const void *a, const void *b)
{
    float ea = ((const ImgEntry *)a)->exposure;
    float eb = ((const ImgEntry *)b)->exposure;
    return (eb > ea) - (eb < ea);   /* descending */
}

uint8_t **load_images(const char **paths, int num_paths,
                      int *w_out, int *h_out,
                      float *B_out)
{
    ImgEntry *entries = (ImgEntry *)malloc((size_t)num_paths * sizeof(ImgEntry));
    if (!entries) return NULL;

    int W = 0, H = 0;
    for (int i = 0; i < num_paths; i++) {
        int w, h, ch;
        /* Force 3 channels (RGB) regardless of source */
        uint8_t *px = stbi_load(paths[i], &w, &h, &ch, 3);
        if (!px) {
            fprintf(stderr, "load_images: failed to load '%s': %s\n",
                    paths[i], stbi_failure_reason());
            for (int j = 0; j < i; j++) free(entries[j].pixels);
            free(entries);
            return NULL;
        }
        if (W == 0) { W = w; H = h; }
        else if (w != W || h != H) {
            fprintf(stderr, "load_images: size mismatch in '%s' (%dx%d vs %dx%d)\n",
                    paths[i], w, h, W, H);
            free(px);
            for (int j = 0; j < i; j++) free(entries[j].pixels);
            free(entries);
            return NULL;
        }
        entries[i].pixels   = px;
        entries[i].exposure = parse_exposure(paths[i]);
        printf("  Loaded '%s'  exposure=%.4f s\n", paths[i], entries[i].exposure);
    }

    /* Sort descending by exposure time (matches Python pipeline) */
    qsort(entries, (size_t)num_paths, sizeof(ImgEntry), cmp_desc_exposure);

    uint8_t **result = (uint8_t **)malloc((size_t)num_paths * sizeof(uint8_t *));
    if (!result) {
        for (int i = 0; i < num_paths; i++) free(entries[i].pixels);
        free(entries);
        return NULL;
    }

    for (int i = 0; i < num_paths; i++) {
        result[i] = entries[i].pixels;
        B_out[i]  = logf(entries[i].exposure);
    }

    *w_out = W;
    *h_out = H;
    free(entries);
    return result;
}

void free_images(uint8_t **images, int num_paths)
{
    for (int i = 0; i < num_paths; i++) free(images[i]);
    free(images);
}
