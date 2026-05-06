#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "load_images.h"
#include "hdr_debevec.h"
#include "compute_irradiance.h"
#include "tonemap.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s <img1> <img2> [<imgN> ...] <output.png>\n"
        "\n"
        "  Input images must be named {num}_{den}.<ext> so that the exposure\n"
        "  time can be parsed from the filename, e.g. 1_4.png = 0.25 s.\n"
        "\n"
        "Options:\n"
        "  --lambda <f>   Smoothness weight for CRF estimation  (default: 50.0)\n"
        "  --num-px <n>   Pixels sampled for CRF estimation     (default: 150)\n"
        "  --alpha  <f>   Reinhard key/exposure parameter       (default: 0.35)\n"
        "  --gamma  <f>   Display gamma exponent                (default: 1/2.2)\n"
        "  --wb     <f>   Gray-world WB blend [0..1]            (default: 0.0)\n"
        "\n"
        "Example:\n"
        "  %s ../images/1_4.png ../images/32_1.png output.png\n",
        prog, prog);
}

int main(int argc, char **argv)
{
    if (argc < 3) { usage(argv[0]); return 1; }

    /* Defaults (match run_hdr_image.py) */
    float lambda_  = 50.0f;
    int   num_px   = 150;
    float alpha    = 0.35f;
    float gamma_   = 1.0f / 2.2f;
    float wb_strength = 0.0f;

    /* Collect positional args and parse optional flags */
    const char **pos = (const char **)malloc((size_t)argc * sizeof(char *));
    int npos = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--lambda") == 0 && i + 1 < argc) {
            lambda_ = (float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--num-px") == 0 && i + 1 < argc) {
            num_px = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--alpha") == 0 && i + 1 < argc) {
            alpha = (float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--gamma") == 0 && i + 1 < argc) {
            gamma_ = (float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--wb") == 0 && i + 1 < argc) {
            wb_strength = (float)atof(argv[++i]);
            if (wb_strength < 0.0f) wb_strength = 0.0f;
            if (wb_strength > 1.0f) wb_strength = 1.0f;
        } else if (argv[i][0] == '-' && argv[i][1] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(argv[0]);
            free(pos);
            return 1;
        } else {
            pos[npos++] = argv[i];
        }
    }

    if (npos < 3) {
        fprintf(stderr, "Error: need at least 2 input images and 1 output path.\n\n");
        usage(argv[0]);
        free(pos);
        return 1;
    }

    const char  *output_path  = pos[npos - 1];
    const char **input_paths  = pos;
    int          num_images   = npos - 1;

    /* ------------------------------------------------------------------ */
    /* 1. Load images                                                       */
    /* ------------------------------------------------------------------ */
    printf("=== Loading %d image(s) ===\n", num_images);
    int W, H;
    float *B = (float *)malloc((size_t)num_images * sizeof(float));
    uint8_t **images = load_images(input_paths, num_images, &W, &H, B);
    if (!images) { free(B); free(pos); return 1; }

    printf("  Dimensions: %d x %d\n", W, H);
    for (int i = 0; i < num_images; i++)
        printf("  [%d] log_exp=%.4f  (t=%.5f s)\n", i, B[i], expf(B[i]));

    /* ------------------------------------------------------------------ */
    /* 2. Estimate CRF                                                      */
    /* ------------------------------------------------------------------ */
    float *crf = (float *)malloc(3 * 256 * sizeof(float));
    float *w   = (float *)malloc(256       * sizeof(float));
    printf("\n=== Estimating Camera Response Function ===\n");
    printf("  lambda=%.1f  num_px=%d\n", lambda_, num_px);
    hdr_debevec(images, num_images, W, H, B, lambda_, num_px, crf, w);

    /* Sanity check: anchor constraint g[128] should be ~0 for all channels */
    printf("  g_R[128]=%.6f  g_G[128]=%.6f  g_B[128]=%.6f  (all should be ~0)\n",
           crf[0*256+128], crf[1*256+128], crf[2*256+128]);

    /* ------------------------------------------------------------------ */
    /* 3. Compute HDR irradiance map                                        */
    /* ------------------------------------------------------------------ */
    printf("\n=== Computing irradiance map ===\n");
    float *irr = (float *)malloc((size_t)W * H * 3 * sizeof(float));
    compute_irradiance(images, num_images, W, H, B, crf, w, irr);

    /* Quick stats */
    float irr_min =  1e30f, irr_max = -1e30f;
    for (int i = 0; i < W * H * 3; i++) {
        if (irr[i] < irr_min) irr_min = irr[i];
        if (irr[i] > irr_max) irr_max = irr[i];
    }
    printf("  Irradiance range: [%.4e, %.4e]\n", irr_min, irr_max);

    /* ------------------------------------------------------------------ */
    /* 4. Tone mapping                                                      */
    /* ------------------------------------------------------------------ */
    printf("\n=== Reinhard tone mapping ===\n");
    printf("  alpha=%.4f  gamma=%.4f  wb=%.3f\n", alpha, gamma_, wb_strength);

    uint8_t *ldr = (uint8_t *)malloc((size_t)W * H * 3);
    reinhard_tonemap(irr, W, H, alpha, gamma_, wb_strength, ldr);

    /* ------------------------------------------------------------------ */
    /* 5. Save output                                                       */
    /* ------------------------------------------------------------------ */
    printf("\n=== Saving output -> %s ===\n", output_path);
    if (!stbi_write_png(output_path, W, H, 3, ldr, W * 3)) {
        fprintf(stderr, "Failed to write '%s'\n", output_path);
        free(ldr); free(irr); free(crf); free(w);
        free(B); free_images(images, num_images); free(pos);
        return 1;
    }
    printf("Done.\n");

    free(ldr);
    free(irr);
    free(crf);
    free(w);
    free(B);
    free_images(images, num_images);
    free(pos);
    return 0;
}
