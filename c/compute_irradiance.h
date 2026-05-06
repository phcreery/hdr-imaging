#pragma once
#include <stdint.h>

/*
 * Compute the HDR radiance (irradiance) map from a multi-exposure image stack.
 *
 * Ports compute_irradiance.py.
 *
 * For each pixel p and colour channel c:
 *
 *   log_E[p,c] = sum_j( w[Z_j] * (g[Z_j] - B[j]) )
 *              / sum_j( w[Z_j] )
 *   E[p,c]     = exp( log_E[p,c] )
 *
 * where Z_j = images[j][p*3+c].
 *
 * Inputs:
 *   images[num_images]   W*H*3 uint8_t RGB buffers, descending exposure order
 *   B[num_images]        log exposure times
 *   crf[768]             camera response: crf[c*256 + z] = g_c(z)
 *   w[256]               weighting function
 *
 * Output:
 *   irr_out[W*H*3]       HDR irradiance map, float, row-major interleaved RGB
 *                        (pre-allocated by caller)
 */
void compute_irradiance(uint8_t **images, int num_images, int W, int H,
                        const float *B, const float *crf, const float *w,
                        float *irr_out);
