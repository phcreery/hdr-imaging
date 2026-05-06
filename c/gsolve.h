#pragma once

/*
 * Solve for the imaging system response function (CRF).
 *
 * Ports gsolve.py: constructs an overdetermined linear system A x = b and
 * solves it via Householder QR least squares.
 *
 * Buffer layout (all flat, GLSL-friendly):
 *   Z[num_px * num_im]   pixel intensities: Z[i*num_im + j] = intensity of
 *                        pixel i in image j  (int, 0-255)
 *   B[num_im]            log exposure times (descending order)
 *   w[256]               per-intensity weighting function
 *   g_out[256]           OUTPUT: log camera response g(z), with g[128]=0
 *   lE_out[num_px]       OUTPUT: log film irradiance at each sampled pixel
 */
void gsolve(const int *Z, int num_px, int num_im,
            const float *B, float lambda_,
            const float *w,
            float *g_out, float *lE_out);
