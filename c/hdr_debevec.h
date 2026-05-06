#pragma once
#include <stdint.h>

/*
 * Estimate the Camera Response Function (CRF) using the Debevec & Malik method.
 *
 * Ports hdr_debevec.py: samples random pixels, builds the Z matrix, then
 * calls gsolve() once per colour channel.
 *
 * Inputs:
 *   images[num_images]   array of pointers to W*H*3 uint8_t RGB buffers
 *                        sorted by exposure descending
 *   B[num_images]        log exposure times (descending)
 *   lambda_              smoothness weight (default: 50.0)
 *   num_px               pixels to sample for CRF estimation (default: 150)
 *
 * Outputs:
 *   crf_out[768]         log response curves: crf_out[c*256 + z] = g_c(z)
 *                        for c in {0=R, 1=G, 2=B}
 *   w_out[256]           triangular weighting function w[z] = min(z, 255-z)
 */
void hdr_debevec(uint8_t **images, int num_images, int W, int H,
                 const float *B, float lambda_, int num_px,
                 float *crf_out, float *w_out);
