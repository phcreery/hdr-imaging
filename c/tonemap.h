#pragma once
#include <stdint.h>

/*
 * Reinhard global tone mapping.
 *
 * Ports reinhard_tonemap() from tonemap.py.
 *
 * Pipeline:
 *   1. Per-channel linear normalisation to [0, 1]
 *   2. Gamma correction:  E = E ^ gamma
 *   3. Luminance:         L = 0.299*R + 0.587*G + 0.114*B  (BT.601)
 *   4. Geometric mean:    L_avg = exp( mean(log(L)) )
 *   5. Key value:         T = (alpha / L_avg) * L
 *   6. Tone curve:        L_tone = T * (1 + T/Tmax^2) / (1 + T)
 *   7. Scaling:           M = L_tone / L
 *   8. Apply + clip:      out[c] = clip(E[c] * M, 0, 1)  -> uint8
 *
 * Inputs:
 *   irr[W*H*3]   HDR irradiance map, float, row-major interleaved RGB
 *   alpha        key/exposure parameter (default: 0.35)
 *   gamma        display gamma exponent (default: 1/2.2 ≈ 0.4545)
 *   wb_strength  gray-world white-balance blend [0..1]
 *                0.0 = exact Python behavior (no WB)
 *                1.0 = full gray-world correction
 *
 * Output:
 *   out[W*H*3]   LDR image, uint8_t [0-255], row-major interleaved RGB
 *                (pre-allocated by caller)
 */
void reinhard_tonemap(const float *irr, int W, int H,
                      float alpha, float gamma, float wb_strength,
                      uint8_t *out);
