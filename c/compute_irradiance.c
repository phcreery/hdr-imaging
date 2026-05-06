#include "compute_irradiance.h"
#include <math.h>

void compute_irradiance(uint8_t **images, int num_images, int W, int H,
                        const float *B, const float *crf, const float *w,
                        float *irr_out)
{
    int N = W * H;
    for (int p = 0; p < N; p++) {
        for (int c = 0; c < 3; c++) {
            const float *g = crf + c * 256;
            float num_sum = 0.0f;
            float den_sum = 0.0f;
            for (int j = 0; j < num_images; j++) {
                int   z  = images[j][p * 3 + c];
                float wz = w[z];
                num_sum += wz * (g[z] - B[j]);
                den_sum += wz;
            }
            irr_out[p * 3 + c] = expf(num_sum / (den_sum + 1e-6f));
        }
    }
}
