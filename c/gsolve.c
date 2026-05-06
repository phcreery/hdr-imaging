#include "gsolve.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#ifdef USE_NML_LS
#include "neat-matrix-library/nml.h"
#endif

/*
 * In-place Householder QR least squares — double precision.
 *
 * A (m x n, row-major) and b (m) are modified in place.
 * After return, b[0..n-1] contains the LS solution x.
 *
 * Using double throughout because the gsolve matrix is sparse and
 * poorly-conditioned; float32 produces visibly wrong CRF curves.
 */
static void qr_lstsq(double *A, double *b, int m, int n)
{
    double *v = (double *)malloc((size_t)m * sizeof(double));

    for (int k = 0; k < n; k++) {
        int len = m - k;

        /* Squared norm of A[k:m, k] */
        double norm2 = 0.0;
        for (int i = 0; i < len; i++) {
            double a = A[(k + i) * n + k];
            norm2 += a * a;
        }
        if (norm2 < 1e-48) continue;
        double norm = sqrt(norm2);

        /* Build Householder vector v = A[k:m,k] + sign(A[k,k])*||...||*e1 */
        for (int i = 0; i < len; i++) v[i] = A[(k + i) * n + k];
        v[0] += (v[0] >= 0.0) ? norm : -norm;

        /* Squared norm of v */
        double vnorm2 = 0.0;
        for (int i = 0; i < len; i++) vnorm2 += v[i] * v[i];
        if (vnorm2 < 1e-48) continue;
        double inv2 = 2.0 / vnorm2;

        /* Apply H = I - 2*v*v^T/||v||^2 to A[k:m, k:n] */
        for (int j = k; j < n; j++) {
            double dot = 0.0;
            for (int i = 0; i < len; i++) dot += v[i] * A[(k + i) * n + j];
            dot *= inv2;
            for (int i = 0; i < len; i++) A[(k + i) * n + j] -= dot * v[i];
        }

        /* Apply same reflector to b */
        {
            double dot = 0.0;
            for (int i = 0; i < len; i++) dot += v[i] * b[k + i];
            dot *= inv2;
            for (int i = 0; i < len; i++) b[k + i] -= dot * v[i];
        }
    }

    free(v);

    /* Back-substitution: solve upper-triangular R x = c (c stored in b[0..n-1]) */
    for (int k = n - 1; k >= 0; k--) {
        double s = b[k];
        for (int j = k + 1; j < n; j++) s -= A[k * n + j] * b[j];
        double diag = A[k * n + k];
        b[k] = (diag * diag > 1e-40) ? (s / diag) : 0.0;
    }
}

static double now_ms(void)
{
    return 1000.0 * (double)clock() / (double)CLOCKS_PER_SEC;
}

#ifdef USE_NML_LS
/*
 * Experimental LS solver path using neat-matrix-library.
 *
 * Solves min ||Ax-b||_2 via normal equations:
 *   (A^T A) x = A^T b
 *
 * This is intended for experimentation only; normal equations can be
 * less numerically stable than QR on ill-conditioned systems.
 */
static int nml_lstsq_normal_eq(const double *A, const double *b, int m, int n, double *x_out)
{
    int ok = 0;
    nml_mat *A_m = NULL, *b_m = NULL, *At = NULL, *AtA = NULL, *Atb = NULL;
    nml_mat_lup *lup = NULL;
    nml_mat *x = NULL;

    A_m = nml_mat_new((unsigned int)m, (unsigned int)n);
    b_m = nml_mat_new((unsigned int)m, 1u);
    if (!A_m || !b_m) goto cleanup;

    for (int i = 0; i < m; i++) {
        b_m->data[i][0] = b[i];
        for (int j = 0; j < n; j++) {
            A_m->data[i][j] = A[(size_t)i * (size_t)n + (size_t)j];
        }
    }

    At  = nml_mat_transp(A_m);
    AtA = nml_mat_dot(At, A_m);
    Atb = nml_mat_dot(At, b_m);
    if (!At || !AtA || !Atb) goto cleanup;

    /*
     * Regularize normal equations to avoid singular/near-singular AtA on
     * underconstrained or ill-conditioned CRF systems (e.g. very few images).
     */
    {
        double trace_abs = 0.0;
        for (int i = 0; i < n; i++) {
            trace_abs += fabs(AtA->data[i][i]);
        }
        double avg_diag = (n > 0) ? (trace_abs / (double)n) : 1.0;
        if (avg_diag < 1e-12) avg_diag = 1.0;

        /* Small Tikhonov term: AtA += ridge * I */
        const double ridge = avg_diag * 1e-8;
        for (int i = 0; i < n; i++) {
            AtA->data[i][i] += ridge;
        }
    }

    lup = nml_mat_lup_solve(AtA);
    if (!lup) goto cleanup;

    x = nml_ls_solve(lup, Atb);
    if (!x) goto cleanup;

    for (int i = 0; i < n; i++) {
        x_out[i] = x->data[i][0];
    }
    ok = 1;

cleanup:
    if (x) nml_mat_free(x);
    if (lup) nml_mat_lup_free(lup);
    if (Atb) nml_mat_free(Atb);
    if (AtA) nml_mat_free(AtA);
    if (At) nml_mat_free(At);
    if (b_m) nml_mat_free(b_m);
    if (A_m) nml_mat_free(A_m);
    return ok;
}
#endif

/*
 * Build and solve the Debevec / Malik linear system for one colour channel.
 *
 * System A x = b  (overdetermined, solved in LS sense):
 *   Unknowns  x  = [ g(0)..g(255) | lE(0)..lE(P-1) ]   size: 256 + P
 *   Rows:
 *     Data equations  (P * N rows): w[Z_ij] * (g[Z_ij] - lE_i) = w[Z_ij]*B_j
 *     Anchor          (1 row):       g[128] = 0
 *     Smoothness      (254 rows):    lambda*w[z+1]*(g[z]-2g[z+1]+g[z+2]) = 0
 *     (one zero-padding row to match the Python +256 allocation)
 */
void gsolve(const int *Z, int num_px, int num_im,
            const float *B, float lambda_,
            const float *w,
            float *g_out, float *lE_out)
{
    const int Zn     = 256;
    const int m_rows = num_px * num_im + Zn;
    const int n_cols = Zn + num_px;

    double *A = (double *)calloc((size_t)m_rows * n_cols, sizeof(double));
    double *b = (double *)calloc((size_t)m_rows,           sizeof(double));
    if (!A || !b) {
        fprintf(stderr, "gsolve: out of memory (%d x %d matrix)\n", m_rows, n_cols);
        free(A); free(b);
        return;
    }

    int k = 0;

    /* Data fitting equations */
    for (int i = 0; i < num_px; i++) {
        for (int j = 0; j < num_im; j++) {
            int    z   = Z[i * num_im + j];
            double wij = (double)w[z];
            A[(size_t)k * n_cols + z]      =  wij;
            A[(size_t)k * n_cols + Zn + i] = -wij;
            b[k] = wij * (double)B[j];
            k++;
        }
    }

    /* Anchor: g[128] = 0 */
    A[(size_t)k * n_cols + Zn / 2] = 1.0;
    /* b[k] = 0 (already calloc'd) */
    k++;

    /* Smoothness: lambda * w[z+1] * (g[z] - 2*g[z+1] + g[z+2]) = 0 */
    for (int z = 0; z < Zn - 2; z++) {
        double lw = (double)lambda_ * (double)w[z + 1];
        A[(size_t)k * n_cols + z]     =  lw;
        A[(size_t)k * n_cols + z + 1] = -2.0 * lw;
        A[(size_t)k * n_cols + z + 2] =  lw;
        /* b[k] = 0 */
        k++;
    }
    /* Row k = num_px*num_im + 255 is the spare zero-padding row — left as zero */

#ifdef USE_NML_LS
    {
        double t0_ms = now_ms();
        double *x = (double *)calloc((size_t)n_cols, sizeof(double));
        if (!x) {
            fprintf(stderr, "gsolve: out of memory (nml solution vector)\n");
            qr_lstsq(A, b, m_rows, n_cols);
            fprintf(stderr, "gsolve: solver=qr_lstsq elapsed=%.3f ms (m=%d n=%d, oom fallback)\n",
                    now_ms() - t0_ms, m_rows, n_cols);
            for (int z = 0; z < Zn; z++) g_out[z] = (float)b[z];
            for (int i = 0; i < num_px; i++) lE_out[i] = (float)b[Zn + i];
            free(A);
            free(b);
            return;
        }

        if (nml_lstsq_normal_eq(A, b, m_rows, n_cols, x)) {
            fprintf(stderr, "gsolve: solver=nml_normal_eq elapsed=%.3f ms (m=%d n=%d)\n",
                    now_ms() - t0_ms, m_rows, n_cols);
            for (int z = 0; z < Zn; z++) g_out[z] = (float)x[z];
            for (int i = 0; i < num_px; i++) lE_out[i] = (float)x[Zn + i];
        } else {
            fprintf(stderr, "gsolve: nml normal-equation solve failed, falling back to qr_lstsq\n");
            qr_lstsq(A, b, m_rows, n_cols);
            fprintf(stderr, "gsolve: solver=qr_lstsq elapsed=%.3f ms (m=%d n=%d, nml fallback)\n",
                    now_ms() - t0_ms, m_rows, n_cols);
            for (int z = 0; z < Zn; z++) g_out[z] = (float)b[z];
            for (int i = 0; i < num_px; i++) lE_out[i] = (float)b[Zn + i];
        }

        free(x);
    }
#else
    {
    double t0_ms = now_ms();
    qr_lstsq(A, b, m_rows, n_cols);
    fprintf(stderr, "gsolve: solver=qr_lstsq elapsed=%.3f ms (m=%d n=%d)\n",
            now_ms() - t0_ms, m_rows, n_cols);
    }

    /* Extract solution */
    for (int z = 0; z < Zn;    z++) g_out[z]  = b[z];
    for (int i = 0; i < num_px; i++) lE_out[i] = b[Zn + i];
#endif

    free(A);
    free(b);
}
