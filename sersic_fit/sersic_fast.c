/*
 * sersic_fast.c — Fast 2D Sérsic model evaluation with analytical Jacobian.
 *
 * Provides three entry points:
 *   sersic_2d_eval()          — model image
 *   sersic_2d_residual()      — data - model
 *   sersic_2d_jacobian()      — d(residual)/d(params), shape (npix, 8)
 *   sersic_2d_cost_and_grad() — scalar cost + gradient (for L-BFGS-B)
 *
 * Build:
 *   gcc -O3 -march=native -ffast-math -fopenmp -shared -fPIC \
 *       -o libsersic_fast.so sersic_fit/sersic_fast.c -lm
 *
 * Parameters (8): Ie, re, n, xc, yc, ellip, theta, bg
 */

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* bn ≈ 1.9992*n - 0.3271 (Ciotti & Bertin 1999) */
static inline double sersic_bn(double n)
{
    return (n > 0.36) ? (1.9992 * n - 0.3271) : 0.01;
}

/* dbn/dn = 1.9992 for n > 0.36, else 0 */
static inline double sersic_dbn_dn(double n)
{
    return (n > 0.36) ? 1.9992 : 0.0;
}

static inline double clamp(double v, double lo, double hi)
{
    return (v < lo) ? lo : ((v > hi) ? hi : v);
}

/*
 * Evaluate 2D Sérsic model at npix pixel positions.
 *
 * params[8] = {Ie, re, n, xc, yc, ellip, theta, bg}
 * x[npix], y[npix] = pixel coordinates (flattened)
 * out[npix] = model values
 */
void sersic_2d_eval(const double *params,
                    const double *x, const double *y,
                    double *out, int npix)
{
    double Ie    = params[0];
    double re    = params[1];
    double n     = params[2];
    double xc    = params[3];
    double yc    = params[4];
    double ellip = params[5];
    double theta = params[6];
    double bg    = params[7];

    if (re < 0.1)  re = 0.1;
    if (n  < 0.1)  n  = 0.1;

    double q = 1.0 - ellip;
    if (q < 0.05) q = 0.05;

    double bn      = sersic_bn(n);
    double inv_n   = 1.0 / n;
    double inv_re  = 1.0 / re;
    double cos_t   = cos(theta);
    double sin_t   = sin(theta);
    double inv_q2  = 1.0 / (q * q);

    #pragma omp parallel for schedule(static) if(npix > 1000)
    for (int i = 0; i < npix; i++) {
        double dx    = x[i] - xc;
        double dy    = y[i] - yc;
        double x_rot = dx * cos_t + dy * sin_t;
        double y_rot = -dx * sin_t + dy * cos_t;

        double r = sqrt(x_rot * x_rot + y_rot * y_rot * inv_q2);

        /* Match Python: 0**(1/n) = 0, so rn=0 at center */
        double rn = (r < 1e-10) ? 0.0 : pow(r * inv_re, inv_n);
        double E  = bn * (rn - 1.0);
        E = clamp(E, -50.0, 50.0);

        out[i] = Ie * exp(-E) + bg;
    }
}

/*
 * Compute residuals: residual[i] = data[i] - model[i].
 */
void sersic_2d_residual(const double *params,
                        const double *data,
                        const double *x, const double *y,
                        double *residual, int npix)
{
    double Ie    = params[0];
    double re    = params[1];
    double n     = params[2];
    double xc    = params[3];
    double yc    = params[4];
    double ellip = params[5];
    double theta = params[6];
    double bg    = params[7];

    if (re < 0.1)  re = 0.1;
    if (n  < 0.1)  n  = 0.1;

    double q = 1.0 - ellip;
    if (q < 0.05) q = 0.05;

    double bn      = sersic_bn(n);
    double inv_n   = 1.0 / n;
    double inv_re  = 1.0 / re;
    double cos_t   = cos(theta);
    double sin_t   = sin(theta);
    double inv_q2  = 1.0 / (q * q);

    #pragma omp parallel for schedule(static) if(npix > 1000)
    for (int i = 0; i < npix; i++) {
        double dx    = x[i] - xc;
        double dy    = y[i] - yc;
        double x_rot = dx * cos_t + dy * sin_t;
        double y_rot = -dx * sin_t + dy * cos_t;

        double r = sqrt(x_rot * x_rot + y_rot * y_rot * inv_q2);

        /* Match Python: 0**(1/n) = 0, so rn=0 at center */
        double rn = (r < 1e-10) ? 0.0 : pow(r * inv_re, inv_n);
        double E  = bn * (rn - 1.0);
        E = clamp(E, -50.0, 50.0);

        double model = Ie * exp(-E) + bg;
        residual[i] = data[i] - model;
    }
}

/*
 * Analytical Jacobian of residuals w.r.t. parameters.
 *
 * jac[i*8 + j] = d(residual_i) / d(param_j)
 *              = -d(model_i) / d(param_j)
 *
 * Parameters: {Ie, re, n, xc, yc, ellip, theta, bg}
 *              0    1   2  3   4   5      6      7
 *
 * Model: f = Ie * exp(-E) + bg
 *   E = bn * ((r/re)^(1/n) - 1)
 *   r = sqrt(x_rot^2 + (y_rot/q)^2),  q = 1-ellip
 *
 * Derivatives (df/d*):
 *   df/dIe    = exp(-E) = profile / Ie
 *   df/dre    = profile * bn * rn / (n * re)
 *   df/dn     = profile * [-dbn_dn*(rn-1) + bn*rn*ln(r/re)/(n*n)]
 *   df/dxc    = profile * (-dE_dr) * dr/dxc
 *   df/dyc    = profile * (-dE_dr) * dr/dyc
 *   df/dellip = profile * (-dE_dr) * dr/dellip
 *   df/dtheta = profile * (-dE_dr) * dr/dtheta
 *   df/dbg    = 1
 *
 *   dE_dr = bn * rn / (n * r)
 *
 * Ellipse radius chain rule:
 *   dr/dxc    = -(x_rot * cos_t + y_rot * sin_t * inv_q2) / r  -- note sign
 *              actually dr/dxc = d/dxc sqrt(x_rot^2 + (y_rot/q)^2)
 *              x_rot = (x-xc)*cos + (y-yc)*sin => dx_rot/dxc = -cos_t
 *              y_rot = -(x-xc)*sin + (y-yc)*cos => dy_rot/dxc = sin_t
 *   dr/dxc = (x_rot*(-cos_t) + y_rot*sin_t*inv_q2) / r
 *   dr/dyc = (x_rot*(-sin_t) + y_rot*(-cos_t)*inv_q2) / r
 *   dr/dellip = y_rot^2 / (q^3 * r)
 *   dr/dtheta: dx_rot/dtheta = -x_rot_perp, dy_rot/dtheta = x_rot_perp_y
 *              x_rot = dx*cos + dy*sin  =>  d/dtheta = -dx*sin + dy*cos = y_rot_orig
 *              Wait, need to be careful:
 *              dx_rot/dtheta = dx*(-sin_t) + dy*cos_t = y_rot (using original def)
 *              dy_rot/dtheta = -dx*cos_t + dy*(-sin_t) = -(dx*cos_t + dy*sin_t) = -x_rot
 *              Wait: y_rot = -dx*sin_t + dy*cos_t, so:
 *              dx_rot/dtheta = -dx*sin_t + dy*cos_t = y_rot
 *              dy_rot/dtheta = -(dx*cos_t + dy*sin_t) = -x_rot
 *   dr/dtheta = (x_rot * y_rot + y_rot * (-x_rot) * inv_q2) / r
 *             = x_rot * y_rot * (1 - inv_q2) / r
 *
 * Jacobian is -df/d* (since residual = data - f).
 */
void sersic_2d_jacobian(const double *params,
                        const double *x, const double *y,
                        double *jac, int npix)
{
    double Ie    = params[0];
    double re    = params[1];
    double n     = params[2];
    double xc    = params[3];
    double yc    = params[4];
    double ellip = params[5];
    double theta = params[6];
    /* bg = params[7]; not needed for derivatives except df/dbg = 1 */

    if (re < 0.1)  re = 0.1;
    if (n  < 0.1)  n  = 0.1;

    double q = 1.0 - ellip;
    if (q < 0.05) q = 0.05;

    double bn       = sersic_bn(n);
    double dbn_dn   = sersic_dbn_dn(n);
    double inv_n    = 1.0 / n;
    double inv_n2   = inv_n * inv_n;
    double inv_re   = 1.0 / re;
    double cos_t    = cos(theta);
    double sin_t    = sin(theta);
    double inv_q2   = 1.0 / (q * q);
    double inv_q3_r = 1.0 / (q * q * q);  /* 1/q^3, will multiply by y_rot^2/r */
    double log_inv_re = log(inv_re);       /* = -log(re) */

    #pragma omp parallel for schedule(static) if(npix > 1000)
    for (int i = 0; i < npix; i++) {
        double dx    = x[i] - xc;
        double dy    = y[i] - yc;
        double x_rot = dx * cos_t + dy * sin_t;
        double y_rot = -dx * sin_t + dy * cos_t;

        double r2 = x_rot * x_rot + y_rot * y_rot * inv_q2;
        double r  = sqrt(r2);

        int at_center = (r < 1e-10);

        /* Match Python: 0**(1/n) = 0, so rn=0 at center */
        double ratio, rn, E, E_clamped;
        if (at_center) {
            ratio = 0.0;
            rn = 0.0;
            E = -bn;  /* bn * (0 - 1) */
        } else {
            ratio = r * inv_re;
            rn = pow(ratio, inv_n);
            E = bn * (rn - 1.0);
        }
        E_clamped = clamp(E, -50.0, 50.0);

        double profile = Ie * exp(-E_clamped);
        double expnE   = exp(-E_clamped);  /* = profile / Ie */

        /* Intermediate: dE/dr = bn * rn / (n * r) */
        double dE_dr = at_center ? 0.0 : (bn * rn * inv_n / r);

        /* dr/d* (ellipse radius derivatives) */
        double dr_dxc, dr_dyc, dr_dellip, dr_dtheta;

        if (at_center) {
            dr_dxc    = 0.0;
            dr_dyc    = 0.0;
            dr_dellip = 0.0;
            dr_dtheta = 0.0;
        } else {
            double inv_r = 1.0 / r;
            dr_dxc    = (x_rot * (-cos_t) + y_rot * sin_t * inv_q2) * inv_r;
            dr_dyc    = (x_rot * (-sin_t) + y_rot * (-cos_t) * inv_q2) * inv_r;
            dr_dellip = (y_rot * y_rot * inv_q3_r) * inv_r;
            dr_dtheta = x_rot * y_rot * (1.0 - inv_q2) * inv_r;
        }

        /* log(r/re) for dn derivative; at center rn=0, use 0 for the term */
        double log_ratio = at_center ? 0.0 : log(ratio);

        /* df/d(param) */
        double df_dIe    = expnE;
        double df_dre    = profile * bn * rn * inv_n * inv_re;
        double df_dn     = profile * (-(dbn_dn * (rn - 1.0))
                           + bn * rn * log_ratio * inv_n2);
        double df_dxc    = profile * (-dE_dr) * dr_dxc;
        double df_dyc    = profile * (-dE_dr) * dr_dyc;
        double df_dellip = profile * (-dE_dr) * dr_dellip;
        double df_dtheta = profile * (-dE_dr) * dr_dtheta;
        double df_dbg    = 1.0;

        /* Jacobian = -df/d(param) because residual = data - f */
        int base = i * 8;
        jac[base + 0] = -df_dIe;
        jac[base + 1] = -df_dre;
        jac[base + 2] = -df_dn;
        jac[base + 3] = -df_dxc;
        jac[base + 4] = -df_dyc;
        jac[base + 5] = -df_dellip;
        jac[base + 6] = -df_dtheta;
        jac[base + 7] = -df_dbg;
    }
}

/*
 * Compute scalar cost (sum of squared residuals) and its gradient in one pass.
 *
 * cost[0] = sum_i (data[i] - model[i])^2
 * grad[j] = d(cost)/d(param_j) = -2 * sum_i (residual_i * df_i/d(param_j))
 *
 * This is the key function for L-BFGS-B optimization: returns only 9 doubles
 * instead of npix residuals, eliminating per-call array allocation overhead.
 */
void sersic_2d_cost_and_grad(const double *params,
                              const double *data,
                              const double *x, const double *y,
                              double *cost, double *grad, int npix)
{
    double Ie    = params[0];
    double re    = params[1];
    double n     = params[2];
    double xc    = params[3];
    double yc    = params[4];
    double ellip = params[5];
    double theta = params[6];
    /* bg = params[7] */

    if (re < 0.1)  re = 0.1;
    if (n  < 0.1)  n  = 0.1;

    double q = 1.0 - ellip;
    if (q < 0.05) q = 0.05;

    double bn       = sersic_bn(n);
    double dbn_dn   = sersic_dbn_dn(n);
    double inv_n    = 1.0 / n;
    double inv_n2   = inv_n * inv_n;
    double inv_re   = 1.0 / re;
    double cos_t    = cos(theta);
    double sin_t    = sin(theta);
    double inv_q2   = 1.0 / (q * q);
    double inv_q3   = 1.0 / (q * q * q);

    double total_cost = 0.0;
    double g[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    #pragma omp parallel for schedule(static) if(npix > 1000) \
        reduction(+:total_cost) reduction(+:g[:8])
    for (int i = 0; i < npix; i++) {
        double dx    = x[i] - xc;
        double dy    = y[i] - yc;
        double x_rot = dx * cos_t + dy * sin_t;
        double y_rot = -dx * sin_t + dy * cos_t;

        double r = sqrt(x_rot * x_rot + y_rot * y_rot * inv_q2);

        int at_center = (r < 1e-10);
        double ratio, rn, E;
        if (at_center) {
            ratio = 0.0;
            rn = 0.0;
            E = -bn;
        } else {
            ratio = r * inv_re;
            rn = pow(ratio, inv_n);
            E = bn * (rn - 1.0);
        }
        double E_clamped = clamp(E, -50.0, 50.0);

        double expnE   = exp(-E_clamped);
        double profile = Ie * expnE;
        double model   = profile + params[7]; /* + bg */
        double residual = data[i] - model;

        total_cost += residual * residual;

        /* df/d(params) — same as in sersic_2d_jacobian */
        double dE_dr = at_center ? 0.0 : (bn * rn * inv_n / r);

        double dr_dxc, dr_dyc, dr_dellip, dr_dtheta;
        if (at_center) {
            dr_dxc = 0.0; dr_dyc = 0.0; dr_dellip = 0.0; dr_dtheta = 0.0;
        } else {
            double inv_r = 1.0 / r;
            dr_dxc    = (x_rot * (-cos_t) + y_rot * sin_t * inv_q2) * inv_r;
            dr_dyc    = (x_rot * (-sin_t) + y_rot * (-cos_t) * inv_q2) * inv_r;
            dr_dellip = (y_rot * y_rot * inv_q3) * inv_r;
            dr_dtheta = x_rot * y_rot * (1.0 - inv_q2) * inv_r;
        }

        double log_ratio = at_center ? 0.0 : log(ratio);

        double neg_dE_dr = -dE_dr;  /* for chain rule */

        double df_dIe    = expnE;
        double df_dre    = profile * bn * rn * inv_n * inv_re;
        double df_dn     = profile * (-(dbn_dn * (rn - 1.0))
                           + bn * rn * log_ratio * inv_n2);
        double df_dxc    = profile * neg_dE_dr * dr_dxc;
        double df_dyc    = profile * neg_dE_dr * dr_dyc;
        double df_dellip = profile * neg_dE_dr * dr_dellip;
        double df_dtheta = profile * neg_dE_dr * dr_dtheta;
        /* df_dbg = 1.0 */

        /* grad[j] = d(cost)/d(param_j) = -2 * sum(residual * df/dparam)
         * because cost = sum((data-model)^2), d(cost)/dp = 2*sum((d-m)*(-dm/dp))
         *                                                = -2*sum(residual * df/dp) */
        double m2r = -2.0 * residual;
        g[0] += m2r * df_dIe;
        g[1] += m2r * df_dre;
        g[2] += m2r * df_dn;
        g[3] += m2r * df_dxc;
        g[4] += m2r * df_dyc;
        g[5] += m2r * df_dellip;
        g[6] += m2r * df_dtheta;
        g[7] += m2r * 1.0; /* df_dbg */
    }

    *cost = total_cost;
    for (int j = 0; j < 8; j++)
        grad[j] = g[j];
}

/*
 * Compute J^T*J (8x8) and J^T*r (8x1) in a single pass over pixels.
 * Also returns cost = sum(r^2).
 * JtJ is stored as a flat array[64], row-major.
 */
static double compute_jtj_jtr(const double *params,
                               const double *data,
                               const double *x, const double *y,
                               double *JtJ, double *Jtr, int npix)
{
    double Ie    = params[0];
    double re    = params[1];
    double n     = params[2];
    double xc    = params[3];
    double yc    = params[4];
    double ellip = params[5];
    double theta = params[6];

    if (re < 0.1)  re = 0.1;
    if (n  < 0.1)  n  = 0.1;

    double q = 1.0 - ellip;
    if (q < 0.05) q = 0.05;

    double bn       = sersic_bn(n);
    double dbn_dn   = sersic_dbn_dn(n);
    double inv_n    = 1.0 / n;
    double inv_n2   = inv_n * inv_n;
    double inv_re   = 1.0 / re;
    double cos_t    = cos(theta);
    double sin_t    = sin(theta);
    double inv_q2   = 1.0 / (q * q);
    double inv_q3   = 1.0 / (q * q * q);

    /* Zero accumulators */
    double jtj[64] = {0};
    double jtr[8] = {0};
    double cost = 0.0;

    #pragma omp parallel for schedule(static) if(npix > 1000) \
        reduction(+:cost) reduction(+:jtj[:64]) reduction(+:jtr[:8])
    for (int i = 0; i < npix; i++) {
        double dx    = x[i] - xc;
        double dy    = y[i] - yc;
        double x_rot = dx * cos_t + dy * sin_t;
        double y_rot = -dx * sin_t + dy * cos_t;

        double r = sqrt(x_rot * x_rot + y_rot * y_rot * inv_q2);
        int at_center = (r < 1e-10);

        double ratio, rn, E;
        if (at_center) { ratio = 0.0; rn = 0.0; E = -bn; }
        else { ratio = r * inv_re; rn = pow(ratio, inv_n); E = bn * (rn - 1.0); }
        double Ec = clamp(E, -50.0, 50.0);

        double expnE   = exp(-Ec);
        double profile = Ie * expnE;
        double model   = profile + params[7];
        double residual = data[i] - model;
        cost += residual * residual;

        /* Jacobian row: -df/d(param) */
        double dE_dr = at_center ? 0.0 : (bn * rn * inv_n / r);
        double dr_dxc, dr_dyc, dr_dellip, dr_dtheta;
        if (at_center) {
            dr_dxc = 0; dr_dyc = 0; dr_dellip = 0; dr_dtheta = 0;
        } else {
            double inv_r = 1.0 / r;
            dr_dxc    = (x_rot * (-cos_t) + y_rot * sin_t * inv_q2) * inv_r;
            dr_dyc    = (x_rot * (-sin_t) + y_rot * (-cos_t) * inv_q2) * inv_r;
            dr_dellip = (y_rot * y_rot * inv_q3) * inv_r;
            dr_dtheta = x_rot * y_rot * (1.0 - inv_q2) * inv_r;
        }
        double log_ratio = at_center ? 0.0 : log(ratio);
        double neg_dE_dr = -dE_dr;

        double J[8];
        J[0] = -expnE;
        J[1] = -(profile * bn * rn * inv_n * inv_re);
        J[2] = -(profile * (-(dbn_dn * (rn - 1.0)) + bn * rn * log_ratio * inv_n2));
        J[3] = -(profile * neg_dE_dr * dr_dxc);
        J[4] = -(profile * neg_dE_dr * dr_dyc);
        J[5] = -(profile * neg_dE_dr * dr_dellip);
        J[6] = -(profile * neg_dE_dr * dr_dtheta);
        J[7] = -1.0;

        /* Accumulate J^T * J and J^T * r */
        for (int a = 0; a < 8; a++) {
            jtr[a] += J[a] * residual;
            for (int b = a; b < 8; b++)
                jtj[a * 8 + b] += J[a] * J[b];
        }
    }

    /* Fill symmetric part */
    for (int a = 0; a < 8; a++)
        for (int b = 0; b < a; b++)
            jtj[a * 8 + b] = jtj[b * 8 + a];

    for (int i = 0; i < 64; i++) JtJ[i] = jtj[i];
    for (int i = 0; i < 8; i++) Jtr[i] = jtr[i];
    return cost;
}

/* Solve 8x8 linear system Ax=b using Cholesky decomposition.
 * A is modified in place (becomes L). Returns 0 on success, -1 on failure. */
static int solve_8x8(double A[64], double b[8], double x[8])
{
    int N = 8;
    /* Cholesky: A = L * L^T */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            double s = A[i * N + j];
            for (int k = 0; k < j; k++)
                s -= A[i * N + k] * A[j * N + k];
            if (i == j) {
                if (s <= 0) return -1; /* not positive definite */
                A[i * N + j] = sqrt(s);
            } else {
                A[i * N + j] = s / A[j * N + j];
            }
        }
    }

    /* Forward solve: L * y = b */
    for (int i = 0; i < N; i++) {
        double s = b[i];
        for (int k = 0; k < i; k++)
            s -= A[i * N + k] * x[k];
        x[i] = s / A[i * N + i];
    }

    /* Backward solve: L^T * x = y */
    for (int i = N - 1; i >= 0; i--) {
        double s = x[i];
        for (int k = i + 1; k < N; k++)
            s -= A[k * N + i] * x[k];
        x[i] = s / A[i * N + i];
    }

    return 0;
}

/*
 * Full Levenberg-Marquardt fitting in C. One call from Python.
 *
 * params_out[8]: fitted parameters (output)
 * Returns: chi2 = cost / (npix - 8)
 *
 * bounds_lo/hi[8]: parameter bounds (clipped after each step)
 * max_iter: maximum iterations
 *
 * Uses Marquardt scaling (tracked diag max) and gain-ratio lambda update.
 */
double sersic_2d_fit(const double *data, const double *x, const double *y,
                     int npix,
                     const double *params_init,
                     const double *bounds_lo, const double *bounds_hi,
                     int max_iter,
                     double *params_out)
{
    double p[8];
    for (int i = 0; i < 8; i++) p[i] = params_init[i];

    double JtJ[64], Jtr[8], A[64], delta[8];
    double cost = compute_jtj_jtr(p, data, x, y, JtJ, Jtr, npix);

    /* Marquardt scaling: D_ii = max over iterations of diag(J^T*J)_ii */
    double D[8];
    for (int i = 0; i < 8; i++) {
        D[i] = JtJ[i * 8 + i];
        if (D[i] < 1e-10) D[i] = 1e-10;
    }

    /* Initialize lambda from typical scale of J^T*J */
    double max_diag = 0;
    for (int i = 0; i < 8; i++)
        if (D[i] > max_diag) max_diag = D[i];
    double lambda = 1e-3 * max_diag;
    if (lambda < 1e-10) lambda = 1e-10;

    double nu = 2.0;  /* Nielsen factor */

    for (int iter = 0; iter < max_iter; iter++) {
        /* Update Marquardt scaling: D_ii = max(D_ii, diag(J^T*J)_ii) */
        for (int i = 0; i < 8; i++) {
            double d = JtJ[i * 8 + i];
            if (d > D[i]) D[i] = d;
        }

        /* Form A = J^T*J + lambda * diag(D) */
        for (int i = 0; i < 64; i++) A[i] = JtJ[i];
        for (int i = 0; i < 8; i++)
            A[i * 8 + i] += lambda * D[i];

        /* Solve: (J^T*J + lambda*D) * delta = -J^T*r */
        double neg_jtr[8];
        for (int i = 0; i < 8; i++) neg_jtr[i] = -Jtr[i];

        if (solve_8x8(A, neg_jtr, delta) != 0) {
            lambda *= nu;
            nu *= 2.0;
            if (lambda > 1e16) break;
            continue;
        }

        /* Trial step with bounds clipping */
        double p_new[8];
        for (int i = 0; i < 8; i++) {
            p_new[i] = p[i] + delta[i];
            if (p_new[i] < bounds_lo[i]) p_new[i] = bounds_lo[i];
            if (p_new[i] > bounds_hi[i]) p_new[i] = bounds_hi[i];
        }

        /* Evaluate at new point */
        double JtJ_new[64], Jtr_new[8];
        double cost_new = compute_jtj_jtr(p_new, data, x, y,
                                           JtJ_new, Jtr_new, npix);

        /* Gain ratio: actual reduction / predicted reduction */
        /* predicted = delta^T * (lambda*D*delta - Jtr) */
        double predicted = 0;
        for (int i = 0; i < 8; i++)
            predicted += delta[i] * (lambda * D[i] * delta[i] - Jtr[i]);

        double rho = (predicted > 0) ?
            (cost - cost_new) / predicted : -1.0;

        if (cost_new < cost && rho > 0) {
            /* Accept step */
            for (int i = 0; i < 8; i++) p[i] = p_new[i];
            for (int i = 0; i < 64; i++) JtJ[i] = JtJ_new[i];
            for (int i = 0; i < 8; i++) Jtr[i] = Jtr_new[i];

            double reduction = (cost - cost_new) / (cost > 0 ? cost : 1.0);
            cost = cost_new;

            /* Nielsen lambda update */
            double tmp = 2.0 * rho - 1.0;
            double scale = 1.0 - tmp * tmp * tmp;
            if (scale < 1.0 / 3.0) scale = 1.0 / 3.0;
            lambda *= scale;
            if (lambda < 1e-15) lambda = 1e-15;
            nu = 2.0;

            /* Convergence: relative cost reduction + gradient norm */
            if (reduction < 1e-10 && iter > 10) {
                double gnorm = 0;
                for (int i = 0; i < 8; i++) gnorm += Jtr[i] * Jtr[i];
                if (gnorm < 1e-10 * cost * cost) break;
            }
        } else {
            /* Reject step, increase damping */
            lambda *= nu;
            nu *= 2.0;
            if (lambda > 1e16) break;
        }
    }

    for (int i = 0; i < 8; i++) params_out[i] = p[i];
    return cost / (npix > 8 ? npix - 8 : 1);
}
