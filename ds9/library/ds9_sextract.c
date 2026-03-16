/*
 * ds9_sextract.c - SExtractor-like source extraction for DS9
 *
 * Uses: sep (C library), cfitsio (FITS I/O)
 *
 * Output parameters (tab-separated):
 *   Photometry:        MAG_AUTO, MAG_ISOCOR, MAG_APER, FLUX_AUTO,
 *                      FLUXERR_AUTO, FLUXERR_APER, MAGERR_AUTO, MAGERR_APER,
 *                      FLUX_APER_2, FLUX_APER_3, FLUX_APER_5,
 *                      FLUXERR_APER_2, FLUXERR_APER_3, FLUXERR_APER_5,
 *                      MAG_APER_2, MAG_APER_3, MAG_APER_5
 *   Morphology & Size: FLUX_RADIUS, A_IMAGE, B_IMAGE, THETA_IMAGE,
 *                      ELLIPTICITY, KRON_RADIUS, FWHM_IMAGE
 *   Surface Brightness: MU_MAX, MU_THRESHOLD
 *   Astrometry:        X_IMAGE, Y_IMAGE, ALPHA_J2000, DELTA_J2000
 *   Classification:    CLASS_STAR, FLAGS
 *
 * Build:
 *   gcc -O2 -o ds9_sextract ds9_sextract.c \
 *       -I/path/to/sep/src -L/path/to/sep -lsep \
 *       -lcfitsio -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "fitsio.h"
#include "sep.h"

/* Default parameters */
#define DEF_DETECT_THRESH   1.5
#define DEF_DETECT_MINAREA  5
#define DEF_DEBLEND_NTHRESH 32
#define DEF_DEBLEND_MINCONT 0.005
#define DEF_PHOT_APERTURE   5.0
#define DEF_PHOT_APERTURE_2 4.0
#define DEF_PHOT_APERTURE_3 6.0
#define DEF_PHOT_APERTURE_5 10.0
#define DEF_MAG_ZEROPOINT   25.0
#define DEF_GAIN            0.0
#define DEF_PIXEL_SCALE     1.0
#define DEF_SEEING_FWHM     3.0
#define DEF_BACK_SIZE       64
#define DEF_BACK_FILTERSIZE 3
#define DEF_SATUR_LEVEL     65535.0

/* Structures */
typedef struct {
    double detect_thresh;
    int    detect_minarea;
    int    deblend_nthresh;
    double deblend_mincont;
    double phot_aperture;
    double phot_aperture_2;   /* 2nd aperture diameter (pixels) */
    double phot_aperture_3;   /* 3rd aperture diameter (pixels) */
    double phot_aperture_5;   /* 4th aperture diameter (pixels) */
    double mag_zeropoint;
    double gain;
    double pixel_scale;
    double seeing_fwhm;
    int    back_size;
    int    back_filtersize;
    double satur_level;
    int    conv_filter;       /* 0=default3x3, 1=gauss5x5, 2=mexhat, 3=tophat */
} config_t;

typedef struct {
    int    number;
    double x_image, y_image;
    double alpha_j2000, delta_j2000;
    double mag_auto, mag_isocor, mag_aper;
    double flux_auto;
    double fluxerr_auto, fluxerr_aper;
    double magerr_auto, magerr_aper;
    double flux_aper_2, flux_aper_3, flux_aper_5;
    double fluxerr_aper_2, fluxerr_aper_3, fluxerr_aper_5;
    double mag_aper_2, mag_aper_3, mag_aper_5;
    double flux_radius;
    double a_image, b_image, theta_image;
    double ellipticity;
    double kron_radius;
    double fwhm_image;
    double mu_max, mu_threshold;
    double class_star;
    int    flags;
    long   npix_iso;
    double iso_radius;
} source_t;

/* Comparison function for sorting by MAG_AUTO (brightest first) */
static int compare_mag(const void *a, const void *b) {
    double ma = ((const source_t *)a)->mag_auto;
    double mb = ((const source_t *)b)->mag_auto;
    if (ma < mb) return -1;
    if (ma > mb) return  1;
    return 0;
}

/* Simple WCS: read CD matrix or CDELT from FITS header */
typedef struct {
    int    valid;
    double crpix1, crpix2;
    double crval1, crval2;
    double cd11, cd12, cd21, cd22;
} simple_wcs_t;

static void read_wcs(fitsfile *fptr, simple_wcs_t *wcs) {
    int status = 0;
    double cdelt1 = 0, cdelt2 = 0, crota2 = 0;
    double pc11 = 1, pc12 = 0, pc21 = 0, pc22 = 1;

    memset(wcs, 0, sizeof(simple_wcs_t));
    wcs->valid = 0;

    /* Read reference pixel */
    fits_read_key(fptr, TDOUBLE, "CRPIX1", &wcs->crpix1, NULL, &status);
    if (status) return;
    status = 0;
    fits_read_key(fptr, TDOUBLE, "CRPIX2", &wcs->crpix2, NULL, &status);
    if (status) return;
    status = 0;
    fits_read_key(fptr, TDOUBLE, "CRVAL1", &wcs->crval1, NULL, &status);
    if (status) return;
    status = 0;
    fits_read_key(fptr, TDOUBLE, "CRVAL2", &wcs->crval2, NULL, &status);
    if (status) return;
    status = 0;

    /* Try CD matrix first */
    fits_read_key(fptr, TDOUBLE, "CD1_1", &wcs->cd11, NULL, &status);
    if (status == 0) {
        status = 0;
        fits_read_key(fptr, TDOUBLE, "CD1_2", &wcs->cd12, NULL, &status);
        if (status) { status = 0; wcs->cd12 = 0; }
        fits_read_key(fptr, TDOUBLE, "CD2_1", &wcs->cd21, NULL, &status);
        if (status) { status = 0; wcs->cd21 = 0; }
        fits_read_key(fptr, TDOUBLE, "CD2_2", &wcs->cd22, NULL, &status);
        if (status) { status = 0; wcs->cd22 = 0; }
        wcs->valid = 1;
        return;
    }

    /* Fall back to CDELT + CROTA2 or PC matrix */
    status = 0;
    fits_read_key(fptr, TDOUBLE, "CDELT1", &cdelt1, NULL, &status);
    if (status) return;
    status = 0;
    fits_read_key(fptr, TDOUBLE, "CDELT2", &cdelt2, NULL, &status);
    if (status) return;
    status = 0;

    /* Try PC matrix */
    fits_read_key(fptr, TDOUBLE, "PC1_1", &pc11, NULL, &status);
    if (status == 0) {
        status = 0;
        fits_read_key(fptr, TDOUBLE, "PC1_2", &pc12, NULL, &status);
        if (status) { status = 0; pc12 = 0; }
        fits_read_key(fptr, TDOUBLE, "PC2_1", &pc21, NULL, &status);
        if (status) { status = 0; pc21 = 0; }
        fits_read_key(fptr, TDOUBLE, "PC2_2", &pc22, NULL, &status);
        if (status) { status = 0; pc22 = 1; }
        wcs->cd11 = cdelt1 * pc11;
        wcs->cd12 = cdelt1 * pc12;
        wcs->cd21 = cdelt2 * pc21;
        wcs->cd22 = cdelt2 * pc22;
    } else {
        /* CDELT + CROTA2 */
        status = 0;
        fits_read_key(fptr, TDOUBLE, "CROTA2", &crota2, NULL, &status);
        if (status) { status = 0; crota2 = 0; }
        double rad = crota2 * M_PI / 180.0;
        wcs->cd11 = cdelt1 * cos(rad);
        wcs->cd12 = -cdelt2 * sin(rad);
        wcs->cd21 = cdelt1 * sin(rad);
        wcs->cd22 = cdelt2 * cos(rad);
    }
    wcs->valid = 1;
}

/* Convert pixel to world (simple TAN projection) */
static void pix2world(const simple_wcs_t *wcs, double xpix, double ypix,
                       double *ra, double *dec) {
    if (!wcs->valid) {
        *ra = 0.0;
        *dec = 0.0;
        return;
    }

    double dx = xpix - wcs->crpix1;
    double dy = ypix - wcs->crpix2;

    /* Intermediate world coordinates (degrees) */
    double xi  = wcs->cd11 * dx + wcs->cd12 * dy;
    double eta = wcs->cd21 * dx + wcs->cd22 * dy;

    /* TAN (gnomonic) deprojection */
    double ra0  = wcs->crval1 * M_PI / 180.0;
    double dec0 = wcs->crval2 * M_PI / 180.0;
    double xi_r  = xi  * M_PI / 180.0;
    double eta_r = eta * M_PI / 180.0;

    double denom = cos(dec0) - eta_r * sin(dec0);
    double a = atan2(xi_r, denom) + ra0;
    double d = atan2((eta_r * cos(dec0) + sin(dec0)) * cos(a - ra0), denom);

    *ra  = a * 180.0 / M_PI;
    *dec = d * 180.0 / M_PI;

    /* Normalize RA to [0, 360) */
    while (*ra < 0)   *ra += 360.0;
    while (*ra >= 360) *ra -= 360.0;
}

/* Convolution filter kernels */
static float conv_default[9] = {
    1, 2, 1,
    2, 4, 2,
    1, 2, 1
};

static float conv_gauss5x5[25] = {
    1,  4,  7,  4, 1,
    4, 16, 26, 16, 4,
    7, 26, 41, 26, 7,
    4, 16, 26, 16, 4,
    1,  4,  7,  4, 1
};

static float conv_mexhat[25] = {
    -1, -1, -1, -1, -1,
    -1,  1,  1,  1, -1,
    -1,  1,  8,  1, -1,
    -1,  1,  1,  1, -1,
    -1, -1, -1, -1, -1
};

static float conv_tophat[25] = {
    0, 1, 1, 1, 0,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    0, 1, 1, 1, 0
};

int main(int argc, char *argv[]) {
    config_t cfg;
    fitsfile *fptr = NULL;
    int status = 0;
    long naxes[2];
    int naxis, bitpix;
    double *data = NULL;
    double *bkg_arr = NULL;
    double *rms_arr = NULL;
    unsigned char *mask_arr = NULL;

    /* Defaults */
    cfg.detect_thresh   = DEF_DETECT_THRESH;
    cfg.detect_minarea  = DEF_DETECT_MINAREA;
    cfg.deblend_nthresh = DEF_DEBLEND_NTHRESH;
    cfg.deblend_mincont = DEF_DEBLEND_MINCONT;
    cfg.phot_aperture   = DEF_PHOT_APERTURE;
    cfg.phot_aperture_2 = DEF_PHOT_APERTURE_2;
    cfg.phot_aperture_3 = DEF_PHOT_APERTURE_3;
    cfg.phot_aperture_5 = DEF_PHOT_APERTURE_5;
    cfg.mag_zeropoint   = DEF_MAG_ZEROPOINT;
    cfg.gain            = DEF_GAIN;
    cfg.pixel_scale     = DEF_PIXEL_SCALE;
    cfg.seeing_fwhm     = DEF_SEEING_FWHM;
    cfg.back_size       = DEF_BACK_SIZE;
    cfg.back_filtersize = DEF_BACK_FILTERSIZE;
    cfg.satur_level     = DEF_SATUR_LEVEL;
    cfg.conv_filter     = 0;

    /* Parse arguments */
    if (argc < 2) {
        fprintf(stderr,
            "Usage: ds9_sextract <fits_file> [options]\n"
            "Options:\n"
            "  --detect-thresh <f>    Detection threshold in sigma (%.1f)\n"
            "  --detect-minarea <n>   Min detection area in pixels (%d)\n"
            "  --deblend-nthresh <n>  Deblending sub-thresholds (%d)\n"
            "  --deblend-mincont <f>  Deblending min contrast (%.4f)\n"
            "  --phot-aperture <f>    Aperture diameter in pixels (%.1f)\n"
            "  --mag-zeropoint <f>    Magnitude zero-point (%.1f)\n"
            "  --gain <f>             Gain in e-/ADU (%.1f)\n"
            "  --pixel-scale <f>      Pixel scale in arcsec (%.1f)\n"
            "  --seeing-fwhm <f>      Seeing FWHM in arcsec (%.1f)\n"
            "  --back-size <n>        Background mesh size (%d)\n"
            "  --back-filtersize <n>  Background filter size (%d)\n"
            "  --phot-aperture-2 <f>  2nd aperture diameter (%.1f)\n"
            "  --phot-aperture-3 <f>  3rd aperture diameter (%.1f)\n"
            "  --phot-aperture-5 <f>  4th aperture diameter (%.1f)\n"
            "  --conv-filter <s>      Convolution filter: default|gauss5x5|mexhat|tophat\n",
            DEF_DETECT_THRESH, DEF_DETECT_MINAREA, DEF_DEBLEND_NTHRESH,
            DEF_DEBLEND_MINCONT, DEF_PHOT_APERTURE, DEF_MAG_ZEROPOINT,
            DEF_GAIN, DEF_PIXEL_SCALE, DEF_SEEING_FWHM,
            DEF_BACK_SIZE, DEF_BACK_FILTERSIZE,
            DEF_PHOT_APERTURE_2, DEF_PHOT_APERTURE_3, DEF_PHOT_APERTURE_5);
        return 1;
    }

    const char *fitsfile = argv[1];
    for (int i = 2; i < argc - 1; i += 2) {
        if      (!strcmp(argv[i], "--detect-thresh"))   cfg.detect_thresh   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "--detect-minarea"))  cfg.detect_minarea  = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "--deblend-nthresh")) cfg.deblend_nthresh = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "--deblend-mincont")) cfg.deblend_mincont = atof(argv[i+1]);
        else if (!strcmp(argv[i], "--phot-aperture"))   cfg.phot_aperture   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "--mag-zeropoint"))   cfg.mag_zeropoint   = atof(argv[i+1]);
        else if (!strcmp(argv[i], "--gain"))            cfg.gain            = atof(argv[i+1]);
        else if (!strcmp(argv[i], "--pixel-scale"))     cfg.pixel_scale     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "--seeing-fwhm"))     cfg.seeing_fwhm     = atof(argv[i+1]);
        else if (!strcmp(argv[i], "--back-size"))       cfg.back_size       = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "--back-filtersize")) cfg.back_filtersize = atoi(argv[i+1]);
        else if (!strcmp(argv[i], "--phot-aperture-2")) cfg.phot_aperture_2 = atof(argv[i+1]);
        else if (!strcmp(argv[i], "--phot-aperture-3")) cfg.phot_aperture_3 = atof(argv[i+1]);
        else if (!strcmp(argv[i], "--phot-aperture-5")) cfg.phot_aperture_5 = atof(argv[i+1]);
        else if (!strcmp(argv[i], "--conv-filter")) {
            if      (!strcmp(argv[i+1], "gauss5x5")) cfg.conv_filter = 1;
            else if (!strcmp(argv[i+1], "mexhat"))    cfg.conv_filter = 2;
            else if (!strcmp(argv[i+1], "tophat"))    cfg.conv_filter = 3;
            else                                     cfg.conv_filter = 0;
        }
    }

    /* ================================================================
     * Step 1: Read FITS image
     * ================================================================ */
    fits_open_file(&fptr, fitsfile, READONLY, &status);
    if (status) {
        char errmsg[80];
        fits_get_errstatus(status, errmsg);
        fprintf(stderr, "ERROR: Cannot open %s: %s\n", fitsfile, errmsg);
        return 1;
    }

    /* Try primary HDU first, then search extensions for 2D image */
    fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status);
    if (status || naxis < 2) {
        /* Primary has no image data; search extensions */
        status = 0;
        int nhdu = 0;
        fits_get_num_hdus(fptr, &nhdu, &status);
        int found = 0;
        for (int hh = 2; hh <= nhdu; hh++) {
            int hdutype = 0;
            fits_movabs_hdu(fptr, hh, &hdutype, &status);
            if (status || hdutype != IMAGE_HDU) { status = 0; continue; }
            fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status);
            if (status == 0 && naxis >= 2 && naxes[0] > 1 && naxes[1] > 1) {
                found = 1;
                fprintf(stderr, "INFO: Using image extension HDU %d (%ldx%ld)\n",
                        hh, naxes[0], naxes[1]);
                break;
            }
            status = 0;
        }
        if (!found) {
            fprintf(stderr, "ERROR: No 2D image data found in any HDU\n");
            fits_close_file(fptr, &status);
            return 1;
        }
    }

    long nx = naxes[0];
    long ny = naxes[1];
    long npix = nx * ny;

    /* Read WCS */
    simple_wcs_t wcs;
    read_wcs(fptr, &wcs);

    /* Get pixel scale from WCS if possible */
    if (wcs.valid && cfg.pixel_scale == DEF_PIXEL_SCALE) {
        double ps = sqrt(fabs(wcs.cd11 * wcs.cd22 - wcs.cd12 * wcs.cd21));
        if (ps > 0) cfg.pixel_scale = ps * 3600.0;  /* deg to arcsec */
    }

    /* Try to get gain from header */
    if (cfg.gain == 0.0) {
        double hgain = 0;
        int s2 = 0;
        fits_read_key(fptr, TDOUBLE, "GAIN", &hgain, NULL, &s2);
        if (s2) { s2 = 0; fits_read_key(fptr, TDOUBLE, "EGAIN", &hgain, NULL, &s2); }
        if (s2 == 0 && hgain > 0) cfg.gain = hgain;
    }

    /* Try to get saturation from header */
    {
        double hsat = 0;
        int s2 = 0;
        fits_read_key(fptr, TDOUBLE, "SATURATE", &hsat, NULL, &s2);
        if (s2 == 0 && hsat > 0) cfg.satur_level = hsat;
    }

    /* Read image data */
    data = (double *)malloc(npix * sizeof(double));
    if (!data) {
        fprintf(stderr, "ERROR: Out of memory\n");
        fits_close_file(fptr, &status);
        return 1;
    }

    long fpixel[2] = {1, 1};
    int anynul = 0;
    fits_read_pix(fptr, TDOUBLE, fpixel, npix, NULL, data, &anynul, &status);
    if (status) {
        fprintf(stderr, "ERROR: Cannot read image data\n");
        free(data);
        fits_close_file(fptr, &status);
        return 1;
    }
    fits_close_file(fptr, &status);

    /* Build mask: mark pixels that are NaN, Inf, or exactly 0 in the
     * original data.  These are outside the actual observation footprint
     * (common in JWST/HST mosaics where unexposed regions are zero-filled). */
    mask_arr = (unsigned char *)calloc(npix, sizeof(unsigned char));
    if (!mask_arr) {
        fprintf(stderr, "ERROR: Out of memory for mask\n");
        free(data);
        fits_close_file(fptr, &status);
        return 1;
    }
    {
        long n_masked = 0;
        for (long i = 0; i < npix; i++) {
            if (isnan(data[i]) || isinf(data[i]) || data[i] == 0.0) {
                mask_arr[i] = 1;
                n_masked++;
            }
        }
        if (n_masked > 0)
            fprintf(stderr, "Masked %ld / %ld pixels (outside observation)\n",
                    n_masked, npix);
    }

    /* Replace NaN/Inf with 0 for safe arithmetic (masked anyway) */
    for (long i = 0; i < npix; i++) {
        if (isnan(data[i]) || isinf(data[i])) data[i] = 0.0;
    }

    /* ================================================================
     * Step 2: Background estimation
     * ================================================================ */
    sep_bkg *bkg = NULL;
    sep_image im;
    memset(&im, 0, sizeof(sep_image));
    im.data = data;
    im.dtype = SEP_TDOUBLE;
    im.w = nx;
    im.h = ny;
    im.noise = NULL;
    im.ndtype = 0;
    im.mask = mask_arr;
    im.mdtype = SEP_TBYTE;
    im.noiseval = 0;
    im.noise_type = SEP_NOISE_NONE;
    im.gain = (cfg.gain > 0) ? cfg.gain : 0.0;
    im.maskthresh = 0.5;

    int sep_status = sep_background(&im,
        cfg.back_size, cfg.back_size,
        cfg.back_filtersize, cfg.back_filtersize,
        0.0, &bkg);
    if (sep_status) {
        char errtxt[512];
        sep_get_errmsg(sep_status, errtxt);
        fprintf(stderr, "ERROR: Background estimation failed: %s\n", errtxt);
        free(data); free(mask_arr);
        return 1;
    }

    /* Subtract background in-place */
    sep_bkg_subarray(bkg, data, SEP_TDOUBLE);

    /* Get background RMS array */
    rms_arr = (double *)malloc(npix * sizeof(double));
    sep_bkg_rmsarray(bkg, rms_arr, SEP_TDOUBLE);

    /* Also keep background array for peak correction */
    bkg_arr = (double *)malloc(npix * sizeof(double));
    sep_bkg_array(bkg, bkg_arr, SEP_TDOUBLE);

    sep_bkg_free(bkg);

    /* ================================================================
     * Step 3: Source detection
     * ================================================================ */
    sep_catalog *cat = NULL;

    /* Set noise in image struct for extraction */
    im.data = data;  /* now background-subtracted */
    im.noise = rms_arr;
    im.ndtype = SEP_TDOUBLE;
    im.noise_type = SEP_NOISE_STDDEV;

    /* Increase pixel stack for large images */
    {
        size_t pstack = (size_t)nx * (size_t)ny;
        if (pstack < 500000) pstack = 500000;
        if (pstack > 50000000) pstack = 50000000;
        sep_set_extract_pixstack(pstack);
    }

    /* Select convolution kernel */
    float *conv_kern;
    int conv_w, conv_h;
    switch (cfg.conv_filter) {
    case 1:  conv_kern = conv_gauss5x5; conv_w = 5; conv_h = 5; break;
    case 2:  conv_kern = conv_mexhat;   conv_w = 5; conv_h = 5; break;
    case 3:  conv_kern = conv_tophat;   conv_w = 5; conv_h = 5; break;
    default: conv_kern = conv_default;  conv_w = 3; conv_h = 3; break;
    }

    /* Retry extraction with increasing threshold on overflow */
    double try_thresh = cfg.detect_thresh;
    int max_retries = 5;
    for (int retry = 0; retry <= max_retries; retry++) {
        sep_status = sep_extract(&im,
            try_thresh,
            SEP_THRESH_REL,
            cfg.detect_minarea,
            conv_kern, conv_w, conv_h,
            SEP_FILTER_CONV,
            cfg.deblend_nthresh,
            cfg.deblend_mincont,
            1, 1.0,
            &cat);

        if (sep_status == 0) break;
        if (cat) { sep_catalog_free(cat); cat = NULL; }

        /* On overflow, increase threshold and retry */
        try_thresh *= 1.5;
        fprintf(stderr, "WARNING: Extraction overflow, retrying with thresh=%.2f\n",
                try_thresh);
    }

    if (sep_status) {
        char errtxt[512];
        sep_get_errmsg(sep_status, errtxt);
        fprintf(stderr, "ERROR: Extraction failed after retries: %s\n", errtxt);
        free(data); free(rms_arr); free(bkg_arr); free(mask_arr);
        return 1;
    }

    int nobj = cat->nobj;

    /* ================================================================
     * Step 4: Measure all parameters
     * ================================================================ */
    source_t *sources = (source_t *)calloc(nobj, sizeof(source_t));
    if (!sources) {
        fprintf(stderr, "ERROR: Out of memory for %d sources\n", nobj);
        sep_catalog_free(cat);
        free(data); free(rms_arr); free(bkg_arr); free(mask_arr);
        return 1;
    }

    double pix_area = cfg.pixel_scale * cfg.pixel_scale;  /* arcsec^2 */

    for (int i = 0; i < nobj; i++) {
        source_t *s = &sources[i];
        s->number = i + 1;

        /* Position (1-indexed FITS convention) */
        s->x_image = cat->x[i] + 1.0;
        s->y_image = cat->y[i] + 1.0;

        /* WCS */
        pix2world(&wcs, s->x_image, s->y_image,
                  &s->alpha_j2000, &s->delta_j2000);

        /* Shape */
        s->a_image = cat->a[i];
        s->b_image = cat->b[i];
        s->theta_image = cat->theta[i] * 180.0 / M_PI;
        s->ellipticity = (s->a_image > 0) ? 1.0 - s->b_image / s->a_image : 0.0;
        if (s->ellipticity < 0) s->ellipticity = 0;
        if (s->ellipticity > 1) s->ellipticity = 1;

        /* Flags */
        s->flags = (int)cat->flag[i];

        /* ---- Kron/AUTO photometry ---- */
        double kronrad = 0;
        short krflag = 0;
        /* sep_kron_radius uses cxx,cyy,cxy ellipse parameterization */
        sep_kron_radius(&im,
            cat->x[i], cat->y[i],
            cat->cxx[i], cat->cyy[i], cat->cxy[i],
            6.0,    /* r */
            0,      /* id (no segmap filtering) */
            &kronrad, &krflag);

        if (kronrad < 3.5) kronrad = 3.5;
        s->kron_radius = kronrad;

        double flux_auto = 0, fluxerr_auto = 0, area_auto = 0;
        short aflag = 0;
        sep_sum_ellipse(&im,
            cat->x[i], cat->y[i],
            cat->a[i], cat->b[i], cat->theta[i],
            2.5 * kronrad,
            0,      /* id */
            5,      /* subpix */
            0,      /* inflags */
            &flux_auto, &fluxerr_auto, &area_auto, &aflag);

        s->flux_auto = flux_auto;
        s->fluxerr_auto = fluxerr_auto;
        s->flags |= aflag;

        /* ---- Isophotal corrected ---- */
        double flux_iso = cat->flux[i];
        long npix_iso = (long)cat->tnpix[i];
        double local_thresh = cat->thresh[i];

        double ati = (flux_iso > 0) ?
            (double)npix_iso * local_thresh / flux_iso : 0;
        double corr = 1.0 - 0.196 * ati - 0.751 * ati * ati;
        if (corr < 0.01) corr = 0.01;
        if (corr > 1.0) corr = 1.0;
        double flux_isocor = (flux_iso > 0) ? flux_iso / corr : flux_iso;

        /* ---- Isophotal radius ---- */
        s->npix_iso = npix_iso;
        /* Semi-major axis of ellipse with same area as isophotal detection */
        /* area = pi * a * b = pi * a * a * (b/a) => a = sqrt(area / (pi * b/a)) */
        if (s->a_image > 0 && s->b_image > 0 && npix_iso > 0) {
            double ba_ratio = s->b_image / s->a_image;
            s->iso_radius = sqrt((double)npix_iso / (M_PI * ba_ratio));
        } else {
            s->iso_radius = sqrt((double)npix_iso / M_PI);
        }

        /* ---- Aperture photometry ---- */
        double flux_aper = 0, fluxerr_aper = 0, area_aper = 0;
        short apflag = 0;
        sep_sum_circle(&im,
            cat->x[i], cat->y[i],
            cfg.phot_aperture / 2.0,  /* radius */
            0,      /* id */
            5,      /* subpix */
            0,      /* inflags */
            &flux_aper, &fluxerr_aper, &area_aper, &apflag);
        s->fluxerr_aper = fluxerr_aper;
        s->flags |= apflag;

        /* ---- Multi-aperture photometry (3 extra apertures) ---- */
        {
            double fa2 = 0, fe2 = 0, aa2 = 0;
            short af2 = 0;
            sep_sum_circle(&im, cat->x[i], cat->y[i],
                cfg.phot_aperture_2 / 2.0, 0, 5, 0,
                &fa2, &fe2, &aa2, &af2);
            s->flux_aper_2 = fa2;
            s->fluxerr_aper_2 = fe2;
            s->flags |= af2;

            double fa3 = 0, fe3 = 0, aa3 = 0;
            short af3 = 0;
            sep_sum_circle(&im, cat->x[i], cat->y[i],
                cfg.phot_aperture_3 / 2.0, 0, 5, 0,
                &fa3, &fe3, &aa3, &af3);
            s->flux_aper_3 = fa3;
            s->fluxerr_aper_3 = fe3;
            s->flags |= af3;

            double fa5 = 0, fe5 = 0, aa5 = 0;
            short af5 = 0;
            sep_sum_circle(&im, cat->x[i], cat->y[i],
                cfg.phot_aperture_5 / 2.0, 0, 5, 0,
                &fa5, &fe5, &aa5, &af5);
            s->flux_aper_5 = fa5;
            s->fluxerr_aper_5 = fe5;
            s->flags |= af5;
        }

        /* ---- Half-light radius ---- */
        double flux_rad = 0;
        short frflag = 0;
        double rmax = 6.0 * cat->a[i];
        if (rmax < 10.0) rmax = 10.0;
        if (rmax > 200.0) rmax = 200.0;
        double frac = 0.5;
        sep_flux_radius(&im,
            cat->x[i], cat->y[i],
            rmax,
            0,      /* id */
            5,      /* subpix */
            0,      /* inflag */
            &flux_auto,  /* fluxtot (pointer) */
            &frac,       /* fluxfrac (pointer) */
            1,           /* n (number of fractions) */
            &flux_rad, &frflag);
        s->flux_radius = flux_rad;

        /* ---- Magnitudes ---- */
        s->mag_auto = (flux_auto > 0) ?
            -2.5 * log10(flux_auto) + cfg.mag_zeropoint : 99.0;
        s->mag_isocor = (flux_isocor > 0) ?
            -2.5 * log10(flux_isocor) + cfg.mag_zeropoint : 99.0;
        s->mag_aper = (flux_aper > 0) ?
            -2.5 * log10(flux_aper) + cfg.mag_zeropoint : 99.0;

        /* ---- Magnitude errors ---- */
        s->magerr_auto = (flux_auto > 0 && fluxerr_auto > 0) ?
            1.0857 * fluxerr_auto / flux_auto : 99.0;
        s->magerr_aper = (flux_aper > 0 && fluxerr_aper > 0) ?
            1.0857 * fluxerr_aper / flux_aper : 99.0;

        /* ---- Multi-aperture magnitudes ---- */
        s->mag_aper_2 = (s->flux_aper_2 > 0) ?
            -2.5 * log10(s->flux_aper_2) + cfg.mag_zeropoint : 99.0;
        s->mag_aper_3 = (s->flux_aper_3 > 0) ?
            -2.5 * log10(s->flux_aper_3) + cfg.mag_zeropoint : 99.0;
        s->mag_aper_5 = (s->flux_aper_5 > 0) ?
            -2.5 * log10(s->flux_aper_5) + cfg.mag_zeropoint : 99.0;

        /* ---- Surface brightness ---- */
        double peak = cat->peak[i];
        s->mu_max = (peak > 0 && pix_area > 0) ?
            -2.5 * log10(peak / pix_area) + cfg.mag_zeropoint : 99.0;
        s->mu_threshold = (local_thresh > 0 && pix_area > 0) ?
            -2.5 * log10(local_thresh / pix_area) + cfg.mag_zeropoint : 99.0;

        /* ---- Star/Galaxy classification ---- */
        double fwhm = 2.0 * sqrt(log(2.0) * (s->a_image * s->a_image +
                                               s->b_image * s->b_image));
        s->fwhm_image = fwhm;
        double seeing_pix = cfg.seeing_fwhm / cfg.pixel_scale;
        double ratio = (seeing_pix > 0) ? fwhm / seeing_pix : 2.0;
        s->class_star = 1.0 / (1.0 + exp(3.0 * (ratio - 1.5)));
        if (s->class_star < 0) s->class_star = 0;
        if (s->class_star > 1) s->class_star = 1;

        /* ---- Saturation flag ---- */
        int ix = (int)(cat->x[i] + 0.5);
        int iy = (int)(cat->y[i] + 0.5);
        if (ix >= 0 && ix < nx && iy >= 0 && iy < ny) {
            double raw_peak = peak + bkg_arr[iy * nx + ix];
            if (raw_peak >= cfg.satur_level) {
                s->flags |= 4;  /* OBJ_SATUR */
            }
        }
    }

    /* ================================================================
     * Step 5: Filter out sources overlapping observation boundary
     * ================================================================ */
    /* For each source, check if any pixel within its ISO_RADIUS ellipse
     * falls on a masked pixel (outside observation area). */
    int nvalid = 0;
    for (int i = 0; i < nobj; i++) {
        source_t *s = &sources[i];
        /* 0-based center */
        double cx = s->x_image - 1.0;
        double cy = s->y_image - 1.0;

        /* Use ISO_RADIUS (isophotal detection footprint) */
        double semi_a = s->iso_radius;
        double semi_b = semi_a;
        if (s->a_image > 0 && s->b_image > 0)
            semi_b = semi_a * s->b_image / s->a_image;
        if (semi_a < 3.0) semi_a = 3.0;
        if (semi_b < 3.0) semi_b = 3.0;
        double theta_rad = s->theta_image * M_PI / 180.0;
        double cos_t = cos(theta_rad);
        double sin_t = sin(theta_rad);

        int has_masked = 0;
        int ir = (int)(semi_a + 1.5);
        for (int dy = -ir; dy <= ir && !has_masked; dy++) {
            for (int dx = -ir; dx <= ir && !has_masked; dx++) {
                /* Check if (dx,dy) is inside the ellipse */
                double u = dx * cos_t + dy * sin_t;
                double v = -dx * sin_t + dy * cos_t;
                if ((u * u) / (semi_a * semi_a) +
                    (v * v) / (semi_b * semi_b) > 1.0) continue;

                int px = (int)(cx + 0.5) + dx;
                int py = (int)(cy + 0.5) + dy;
                if (px < 0 || px >= nx || py < 0 || py >= ny) {
                    has_masked = 1;
                    break;
                }
                if (mask_arr[py * nx + px]) {
                    has_masked = 1;
                }
            }
        }
        if (has_masked) {
            s->number = -1;  /* Mark for removal */
        } else {
            nvalid++;
        }
    }
    fprintf(stderr, "Filtered: %d / %d sources outside observation boundary\n",
            nobj - nvalid, nobj);

    /* ================================================================
     * Step 6: Sort by brightness and output
     * ================================================================ */
    qsort(sources, nobj, sizeof(source_t), compare_mag);

    /* Header */
    printf("NUMBER\tX_IMAGE\tY_IMAGE\tALPHA_J2000\tDELTA_J2000\t"
           "MAG_AUTO\tMAG_ISOCOR\tMAG_APER\tFLUX_AUTO\t"
           "FLUXERR_AUTO\tFLUXERR_APER\tMAGERR_AUTO\tMAGERR_APER\t"
           "FLUX_APER_2\tFLUX_APER_3\tFLUX_APER_5\t"
           "FLUXERR_APER_2\tFLUXERR_APER_3\tFLUXERR_APER_5\t"
           "MAG_APER_2\tMAG_APER_3\tMAG_APER_5\t"
           "FLUX_RADIUS\tA_IMAGE\tB_IMAGE\tTHETA_IMAGE\tELLIPTICITY\t"
           "KRON_RADIUS\tFWHM_IMAGE\tISO_RADIUS\tNPIX_ISO\t"
           "MU_MAX\tMU_THRESHOLD\tCLASS_STAR\tFLAGS\n");

    /* Data rows — skip sources marked for removal */
    int out_num = 0;
    for (int i = 0; i < nobj; i++) {
        source_t *s = &sources[i];
        if (s->number < 0) continue;  /* Filtered out */
        out_num++;
        s->number = out_num;  /* Renumber sequentially */
        printf("%d\t%.2f\t%.2f\t%.6f\t%.6f\t"
               "%.3f\t%.3f\t%.3f\t%.2f\t"
               "%.2f\t%.2f\t%.4f\t%.4f\t"
               "%.2f\t%.2f\t%.2f\t"
               "%.2f\t%.2f\t%.2f\t"
               "%.3f\t%.3f\t%.3f\t"
               "%.2f\t%.2f\t%.2f\t%.1f\t%.3f\t"
               "%.2f\t%.2f\t%.2f\t%ld\t"
               "%.2f\t%.2f\t%.3f\t%d\n",
               s->number, s->x_image, s->y_image,
               s->alpha_j2000, s->delta_j2000,
               s->mag_auto, s->mag_isocor, s->mag_aper, s->flux_auto,
               s->fluxerr_auto, s->fluxerr_aper,
               s->magerr_auto, s->magerr_aper,
               s->flux_aper_2, s->flux_aper_3, s->flux_aper_5,
               s->fluxerr_aper_2, s->fluxerr_aper_3, s->fluxerr_aper_5,
               s->mag_aper_2, s->mag_aper_3, s->mag_aper_5,
               s->flux_radius,
               s->a_image, s->b_image, s->theta_image, s->ellipticity,
               s->kron_radius, s->fwhm_image, s->iso_radius, s->npix_iso,
               s->mu_max, s->mu_threshold,
               s->class_star, s->flags);
    }

    /* Cleanup */
    free(sources);
    sep_catalog_free(cat);
    free(data);
    free(rms_arr);
    free(bkg_arr);
    free(mask_arr);

    return 0;
}
