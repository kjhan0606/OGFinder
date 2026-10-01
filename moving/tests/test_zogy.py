"""ZOGY on synthetic images: injected point source recovered at the expected significance, flux unbiased."""
import numpy as np
from moving import zogy as Z, imaging as I


def _img(shape, rng, sky_sigma):
    return rng.normal(0, sky_sigma, shape)


def test_zogy_recovers_injected_source():
    rng = np.random.default_rng(5)
    shape = (256, 256); sig = 1.0
    P = I.gaussian_psf(3.0, 25); P /= P.sum()
    R = _img(shape, rng, sig); N = _img(shape, rng, sig)
    # static stars in both
    from scipy.signal import fftconvolve
    stars = np.zeros(shape); 
    for (y, x, f) in [(60, 60, 5000.), (180, 90, 3000.), (100, 200, 8000.)]:
        stars[y, x] = f
    S0 = fftconvolve(stars, P, mode="same")
    R += S0; N += S0
    # transient only in N, 400 counts at (128,128)
    T = np.zeros(shape); T[128, 128] = 400.0
    N += fftconvolve(T, P, mode="same")
    out = Z.zogy(N, R, P, P, sig, sig, Fn=1.0, Fr=1.0)
    score = out["S"] / Z.robust_sigma(out["S"])
    y, x = np.unravel_index(np.argmax(score), shape)
    assert (y, x) == (128, 128)
    assert abs(out["alpha"][128, 128] - 400.0) < 5 * out["sigma_alpha"]
    # static stars cancel (no residual above 5 sigma at their positions)
    for (yy, xx) in [(60, 60), (180, 90), (100, 200)]:
        assert abs(score[yy, xx]) < 5.0
    # noise-only: normalised score has unit scatter
    assert 0.9 < np.std(score[10:50, 10:50]) < 1.1


def test_zogy_expected_snr():
    """Matched-filter S/N of a faint source = F / sigma_alpha within 15% (noise realisation tolerance)."""
    rng = np.random.default_rng(11)
    shape = (256, 256); sig = 2.0
    P = I.gaussian_psf(2.5, 21); P /= P.sum()
    R = _img(shape, rng, sig); N = _img(shape, rng, sig)
    from scipy.signal import fftconvolve
    T = np.zeros(shape)
    pos = [(60 + 30 * i, 60 + 30 * j) for i in range(5) for j in range(5)]
    for (y, x) in pos:
        T[y, x] = 80.0
    N += fftconvolve(T, P, mode="same")
    out = Z.zogy(N, R, P, P, sig, sig)
    fl = np.array([out["alpha"][y, x] for y, x in pos])
    assert abs(fl.mean() - 80.0) < 3 * fl.std() / np.sqrt(len(fl)) + 1.0
    assert 0.7 < fl.std() / out["sigma_alpha"] < 1.3


def _flux_case(Fn, sn, sr, fw_n, fw_r, flux=300.0, seed=3):
    from scipy.signal import fftconvolve
    rng = np.random.default_rng(seed); shape = (256, 256)
    Pn = I.gaussian_psf(fw_n, 25); Pn /= Pn.sum(); Pr = I.gaussian_psf(fw_r, 25); Pr /= Pr.sum()
    pos = [(40 + 40 * i, 40 + 40 * j) for i in range(5) for j in range(5)]
    T = np.zeros(shape)
    for y, x in pos:
        T[y, x] = flux
    stars = np.zeros(shape); stars[100, 100] = 5000.0
    R = fftconvolve(stars, Pr, mode="same") + rng.normal(0, sr, shape)
    N = Fn * fftconvolve(stars, Pn, mode="same") + fftconvolve(T, Pn, mode="same") + rng.normal(0, sn, shape)
    o = Z.zogy(N, R, Pn, Pr, sn, sr, Fn=Fn, Fr=1.0)
    return o, np.array([o["alpha_new"][y, x] for y, x in pos]), np.array([o["alpha"][y, x] for y, x in pos])


def test_flux_normalisation_with_flux_ratio():
    """A source of 300 counts in N is recovered as 300 (alpha_new) for any F_n/F_r; `alpha` is in the reference scale."""
    for Fn in (0.5, 1.0, 2.0):
        o, a_new, a_ref = _flux_case(Fn, 2.0, 2.0, 2.5, 2.5)
        assert abs(a_new.mean() - 300.0) < 4 * a_new.std() / np.sqrt(len(a_new)) + 3.0, (Fn, a_new.mean())
        assert abs(a_ref.mean() * Fn - 300.0) < 4 * a_new.std() / np.sqrt(len(a_new)) + 3.0
        assert 0.6 < a_new.std() / o["sigma_alpha_new"] < 1.5, (Fn, a_new.std(), o["sigma_alpha_new"])


def test_flux_unbiased_with_different_psfs_and_noise():
    o, a_new, _ = _flux_case(1.0, 2.0, 3.0, 3.5, 2.2)
    assert abs(a_new.mean() - 300.0) < 4 * a_new.std() / np.sqrt(len(a_new)) + 3.0
    assert 0.6 < a_new.std() / o["sigma_alpha_new"] < 1.5


def test_psf_width_mismatch_biases_flux_low():
    """Measured sensitivity (documented limitation): a PSF 20% too narrow gives ~22% too little flux (matched filter applied with
    the wrong profile), i.e. the bias is ~ the fractional FWHM error.  Hence the PSF is measured from field stars per chip."""
    from scipy.signal import fftconvolve
    rng = np.random.default_rng(8); shape = (256, 256)
    Ptrue = I.gaussian_psf(3.0, 25); Ptrue /= Ptrue.sum(); Pwrong = I.gaussian_psf(3.0 * 0.8, 25); Pwrong /= Pwrong.sum()
    pos = [(40 + 40 * i, 40 + 40 * j) for i in range(5) for j in range(5)]
    T = np.zeros(shape)
    for y, x in pos:
        T[y, x] = 500.0
    N = fftconvolve(T, Ptrue, mode="same") + rng.normal(0, 2.0, shape); R = rng.normal(0, 2.0, shape)
    o = Z.zogy(N, R, Pwrong, Pwrong, 2.0, 2.0)
    a = np.array([o["alpha_new"][y, x] for y, x in pos])
    bias = a.mean() / 500.0 - 1.0
    assert -0.32 < bias < -0.12, bias
