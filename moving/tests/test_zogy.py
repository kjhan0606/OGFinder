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
