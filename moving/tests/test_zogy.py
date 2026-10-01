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


# ------------------------------------------------------------------ source-noise, astrometric terms, spatially varying PSF
def _star_field(seed, shape=(256, 256), nstars=25, sky=5.0, shift=(0.0, 0.0), psf_fn=None, tflux=300.0, ntrans=16):
    """Poisson images of bright stars (static) + transients only in N.  psf_fn(x, y) -> PSF; returns N, R, transient pixel list."""
    from scipy import ndimage as ndi
    rng = np.random.default_rng(seed)
    psf_fn = psf_fn or (lambda x, y: I.gaussian_psf(2.6, 25))
    pos = [(rng.uniform(30, shape[0] - 30), rng.uniform(30, shape[1] - 30)) for _ in range(nstars)]
    fl = 10 ** rng.uniform(3.5, 5.0, nstars)
    def render(positions, fluxes, sh):
        img = np.zeros(shape)
        for (y, x), f in zip(positions, fluxes):
            P = psf_fn(x, y); r = P.shape[0] // 2; yi, xi = int(y), int(x)
            img[yi - r:yi + r + 1, xi - r:xi + r + 1] += f * ndi.shift(P, (sh[1] + y - yi, sh[0] + x - xi), order=3, mode="constant")
        return img
    Rt = render(pos, fl, (0, 0)); Nt = render(pos, fl, shift)
    tp = [(int(rng.uniform(40, shape[0] - 40)), int(rng.uniform(40, shape[1] - 40))) for _ in range(ntrans)]
    for y, x in tp:
        P = psf_fn(x, y); r = 12; Nt[y - r:y + r + 1, x - r:x + r + 1] += tflux * P
    N = rng.poisson(np.clip(Nt + sky, 0, None)).astype(float) - sky
    R = rng.poisson(np.clip(Rt + sky, 0, None)).astype(float) - sky
    return N, R, tp, np.sqrt(sky)


def _false_positives(score, tp, thr=5.0, rad=3, border=30):
    from scipy import ndimage as ndi
    mask = np.zeros(score.shape, bool)
    for y, x in tp:
        mask[y - rad:y + rad + 1, x - rad:x + rad + 1] = True
    s = score.copy(); s[:border] = 0; s[-border:] = 0; s[:, :border] = 0; s[:, -border:] = 0
    return ndi.label((np.abs(s) > thr) & ~mask)[1]


def test_new_keywords_are_optional_and_default_output_unchanged():
    rng = np.random.default_rng(1); P = I.gaussian_psf(2.5, 21)
    N = rng.normal(0, 2, (128, 128)); R = rng.normal(0, 2, (128, 128))
    o = Z.zogy(N, R, P, P, 2.0, 2.0)
    assert set(o) == {"D", "S", "alpha", "sigma_alpha", "alpha_new", "sigma_alpha_new", "PD", "FD", "sumP2"}
    o2 = Z.zogy(N, R, P, P, 2.0, 2.0, Vn=np.zeros_like(N), Vr=np.zeros_like(N))
    # zero source variance: the corrected score is the classic score / (F_D sqrt(sum P_D^2))
    ref = o["S"] / (o["FD"] * np.sqrt(o["sumP2"]))
    assert np.allclose(o2["S_corr"], ref, rtol=1e-6, atol=1e-8)


def test_source_noise_term_removes_bright_star_false_positives():
    """Poisson noise of bright stars gives hundreds of >5 sigma residuals with the sky-only normalisation; with Vn/Vr they vanish while the
    injected transients stay (>= 95 % recovered at 5 sigma)."""
    from scipy import ndimage as ndi
    fp_c = []; fp_s = []; rec = []
    for seed in range(3):
        N, R, tp, sg = _star_field(seed)
        P = I.gaussian_psf(2.6, 25)
        o = Z.zogy(N, R, P, P, sg, sg)
        s0 = o["S"] / Z.robust_sigma(o["S"]); fp_c.append(_false_positives(s0, tp))
        Vn = ndi.gaussian_filter(np.clip(N, 0, None), 1.0); Vr = ndi.gaussian_filter(np.clip(R, 0, None), 1.0)
        o2 = Z.zogy(N, R, P, P, sg, sg, Vn=Vn, Vr=Vr)
        s2 = o2["S_corr"] / Z.robust_sigma(o2["S_corr"]); fp_s.append(_false_positives(s2, tp))
        rec.append(np.mean([s2[y, x] > 5 for y, x in tp]))
    assert np.mean(fp_c) > 10 * max(np.mean(fp_s), 1.0)
    assert np.mean(rec) > 0.95


def test_flux_bias_unchanged_by_noise_terms():
    """The noise terms only change the variance map: alpha (flux) is identical, with a flux bias below 3 % for 300-count transients."""
    from scipy import ndimage as ndi
    N, R, tp, sg = _star_field(2)
    P = I.gaussian_psf(2.6, 25)
    o = Z.zogy(N, R, P, P, sg, sg)
    o2 = Z.zogy(N, R, P, P, sg, sg, Vn=ndi.gaussian_filter(np.clip(N, 0, None), 1.0), Vr=ndi.gaussian_filter(np.clip(R, 0, None), 1.0))
    assert np.allclose(o["alpha"], o2["alpha"])
    ratio = np.array([o["alpha"][y, x] for y, x in tp]) / 300.0
    assert abs(ratio.mean() - 1.0) < 0.05 and ratio.std() < 0.2


def test_astrometric_term_downweights_misregistered_star():
    """A star pair shifted by 0.3 px between N and R leaves a large dipole in S; the astrometric variance term (sigma = 0.3 px) lowers |S_corr| there."""
    from scipy import ndimage as ndi
    rng = np.random.default_rng(3); shape = (192, 192); P = I.gaussian_psf(2.6, 25)
    star = np.zeros(shape); star[96, 96] = 2.0e5
    Rt = ndi.gaussian_filter(star, 2.6 / 2.355) ; Nt = ndi.shift(Rt, (0.0, 0.3), order=3)
    sky = 4.0
    N = Nt + rng.normal(0, np.sqrt(sky), shape); R = Rt + rng.normal(0, np.sqrt(sky), shape)
    sg = np.sqrt(sky)
    o = Z.zogy(N, R, P, P, sg, sg)
    s_plain = np.abs(o["S"] / o["S"].std())[90:103, 90:103].max()
    o2 = Z.zogy(N, R, P, P, sg, sg, astrom_n=(0.3, 0.3), astrom_r=(0.0, 0.0))
    s_ast = np.abs(o2["S_corr"])[90:103, 90:103].max()
    s_plain_corr = np.abs(Z.zogy(N, R, P, P, sg, sg, Vn=np.zeros(shape), Vr=np.zeros(shape))["S_corr"])[90:103, 90:103].max()
    assert s_plain_corr > 50 and s_ast < 0.2 * s_plain_corr, (s_plain_corr, s_ast)
    # the term is local: far from the star the score is unchanged
    far = np.abs(o2["S_corr"] - Z.zogy(N, R, P, P, sg, sg, Vn=np.zeros(shape), Vr=np.zeros(shape))["S_corr"])[:40, :40].max()
    assert far < 0.2


def test_tiled_zogy_with_constant_psf_matches_global_zogy_in_the_interior():
    N, R, tp, sg = _star_field(4, ntrans=0)
    P = I.gaussian_psf(2.6, 25)
    g = Z.zogy(N, R, P, P, sg, sg)
    t = Z.zogy_tiled(N, R, P, P, sg, sg, tile=(128, 128), margin=32)
    sl = (slice(40, -40), slice(40, -40))
    assert np.allclose(g["alpha"][sl], t["alpha"][sl], atol=0.1 * g["sigma_alpha"] * 10)
    assert np.corrcoef(g["alpha"][sl].ravel(), t["alpha"][sl].ravel())[0, 1] > 0.99


def test_spatially_varying_psf_reduces_flux_bias_gradient():
    """FWHM 1.6 -> 3.6 px across the image: a single constant PSF biases the flux by ~ +5 % on the narrow side and ~ -8 % on the wide side
    (measured: +0.049 / -0.078); per-tile PSFs bring both below 4 %."""
    shape = (256, 512)
    def psf_fn(x, y):
        return I.gaussian_psf(1.6 + 2.0 * x / shape[1], 25)
    n = 4; t = 128
    arr = np.empty((2, n), object)
    for i in range(2):
        for j in range(n):
            arr[i, j] = psf_fn(t * j + t // 2, t * i + t // 2)
    F = I.PSFField(shape, (t, t), arr, np.zeros((2, n), int), np.ones((2, n), bool), np.full((2, n), 2.6), arr[0, 0], 2.6)
    left = {"const": [], "tiled": []}; right = {"const": [], "tiled": []}
    for seed in range(6):
        N, R, tp, sg = _star_field(seed, shape=shape, nstars=10, psf_fn=psf_fn, ntrans=24)
        Pc = psf_fn(shape[1] / 2, 0)
        a_c = Z.zogy(N, R, Pc, Pc, sg, sg)["alpha"]; a_t = Z.zogy_tiled(N, R, F, F, sg, sg, tile=(t, t), margin=24)["alpha"]
        for y, x in tp:
            side = left if x < 170 else right if x > 340 else None
            if side is not None:
                side["const"].append(a_c[y, x] / 300 - 1); side["tiled"].append(a_t[y, x] / 300 - 1)
    bc = max(abs(np.mean(left["const"])), abs(np.mean(right["const"])))
    bt = max(abs(np.mean(left["tiled"])), abs(np.mean(right["tiled"])))
    assert bt < 0.05 and bt < 0.7 * bc, (bc, bt)


def test_psf_field_measured_from_stars_falls_back_when_few_stars():
    from moving import imaging as Im
    P = I.gaussian_psf(2.5, 25)
    F = Im.constant_psf_field((300, 300), P, 2.5, tile=(100, 100))
    assert F.psfs.shape == (3, 3) and not F.varies and np.allclose(F.at(250, 250), P)
