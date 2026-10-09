"""Desktop catalog columns, image info, and forced photometry. No C extractor."""
import os
import sys

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "2")

from ogfmeas.sextract import detection_text, forced_text, info_text, main


def _gauss(shape, x0, y0, amp, sig, sky=0.0):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return sky + amp * np.exp(-0.5 * ((xx - x0) ** 2 + (yy - y0) ** 2) / sig ** 2)


def _write(path, image, **cards):
    from astropy.io import fits
    hdu = fits.PrimaryHDU(np.asarray(image, dtype=np.float32))
    for key, value in cards.items():
        hdu.header[key] = value
    hdu.writeto(path, overwrite=True)


def _table(text):
    lines = [line.split("\t") for line in text.strip().splitlines()]
    header = lines[0]
    return header, [dict(zip(header, row)) for row in lines[1:]]


def test_import_does_not_load_sep():
    before = "sep" in sys.modules
    import ogfmeas.sextract as sextract
    assert ("sep" in sys.modules) is before
    assert "ds9_sextract" not in sextract.__file__


def test_desktop_columns_and_position(tmp_path):
    image = _gauss((48, 48), 20.0, 22.0, 120.0, 1.8, sky=3.0)
    path = tmp_path / "one.fits"
    out = tmp_path / "cat.tsv"
    _write(path, image)
    assert main([str(path), "--out", str(out), "--detect-thresh", "8", "--back-size", "48",
                 "--mag-zeropoint", "25", "--conv-filter", "gauss5x5"]) == 0
    header, rows = _table(out.read_text(encoding="utf-8"))
    assert header[1] == "X_IMAGE" and "MAG_AUTO" in header and "CLASS_STAR" in header
    assert len(rows) == 1
    assert abs(float(rows[0]["X_IMAGE"]) - 21.0) < 0.6
    assert abs(float(rows[0]["Y_IMAGE"]) - 23.0) < 0.6
    assert float(rows[0]["MAG_AUTO"]) < 20.0
    assert rows[0]["CLASS_STAR"] == "0.000"
    direct = detection_text(image, thresh=8, bw=48, zp=25)
    assert "MAG_AUTO" in direct.splitlines()[0]


def test_info_prints_the_band_registry_keys(tmp_path):
    path = tmp_path / "wcs.fits"
    _write(path, np.zeros((8, 10)), FILTER="F814W", CRPIX1=5.0, CRPIX2=4.0,
           CRVAL1=150.0, CRVAL2=2.0, CD1_1=-1.0 / 3600.0, CD1_2=0.0, CD2_1=0.0, CD2_2=1.0 / 3600.0)
    text = info_text(str(path))
    values = dict(line.split("=", 1) for line in text.strip().splitlines())
    assert values["NAXIS1"] == "10" and values["NAXIS2"] == "8"
    assert values["FILTER"] == "F814W" and values["WCSVALID"] == "1"
    assert abs(float(values["PIXSCALE"]) - 1.0) < 1e-4
    assert "CRPIX1" in values and "CD1_1" in values
    assert main(["--info", str(path), "--out", str(tmp_path / "info.txt")]) == 0


def test_forced_photometry_keeps_the_source_and_rejects_empty_sky(tmp_path):
    image = _gauss((64, 64), 30.0, 28.0, 200.0, 2.0, sky=4.0)
    image = image + np.random.default_rng(1).normal(0.0, 0.4, image.shape)
    fits = tmp_path / "im.fits"
    cat = tmp_path / "det.tsv"
    _write(fits, image)
    cat.write_text(
        "NUMBER\tX_IMAGE\tY_IMAGE\tA_IMAGE\tB_IMAGE\tTHETA_IMAGE\tKRON_RADIUS\n"
        "1\t31.0000\t29.0000\t2.0\t2.0\t0\t3.5\n"
        "2\t4.0000\t4.0000\t2.0\t2.0\t0\t3.5\n",
        encoding="utf-8")
    text = forced_text(str(fits), str(cat), str(fits), band="R", zp=25.0, snr_min=8.0)
    _header, rows = _table(text)
    by_number = {row["NUMBER"]: row for row in rows}
    assert "MAG_AUTO_R" in rows[0] and "MAGERR_AUTO_R" in rows[0]
    assert float(by_number["1"]["MAG_AUTO_R"]) < 20.0
    assert float(by_number["2"]["MAG_AUTO_R"]) == 99.0
    assert main([str(fits), "--forced-catalog", str(cat), "--measure-image", str(fits),
                 "--band", "R", "--snr-min", "8", "--mag-zeropoint", "25",
                 "--out", str(tmp_path / "forced.tsv")]) == 0
