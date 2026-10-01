from moving import transients as TR


def _det(i, ra, dec, f, ex):
    return dict(ra=ra, dec=dec, sign=1, cls="point", snr=20.0, file=f, ex=ex)


def test_static_grouping_host_association_and_cr_rejection(tmp_path):
    ra, de = 150.1000, 2.2000
    dets = [_det(0, ra, de, "a", 0), _det(1, ra + 0.1 / 3600, de, "b", 1),          # repeats on 2 files: kept
            _det(2, ra + 0.01, de, "a", 0)]                                         # single file: rejected (CR-like)
    g = TR.group_static(dets, tol_arcsec=0.4, min_epochs=2, snr_min=6)
    assert len(g) == 1 and sorted(g[0]) == [0, 1]
    cat = tmp_path / "c.tsv"
    cat.write_text("NUMBER\tALPHA_J2000\tDELTA_J2000\tFLUX_RADIUS\tz_phot\n7\t%.7f\t%.7f\t10\t0.5\n" % (ra, de + 1.0 / 3600))
    rows, cmap = TR.read_catalog(str(cat))
    h = TR.associate_host(ra, de, rows, cmap, 0.05)
    assert abs(h["sep_arcsec"] - 1.0) < 0.01 and h["z"] == 0.5 and str(h["id"]) == "7"
    cls, why = TR.heuristic_class(h, None)
    assert cls.startswith("SN-like")
    assert TR.heuristic_class(None, None)[0].startswith("hostless")
